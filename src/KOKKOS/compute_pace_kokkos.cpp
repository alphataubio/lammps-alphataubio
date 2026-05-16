// clang-format off
/* ----------------------------------------------------------------------
   compute pace/kk — Kokkos OpenMP descriptor compute for FitSNAP.
   See compute_pace_kokkos.h for design rationale.

   Key facts learned from compute_pace.cpp:
   - element_type_mapping(ik) = ik-1  (identity: LAMMPS type i → element i-1)
   - type_offsets is 1-indexed (type_offsets[itype], itype=1..ntypes)
   - projections[func_ind] = sum_ms sum_p ctildes[ms,p] * Re[prod(A)]
------------------------------------------------------------------------- */

#include "compute_pace_kokkos.h"

#include "ace-evaluator/ace_radial.h"
#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>

// ACECimpl redeclared here (body in compute_pace.cpp, same layout)
namespace LAMMPS_NS {
struct ACECimpl {
  ACECimpl() : basis_set(nullptr), ace(nullptr) {}
  ~ACECimpl() { delete basis_set; delete ace; }
  ACECTildeBasisSet  *basis_set;
  ACECTildeEvaluator *ace;
};
}

using namespace LAMMPS_NS;
using namespace MathConst;
// Y00 / sq3 / sq3o2 are class static constexpr in the header —
// do NOT redefine at file scope (ace_spherical_cart.h defines conflicting globals).

// ═════════════════════════════════════════════════════════════════════════════
// Constructor / Destructor
// ═════════════════════════════════════════════════════════════════════════════

template<class DeviceType>
ComputePACEKokkos<DeviceType>::ComputePACEKokkos(LAMMPS *lmp, int narg, char **arg)
    : ComputePACE(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read   = X_MASK | TYPE_MASK | TAG_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;

  if (dgradflag)
    error->all(FLERR, "compute pace/kk does not support dgradflag=1");
}

template<class DeviceType>
ComputePACEKokkos<DeviceType>::~ComputePACEKokkos()
{
  deallocate_views_of_views();
}

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::deallocate_views_of_views()
{
  if (!k_splines_gk.view_host().data()) return;
  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++) {
      k_splines_gk.view_host()(i,j).deallocate();
      k_splines_rnl.view_host()(i,j).deallocate();
      k_splines_hc.view_host()(i,j).deallocate();
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// init
// ═════════════════════════════════════════════════════════════════════════════

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::init()
{
  ComputePACE::init();

  auto basis_set = acecimpl->basis_set;
  nelements = basis_set->nelements;
  lmax      = basis_set->lmax;
  nradmax   = basis_set->nradmax;
  nradbase  = basis_set->nradbase;
  nfuncs_total = ncoeff;  // ncoeff set by ComputePACE constructor

  // element_type_mapping(ik) = ik-1 (identity mapping set in compute_pace.cpp)
  // so element index mu = LAMMPS atom type - 1; no d_map needed.

  // cutsq (conservative upper bound for all type pairs)
  int n = atom->ntypes + 1;
  MemKK::realloc_kokkos(k_cutsq, "pace:cutsq", n, n);
  d_cutsq = k_cutsq.template view<DeviceType>();
  {
    auto h = k_cutsq.view_host();
    double rcut = basis_set->cutoffmax;
    for (int i = 1; i < n; i++)
      for (int j = 1; j < n; j++)
        h(i,j) = rcut*rcut;
  }
  k_cutsq.modify_host();
  k_cutsq.sync_device();

  // type_offsets[itype] is 1-indexed in ComputePACE → d_type_offsets(mu) = type_offsets[mu+1]
  MemKK::realloc_kokkos(d_type_offsets, "pace:type_offsets", nelements);
  {
    auto h = Kokkos::create_mirror_view(d_type_offsets);
    for (int mu = 0; mu < nelements; mu++)
      h(mu) = type_offsets.at(mu + 1);
    Kokkos::deep_copy(d_type_offsets, h);
  }

  // Spherical harmonic coefficients
  MemKK::realloc_kokkos(d_idx_sph, "pace:idx_sph", (lmax+1)*(lmax+1));
  MemKK::realloc_kokkos(alm, "pace:alm", (lmax+1)*(lmax+1));
  MemKK::realloc_kokkos(blm, "pace:blm", (lmax+1)*(lmax+1));
  MemKK::realloc_kokkos(cl,  "pace:cl",  lmax+1);
  MemKK::realloc_kokkos(dl,  "pace:dl",  lmax+1);
  pre_compute_harmonics(lmax);

  copy_pertype();
  copy_splines();
  copy_tilde();

  // Request Kokkos-aware neighbour list
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(
      std::is_same_v<DeviceType,LMPHostType> &&
      !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(
      std::is_same_v<DeviceType,LMPDeviceType>);
}

// ─────────────────────────────────────────────────────────────────────────────
// grow
// ─────────────────────────────────────────────────────────────────────────────

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::grow(int natom, int maxn)
{
  auto basis_set = acecimpl->basis_set;

  if ((int)A.extent(0) < natom) {
    MemKK::realloc_kokkos(A_sph,         "pace:A_sph",         natom, nelements, idx_sph_max, nradmax+1);
    MemKK::realloc_kokkos(A,             "pace:A",             natom, nelements, (lmax+1)*(lmax+1), nradmax+1);
    MemKK::realloc_kokkos(A_rank1,       "pace:A_rank1",       natom, nelements, nradbase);
    MemKK::realloc_kokkos(A_list,        "pace:A_list",        natom, idx_ms_combs_max, basis_set->rankmax);
    MemKK::realloc_kokkos(A_forward_prod,"pace:A_forward_prod",natom, idx_ms_combs_max, basis_set->rankmax+1);
    MemKK::realloc_kokkos(dB_flatten,    "pace:dB_flatten",    natom, idx_ms_combs_max, basis_set->rankmax);
    MemKK::realloc_kokkos(rho_core,      "pace:rho_core",      natom);
  }

  if ((int)fr.extent(0) < natom || (int)fr.extent(1) < maxn) {
    MemKK::realloc_kokkos(fr,  "pace:fr",  natom, maxn, lmax+1, nradmax);
    MemKK::realloc_kokkos(dfr, "pace:dfr", natom, maxn, lmax+1, nradmax);
    MemKK::realloc_kokkos(gr,  "pace:gr",  natom, maxn, nradbase);
    MemKK::realloc_kokkos(dgr, "pace:dgr", natom, maxn, nradbase);
    const int maxfun = MAX(nradbase, nradmax*(lmax+1));
    MemKK::realloc_kokkos(d_values,      "pace:d_values",      natom, maxn, maxfun);
    MemKK::realloc_kokkos(d_derivatives, "pace:d_derivatives", natom, maxn, maxfun);
    // cr/dcr are 2D (natom × maxneigh) — hard-core repulsion
    MemKK::realloc_kokkos(cr,       "pace:cr",      natom, maxn);
    MemKK::realloc_kokkos(dcr,      "pace:dcr",     natom, maxn);
    MemKK::realloc_kokkos(d_ncount, "pace:ncount",  natom);
    MemKK::realloc_kokkos(d_mu,     "pace:mu",      natom, maxn);
    MemKK::realloc_kokkos(d_rhats,  "pace:rhats",   natom, maxn);
    MemKK::realloc_kokkos(d_rnorms, "pace:rnorms",  natom, maxn);
    MemKK::realloc_kokkos(d_nearest,"pace:nearest", natom, maxn);
  }

  if ((int)d_descriptors.extent(0) < natom ||
      (int)d_descriptors.extent(1) < nfuncs_total)
    MemKK::realloc_kokkos(d_descriptors, "pace:descriptors", natom, nfuncs_total);

  if ((int)d_desc_force_ij.extent(0) < natom ||
      (int)d_desc_force_ij.extent(1) < maxn)
    MemKK::realloc_kokkos(d_desc_force_ij, "pace:desc_force_ij",
                          natom, maxn, nfuncs_total, 3);
}

// ─────────────────────────────────────────────────────────────────────────────
// copy_pertype, copy_splines, copy_tilde, pre_compute_harmonics
// ─────────────────────────────────────────────────────────────────────────────

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::copy_pertype()
{
  auto basis_set = acecimpl->basis_set;
  MemKK::realloc_kokkos(d_ndensity, "pace:ndensity", nelements);
  auto h = Kokkos::create_mirror_view(d_ndensity);
  for (int n = 0; n < nelements; n++)
    h(n) = basis_set->map_embedding_specifications.at(n).ndensity;
  Kokkos::deep_copy(d_ndensity, h);
}

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::copy_splines()
{
  auto basis_set = acecimpl->basis_set;
  deallocate_views_of_views();

  k_splines_gk  = Kokkos::DualView<SplineInterpolatorKokkos**,DeviceType>("pace:splines_gk",  nelements, nelements);
  k_splines_rnl = Kokkos::DualView<SplineInterpolatorKokkos**,DeviceType>("pace:splines_rnl", nelements, nelements);
  k_splines_hc  = Kokkos::DualView<SplineInterpolatorKokkos**,DeviceType>("pace:splines_hc",  nelements, nelements);

  // Kokkos allocates DualView memory without calling constructors (device-design).
  // The SplineInterpolatorKokkos structs contain a Kokkos::View member whose
  // internal reference-count pointer is garbage.  operator=(SplineInterpolator&)
  // calls View::operator= which inspects the OLD reference count; garbage bits
  // → undefined behaviour → lookupTable.data() ends up null → crash in calcSplines.
  // Fix: zero the struct bytes first so m_record_bits == 0 and the decrement is skipped.
  {
    auto h_gk  = k_splines_gk.view_host();
    auto h_rnl = k_splines_rnl.view_host();
    auto h_hc  = k_splines_hc.view_host();
    for (int i = 0; i < nelements; i++)
      for (int j = 0; j < nelements; j++) {
        memset(&h_gk(i,j),  0, sizeof(SplineInterpolatorKokkos));
        memset(&h_rnl(i,j), 0, sizeof(SplineInterpolatorKokkos));
        memset(&h_hc(i,j),  0, sizeof(SplineInterpolatorKokkos));
      }
  }

  auto *rf = dynamic_cast<ACERadialFunctions*>(basis_set->radial_functions);
  if (!rf) error->all(FLERR,"compute pace/kk: radial basis style not supported");

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++) {
      k_splines_gk.view_host()(i,j)  = rf->splines_gk(i,j);
      k_splines_rnl.view_host()(i,j) = rf->splines_rnl(i,j);
      k_splines_hc.view_host()(i,j)  = rf->splines_hc(i,j);
    }
  k_splines_gk.modify_host();  k_splines_rnl.modify_host();  k_splines_hc.modify_host();
  k_splines_gk.sync_device();  k_splines_rnl.sync_device();  k_splines_hc.sync_device();
}

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::copy_tilde()
{
  auto basis_set = acecimpl->basis_set;
  idx_ms_combs_max = 0;
  int total_basis_size_max = 0;

  MemKK::realloc_kokkos(d_idx_ms_combs_count,"pace:idx_ms_combs_count",nelements);
  auto h_cnt = Kokkos::create_mirror_view(d_idx_ms_combs_count);

  for (int mu = 0; mu < nelements; mu++) {
    int cnt = 0;
    int sz1 = basis_set->total_basis_size_rank1[mu];
    int sz  = basis_set->total_basis_size[mu];
    ACECTildeBasisFunction *basis = basis_set->basis[mu];
    for (int f = 0; f < sz1; ++f) cnt++;
    for (int f = 0; f < sz;  ++f)
      for (int ms = 0; ms < basis[f].num_ms_combs; ++ms) cnt++;
    h_cnt(mu) = cnt;
    idx_ms_combs_max    = MAX(idx_ms_combs_max, cnt);
    total_basis_size_max = MAX(total_basis_size_max, sz1 + sz);
  }
  Kokkos::deep_copy(d_idx_ms_combs_count, h_cnt);

  MemKK::realloc_kokkos(d_rank,        "pace:rank",        nelements, total_basis_size_max);
  MemKK::realloc_kokkos(d_num_ms_combs,"pace:num_ms_combs",nelements, total_basis_size_max);
  MemKK::realloc_kokkos(d_idx_funcs,   "pace:idx_func",    nelements, idx_ms_combs_max);
  MemKK::realloc_kokkos(d_mus,  "pace:mus", nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ns,   "pace:ns",  nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ls,   "pace:ls",  nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ms_combs,"pace:ms_combs",nelements,idx_ms_combs_max,basis_set->rankmax);
  MemKK::realloc_kokkos(d_ctildes, "pace:ctildes", nelements,idx_ms_combs_max,basis_set->ndensitymax);

  auto h_rank         = Kokkos::create_mirror_view(d_rank);
  auto h_num_ms_combs = Kokkos::create_mirror_view(d_num_ms_combs);
  auto h_idx_funcs    = Kokkos::create_mirror_view(d_idx_funcs);
  auto h_mus          = Kokkos::create_mirror_view(d_mus);
  auto h_ns           = Kokkos::create_mirror_view(d_ns);
  auto h_ls           = Kokkos::create_mirror_view(d_ls);
  auto h_ms_combs     = Kokkos::create_mirror_view(d_ms_combs);
  auto h_ctildes      = Kokkos::create_mirror_view(d_ctildes);

  for (int mu = 0; mu < nelements; mu++) {
    int sz1 = basis_set->total_basis_size_rank1[mu];
    int sz  = basis_set->total_basis_size[mu];
    int ndensity = basis_set->map_embedding_specifications.at(mu).ndensity;
    ACECTildeBasisFunction *br1 = basis_set->basis_rank1[mu];
    ACECTildeBasisFunction *b   = basis_set->basis[mu];
    int mc = 0;

    for (int f = 0; f < sz1; ++f) {
      h_rank(mu,f) = 1;
      h_mus(mu,f,0) = br1[f].mus[0];
      h_ns(mu,f,0)  = br1[f].ns[0];
      for (int p = 0; p < ndensity; ++p) h_ctildes(mu,mc,p) = br1[f].ctildes[p];
      h_idx_funcs(mu,mc) = f;
      mc++;
    }
    for (int f = 0; f < sz; ++f) {
      int ft = sz1 + f;
      int rank = h_rank(mu,ft) = b[f].rank;
      h_num_ms_combs(mu,ft) = b[f].num_ms_combs;
      for (int t = 0; t < rank; t++) {
        h_mus(mu,ft,t) = b[f].mus[t];
        h_ns(mu,ft,t)  = b[f].ns[t];
        h_ls(mu,ft,t)  = b[f].ls[t];
      }
      for (int ms = 0; ms < b[f].num_ms_combs; ++ms) {
        auto ms_ptr = &b[f].ms_combs[ms*rank];
        for (int t = 0; t < rank; t++) h_ms_combs(mu,mc,t) = ms_ptr[t];
        for (int p = 0; p < ndensity; ++p)
          h_ctildes(mu,mc,p) = b[f].ctildes[ms*ndensity+p];
        h_idx_funcs(mu,mc) = ft;
        mc++;
      }
    }
  }

  Kokkos::deep_copy(d_rank, h_rank);
  Kokkos::deep_copy(d_num_ms_combs, h_num_ms_combs);
  Kokkos::deep_copy(d_idx_funcs,    h_idx_funcs);
  Kokkos::deep_copy(d_mus,      h_mus);
  Kokkos::deep_copy(d_ns,       h_ns);
  Kokkos::deep_copy(d_ls,       h_ls);
  Kokkos::deep_copy(d_ms_combs, h_ms_combs);
  Kokkos::deep_copy(d_ctildes,  h_ctildes);
}

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::pre_compute_harmonics(int lmax_)
{
  auto h_idx_sph = Kokkos::create_mirror_view(d_idx_sph);
  auto h_alm     = Kokkos::create_mirror_view(alm);
  auto h_blm     = Kokkos::create_mirror_view(blm);
  auto h_cl      = Kokkos::create_mirror_view(cl);
  auto h_dl      = Kokkos::create_mirror_view(dl);
  Kokkos::deep_copy(h_idx_sph, -1);

  int idx_sph = 0;
  for (int m = 0; m <= lmax_; m++) {
    const double msq = m*m;
    for (int l = m; l <= lmax_; l++) {
      h_idx_sph(l*(l+1)+m) = idx_sph;
      double a = 0.0, b = 0.0;
      if (l > 1 && l != m) {
        const double lsq = l*l, ld = 2*l;
        a = sqrt(double(4*lsq-1)/double(lsq-msq));
        b = -sqrt(double(lsq-ld+1-msq)/double(4*(lsq-ld+1)-1));
      }
      h_alm(idx_sph) = a;
      h_blm(idx_sph) = b;
      idx_sph++;
    }
  }
  idx_sph_max = idx_sph;
  for (int l = 1; l <= lmax_; l++) {
    h_cl(l) = -sqrt(1.0+0.5/double(l));
    h_dl(l) =  sqrt(double(2*(l-1)+3));
  }
  Kokkos::deep_copy(d_idx_sph, h_idx_sph);
  Kokkos::deep_copy(alm, h_alm);
  Kokkos::deep_copy(blm, h_blm);
  Kokkos::deep_copy(cl,  h_cl);
  Kokkos::deep_copy(dl,  h_dl);
}

// ═════════════════════════════════════════════════════════════════════════════
// compute_array
// ═════════════════════════════════════════════════════════════════════════════

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::compute_array()
{
  invoked_array = update->ntimestep;

  for (int r = 0; r < size_array_rows; r++)
    for (int c = 0; c < size_array_cols; c++)
      pace[r][c] = 0.0;

  atomKK->sync(HostKK, datamask_read);
  neighbor->build_one(list);

  auto *k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh  = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist     = k_list->d_ilist;
  inum = list->inum;

  x    = atomKK->k_x.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  maxneigh = 0;
  for (int ii = 0; ii < inum; ii++)
    maxneigh = MAX(maxneigh, d_numneigh(d_ilist(ii)));

  chunk_size   = MIN(512, inum);
  chunk_offset = 0;
  grow(chunk_size, maxneigh);

  while (chunk_offset < inum) {
    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    Kokkos::deep_copy(A_sph,          complex(0,0));
    Kokkos::deep_copy(A_rank1,        0.0);
    Kokkos::deep_copy(rho_core,       0.0);
    Kokkos::deep_copy(d_descriptors,  0.0);
    Kokkos::deep_copy(d_desc_force_ij,0.0);

    int vector_length = 1, team_size = 1;

    // ComputeNeigh
    {
      check_team_size_for<TagComputePACEComputeNeigh>(chunk_size, team_size, vector_length);
      int scratch = scratch_size_helper<int>(team_size * maxneigh);
      Kokkos::TeamPolicy<DeviceType,TagComputePACEComputeNeigh> pol(chunk_size,team_size,vector_length);
      Kokkos::parallel_for("ComputeNeigh", pol.set_scratch_size(0,Kokkos::PerTeam(scratch)), *this);
    }
    // ComputeRadial
    {
      int nt = ((chunk_size+team_size-1)/team_size)*maxneigh;
      check_team_size_for<TagComputePACEComputeRadial>(nt, team_size, vector_length);
      Kokkos::parallel_for("ComputeRadial",
          Kokkos::TeamPolicy<DeviceType,TagComputePACEComputeRadial>(nt,team_size,vector_length), *this);
    }
    // ComputeAi
    {
      int nt = ((chunk_size+team_size-1)/team_size)*maxneigh;
      check_team_size_for<TagComputePACEComputeAi>(nt, team_size, vector_length);
      Kokkos::parallel_for("ComputeAi",
          Kokkos::TeamPolicy<DeviceType,TagComputePACEComputeAi>(nt,team_size,vector_length), *this);
    }
    // ConjugateAi
    Kokkos::parallel_for("ConjugateAi",
        Kokkos::RangePolicy<DeviceType,TagComputePACEConjugateAi>(0,chunk_size), *this);

    // ComputeB
    Kokkos::parallel_for("ComputeB",
        Kokkos::RangePolicy<DeviceType,TagComputePACEComputeB>(0, chunk_size * idx_ms_combs_max), *this);

    // DescGrad
    Kokkos::parallel_for("DescGrad",
        Kokkos::RangePolicy<DeviceType,TagComputePACEDescGrad>(0, chunk_size * maxneigh), *this);

    Kokkos::fence();

    // Copy chunk to CPU pace[][] array
    auto h_desc   = Kokkos::create_mirror_view(d_descriptors);
    auto h_dforce = Kokkos::create_mirror_view(d_desc_force_ij);
    auto h_nearest = Kokkos::create_mirror_view(d_nearest);
    auto h_ncount  = Kokkos::create_mirror_view(d_ncount);
    Kokkos::deep_copy(h_desc,    d_descriptors);
    Kokkos::deep_copy(h_dforce,  d_desc_force_ij);
    Kokkos::deep_copy(h_nearest, d_nearest);
    Kokkos::deep_copy(h_ncount,  d_ncount);

    tagint *tag  = atom->tag;
    double **xa  = atom->x;

    for (int ii = 0; ii < chunk_size; ii++) {
      const int i      = d_ilist(ii + chunk_offset);
      const int ncount = h_ncount(ii);
      const int irow   = bik_rows + 3*(tag[i]-1);
      const int vrow0  = size_array_rows - 6;

      // Energy row
      for (int col = 0; col < nfuncs_total; col++)
        pace[0][col] += h_desc(ii, col);

      // Force and virial rows
      for (int jj = 0; jj < ncount; jj++) {
        const int j    = h_nearest(ii, jj);
        const int jrow = bik_rows + 3*(tag[j]-1);
        const double dx = xa[i][0] - xa[j][0];
        const double dy = xa[i][1] - xa[j][1];
        const double dz = xa[i][2] - xa[j][2];

        for (int col = 0; col < nfuncs_total; col++) {
          const double fx = h_dforce(ii, jj, col, 0);
          const double fy = h_dforce(ii, jj, col, 1);
          const double fz = h_dforce(ii, jj, col, 2);
          pace[irow+0][col] += fx;
          pace[irow+1][col] += fy;
          pace[irow+2][col] += fz;
          pace[jrow+0][col] -= fx;
          pace[jrow+1][col] -= fy;
          pace[jrow+2][col] -= fz;
          pace[vrow0+0][col] += fx*dx;
          pace[vrow0+1][col] += fy*dy;
          pace[vrow0+2][col] += fz*dz;
          pace[vrow0+3][col] += fz*dy;
          pace[vrow0+4][col] += fz*dx;
          pace[vrow0+5][col] += fy*dx;
        }
      }
    }
    chunk_offset += chunk_size;
  }

  // Normalise energy row by N_atoms
  for (int col = 0; col < nfuncs_total; col++)
    pace[0][col] /= natoms;

  MPI_Allreduce(&pace[0][0], &paceall[0][0],
                size_array_rows * size_array_cols,
                MPI_DOUBLE, MPI_SUM, world);

  for (int r = 0; r < bik_rows; r++) paceall[r][lastcol] = 0.0;
  paceall[0][lastcol] = c_pe->compute_scalar();

  double **f   = atom->f;
  tagint *tag2 = atom->tag;
  for (int i = 0; i < atom->nlocal; i++) {
    int irow = bik_rows + 3*(tag2[i]-1);
    paceall[irow++][lastcol] = f[i][0];
    paceall[irow++][lastcol] = f[i][1];
    paceall[irow  ][lastcol] = f[i][2];
  }

  c_virial->compute_vector();
  int vrow = 3*natoms + bik_rows;
  paceall[vrow++][lastcol] = c_virial->vector[0];
  paceall[vrow++][lastcol] = c_virial->vector[1];
  paceall[vrow++][lastcol] = c_virial->vector[2];
  paceall[vrow++][lastcol] = c_virial->vector[5];
  paceall[vrow++][lastcol] = c_virial->vector[4];
  paceall[vrow  ][lastcol] = c_virial->vector[3];
}

// ═════════════════════════════════════════════════════════════════════════════
// Kernels: ComputeNeigh, ComputeRadial, ComputeAi, ConjugateAi
// Verbatim from pair_pace_kokkos; only d_map(type(x)) → (type(x)-1)
// ═════════════════════════════════════════════════════════════════════════════

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::operator()(
    TagComputePACEComputeNeigh,
    const typename Kokkos::TeamPolicy<DeviceType,TagComputePACEComputeNeigh>::member_type& team) const
{
  const int ii   = team.league_rank();
  const int i    = d_ilist[ii + chunk_offset];
  const int jnum = d_numneigh[i];
  const KK_FLOAT xtmp=x(i,0), ytmp=x(i,1), ztmp=x(i,2);

  const int team_rank   = team.team_rank();
  int *inside = (int*)team.team_shmem().get_shmem(team.team_size()*maxneigh*sizeof(int),0)
                + team_rank*maxneigh;

  int ncount = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,jnum),
      [&](const int jj, int& cnt){
        int j = d_neighbors(i,jj) & NEIGHMASK;
        KK_FLOAT dx=xtmp-x(j,0), dy=ytmp-x(j,1), dz=ztmp-x(j,2);
        inside[jj] = (dx*dx+dy*dy+dz*dz < d_cutsq(type(i),type(j))) ? 1 : -1;
        if (inside[jj]>0) cnt++;
      }, ncount);
  d_ncount(ii) = ncount;

  Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,jnum),
      [&](const int jj, int& offset, bool final){
        if (inside[jj]<0) return;
        if (final) {
          int j = d_neighbors(i,jj) & NEIGHMASK;
          KK_FLOAT dx=xtmp-x(j,0), dy=ytmp-x(j,1), dz=ztmp-x(j,2);
          KK_FLOAT r=sqrt(dx*dx+dy*dy+dz*dz), ri=1.0/r;
          d_mu(ii,offset)      = type(j) - 1;   // element index = LAMMPS type - 1
          d_rnorms(ii,offset)  = r;
          d_rhats(ii,offset,0) = -dx*ri;
          d_rhats(ii,offset,1) = -dy*ri;
          d_rhats(ii,offset,2) = -dz*ri;
          d_nearest(ii,offset) = j;
        }
        offset++;
      });
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::operator()(
    TagComputePACEComputeRadial,
    const typename Kokkos::TeamPolicy<DeviceType,TagComputePACEComputeRadial>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size()*(team.league_rank()%
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int jj = team.league_rank()/((chunk_size+team.team_size()-1)/team.team_size());
  if (jj >= d_ncount(ii)) return;
  const int i    = d_ilist[ii+chunk_offset];
  const int mu_i = type(i) - 1;
  evaluate_splines(ii, jj, d_rnorms(ii,jj), nradbase, nradmax, mu_i, d_mu(ii,jj));
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::operator()(
    TagComputePACEComputeAi,
    const typename Kokkos::TeamPolicy<DeviceType,TagComputePACEComputeAi>::member_type& team) const
{
  int ii = team.team_rank() + team.team_size()*(team.league_rank()%
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;
  const int jj   = team.league_rank()/((chunk_size+team.team_size()-1)/team.team_size());
  if (jj >= d_ncount(ii)) return;
  const int mu_j = d_mu(ii,jj);

  for (int n = 0; n < nradbase; n++)
    Kokkos::atomic_add(&A_rank1(ii,mu_j,n), gr(ii,jj,n) * Y00);

  complex ylm, phase;
  const KK_FLOAT rx=d_rhats(ii,jj,0), ry=d_rhats(ii,jj,1), rz=d_rhats(ii,jj,2);
  phase.re=rx; phase.im=ry;
  KK_FLOAT plm=0,plm1=0,plm2=0;
  int idx_sph=0;

  // m=0
  for (int l=0; l<=lmax; l++) {
    if      (l==0) plm=Y00;
    else if (l==1) plm=Y00*sq3*rz;
    else           plm=alm(idx_sph)*(rz*plm1+blm(idx_sph)*plm2);
    ylm.re=plm; ylm.im=0;
    for (int n=0; n<nradmax; n++) {
      Kokkos::atomic_add(&A_sph(ii,mu_j,idx_sph,n).re, fr(ii,jj,l,n)*ylm.re);
      Kokkos::atomic_add(&A_sph(ii,mu_j,idx_sph,n).im, fr(ii,jj,l,n)*ylm.im);
    }
    plm2=plm1; plm1=plm; idx_sph++;
  }
  plm=plm1=plm2=0;
  // m=1
  for (int l=1; l<=lmax; l++) {
    if      (l==1) plm=-sq3o2*Y00;
    else if (l==2) plm=dl(l)*plm1*rz;
    else           plm=alm(idx_sph)*(rz*plm1+blm(idx_sph)*plm2);
    ylm=phase*plm;
    for (int n=0; n<nradmax; n++) {
      Kokkos::atomic_add(&A_sph(ii,mu_j,idx_sph,n).re, fr(ii,jj,l,n)*ylm.re);
      Kokkos::atomic_add(&A_sph(ii,mu_j,idx_sph,n).im, fr(ii,jj,l,n)*ylm.im);
    }
    plm2=plm1; plm1=plm; idx_sph++;
  }
  plm=plm1=plm2=0;
  KK_FLOAT plm_mm1=-sq3o2*Y00;
  complex phasem=phase;
  // m>1
  for (int m=2; m<=lmax; m++) {
    phasem=phasem*phase;
    for (int l=m; l<=lmax; l++) {
      if      (l==m)   { plm=cl(l)*plm_mm1; plm_mm1=plm; }
      else if (l==m+1) plm=dl(l)*plm_mm1*rz;
      else             plm=alm(idx_sph)*(rz*plm1+blm(idx_sph)*plm2);
      ylm.re=phasem.re*plm; ylm.im=phasem.im*plm;
      for (int n=0; n<nradmax; n++) {
        Kokkos::atomic_add(&A_sph(ii,mu_j,idx_sph,n).re, fr(ii,jj,l,n)*ylm.re);
        Kokkos::atomic_add(&A_sph(ii,mu_j,idx_sph,n).im, fr(ii,jj,l,n)*ylm.im);
      }
      plm2=plm1; plm1=plm; idx_sph++;
    }
  }
  Kokkos::atomic_add(&rho_core(ii), cr(ii,jj));
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::operator()(
    TagComputePACEConjugateAi, const int& ii) const
{
  for (int mu_j=0; mu_j<nelements; mu_j++) {
    int idx_sph=0;
    for (int m=0; m<=lmax; m++)
      for (int l=m; l<=lmax; l++) {
        for (int n=0; n<nradmax; n++) A(ii,mu_j,l*(l+1)+m,n)=A_sph(ii,mu_j,idx_sph,n);
        idx_sph++;
      }
    for (int l=0; l<=lmax; l++)
      for (int m=1; m<=l; m++) {
        const int is  = d_idx_sph(l*(l+1)+m);
        const int fac = (m%2==0)?1:-1;
        for (int n=0; n<nradmax; n++)
          A(ii,mu_j,l*(l+1)-m,n) = A_sph(ii,mu_j,is,n).conj()*(KK_FLOAT)fac;
      }
  }
}

// ═════════════════════════════════════════════════════════════════════════════
// ComputeB — accumulates B_k = sum_ms sum_p ctildes[p] * Re[prod(A)]
// into d_descriptors(ii, type_offset + idx_func).
// Also fills A_forward_prod and dB_flatten for DescGrad.
// ═════════════════════════════════════════════════════════════════════════════

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::operator()(
    TagComputePACEComputeB, const int& iter) const
{
  const int mc   = iter / chunk_size;
  const int ii   = iter % chunk_size;
  const int i    = d_ilist[ii + chunk_offset];
  const int mu_i = type(i) - 1;

  if (mc >= d_idx_ms_combs_count(mu_i)) return;

  const int idx_func  = d_idx_funcs(mu_i, mc);
  const int rank      = d_rank(mu_i, idx_func);
  const int col       = d_type_offsets(mu_i) + idx_func;
  const int ndensity  = d_ndensity(mu_i);

  if (rank == 1) {
    const int mu = d_mus(mu_i, idx_func, 0);
    const int n  = d_ns(mu_i,  idx_func, 0);
    const KK_FLOAT Aval = A_rank1(ii, mu, n-1);
    KK_FLOAT contrib = 0;
    for (int p = 0; p < ndensity; p++) contrib += d_ctildes(mu_i, mc, p) * Aval;
    Kokkos::atomic_add(&d_descriptors(ii, col), contrib);
  } else {
    // Forward product
    A_forward_prod(ii, mc, 0) = complex::one();
    for (int t = 0; t < rank; t++) {
      const int mu = d_mus(mu_i, idx_func, t);
      const int n  = d_ns(mu_i,  idx_func, t);
      const int l  = d_ls(mu_i,  idx_func, t);
      const int m  = d_ms_combs(mu_i, mc, t);
      A_list(ii, mc, t) = A(ii, mu, l*(l+1)+m, n-1);
      A_forward_prod(ii, mc, t+1) = A_forward_prod(ii, mc, t) * A_list(ii, mc, t);
    }
    // Backward product (fills dB_flatten for DescGrad)
    complex A_bwd = complex::one();
    for (int t = rank-1; t >= 1; t--) {
      dB_flatten(ii, mc, t) = A_forward_prod(ii, mc, t) * A_bwd;
      A_bwd = A_bwd * A_list(ii, mc, t);
    }
    dB_flatten(ii, mc, 0) = A_forward_prod(ii, mc, 0) * A_bwd;

    const complex B = A_forward_prod(ii, mc, rank);
    KK_FLOAT contrib = 0;
    for (int p = 0; p < ndensity; p++)
      contrib += B.real_part_product(d_ctildes(mu_i, mc, p));
    Kokkos::atomic_add(&d_descriptors(ii, col), contrib);
  }
}

// ═════════════════════════════════════════════════════════════════════════════
// DescGrad — gradient of B_k w.r.t. neighbour displacement → d_desc_force_ij
// ═════════════════════════════════════════════════════════════════════════════

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::operator()(
    TagComputePACEDescGrad, const int& iter) const
{
  const int ii = iter % chunk_size;
  const int jj = iter / chunk_size;
  if (jj >= d_ncount(ii)) return;

  const int i    = d_ilist[ii + chunk_offset];
  const int mu_i = type(i) - 1;
  const int mu_j = d_mu(ii, jj);
  const KK_FLOAT r    = d_rnorms(ii, jj);
  const KK_FLOAT rinv = 1.0/r;
  const KK_FLOAT rx = d_rhats(ii,jj,0);
  const KK_FLOAT ry = d_rhats(ii,jj,1);
  const KK_FLOAT rz = d_rhats(ii,jj,2);
  const int NRAD = nradmax;

  // Local stack: grad_phi for each (idx_sph, n) — host-only (no GPU VLA)
  KK_FLOAT gx[512]={}, gy[512]={}, gz[512]={};

  // Build grad_phi[is*NRAD+n] for all (l,m,n) — same as ComputeDerivative
  // in pair_pace_kokkos, but stored locally for the ms_combs loop below.
  complex ylm, dylm[3], phase, phasem, dyx, dyy, dyz, rdy;
  phase.re=rx; phase.im=ry;
  KK_FLOAT plm=0,plm1=0,plm2=0,dplm=0,dplm1=0,dplm2=0;
  int idx_sph=0;

  // m=0
  for (int l=0; l<=lmax; l++) {
    if      (l==0) { plm=Y00;dplm=0; }
    else if (l==1) { plm=Y00*sq3*rz; dplm=Y00*sq3; }
    else { plm=alm(idx_sph)*(rz*plm1+blm(idx_sph)*plm2);
           dplm=alm(idx_sph)*(plm1+rz*dplm1+blm(idx_sph)*dplm2); }
    ylm.re=plm; ylm.im=0;
    dyz.re=dplm; rdy.re=dyz.re*rz;
    dylm[0].re=-rdy.re*rx; dylm[0].im=0;
    dylm[1].re=-rdy.re*ry; dylm[1].im=0;
    dylm[2].re=dyz.re-rdy.re*rz; dylm[2].im=0;
    for (int n=0; n<NRAD; n++) {
      KK_FLOAT Ror=fr(ii,jj,l,n)*rinv, DR=dfr(ii,jj,l,n);
      int gi=idx_sph*NRAD+n;
      gx[gi]=ylm.re*DR*rx+dylm[0].re*Ror;
      gy[gi]=ylm.re*DR*ry+dylm[1].re*Ror;
      gz[gi]=ylm.re*DR*rz+dylm[2].re*Ror;
    }
    plm2=plm1; plm1=plm; dplm2=dplm1; dplm1=dplm; idx_sph++;
  }
  plm=plm1=plm2=dplm=dplm1=dplm2=0;
  // m=1
  for (int l=1; l<=lmax; l++) {
    if      (l==1) { plm=-sq3o2*Y00; dplm=0; }
    else if (l==2) { KK_FLOAT t=dl(l)*plm1; plm=t*rz; dplm=t; }
    else { plm=alm(idx_sph)*(rz*plm1+blm(idx_sph)*plm2);
           dplm=alm(idx_sph)*(plm1+rz*dplm1+blm(idx_sph)*dplm2); }
    ylm=phase*plm;
    dyx.re=plm; dyx.im=0; dyy.re=0; dyy.im=plm; dyz=phase*dplm;
    rdy.re=rx*dyx.re+rz*dyz.re; rdy.im=ry*dyy.im+rz*dyz.im;
    dylm[0].re=dyx.re-rdy.re*rx; dylm[0].im=-rdy.im*rx;
    dylm[1].re=-rdy.re*ry;       dylm[1].im=dyy.im-rdy.im*ry;
    dylm[2].re=dyz.re-rdy.re*rz; dylm[2].im=dyz.im-rdy.im*rz;
    for (int n=0; n<NRAD; n++) {
      KK_FLOAT Ror=fr(ii,jj,l,n)*rinv*2.0, DR=dfr(ii,jj,l,n)*2.0;
      complex Y_DR=ylm*DR;
      int gi=idx_sph*NRAD+n;
      gx[gi]=Y_DR.re*rx+dylm[0].re*Ror;
      gy[gi]=Y_DR.re*ry+dylm[1].re*Ror;
      gz[gi]=Y_DR.re*rz+dylm[2].re*Ror;
    }
    plm2=plm1; plm1=plm; dplm2=dplm1; dplm1=dplm; idx_sph++;
  }
  plm=plm1=plm2=dplm=dplm1=dplm2=0;
  KK_FLOAT plm_mm1=-sq3o2*Y00;
  phasem=phase;
  // m>1
  for (int m=2; m<=lmax; m++) {
    phasem=phasem*phase;
    for (int l=m; l<=lmax; l++) {
      if      (l==m)   { plm=cl(l)*plm_mm1; dplm=0; plm_mm1=plm; }
      else if (l==m+1) { KK_FLOAT t=dl(l)*plm_mm1; plm=t*rz; dplm=t; }
      else { plm=alm(idx_sph)*(rz*plm1+blm(idx_sph)*plm2);
             dplm=alm(idx_sph)*(plm1+rz*dplm1+blm(idx_sph)*dplm2); }
      ylm.re=phasem.re*plm; ylm.im=phasem.im*plm;
      dyx=complex{phasem.re*KK_FLOAT(m)*plm, phasem.im*KK_FLOAT(m)*plm};
      dyy.re=-dyx.im; dyy.im=dyx.re;
      dyz=phasem*dplm;
      rdy.re=rx*dyx.re+ry*dyy.re+rz*dyz.re;
      rdy.im=rx*dyx.im+ry*dyy.im+rz*dyz.im;
      dylm[0].re=dyx.re-rdy.re*rx; dylm[0].im=dyx.im-rdy.im*rx;
      dylm[1].re=dyy.re-rdy.re*ry; dylm[1].im=dyy.im-rdy.im*ry;
      dylm[2].re=dyz.re-rdy.re*rz; dylm[2].im=dyz.im-rdy.im*rz;
      for (int n=0; n<NRAD; n++) {
        KK_FLOAT Ror=fr(ii,jj,l,n)*rinv*2.0, DR=dfr(ii,jj,l,n)*2.0;
        complex Y_DR=ylm*DR;
        int gi=idx_sph*NRAD+n;
        gx[gi]=Y_DR.re*rx+dylm[0].re*Ror;
        gy[gi]=Y_DR.re*ry+dylm[1].re*Ror;
        gz[gi]=Y_DR.re*rz+dylm[2].re*Ror;
      }
      plm2=plm1; plm1=plm; dplm2=dplm1; dplm1=dplm; idx_sph++;
    }
  }

  // Loop over ms_combs and accumulate per-function gradient
  const int ndensity = d_ndensity(mu_i);
  for (int mc = 0; mc < d_idx_ms_combs_count(mu_i); mc++) {
    const int idx_func = d_idx_funcs(mu_i, mc);
    const int rank     = d_rank(mu_i, idx_func);
    const int col      = d_type_offsets(mu_i) + idx_func;

    if (rank == 1) {
      if (d_mus(mu_i, idx_func, 0) != mu_j) continue;
      const int n  = d_ns(mu_i, idx_func, 0);
      KK_FLOAT ctot = 0;
      for (int p = 0; p < ndensity; p++) ctot += d_ctildes(mu_i, mc, p);
      const KK_FLOAT DGR = ctot * dgr(ii,jj,n-1) * Y00;
      d_desc_force_ij(ii,jj,col,0) += DGR * rx;
      d_desc_force_ij(ii,jj,col,1) += DGR * ry;
      d_desc_force_ij(ii,jj,col,2) += DGR * rz;
    } else {
      for (int t = 0; t < rank; t++) {
        if (d_mus(mu_i, idx_func, t) != mu_j) continue;
        const int n_t = d_ns(mu_i, idx_func, t);
        const int l_t = d_ls(mu_i, idx_func, t);
        const int m_t = d_ms_combs(mu_i, mc, t);
        const int is  = d_idx_sph(l_t*(l_t+1)+m_t);
        if (is < 0) continue;
        const int gi = is * NRAD + (n_t-1);
        const complex dB = dB_flatten(ii, mc, t);
        KK_FLOAT w = 0;
        for (int p = 0; p < ndensity; p++) w += d_ctildes(mu_i, mc, p) * dB.re;
        d_desc_force_ij(ii,jj,col,0) += w * gx[gi];
        d_desc_force_ij(ii,jj,col,1) += w * gy[gi];
        d_desc_force_ij(ii,jj,col,2) += w * gz[gi];
      }
    }
  }
}

// ═════════════════════════════════════════════════════════════════════════════
// SplineInterpolatorKokkos
// ═════════════════════════════════════════════════════════════════════════════

template<class DeviceType>
void ComputePACEKokkos<DeviceType>::SplineInterpolatorKokkos::operator=(
    const SplineInterpolator &spline)
{
  cutoff=spline.cutoff; deltaSplineBins=spline.deltaSplineBins;
  ntot=spline.ntot; nlut=spline.nlut;
  invrscalelookup=spline.invrscalelookup; rscalelookup=spline.rscalelookup;
  num_of_functions=spline.num_of_functions;
  // lookupTable shape: (ntot+1, num_of_functions) with compile-time last dim [4]
  lookupTable = t_ace_3d4_lr("lookupTable", ntot+1, num_of_functions);
  auto h = Kokkos::create_mirror_view(lookupTable);
  for (int i=0; i<ntot+1; i++)
    for (int j=0; j<num_of_functions; j++)
      for (int k=0; k<4; k++)
        h(i,j,k) = spline.lookupTable(i,j,k);
  Kokkos::deep_copy(lookupTable, h);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::SplineInterpolatorKokkos::calcSplines(
    int ii, int jj, KK_FLOAT r,
    const t_ace_3d &vals, const t_ace_3d &derivs) const
{
  KK_FLOAT x = r * rscalelookup;
  int nl = static_cast<int>(floor(x));
  if (nl <= 0) Kokkos::abort("Very small distance in pace/kk");
  if (nl < nlut) {
    KK_FLOAT wl=x-KK_FLOAT(nl), wl2=wl*wl, wl3=wl2*wl;
    for (int f=0; f<num_of_functions; f++) {
      KK_FLOAT c0=lookupTable(nl,f,0), c1=lookupTable(nl,f,1),
               c2=lookupTable(nl,f,2), c3=lookupTable(nl,f,3);
      vals(ii,jj,f)   = c0+c1*wl+c2*wl2+c3*wl3;
      derivs(ii,jj,f) = (c1+c2*2*wl+c3*3*wl2)*rscalelookup;
    }
  } else {
    for (int f=0; f<num_of_functions; f++) { vals(ii,jj,f)=0; derivs(ii,jj,f)=0; }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// evaluate_splines, cutoff_func_poly, check_team_size_for, scratch_size_helper
// ─────────────────────────────────────────────────────────────────────────────

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::evaluate_splines(
    int ii, int jj, KK_FLOAT r, int, int, int mu_i, int mu_j) const
{
  k_splines_gk.template view<DeviceType>()(mu_i,mu_j).calcSplines(ii,jj,r,gr,dgr);
  k_splines_rnl.template view<DeviceType>()(mu_i,mu_j).calcSplines(ii,jj,r,d_values,d_derivatives);
  for (int l=0; l<(int)fr.extent(2); l++)
    for (int k=0; k<(int)fr.extent(3); k++) {
      const int flat=k*fr.extent(2)+l;
      fr(ii,jj,l,k)  = d_values(ii,jj,flat);
      dfr(ii,jj,l,k) = d_derivatives(ii,jj,flat);
    }
  k_splines_hc.template view<DeviceType>()(mu_i,mu_j).calcSplines(ii,jj,r,d_values,d_derivatives);
  cr(ii,jj)  = d_values(ii,jj,0);
  dcr(ii,jj) = d_derivatives(ii,jj,0);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputePACEKokkos<DeviceType>::cutoff_func_poly(
    KK_FLOAT r, KK_FLOAT r_in, KK_FLOAT delta, KK_FLOAT &fc, KK_FLOAT &dfc) const
{
  if (r <= r_in-delta) { fc=1; dfc=0; }
  else if (r >= r_in)  { fc=0; dfc=0; }
  else {
    KK_FLOAT x=1-2*(1+(r-r_in)/delta);
    fc  = 0.5+7.5/2.*(x/4.-x*x*x/6.+x*x*x*x*x/20.);
    dfc = -7.5/delta*(0.25-x*x/2.+x*x*x*x/4.);
  }
}

template<class DeviceType>
template<class TagStyle>
void ComputePACEKokkos<DeviceType>::check_team_size_for(
    int inum, int &team_size, int vector_length)
{
  int mx = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO)
               .team_size_max(*this, Kokkos::ParallelForTag());
  if (team_size*vector_length > mx) team_size = mx/vector_length;
}

template<class DeviceType>
template<typename scratch_type>
int ComputePACEKokkos<DeviceType>::scratch_size_helper(int values_per_team)
{
  typedef Kokkos::View<scratch_type*,
      Kokkos::DefaultExecutionSpace::scratch_memory_space,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>> ScratchView;
  return ScratchView::shmem_size(values_per_team);
}

// ─────────────────────────────────────────────────────────────────────────────

namespace LAMMPS_NS {
  template class ComputePACEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
  template class ComputePACEKokkos<LMPHostType>;
#endif
}
