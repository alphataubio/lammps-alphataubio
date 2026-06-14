// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors (pair style dispersion/d3):
      Sonia Salomoni, Arthur France-Lanord
   KOKKOS port:
      Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "pair_dispersion_d3_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

// global ad hoc parameters (same as the host pair style)

static constexpr double K1 = 16.0;
static constexpr double K3 = -4.0;
static constexpr double AUTOANG = 0.52917725;    // Bohr -> Angstrom
static constexpr double AUTOEV = 27.21140795;    // Hartree -> eV

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDispersionD3Kokkos<DeviceType>::PairDispersionD3Kokkos(LAMMPS *lmp) : PairDispersionD3(lmp)
{
  respa_enable = 0;
  single_enable = 0;

  kokkosable = 1;
  reverse_comm_device = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairDispersionD3Kokkos<DeviceType>::~PairDispersionD3Kokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ----------------------------------------------------------------------
   build device copies of the per-type and per-type-pair coefficient
   tables that the host coeff()/read_*() routines have already filled
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::allocate_tables_kokkos()
{
  const int n = atom->ntypes;

  k_r2r4 = DAT::tdual_kkfloat_1d("pair:r2r4",n+1);
  k_rcov = DAT::tdual_kkfloat_1d("pair:rcov",n+1);
  k_mxci = DAT::tdual_int_1d("pair:mxci",n+1);
  k_r0ab = DAT::tdual_kkfloat_2d("pair:r0ab",n+1,n+1);
  k_c6ab = tdual_c6_5d("pair:c6ab",n+1,n+1,5,5,3);

  auto hv_r2r4 = k_r2r4.view_host();
  auto hv_rcov = k_rcov.view_host();
  auto hv_mxci = k_mxci.view_host();
  auto hv_r0ab = k_r0ab.view_host();
  auto hv_c6ab = k_c6ab.view_host();

  for (int i = 0; i <= n; i++) {
    hv_r2r4(i) = static_cast<KK_FLOAT>(r2r4[i]);
    hv_rcov(i) = static_cast<KK_FLOAT>(rcov[i]);
    hv_mxci(i) = mxci[i];
    for (int j = 0; j <= n; j++) {
      hv_r0ab(i,j) = static_cast<KK_FLOAT>(r0ab[i][j]);
      for (int ci = 0; ci < 5; ci++)
        for (int cj = 0; cj < 5; cj++)
          for (int k = 0; k < 3; k++)
            hv_c6ab(i,j,ci,cj,k) = static_cast<KK_FLOAT>(c6ab[i][j][ci][cj][k]);
    }
  }

  k_r2r4.modify_host(); k_r2r4.template sync<DeviceType>(); d_r2r4 = k_r2r4.template view<DeviceType>();
  k_rcov.modify_host(); k_rcov.template sync<DeviceType>(); d_rcov = k_rcov.template view<DeviceType>();
  k_mxci.modify_host(); k_mxci.template sync<DeviceType>(); d_mxci = k_mxci.template view<DeviceType>();
  k_r0ab.modify_host(); k_r0ab.template sync<DeviceType>(); d_r0ab = k_r0ab.template view<DeviceType>();
  k_c6ab.modify_host(); k_c6ab.template sync<DeviceType>(); d_c6ab = k_c6ab.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::init_style()
{
  // checks atom IDs, sets metal units, requests a (half) neighbor list

  PairDispersionD3::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  // the host coeff() has already filled the coefficient tables: mirror to device

  allocate_tables_kokkos();
}

/* ----------------------------------------------------------------------
   Compute : energy, force and stress.

   Mirrors the three-stage host algorithm, where cn is the analog of the
   EAM density rho and dc6 is the analog of the EAM embedding derivative fp:
     A) accumulate coordination number cn, reverse+forward communicate
     B) dispersion energy/forces; accumulate dc6 = dE/dCN, reverse+forward comm
     C) coordination-number-derivative contribution to the forces
------------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  // grow cn / dc6 arrays if necessary (need to be atom->nmax in length)

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    k_cn = DAT::tdual_kkfloat_1d("pair:cn",nmax);
    k_dc6 = DAT::tdual_kkfloat_1d("pair:dc6",nmax);
    d_cn = k_cn.template view<DeviceType>();
    d_dc6 = k_dc6.template view<DeviceType>();
    h_cn = k_cn.view_host();
    h_dc6 = k_dc6.view_host();
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  // duplicated/atomic scatter views for cn, dc6, forces and per-atom e/virial

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_cn    = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterDuplicated>(d_cn);
    dup_dc6   = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterDuplicated>(d_dc6);
    dup_f     = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterDuplicated>(d_vatom);
  } else {
    ndup_cn    = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterNonDuplicated>(d_cn);
    ndup_dc6   = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterNonDuplicated>(d_dc6);
    ndup_f     = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<KKScatterSum, KKScatterNonDuplicated>(d_vatom);
  }

  // compile/run-time constants and parameters, cast to the working precision

  const KK_FLOAT k1       = static_cast<KK_FLOAT>(K1);
  const KK_FLOAT k3       = static_cast<KK_FLOAT>(K3);
  const KK_FLOAT autoang  = static_cast<KK_FLOAT>(AUTOANG);
  const KK_FLOAT autoang2 = autoang*autoang;
  const KK_FLOAT c6_conv  = static_cast<KK_FLOAT>(AUTOEV)*autoang2*autoang2*autoang2;    // autoev*autoang^6
  const KK_FLOAT l_rthr   = static_cast<KK_FLOAT>(rthr);
  const KK_FLOAT l_cn_thr = static_cast<KK_FLOAT>(cn_thr);
  const int      l_damping = dampingCode;
  const KK_FLOAT l_s6     = static_cast<KK_FLOAT>(s6);
  const KK_FLOAT l_s8     = static_cast<KK_FLOAT>(s8);
  const KK_FLOAT l_rs6    = static_cast<KK_FLOAT>(rs6);
  const KK_FLOAT l_rs8    = static_cast<KK_FLOAT>(rs8);
  const KK_FLOAT l_a1     = static_cast<KK_FLOAT>(a1);
  const KK_FLOAT l_a2     = static_cast<KK_FLOAT>(a2);
  const KK_FLOAT l_alpha  = static_cast<KK_FLOAT>(alpha);

  const int l_nlocal       = nlocal;
  const int l_eflag        = eflag;
  const int l_eflag_global = eflag_global;
  const int l_eflag_atom   = eflag_atom;
  const int l_vflag_either = vflag_either;
  const int l_vflag_global = vflag_global;
  const int l_vflag_atom   = vflag_atom;

  copymode = 1;

  // zero the base copies of the cn / dc6 scatter views

  Kokkos::deep_copy(d_cn, KK_ZERO);
  Kokkos::deep_copy(d_dc6, KK_ZERO);

  // ===================================================================
  // Pass A : coordination number
  // ===================================================================

  auto launch_cn = [&]<int NEIGHFLAG, int NEWTON_PAIR>() {
    auto l_x = x; auto l_type = type;
    auto l_ilist = d_ilist; auto l_numneigh = d_numneigh; auto l_neighbors = d_neighbors;
    auto l_rcov = d_rcov;

    auto v_cn = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_cn),decltype(ndup_cn)>::get(dup_cn,ndup_cn);
    auto a_cn = v_cn.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    Kokkos::parallel_for("PairDispersionD3:cn",
      Kokkos::RangePolicy<DeviceType>(0,inum),
      KOKKOS_LAMBDA(const int ii) {
        const int i = l_ilist(ii);
        const KK_FLOAT xtmp = l_x(i,0);
        const KK_FLOAT ytmp = l_x(i,1);
        const KK_FLOAT ztmp = l_x(i,2);
        const int itype = l_type(i);
        const int jnum = l_numneigh(i);

        KK_ACC_FLOAT cnitmp = 0.0;

        for (int jj = 0; jj < jnum; jj++) {
          int j = l_neighbors(i,jj);
          j &= NEIGHMASK;
          const KK_FLOAT delx = xtmp - l_x(j,0);
          const KK_FLOAT dely = ytmp - l_x(j,1);
          const KK_FLOAT delz = ztmp - l_x(j,2);
          const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < l_cn_thr) {
            const int jtype = l_type(j);
            const KK_FLOAT rr = sqrt(rsq);
            const KK_FLOAT rcovij = (l_rcov(itype) + l_rcov(jtype))*autoang;
            const KK_FLOAT cn_ij = static_cast<KK_FLOAT>(1.0) /
              (static_cast<KK_FLOAT>(1.0) + exp(-k1*((rcovij/rr) - static_cast<KK_FLOAT>(1.0))));
            cnitmp += cn_ij;
            if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < l_nlocal))
              a_cn(j) += cn_ij;
          }
        }
        a_cn(i) += static_cast<KK_FLOAT>(cnitmp);
      });
  };

  if (neighflag == HALF) {
    if (newton_pair) launch_cn.template operator()<HALF,1>();
    else             launch_cn.template operator()<HALF,0>();
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair) launch_cn.template operator()<HALFTHREAD,1>();
    else             launch_cn.template operator()<HALFTHREAD,0>();
  } else {
    if (newton_pair) launch_cn.template operator()<FULL,1>();
    else             launch_cn.template operator()<FULL,0>();
  }

  if (need_dup) Kokkos::Experimental::contribute(d_cn, dup_cn);

  // communicate coordination number (reverse to sum, forward to distribute)

  communicationStage = 1;
  if (newton_pair && (neighflag == HALF || neighflag == HALFTHREAD)) {
    k_cn.template modify<DeviceType>();
    comm->reverse_comm(this);
    k_cn.template sync<DeviceType>();
  }
  k_cn.template modify<DeviceType>();
  comm->forward_comm(this);
  k_cn.template sync<DeviceType>();

  // ===================================================================
  // Pass B : dispersion energy / forces and dc6 = dE/dCN accumulation
  // ===================================================================

  auto launch_force = [&]<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>(EV_FLOAT &ev_out) {
    auto l_x = x; auto l_type = type;
    auto l_ilist = d_ilist; auto l_numneigh = d_numneigh; auto l_neighbors = d_neighbors;
    auto l_cn = d_cn;
    auto l_r2r4 = d_r2r4; auto l_r0ab = d_r0ab; auto l_mxci = d_mxci; auto l_c6ab = d_c6ab;

    auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
    auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
    auto v_dc6 = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_dc6),decltype(ndup_dc6)>::get(dup_dc6,ndup_dc6);
    auto a_dc6 = v_dc6.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
    auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
    auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
    auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
    auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    const KK_FLOAT ls0 = special_lj[0], ls1 = special_lj[1], ls2 = special_lj[2], ls3 = special_lj[3];

    Kokkos::parallel_reduce("PairDispersionD3:force",
      Kokkos::RangePolicy<DeviceType>(0,inum),
      KOKKOS_LAMBDA(const int ii, EV_FLOAT &ev) {
        const int i = l_ilist(ii);
        const KK_FLOAT xtmp = l_x(i,0);
        const KK_FLOAT ytmp = l_x(i,1);
        const KK_FLOAT ztmp = l_x(i,2);
        const int itype = l_type(i);
        const KK_FLOAT cni = l_cn(i);
        const KK_FLOAT r2r4i = l_r2r4(itype);
        const int jnum = l_numneigh(i);
        const int mxi = l_mxci(itype);

        KK_ACC_FLOAT fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;
        KK_ACC_FLOAT dc6itmp = 0.0;

        for (int jj = 0; jj < jnum; jj++) {
          int jraw = l_neighbors(i,jj);
          const int sb = (jraw >> SBBITS) & 3;
          const KK_FLOAT factor_lj = (sb==0) ? ls0 : ((sb==1) ? ls1 : ((sb==2) ? ls2 : ls3));
          int j = jraw & NEIGHMASK;
          const KK_FLOAT delx = xtmp - l_x(j,0);
          const KK_FLOAT dely = ytmp - l_x(j,1);
          const KK_FLOAT delz = ztmp - l_x(j,2);
          const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < l_rthr) {
            const int jtype = l_type(j);
            const KK_FLOAT cnj = l_cn(j);
            const KK_FLOAT r = sqrt(rsq);
            const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
            const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
            const KK_FLOAT r8inv = r6inv*r2inv;
            const KK_FLOAT r10inv = r8inv*r2inv;

            // --- get_dC6: Gaussian-weighted interpolation over the C6 grid ---

            KK_FLOAT C6 = 0.0, dC6i = 0.0, dC6j = 0.0;
            {
              const int mxj = l_mxci(jtype);
              KK_FLOAT c6mem = static_cast<KK_FLOAT>(-1.0e20);
              KK_FLOAT r_save = static_cast<KK_FLOAT>(1.0e20);
              KK_FLOAT num = 0.0, den = 0.0;
              KK_FLOAT d_num_i = 0.0, d_num_j = 0.0, d_den_i = 0.0, d_den_j = 0.0;

              for (int ci = 0; ci <= mxi; ci++) {
                for (int cj = 0; cj <= mxj; cj++) {
                  const KK_FLOAT c6_ref = l_c6ab(itype,jtype,ci,cj,0)*c6_conv;
                  if (c6_ref > static_cast<KK_FLOAT>(0.0)) {
                    const KK_FLOAT cni_ref = l_c6ab(itype,jtype,ci,cj,1);
                    const KK_FLOAT cnj_ref = l_c6ab(itype,jtype,ci,cj,2);
                    const KK_FLOAT dcni = cni - cni_ref;
                    const KK_FLOAT dcnj = cnj - cnj_ref;
                    const KK_FLOAT rr = dcni*dcni + dcnj*dcnj;
                    if (rr < r_save) { r_save = rr; c6mem = c6_ref; }
                    KK_FLOAT expterm = exp(k3*rr);
                    num += c6_ref*expterm;
                    den += expterm;
                    expterm *= static_cast<KK_FLOAT>(2.0)*k3;
                    KK_FLOAT term = expterm*dcni;
                    d_num_i += c6_ref*term;
                    d_den_i += term;
                    term = expterm*dcnj;
                    d_num_j += c6_ref*term;
                    d_den_j += term;
                  }
                }
              }
              if (den > static_cast<KK_FLOAT>(1.0e-99)) {
                const KK_FLOAT deninv = static_cast<KK_FLOAT>(1.0)/den;
                C6 = num*deninv;
                dC6i = (d_num_i*den - d_den_i*num)*deninv*deninv;
                dC6j = (d_num_j*den - d_den_j*num)*deninv*deninv;
              } else {
                C6 = c6mem; dC6i = 0.0; dC6j = 0.0;
              }
            }

            const KK_FLOAT C8 = static_cast<KK_FLOAT>(3.0)*C6*r2r4i*l_r2r4(jtype)*autoang2;
            const KK_FLOAT alpha6 = l_alpha;
            const KK_FLOAT alpha8 = l_alpha + static_cast<KK_FLOAT>(2.0);

            KK_FLOAT e6 = 0.0, e8 = 0.0, fpair = 0.0;

            if (l_damping == 1) {                          // original / zero
              const KK_FLOAT r0 = r/l_r0ab(itype,jtype);
              const KK_FLOAT t6 = pow(l_rs6/r0, alpha6);
              const KK_FLOAT damp6 = static_cast<KK_FLOAT>(1.0)/(static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(6.0)*t6);
              const KK_FLOAT t8 = pow(l_rs8/r0, alpha8);
              const KK_FLOAT damp8 = static_cast<KK_FLOAT>(1.0)/(static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(6.0)*t8);
              e6 = C6*damp6*r6inv;
              e8 = C8*damp8*r8inv;
              const KK_FLOAT tmp6 = static_cast<KK_FLOAT>(6.0)*l_s6*C6*r8inv*damp6;
              const KK_FLOAT tmp8 = static_cast<KK_FLOAT>(8.0)*l_s8*C8*r10inv*damp8;
              const KK_FLOAT fp1 = -tmp6 - tmp8;
              const KK_FLOAT fp2 = tmp6*alpha6*t6*damp6 + static_cast<KK_FLOAT>(0.75)*tmp8*alpha8*t8*damp8;
              fpair = (fp1 + fp2)*factor_lj;
            } else if (l_damping == 2) {                   // zerom
              const KK_FLOAT r0 = l_r0ab(itype,jtype);
              const KK_FLOAT t6 = pow((r/(l_rs6*r0)) + l_rs8*r0, -alpha6);
              const KK_FLOAT damp6 = static_cast<KK_FLOAT>(1.0)/(static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(6.0)*t6);
              const KK_FLOAT t8 = pow((r/r0) + l_rs8*r0, -alpha8);
              const KK_FLOAT damp8 = static_cast<KK_FLOAT>(1.0)/(static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(6.0)*t8);
              e6 = C6*damp6*r6inv;
              e8 = C8*damp8*r8inv;
              const KK_FLOAT tmp6 = static_cast<KK_FLOAT>(6.0)*l_s6*C6*r8inv*damp6;
              const KK_FLOAT tmp8 = static_cast<KK_FLOAT>(8.0)*l_s8*C8*r10inv*damp8;
              const KK_FLOAT fp1 = -tmp6 - tmp8;
              const KK_FLOAT fp26 = tmp6*alpha6*t6*damp6*r/(r + l_rs6*l_rs8*r0*r0);
              const KK_FLOAT fp28 = tmp8*alpha8*t8*damp8*r/(r + l_rs8*r0*r0);
              const KK_FLOAT fp2 = fp26 + static_cast<KK_FLOAT>(0.75)*fp28;
              fpair = (fp1 + fp2)*factor_lj;
            } else {                                       // bj (3) and bjm (4)
              const KK_FLOAT r0 = sqrt(C8/C6);
              const KK_FLOAT r4 = rsq*rsq;
              const KK_FLOAT r6 = r4*rsq;
              const KK_FLOAT r8 = r6*rsq;
              const KK_FLOAT b = l_a1*r0 + l_a2;
              const KK_FLOAT b2 = b*b;
              const KK_FLOAT b6 = b2*b2*b2;
              const KK_FLOAT b8 = b6*b2;
              const KK_FLOAT t6 = r6 + b6;
              const KK_FLOAT t8 = r8 + b8;
              e6 = C6/t6;
              e8 = C8/t8;
              const KK_FLOAT tmp6 = static_cast<KK_FLOAT>(6.0)*l_s6*C6*r4/(t6*t6);
              const KK_FLOAT tmp8 = static_cast<KK_FLOAT>(8.0)*l_s8*C8*r6/(t8*t8);
              fpair = -(tmp6 + tmp8)*factor_lj;
            }

            const KK_FLOAT ed = l_s6*e6 + l_s8*e8;
            const KK_FLOAT evdwl = -ed*factor_lj;
            const KK_FLOAT rest = ed/C6;

            // dispersion pair force

            fxtmp += static_cast<KK_ACC_FLOAT>(delx*fpair);
            fytmp += static_cast<KK_ACC_FLOAT>(dely*fpair);
            fztmp += static_cast<KK_ACC_FLOAT>(delz*fpair);

            // accumulate dE/dCN (dc6): the i-side is always added; the
            // j-side only for half-list pairs that this rank owns

            dc6itmp += static_cast<KK_ACC_FLOAT>(rest*dC6i);

            if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < l_nlocal)) {
              a_f(j,0) -= static_cast<KK_ACC_FLOAT>(delx*fpair);
              a_f(j,1) -= static_cast<KK_ACC_FLOAT>(dely*fpair);
              a_f(j,2) -= static_cast<KK_ACC_FLOAT>(delz*fpair);
              a_dc6(j) += static_cast<KK_FLOAT>(rest*dC6j);
            }

            if (EVFLAG) {
              const int fullpair = ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) &&
                                    (NEWTON_PAIR || j < l_nlocal));
              if (l_eflag) {
                const KK_FLOAT ewt = fullpair ? static_cast<KK_FLOAT>(1.0) : static_cast<KK_FLOAT>(0.5);
                if (l_eflag_global) ev.evdwl += static_cast<KK_ACC_FLOAT>(ewt*evdwl);
                if (l_eflag_atom) {
                  const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*evdwl);
                  if (NEIGHFLAG != FULL) {
                    if (NEWTON_PAIR || i < l_nlocal) a_eatom(i) += epairhalf;
                    if (NEWTON_PAIR || j < l_nlocal) a_eatom(j) += epairhalf;
                  } else a_eatom(i) += epairhalf;
                }
              }
              if (l_vflag_either) {
                const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delx*fpair);
                const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*dely*fpair);
                const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delz*delz*fpair);
                const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*dely*fpair);
                const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delz*fpair);
                const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*delz*fpair);
                if (l_vflag_global) {
                  if (NEIGHFLAG != FULL) {
                    if (NEWTON_PAIR || i < l_nlocal) { ev.v[0]+=v0; ev.v[1]+=v1; ev.v[2]+=v2; ev.v[3]+=v3; ev.v[4]+=v4; ev.v[5]+=v5; }
                    if (NEWTON_PAIR || j < l_nlocal) { ev.v[0]+=v0; ev.v[1]+=v1; ev.v[2]+=v2; ev.v[3]+=v3; ev.v[4]+=v4; ev.v[5]+=v5; }
                  } else { ev.v[0]+=v0; ev.v[1]+=v1; ev.v[2]+=v2; ev.v[3]+=v3; ev.v[4]+=v4; ev.v[5]+=v5; }
                }
                if (l_vflag_atom) {
                  if (NEIGHFLAG != FULL) {
                    if (NEWTON_PAIR || i < l_nlocal) { a_vatom(i,0)+=v0; a_vatom(i,1)+=v1; a_vatom(i,2)+=v2; a_vatom(i,3)+=v3; a_vatom(i,4)+=v4; a_vatom(i,5)+=v5; }
                    if (NEWTON_PAIR || j < l_nlocal) { a_vatom(j,0)+=v0; a_vatom(j,1)+=v1; a_vatom(j,2)+=v2; a_vatom(j,3)+=v3; a_vatom(j,4)+=v4; a_vatom(j,5)+=v5; }
                  } else { a_vatom(i,0)+=v0; a_vatom(i,1)+=v1; a_vatom(i,2)+=v2; a_vatom(i,3)+=v3; a_vatom(i,4)+=v4; a_vatom(i,5)+=v5; }
                }
              }
            }
          }
        }

        a_f(i,0) += fxtmp;
        a_f(i,1) += fytmp;
        a_f(i,2) += fztmp;
        a_dc6(i) += static_cast<KK_FLOAT>(dc6itmp);
      }, ev_out);
  };

  EV_FLOAT ev;
  {
    EV_FLOAT evB;
    if (neighflag == HALF) {
      if (newton_pair) { if (evflag) launch_force.template operator()<HALF,1,1>(evB); else launch_force.template operator()<HALF,1,0>(evB); }
      else             { if (evflag) launch_force.template operator()<HALF,0,1>(evB); else launch_force.template operator()<HALF,0,0>(evB); }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) { if (evflag) launch_force.template operator()<HALFTHREAD,1,1>(evB); else launch_force.template operator()<HALFTHREAD,1,0>(evB); }
      else             { if (evflag) launch_force.template operator()<HALFTHREAD,0,1>(evB); else launch_force.template operator()<HALFTHREAD,0,0>(evB); }
    } else {
      if (newton_pair) { if (evflag) launch_force.template operator()<FULL,1,1>(evB); else launch_force.template operator()<FULL,1,0>(evB); }
      else             { if (evflag) launch_force.template operator()<FULL,0,1>(evB); else launch_force.template operator()<FULL,0,0>(evB); }
    }
    ev += evB;
  }

  if (need_dup) Kokkos::Experimental::contribute(d_dc6, dup_dc6);

  // communicate dc6 (reverse to sum, forward to distribute to ghosts)

  communicationStage = 2;
  if (newton_pair && (neighflag == HALF || neighflag == HALFTHREAD)) {
    k_dc6.template modify<DeviceType>();
    comm->reverse_comm(this);
    k_dc6.template sync<DeviceType>();
  }
  k_dc6.template modify<DeviceType>();
  comm->forward_comm(this);
  k_dc6.template sync<DeviceType>();

  // ===================================================================
  // Pass C : coordination-number-derivative contribution to the forces
  // ===================================================================

  auto launch_cnforce = [&]<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>(EV_FLOAT &ev_out) {
    auto l_x = x; auto l_type = type;
    auto l_ilist = d_ilist; auto l_numneigh = d_numneigh; auto l_neighbors = d_neighbors;
    auto l_dc6 = d_dc6; auto l_rcov = d_rcov;

    auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
    auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
    auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
    auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    const KK_FLOAT ls0 = special_lj[0], ls1 = special_lj[1], ls2 = special_lj[2], ls3 = special_lj[3];

    Kokkos::parallel_reduce("PairDispersionD3:cnforce",
      Kokkos::RangePolicy<DeviceType>(0,inum),
      KOKKOS_LAMBDA(const int ii, EV_FLOAT &ev) {
        const int i = l_ilist(ii);
        const KK_FLOAT xtmp = l_x(i,0);
        const KK_FLOAT ytmp = l_x(i,1);
        const KK_FLOAT ztmp = l_x(i,2);
        const int itype = l_type(i);
        const KK_FLOAT dc6i = l_dc6(i);
        const int jnum = l_numneigh(i);

        KK_ACC_FLOAT fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

        for (int jj = 0; jj < jnum; jj++) {
          int jraw = l_neighbors(i,jj);
          const int sb = (jraw >> SBBITS) & 3;
          const KK_FLOAT factor_lj = (sb==0) ? ls0 : ((sb==1) ? ls1 : ((sb==2) ? ls2 : ls3));
          int j = jraw & NEIGHMASK;
          const KK_FLOAT delx = xtmp - l_x(j,0);
          const KK_FLOAT dely = ytmp - l_x(j,1);
          const KK_FLOAT delz = ztmp - l_x(j,2);
          const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < l_rthr) {
            const int jtype = l_type(j);
            const KK_FLOAT r = sqrt(rsq);

            // dcn = dCN_i/dr = dCN_j/dr (only within the CN cutoff)

            KK_FLOAT dcn = 0.0;
            if (rsq < l_cn_thr) {
              const KK_FLOAT rcovij = (l_rcov(itype) + l_rcov(jtype))*autoang;
              const KK_FLOAT expterm = exp(-k1*(rcovij/r - static_cast<KK_FLOAT>(1.0)));
              const KK_FLOAT denom = (expterm + static_cast<KK_FLOAT>(1.0));
              dcn = -k1*rcovij*expterm/(rsq*denom*denom);
            }

            const KK_FLOAT fpair = dcn*(dc6i + l_dc6(j))/r*factor_lj;

            fxtmp += static_cast<KK_ACC_FLOAT>(delx*fpair);
            fytmp += static_cast<KK_ACC_FLOAT>(dely*fpair);
            fztmp += static_cast<KK_ACC_FLOAT>(delz*fpair);

            if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < l_nlocal)) {
              a_f(j,0) -= static_cast<KK_ACC_FLOAT>(delx*fpair);
              a_f(j,1) -= static_cast<KK_ACC_FLOAT>(dely*fpair);
              a_f(j,2) -= static_cast<KK_ACC_FLOAT>(delz*fpair);
            }

            if (EVFLAG && l_vflag_either) {    // energy contribution here is zero
              const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delx*fpair);
              const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*dely*fpair);
              const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delz*delz*fpair);
              const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*dely*fpair);
              const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delz*fpair);
              const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*delz*fpair);
              if (l_vflag_global) {
                if (NEIGHFLAG != FULL) {
                  if (NEWTON_PAIR || i < l_nlocal) { ev.v[0]+=v0; ev.v[1]+=v1; ev.v[2]+=v2; ev.v[3]+=v3; ev.v[4]+=v4; ev.v[5]+=v5; }
                  if (NEWTON_PAIR || j < l_nlocal) { ev.v[0]+=v0; ev.v[1]+=v1; ev.v[2]+=v2; ev.v[3]+=v3; ev.v[4]+=v4; ev.v[5]+=v5; }
                } else { ev.v[0]+=v0; ev.v[1]+=v1; ev.v[2]+=v2; ev.v[3]+=v3; ev.v[4]+=v4; ev.v[5]+=v5; }
              }
              if (l_vflag_atom) {
                if (NEIGHFLAG != FULL) {
                  if (NEWTON_PAIR || i < l_nlocal) { a_vatom(i,0)+=v0; a_vatom(i,1)+=v1; a_vatom(i,2)+=v2; a_vatom(i,3)+=v3; a_vatom(i,4)+=v4; a_vatom(i,5)+=v5; }
                  if (NEWTON_PAIR || j < l_nlocal) { a_vatom(j,0)+=v0; a_vatom(j,1)+=v1; a_vatom(j,2)+=v2; a_vatom(j,3)+=v3; a_vatom(j,4)+=v4; a_vatom(j,5)+=v5; }
                } else { a_vatom(i,0)+=v0; a_vatom(i,1)+=v1; a_vatom(i,2)+=v2; a_vatom(i,3)+=v3; a_vatom(i,4)+=v4; a_vatom(i,5)+=v5; }
              }
            }
          }
        }

        a_f(i,0) += fxtmp;
        a_f(i,1) += fytmp;
        a_f(i,2) += fztmp;
      }, ev_out);
  };

  {
    EV_FLOAT evC;
    if (neighflag == HALF) {
      if (newton_pair) { if (evflag) launch_cnforce.template operator()<HALF,1,1>(evC); else launch_cnforce.template operator()<HALF,1,0>(evC); }
      else             { if (evflag) launch_cnforce.template operator()<HALF,0,1>(evC); else launch_cnforce.template operator()<HALF,0,0>(evC); }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) { if (evflag) launch_cnforce.template operator()<HALFTHREAD,1,1>(evC); else launch_cnforce.template operator()<HALFTHREAD,1,0>(evC); }
      else             { if (evflag) launch_cnforce.template operator()<HALFTHREAD,0,1>(evC); else launch_cnforce.template operator()<HALFTHREAD,0,0>(evC); }
    } else {
      if (newton_pair) { if (evflag) launch_cnforce.template operator()<FULL,1,1>(evC); else launch_cnforce.template operator()<FULL,1,0>(evC); }
      else             { if (evflag) launch_cnforce.template operator()<FULL,0,1>(evC); else launch_cnforce.template operator()<FULL,0,0>(evC); }
    }
    ev += evC;
  }

  if (need_dup) Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += static_cast<double>(ev.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    if (need_dup) Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }
  if (vflag_atom) {
    if (need_dup) Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  copymode = 0;

  // free duplicated memory

  if (need_dup) {
    dup_cn    = {};
    dup_dc6   = {};
    dup_f     = {};
    dup_eatom = {};
    dup_vatom = {};
  }
}

/* ----------------------------------------------------------------------
   Kokkos-resident communication
     communicationStage == 1 : coordination number cn
     communicationStage == 2 : derivative dc6
------------------------------------------------------------------------- */

template<class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                                 DAT::tdual_double_1d &k_buf,
                                                                 int /*pbc_flag*/, int * /*pbc*/)
{
  auto l_sendlist = k_sendlist.template view<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();
  if (communicationStage == 1) {
    auto l_cn = d_cn;
    Kokkos::parallel_for("PairDispersionD3:pack_forward_cn",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { l_buf[i] = static_cast<double>(l_cn(l_sendlist(i))); });
  } else {
    auto l_dc6 = d_dc6;
    Kokkos::parallel_for("PairDispersionD3:pack_forward_dc6",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { l_buf[i] = static_cast<double>(l_dc6(l_sendlist(i))); });
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first,
                                                                    DAT::tdual_double_1d &k_buf)
{
  auto l_buf = k_buf.template view<DeviceType>();
  const int l_first = first;
  if (communicationStage == 1) {
    auto l_cn = d_cn;
    Kokkos::parallel_for("PairDispersionD3:unpack_forward_cn",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { l_cn(i+l_first) = static_cast<KK_FLOAT>(l_buf[i]); });
  } else {
    auto l_dc6 = d_dc6;
    Kokkos::parallel_for("PairDispersionD3:unpack_forward_dc6",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { l_dc6(i+l_first) = static_cast<KK_FLOAT>(l_buf[i]); });
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_reverse_comm_kokkos(int n, int first,
                                                                 DAT::tdual_double_1d &k_buf)
{
  auto l_buf = k_buf.template view<DeviceType>();
  const int l_first = first;
  if (communicationStage == 1) {
    auto l_cn = d_cn;
    Kokkos::parallel_for("PairDispersionD3:pack_reverse_cn",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { l_buf[i] = static_cast<double>(l_cn(l_first+i)); });
  } else {
    auto l_dc6 = d_dc6;
    Kokkos::parallel_for("PairDispersionD3:pack_reverse_dc6",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { l_buf[i] = static_cast<double>(l_dc6(l_first+i)); });
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_reverse_comm_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                                    DAT::tdual_double_1d &k_buf)
{
  auto l_sendlist = k_sendlist.template view<DeviceType>();
  auto l_buf = k_buf.template view<DeviceType>();
  if (communicationStage == 1) {
    auto l_cn = d_cn;
    Kokkos::parallel_for("PairDispersionD3:unpack_reverse_cn",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { Kokkos::atomic_add(&l_cn(l_sendlist(i)), static_cast<KK_FLOAT>(l_buf[i])); });
  } else {
    auto l_dc6 = d_dc6;
    Kokkos::parallel_for("PairDispersionD3:unpack_reverse_dc6",
      Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i) { Kokkos::atomic_add(&l_dc6(l_sendlist(i)), static_cast<KK_FLOAT>(l_buf[i])); });
  }
}

/* ----------------------------------------------------------------------
   host fallback communication (used by the /kk/host variant and legacy comm)
------------------------------------------------------------------------- */

template<class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf,
                                                          int /*pbc_flag*/, int * /*pbc*/)
{
  if (communicationStage == 1) {
    k_cn.sync_host();
    for (int i = 0; i < n; i++) buf[i] = static_cast<double>(h_cn[list[i]]);
  } else {
    k_dc6.sync_host();
    for (int i = 0; i < n; i++) buf[i] = static_cast<double>(h_dc6[list[i]]);
  }
  return n;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  if (communicationStage == 1) {
    k_cn.sync_host();
    for (int i = 0; i < n; i++) h_cn[i+first] = static_cast<KK_FLOAT>(buf[i]);
    k_cn.modify_host();
  } else {
    k_dc6.sync_host();
    for (int i = 0; i < n; i++) h_dc6[i+first] = static_cast<KK_FLOAT>(buf[i]);
    k_dc6.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairDispersionD3Kokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  if (communicationStage == 1) {
    k_cn.sync_host();
    for (int i = first; i < last; i++) buf[m++] = static_cast<double>(h_cn[i]);
  } else {
    k_dc6.sync_host();
    for (int i = first; i < last; i++) buf[m++] = static_cast<double>(h_dc6[i]);
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairDispersionD3Kokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  if (communicationStage == 1) {
    k_cn.sync_host();
    for (int i = 0; i < n; i++) h_cn[list[i]] += static_cast<KK_FLOAT>(buf[m++]);
    k_cn.modify_host();
  } else {
    k_dc6.sync_host();
    for (int i = 0; i < n; i++) h_dc6[list[i]] += static_cast<KK_FLOAT>(buf[m++]);
    k_dc6.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairDispersionD3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairDispersionD3Kokkos<LMPHostType>;
#endif
}
