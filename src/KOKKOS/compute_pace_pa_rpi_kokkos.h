/* -*- c++ -*- ----------------------------------------------------------
   compute pace/kk  —  Kokkos OpenMP descriptor compute for FitSNAP.

   Inherits ComputePACE only (NOT PairPACEKokkos) to avoid diamond via
   Pointers and to stay in modify->compute rather than force->pair.

   dgradflag is NOT supported — error at construction.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pace/kk,        ComputePACEKokkos<LMPDeviceType>);
ComputeStyle(pace/kk/device, ComputePACEKokkos<LMPDeviceType>);
ComputeStyle(pace/kk/host,   ComputePACEKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_PACE_KOKKOS_H
#define LMP_COMPUTE_PACE_KOKKOS_H

#include "compute_pace.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"

// Forward-declare so SplineInterpolatorKokkos::operator= can be declared
// without pulling ace_radial.h into every TU that includes this header.
class SplineInterpolator;

namespace LAMMPS_NS {

template<class DeviceType>
class ComputePACEKokkos : public ComputePACE {

 public:
  // Tag structs
  struct TagComputePACEComputeNeigh{};
  struct TagComputePACEComputeRadial{};
  struct TagComputePACEComputeAi{};
  struct TagComputePACEConjugateAi{};
  struct TagComputePACEComputeB{};
  struct TagComputePACEDescGrad{};

  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  using KKDeviceType = typename KKDevice<DeviceType>::value;
  using complex = SNAComplex<KK_FLOAT>;

  // Class-scope constants to shadow the global Y00=1 in ace_spherical_cart.h
  // (same trick used in pair_pace_kokkos.h)
  static constexpr KK_FLOAT Y00   = 0.2820947917738782;
  static constexpr KK_FLOAT sq3   = 1.7320508075688772;
  static constexpr KK_FLOAT sq3o2 = 1.2247448713915890;

  ComputePACEKokkos(class LAMMPS *, int, char **);
  ~ComputePACEKokkos() override;

  void init() override;
  void compute_array() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePACEComputeNeigh,
                  const typename Kokkos::TeamPolicy<DeviceType,
                  TagComputePACEComputeNeigh>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePACEComputeRadial,
                  const typename Kokkos::TeamPolicy<DeviceType,
                  TagComputePACEComputeRadial>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePACEComputeAi,
                  const typename Kokkos::TeamPolicy<DeviceType,
                  TagComputePACEComputeAi>::member_type& team) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePACEConjugateAi, const int& ii) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePACEComputeB, const int& iter) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePACEDescGrad, const int& iter) const;

 protected:
  int inum, maxneigh, chunk_size, chunk_offset;
  int idx_ms_combs_max, idx_sph_max;
  int nelements, lmax, nradmax, nradbase;
  int nfuncs_total;

  class AtomKokkos *atomKK;
  int execution_space;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread type;

  typedef Kokkos::DualView<KK_FLOAT**, DeviceType> tdual_fparams;
  tdual_fparams k_cutsq;
  typedef Kokkos::View<KK_FLOAT**, DeviceType> t_fparams;
  t_fparams d_cutsq;

  // View type aliases
  typedef Kokkos::View<int*,                                   DeviceType> t_ace_1i;
  typedef Kokkos::View<int**,                                  DeviceType> t_ace_2i;
  typedef Kokkos::View<int**, Kokkos::LayoutRight,             DeviceType> t_ace_2i_lr;
  typedef Kokkos::View<int***,                                 DeviceType> t_ace_3i;
  typedef Kokkos::View<int***, Kokkos::LayoutRight,            DeviceType> t_ace_3i_lr;
  typedef Kokkos::View<KK_FLOAT*,                              DeviceType> t_ace_1d;
  typedef Kokkos::View<KK_FLOAT**,                             DeviceType> t_ace_2d;
  typedef Kokkos::View<KK_FLOAT**, Kokkos::LayoutRight,        DeviceType> t_ace_2d_lr;
  typedef Kokkos::View<KK_FLOAT*[3],                           DeviceType> t_ace_2d3;
  typedef Kokkos::View<KK_FLOAT***,                            DeviceType> t_ace_3d;
  typedef Kokkos::View<KK_FLOAT**[3],                          DeviceType> t_ace_3d3;
  typedef Kokkos::View<KK_FLOAT****,                           DeviceType> t_ace_4d;
  // compile-time last dim = 4 (used for spline lookup table)
  typedef Kokkos::View<KK_FLOAT**[4], Kokkos::LayoutRight,     DeviceType> t_ace_3d4_lr;
  typedef Kokkos::View<KK_FLOAT****, Kokkos::LayoutRight,      DeviceType> t_ace_4d_lr;
  typedef Kokkos::View<complex*,                               DeviceType> t_ace_1c;
  typedef Kokkos::View<complex**,                              DeviceType> t_ace_2c;
  typedef Kokkos::View<complex***,                             DeviceType> t_ace_3c;
  typedef Kokkos::View<complex****,                            DeviceType> t_ace_4c;
  typedef Kokkos::View<complex***[3],                          DeviceType> t_ace_3c3;

  // A-matrix and product Views
  t_ace_3d  A_rank1;           // (natom, nelements, nradbase)
  t_ace_4c  A;                 // (natom, nelements, (lmax+1)^2, nradmax)
  t_ace_4c  A_sph;             // (natom, nelements, idx_sph_max, nradmax)
  t_ace_3c  A_list;            // (natom, idx_ms_combs_max, rankmax)
  t_ace_3c  A_forward_prod;    // (natom, idx_ms_combs_max, rankmax+1)
  t_ace_3c  dB_flatten;        // (natom, idx_ms_combs_max, rankmax)

  // Radial Views
  t_ace_4d  fr, dfr;           // (natom, maxneigh, lmax+1, nradmax)
  t_ace_3d  gr, dgr;           // (natom, maxneigh, nradbase)
  t_ace_3d  d_values, d_derivatives;

  // Neighbour geometry
  t_ace_1i  d_ncount;
  t_ace_2d  d_mu;
  t_ace_3d3 d_rhats;
  t_ace_2d  d_rnorms;
  t_ace_2i  d_nearest;

  // Hard-core repulsion (2D: natom × maxneigh)
  t_ace_2d  cr, dcr;
  t_ace_1d  rho_core;

  // Spherical harmonic helpers
  t_ace_1d  d_idx_sph;
  t_ace_1d  alm, blm, cl, dl;

  // Per-element basis metadata
  t_ace_1i   d_idx_ms_combs_count;
  t_ace_2i_lr d_rank, d_num_ms_combs, d_idx_funcs;
  t_ace_3i_lr d_mus, d_ns, d_ls, d_ms_combs;
  t_ace_3d    d_ctildes;       // (nelements, idx_ms_combs_max, ndensitymax)
  t_ace_1i    d_type_offsets;  // column offset per element (0-indexed)
  t_ace_1i    d_ndensity;      // ndensity per element

  // Descriptor output
  t_ace_2d  d_descriptors;     // (chunk_size, nfuncs_total)
  t_ace_4d  d_desc_force_ij;   // (chunk_size, maxneigh, nfuncs_total, 3)

  // Spline interpolator (identical layout to PairPACEKokkos)
  struct SplineInterpolatorKokkos {
    int ntot, nlut, num_of_functions;
    KK_FLOAT cutoff, deltaSplineBins, invrscalelookup, rscalelookup;
    t_ace_3d4_lr lookupTable;  // (ntot+1, num_of_functions, [4])

    void operator=(const SplineInterpolator &spline);
    void deallocate() { lookupTable = t_ace_3d4_lr(); }
    KK_FLOAT memory_usage() {
      return lookupTable.span() *
             sizeof(typename decltype(lookupTable)::value_type);
    }
    KOKKOS_INLINE_FUNCTION
    void calcSplines(int ii, int jj, KK_FLOAT r,
                     const t_ace_3d &vals, const t_ace_3d &derivs) const;
  };

  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_gk;
  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_rnl;
  Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType> k_splines_hc;

  void grow(int natom, int maxneigh);
  void copy_pertype();
  void copy_splines();
  void copy_tilde();
  void pre_compute_harmonics(int lmax);
  void deallocate_views_of_views();

  KOKKOS_INLINE_FUNCTION
  void evaluate_splines(int ii, int jj, KK_FLOAT r,
                        int nradbase, int nradmax, int mu_i, int mu_j) const;

  KOKKOS_INLINE_FUNCTION
  void cutoff_func_poly(KK_FLOAT r, KK_FLOAT r_in, KK_FLOAT delta,
                        KK_FLOAT &fc, KK_FLOAT &dfc) const;

  template<class TagStyle>
  void check_team_size_for(int inum, int &team_size, int vector_length);

  template<typename scratch_type>
  int scratch_size_helper(int values_per_team);
};

}  // namespace LAMMPS_NS
#endif
#endif
