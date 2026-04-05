/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing author: Ajinkya Hire (Univ. of Florida),
                        Hendrik Krass (Univ. of Constance),
                        Matthias Rupp (Luxembourg Institute of Science and Technology),
                        Richard Hennig (Univ of Florida)
---------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(uf3/kk,PairUF3Kokkos<LMPDeviceType>)
PairStyle(uf3/kk/device,PairUF3Kokkos<LMPDeviceType>)
PairStyle(uf3/kk/host,PairUF3Kokkos<LMPHostType>)
// clang-format on
#else

#ifndef LMP_PAIR_UF3_KOKKOS_H
#define LMP_PAIR_UF3_KOKKOS_H

#include "kokkos.h"
#include "pair_kokkos.h"
#include "pair_uf3.h"

template <int NEIGHFLAG, int EVFLAG> struct TagPairUF3ComputeFullA {};
struct TagPairUF3ComputeShortNeigh {};

namespace LAMMPS_NS {

template <class DeviceType> class PairUF3Kokkos : public PairUF3 {
 public:
  PairUF3Kokkos(class LAMMPS *);
  ~PairUF3Kokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void allocate() override;
  void init_style() override;
  void init_list(int, class NeighList *) override;    // needed for ptr to full neigh list
  double init_one(int, int) override;                 // needed for cutoff radius for neighbour list
  double single(int, int, int, int, double, double, double, double &) override;

  template <typename T, typename V> void copy_2d(V &d, T **h, int m, int n);
  template <typename T, typename V> void copy_3d(V &d, T ***h, int m, int n, int o);

  template <int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION void operator()(TagPairUF3ComputeFullA<NEIGHFLAG, EVFLAG>, const int &,
                                         EV_FLOAT &) const;

  template <int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION void operator()(TagPairUF3ComputeFullA<NEIGHFLAG, EVFLAG>,
                                         const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairUF3ComputeShortNeigh, const int &) const;

  enum { EnabledNeighFlags = FULL };
  enum { COUL_FLAG = 0 };
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

 protected:
  DAT::ttransform_kkfloat_2d k_cutsq;//Create a DualView, defination of tdual_kkfloat_2d in kokkos_type.h
  typename AT::t_kkfloat_2d d_cutsq; //t_kkfloat_2d = t_dev ==> Creates a new View d_cutsq
  //the type of d_cutsq is decided by the Device(not host) type for the DualView k_cutsq
  //Meaning the memory location of d_cutsq is the same as the Device(not host) memory location of
  //k_cutsq
  DAT::ttransform_kkfloat_3d k_cut_3b;
  DAT::ttransform_kkfloat_4d k_min_cut_3b;
  template <typename TYPE> void destroy_3d(TYPE data, typename TYPE::value_type*** &array);
  template <typename TYPE> void destroy_4d(TYPE data, typename TYPE::value_type**** &array);
  Kokkos::View<KK_FLOAT **, LMPDeviceLayout, LMPDeviceType> /*d_cutsq,*/ d_cut_3b_list;

  typename AT::t_int_2d map2b;
  typename AT::t_int_3d map3b;
  typename AT::t_kkfloat_2d d_n2b_knot;
  typename AT::t_kkfloat_3d constants_2b, dnconstants_2b, d_cut_3b, d_n3b_knot_matrix;
  typename AT::t_kkfloat_4d d_min_cut_3b, constants_3b, dnconstants_3b, d_coefficients_3b;
  Kokkos::View<KK_FLOAT*****, LMPDeviceLayout, LMPDeviceType> d_dncoefficients_3b;

  int coefficients_created = 0;
  void create_coefficients();

  template <int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void twobody(const int itype, const int jtype, const KK_FLOAT r, KK_FLOAT &evdwl,
               KK_FLOAT &fpair) const;
  template <int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void threebody(const int itype, const int jtype, const int ktype, const KK_FLOAT value_rij,
                 const KK_FLOAT value_rik, const KK_FLOAT value_rjk, KK_FLOAT &evdwl3,
                 KK_FLOAT (&fforce)[3]) const;

  template <int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION void
  ev_tally(EV_FLOAT &ev, const int &i, const int &j, const KK_FLOAT &epair, const KK_FLOAT &fpair,
           const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template <int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION void ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
                                        const KK_FLOAT &evdwl, const KK_FLOAT &ecoul, KK_ACC_FLOAT *fj,
                                        KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_tagint_1d tag;
  typename AT::t_int_1d_randomread type;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  DAT::ttransform_kkacc_1d_9 k_cvatom;
  
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
  typename AT::t_kkacc_1d_9 d_cvatom;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  using ScatterFType = Kokkos::Experimental::ScatterView<KK_ACC_FLOAT *[3], DAT::t_kkacc_1d_3::array_layout,
                                                         KKDeviceType>;
  ScatterFType fscatter;
  using ScatterVType = Kokkos::Experimental::ScatterView<KK_ACC_FLOAT *[6], LMPDeviceLayout,
                                                         KKDeviceType>;
  ScatterVType vscatter;
  using ScatterCVType = Kokkos::Experimental::ScatterView<KK_ACC_FLOAT *[9], LMPDeviceLayout,
                                                          KKDeviceType>;
  ScatterCVType cvscatter;
  using ScatterEType = Kokkos::Experimental::ScatterView<KK_ACC_FLOAT *, Kokkos::LayoutRight,
                                                         KKDeviceType>;
  ScatterEType escatter;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int neighflag, newton_pair;
  int nlocal, nall, eflag, vflag;

  int inum;
  Kokkos::View<int **, DeviceType> d_neighbors_short;
  Kokkos::View<int *, DeviceType> d_numneigh_short;

  friend void pair_virial_fdotr_compute<PairUF3Kokkos>(PairUF3Kokkos *);
};

}    // namespace LAMMPS_NS

#endif // !LMP_PAIR_UF3_KOKKOS_H
#endif

