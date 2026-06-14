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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(dispersion/d3/kk,PairDispersionD3Kokkos<LMPDeviceType>);
PairStyle(dispersion/d3/kk/device,PairDispersionD3Kokkos<LMPDeviceType>);
PairStyle(dispersion/d3/kk/host,PairDispersionD3Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_DISPERSION_D3_KOKKOS_H
#define LMP_PAIR_DISPERSION_D3_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_dispersion_d3.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairDispersionD3Kokkos : public PairDispersionD3, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairDispersionD3Kokkos(class LAMMPS *);
  ~PairDispersionD3Kokkos() override;

  void compute(int, int) override;
  void init_style() override;

  // device (Kokkos-resident) communication
  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &, int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  int pack_reverse_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &) override;

  // host fallback communication (used for the /kk/host variant and legacy comm)
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  void allocate_tables_kokkos();

  // atom data

  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d type;

  // per-atom energy / virial

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  // per-atom coordination number (cn) and dE_disp/dCN (dc6)
  //   cn  : analog of the EAM electron density rho
  //   dc6 : analog of the EAM embedding-energy derivative fp

  DAT::tdual_kkfloat_1d k_cn, k_dc6;
  typename AT::t_kkfloat_1d d_cn, d_dc6;
  HAT::t_kkfloat_1d h_cn, h_dc6;

  // per-type and per-type-pair coefficient tables on device

  DAT::tdual_kkfloat_1d k_r2r4, k_rcov;
  DAT::tdual_int_1d k_mxci;
  DAT::tdual_kkfloat_2d k_r0ab;
  typename AT::t_kkfloat_1d d_r2r4, d_rcov;
  typename AT::t_int_1d d_mxci;
  typename AT::t_kkfloat_2d d_r0ab;

  // c6ab(itype,jtype,ci,cj,k): k=0 -> C6 reference, k=1/2 -> reference CNs
  typedef Kokkos::DualView<KK_FLOAT*****, Kokkos::LayoutRight, DeviceType> tdual_c6_5d;
  tdual_c6_5d k_c6ab;
  typename tdual_c6_5d::t_dev d_c6ab;

  // neighbor list

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;

  // scatter views (duplicated for threaded host, atomic for GPU, plain for serial)

  int need_dup;
  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> dup_cn;
  DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> dup_dc6;
  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> ndup_cn;
  NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> ndup_dc6;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  int neighflag, newton_pair;
  int nlocal, nall, eflag, vflag, inum;
  KK_FLOAT special_lj[4];

  friend void pair_virial_fdotr_compute<PairDispersionD3Kokkos>(PairDispersionD3Kokkos*);
};

}    // namespace LAMMPS_NS
#endif
#endif
