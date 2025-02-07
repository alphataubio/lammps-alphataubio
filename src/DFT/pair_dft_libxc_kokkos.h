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
PairStyle(dft/kk,PairDftKokkos<LMPDeviceType>);
PairStyle(dft/kk/device,PairDftKokkos<LMPDeviceType>);
PairStyle(dft/kk/host,PairDftKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_DFT_KOKKOS_H
#define LMP_PAIR_DFT_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_dft.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairDftKokkos : public PairEAM, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairDftKokkos(class LAMMPS *);
  ~PairDftKokkos() override;
  void compute(int, int) override;
  void init_style() override;


 protected:

};

}
#endif
#endif

