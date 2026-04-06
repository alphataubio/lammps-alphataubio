/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(uf3,ComputeUF3);
// clang-format on
#else

#ifndef LMP_COMPUTE_UF3_H
#define LMP_COMPUTE_UF3_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeUF3 : public Compute {
 public:
  ComputeUF3(class LAMMPS *, int, char **);
  ~ComputeUF3() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:

  class NeighList *list;
  class UF3Potential *uf3_potential;
  int *neighshort, maxshort;    // short neighbor list array for 3body interaction

  double **cutsq;     // cutoff sq for each atom pair
  int **setflag;      // 0/1 = whether each i,j has been set

  bool pot_3b;
  int lastcol, virial_flag;
  double **array_local;

  Compute *c_pe, *c_virial;
  std::string id_virial;

};

}    // namespace LAMMPS_NS

#endif // !LMP_COMPUTE_UF3_H
#endif
