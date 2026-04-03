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
PairStyle(uf3,PairUF3);
// clang-format on
#else

#ifndef LMP_PAIR_UF3_H
#define LMP_PAIR_UF3_H

#include "pair.h"

namespace LAMMPS_NS {

class PairUF3 : public Pair {
 public:
  PairUF3(class LAMMPS *);
  ~PairUF3() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;    // needed for ptr to full neigh list
  double init_one(int, int) override;                 // needed for cutoff radius for neighbour list
  double single(int, int, int, int, double, double, double, double &) override;
  double memory_usage() override;

 protected:

  class UF3Potential *uf3_potential;

  int *neighshort, maxshort;    // short neighbor list array for 3body interaction

  bool pot_3b;

  virtual void allocate();

  int nbody_flag;
  int max_num_knots_2b;
  int max_num_coeff_2b;
  int max_num_knots_3b;
  int max_num_coeff_3b;
  int tot_interaction_count_3b;
};

}    // namespace LAMMPS_NS

#endif
#endif
