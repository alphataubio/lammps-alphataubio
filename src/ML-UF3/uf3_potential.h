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

#ifndef UF3_POTENTIAL_H
#define UF3_POTENTIAL_H

#include "pointers.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace LAMMPS_NS {

class UF3Potential : protected Pointers {
 public:
  UF3Potential(class LAMMPS *, const std::string &, double **cutsq_, int **, std::vector<std::string> &, bool);
  ~UF3Potential() override;

  double memory_usage();

  int get_starting_index_2b(int i, int j, double r);
  int get_starting_index_3b(int i, int j, int k, double r, int knot_dim);

  int allocated;
  void allocate();

  int **n2b_knots_array_size, **n2b_coeff_array_size;
  double **cut_2b, ***n2b_knots_array, ***n2b_coeff_array;

  int ***map_3b, **n3b_knots_array_size, **n3b_coeff_array_size;
  double ***cut_3b, **cut_3b_list, ****min_cut_3b, ***n3b_knots_array, ****n3b_coeff_array;
  double ****cached_constants_2b, ****cached_constants_2b_deri;
  double ****cached_constants_3b, ****cached_constants_3b_deri;
  double ****coeff_for_der_jk, ****coeff_for_der_ik, ****coeff_for_der_ij;

  int max_num_knots_2b, max_num_coeff_2b, max_num_knots_3b, max_num_coeff_3b, tot_interaction_count_3b;

 protected:

  double **cutsq;     // cutoff sq for each atom pair
  int **setflag;      // 0/1 = whether each i,j has been set

  std::vector<std::string> elements;
  std::unordered_map<std::string, int> element_to_type;

  int ***setflag_3b, **knot_spacing_type_2b, ***knot_spacing_type_3b;
  double **knot_spacing_2b, ****knot_spacing_3b;

  void uf3_read_unified_pot_file(const std::string &potf_name);
  void communicate();
  int bsplines_created;
  bool pot_3b;

  void create_bsplines();
  void create_cached_constants_2b();
  void create_cached_constants_3b();

  int (UF3Potential::*get_starting_index_2b_ptr)(int i, int j, double r);
  int (UF3Potential::*get_starting_index_3b_ptr)(int i, int j, int k, double r, int knot_dim);
  int get_starting_index_uniform_2b(int i, int j, double r);
  int get_starting_index_uniform_3b(int i, int j, int k, double r, int knot_dim);
  int get_starting_index_nonuniform_2b(int i, int j, double r);
  int get_starting_index_nonuniform_3b(int i, int j, int k, double r, int knot_dim);

};

}    // namespace LAMMPS_NS

#endif // !UF3_POTENTIAL_H

