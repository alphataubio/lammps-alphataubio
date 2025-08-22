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
PairStyle(libxc,PairLibXC);
// clang-format on
#else

#ifndef LMP_PAIR_LIBXC_H
#define LMP_PAIR_LIBXC_H

#include "pair.h"
#include <xc.h>
#include <vector>
#include <map>

namespace LAMMPS_NS {

class PairLibXC : public Pair {
 public:
  PairLibXC(class LAMMPS *);
  ~PairLibXC() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double cut_global;
  double **cut;
  double **offset;
  
  // LibXC specific members
  xc_func_type *xc_func_x;      // Exchange functional
  xc_func_type *xc_func_c;      // Correlation functional
  xc_func_type *xc_func_xc;     // Combined XC functional
  
  int xc_functional_x;           // Exchange functional ID
  int xc_functional_c;           // Correlation functional ID  
  int xc_functional_xc;          // Combined XC functional ID
  
  bool use_combined_xc;          // Use combined XC functional instead of separate X and C
  
  // Electron density parameters
  double **rho0;                 // Reference electron density for each atom type
  double **decay_length;         // Decay length for electron density
  double **atomic_volume;        // Atomic volume for each type
  
  // Grid parameters for density computation
  int ngrid;                     // Number of grid points for numerical integration
  double grid_spacing;           // Grid spacing for numerical integration
  
  // Energy and force calculation methods
  void compute_density(int, int, double *, double *, double *);
  void compute_xc_energy_force(int, int, double *, double *, double *);
  double compute_embedding_energy(double, int);
  void compute_embedding_force(double, double *, int);
  
  // Electron density functions
  double electron_density(double, int);
  double electron_density_derivative(double, int);
  
  // Utility functions
  void allocate();
  void parse_functional_name(const char *, int &);
  double interpolate_density(double, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
