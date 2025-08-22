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
PairStyle(libxc/gga,PairLibXCGGA);
// clang-format on
#else

#ifndef LMP_PAIR_LIBXC_GGA_H
#define LMP_PAIR_LIBXC_GGA_H

#include "pair_libxc.h"
#include <xc.h>
#include <vector>

namespace LAMMPS_NS {

class PairLibXCGGA : public PairLibXC {
 public:
  PairLibXCGGA(class LAMMPS *);
  ~PairLibXCGGA() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  
 protected:
  // GGA-specific members for gradient computation
  double ***density_gradient;    // Gradient of electron density
  double **laplacian_density;    // Laplacian of electron density
  
  // Additional parameters for GGA
  double gradient_cutoff;         // Cutoff for gradient contributions
  double enhancement_factor;      // Enhancement factor for GGA corrections
  
  // GGA-specific methods
  void compute_density_gradient(int, int, double *, double **, double *);
  void compute_gga_xc_energy_force(int, int, double *, double **, double *, double *, double **);
  void compute_gradient_vector(int, double *, double *, double *);
  double compute_gradient_magnitude(double *);
  
  // Meta-GGA extensions
  bool use_meta_gga;              // Flag for meta-GGA functionals
  double **kinetic_energy_density; // Kinetic energy density tau
  void compute_kinetic_energy_density(int, int, double *);
  
  // Utility functions for GGA
  void allocate_gga();
  void deallocate_gga();
  double gradient_correction_factor(double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
