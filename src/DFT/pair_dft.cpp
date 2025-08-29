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

#include "pair_dft.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>

using namespace LAMMPS_NS;

// Include the implementations from the .hpp files
#include "pair_dft_basis.hpp"
#include "pair_dft_libxc.hpp"
#include "pair_dft_libint2.hpp"
#include "pair_dft_scf.hpp"

static constexpr double ANGSTROM_TO_BOHR = 1/0.52917721067;

/* ---------------------------------------------------------------------- */

PairDFT::PairDFT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  manybody_flag = 1;
  one_coeff = 0;
  
  // Initialize pointers
  xc_func_x = nullptr;
  xc_func_c = nullptr;
  xc_func_xc = nullptr;
  
  // SCF parameters
  energy_tolerance = 1.0e-8;
  density_tolerance = 1.0e-6;
  max_scf_iterations = 100;
  current_iteration = 0;
  scf_converged = false;
  
  // Initialize energy components
  total_dft_energy = 0.0;
  kinetic_energy = 0.0;
  nuclear_repulsion = 0.0;
  xc_energy = 0.0;
  
  // Grid parameters
  grid_size = 50000;
  grid_type = "Lebedev";
  
  // Basis function information
  n_basis_functions = 0;
  n_occupied_orbitals = 0;
  n_electrons = 0;
  
  // Initialize libint2
  libint2::initialize();
}

/* ---------------------------------------------------------------------- */

PairDFT::~PairDFT()
{
  // Cleanup LibXC
  cleanup_libxc();
  
  // Cleanup libint2
  cleanup_libint();
  libint2::finalize();
  
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  
  // Initialize basis set and integral engines if needed
  if (n_basis_functions == 0) initialize_basis_set();
  
  // Update basis set positions with current atom positions
  std::vector<std::vector<double>> positions;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    positions.push_back({
      x[i][0] * ANGSTROM_TO_BOHR,
      x[i][1] * ANGSTROM_TO_BOHR,
      x[i][2] * ANGSTROM_TO_BOHR
    });
  }
  set_atom_positions(positions);
  
  // Perform SCF calculation for the current configuration
  perform_scf();
  
  // Compute DFT forces
  compute_hellmann_feynman_forces();
  compute_pulay_forces();
  
  // Store energy
  if (eflag) eng_vdwl = total_dft_energy;
  
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairDFT::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal pair_style dft command");
  
  // Parse: pair_style dft <functional> <basis.json> [options]
  functional_name = std::string(arg[0]);
  basis_file = std::string(arg[1]);
  
  // Parse the functional name and initialize LibXC
  parse_functional_name(functional_name.c_str());
  
  // Check if basis file exists
  std::ifstream file(basis_file);
  if (!file.good()) {
    error->all(FLERR, "Cannot open basis set file");
  }
  file.close();
  
  // Parse additional options
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "grid") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing grid size");
      grid_size = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "tol") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing tolerance");
      energy_tolerance = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      density_tolerance = energy_tolerance * 100.0;
      iarg += 2;
    } else if (strcmp(arg[iarg], "maxiter") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing max iterations");
      max_scf_iterations = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown pair_style dft option");
    }
  }
  
  // Set cutoff - for DFT we typically use a large cutoff
  cut_global = 20.0;  // Angstroms
}

/* ---------------------------------------------------------------------- */

void PairDFT::coeff(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDFT::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDFT::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style dft requires atom IDs");
  if (atom->q_flag == 0) error->all(FLERR, "Pair style dft requires atom attribute q");
  
  // Request standard neighbor list
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDFT::init_one(int i, int j)
{
  return cut_global;
}
