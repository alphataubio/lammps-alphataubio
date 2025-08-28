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

#include <libint2.hpp>



using namespace LAMMPS_NS;
using json = nlohmann_lmp::json;

/* ---------------------------------------------------------------------- */

PairDFT::PairDFT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 0;
  
  // Initialize LibXC pointers
  xc_func_x = nullptr;
  xc_func_c = nullptr;
  xc_func_xc = nullptr;
  
  xc_functional_x = -1;
  xc_functional_c = -1;
  xc_functional_xc = -1;
  
  use_combined_xc = false;
  is_hybrid = false;
  is_meta_gga = false;
  is_range_separated = false;
  
  hybrid_coeff = 0.0;
  range_separation_param = 0.0;
  
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
  electron_nuclear_energy = 0.0;
  electron_electron_energy = 0.0;
  xc_energy = 0.0;
  exact_exchange_energy = 0.0;
  
  // Grid parameters
  grid_size = 50000;  // Medium grid
  grid_tolerance = 1.0e-10;
  grid_type = "Lebedev";
  
  // Basis function information
  n_basis_functions = 0;
  n_occupied_orbitals = 0;
  n_electrons = 0;
  
  // Dispersion correction
  use_dispersion = false;
  dispersion_type = "D3BJ";
  dispersion_energy = 0.0;
  
  atomic_charges = nullptr;
  vdw_radii = nullptr;
  
  // Initialize libint2
  libint2::initialize();
}

/* ---------------------------------------------------------------------- */

PairDFT::~PairDFT()
{
  cleanup_libxc();
  cleanup_integrals();
  
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(offset);
    memory->destroy(atomic_charges);
    memory->destroy(vdw_radii);
  }
  
  // Cleanup libint2
  libint2::finalize();
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Initialize basis set and integral engines if needed
  if (n_basis_functions == 0) {
    initialize_basis_set();
  }
  
  // Perform SCF calculation for the current configuration
  perform_scf();
  
  // Compute DFT forces
  compute_forces(eflag, vflag);
  
  // Add dispersion correction if requested
  if (use_dispersion) {
    compute_dispersion_correction();
  }
  
  // Add classical repulsion at short range (optional)
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        
        // Add short-range repulsion to prevent core overlap
        if (r < 0.5) {  // Very short range
          double repulsion = 100.0 * exp(-10.0 * r);
          fpair = 1000.0 * exp(-10.0 * r) / r;
          
          f[i][0] += delx * fpair;
          f[i][1] += dely * fpair;
          f[i][2] += delz * fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx * fpair;
            f[j][1] -= dely * fpair;
            f[j][2] -= delz * fpair;
          }
          
          if (eflag) {
            evdwl = repulsion;
          }
          
          if (evflag) ev_tally(i, j, nlocal, newton_pair,
                              evdwl, 0.0, fpair, delx, dely, delz);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   Initialize basis set from JSON file
------------------------------------------------------------------------- */

void PairDFT::initialize_basis_set()
{
  basis_manager = std::make_unique<BasisSetManager>();
  basis_manager->load_from_json(basis_file);
  
  n_basis_functions = basis_manager->get_n_basis();
  
  // Initialize integral engine
  integral_engine = std::make_unique<IntegralEngine>(basis_manager.get());
  
  // Initialize density matrix
  density_matrix = std::make_unique<DensityMatrix>(n_basis_functions);
  
  // Initialize XC functional
  xc_functional = std::make_unique<XCFunctional>(functional_name);
  
  // Initialize grid integrator
  grid_integrator = std::make_unique<GridIntegrator>(grid_size, grid_type);
  
  // Allocate matrices
  overlap_matrix.resize(n_basis_functions, n_basis_functions);
  kinetic_matrix.resize(n_basis_functions, n_basis_functions);
  nuclear_matrix.resize(n_basis_functions, n_basis_functions);
  coulomb_matrix.resize(n_basis_functions, n_basis_functions);
  exchange_matrix.resize(n_basis_functions, n_basis_functions);
  fock_matrix.resize(n_basis_functions, n_basis_functions);
  density_matrix_eigen.resize(n_basis_functions, n_basis_functions);
  mo_coefficients.resize(n_basis_functions, n_basis_functions);
  mo_energies.resize(n_basis_functions);
  
  // Compute one-electron integrals
  compute_one_electron_integrals();
  
  // Initialize density guess
  initialize_density_guess();
}

/* ----------------------------------------------------------------------
   Perform SCF calculation
------------------------------------------------------------------------- */

void PairDFT::perform_scf()
{
  print_scf_header();
  
  double prev_energy = 0.0;
  scf_converged = false;
  
  for (current_iteration = 1; current_iteration <= max_scf_iterations; current_iteration++) {
    // Build Fock matrix
    build_fock_matrix();
    
    // Solve Roothaan-Hall equations
    solve_roothaan_hall();
    
    // Update density matrix
    compute_density_matrix();
    
    // Mix with previous density for better convergence
    if (current_iteration > 1) {
      density_matrix->mix_with_previous(0.5);
    }
    
    // Compute energy
    compute_energy();
    
    // Check convergence
    double energy_change = std::abs(total_dft_energy - prev_energy);
    double density_change = density_matrix->get_change();
    
    print_scf_iteration();
    
    if (energy_change < energy_tolerance && density_change < density_tolerance) {
      scf_converged = true;
      break;
    }
    
    prev_energy = total_dft_energy;
  }
  
  if (!scf_converged) {
    if (comm->me == 0) {
      error->warning(FLERR, "SCF did not converge within maximum iterations");
    }
  }
  
  print_scf_summary();
}

/* ----------------------------------------------------------------------
   Build Fock matrix
------------------------------------------------------------------------- */

void PairDFT::build_fock_matrix()
{
  // Core Hamiltonian
  fock_matrix = kinetic_matrix + nuclear_matrix;
  
  // Two-electron part: Coulomb and Exchange
  integral_engine->compute_eri_with_density(density_matrix_eigen, 
                                           coulomb_matrix, exchange_matrix);
  
  // Add Coulomb contribution
  fock_matrix += 2.0 * coulomb_matrix;
  
  // Handle exchange based on functional type
  if (xc_functional->is_hybrid()) {
    double alpha = xc_functional->get_exact_exchange_fraction();
    fock_matrix -= alpha * exchange_matrix;
    
    // Add XC contribution (1-alpha) for exchange + full correlation
    Eigen::MatrixXd vxc_matrix(n_basis_functions, n_basis_functions);
    grid_integrator->integrate_xc(xc_functional.get(), density_matrix_eigen,
                                  basis_manager.get(), xc_energy, vxc_matrix);
    fock_matrix += (1.0 - alpha) * vxc_matrix;
  } else {
    // Pure DFT functional - use grid integration for XC
    Eigen::MatrixXd vxc_matrix(n_basis_functions, n_basis_functions);
    grid_integrator->integrate_xc(xc_functional.get(), density_matrix_eigen,
                                  basis_manager.get(), xc_energy, vxc_matrix);
    fock_matrix += vxc_matrix;
  }
  
  // Special handling for range-separated functionals
  if (xc_functional->is_range_separated()) {
    compute_exact_exchange();
  }
}

/* ----------------------------------------------------------------------
   Solve Roothaan-Hall equations
------------------------------------------------------------------------- */

void PairDFT::solve_roothaan_hall()
{
  // Transform Fock matrix to orthogonal basis
  // F' = S^(-1/2) * F * S^(-1/2)
  
  // Compute S^(-1/2) using eigendecomposition
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(overlap_matrix);
  Eigen::MatrixXd S_sqrt = es.operatorSqrt();
  Eigen::MatrixXd S_sqrt_inv = es.operatorInverseSqrt();
  
  // Transform Fock matrix
  Eigen::MatrixXd F_prime = S_sqrt_inv * fock_matrix * S_sqrt_inv;
  
  // Diagonalize F'
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_fock(F_prime);
  mo_energies = es_fock.eigenvalues();
  Eigen::MatrixXd C_prime = es_fock.eigenvectors();
  
  // Back-transform eigenvectors to non-orthogonal basis
  mo_coefficients = S_sqrt_inv * C_prime;
}

/* ----------------------------------------------------------------------
   Compute density matrix from MO coefficients
------------------------------------------------------------------------- */

void PairDFT::compute_density_matrix()
{
  density_matrix->update_from_mo(mo_coefficients, n_occupied_orbitals);
  density_matrix_eigen = density_matrix->get_current();
}

/* ----------------------------------------------------------------------
   Compute total energy
------------------------------------------------------------------------- */

void PairDFT::compute_energy()
{
  // Nuclear repulsion
  nuclear_repulsion = compute_nuclear_repulsion_energy();
  
  // One-electron energy
  Eigen::MatrixXd H_core = kinetic_matrix + nuclear_matrix;
  double one_electron = (density_matrix_eigen.cwiseProduct(H_core)).sum();
  
  // Two-electron energy
  double j_energy = (density_matrix_eigen.cwiseProduct(coulomb_matrix)).sum();
  double k_energy = 0.0;
  
  if (xc_functional->is_hybrid()) {
    double alpha = xc_functional->get_exact_exchange_fraction();
    k_energy = alpha * (density_matrix_eigen.cwiseProduct(exchange_matrix)).sum();
  }
  
  // Total energy
  kinetic_energy = (density_matrix_eigen.cwiseProduct(kinetic_matrix)).sum();
  electron_nuclear_energy = (density_matrix_eigen.cwiseProduct(nuclear_matrix)).sum();
  electron_electron_energy = 2.0 * j_energy - k_energy;
  
  total_dft_energy = nuclear_repulsion + one_electron + j_energy - k_energy + xc_energy;
  
  // Add dispersion if requested
  if (use_dispersion) {
    total_dft_energy += dispersion_energy;
  }
}

/* ----------------------------------------------------------------------
   Compute nuclear repulsion energy
------------------------------------------------------------------------- */

double PairDFT::compute_nuclear_repulsion_energy()
{
  double energy = 0.0;
  
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    for (int j = i + 1; j < nlocal; j++) {
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      double r = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (r > 1e-10) {
        double Zi = atomic_charges[type[i]][type[i]];
        double Zj = atomic_charges[type[j]][type[j]];
        energy += Zi * Zj / r;
      }
    }
  }
  
  return energy;
}

/* ----------------------------------------------------------------------
   Initialize density guess
------------------------------------------------------------------------- */

void PairDFT::initialize_density_guess()
{
  // Use core Hamiltonian guess
  fock_matrix = kinetic_matrix + nuclear_matrix;
  solve_roothaan_hall();
  
  // Determine number of occupied orbitals
  n_electrons = 0;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    n_electrons += static_cast<int>(atomic_charges[type[i]][type[i]]);
  }
  
  n_occupied_orbitals = n_electrons / 2;  // Assuming closed-shell
  
  compute_density_matrix();
}

/* ----------------------------------------------------------------------
   Compute DFT forces
------------------------------------------------------------------------- */

void PairDFT::compute_forces(int eflag, int vflag)
{
  // Hellmann-Feynman forces
  compute_hellmann_feynman_forces();
  
  // Pulay forces (for non-orthogonal basis)
  compute_pulay_forces();
}

/* ----------------------------------------------------------------------
   Compute Hellmann-Feynman forces
------------------------------------------------------------------------- */

void PairDFT::compute_hellmann_feynman_forces()
{
  double **f = atom->f;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  // Nuclear-nuclear repulsion gradient
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < nlocal; j++) {
      if (i != j) {
        double dx = x[i][0] - x[j][0];
        double dy = x[i][1] - x[j][1];
        double dz = x[i][2] - x[j][2];
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        
        if (r > 1e-10) {
          double Zi = atomic_charges[type[i]][type[i]];
          double Zj = atomic_charges[type[j]][type[j]];
          double force_mag = Zi * Zj / (r * r * r);
          
          f[i][0] += force_mag * dx;
          f[i][1] += force_mag * dy;
          f[i][2] += force_mag * dz;
        }
      }
    }
  }
  
  // Electronic contribution to forces
  // This requires gradient integrals
  std::vector<Eigen::MatrixXd> dV(3 * nlocal);
  integral_engine->compute_nuclear_gradient(dV);
  
  for (int i = 0; i < nlocal; i++) {
    for (int k = 0; k < 3; k++) {
      double force_component = (density_matrix_eigen.cwiseProduct(dV[3*i + k])).sum();
      f[i][k] -= 2.0 * force_component;  // Factor of 2 for closed shell
    }
  }
}

/* ----------------------------------------------------------------------
   Compute Pulay forces
------------------------------------------------------------------------- */

void PairDFT::compute_pulay_forces()
{
  double **f = atom->f;
  int nlocal = atom->nlocal;
  
  // Energy-weighted density matrix
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(n_basis_functions, n_basis_functions);
  
  for (int i = 0; i < n_occupied_orbitals; i++) {
    W += 2.0 * mo_energies(i) * mo_coefficients.col(i) * mo_coefficients.col(i).transpose();
  }
  
  // Compute overlap gradient
  std::vector<Eigen::MatrixXd> dS(3 * nlocal);
  integral_engine->compute_overlap_gradient(dS);
  
  // Pulay force contribution
  for (int i = 0; i < nlocal; i++) {
    for (int k = 0; k < 3; k++) {
      double force_component = -(W.cwiseProduct(dS[3*i + k])).sum();
      f[i][k] += force_component;
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDFT::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal pair_style dft command");
  
  // Parse: pair_style dft <functional> <basis.json> [options]
  functional_name = std::string(arg[0]);
  basis_file = std::string(arg[1]);
  
  // Check if basis file exists
  std::ifstream file(basis_file);
  if (!file.good()) {
    error->all(FLERR, "Cannot open basis set file");
  }
  file.close();
  
  // Parse functional name
  parse_functional_name(functional_name.c_str());
  
  // Parse additional options
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "grid") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing grid size");
      grid_size = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "grid_type") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing grid type");
      grid_type = std::string(arg[iarg + 1]);
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
    } else if (strcmp(arg[iarg], "dispersion") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing dispersion type");
      use_dispersion = true;
      dispersion_type = std::string(arg[iarg + 1]);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown pair_style dft option");
    }
  }
  
  // Set cutoff - for DFT we typically use a large cutoff
  cut_global = 20.0;  // Angstroms
  
  // Reset cutoffs that have been explicitly set
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}


/* ----------------------------------------------------------------------
   set coefficients for one or more type pairs
------------------------------------------------------------------------- */

void PairDFT::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) 
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double charge = utils::numeric(FLERR, arg[2], false, lmp);
  double vdw_radius = utils::numeric(FLERR, arg[3], false, lmp);
  
  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      atomic_charges[i][j] = charge;
      vdw_radii[i][j] = vdw_radius;
      cut[i][j] = cut_one;
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
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(offset, n + 1, n + 1, "pair:offset");
  memory->create(atomic_charges, n + 1, n + 1, "pair:atomic_charges");
  memory->create(vdw_radii, n + 1, n + 1, "pair:vdw_radii");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDFT::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style dft requires atom IDs");
  
  // Request standard neighbor list
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDFT::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    atomic_charges[i][j] = mix_energy(atomic_charges[i][i], atomic_charges[j][j],
                                      vdw_radii[i][i], vdw_radii[j][j]);
    vdw_radii[i][j] = mix_distance(vdw_radii[i][i], vdw_radii[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  atomic_charges[j][i] = atomic_charges[i][j];
  vdw_radii[j][i] = vdw_radii[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   cleanup LibXC functionals
------------------------------------------------------------------------- */

void PairDFT::cleanup_libxc()
{
  if (xc_func_x) {
    xc_func_end(xc_func_x);
    delete xc_func_x;
    xc_func_x = nullptr;
  }
  if (xc_func_c) {
    xc_func_end(xc_func_c);
    delete xc_func_c;
    xc_func_c = nullptr;
  }
  if (xc_func_xc) {
    xc_func_end(xc_func_xc);
    delete xc_func_xc;
    xc_func_xc = nullptr;
  }
}



/* ----------------------------------------------------------------------
   compute one-electron integrals
------------------------------------------------------------------------- */

void PairDFT::compute_one_electron_integrals()
{
  // Overlap integrals
  integral_engine->compute_overlap(overlap_matrix);
  
  // Kinetic energy integrals
  integral_engine->compute_kinetic(kinetic_matrix);
  
  // Nuclear attraction integrals
  std::vector<double> charges;
  std::vector<std::vector<double>> positions;
  
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    charges.push_back(atomic_charges[type[i]][type[i]]);
    positions.push_back({x[i][0], x[i][1], x[i][2]});
  }
  
  integral_engine->compute_nuclear(nuclear_matrix, charges, positions);
}

/* ----------------------------------------------------------------------
   compute exact exchange for range-separated functionals
------------------------------------------------------------------------- */

void PairDFT::compute_exact_exchange()
{
  // For range-separated functionals, we need to compute
  // both short-range and long-range exchange
  
  // This requires specialized ERIs with error function
  // erf(omega * r12) / r12 for long-range
  // erfc(omega * r12) / r12 for short-range
  
  // Placeholder for now
  exact_exchange_energy = 0.0;
}

/* ----------------------------------------------------------------------
   SCF output methods
------------------------------------------------------------------------- */

void PairDFT::print_scf_header()
{
  if (comm->me == 0) {
    utils::logmesg(lmp, "\n");
    utils::logmesg(lmp, "================================================\n");
    utils::logmesg(lmp, "            DFT SCF CALCULATION\n");
    utils::logmesg(lmp, "================================================\n");
    utils::logmesg(lmp, fmt::format("Functional: {}\n", functional_name));
    utils::logmesg(lmp, fmt::format("Basis functions: {}\n", n_basis_functions));
    utils::logmesg(lmp, fmt::format("Electrons: {}\n", n_electrons));
    utils::logmesg(lmp, fmt::format("Grid points: {}\n", grid_size));
    utils::logmesg(lmp, "------------------------------------------------\n");
    utils::logmesg(lmp, " Iter    Energy         Delta E      RMS Dens\n");
    utils::logmesg(lmp, "------------------------------------------------\n");
  }
}

void PairDFT::print_scf_iteration()
{
  if (comm->me == 0) {
    double energy_change = 0.0;
    double density_change = density_matrix->get_change();
    
    utils::logmesg(lmp, fmt::format("{:4d} {:15.8f} {:12.5e} {:12.5e}\n",
                                    current_iteration, total_dft_energy,
                                    energy_change, density_change));
  }
}

void PairDFT::print_scf_summary()
{
  if (comm->me == 0) {
    utils::logmesg(lmp, "------------------------------------------------\n");
    if (scf_converged) {
      utils::logmesg(lmp, "SCF CONVERGED\n");
    } else {
      utils::logmesg(lmp, "SCF NOT CONVERGED\n");
    }
    utils::logmesg(lmp, "\nEnergy Components:\n");
    utils::logmesg(lmp, fmt::format("  Kinetic:         {:15.8f}\n", kinetic_energy));
    utils::logmesg(lmp, fmt::format("  Nuclear:         {:15.8f}\n", nuclear_repulsion));
    utils::logmesg(lmp, fmt::format("  Electron-Nuclear:{:15.8f}\n", electron_nuclear_energy));
    utils::logmesg(lmp, fmt::format("  Coulomb:         {:15.8f}\n", electron_electron_energy));
    utils::logmesg(lmp, fmt::format("  XC:              {:15.8f}\n", xc_energy));
    if (use_dispersion) {
      utils::logmesg(lmp, fmt::format("  Dispersion:      {:15.8f}\n", dispersion_energy));
    }
    utils::logmesg(lmp, fmt::format("  TOTAL:           {:15.8f}\n", total_dft_energy));
    utils::logmesg(lmp, "================================================\n\n");
  }
}

/* ----------------------------------------------------------------------
   DensityMatrix implementation
------------------------------------------------------------------------- */

DensityMatrix::DensityMatrix(int nbasis) : n_basis(nbasis)
{
  current.setZero(n_basis, n_basis);
  previous.setZero(n_basis, n_basis);
  difference.setZero(n_basis, n_basis);
  
  use_diis = true;
  diis_size = 6;
  
  diis_fock.clear();
  diis_error.clear();
}

DensityMatrix::~DensityMatrix()
{
}

void DensityMatrix::initialize_guess(const Eigen::MatrixXd &S)
{
  // Simple guess: set diagonal elements proportional to overlap
  current.setZero();
  
  for (int i = 0; i < n_basis; i++) {
    current(i, i) = 1.0 / sqrt(S(i, i));
  }
  
  // Normalize
  double trace = (current * S).trace();
  if (trace > 0) {
    current /= trace;
  }
  
  previous = current;
}

void DensityMatrix::update_from_mo(const Eigen::MatrixXd &C, int nocc)
{
  // Save previous density
  previous = current;
  
  // Build new density matrix from occupied orbitals
  // D = 2 * C_occ * C_occ^T for closed shell
  current.setZero();
  
  for (int i = 0; i < nocc; i++) {
    current += 2.0 * C.col(i) * C.col(i).transpose();
  }
  
  // Compute difference for convergence check
  difference = current - previous;
}

void DensityMatrix::mix_with_previous(double mixing_param)
{
  // Simple linear mixing for stability
  // D_new = (1-alpha)*D_old + alpha*D_current
  current = mixing_param * current + (1.0 - mixing_param) * previous;
}

double DensityMatrix::get_change() const
{
  // RMS change in density matrix
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < n_basis; i++) {
    for (int j = 0; j < n_basis; j++) {
      sum += difference(i, j) * difference(i, j);
      count++;
    }
  }
  
  return sqrt(sum / count);
}

std::vector<double> DensityMatrix::compute_mulliken_charges(const Eigen::MatrixXd &S)
{
  // Mulliken population analysis
  // q_A = Z_A - sum_mu(P_mu,mu * S_mu,mu) for mu on atom A
  
  // For now, return empty vector (need atom mapping)
  std::vector<double> charges;
  
  // TODO: Implement with proper basis-to-atom mapping
  
  return charges;
}

std::vector<double> DensityMatrix::compute_lowdin_charges(const Eigen::MatrixXd &S)
{
  // Löwdin population analysis
  // Uses S^(1/2) transformation
  
  // Compute S^(1/2)
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  Eigen::MatrixXd S_sqrt = es.operatorSqrt();
  
  // Transform density matrix
  Eigen::MatrixXd P_lowdin = S_sqrt * current * S_sqrt;
  
  // For now, return empty vector (need atom mapping)
  std::vector<double> charges;
  
  // TODO: Implement with proper basis-to-atom mapping
  
  return charges;
}

void DensityMatrix::apply_diis(Eigen::MatrixXd &F)
{
  if (!use_diis) return;
  
  // DIIS (Direct Inversion of Iterative Subspace)
  // Accelerates SCF convergence
  
  // Store current Fock matrix and error
  diis_fock.push_back(F);
  
  // Compute error matrix: e = FDS - SDF
  // For simplicity, using gradient of energy wrt density
  Eigen::MatrixXd error = F * current - current * F;
  diis_error.push_back(error);
  
  // Keep only last diis_size iterations
  if (diis_fock.size() > static_cast<size_t>(diis_size)) {
    diis_fock.erase(diis_fock.begin());
    diis_error.erase(diis_error.begin());
  }
  
  int n = diis_fock.size();
  if (n < 2) return;  // Need at least 2 iterations
  
  // Build B matrix
  Eigen::MatrixXd B(n + 1, n + 1);
  B.setZero();
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      double val = (diis_error[i].cwiseProduct(diis_error[j])).sum();
      B(i, j) = val;
      B(j, i) = val;
    }
    B(i, n) = -1.0;
    B(n, i) = -1.0;
  }
  B(n, n) = 0.0;
  
  // Solve for coefficients
  Eigen::VectorXd rhs(n + 1);
  rhs.setZero();
  rhs(n) = -1.0;
  
  Eigen::VectorXd c = B.colPivHouseholderQr().solve(rhs);
  
  // Build extrapolated Fock matrix
  F.setZero();
  for (int i = 0; i < n; i++) {
    F += c(i) * diis_fock[i];
  }
}

/* ----------------------------------------------------------------------
   GridIntegrator implementation
------------------------------------------------------------------------- */

GridIntegrator::GridIntegrator(int grid_size, const std::string &type) 
  : target_grid_size(grid_size), grid_type(type)
{
  grid_points.clear();
  grid_weights.clear();
  becke_weights.clear();
}

GridIntegrator::~GridIntegrator()
{
}

void GridIntegrator::generate_grid(const std::vector<std::vector<double>> &atom_positions,
                                   const std::vector<int> &atomic_numbers)
{
  grid_points.clear();
  grid_weights.clear();
  
  int n_atoms = atom_positions.size();
  if (n_atoms == 0) return;
  
  // Generate atomic grids
  int points_per_atom = target_grid_size / n_atoms;
  
  for (int atom = 0; atom < n_atoms; atom++) {
    // Get radial and angular grid sizes
    int n_radial = 50;  // Typical value
    int n_angular = points_per_atom / n_radial;
    
    // Generate radial grid (Chebyshev-Gauss)
    std::vector<double> r_points(n_radial);
    std::vector<double> r_weights(n_radial);
    generate_radial_grid(n_radial, atomic_numbers[atom], r_points, r_weights);
    
    // Generate angular grid (Lebedev)
    std::vector<std::vector<double>> angular_points;
    std::vector<double> angular_weights;
    generate_lebedev_grid(n_angular, angular_points, angular_weights);
    
    // Combine radial and angular grids
    for (int i_r = 0; i_r < n_radial; i_r++) {
      double r = r_points[i_r];
      double w_r = r_weights[i_r];
      
      for (size_t i_ang = 0; i_ang < angular_points.size(); i_ang++) {
        std::vector<double> point(3);
        point[0] = atom_positions[atom][0] + r * angular_points[i_ang][0];
        point[1] = atom_positions[atom][1] + r * angular_points[i_ang][1];
        point[2] = atom_positions[atom][2] + r * angular_points[i_ang][2];
        
        grid_points.push_back(point);
        grid_weights.push_back(w_r * angular_weights[i_ang] * r * r);
      }
    }
  }
  
  // Compute Becke partitioning weights
  compute_becke_weights(atom_positions);
}

void GridIntegrator::generate_radial_grid(int n_points, double Z,
                                          std::vector<double> &r,
                                          std::vector<double> &w)
{
  // Chebyshev-Gauss radial grid
  // Transform from [-1,1] to [0,inf) using appropriate mapping
  
  r.resize(n_points);
  w.resize(n_points);
  
  // Bragg radius for scaling
  double R_bragg = 1.0;  // Default, should depend on Z
  if (Z > 0) {
    R_bragg = 0.5 * (3.0 - 0.01 * Z);  // Simple approximation
  }
  
  for (int i = 0; i < n_points; i++) {
    // Chebyshev nodes
    double xi = cos(M_PI * (i + 0.5) / n_points);
    
    // Becke transformation
    double x = (1.0 + xi) / (1.0 - xi);
    r[i] = R_bragg * x;
    
    // Weight includes Jacobian
    w[i] = (M_PI / n_points) * 2.0 * R_bragg / ((1.0 - xi) * (1.0 - xi));
  }
}

void GridIntegrator::generate_lebedev_grid(int n_points,
                                           std::vector<std::vector<double>> &points,
                                           std::vector<double> &weights)
{
  // Simplified Lebedev grid generation
  // Real implementation would use pre-computed Lebedev grids
  
  points.clear();
  weights.clear();
  
  // Find closest available Lebedev grid
  // Available sizes: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, ...
  int actual_points = 50;  // Default to 50-point grid
  
  if (n_points <= 6) actual_points = 6;
  else if (n_points <= 14) actual_points = 14;
  else if (n_points <= 26) actual_points = 26;
  else if (n_points <= 38) actual_points = 38;
  else if (n_points <= 50) actual_points = 50;
  else if (n_points <= 74) actual_points = 74;
  else if (n_points <= 110) actual_points = 110;
  else if (n_points <= 170) actual_points = 170;
  else actual_points = 194;
  
  // Generate points on unit sphere
  // This is a simplified version - real Lebedev grids are more complex
  
  if (actual_points == 6) {
    // Octahedron vertices
    points.push_back({1, 0, 0});
    points.push_back({-1, 0, 0});
    points.push_back({0, 1, 0});
    points.push_back({0, -1, 0});
    points.push_back({0, 0, 1});
    points.push_back({0, 0, -1});
    
    double w = 4.0 * M_PI / 6.0;
    for (int i = 0; i < 6; i++) {
      weights.push_back(w);
    }
  } else {
    // Use spherical Fibonacci grid as approximation
    double phi = (1.0 + sqrt(5.0)) / 2.0;  // Golden ratio
    
    for (int i = 0; i < actual_points; i++) {
      double y = 1.0 - 2.0 * i / (actual_points - 1.0);
      double radius = sqrt(1.0 - y * y);
      double theta = 2.0 * M_PI * i / phi;
      
      std::vector<double> point(3);
      point[0] = radius * cos(theta);
      point[1] = radius * sin(theta);
      point[2] = y;
      
      points.push_back(point);
      weights.push_back(4.0 * M_PI / actual_points);
    }
  }
}

void GridIntegrator::compute_becke_weights(const std::vector<std::vector<double>> &atoms)
{
  // Becke partitioning for multi-center integration
  int n_points = grid_points.size();
  int n_atoms = atoms.size();
  
  becke_weights.resize(n_points);
  
  for (int i_point = 0; i_point < n_points; i_point++) {
    std::vector<double> P(n_atoms, 1.0);
    
    // Compute partition function for each atom
    for (int i_atom = 0; i_atom < n_atoms; i_atom++) {
      for (int j_atom = 0; j_atom < n_atoms; j_atom++) {
        if (i_atom == j_atom) continue;
        
        // Distance from point to atoms i and j
        double r_i = 0.0, r_j = 0.0;
        for (int k = 0; k < 3; k++) {
          double d_i = grid_points[i_point][k] - atoms[i_atom][k];
          double d_j = grid_points[i_point][k] - atoms[j_atom][k];
          r_i += d_i * d_i;
          r_j += d_j * d_j;
        }
        r_i = sqrt(r_i);
        r_j = sqrt(r_j);
        
        // Distance between atoms i and j
        double R_ij = 0.0;
        for (int k = 0; k < 3; k++) {
          double d = atoms[i_atom][k] - atoms[j_atom][k];
          R_ij += d * d;
        }
        R_ij = sqrt(R_ij);
        
        if (R_ij < 1e-10) continue;
        
        // Confocal elliptical coordinate
        double mu = (r_i - r_j) / R_ij;
        
        // Becke's step function with smoothing
        double f = mu;
        for (int iter = 0; iter < 3; iter++) {
          f = 0.5 * f * (3.0 - f * f);
        }
        double s = 0.5 * (1.0 - f);
        
        P[i_atom] *= s;
      }
    }
    
    // Normalize partition functions
    double sum = 0.0;
    for (int i_atom = 0; i_atom < n_atoms; i_atom++) {
      sum += P[i_atom];
    }
    
    if (sum > 1e-15) {
      // For simplicity, assign full weight to nearest atom
      // Real implementation would properly partition weights
      becke_weights[i_point] = 1.0;
    } else {
      becke_weights[i_point] = 0.0;
    }
  }
}

void GridIntegrator::integrate_xc(XCFunctional *xc_func,
                                  const Eigen::MatrixXd &density_matrix,
                                  BasisSetManager *basis,
                                  double &exc_energy,
                                  Eigen::MatrixXd &vxc_matrix)
{
  int n_basis = basis->get_n_basis();
  int n_points = grid_points.size();
  
  vxc_matrix.setZero(n_basis, n_basis);
  exc_energy = 0.0;
  
  if (n_points == 0) return;
  
  // Evaluate basis functions at grid points
  std::vector<std::vector<double>> basis_values(n_points, std::vector<double>(n_basis));
  std::vector<std::vector<std::vector<double>>> basis_gradients;
  
  if (xc_func->is_gga() || xc_func->is_meta()) {
    basis_gradients.resize(n_points, 
                          std::vector<std::vector<double>>(n_basis, 
                                                          std::vector<double>(3)));
  }
  
  evaluate_basis_at_points(basis, basis_values, basis_gradients);
  
  // Compute density and gradients at grid points
  std::vector<double> rho(n_points, 0.0);
  std::vector<double> sigma(n_points, 0.0);  // |grad rho|^2
  std::vector<double> lapl(n_points, 0.0);   // Laplacian
  std::vector<double> tau(n_points, 0.0);     // Kinetic energy density
  
  for (int i_point = 0; i_point < n_points; i_point++) {
    // Density
    for (int i = 0; i < n_basis; i++) {
      for (int j = 0; j < n_basis; j++) {
        rho[i_point] += density_matrix(i, j) * 
                        basis_values[i_point][i] * 
                        basis_values[i_point][j];
      }
    }
    
    // Gradient and kinetic energy density for GGA/meta-GGA
    if (xc_func->is_gga() || xc_func->is_meta()) {
      std::vector<double> grad_rho(3, 0.0);
      
      for (int i = 0; i < n_basis; i++) {
        for (int j = 0; j < n_basis; j++) {
          double P_ij = density_matrix(i, j);
          
          for (int k = 0; k < 3; k++) {
            grad_rho[k] += P_ij * (basis_gradients[i_point][i][k] * basis_values[i_point][j] +
                                   basis_values[i_point][i] * basis_gradients[i_point][j][k]);
          }
          
          if (xc_func->is_meta()) {
            for (int k = 0; k < 3; k++) {
              tau[i_point] += 0.5 * P_ij * 
                             basis_gradients[i_point][i][k] * 
                             basis_gradients[i_point][j][k];
            }
          }
        }
      }
      
      sigma[i_point] = grad_rho[0] * grad_rho[0] + 
                       grad_rho[1] * grad_rho[1] + 
                       grad_rho[2] * grad_rho[2];
    }
  }
  
  // Evaluate XC functional
  std::vector<double> exc(n_points);
  std::vector<double> vrho(n_points);
  std::vector<double> vsigma(n_points);
  std::vector<double> vlapl(n_points);
  std::vector<double> vtau(n_points);
  
  xc_func->evaluate(rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);
  
  // Integrate XC energy
  for (int i_point = 0; i_point < n_points; i_point++) {
    exc_energy += exc[i_point] * rho[i_point] * 
                  grid_weights[i_point] * becke_weights[i_point];
  }
  
  // Build XC potential matrix
  for (int i_point = 0; i_point < n_points; i_point++) {
    double w = grid_weights[i_point] * becke_weights[i_point];
    
    // LDA contribution
    for (int i = 0; i < n_basis; i++) {
      for (int j = 0; j <= i; j++) {
        double val = vrho[i_point] * basis_values[i_point][i] * 
                    basis_values[i_point][j] * w;
        vxc_matrix(i, j) += val;
        if (i != j) vxc_matrix(j, i) += val;
      }
    }
    
    // GGA contribution
    if (xc_func->is_gga() && vsigma[i_point] != 0.0) {
      // Simplified - full implementation would include gradient contributions
      // This requires second derivatives of basis functions
    }
    
    // Meta-GGA contribution
    if (xc_func->is_meta() && vtau[i_point] != 0.0) {
      // Simplified - full implementation would include tau contributions
    }
  }
}

void GridIntegrator::evaluate_basis_at_points(BasisSetManager *basis,
                                              std::vector<std::vector<double>> &basis_values,
                                              std::vector<std::vector<std::vector<double>>> &basis_gradients)
{
  int n_points = grid_points.size();
  int n_basis = basis->get_n_basis();
  
  // For each grid point, evaluate all basis functions
  for (int i_point = 0; i_point < n_points; i_point++) {
    basis->evaluate_basis(grid_points[i_point][0], 
                         grid_points[i_point][1], 
                         grid_points[i_point][2],
                         basis_values[i_point]);
    
    if (!basis_gradients.empty()) {
      std::vector<double> dummy_values;
      basis->evaluate_basis_gradient(grid_points[i_point][0], 
                                     grid_points[i_point][1], 
                                     grid_points[i_point][2],
                                     dummy_values,
                                     basis_gradients[i_point]);
    }
  }
}
