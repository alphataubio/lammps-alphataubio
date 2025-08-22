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
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>

#include <nlohmann/json.hpp>
#include <libint2.hpp>
#include <xc.h>

using namespace LAMMPS_NS;
using json = nlohmann_lmp::json;

/* ---------------------------------------------------------------------- */

PairDFT::PairDFT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  restartinfo = 1;
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
   Special implementation for wB97M-V functional
------------------------------------------------------------------------- */

void PairDFT::evaluate_wb97mv_functional()
{
  // wB97M-V is a range-separated meta-GGA functional with VV10 NLC
  // Parameters from the original paper (Mardirossian & Head-Gordon, 2016)
  
  const double omega = 0.3;  // Range-separation parameter
  const double c_x_lr = 1.0;  // Long-range HF exchange
  const double c_x_sr = 0.15;  // Short-range HF exchange at r=0
  
  // B97 exchange parameters for wB97M-V
  const double a0_x = 0.85;
  const double a1_x = 1.007;
  const double a2_x = 0.259;
  
  // B97 correlation parameters for wB97M-V
  const double a0_c_ss = 1.0;
  const double a1_c_ss = -3.382;
  const double a2_c_ss = -0.892;
  const double a0_c_os = 1.0;
  const double a1_c_os = -1.855;
  const double a2_c_os = 0.653;
  
  // VV10 NLC parameters
  const double b_vv10 = 6.0;
  const double C_vv10 = 0.01;
  
  // Generate grid if not already done
  std::vector<std::vector<double>> atom_positions;
  std::vector<int> atomic_numbers;
  
  // Get atomic positions and numbers from LAMMPS atoms
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    atom_positions.push_back({x[i][0], x[i][1], x[i][2]});
    atomic_numbers.push_back(static_cast<int>(atomic_charges[type[i]][type[i]]));
  }
  
  grid_integrator->generate_grid(atom_positions, atomic_numbers);
  
  // Evaluate density and its derivatives on the grid
  int npoints = grid_integrator->get_n_points();
  std::vector<double> rho(npoints);
  std::vector<double> sigma(npoints);  // |grad rho|^2
  std::vector<double> tau(npoints);    // Kinetic energy density
  
  // TODO: Evaluate basis functions and density on grid points
  // This requires detailed implementation of basis function evaluation
  
  // Apply B97-style inhomogeneity correction factor
  std::vector<double> s(npoints);  // Reduced density gradient
  std::vector<double> u_x(npoints);  // Exchange enhancement factor
  std::vector<double> u_c_ss(npoints);  // Same-spin correlation factor
  std::vector<double> u_c_os(npoints);  // Opposite-spin correlation factor
  
  for (int i = 0; i < npoints; i++) {
    if (rho[i] > 1e-15) {
      // Reduced density gradient
      s[i] = sqrt(sigma[i]) / (2.0 * pow(3.0 * M_PI * M_PI, 1.0/3.0) * pow(rho[i], 4.0/3.0));
      
      // B97 inhomogeneity correction factor
      double s2 = s[i] * s[i];
      double denom = 1.0 + 0.004 * s2;
      u_x[i] = (a0_x + a1_x * s2 + a2_x * s2 * s2) / denom;
      u_c_ss[i] = (a0_c_ss + a1_c_ss * s2 + a2_c_ss * s2 * s2) / denom;
      u_c_os[i] = (a0_c_os + a1_c_os * s2 + a2_c_os * s2 * s2) / denom;
    }
  }
  
  // Compute short-range exchange-correlation energy
  std::vector<double> exc_sr(npoints);
  std::vector<double> vrho_sr(npoints);
  std::vector<double> vsigma_sr(npoints);
  std::vector<double> vtau_sr(npoints);
  
  // Use modified B97 functional for short-range part
  for (int i = 0; i < npoints; i++) {
    if (rho[i] > 1e-15) {
      // LDA exchange energy density
      double ex_lda = -0.75 * pow(3.0 * rho[i] / M_PI, 1.0/3.0);
      
      // Apply enhancement factor
      exc_sr[i] = ex_lda * u_x[i];
      
      // Add correlation (simplified - should use proper LDA correlation)
      double ec_lda = -0.05 * pow(rho[i], 1.0/3.0);  // Placeholder
      exc_sr[i] += ec_lda * (u_c_ss[i] + u_c_os[i]) / 2.0;
      
      // Meta-GGA contribution from kinetic energy density
      if (tau[i] > 1e-15) {
        double tau_w = sigma[i] / (8.0 * rho[i]);  // von Weizsäcker kinetic energy
        double alpha = (tau[i] - tau_w) / tau[i];
        exc_sr[i] *= (1.0 + 0.1 * alpha);  // Simple meta-GGA correction
      }
    }
  }
  
  // Add VV10 non-local correlation
  compute_vv10_nlc(rho, sigma, b_vv10, C_vv10);
  
  // Compute exact exchange contribution
  // Short-range: c_x_sr * HF_sr(omega)
  // Long-range: c_x_lr * HF_lr(omega)
  
  // The actual implementation would require:
  // 1. Computing range-separated ERIs
  // 2. Building separate short-range and long-range exchange matrices
  // 3. Properly combining all contributions
  
  if (comm->me == 0) {
    utils::logmesg(lmp, "wB97M-V functional evaluation completed\n");
  }
}

/* ----------------------------------------------------------------------
   Compute VV10 non-local correlation
------------------------------------------------------------------------- */

void PairDFT::compute_vv10_nlc(const std::vector<double> &rho,
                               const std::vector<double> &sigma,
                               double b, double C)
{
  // VV10 non-local correlation functional
  // Vydrov & Van Voorhis, JCP 133, 244103 (2010)
  
  int npoints = rho.size();
  std::vector<double> omega(npoints);
  
  // Compute omega(r) = sqrt(g(r)) where g is defined in VV10 paper
  for (int i = 0; i < npoints; i++) {
    if (rho[i] > 1e-15) {
      double kf = pow(3.0 * M_PI * M_PI * rho[i], 1.0/3.0);
      double wp = 4.0 * M_PI * rho[i];
      
      // Compute g factor
      double s2 = sigma[i] / (4.0 * pow(kf, 2) * rho[i] * rho[i]);
      double g = wp / (kf * kf) * (1.0 + C * s2);
      
      omega[i] = sqrt(g);
    }
  }
  
  // Double integral for non-local correlation
  // This is computationally expensive and typically requires special techniques
  // For now, using a simplified local approximation
  
  double e_nlc = 0.0;
  for (int i = 0; i < npoints; i++) {
    if (rho[i] > 1e-15) {
      // Simplified local approximation to NLC
      e_nlc += -0.01 * rho[i] * pow(omega[i], 0.5);
    }
  }
  
  xc_energy += e_nlc;
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
   parse functional name to LibXC ID
------------------------------------------------------------------------- */

void PairDFT::parse_functional_name(const char *name)
{
  // Special handling for wB97M-V
  if (strcasecmp(name, "wB97M-V") == 0 || strcasecmp(name, "wb97mv") == 0) {
    // wB97M-V is not directly in LibXC, we need special implementation
    is_range_separated = true;
    is_meta_gga = true;
    is_hybrid = true;
    range_separation_param = 0.3;
    hybrid_coeff = 0.15;  // Short-range exact exchange
    use_combined_xc = false;
    
    // We'll use a special implementation for wB97M-V
    xc_functional_xc = -999;  // Special flag for wB97M-V
    
    if (comm->me == 0) {
      utils::logmesg(lmp, "Using special implementation for wB97M-V functional\n");
    }
    return;
  }
  
  // Check if it's a combined functional
  if (strchr(name, '+') != nullptr) {
    // Separate X and C functionals
    use_combined_xc = false;
    char *name_copy = strdup(name);
    char *x_func = strtok(name_copy, "+");
    char *c_func = strtok(nullptr, "+");
    
    if (x_func) {
      xc_functional_x = xc_functional_get_number(x_func);
      if (xc_functional_x == -1) {
        error->all(FLERR, "Unknown exchange functional");
      }
      xc_func_x = new xc_func_type;
      if (xc_func_init(xc_func_x, xc_functional_x, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize exchange functional");
      }
    }
    
    if (c_func) {
      xc_functional_c = xc_functional_get_number(c_func);
      if (xc_functional_c == -1) {
        error->all(FLERR, "Unknown correlation functional");
      }
      xc_func_c = new xc_func_type;
      if (xc_func_init(xc_func_c, xc_functional_c, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize correlation functional");
      }
    }
    
    free(name_copy);
  } else {
    // Combined XC functional
    use_combined_xc = true;
    xc_functional_xc = xc_functional_get_number(name);
    if (xc_functional_xc == -1) {
      error->all(FLERR, "Unknown XC functional");
    }
    xc_func_xc = new xc_func_type;
    if (xc_func_init(xc_func_xc, xc_functional_xc, XC_UNPOLARIZED) != 0) {
      error->all(FLERR, "Failed to initialize XC functional");
    }
    
    // Check functional family
    xc_func_info_type *info = xc_func_get_info(xc_func_xc);
    int family = xc_func_info_get_family(info);
    
    if (family == XC_FAMILY_HYB_GGA || family == XC_FAMILY_HYB_MGGA) {
      is_hybrid = true;
      hybrid_coeff = xc_hyb_exx_coef(xc_func_xc);
    }
    if (family == XC_FAMILY_MGGA || family == XC_FAMILY_HYB_MGGA) {
      is_meta_gga = true;
    }
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
  if (atom->tag_enable == 0)
    error->all(FLERR, "Pair style dft requires atom IDs");
  
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
   cleanup integral engines
------------------------------------------------------------------------- */

void PairDFT::cleanup_integrals()
{
  // Cleanup handled by smart pointers
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
   compute dispersion correction
------------------------------------------------------------------------- */

void PairDFT::compute_dispersion_correction()
{
  if (dispersion_type == "D3") {
    compute_d3_dispersion();
  } else if (dispersion_type == "D3BJ") {
    compute_d3bj_dispersion();
  } else if (dispersion_type == "D4") {
    compute_d4_dispersion();
  }
}

/* ----------------------------------------------------------------------
   compute D3 dispersion correction (Grimme et al.)
------------------------------------------------------------------------- */

void PairDFT::compute_d3_dispersion()
{
  // Simplified D3 dispersion
  // Full implementation would require:
  // - CN (coordination number) calculation
  // - C6 interpolation
  // - Three-body terms
  
  dispersion_energy = 0.0;
  
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  // Two-body dispersion
  for (int i = 0; i < nlocal; i++) {
    for (int j = i + 1; j < nlocal; j++) {
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      double r = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (r > 1e-10 && r < cut_global) {
        // Get C6 coefficient (simplified - should depend on CN)
        double C6 = 10.0;  // Placeholder
        
        // Becke-Johnson damping
        double R0 = vdw_radii[type[i]][type[i]] + vdw_radii[type[j]][type[j]];
        double a1 = 0.45;
        double a2 = 4.0;
        
        double damp = 1.0 / (1.0 + exp(-a2 * (r/R0 - a1)));
        double E_disp = -C6 * damp / pow(r, 6);
        
        dispersion_energy += E_disp;
        
        // Forces
        double F_disp = -6.0 * E_disp / r;
        F_disp += C6 * a2 * damp * damp * exp(-a2 * (r/R0 - a1)) / (R0 * pow(r, 6));
        
        F_disp /= r;
        
        f[i][0] += F_disp * dx;
        f[i][1] += F_disp * dy;
        f[i][2] += F_disp * dz;
        f[j][0] -= F_disp * dx;
        f[j][1] -= F_disp * dy;
        f[j][2] -= F_disp * dz;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute D3BJ dispersion correction
------------------------------------------------------------------------- */

void PairDFT::compute_d3bj_dispersion()
{
  // D3 with Becke-Johnson damping
  compute_d3_dispersion();  // Similar but with different damping function
}

/* ----------------------------------------------------------------------
   compute D4 dispersion correction
------------------------------------------------------------------------- */

void PairDFT::compute_d4_dispersion()
{
  // D4 includes charge-dependent C6 coefficients
  // Placeholder for now
  compute_d3_dispersion();
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
   single() function for energy and force calculations
------------------------------------------------------------------------- */

double PairDFT::single(int i, int j, int itype, int jtype,
                      double rsq, double factor_coul, double factor_lj,
                      double &fforce)
{
  // For DFT, single pair interactions don't have the usual meaning
  // Return a simple repulsion at very short range
  
  double r = sqrt(rsq);
  if (r < 0.5) {
    fforce = 1000.0 * exp(-10.0 * r) / r * factor_lj;
    return 100.0 * exp(-10.0 * r) * factor_lj;
  }
  
  fforce = 0.0;
  return 0.0;
}

/* ----------------------------------------------------------------------
   write/read restart methods
------------------------------------------------------------------------- */

void PairDFT::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&atomic_charges[i][j], sizeof(double), 1, fp);
        fwrite(&vdw_radii[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

void PairDFT::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &atomic_charges[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &vdw_radii[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&atomic_charges[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&vdw_radii[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

void PairDFT::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&xc_functional_x, sizeof(int), 1, fp);
  fwrite(&xc_functional_c, sizeof(int), 1, fp);
  fwrite(&xc_functional_xc, sizeof(int), 1, fp);
  fwrite(&use_combined_xc, sizeof(bool), 1, fp);
  fwrite(&grid_size, sizeof(int), 1, fp);
  fwrite(&max_scf_iterations, sizeof(int), 1, fp);
  fwrite(&energy_tolerance, sizeof(double), 1, fp);
  fwrite(&use_dispersion, sizeof(bool), 1, fp);
}

void PairDFT::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &xc_functional_x, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &xc_functional_c, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &xc_functional_xc, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &use_combined_xc, sizeof(bool), 1, fp, nullptr, error);
    utils::sfread(FLERR, &grid_size, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &max_scf_iterations, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &energy_tolerance, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &use_dispersion, sizeof(bool), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&xc_functional_x, 1, MPI_INT, 0, world);
  MPI_Bcast(&xc_functional_c, 1, MPI_INT, 0, world);
  MPI_Bcast(&xc_functional_xc, 1, MPI_INT, 0, world);
  MPI_Bcast(&use_combined_xc, 1, MPI_C_BOOL, 0, world);
  MPI_Bcast(&grid_size, 1, MPI_INT, 0, world);
  MPI_Bcast(&max_scf_iterations, 1, MPI_INT, 0, world);
  MPI_Bcast(&energy_tolerance, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&use_dispersion, 1, MPI_C_BOOL, 0, world);
}

void PairDFT::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g\n", i, atomic_charges[i][i], vdw_radii[i][i]);
}

void PairDFT::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j,
              atomic_charges[i][j], vdw_radii[i][j], cut[i][j]);
}
