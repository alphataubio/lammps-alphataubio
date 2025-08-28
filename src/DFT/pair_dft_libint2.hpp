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

#include <libint2.hpp>
#include <libint2/basis.h>
#include <libint2/shell.h>
#include <libint2/engine.h>
#include <cmath>
#include <vector>
#include <algorithm>

// Type alias for libint2's small_vector
template<typename T>
using svector = libint2::svector<T>;

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
  }
  
  // Add XC contribution via grid integration
  Eigen::MatrixXd vxc_matrix(n_basis_functions, n_basis_functions);
  
  // Get atomic numbers from atom->q
  std::vector<int> atomic_numbers;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    atomic_numbers.push_back(static_cast<int>(q[i]));
  }
  
  // Generate grid and integrate XC
  std::vector<std::vector<double>> positions;
  double **x = atom->x;
  for (int i = 0; i < nlocal; i++) {
    positions.push_back({x[i][0], x[i][1], x[i][2]});
  }
  
  grid_integrator->generate_grid(positions, atomic_numbers);
  grid_integrator->integrate_xc(xc_functional.get(), density_matrix_eigen,
                                basis_manager.get(), xc_energy, vxc_matrix);
  
  fock_matrix += vxc_matrix;
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
  
  // Kinetic energy component
  kinetic_energy = (density_matrix_eigen.cwiseProduct(kinetic_matrix)).sum();
  
  // Total energy
  total_dft_energy = nuclear_repulsion + one_electron + j_energy - k_energy + xc_energy;
}

/* ----------------------------------------------------------------------
   Compute nuclear repulsion energy
------------------------------------------------------------------------- */

double PairDFT::compute_nuclear_repulsion_energy()
{
  double energy = 0.0;
  
  double **x = atom->x;
  double *q = atom->q;  // Use atom->q for nuclear charges
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    for (int j = i + 1; j < nlocal; j++) {
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      double r = sqrt(dx*dx + dy*dy + dz*dz);
      if (r > 1e-10) energy += q[i] * q[j] / r;
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
  double *q = atom->q;  // Use atom->q
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) n_electrons += static_cast<int>(q[i]);
  
  n_occupied_orbitals = n_electrons / 2;  // Assuming closed-shell
  
  compute_density_matrix();
}

/* ----------------------------------------------------------------------
   Compute one-electron integrals
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
  double *q = atom->q;  // Use atom->q
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    charges.push_back(q[i]);
    positions.push_back({x[i][0], x[i][1], x[i][2]});
  }
  
  integral_engine->compute_nuclear(nuclear_matrix, charges, positions);
}


/* ----------------------------------------------------------------------
   Compute Hellmann-Feynman forces
------------------------------------------------------------------------- */

void PairDFT::compute_hellmann_feynman_forces()
{
  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;  // Use atom->q
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
          double force_mag = q[i] * q[j] / (r * r * r);
          f[i][0] += force_mag * dx;
          f[i][1] += force_mag * dy;
          f[i][2] += force_mag * dz;
        }
      }
    }
  }
  
  // Electronic contribution to forces
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
   Solve Roothaan-Hall equations
------------------------------------------------------------------------- */

void PairDFT::solve_roothaan_hall()
{
  // Transform Fock matrix to orthogonal basis
  // F' = S^(-1/2) * F * S^(-1/2)
  
  // Compute S^(-1/2) using eigendecomposition
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(overlap_matrix);
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


/* ---------------------------------------------------------------------- */

void PairDFT::initialize_libint()
{
  // Create libint2 basis set from our basis manager
  std::vector<libint2::Shell> shells;
  
  int n_shells = basis_set->get_n_shells();
  auto shell_to_atom = basis_set->get_shell_to_atom_map();
  
  for (int s = 0; s < n_shells; s++) {
    int l = basis_set->get_angular_momentum(s);
    auto exponents = basis_set->get_exponents(s);
    auto coefficients = basis_set->get_coefficients(s);
    
    // Convert std::vector to libint2's svector (small_vector)
    libint2::svector<double> alpha_svec;
    for (auto exp : exponents) {
      alpha_svec.push_back(exp);
    }
    
    // Create contraction with proper svector type
    libint2::Shell::Contraction contr;
    contr.l = l;
    contr.pure = false;  // Use Cartesian Gaussians
    
    // Convert coefficients to svector
    for (auto coeff : coefficients) {
      contr.coeff.push_back(coeff);
    }
    
    // Create svector of contractions
    libint2::svector<libint2::Shell::Contraction> contr_svec;
    contr_svec.push_back(contr);
    
    // Create shell using default constructor and then set members
    libint2::Shell shell;
    shell.alpha = alpha_svec;
    shell.contr = contr_svec;
    shell.O = {{0.0, 0.0, 0.0}};  // Origin - placeholder, should come from atom positions
    
    shells.push_back(shell);
  }
  
  libint_basis = std::make_unique<libint2::BasisSet>(shells);
  
  // Initialize engines for different integral types
  engines.resize(4);
  
  // Get max angular momentum and number of primitives
  int max_l = 0;
  size_t max_nprim = 0;
  for (const auto& shell : shells) {
    max_l = std::max(max_l, shell.contr[0].l);
    max_nprim = std::max(max_nprim, shell.alpha.size());
  }
  
  engines[0] = std::make_unique<libint2::Engine>(libint2::Operator::overlap, 
                                                 max_nprim, max_l);
  engines[1] = std::make_unique<libint2::Engine>(libint2::Operator::kinetic,
                                                 max_nprim, max_l);
  engines[2] = std::make_unique<libint2::Engine>(libint2::Operator::nuclear,
                                                 max_nprim, max_l);
  engines[3] = std::make_unique<libint2::Engine>(libint2::Operator::coulomb,
                                                 max_nprim, max_l);
}

/* ---------------------------------------------------------------------- */

void PairDFT::cleanup_libint()
{
  engines.clear();
  libint_basis.reset();
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_overlap(Eigen::MatrixXd &S)
{
  int n = basis_set->get_n_basis();
  S.setZero(n, n);
  
  auto& engine = *engines[0];
  const auto& shells = libint_basis->shells();
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int n1 = shells[s1].size();  // Number of basis functions in shell
    int bf2 = 0;
    
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      int n2 = shells[s2].size();
      
      engine.compute(shells[s1], shells[s2]);
      const auto* buf = engine.results()[0];
      
      if (buf != nullptr) {
        for (int f1 = 0; f1 < n1; ++f1) {
          for (int f2 = 0; f2 < n2; ++f2) {
            S(bf1 + f1, bf2 + f2) = buf[f1 * n2 + f2];
            if (s1 != s2) {
              S(bf2 + f2, bf1 + f1) = buf[f1 * n2 + f2];
            }
          }
        }
      }
      
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_kinetic(Eigen::MatrixXd &T)
{
  int n = basis_set->get_n_basis();
  T.setZero(n, n);
  
  auto& engine = *engines[1];
  const auto& shells = libint_basis->shells();
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int n1 = shells[s1].size();
    int bf2 = 0;
    
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      int n2 = shells[s2].size();
      
      engine.compute(shells[s1], shells[s2]);
      const auto* buf = engine.results()[0];
      
      if (buf != nullptr) {
        for (int f1 = 0; f1 < n1; ++f1) {
          for (int f2 = 0; f2 < n2; ++f2) {
            T(bf1 + f1, bf2 + f2) = buf[f1 * n2 + f2];
            if (s1 != s2) {
              T(bf2 + f2, bf1 + f1) = buf[f1 * n2 + f2];
            }
          }
        }
      }
      
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_nuclear(Eigen::MatrixXd &V, 
                                     const std::vector<double> &charges,
                                     const std::vector<std::vector<double>> &positions)
{
  int n = basis_set->get_n_basis();
  V.setZero(n, n);
  
  auto& engine = *engines[2];
  const auto& shells = libint_basis->shells();
  
  // Set nuclear charges and positions
  std::vector<std::pair<double, std::array<double, 3>>> charges_pos;
  for (size_t i = 0; i < charges.size(); ++i) {
    charges_pos.push_back({charges[i], 
                          {positions[i][0], positions[i][1], positions[i][2]}});
  }
  engine.set_params(charges_pos);
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int n1 = shells[s1].size();
    int bf2 = 0;
    
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      int n2 = shells[s2].size();
      
      engine.compute(shells[s1], shells[s2]);
      const auto* buf = engine.results()[0];
      
      if (buf != nullptr) {
        for (int f1 = 0; f1 < n1; ++f1) {
          for (int f2 = 0; f2 < n2; ++f2) {
            V(bf1 + f1, bf2 + f2) += buf[f1 * n2 + f2];
            if (s1 != s2) {
              V(bf2 + f2, bf1 + f1) += buf[f1 * n2 + f2];
            }
          }
        }
      }
      
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_eri(std::vector<double> &eri_tensor)
{
  int n = basis_set->get_n_basis();
  eri_tensor.resize(n * n * n * n);
  std::fill(eri_tensor.begin(), eri_tensor.end(), 0.0);
  
  auto& engine = *engines[3];
  const auto& shells = libint_basis->shells();
  
  // This is expensive - compute all unique shell quartets
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int n1 = shells[s1].size();
    int bf2 = 0;
    
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      int n2 = shells[s2].size();
      int bf3 = 0;
      
      for (size_t s3 = 0; s3 <= s1; ++s3) {
        int n3 = shells[s3].size();
        int bf4 = 0;
        
        for (size_t s4 = 0; s4 <= ((s1 == s3) ? s2 : s3); ++s4) {
          int n4 = shells[s4].size();
          
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf = engine.results()[0];
          
          if (buf != nullptr) {
            for (int f1 = 0; f1 < n1; ++f1) {
              for (int f2 = 0; f2 < n2; ++f2) {
                for (int f3 = 0; f3 < n3; ++f3) {
                  for (int f4 = 0; f4 < n4; ++f4) {
                    double val = buf[f1*n2*n3*n4 + f2*n3*n4 + f3*n4 + f4];
                    
                    int i = bf1 + f1;
                    int j = bf2 + f2;
                    int k = bf3 + f3;
                    int l = bf4 + f4;
                    
                    // Store with 8-fold permutation symmetry
                    eri_tensor[i*n*n*n + j*n*n + k*n + l] = val;
                    eri_tensor[j*n*n*n + i*n*n + k*n + l] = val;
                    eri_tensor[i*n*n*n + j*n*n + l*n + k] = val;
                    eri_tensor[j*n*n*n + i*n*n + l*n + k] = val;
                    eri_tensor[k*n*n*n + l*n*n + i*n + j] = val;
                    eri_tensor[l*n*n*n + k*n*n + i*n + j] = val;
                    eri_tensor[k*n*n*n + l*n*n + j*n + i] = val;
                    eri_tensor[l*n*n*n + k*n*n + j*n + i] = val;
                  }
                }
              }
            }
          }
          
          bf4 += n4;
        }
        bf3 += n3;
      }
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_eri_with_density(const Eigen::MatrixXd &D,
                                              Eigen::MatrixXd &J,
                                              Eigen::MatrixXd &K)
{
  int n = basis_set->get_n_basis();
  J.setZero(n, n);
  K.setZero(n, n);
  
  auto& engine = *engines[3];
  const auto& shells = libint_basis->shells();
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int n1 = shells[s1].size();
    int bf2 = 0;
    
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      int n2 = shells[s2].size();
      int bf3 = 0;
      
      for (size_t s3 = 0; s3 < shells.size(); ++s3) {
        int n3 = shells[s3].size();
        int bf4 = 0;
        
        for (size_t s4 = 0; s4 <= s3; ++s4) {
          int n4 = shells[s4].size();
          
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          const auto* buf = engine.results()[0];
          
          if (buf != nullptr) {
            for (int f1 = 0; f1 < n1; ++f1) {
              for (int f2 = 0; f2 < n2; ++f2) {
                for (int f3 = 0; f3 < n3; ++f3) {
                  for (int f4 = 0; f4 < n4; ++f4) {
                    double val = buf[f1*n2*n3*n4 + f2*n3*n4 + f3*n4 + f4];
                    
                    int i = bf1 + f1;
                    int j = bf2 + f2;
                    int k = bf3 + f3;
                    int l = bf4 + f4;
                    
                    // Coulomb contribution
                    J(i, j) += D(k, l) * val;
                    if (s3 != s4) J(i, j) += D(l, k) * val;
                    if (s1 != s2) {
                      J(j, i) += D(k, l) * val;
                      if (s3 != s4) J(j, i) += D(l, k) * val;
                    }
                    
                    // Exchange contribution
                    K(i, l) += D(j, k) * val;
                    K(i, k) += D(j, l) * val;
                    K(j, l) += D(i, k) * val;
                    K(j, k) += D(i, l) * val;
                  }
                }
              }
            }
          }
          
          bf4 += n4;
        }
        bf3 += n3;
      }
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_overlap_gradient(std::vector<Eigen::MatrixXd> &dS)
{
  // Gradient of overlap integrals
  // This would require derivative integrals from libint2
  // Placeholder implementation
  
  int n = basis_set->get_n_basis();
  int n_atoms = basis_set->get_shell_to_atom_map().size();
  
  dS.resize(3 * n_atoms);
  for (auto& mat : dS) {
    mat.setZero(n, n);
  }
  
  // TODO: Implement using libint2 derivative engines
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_kinetic_gradient(std::vector<Eigen::MatrixXd> &dT)
{
  // Gradient of kinetic integrals
  int n = basis_set->get_n_basis();
  int n_atoms = basis_set->get_shell_to_atom_map().size();
  
  dT.resize(3 * n_atoms);
  for (auto& mat : dT) {
    mat.setZero(n, n);
  }
  
  // TODO: Implement using libint2 derivative engines
}

/* ---------------------------------------------------------------------- */

void PairDFT::compute_nuclear_gradient(std::vector<Eigen::MatrixXd> &dV)
{
  // Gradient of nuclear attraction integrals
  int n = basis_set->get_n_basis();
  int n_atoms = basis_set->get_shell_to_atom_map().size();
  
  dV.resize(3 * n_atoms);
  for (auto& mat : dV) {
    mat.setZero(n, n);
  }
  
  // TODO: Implement using libint2 derivative engines
}
