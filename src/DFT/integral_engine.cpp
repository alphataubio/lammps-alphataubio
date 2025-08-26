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
#include <libint2.hpp>
#include <cmath>
#include <vector>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

IntegralEngine::IntegralEngine(BasisSetManager *basis) : basis_set(basis)
{
  initialize_libint();
}

/* ---------------------------------------------------------------------- */

IntegralEngine::~IntegralEngine()
{
  cleanup_libint();
}

/* ---------------------------------------------------------------------- */

void IntegralEngine::initialize_libint()
{
  // Create libint2 basis set from our basis manager
  std::vector<libint2::Shell> shells;
  
  int n_shells = basis_set->get_n_shells();
  auto shell_to_atom = basis_set->get_shell_to_atom_map();
  
  for (int s = 0; s < n_shells; s++) {
    int l = basis_set->get_angular_momentum(s);
    auto exponents = basis_set->get_exponents(s);
    auto coefficients = basis_set->get_coefficients(s);
    
    // Create contracted Gaussian shell
    libint2::Shell shell;
    shell.alpha = exponents;
    shell.contr = {coefficients};
    shell.l = {l};
    
    // Set origin (would come from atom positions in real use)
    shell.O = {0.0, 0.0, 0.0};  // Placeholder
    
    shells.push_back(shell);
  }
  
  libint_basis = std::make_unique<libint2::BasisSet>(shells);
  
  // Initialize engines for different integral types
  engines.resize(4);
  engines[0] = std::make_unique<libint2::Engine>(libint2::Operator::overlap, 
                                                 libint_basis->max_nprim(), 
                                                 libint_basis->max_l());
  engines[1] = std::make_unique<libint2::Engine>(libint2::Operator::kinetic,
                                                 libint_basis->max_nprim(),
                                                 libint_basis->max_l());
  engines[2] = std::make_unique<libint2::Engine>(libint2::Operator::nuclear,
                                                 libint_basis->max_nprim(),
                                                 libint_basis->max_l());
  engines[3] = std::make_unique<libint2::Engine>(libint2::Operator::coulomb,
                                                 libint_basis->max_nprim(),
                                                 libint_basis->max_l());
}

/* ---------------------------------------------------------------------- */

void IntegralEngine::cleanup_libint()
{
  engines.clear();
  libint_basis.reset();
}

/* ---------------------------------------------------------------------- */

void IntegralEngine::compute_overlap(Eigen::MatrixXd &S)
{
  int n = basis_set->get_n_basis();
  S.setZero(n, n);
  
  auto& engine = *engines[0];
  const auto& shells = libint_basis->shells();
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int bf2 = 0;
    int n1 = shells[s1].size();
    
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

void IntegralEngine::compute_kinetic(Eigen::MatrixXd &T)
{
  int n = basis_set->get_n_basis();
  T.setZero(n, n);
  
  auto& engine = *engines[1];
  const auto& shells = libint_basis->shells();
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int bf2 = 0;
    int n1 = shells[s1].size();
    
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

void IntegralEngine::compute_nuclear(Eigen::MatrixXd &V, 
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
    int bf2 = 0;
    int n1 = shells[s1].size();
    
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

void IntegralEngine::compute_eri(std::vector<double> &eri_tensor)
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

void IntegralEngine::compute_eri_with_density(const Eigen::MatrixXd &D,
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

void IntegralEngine::compute_overlap_gradient(std::vector<Eigen::MatrixXd> &dS)
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

void IntegralEngine::compute_kinetic_gradient(std::vector<Eigen::MatrixXd> &dT)
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

void IntegralEngine::compute_nuclear_gradient(std::vector<Eigen::MatrixXd> &dV)
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
