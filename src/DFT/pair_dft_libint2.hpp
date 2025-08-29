/* ----------------------------------------------------------------------
   Libint2 integral methods for PairDFT
------------------------------------------------------------------------- */

#include <libint2/basis.h>
#include <libint2/shell.h>
#include <libint2/engine.h>

/* ----------------------------------------------------------------------
   Initialize libint2 engines
------------------------------------------------------------------------- */

void PairDFT::initialize_libint()
{
  // Create libint2 basis set from our basis data
  std::vector<libint2::Shell> shells;
  
  for (int s = 0; s < n_shells; s++) {
    int l = angular_momentum[s];
    auto shell_exponents = get_exponents(s);
    auto shell_coefficients = get_coefficients(s);
    
    // Convert to libint2's svector (small_vector)
    libint2::svector<double> alpha_svec;
    for (auto exp : shell_exponents) alpha_svec.push_back(exp);
    
    // Create contraction
    libint2::Shell::Contraction contr;
    contr.l = l;
    contr.pure = false;  // Use Cartesian Gaussians
    
    // Convert coefficients to svector
    for (auto coeff : shell_coefficients) contr.coeff.push_back(coeff);
    
    // Create svector of contractions
    libint2::svector<libint2::Shell::Contraction> contr_svec;
    contr_svec.push_back(contr);
    
    // Create shell
    libint2::Shell shell;
    shell.alpha = alpha_svec;
    shell.contr = contr_svec;
    shell.O = {{0.0, 0.0, 0.0}};  // Origin - placeholder, updated later with atom
    
    shells.push_back(shell);
  }
  
  libint_basis = std::make_unique<libint2::BasisSet>(shells);
  
  // Initialize engines for different integral types
  engines.resize(6); // Added space for gradient engines
  
  // Get max angular momentum and number of primitives
  int max_l = 0;
  size_t max_nprim = 0;
  for (const auto& shell : shells) {
    max_l = std::max(max_l, shell.contr[0].l);
    max_nprim = std::max(max_nprim, shell.alpha.size());
  }
  
  engines[0] = std::make_unique<libint2::Engine>(libint2::Operator::overlap, max_nprim, max_l);
  engines[1] = std::make_unique<libint2::Engine>(libint2::Operator::kinetic, max_nprim, max_l);
  engines[2] = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_l);
  engines[3] = std::make_unique<libint2::Engine>(libint2::Operator::coulomb, max_nprim, max_l);
  
  // Initialize derivative engines for gradient calculations
  engines[4] = std::make_unique<libint2::Engine>(libint2::Operator::overlap, max_nprim, max_l, 1); // overlap gradient
  engines[5] = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_l, 1); // nuclear gradient
  
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("Initialized libint2 with max_l={}, max_nprim={}\n", 
                                    max_l, max_nprim));
  }
}

/* ----------------------------------------------------------------------
   Cleanup libint2
------------------------------------------------------------------------- */

void PairDFT::cleanup_libint()
{
  engines.clear();
  libint_basis.reset();
}

/* ----------------------------------------------------------------------
   Compute one-electron integrals
------------------------------------------------------------------------- */

void PairDFT::compute_one_electron_integrals()
{
  compute_overlap_integrals();
  compute_kinetic_integrals();
  compute_nuclear_integrals();
}

/* ----------------------------------------------------------------------
   Compute overlap integrals
------------------------------------------------------------------------- */

void PairDFT::compute_overlap_integrals()
{
  overlap_matrix.setZero();
  
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
            overlap_matrix(bf1 + f1, bf2 + f2) = buf[f1 * n2 + f2];
            if (s1 != s2) overlap_matrix(bf2 + f2, bf1 + f1) = buf[f1 * n2 + f2];
          }
        }
      }
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ----------------------------------------------------------------------
   Compute kinetic energy integrals
------------------------------------------------------------------------- */

void PairDFT::compute_kinetic_integrals()
{
  kinetic_matrix.setZero();
  
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
            kinetic_matrix(bf1 + f1, bf2 + f2) = buf[f1 * n2 + f2];
            if (s1 != s2) kinetic_matrix(bf2 + f2, bf1 + f1) = buf[f1 * n2 + f2];
          }
        }
      }
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ----------------------------------------------------------------------
   Compute nuclear attraction integrals
------------------------------------------------------------------------- */

void PairDFT::compute_nuclear_integrals()
{
  nuclear_matrix.setZero();
  
  // Get nuclear charges and positions
  std::vector<std::pair<double, std::array<double, 3>>> charges_pos;
  double **x = atom->x;
  double *q = atom->q;  // Use atom->q for nuclear charges
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    charges_pos.push_back({q[i], {x[i][0], x[i][1], x[i][2]}});
  }
  
  auto& engine = *engines[2];
  engine.set_params(charges_pos);
  
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
            nuclear_matrix(bf1 + f1, bf2 + f2) += buf[f1 * n2 + f2];
            if (s1 != s2) {
              nuclear_matrix(bf2 + f2, bf1 + f1) += buf[f1 * n2 + f2];
            }
          }
        }
      }
      
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ----------------------------------------------------------------------
   Compute electron repulsion integrals with density
------------------------------------------------------------------------- */

void PairDFT::compute_eri_with_density(const Eigen::MatrixXd &D,
                                       Eigen::MatrixXd &J,
                                       Eigen::MatrixXd &K)
{
  J.setZero();
  K.setZero();
  
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
                    double eri = buf[f1*n2*n3*n4 + f2*n3*n4 + f3*n4 + f4];
                    
                    int mu = bf1 + f1;
                    int nu = bf2 + f2;
                    int lambda = bf3 + f3;
                    int sigma = bf4 + f4;
                    
                    // Coulomb matrix: J[μν] = sum_λσ D[λσ] * (μν|λσ)
                    double coulomb_contrib = D(lambda, sigma) * eri;
                    J(mu, nu) += coulomb_contrib;
                    if (s1 != s2) J(nu, mu) += coulomb_contrib;
                    if (s3 != s4) {
                      double coulomb_contrib2 = D(sigma, lambda) * eri;
                      J(mu, nu) += coulomb_contrib2;
                      if (s1 != s2) J(nu, mu) += coulomb_contrib2;
                    }
                    
                    // Exchange matrix: K[μν] = sum_λσ D[μλ] * (μν|λσ)
                    // We need to handle the 8-fold permutation symmetry properly
                    
                    // (μν|λσ) = (νμ|λσ) = (μν|σλ) = (νμ|σλ) = (λσ|μν) = (σλ|μν) = (λσ|νμ) = (σλ|νμ)
                    
                    // K[μσ] -= 0.5 * D[νλ] * (μν|λσ)
                    K(mu, sigma) -= 0.5 * D(nu, lambda) * eri;
                    K(sigma, mu) -= 0.5 * D(nu, lambda) * eri;
                    
                    if (s1 != s2) {
                      K(nu, sigma) -= 0.5 * D(mu, lambda) * eri;
                      K(sigma, nu) -= 0.5 * D(mu, lambda) * eri;
                    }
                    
                    if (s3 != s4) {
                      K(mu, lambda) -= 0.5 * D(nu, sigma) * eri;
                      K(lambda, mu) -= 0.5 * D(nu, sigma) * eri;
                      
                      if (s1 != s2) {
                        K(nu, lambda) -= 0.5 * D(mu, sigma) * eri;
                        K(lambda, nu) -= 0.5 * D(mu, sigma) * eri;
                      }
                    }
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

/* ----------------------------------------------------------------------
   Compute overlap gradient for force calculation
------------------------------------------------------------------------- */

void PairDFT::compute_overlap_gradient(std::vector<Eigen::MatrixXd> &dS)
{
  int nlocal = atom->nlocal;
  
  dS.resize(3 * nlocal);
  for (auto& mat : dS) {
    mat.setZero(n_basis_functions, n_basis_functions);
  }
  
  auto& engine = *engines[4]; // overlap gradient engine
  const auto& shells = libint_basis->shells();
  
  // Set up derivative calculations
  engine.set(libint2::BraKet::xx_xx); // derivatives w.r.t. bra and ket centers
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int n1 = shells[s1].size();
    int bf2 = 0;
    int atom1 = shell_to_atom[s1];
    
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      int n2 = shells[s2].size();
      int atom2 = shell_to_atom[s2];
      
      engine.compute(shells[s1], shells[s2]);
      
      // Results are organized as [x1, y1, z1, x2, y2, z2] for derivatives
      // w.r.t. center 1 (x,y,z) and center 2 (x,y,z)
      for (int deriv = 0; deriv < 6; ++deriv) {
        const auto* buf = engine.results()[deriv];
        if (buf == nullptr) continue;
        
        int atom_idx = (deriv < 3) ? atom1 : atom2;
        int cart_comp = deriv % 3;
        int mat_idx = 3 * atom_idx + cart_comp;
        
        if (mat_idx >= 3 * nlocal) continue;
        
        for (int f1 = 0; f1 < n1; ++f1) {
          for (int f2 = 0; f2 < n2; ++f2) {
            double value = buf[f1 * n2 + f2];
            dS[mat_idx](bf1 + f1, bf2 + f2) += value;
            
            if (s1 != s2) {
              dS[mat_idx](bf2 + f2, bf1 + f1) += value;
            }
          }
        }
      }
      
      bf2 += n2;
    }
    bf1 += n1;
  }
}

/* ----------------------------------------------------------------------
   Compute nuclear gradient for force calculation
------------------------------------------------------------------------- */

void PairDFT::compute_nuclear_gradient(std::vector<Eigen::MatrixXd> &dV)
{
  int nlocal = atom->nlocal;
  
  dV.resize(3 * nlocal);
  for (auto& mat : dV) {
    mat.setZero(n_basis_functions, n_basis_functions);
  }
  
  // Get nuclear charges and positions
  std::vector<std::pair<double, std::array<double, 3>>> charges_pos;
  double **x = atom->x;
  double *q = atom->q;
  
  for (int i = 0; i < nlocal; i++) {
    charges_pos.push_back({q[i], {x[i][0], x[i][1], x[i][2]}});
  }
  
  auto& engine = *engines[5]; // nuclear gradient engine
  engine.set_params(charges_pos);
  engine.set(libint2::BraKet::xx_xx); // derivatives w.r.t. bra and ket centers
  
  const auto& shells = libint_basis->shells();
  
  int bf1 = 0;
  for (size_t s1 = 0; s1 < shells.size(); ++s1) {
    int n1 = shells[s1].size();
    int bf2 = 0;
    int atom1 = shell_to_atom[s1];
    
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      int n2 = shells[s2].size();
      int atom2 = shell_to_atom[s2];
      
      engine.compute(shells[s1], shells[s2]);
      
      // Process derivatives w.r.t. basis function centers
      for (int deriv = 0; deriv < 6; ++deriv) {
        const auto* buf = engine.results()[deriv];
        if (buf == nullptr) continue;
        
        int atom_idx = (deriv < 3) ? atom1 : atom2;
        int cart_comp = deriv % 3;
        int mat_idx = 3 * atom_idx + cart_comp;
        
        if (mat_idx >= 3 * nlocal) continue;
        
        for (int f1 = 0; f1 < n1; ++f1) {
          for (int f2 = 0; f2 < n2; ++f2) {
            double value = buf[f1 * n2 + f2];
            dV[mat_idx](bf1 + f1, bf2 + f2) += value;
            
            if (s1 != s2) {
              dV[mat_idx](bf2 + f2, bf1 + f1) += value;
            }
          }
        }
      }
      
      // Also need to handle derivatives w.r.t. nuclear positions
      // This would involve additional derivative engines and more complex logic
      // For now, we handle only the basis function center derivatives
      
      bf2 += n2;
    }
    bf1 += n1;
  }
}

