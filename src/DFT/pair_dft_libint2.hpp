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
            if (s1 != s2) {
              overlap_matrix(bf2 + f2, bf1 + f1) = buf[f1 * n2 + f2];
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
            if (s1 != s2) {
              kinetic_matrix(bf2 + f2, bf1 + f1) = buf[f1 * n2 + f2];
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
                    
                    // FIXME Exchange contribution (simplified - needs proper symmetry handling)
                    K(i, l) += 0.5 * D(j, k) * val;
                    K(i, k) += 0.5 * D(j, l) * val;
                    K(j, l) += 0.5 * D(i, k) * val;
                    K(j, k) += 0.5 * D(i, l) * val;
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
  // Gradient of overlap integrals
  int nlocal = atom->nlocal;
  
  dS.resize(3 * nlocal);
  for (auto& mat : dS) {
    mat.setZero(n_basis_functions, n_basis_functions);
  }
  
  // FIXME Implement using libint2 derivative engines
  // This requires setting up gradient engines and computing
  // derivative integrals properly
}

/* ----------------------------------------------------------------------
   Compute nuclear gradient for force calculation
------------------------------------------------------------------------- */

void PairDFT::compute_nuclear_gradient(std::vector<Eigen::MatrixXd> &dV)
{
  // Gradient of nuclear attraction integrals
  int nlocal = atom->nlocal;
  
  dV.resize(3 * nlocal);
  for (auto& mat : dV) {
    mat.setZero(n_basis_functions, n_basis_functions);
  }
  
  // FIXME Implement using libint2 derivative engines
  // This requires setting up gradient engines and computing
  // derivative integrals properly
}
