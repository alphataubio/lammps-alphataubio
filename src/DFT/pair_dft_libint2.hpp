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
  // Initialize libint2 library
  libint2::initialize();
  
  // Create libint2 basis set from our basis data
  std::vector<libint2::Shell> shells;
  
  // Get atom positions for shell origins
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  for (int s = 0; s < n_shells; s++) {
    int l = angular_momentum[s];
    auto shell_exponents = exponents[s];
    auto shell_coefficients = coefficients[s];
    
    if (comm->me == 0) {
      utils::logmesg(lmp, fmt::format("Shell {}: l={}, n_exp={}, n_coeff={}\n", 
                                      s, l, shell_exponents.size(), shell_coefficients.size()));
    }
    
    // Skip shells with no exponents
    if (shell_exponents.empty()) {
      if (comm->me == 0) {
        utils::logmesg(lmp, fmt::format("Warning: Shell {} has no exponents, skipping\n", s));
      }
      continue;
    }
    
    // Convert to libint2's svector (small_vector)
    libint2::svector<double> alpha_svec;
    for (auto exp : shell_exponents) alpha_svec.push_back(exp);
    
    // Create contraction - coefficients must match exponents size
    libint2::Shell::Contraction contr;
    contr.l = l;
    contr.pure = false;  // Use Cartesian Gaussians
    
    // Make sure we have the right number of coefficients
    if (shell_coefficients.size() != shell_exponents.size()) {
      error->all(FLERR, fmt::format("Shell {} has {} exponents but {} coefficients", 
                                    s, shell_exponents.size(), shell_coefficients.size()));
    }
    
    // Convert coefficients to svector
    // Libint2 uses normalized primitives by default
    // Check if BSE coefficients need renormalization
    std::vector<double> normalized_coeff;
    
    // First compute the self-overlap of the contraction
    double S_contr = 0.0;
    for (size_t i = 0; i < shell_coefficients.size(); i++) {
      for (size_t j = 0; j < shell_coefficients.size(); j++) {
        double alpha_i = shell_exponents[i];
        double alpha_j = shell_exponents[j];
        double ci = shell_coefficients[i];
        double cj = shell_coefficients[j];
        
        // Overlap integral for primitives with same center
        // S_ij = (pi/(alpha_i + alpha_j))^(3/2)
        double S_ij = pow(M_PI / (alpha_i + alpha_j), 1.5);
        
        // Add angular momentum factor for l > 0
        double ang_factor = 1.0;
        for (int k = 0; k < l; k++) {
          ang_factor *= (2*k + 1) / (2.0 * (alpha_i + alpha_j));
        }
        
        S_contr += ci * cj * S_ij * ang_factor;
      }
    }
    
    // Normalize the contraction if needed
    double norm_factor = 1.0 / sqrt(S_contr);
    
    for (auto coeff : shell_coefficients) {
      contr.coeff.push_back(coeff * norm_factor);
    }
    
    // Create svector of contractions
    libint2::svector<libint2::Shell::Contraction> contr_svec;
    contr_svec.push_back(contr);
    
    // Create shell using the constructor
    // Set origin based on atom position
    int atom_idx = (s < shell_to_atom.size()) ? shell_to_atom[s] : 0;
    std::array<double, 3> origin = {{0.0, 0.0, 0.0}};
    if (atom_idx < nlocal && nlocal > 0) {
      origin = {{x[atom_idx][0] * ANGSTROM_TO_BOHR, 
                 x[atom_idx][1] * ANGSTROM_TO_BOHR, 
                 x[atom_idx][2] * ANGSTROM_TO_BOHR}};
    }
    
    // Use the constructor instead of setting fields directly
    libint2::Shell shell(alpha_svec, contr_svec, origin);
    
    shells.push_back(shell);
  }
  
  // Validate shells before creating basis set
  if (shells.empty()) {
    error->all(FLERR, "No shells created from basis set data");
  }
  
  for (size_t i = 0; i < shells.size(); ++i) {
    if (shells[i].alpha.empty()) {
      error->all(FLERR, fmt::format("Shell {} has no exponents", i));
    }
    if (shells[i].contr.empty()) {
      error->all(FLERR, fmt::format("Shell {} has no contractions", i));
    }
    for (const auto& contr : shells[i].contr) {
      if (contr.coeff.empty()) {
        error->all(FLERR, fmt::format("Shell {} has contraction with no coefficients", i));
      }
    }
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
  
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("Creating engines with max_l={}, max_nprim={}\n", 
                                    max_l, max_nprim));
  }
  
  // Initialize standard integral engines
  try {
    engines[0] = std::make_unique<libint2::Engine>(libint2::Operator::overlap, max_nprim, max_l);
    if (comm->me == 0) utils::logmesg(lmp, "Overlap engine created\n");
    
    engines[1] = std::make_unique<libint2::Engine>(libint2::Operator::kinetic, max_nprim, max_l);
    if (comm->me == 0) utils::logmesg(lmp, "Kinetic engine created\n");
    
    engines[2] = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_l);
    if (comm->me == 0) utils::logmesg(lmp, "Nuclear engine created\n");
    
    engines[3] = std::make_unique<libint2::Engine>(libint2::Operator::coulomb, max_nprim, max_l);
    if (comm->me == 0) utils::logmesg(lmp, "Coulomb engine created\n");
  } catch (const std::exception& e) {
    error->all(FLERR, fmt::format("Failed to create standard integral engines: {}", e.what()));
  }
  
  // Initialize derivative engines for gradient calculations
  // According to libint2 documentation, the 4th parameter is deriv_order
  // Default max_l and max_nprim should work, but the issue might be:
  // 1. Shells not properly constructed
  // 2. max_nprim or max_l being 0 or invalid
  
  if (max_l < 0 || max_nprim < 1) {
    engines[4] = nullptr;
    engines[5] = nullptr;
    if (comm->me == 0) {
      utils::logmesg(lmp, fmt::format("Warning: Invalid max_l={} or max_nprim={}, gradient engines disabled\n", 
                                      max_l, max_nprim));
    }
  } else {
    try {
      // Create derivative engines with deriv_order=1
      //engines[4] = std::make_unique<libint2::Engine>(libint2::Operator::overlap, max_nprim, max_l, 1);
      //engines[5] = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_l, 1);
      
      if (comm->me == 0) {
        utils::logmesg(lmp, "FIXME Gradient engines initialized successfully\n");
      }
    } catch (const std::exception& e) {
      engines[4] = nullptr;
      engines[5] = nullptr;
      if (comm->me == 0) {
        utils::logmesg(lmp, fmt::format("Failed to initialize gradient engines: {}\n", e.what()));
        utils::logmesg(lmp, "Gradient calculations will be disabled\n");
      }
    }
  }
  
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
  libint2::finalize();
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
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    // For hydrogen atoms, nuclear charge Z = 1
    // TODO: Get actual atomic numbers from atom types
    double nuclear_charge = 1.0;
    charges_pos.push_back({nuclear_charge, {x[i][0] * ANGSTROM_TO_BOHR, x[i][1] * ANGSTROM_TO_BOHR, x[i][2] * ANGSTROM_TO_BOHR}});
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
                    J(mu, nu) += D(lambda, sigma) * eri;
                    if (s1 != s2) J(nu, mu) += D(lambda, sigma) * eri;
                    
                    // Exchange matrix: K[μν] = sum_λσ D[μσ] * (μλ|νσ)
                    // Using permutation symmetry: (μν|λσ) = (λσ|μν)
                    K(mu, lambda) += D(nu, sigma) * eri;
                    if (mu != lambda) K(lambda, mu) += D(sigma, nu) * eri;
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
  
  // Check if gradient engine is available
  if (!engines[4]) {
    if (comm->me == 0) {
      error->warning(FLERR, "Overlap gradient engine not available - gradients will be zero");
    }
    return;
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
  
  // Check if gradient engine is available
  if (!engines[5]) {
    if (comm->me == 0) {
      error->warning(FLERR, "Nuclear gradient engine not available - gradients will be zero");
    }
    return;
  }
  
  // Get nuclear charges and positions
  std::vector<std::pair<double, std::array<double, 3>>> charges_pos;
  double **x = atom->x;
  
  for (int i = 0; i < nlocal; i++) {
    // For hydrogen atoms, nuclear charge Z = 1
    // TODO: Get actual atomic numbers from atom types
    double nuclear_charge = 1.0;
    charges_pos.push_back({nuclear_charge, {x[i][0] * ANGSTROM_TO_BOHR, x[i][1] * ANGSTROM_TO_BOHR, x[i][2] * ANGSTROM_TO_BOHR}});
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

