/* ----------------------------------------------------------------------
   SCF and grid integration methods for PairDFT
------------------------------------------------------------------------- */

#include <cmath>
#include <vector>
#include <algorithm>

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
    update_density_matrix();
    
    // Mix with previous density for better convergence
    if (current_iteration > 1) {
      mix_density_matrices(0.5);
    }
    
    // Compute energy
    compute_energy();
    
    // Check convergence
    double energy_change = std::abs(total_dft_energy - prev_energy);
    double density_change = compute_density_change();
    
    print_scf_iteration();
    
    if (energy_change < energy_tolerance && density_change < density_tolerance) {
      scf_converged = true;
      break;
    }
    
    prev_energy = total_dft_energy;
  }
  
  if (!scf_converged && comm->me == 0) {
    error->warning(FLERR, "SCF did not converge within maximum iterations");
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
  compute_eri_with_density(density_matrix, coulomb_matrix, exchange_matrix);
  
  // Add Coulomb contribution (factor of 2 for closed-shell)
  fock_matrix += 2.0 * coulomb_matrix;
  
  // Add exchange contribution (factor of -1 for exchange, no factor of 2)
  // Note: exchange_matrix already includes appropriate factors from ERI computation
  // For pure DFT (non-hybrid), exchange is handled via XC functional
  // For hybrid functionals, we include a fraction of exact exchange
  if (is_hybrid) {
    fock_matrix -= hybrid_coeff * exchange_matrix;
  }
  
  // Add XC contribution via grid integration
  Eigen::MatrixXd vxc_matrix(n_basis_functions, n_basis_functions);
  vxc_matrix.setZero();
  
  // Generate grid and integrate XC
  generate_integration_grid();
  integrate_xc_on_grid(xc_energy, vxc_matrix);
  
  fock_matrix += vxc_matrix;
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

/* ----------------------------------------------------------------------
   Update density matrix from MO coefficients
------------------------------------------------------------------------- */

void PairDFT::update_density_matrix()
{
  // Save previous density
  density_matrix_prev = density_matrix;
  
  // Build new density matrix from occupied orbitals
  // D = 2 * C_occ * C_occ^T for closed shell
  density_matrix.setZero();
  
  for (int i = 0; i < n_occupied_orbitals; i++) {
    density_matrix += 2.0 * mo_coefficients.col(i) * mo_coefficients.col(i).transpose();
  }
}

/* ----------------------------------------------------------------------
   Mix density matrices for convergence
------------------------------------------------------------------------- */

void PairDFT::mix_density_matrices(double mixing_param)
{
  // Simple linear mixing for stability
  // D_new = (1-alpha)*D_old + alpha*D_current
  density_matrix = mixing_param * density_matrix + (1.0 - mixing_param) * density_matrix_prev;
}

/* ----------------------------------------------------------------------
   Compute RMS density change
------------------------------------------------------------------------- */

double PairDFT::compute_density_change()
{
  Eigen::MatrixXd diff = density_matrix - density_matrix_prev;
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < n_basis_functions; i++) {
    for (int j = 0; j < n_basis_functions; j++) {
      sum += diff(i, j) * diff(i, j);
      count++;
    }
  }
  
  return sqrt(sum / count);
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
  double one_electron = (density_matrix.cwiseProduct(H_core)).sum();
  
  // Two-electron energy (J and K contributions)
  // For closed-shell: E = Tr(P*H) + 0.5*Tr(P*G)
  // where G = 2*J - K for RHF or 2*J - alpha*K for hybrid DFT
  // and P is the total density matrix (factor of 2 already included)
  double j_energy = (density_matrix.cwiseProduct(coulomb_matrix)).sum();
  double k_energy = 0.0;
  
  if (is_hybrid) {
    k_energy = hybrid_coeff * 0.5 * (density_matrix.cwiseProduct(exchange_matrix)).sum();
  }
  
  // Kinetic energy component (for output)
  kinetic_energy = (density_matrix.cwiseProduct(kinetic_matrix)).sum();
  
  // Total energy: E = E_nuc + Tr(P*H_core) + 0.5*Tr(P*(2J-K)) + E_xc
  // For closed-shell: Tr(P*H_core) already includes factor of 2 from density
  // The 0.5 factor accounts for double-counting in two-electron terms
  total_dft_energy = nuclear_repulsion + one_electron + 0.5 * j_energy - k_energy + xc_energy;
}

/* ----------------------------------------------------------------------
   Compute nuclear repulsion energy
------------------------------------------------------------------------- */

double PairDFT::compute_nuclear_repulsion_energy()
{
  double energy = 0.0;
  
  double **x = atom->x;
  // Note: We don't use atom->q anymore, nuclear charges are determined from atom types
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    for (int j = i + 1; j < nlocal; j++) {
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      double r = sqrt(dx*dx + dy*dy + dz*dz) * ANGSTROM_TO_BOHR;
      
      // Use atomic number Z for nuclear charges (for H2, Z=1)
      // TODO: Get actual atomic numbers from atom types
      const double Zi = 1.0, Zj = 1.0;  // Hydrogen atoms

      if (r > 1e-10) energy += Zi * Zj / r;
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
  int nlocal = atom->nlocal;
  
  for (int i = 0; i < nlocal; i++) {
    // For H2 molecule, each H atom contributes 1 electron
    // TODO: Get actual atomic numbers from atom types
    int nuclear_charge = 1;  // Hydrogen
    n_electrons += nuclear_charge;
  }
  
  n_occupied_orbitals = n_electrons / 2;  // Assuming closed-shell
  
  update_density_matrix();
}

/* ----------------------------------------------------------------------
   Generate integration grid
------------------------------------------------------------------------- */

void PairDFT::generate_integration_grid()
{
  grid_points.clear();
  grid_weights.clear();
  grid_atom_owners.clear();
  
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  if (nlocal == 0) return;
  
  // Generate atomic grids
  int points_per_atom = grid_size / nlocal;
  
  for (int atom_idx = 0; atom_idx < nlocal; atom_idx++) {
    // Get radial and angular grid sizes
    int n_radial = 50;  // Typical value
    int n_angular = points_per_atom / n_radial;
    
    // Generate radial grid
    std::vector<double> r_points(n_radial);
    std::vector<double> r_weights(n_radial);
    // Use Z=1 for hydrogen atoms
    double nuclear_charge = 1.0;
    generate_radial_grid(n_radial, nuclear_charge, r_points, r_weights);
    
    // Generate angular grid
    std::vector<std::vector<double>> angular_points;
    std::vector<double> angular_weights;
    generate_lebedev_grid(n_angular, angular_points, angular_weights);
    
    // Combine radial and angular grids
    for (int i_r = 0; i_r < n_radial; i_r++) {
      double r = r_points[i_r];
      double w_r = r_weights[i_r];
      
      for (size_t i_ang = 0; i_ang < angular_points.size(); i_ang++) {
        std::vector<double> point(3);
        point[0] = x[atom_idx][0] * ANGSTROM_TO_BOHR + r * angular_points[i_ang][0];
        point[1] = x[atom_idx][1] * ANGSTROM_TO_BOHR + r * angular_points[i_ang][1];
        point[2] = x[atom_idx][2] * ANGSTROM_TO_BOHR + r * angular_points[i_ang][2];
        
        grid_points.push_back(point);
        grid_weights.push_back(w_r * angular_weights[i_ang] * r * r);
        grid_atom_owners.push_back(atom_idx);
      }
    }
  }
  
  // Compute Becke partitioning weights
  std::vector<std::vector<double>> atom_positions;
  for (int i = 0; i < nlocal; i++) {
    atom_positions.push_back({x[i][0] * ANGSTROM_TO_BOHR, x[i][1] * ANGSTROM_TO_BOHR, x[i][2] * ANGSTROM_TO_BOHR});
  }
  compute_becke_weights(atom_positions);
}

/* ----------------------------------------------------------------------
   Generate Lebedev angular grid
------------------------------------------------------------------------- */

void PairDFT::generate_lebedev_grid(int n_points,
                                    std::vector<std::vector<double>> &points,
                                    std::vector<double> &weights)
{
  points.clear();
  weights.clear();
  
  // Find closest available Lebedev grid
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

/* ----------------------------------------------------------------------
   Generate radial grid
------------------------------------------------------------------------- */

void PairDFT::generate_radial_grid(int n_points, double Z,
                                   std::vector<double> &r,
                                   std::vector<double> &w)
{
  r.resize(n_points);
  w.resize(n_points);
  
  // Bragg radius for scaling
  double R_bragg = 1.0;
  if (Z > 0) R_bragg = 0.5 * (3.0 - 0.01 * Z);
  
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

/* ----------------------------------------------------------------------
   Compute Becke partitioning weights
------------------------------------------------------------------------- */

void PairDFT::compute_becke_weights(const std::vector<std::vector<double>> &atoms)
{
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
        for (int iter = 0; iter < 3; iter++) f = 0.5 * f * (3.0 - f * f);
        double s = 0.5 * (1.0 - f);
        
        P[i_atom] *= s;
      }
    }
    
    // Normalize partition functions to get actual Becke weights
    double sum = 0.0;
    for (int i_atom = 0; i_atom < n_atoms; i_atom++) sum += P[i_atom];
    
    // Get the weight for the atom that owns this grid point
    if (sum > 1e-15) {
      int owner_atom = grid_atom_owners[i_point];
      becke_weights[i_point] = P[owner_atom] / sum;
    } else {
      becke_weights[i_point] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   Integrate XC on grid
------------------------------------------------------------------------- */

void PairDFT::integrate_xc_on_grid(double &exc_energy, Eigen::MatrixXd &vxc_matrix)
{
  int n_points = grid_points.size();
  
  vxc_matrix.setZero();
  exc_energy = 0.0;
  
  if (n_points == 0) return;
  
  // Evaluate basis functions at grid points
  std::vector<std::vector<double>> basis_values;
  std::vector<std::vector<std::vector<double>>> basis_gradients;
  
  // For GGA and meta-GGA functionals, we need gradients
  // Initialize basis_gradients to request gradient computation
  // PBE is a GGA functional, so we always need gradients
  basis_gradients.resize(n_points);
  
  evaluate_basis_at_points(grid_points, basis_values, basis_gradients);
  
  // Compute density and gradients at grid points
  std::vector<double> rho(n_points, 0.0);
  std::vector<double> sigma(n_points, 0.0);  // |grad rho|^2
  std::vector<double> lapl(n_points, 0.0);   // Laplacian
  std::vector<double> tau(n_points, 0.0);     // Kinetic energy density
  
  for (int i_point = 0; i_point < n_points; i_point++) {
    // Density
    for (int i = 0; i < n_basis_functions; i++) {
      for (int j = 0; j < n_basis_functions; j++) {
        rho[i_point] += density_matrix(i, j) * basis_values[i_point][i] * basis_values[i_point][j];
      }
    }
    
    // Gradient for GGA/meta-GGA
    if (!basis_gradients.empty()) {
      std::vector<double> grad_rho(3, 0.0);
      
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j < n_basis_functions; j++) {
          double P_ij = density_matrix(i, j);
          
          for (int k = 0; k < 3; k++) {
            grad_rho[k] += P_ij * (basis_gradients[i_point][i][k] * basis_values[i_point][j] +
                                   basis_values[i_point][i] * basis_gradients[i_point][j][k]);
          }
          
          if (is_meta_gga) {
            for (int k = 0; k < 3; k++) {
              tau[i_point] += 0.5 * P_ij * basis_gradients[i_point][i][k] * basis_gradients[i_point][j][k];
            }
          }
        }
      }
      
      sigma[i_point] = grad_rho[0]*grad_rho[0] + grad_rho[1]*grad_rho[1] + grad_rho[2]*grad_rho[2];
    }
  }
  
  // Evaluate XC functional
  std::vector<double> exc(n_points);
  std::vector<double> vrho(n_points);
  std::vector<double> vsigma(n_points);
  std::vector<double> vlapl(n_points);
  std::vector<double> vtau(n_points);
  
  evaluate_xc_functional(rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);
  
  // Integrate XC energy
  for (int i_point = 0; i_point < n_points; i_point++) {
    exc_energy += exc[i_point] * rho[i_point] * grid_weights[i_point] * becke_weights[i_point];
  }
  
  // Build XC potential matrix
  for (int i_point = 0; i_point < n_points; i_point++) {
    double w = grid_weights[i_point] * becke_weights[i_point];
    
    // LDA contribution: v_xc[rho] * phi_i * phi_j
    for (int i = 0; i < n_basis_functions; i++) {
      for (int j = 0; j <= i; j++) {
        double val = vrho[i_point] * basis_values[i_point][i] * basis_values[i_point][j] * w;
        vxc_matrix(i, j) += val;
        if (i != j) vxc_matrix(j, i) += val;
      }
    }
    
    // GGA contribution: 2 * v_xc[sigma] * grad(rho) . grad(phi)
    if (!basis_gradients.empty() && vsigma[i_point] != 0.0) {
      // First compute gradient of density at this point
      std::vector<double> grad_rho(3, 0.0);
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j < n_basis_functions; j++) {
          double P_ij = density_matrix(i, j);
          for (int k = 0; k < 3; k++) {
            grad_rho[k] += P_ij * (basis_gradients[i_point][i][k] * basis_values[i_point][j] +
                                   basis_values[i_point][i] * basis_gradients[i_point][j][k]);
          }
        }
      }
      
      // Add GGA contribution to matrix elements
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j <= i; j++) {
          double gga_contrib = 0.0;
          
          // Contribution from derivative of functional w.r.t. |grad rho|^2
          for (int k = 0; k < 3; k++) {
            gga_contrib += 2.0 * vsigma[i_point] * w * (
              grad_rho[k] * (basis_gradients[i_point][i][k] * basis_values[i_point][j] +
                            basis_values[i_point][i] * basis_gradients[i_point][j][k])
            );
          }
          
          vxc_matrix(i, j) += gga_contrib;
          if (i != j) vxc_matrix(j, i) += gga_contrib;
        }
      }
    }
    
    // Meta-GGA contribution: v_xc[tau] * grad(phi_i) . grad(phi_j)
    if (is_meta_gga && vtau[i_point] != 0.0 && !basis_gradients.empty()) {
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j <= i; j++) {
          double tau_contrib = 0.0;
          
          // Kinetic energy density contribution
          for (int k = 0; k < 3; k++) {
            tau_contrib += 0.5 * vtau[i_point] * w * 
                          basis_gradients[i_point][i][k] * basis_gradients[i_point][j][k];
          }
          
          vxc_matrix(i, j) += tau_contrib;
          if (i != j) vxc_matrix(j, i) += tau_contrib;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Compute Hellmann-Feynman forces
------------------------------------------------------------------------- */

void PairDFT::compute_hellmann_feynman_forces()
{
  double **f = atom->f;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  // Nuclear-nuclear repulsion gradient
  for (int i = 0; i < nlocal; i++) {
    for (int j = 0; j < nlocal; j++) {
      if (i != j) {
        double dx = x[i][0] - x[j][0];
        double dy = x[i][1] - x[j][1];
        double dz = x[i][2] - x[j][2];
        double r = sqrt(dx*dx + dy*dy + dz*dz) * ANGSTROM_TO_BOHR;
        
        // Nuclear charges (Z=1 for hydrogen)
        double Zi = 1.0, Zj = 1.0;
        
        if (r > 1e-10) {
          double force_mag = Zi * Zj / (r * r * r);
          f[i][0] += force_mag * dx;
          f[i][1] += force_mag * dy;
          f[i][2] += force_mag * dz;
        }
      }
    }
  }
  
  // Electronic contribution to forces
  std::vector<Eigen::MatrixXd> dV(3 * nlocal);
  compute_nuclear_gradient(dV);
  
  for (int i = 0; i < nlocal; i++) {
    for (int k = 0; k < 3; k++) {
      double force_component = (density_matrix.cwiseProduct(dV[3*i + k])).sum();
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
  compute_overlap_gradient(dS);
  
  // Pulay force contribution
  for (int i = 0; i < nlocal; i++) {
    for (int k = 0; k < 3; k++) {
      double force_component = -(W.cwiseProduct(dS[3*i + k])).sum();
      f[i][k] += force_component;
    }
  }
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
    utils::logmesg(lmp, " Iter    Energy (Eh)    Delta E      RMS Dens\n");
    utils::logmesg(lmp, "------------------------------------------------\n");
  }
}

void PairDFT::print_scf_iteration()
{
  if (comm->me == 0) {
    double density_change = compute_density_change();
    // Note: energy_change is already computed in perform_scf() from prev_energy
    // This local calculation was incorrect
    static double prev_print_energy = 0.0;
    double energy_change = (current_iteration == 1) ? 0.0 : 
                          total_dft_energy - prev_print_energy;
    prev_print_energy = total_dft_energy;
    
    utils::logmesg(lmp, fmt::format("{:4d} {:15.8f} {:12.5e} {:12.5e}\n",
                                    current_iteration, total_dft_energy,
                                    energy_change, density_change));
  }
}

void PairDFT::print_scf_summary()
{
  if (comm->me == 0) {
  
    double energy_conversion = 1.0;
    std::string energy_units;
    
    if (strcmp(update->unit_style,"metal") == 0) {
      energy_conversion = 27.211386245981;
      energy_units = "eV";
    } else if (strcmp(update->unit_style,"real") == 0) {
      energy_conversion = 627.5094740631;
      energy_units = "kcal/mol";
    }
    utils::logmesg(lmp, "------------------------------------------------\n");
    if (scf_converged) utils::logmesg(lmp, "SCF CONVERGED\n");
    else utils::logmesg(lmp, "SCF NOT CONVERGED\n");
    utils::logmesg(lmp, "\nEnergy Components:\n");
    
    utils::logmesg(lmp, fmt::format("  Kinetic: {:15.8f} Eh {:15.8f} {}\n",
      kinetic_energy, kinetic_energy*energy_conversion, energy_units));
    
    utils::logmesg(lmp, fmt::format("  Nuclear: {:15.8f} Eh {:15.8f} {}\n",
      nuclear_repulsion, nuclear_repulsion*energy_conversion, energy_units));
    
    utils::logmesg(lmp, fmt::format("  XC:      {:15.8f} Eh {:15.8f} {}\n",
      xc_energy, xc_energy*energy_conversion, energy_units));
    
    utils::logmesg(lmp, fmt::format("  TOTAL:   {:15.8f} Eh {:15.8f} {}\n",
      total_dft_energy, total_dft_energy*energy_conversion, energy_units));
    
    utils::logmesg(lmp, "================================================\n\n");
  }
}
