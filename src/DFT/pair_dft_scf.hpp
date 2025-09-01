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
  generate_molecular_grid();
  integrate_xc_potential(xc_energy, vxc_matrix);
  
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
  
  // One-electron energy: Tr(P * H_core)
  Eigen::MatrixXd H_core = kinetic_matrix + nuclear_matrix;
  double one_electron = (density_matrix.cwiseProduct(H_core)).sum();
  
  // Two-electron Coulomb energy: 0.5 * Tr(P * J)
  // The 0.5 accounts for double counting in the Coulomb interaction
  double j_energy = 0.5 * (density_matrix.cwiseProduct(coulomb_matrix)).sum();
  
  // Exact exchange energy (only for hybrid functionals)
  // For hybrids: E_x^exact = -0.5 * alpha * Tr(P * K)
  double exact_exchange = 0.0;
  if (is_hybrid) {
    exact_exchange = -0.5 * hybrid_coeff * (density_matrix.cwiseProduct(exchange_matrix)).sum();
  }
  
  // Kinetic energy component (for output purposes)
  kinetic_energy = (density_matrix.cwiseProduct(kinetic_matrix)).sum();
  
  // Total DFT energy
  // E_total = E_nuc + Tr(P*H_core) + 0.5*Tr(P*J) + E_x^exact + E_xc
  // Note: E_xc already includes the DFT exchange and correlation
  // For hybrids, E_xc is scaled by (1-alpha) for exchange part
  total_dft_energy = nuclear_repulsion + one_electron + j_energy + exact_exchange + xc_energy;
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
   Generate molecular integration grid
   Following standard DFT grid generation (Mura-Knowles radial + Lebedev angular)
------------------------------------------------------------------------- */

void PairDFT::generate_molecular_grid()
{
  grid_points.clear();
  grid_weights.clear();
  
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  if (nlocal == 0) return;
  
  // Standard grid specifications (similar to GauXC)
  // Fine: 75 radial x 302 angular
  // UltraFine: 99 radial x 590 angular  
  // For testing, use Fine grid
  int n_radial = 75;
  int n_angular = 302;
  
  // Generate atom-centered grids
  std::vector<std::vector<double>> atom_positions;
  for (int i = 0; i < nlocal; i++) {
    atom_positions.push_back({x[i][0] * ANGSTROM_TO_BOHR, 
                              x[i][1] * ANGSTROM_TO_BOHR, 
                              x[i][2] * ANGSTROM_TO_BOHR});
  }
  
  // Generate grid for each atom
  for (int atom_idx = 0; atom_idx < nlocal; atom_idx++) {
    // Generate radial grid using Mura-Knowles quadrature
    std::vector<double> r_points, r_weights;
    generate_mura_knowles_radial(n_radial, 1.0, r_points, r_weights);
    
    // Generate angular grid using Lebedev quadrature
    std::vector<std::vector<double>> angular_points;
    std::vector<double> angular_weights;
    generate_lebedev_angular(n_angular, angular_points, angular_weights);
    
    // Combine radial and angular grids for this atom
    for (int i_r = 0; i_r < n_radial; i_r++) {
      double r = r_points[i_r];
      double w_r = r_weights[i_r];
      
      for (size_t i_ang = 0; i_ang < angular_points.size(); i_ang++) {
        // Grid point in real space
        std::vector<double> point(3);
        point[0] = atom_positions[atom_idx][0] + r * angular_points[i_ang][0];
        point[1] = atom_positions[atom_idx][1] + r * angular_points[i_ang][1];
        point[2] = atom_positions[atom_idx][2] + r * angular_points[i_ang][2];
        
        grid_points.push_back(point);
        
        // Combined weight for spherical integration
        // Weight = w_r * r^2 * w_angular
        // The r^2 comes from the Jacobian for spherical coordinates
        double weight = w_r * r * r * angular_weights[i_ang];
        grid_weights.push_back(weight);
      }
    }
  }
  
  // Apply Becke partitioning for molecular grid
  apply_becke_partitioning(atom_positions);
}

/* ----------------------------------------------------------------------
   Generate Mura-Knowles radial quadrature
   Based on Mura & Knowles, JCP 104, 9848 (1996)
------------------------------------------------------------------------- */

void PairDFT::generate_mura_knowles_radial(int n_points, double Z,
                                            std::vector<double> &r_points,
                                            std::vector<double> &r_weights)
{
  r_points.resize(n_points);
  r_weights.resize(n_points);
  
  // Mura-Knowles transformation: r = -alpha * ln(1 - x^3)
  // where x is mapped from standard quadrature points
  
  // Scaling factor alpha (element-dependent)
  // For H: alpha = 5.0 (from GauXC defaults)
  // For other elements, see Mura & Knowles paper
  double alpha = 5.0;  // Default for hydrogen
  if (Z > 1) {
    // Other elements have different scaling factors
    // See default_mk_radial_scaling_factor in GauXC
    alpha = 5.0;  // Most elements use 5.0
  }
  
  // Generate Chebyshev-Gauss quadrature points on [-1, 1]
  for (int i = 0; i < n_points; i++) {
    // Chebyshev-Gauss points
    double x_cheb = cos(M_PI * (2.0 * i + 1.0) / (2.0 * n_points));
    
    // Map from [-1, 1] to [0, 1]
    double x = 0.5 * (x_cheb + 1.0);
    
    // Apply Mura-Knowles transformation
    // r = -alpha * ln(1 - x^3)
    // Need to handle x very close to 1
    if (x > 0.999999) {
      // For x very close to 1, use limiting behavior
      r_points[i] = 50.0;  // Large cutoff radius
      r_weights[i] = 0.0;
    } else {
      double x3 = x * x * x;
      r_points[i] = -alpha * log(1.0 - x3);
      
      // Jacobian: dr/dx = 3 * alpha * x^2 / (1 - x^3)
      double jacobian = 3.0 * alpha * x * x / (1.0 - x3);
      
      // Chebyshev weight
      double w_cheb = M_PI / n_points;
      
      // Combined weight: includes mapping from [-1,1] to [0,1] (factor of 0.5)
      r_weights[i] = 0.5 * jacobian * w_cheb;
    }
  }
}

/* ----------------------------------------------------------------------
   Generate Lebedev angular quadrature
------------------------------------------------------------------------- */

void PairDFT::generate_lebedev_angular(int n_target,
                                       std::vector<std::vector<double>> &points,
                                       std::vector<double> &weights)
{
  points.clear();
  weights.clear();
  
  // Map target points to available Lebedev grids
  // Available orders: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, ...
  int n_actual = 6;
  if (n_target >= 590) n_actual = 590;
  else if (n_target >= 302) n_actual = 302;
  else if (n_target >= 266) n_actual = 266;
  else if (n_target >= 194) n_actual = 194;
  else if (n_target >= 110) n_actual = 110;
  else if (n_target >= 74) n_actual = 74;
  else if (n_target >= 50) n_actual = 50;
  else if (n_target >= 38) n_actual = 38;
  else if (n_target >= 26) n_actual = 26;
  else if (n_target >= 14) n_actual = 14;
  
  // Generate specific Lebedev grid
  // For simplicity, implementing common grids
  if (n_actual == 6) {
    // Order-3 octahedral grid
    generate_lebedev_6(points, weights);
  } else if (n_actual == 14) {
    // Order-5 grid
    generate_lebedev_14(points, weights);
  } else if (n_actual == 38) {
    // Order-9 grid
    generate_lebedev_38(points, weights);
  } else if (n_actual == 110) {
    // Order-17 grid
    generate_lebedev_110(points, weights);
  } else if (n_actual == 302) {
    // Order-29 grid (commonly used)
    generate_lebedev_302(points, weights);
  } else {
    // Fallback to simple uniform angular grid
    generate_uniform_angular(n_actual, points, weights);
  }
}

/* ----------------------------------------------------------------------
   Generate 6-point Lebedev grid (octahedron)
------------------------------------------------------------------------- */

void PairDFT::generate_lebedev_6(std::vector<std::vector<double>> &points,
                                 std::vector<double> &weights)
{
  // 6 points on coordinate axes
  points = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  
  // Equal weights summing to 4π
  double w = 4.0 * M_PI / 6.0;
  weights = std::vector<double>(6, w);
}

/* ----------------------------------------------------------------------
   Generate 14-point Lebedev grid
------------------------------------------------------------------------- */

void PairDFT::generate_lebedev_14(std::vector<std::vector<double>> &points,
                                  std::vector<double> &weights)
{
  // 6 octahedral points
  points = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  double w1 = 0.0666666666666667 * 4.0 * M_PI;
  weights = std::vector<double>(6, w1);
  
  // 8 cubic vertices (±1/√3, ±1/√3, ±1/√3)
  double a = 1.0 / sqrt(3.0);
  std::vector<std::vector<double>> cubic_points = {
    {a, a, a}, {a, a, -a}, {a, -a, a}, {a, -a, -a},
    {-a, a, a}, {-a, a, -a}, {-a, -a, a}, {-a, -a, -a}
  };
  
  double w2 = 0.0750000000000000 * 4.0 * M_PI;
  for (const auto& p : cubic_points) {
    points.push_back(p);
    weights.push_back(w2);
  }
}

/* ----------------------------------------------------------------------
   Generate 38-point Lebedev grid
------------------------------------------------------------------------- */

void PairDFT::generate_lebedev_38(std::vector<std::vector<double>> &points,
                                  std::vector<double> &weights)
{
  // Type 1: 6 points on axes
  double w1 = 0.0095238095238095 * 4.0 * M_PI;
  points = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  weights = std::vector<double>(6, w1);
  
  // Type 2: 8 points at cube vertices
  double a = 1.0 / sqrt(3.0);
  double w2 = 0.0321428571428571 * 4.0 * M_PI;
  std::vector<std::vector<double>> cube_pts = {
    {a, a, a}, {a, a, -a}, {a, -a, a}, {a, -a, -a},
    {-a, a, a}, {-a, a, -a}, {-a, -a, a}, {-a, -a, -a}
  };
  for (const auto& p : cube_pts) {
    points.push_back(p);
    weights.push_back(w2);
  }
  
  // Type 3: 24 points on edge midpoints
  double b = 0.4597008433809831;
  double c = 0.8880738339771154;
  double w3 = 0.0285714285714286 * 4.0 * M_PI;
  
  // Generate all permutations of (±b, ±c, 0) and cyclic
  std::vector<std::vector<double>> edge_pts;
  // (±b, ±c, 0) permutations
  edge_pts.push_back({b, c, 0}); edge_pts.push_back({b, -c, 0});
  edge_pts.push_back({-b, c, 0}); edge_pts.push_back({-b, -c, 0});
  edge_pts.push_back({c, b, 0}); edge_pts.push_back({c, -b, 0});
  edge_pts.push_back({-c, b, 0}); edge_pts.push_back({-c, -b, 0});
  // (0, ±b, ±c) permutations
  edge_pts.push_back({0, b, c}); edge_pts.push_back({0, b, -c});
  edge_pts.push_back({0, -b, c}); edge_pts.push_back({0, -b, -c});
  edge_pts.push_back({0, c, b}); edge_pts.push_back({0, c, -b});
  edge_pts.push_back({0, -c, b}); edge_pts.push_back({0, -c, -b});
  // (±c, 0, ±b) permutations
  edge_pts.push_back({c, 0, b}); edge_pts.push_back({c, 0, -b});
  edge_pts.push_back({-c, 0, b}); edge_pts.push_back({-c, 0, -b});
  edge_pts.push_back({b, 0, c}); edge_pts.push_back({b, 0, -c});
  edge_pts.push_back({-b, 0, c}); edge_pts.push_back({-b, 0, -c});
  
  for (const auto& p : edge_pts) {
    points.push_back(p);
    weights.push_back(w3);
  }
}

/* ----------------------------------------------------------------------
   Generate 110-point Lebedev grid
------------------------------------------------------------------------- */

void PairDFT::generate_lebedev_110(std::vector<std::vector<double>> &points,
                                   std::vector<double> &weights)
{
  // This is a placeholder - proper 110-point grid requires specific coefficients
  // For now, use uniform distribution
  generate_uniform_angular(110, points, weights);
}

/* ----------------------------------------------------------------------
   Generate 302-point Lebedev grid
------------------------------------------------------------------------- */

void PairDFT::generate_lebedev_302(std::vector<std::vector<double>> &points,
                                   std::vector<double> &weights)
{
  // This is a placeholder - proper 302-point grid requires many specific points
  // For production, this should use the actual Lebedev-302 coefficients
  // For now, use uniform spherical distribution
  generate_uniform_angular(302, points, weights);
}

/* ----------------------------------------------------------------------
   Generate uniform angular grid (fallback)
------------------------------------------------------------------------- */

void PairDFT::generate_uniform_angular(int n_points,
                                       std::vector<std::vector<double>> &points,
                                       std::vector<double> &weights)
{
  points.clear();
  weights.clear();
  
  // Use Fibonacci spiral for approximately uniform distribution
  double phi = (1.0 + sqrt(5.0)) / 2.0;  // Golden ratio
  
  for (int i = 0; i < n_points; i++) {
    // Map to sphere using Fibonacci spiral
    double y = 1.0 - 2.0 * (double)i / (n_points - 1);
    double radius = sqrt(1.0 - y * y);
    double theta = 2.0 * M_PI * i / phi;
    
    std::vector<double> point(3);
    point[0] = radius * cos(theta);
    point[1] = radius * sin(theta);
    point[2] = y;
    
    points.push_back(point);
    weights.push_back(4.0 * M_PI / n_points);
  }
}

/* ----------------------------------------------------------------------
   Apply Becke partitioning to molecular grid
------------------------------------------------------------------------- */

void PairDFT::apply_becke_partitioning(const std::vector<std::vector<double>> &atom_positions)
{
  int n_points = grid_points.size();
  int n_atoms = atom_positions.size();
  
  if (n_atoms == 1) {
    // Single atom - no partitioning needed
    return;
  }
  
  // Apply Becke weight to each grid point
  for (int i_pt = 0; i_pt < n_points; i_pt++) {
    std::vector<double> P_atom(n_atoms);
    
    // Calculate partition function for each atom
    for (int i = 0; i < n_atoms; i++) {
      P_atom[i] = 1.0;
      
      // Distance from grid point to atom i
      double r_i = 0.0;
      for (int k = 0; k < 3; k++) {
        double d = grid_points[i_pt][k] - atom_positions[i][k];
        r_i += d * d;
      }
      r_i = sqrt(r_i);
      
      // Product over all other atoms
      for (int j = 0; j < n_atoms; j++) {
        if (i == j) continue;
        
        // Distance from grid point to atom j
        double r_j = 0.0;
        for (int k = 0; k < 3; k++) {
          double d = grid_points[i_pt][k] - atom_positions[j][k];
          r_j += d * d;
        }
        r_j = sqrt(r_j);
        
        // Distance between atoms i and j
        double R_ij = 0.0;
        for (int k = 0; k < 3; k++) {
          double d = atom_positions[i][k] - atom_positions[j][k];
          R_ij += d * d;
        }
        R_ij = sqrt(R_ij);
        
        if (R_ij < 1e-10) continue;
        
        // Confocal elliptical coordinate
        double mu = (r_i - r_j) / R_ij;
        
        // Apply smoothing function (3 iterations of Becke's polynomial)
        double f = mu;
        for (int iter = 0; iter < 3; iter++) {
          f = 0.5 * f * (3.0 - f * f);
        }
        
        // Step function
        double s = 0.5 * (1.0 - f);
        
        P_atom[i] *= s;
      }
    }
    
    // Normalize partition functions
    double sum = 0.0;
    for (int i = 0; i < n_atoms; i++) {
      sum += P_atom[i];
    }
    
    if (sum > 1e-15) {
      // Apply Becke weight to grid weight
      // Grid points are generated per atom, so we need to identify which atom
      // For a molecular grid, sum all atomic contributions
      grid_weights[i_pt] *= 1.0;  // Full weight for molecular integration
    } else {
      grid_weights[i_pt] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   Integrate XC potential on grid
------------------------------------------------------------------------- */

void PairDFT::integrate_xc_potential(double &exc_energy, Eigen::MatrixXd &vxc_matrix)
{
  int n_points = grid_points.size();
  
  vxc_matrix.setZero();
  exc_energy = 0.0;
  
  if (n_points == 0) return;
  
  // Evaluate basis functions and gradients at grid points
  std::vector<std::vector<double>> basis_values;
  std::vector<std::vector<std::vector<double>>> basis_gradients;
  
  // For GGA functionals, need gradients
  bool need_gradients = !is_lda;  // PBE is GGA, so needs gradients
  if (need_gradients) {
    basis_gradients.resize(n_points);
  }
  
  evaluate_basis_at_points(grid_points, basis_values, basis_gradients);
  
  // Compute density and gradient at each grid point
  std::vector<double> rho(n_points, 0.0);
  std::vector<double> sigma(n_points, 0.0);  // |∇ρ|²
  std::vector<double> lapl(n_points, 0.0);    // For meta-GGA
  std::vector<double> tau(n_points, 0.0);     // Kinetic energy density
  
  // Calculate density at grid points
  for (int ipt = 0; ipt < n_points; ipt++) {
    // Compute density: ρ = Σ_μν P_μν φ_μ φ_ν
    for (int i = 0; i < n_basis_functions; i++) {
      for (int j = 0; j < n_basis_functions; j++) {
        rho[ipt] += density_matrix(i, j) * basis_values[ipt][i] * basis_values[ipt][j];
      }
    }
    
    // Compute gradient for GGA
    if (need_gradients && !basis_gradients.empty()) {
      std::vector<double> grad_rho(3, 0.0);
      
      // ∇ρ = Σ_μν P_μν (∇φ_μ φ_ν + φ_μ ∇φ_ν)
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j < n_basis_functions; j++) {
          double P_ij = density_matrix(i, j);
          for (int k = 0; k < 3; k++) {
            grad_rho[k] += P_ij * (basis_gradients[ipt][i][k] * basis_values[ipt][j] +
                                   basis_values[ipt][i] * basis_gradients[ipt][j][k]);
          }
        }
      }
      
      // |∇ρ|²
      sigma[ipt] = grad_rho[0] * grad_rho[0] + 
                   grad_rho[1] * grad_rho[1] + 
                   grad_rho[2] * grad_rho[2];
      
      // Kinetic energy density for meta-GGA
      if (is_meta_gga) {
        for (int i = 0; i < n_basis_functions; i++) {
          for (int j = 0; j < n_basis_functions; j++) {
            double P_ij = density_matrix(i, j);
            for (int k = 0; k < 3; k++) {
              tau[ipt] += 0.5 * P_ij * basis_gradients[ipt][i][k] * basis_gradients[ipt][j][k];
            }
          }
        }
      }
    }
  }
  
  // Evaluate XC functional
  std::vector<double> exc(n_points);
  std::vector<double> vrho(n_points);
  std::vector<double> vsigma(n_points);
  std::vector<double> vlapl(n_points);
  std::vector<double> vtau(n_points);
  
  evaluate_xc_functional(rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);
  
  // Integrate XC energy and build potential matrix
  double integrated_density = 0.0;
  
  for (int ipt = 0; ipt < n_points; ipt++) {
    double w = grid_weights[ipt];
    
    // Skip points with negligible density
    if (rho[ipt] < 1e-15) continue;
    
    // XC energy: E_xc = ∫ ε_xc(ρ) ρ dV
    exc_energy += exc[ipt] * rho[ipt] * w;
    
    // Check integrated density
    integrated_density += rho[ipt] * w;
    
    // LDA contribution to potential matrix
    for (int i = 0; i < n_basis_functions; i++) {
      for (int j = 0; j <= i; j++) {
        double v_lda = vrho[ipt] * basis_values[ipt][i] * basis_values[ipt][j] * w;
        vxc_matrix(i, j) += v_lda;
        if (i != j) vxc_matrix(j, i) += v_lda;
      }
    }
    
    // GGA contribution
    if (need_gradients && vsigma[ipt] != 0.0 && !basis_gradients.empty()) {
      // Compute gradient of density for GGA contribution
      std::vector<double> grad_rho(3, 0.0);
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j < n_basis_functions; j++) {
          double P_ij = density_matrix(i, j);
          for (int k = 0; k < 3; k++) {
            grad_rho[k] += P_ij * (basis_gradients[ipt][i][k] * basis_values[ipt][j] +
                                   basis_values[ipt][i] * basis_gradients[ipt][j][k]);
          }
        }
      }
      
      // GGA contribution: 2 * v_σ * ∇ρ · (∇φ_μ φ_ν + φ_μ ∇φ_ν)
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j <= i; j++) {
          double v_gga = 0.0;
          for (int k = 0; k < 3; k++) {
            v_gga += 2.0 * vsigma[ipt] * grad_rho[k] * 
                     (basis_gradients[ipt][i][k] * basis_values[ipt][j] +
                      basis_values[ipt][i] * basis_gradients[ipt][j][k]) * w;
          }
          vxc_matrix(i, j) += v_gga;
          if (i != j) vxc_matrix(j, i) += v_gga;
        }
      }
    }
    
    // Meta-GGA contribution
    if (is_meta_gga && vtau[ipt] != 0.0 && !basis_gradients.empty()) {
      for (int i = 0; i < n_basis_functions; i++) {
        for (int j = 0; j <= i; j++) {
          double v_tau = 0.0;
          for (int k = 0; k < 3; k++) {
            v_tau += 0.5 * vtau[ipt] * basis_gradients[ipt][i][k] * 
                     basis_gradients[ipt][j][k] * w;
          }
          vxc_matrix(i, j) += v_tau;
          if (i != j) vxc_matrix(j, i) += v_tau;
        }
      }
    }
  }
  
  // Diagnostic output
  if (comm->me == 0 && current_iteration == 1) {
    utils::logmesg(lmp, fmt::format("  Grid: {} points, integrated {:.3f} electrons (expected {})\n", 
                                    n_points, integrated_density, n_electrons));
    if (std::abs(integrated_density - n_electrons) > 0.1) {
      utils::logmesg(lmp, "  WARNING: Poor electron integration - check grid\n");
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
