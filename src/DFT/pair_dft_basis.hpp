/* ----------------------------------------------------------------------
   Basis set methods for PairDFT
------------------------------------------------------------------------- */

#include "json.h"

#include <fstream>
#include <iostream>
#include <cmath>


/* ----------------------------------------------------------------------
   Initialize basis set from JSON file
------------------------------------------------------------------------- */

void PairDFT::initialize_basis_set()
{
  // Load basis from JSON
  load_basis_from_json(basis_file);
  
  // Initialize libint2
  initialize_libint();
  
  // Allocate matrices
  overlap_matrix.resize(n_basis_functions, n_basis_functions);
  kinetic_matrix.resize(n_basis_functions, n_basis_functions);
  nuclear_matrix.resize(n_basis_functions, n_basis_functions);
  coulomb_matrix.resize(n_basis_functions, n_basis_functions);
  exchange_matrix.resize(n_basis_functions, n_basis_functions);
  fock_matrix.resize(n_basis_functions, n_basis_functions);
  density_matrix.resize(n_basis_functions, n_basis_functions);
  density_matrix_prev.resize(n_basis_functions, n_basis_functions);
  mo_coefficients.resize(n_basis_functions, n_basis_functions);
  mo_energies.resize(n_basis_functions);
  
  // Initialize all matrices to zero
  overlap_matrix.setZero();
  kinetic_matrix.setZero();
  nuclear_matrix.setZero();
  coulomb_matrix.setZero();
  exchange_matrix.setZero();
  fock_matrix.setZero();
  density_matrix.setZero();
  density_matrix_prev.setZero();
  mo_coefficients.setZero();
  mo_energies.setZero();
  
  // Compute one-electron integrals
  compute_one_electron_integrals();
  
  // Initialize density guess
  initialize_density_guess();
}

/* ----------------------------------------------------------------------
   Load basis set from JSON file
------------------------------------------------------------------------- */

void PairDFT::load_basis_from_json(const std::string &filename) 
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    error->all(FLERR, "Cannot open basis set file: " + filename);
  }
  
  json j;
  file >> j;
  file.close();
  
  // Parse basis set from JSON
  n_shells = 0;
  n_basis_functions = 0;
  
  // Clear existing data
  shell_to_atom.clear();
  angular_momentum.clear();
  exponents.clear();
  coefficients.clear();
  normalized_coefficients.clear();
  
  // Parse elements
  if (j.contains("elements")) {
    for (auto& [element_str, element_data] : j["elements"].items()) {
      int atomic_number = std::stoi(element_str);
      
      if (element_data.contains("electron_shells")) {
        // Store shells for this element (will be assigned to atoms later)
        std::vector<int> element_am;
        std::vector<std::vector<double>> element_exp;
        std::vector<std::vector<double>> element_coeff;
        
        for (auto& shell : element_data["electron_shells"]) {
          int am = -1;
          
          // Get angular momentum
          if (shell.contains("angular_momentum")) {
            auto am_array = shell["angular_momentum"];
            if (am_array.is_array() && am_array.size() > 0) {
              am = am_array[0].get<int>();
            }
          }
          
          if (am < 0) continue;
          
          element_am.push_back(am);
          
          // Get exponents
          std::vector<double> shell_exponents;
          if (shell.contains("exponents")) {
            for (auto& exp_str : shell["exponents"]) {
              double exp_val = std::stod(exp_str.get<std::string>());
              shell_exponents.push_back(exp_val);
            }
          }
          element_exp.push_back(shell_exponents);
          
          // Get coefficients
          std::vector<double> shell_coefficients;
          if (shell.contains("coefficients")) {
            auto coeff_array = shell["coefficients"];
            if (coeff_array.is_array() && coeff_array.size() > 0) {
              for (auto& coeff_str : coeff_array[0]) {
                double coeff_val = std::stod(coeff_str.get<std::string>());
                shell_coefficients.push_back(coeff_val);
              }
            }
          }
          element_coeff.push_back(shell_coefficients);
        }
        
        // Now assign shells to each atom of this element type
        // For H2, both atoms are hydrogen (element 1)
        int nlocal = atom->nlocal;
        for (int atom_idx = 0; atom_idx < nlocal; atom_idx++) {
          // TODO: Check atom type to match with element
          // For now, assume all atoms are of this element type
          
          // Add all shells for this element to this atom
          for (size_t s = 0; s < element_am.size(); s++) {
            angular_momentum.push_back(element_am[s]);
            exponents.push_back(element_exp[s]);
            coefficients.push_back(element_coeff[s]);
            shell_to_atom.push_back(atom_idx);
            n_shells++;
            
            // Count basis functions
            // s: 1, p: 3, d: 6, f: 10, g: 15
            int n_funcs = (element_am[s] + 1) * (element_am[s] + 2) / 2;
            n_basis_functions += n_funcs;
            
            if (comm->me == 0) {
              utils::logmesg(lmp, fmt::format("Added shell {} for atom {}: l={}, n_exp={}, n_coeff={}, n_funcs={}\n",
                                             n_shells-1, atom_idx, element_am[s], 
                                             element_exp[s].size(), element_coeff[s].size(), n_funcs));
            }
          }
        }
      }
    }
  }
  
  // Skip normalization - basis sets from BSE are already normalized
  // normalize_basis_functions();
  // Don't set normalized_coefficients, just use coefficients directly
  
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("Loaded basis set with {} shells and {} basis functions\n", 
                                    n_shells, n_basis_functions));
    utils::logmesg(lmp, fmt::format("  exponents.size() = {}\n", exponents.size()));
    utils::logmesg(lmp, fmt::format("  coefficients.size() = {}\n", coefficients.size()));
    utils::logmesg(lmp, fmt::format("  angular_momentum.size() = {}\n", angular_momentum.size()));
    utils::logmesg(lmp, fmt::format("  shell_to_atom.size() = {}\n", shell_to_atom.size()));
  }
}

/* ----------------------------------------------------------------------
   Normalize basis functions
------------------------------------------------------------------------- */

void PairDFT::normalize_basis_functions()
{
  normalized_coefficients.clear();
  
  for (int shell = 0; shell < n_shells; shell++) {
    int l = angular_momentum[shell];
    std::vector<double> norm_coeffs;
    
    for (size_t i = 0; i < exponents[shell].size(); i++) {
      double alpha = exponents[shell][i];
      double norm = compute_normalization(l, alpha);
      double coeff = coefficients[shell][i];
      norm_coeffs.push_back(coeff * norm);
    }
    
    normalized_coefficients.push_back(norm_coeffs);
  }
}

/* ----------------------------------------------------------------------
   Compute normalization constant for Gaussian basis function
------------------------------------------------------------------------- */

double PairDFT::compute_normalization(int l, double exponent)
{
  // Normalization constant for Gaussian basis function
  // N = (2*alpha/pi)^(3/4) * sqrt((8*alpha)^l * (2l-1)!! / (2l+1)!!)
  
  double pi = M_PI;
  double norm = pow(2.0 * exponent / pi, 0.75);
  
  // Add angular momentum dependent part
  if (l > 0) {
    double factor = pow(2.0, l) * pow(exponent, l/2.0);
    
    // Compute double factorial
    int double_fact_numerator = 1;
    for (int i = 2*l - 1; i > 0; i -= 2) double_fact_numerator *= i;
    
    int double_fact_denominator = 1;
    for (int i = 2*l + 1; i > 0; i -= 2) double_fact_denominator *= i;
    
    factor *= sqrt(double(double_fact_numerator) / double(double_fact_denominator));
    norm *= factor;
  }
  
  return norm;
}

/* ----------------------------------------------------------------------
   Get angular momentum for shell
------------------------------------------------------------------------- */

int PairDFT::get_angular_momentum(int shell) const
{
  if (shell < 0 || shell >= n_shells) return -1;
  return angular_momentum[shell];
}

/* ----------------------------------------------------------------------
   Get exponents for shell
------------------------------------------------------------------------- */

std::vector<double> PairDFT::get_exponents(int shell) const
{
  if (shell < 0 || shell >= n_shells || shell >= exponents.size()) 
    return std::vector<double>();
  return exponents[shell];
}

/* ----------------------------------------------------------------------
   Get coefficients for shell
------------------------------------------------------------------------- */

std::vector<double> PairDFT::get_coefficients(int shell) const
{
  if (shell < 0 || shell >= n_shells || shell >= coefficients.size()) 
    return std::vector<double>();
  // Always return from coefficients since we're not normalizing
  return coefficients[shell];
}

/* ----------------------------------------------------------------------
   Get Cartesian Gaussian function indices for angular momentum
------------------------------------------------------------------------- */

void PairDFT::get_cartesian_indices(int l, std::vector<std::vector<int>> &indices)
{
  indices.clear();
  
  // Generate all combinations of (nx, ny, nz) where nx + ny + nz = l
  for (int nx = 0; nx <= l; nx++) {
    for (int ny = 0; ny <= l - nx; ny++) {
      int nz = l - nx - ny;
      indices.push_back({nx, ny, nz});
    }
  }
}

/* ----------------------------------------------------------------------
   Evaluate basis functions at grid points
------------------------------------------------------------------------- */

void PairDFT::evaluate_basis_at_points(const std::vector<std::vector<double>> &points,
                                       std::vector<std::vector<double>> &basis_values,
                                       std::vector<std::vector<std::vector<double>>> &basis_gradients)
{
  int n_points = points.size();
  basis_values.resize(n_points, std::vector<double>(n_basis_functions));
  
  // Debug check
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("evaluate_basis_at_points: n_shells={}, exponents.size()={}, coefficients.size()={}\n",
                                    n_shells, exponents.size(), coefficients.size()));
  }
  
  bool need_gradients = !basis_gradients.empty();
  if (need_gradients) {
    basis_gradients.resize(n_points, 
                          std::vector<std::vector<double>>(n_basis_functions, 
                                                          std::vector<double>(3)));
  }
  
  // Get current atom positions from LAMMPS atom structure
  double **x = atom->x;
  int nlocal = atom->nlocal;
  
  // For each grid point, evaluate all basis functions
  for (int i_point = 0; i_point < n_points; i_point++) {
    int bf_idx = 0;
    
    for (int shell = 0; shell < n_shells; shell++) {
      int l = angular_momentum[shell];
      
      // Bounds checking for exponents and coefficients
      if (shell >= exponents.size() || shell >= coefficients.size()) {
        if (comm->me == 0) {
          error->warning(FLERR, fmt::format("Shell {} out of bounds (exp size={}, coeff size={})", 
                                           shell, exponents.size(), coefficients.size()));
        }
        continue;
      }
      
      auto shell_exponents = exponents[shell];
      auto shell_coefficients = coefficients[shell];
      int atom_id = shell_to_atom[shell];
      
      // Skip if exponents or coefficients are empty
      if (shell_exponents.empty() || shell_coefficients.empty()) {
        if (comm->me == 0) {
          error->warning(FLERR, fmt::format("Shell {} has empty exponents or coefficients", shell));
        }
        // Still need to advance bf_idx for the basis functions of this shell
        std::vector<std::vector<int>> cart_indices;
        get_cartesian_indices(l, cart_indices);
        bf_idx += cart_indices.size();
        continue;
      }
      
      // Get atom position directly from LAMMPS
      std::vector<double> atom_pos(3, 0.0);
      if (atom_id >= 0 && atom_id < nlocal) {
        atom_pos[0] = x[atom_id][0] * ANGSTROM_TO_BOHR;
        atom_pos[1] = x[atom_id][1] * ANGSTROM_TO_BOHR;
        atom_pos[2] = x[atom_id][2] * ANGSTROM_TO_BOHR;
      }
      
      // Calculate distance vector from atom to grid point
      double dx = points[i_point][0] - atom_pos[0];
      double dy = points[i_point][1] - atom_pos[1];
      double dz = points[i_point][2] - atom_pos[2];
      double r2 = dx*dx + dy*dy + dz*dz;
      
      // Get Cartesian indices for this angular momentum
      std::vector<std::vector<int>> cart_indices;
      get_cartesian_indices(l, cart_indices);
      
      // Evaluate each Cartesian Gaussian in the shell
      for (size_t func = 0; func < cart_indices.size(); func++) {
        int nx = cart_indices[func][0];
        int ny = cart_indices[func][1];
        int nz = cart_indices[func][2];
        
        double value = 0.0;
        std::vector<double> gradient = {0.0, 0.0, 0.0};
        
        // Contract over primitives
        for (size_t prim = 0; prim < shell_exponents.size(); prim++) {
          double alpha = shell_exponents[prim];
          double coeff = shell_coefficients[prim];
          
          // Gaussian exponential part
          double gauss_exp = exp(-alpha * r2);
          
          // Cartesian polynomial part
          double poly_part = pow(dx, nx) * pow(dy, ny) * pow(dz, nz);
          
          // Full basis function value
          double prim_value = coeff * gauss_exp * poly_part;
          value += prim_value;
          
          // Compute gradients if needed
          if (need_gradients) {
            // Gradient of N * exp(-alpha*r^2) * x^nx * y^ny * z^nz
            // d/dx = N * exp(-alpha*r^2) * y^ny * z^nz * (nx * x^(nx-1) - 2*alpha*x * x^nx)
            
            double grad_x, grad_y, grad_z;
            
            // x-component
            if (nx > 0) {
              grad_x = coeff * gauss_exp * pow(dy, ny) * pow(dz, nz) * 
                      pow(dx, nx-1) * (nx - 2*alpha*dx*dx);
            } else {
              grad_x = coeff * gauss_exp * pow(dy, ny) * pow(dz, nz) * 
                      (-2*alpha*dx);
            }
            
            // y-component
            if (ny > 0) {
              grad_y = coeff * gauss_exp * pow(dx, nx) * pow(dz, nz) *
                      pow(dy, ny-1) * (ny - 2*alpha*dy*dy);
            } else {
              grad_y = coeff * gauss_exp * pow(dx, nx) * pow(dz, nz) *
                      (-2*alpha*dy);
            }
            
            // z-component
            if (nz > 0) {
              grad_z = coeff * gauss_exp * pow(dx, nx) * pow(dy, ny) *
                      pow(dz, nz-1) * (nz - 2*alpha*dz*dz);
            } else {
              grad_z = coeff * gauss_exp * pow(dx, nx) * pow(dy, ny) *
                      (-2*alpha*dz);
            }
            
            gradient[0] += grad_x;
            gradient[1] += grad_y; 
            gradient[2] += grad_z;
          }
        }
        
        basis_values[i_point][bf_idx] = value;
        
        if (need_gradients) {
          basis_gradients[i_point][bf_idx][0] = gradient[0];
          basis_gradients[i_point][bf_idx][1] = gradient[1];
          basis_gradients[i_point][bf_idx][2] = gradient[2];
        }
        
        bf_idx++;
      }
    }
  }
}

