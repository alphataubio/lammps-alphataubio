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
          
          angular_momentum.push_back(am);
          n_shells++;
          
          // Get exponents
          std::vector<double> shell_exponents;
          if (shell.contains("exponents")) {
            for (auto& exp_str : shell["exponents"]) {
              double exp_val = std::stod(exp_str.get<std::string>());
              shell_exponents.push_back(exp_val);
            }
          }
          exponents.push_back(shell_exponents);
          
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
          coefficients.push_back(shell_coefficients);
          
          // Count basis functions
          // s: 1, p: 3, d: 6, f: 10, g: 15
          int n_funcs = (am + 1) * (am + 2) / 2;
          n_basis_functions += n_funcs;
          
          // Map shells to atoms - for now, distribute shells among atoms
          // In a real implementation, this would need proper atom-type mapping
          // based on the actual atom types in the simulation
          int atom_index = 0;  // Default to first atom
          shell_to_atom.push_back(atom_index);
        }
      }
    }
  }
  
  // Normalize basis functions
  normalize_basis_functions();
  
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("Loaded basis set with {} shells and {} basis functions\n", 
                                    n_shells, n_basis_functions));
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
   Set atom positions for basis functions
------------------------------------------------------------------------- */

void PairDFT::set_atom_positions(const std::vector<std::vector<double>> &positions)
{
  atom_positions = positions;
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
  if (shell < 0 || shell >= n_shells) return std::vector<double>();
  return exponents[shell];
}

/* ----------------------------------------------------------------------
   Get coefficients for shell
------------------------------------------------------------------------- */

std::vector<double> PairDFT::get_coefficients(int shell) const
{
  if (shell < 0 || shell >= n_shells) return std::vector<double>();
  return normalized_coefficients.empty() ? coefficients[shell] : normalized_coefficients[shell];
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
  
  bool need_gradients = !basis_gradients.empty();
  if (need_gradients) {
    basis_gradients.resize(n_points, 
                          std::vector<std::vector<double>>(n_basis_functions, 
                                                          std::vector<double>(3)));
  }
  
  // For each grid point, evaluate all basis functions
  for (int i_point = 0; i_point < n_points; i_point++) {
    int bf_idx = 0;
    
    for (int shell = 0; shell < n_shells; shell++) {
      int l = angular_momentum[shell];
      auto shell_exponents = get_exponents(shell);
      auto shell_coefficients = get_coefficients(shell);
      int atom_id = shell_to_atom[shell];
      
      // Get atom position (default to origin if not set)
      std::vector<double> atom_pos = {0.0, 0.0, 0.0};
      if (atom_id < atom_positions.size()) {
        atom_pos = atom_positions[atom_id];
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
            // Gradient of exp(-alpha*r^2) * x^nx * y^ny * z^nz
            
            // x-component: d/dx [exp(-alpha*r^2) * x^nx * y^ny * z^nz]
            double grad_x = coeff * gauss_exp * pow(dy, ny) * pow(dz, nz) * 
                           (nx * pow(dx, std::max(0, nx-1)) - 2*alpha * dx * pow(dx, nx));
            
            // y-component: d/dy [exp(-alpha*r^2) * x^nx * y^ny * z^nz]  
            double grad_y = coeff * gauss_exp * pow(dx, nx) * pow(dz, nz) *
                           (ny * pow(dy, std::max(0, ny-1)) - 2*alpha * dy * pow(dy, ny));
            
            // z-component: d/dz [exp(-alpha*r^2) * x^nx * y^ny * z^nz]
            double grad_z = coeff * gauss_exp * pow(dx, nx) * pow(dy, ny) *
                           (nz * pow(dz, std::max(0, nz-1)) - 2*alpha * dz * pow(dz, nz));
            
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

