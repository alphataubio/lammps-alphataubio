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

#include <fstream>
#include <iostream>
#include <cmath>

#include "json.h"

/* ---------------------------------------------------------------------- */

void PairDFT::load_basis_from_json(const std::string &filename) {

  if (comm->me == 0) {
    fp = fopen(filename.c_str(), "r");
    if (fp == nullptr)
      error->one(FLERR, fileiarg, "Cannot open basis set file {}: {}", filename, utils::getsyserror());
    try {
      // try to parse as a JSON file. parser throws an exception on errors
      // if successful, temporarily serialize to bytearray for communication
      moldata = json::parse(fp);
      jsondata = json::to_ubjson(moldata);
      jsondata_size = jsondata.size();
      fclose(fp);
    } catch (std::exception &e) {
      fclose(fp);
      error->one(FLERR, fileiarg, "Error parsing JSON file {}: {}", filename, e.what());
    }
  }

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open basis set file: " + filename);
  }
  
  json j;
  file >> j;
  file.close();
  
  // Parse basis set from JSON
  // Expected format from basis set exchange
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
          
          // For now, assume one atom type - would need atom mapping in real use
          shell_to_atom.push_back(0);
        }
      }
    }
  }
  
  // Normalize basis functions
  normalize_basis_functions();
}

/* ----------------------------------------------------------------------
   Initialize basis set from JSON file
------------------------------------------------------------------------- */

void PairDFT::initialize_basis_set()
{
  // Create modular components
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

/* ---------------------------------------------------------------------- */

std::vector<double> PairDFT::get_exponents(int shell) const
{
  if (shell < 0 || shell >= n_shells) return std::vector<double>();
  return exponents[shell];
}

/* ---------------------------------------------------------------------- */

std::vector<double> PairDFT::get_coefficients(int shell) const
{
  if (shell < 0 || shell >= n_shells) {
    return std::vector<double>();
  }
  return normalized_coefficients.empty() ? coefficients[shell] : normalized_coefficients[shell];
}

/* ---------------------------------------------------------------------- */

int PairDFT::get_angular_momentum(int shell) const
{
  if (shell < 0 || shell >= n_shells) {
    return -1;
  }
  return angular_momentum[shell];
}

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

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
    for (int i = 2*l - 1; i > 0; i -= 2) {
      double_fact_numerator *= i;
    }
    
    int double_fact_denominator = 1;
    for (int i = 2*l + 1; i > 0; i -= 2) {
      double_fact_denominator *= i;
    }
    
    factor *= sqrt(double(double_fact_numerator) / double(double_fact_denominator));
    norm *= factor;
  }
  
  return norm;
}
