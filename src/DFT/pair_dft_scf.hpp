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
#include <cmath>
#include <vector>
#include <algorithm>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridIntegrator::GridIntegrator(int grid_size, const std::string &type) 
  : target_grid_size(grid_size), grid_type(type)
{
  grid_points.clear();
  grid_weights.clear();
  becke_weights.clear();
}

/* ---------------------------------------------------------------------- */

GridIntegrator::~GridIntegrator()
{
}

/* ---------------------------------------------------------------------- */

void GridIntegrator::generate_grid(const std::vector<std::vector<double>> &atom_positions,
                                   const std::vector<int> &atomic_numbers)
{
  grid_points.clear();
  grid_weights.clear();
  
  int n_atoms = atom_positions.size();
  if (n_atoms == 0) return;
  
  // Generate atomic grids
  int points_per_atom = target_grid_size / n_atoms;
  
  for (int atom = 0; atom < n_atoms; atom++) {
    // Get radial and angular grid sizes
    int n_radial = 50;  // Typical value
    int n_angular = points_per_atom / n_radial;
    
    // Generate radial grid (Chebyshev-Gauss)
    std::vector<double> r_points(n_radial);
    std::vector<double> r_weights(n_radial);
    generate_radial_grid(n_radial, atomic_numbers[atom], r_points, r_weights);
    
    // Generate angular grid (Lebedev)
    std::vector<std::vector<double>> angular_points;
    std::vector<double> angular_weights;
    generate_lebedev_grid(n_angular, angular_points, angular_weights);
    
    // Combine radial and angular grids
    for (int i_r = 0; i_r < n_radial; i_r++) {
      double r = r_points[i_r];
      double w_r = r_weights[i_r];
      
      for (size_t i_ang = 0; i_ang < angular_points.size(); i_ang++) {
        std::vector<double> point(3);
        point[0] = atom_positions[atom][0] + r * angular_points[i_ang][0];
        point[1] = atom_positions[atom][1] + r * angular_points[i_ang][1];
        point[2] = atom_positions[atom][2] + r * angular_points[i_ang][2];
        
        grid_points.push_back(point);
        grid_weights.push_back(w_r * angular_weights[i_ang] * r * r);
      }
    }
  }
  
  // Compute Becke partitioning weights
  compute_becke_weights(atom_positions);
}

/* ---------------------------------------------------------------------- */

void GridIntegrator::generate_radial_grid(int n_points, double Z,
                                          std::vector<double> &r,
                                          std::vector<double> &w)
{
  // Chebyshev-Gauss radial grid
  // Transform from [-1,1] to [0,inf) using appropriate mapping
  
  r.resize(n_points);
  w.resize(n_points);
  
  // Bragg radius for scaling
  double R_bragg = 1.0;  // Default, should depend on Z
  if (Z > 0) {
    R_bragg = 0.5 * (3.0 - 0.01 * Z);  // Simple approximation
  }
  
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

/* ---------------------------------------------------------------------- */

void GridIntegrator::generate_lebedev_grid(int n_points,
                                           std::vector<std::vector<double>> &points,
                                           std::vector<double> &weights)
{
  // Simplified Lebedev grid generation
  // Real implementation would use pre-computed Lebedev grids
  
  points.clear();
  weights.clear();
  
  // Find closest available Lebedev grid
  // Available sizes: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, ...
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
  // This is a simplified version - real Lebedev grids are more complex
  
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

/* ---------------------------------------------------------------------- */

void GridIntegrator::compute_becke_weights(const std::vector<std::vector<double>> &atoms)
{
  // Becke partitioning for multi-center integration
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
        for (int iter = 0; iter < 3; iter++) {
          f = 0.5 * f * (3.0 - f * f);
        }
        double s = 0.5 * (1.0 - f);
        
        P[i_atom] *= s;
      }
    }
    
    // Normalize partition functions
    double sum = 0.0;
    for (int i_atom = 0; i_atom < n_atoms; i_atom++) {
      sum += P[i_atom];
    }
    
    if (sum > 1e-15) {
      // For simplicity, assign full weight to nearest atom
      // Real implementation would properly partition weights
      becke_weights[i_point] = 1.0;
    } else {
      becke_weights[i_point] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridIntegrator::integrate_xc(XCFunctional *xc_func,
                                  const Eigen::MatrixXd &density_matrix,
                                  BasisSetManager *basis,
                                  double &exc_energy,
                                  Eigen::MatrixXd &vxc_matrix)
{
  int n_basis = basis->get_n_basis();
  int n_points = grid_points.size();
  
  vxc_matrix.setZero(n_basis, n_basis);
  exc_energy = 0.0;
  
  if (n_points == 0) return;
  
  // Evaluate basis functions at grid points
  std::vector<std::vector<double>> basis_values(n_points, std::vector<double>(n_basis));
  std::vector<std::vector<std::vector<double>>> basis_gradients;
  
  if (xc_func->is_gga() || xc_func->is_meta()) {
    basis_gradients.resize(n_points, 
                          std::vector<std::vector<double>>(n_basis, 
                                                          std::vector<double>(3)));
  }
  
  evaluate_basis_at_points(basis, basis_values, basis_gradients);
  
  // Compute density and gradients at grid points
  std::vector<double> rho(n_points, 0.0);
  std::vector<double> sigma(n_points, 0.0);  // |grad rho|^2
  std::vector<double> lapl(n_points, 0.0);   // Laplacian
  std::vector<double> tau(n_points, 0.0);     // Kinetic energy density
  
  for (int i_point = 0; i_point < n_points; i_point++) {
    // Density
    for (int i = 0; i < n_basis; i++) {
      for (int j = 0; j < n_basis; j++) {
        rho[i_point] += density_matrix(i, j) * 
                        basis_values[i_point][i] * 
                        basis_values[i_point][j];
      }
    }
    
    // Gradient and kinetic energy density for GGA/meta-GGA
    if (xc_func->is_gga() || xc_func->is_meta()) {
      std::vector<double> grad_rho(3, 0.0);
      
      for (int i = 0; i < n_basis; i++) {
        for (int j = 0; j < n_basis; j++) {
          double P_ij = density_matrix(i, j);
          
          for (int k = 0; k < 3; k++) {
            grad_rho[k] += P_ij * (basis_gradients[i_point][i][k] * basis_values[i_point][j] +
                                   basis_values[i_point][i] * basis_gradients[i_point][j][k]);
          }
          
          if (xc_func->is_meta()) {
            for (int k = 0; k < 3; k++) {
              tau[i_point] += 0.5 * P_ij * 
                             basis_gradients[i_point][i][k] * 
                             basis_gradients[i_point][j][k];
            }
          }
        }
      }
      
      sigma[i_point] = grad_rho[0] * grad_rho[0] + 
                       grad_rho[1] * grad_rho[1] + 
                       grad_rho[2] * grad_rho[2];
    }
  }
  
  // Evaluate XC functional
  std::vector<double> exc(n_points);
  std::vector<double> vrho(n_points);
  std::vector<double> vsigma(n_points);
  std::vector<double> vlapl(n_points);
  std::vector<double> vtau(n_points);
  
  xc_func->evaluate(rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau);
  
  // Integrate XC energy
  for (int i_point = 0; i_point < n_points; i_point++) {
    exc_energy += exc[i_point] * rho[i_point] * 
                  grid_weights[i_point] * becke_weights[i_point];
  }
  
  // Build XC potential matrix
  for (int i_point = 0; i_point < n_points; i_point++) {
    double w = grid_weights[i_point] * becke_weights[i_point];
    
    // LDA contribution
    for (int i = 0; i < n_basis; i++) {
      for (int j = 0; j <= i; j++) {
        double val = vrho[i_point] * basis_values[i_point][i] * 
                    basis_values[i_point][j] * w;
        vxc_matrix(i, j) += val;
        if (i != j) vxc_matrix(j, i) += val;
      }
    }
    
    // GGA contribution
    if (xc_func->is_gga() && vsigma[i_point] != 0.0) {
      // Simplified - full implementation would include gradient contributions
      // This requires second derivatives of basis functions
    }
    
    // Meta-GGA contribution
    if (xc_func->is_meta() && vtau[i_point] != 0.0) {
      // Simplified - full implementation would include tau contributions
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridIntegrator::evaluate_basis_at_points(BasisSetManager *basis,
                                              std::vector<std::vector<double>> &basis_values,
                                              std::vector<std::vector<std::vector<double>>> &basis_gradients)
{
  int n_points = grid_points.size();
  int n_basis = basis->get_n_basis();
  int n_shells = basis->get_n_shells();
  
  // For each grid point, evaluate all basis functions
  for (int i_point = 0; i_point < n_points; i_point++) {
    int bf_idx = 0;
    
    for (int shell = 0; shell < n_shells; shell++) {
      int l = basis->get_angular_momentum(shell);
      auto exponents = basis->get_exponents(shell);
      auto coefficients = basis->get_coefficients(shell);
      
      // Number of functions in this shell
      int n_funcs = (l + 1) * (l + 2) / 2;
      
      // Evaluate contracted Gaussian
      for (int func = 0; func < n_funcs; func++) {
        double value = 0.0;
        
        // Contract over primitives
        for (size_t prim = 0; prim < exponents.size(); prim++) {
          double alpha = exponents[prim];
          double coeff = coefficients[prim];
          
          // Distance from basis function center (assuming at origin for now)
          double r2 = 0.0;
          for (int k = 0; k < 3; k++) {
            r2 += grid_points[i_point][k] * grid_points[i_point][k];
          }
          
          // Gaussian value
          double gauss = coeff * exp(-alpha * r2);
          
          // Add angular part (simplified - should use spherical harmonics)
          if (l == 0) {
            // s-orbital
            value += gauss;
          } else if (l == 1) {
            // p-orbitals (px, py, pz)
            if (func < 3) {
              value += gauss * grid_points[i_point][func];
            }
          } else if (l == 2) {
            // d-orbitals (simplified)
            value += gauss * pow(grid_points[i_point][func % 3], 2);
          }
        }
        
        basis_values[i_point][bf_idx] = value;
        
        // Compute gradients if needed
        if (!basis_gradients.empty()) {
          // Simplified gradient - full implementation would be more complex
          for (int k = 0; k < 3; k++) {
            basis_gradients[i_point][bf_idx][k] = 
              -2.0 * exponents[0] * grid_points[i_point][k] * value;
          }
        }
        
        bf_idx++;
      }
    }
  }
}
