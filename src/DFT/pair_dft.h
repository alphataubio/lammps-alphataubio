/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(dft,PairDFT);
// clang-format on
#else

#ifndef LMP_PAIR_DFT_H
#define LMP_PAIR_DFT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairDFT : public Pair {
 public:
  PairDFT(class LAMMPS *);
  ~PairDFT() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;

 protected:

};


 
 // *** FIXME EVERYTHING BELOW THIS NEEDS TO BE REFACTORED TO THE PROPER .HPP ***
 
/*
 #include <vector>
#include <map>
#include <memory>
#include <string>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


// Forward declarations
namespace libint2 {
  class BasisSet;
  class Engine;
}


  // DFT-specific public methods
  void set_convergence_criteria(double energy_tol, double density_tol, int max_iter);
  double get_total_energy() const { return total_dft_energy; }
  double get_xc_energy() const { return xc_energy; }
  double get_kinetic_energy() const { return kinetic_energy; }
  double get_nuclear_repulsion() const { return nuclear_repulsion; }

  double cut_global;
  double **cut;
  double **offset;
  
  // DFT functional and basis set
  std::string functional_name;
  std::string basis_file;
  
  // Modular components
  std::unique_ptr<BasisSetManager> basis_manager;
  std::unique_ptr<IntegralEngine> integral_engine;
  std::unique_ptr<DensityMatrix> density_matrix;
  std::unique_ptr<XCFunctional> xc_functional;
  std::unique_ptr<GridIntegrator> grid_integrator;
  
  // LibXC functionals
  xc_func_type *xc_func_x;      // Exchange functional
  xc_func_type *xc_func_c;      // Correlation functional
  xc_func_type *xc_func_xc;     // Combined XC functional
  
  int xc_functional_x;           // Exchange functional ID
  int xc_functional_c;           // Correlation functional ID  
  int xc_functional_xc;          // Combined XC functional ID
  
  bool use_combined_xc;          // Use combined XC functional
  bool is_hybrid;                // Is this a hybrid functional
  bool is_meta_gga;              // Is this a meta-GGA functional
  bool is_range_separated;       // Is this a range-separated functional
  
  double hybrid_coeff;           // Exact exchange mixing coefficient
  double range_separation_param; // Range separation parameter (omega)
  
  // SCF parameters
  double energy_tolerance;
  double density_tolerance;
  int max_scf_iterations;
  int current_iteration;
  bool scf_converged;
  
  // Energy components
  double total_dft_energy;
  double kinetic_energy;
  double nuclear_repulsion;
  double electron_nuclear_energy;
  double electron_electron_energy;
  double xc_energy;
  double exact_exchange_energy;
  
  // Grid parameters for numerical integration
  int grid_size;
  double grid_tolerance;
  std::string grid_type;  // "Lebedev", "Becke", "MultiExp"
  
  // Density and potential arrays
  std::vector<double> electron_density;
  std::vector<double> density_gradient;
  std::vector<double> density_laplacian;
  std::vector<double> kinetic_density;  // For meta-GGA
  std::vector<double> xc_potential;
  
  // Matrices using Eigen for efficient linear algebra
  Eigen::MatrixXd overlap_matrix;
  Eigen::MatrixXd kinetic_matrix;
  Eigen::MatrixXd nuclear_matrix;
  Eigen::MatrixXd coulomb_matrix;
  Eigen::MatrixXd exchange_matrix;
  Eigen::MatrixXd fock_matrix;
  Eigen::MatrixXd density_matrix_eigen;
  Eigen::MatrixXd mo_coefficients;
  Eigen::VectorXd mo_energies;
  
  // Basis function information
  int n_basis_functions;
  int n_occupied_orbitals;
  int n_electrons;
  
  // Atom-specific parameters
  double **atomic_charges;       // Nuclear charges for each atom type
  double **vdw_radii;            // van der Waals radii for dispersion
    
  // Methods for DFT calculations
  void initialize_basis_set();
  void compute_one_electron_integrals();
  void compute_two_electron_integrals();
  void build_fock_matrix();
  void solve_roothaan_hall();
  void compute_density_matrix();
  void compute_energy();
  void perform_scf();
  
  // XC functional methods
  void evaluate_xc_functional();
  void compute_xc_potential();
  void compute_exact_exchange();
  
  // Grid integration methods
  void generate_integration_grid();
  void integrate_xc_on_grid();
  
  // Integral evaluation
  void compute_overlap_integrals();
  void compute_kinetic_integrals();
  void compute_nuclear_attraction_integrals();
  void compute_electron_repulsion_integrals();
  
  // Force calculation
  void compute_forces(int, int);
  void compute_pulay_forces();
  void compute_hellmann_feynman_forces();
  
  // Utility functions
  void allocate();
  void parse_functional_name(const char *);
  void load_basis_set(const char *);
  double compute_nuclear_repulsion_energy();
  void initialize_density_guess();
  bool check_convergence();
  void mix_density(double mixing_param);
  
  // Output and debugging
  void print_scf_header();
  void print_scf_iteration();
  void print_scf_summary();
  void write_molecular_orbitals(const char *filename);
  void write_density_cube(const char *filename);
  
};

// Density Matrix Manager
class DensityMatrix {
 public:
  DensityMatrix(int nbasis);
  ~DensityMatrix();
  
  void initialize_guess(const Eigen::MatrixXd &S);
  void update_from_mo(const Eigen::MatrixXd &C, int nocc);
  void mix_with_previous(double mixing_param);
  
  const Eigen::MatrixXd& get_current() const { return current; }
  const Eigen::MatrixXd& get_previous() const { return previous; }
  double get_change() const;
    
 private:
  int n_basis;
  Eigen::MatrixXd current;
  Eigen::MatrixXd previous;
  Eigen::MatrixXd difference;
  
  // DIIS acceleration
  bool use_diis;
  int diis_size;
  std::vector<Eigen::MatrixXd> diis_fock;
  std::vector<Eigen::MatrixXd> diis_error;
};

// Grid Integration Class
class GridIntegrator {
 public:
  GridIntegrator(int grid_size, const std::string &grid_type);
  ~GridIntegrator();
  
  void generate_grid(const std::vector<std::vector<double>> &atom_positions,
                     const std::vector<int> &atomic_numbers);
  
  void integrate_xc(XCFunctional *xc_func,
                    const Eigen::MatrixXd &density_matrix,
                    BasisSetManager *basis,
                    double &exc_energy,
                    Eigen::MatrixXd &vxc_matrix);
  
  int get_n_points() const { return grid_points.size(); }
  
 private:
  int target_grid_size;
  std::string grid_type;
  
  std::vector<std::vector<double>> grid_points;
  std::vector<double> grid_weights;
  
  // Becke partitioning
  std::vector<double> becke_weights;
  void compute_becke_weights(const std::vector<std::vector<double>> &atoms);
  
  // Radial and angular grids
  void generate_lebedev_grid(int n_points, std::vector<std::vector<double>> &points,
                             std::vector<double> &weights);
  void generate_radial_grid(int n_points, double Z,
                           std::vector<double> &r, std::vector<double> &w);
  
  // Basis function evaluation
  void evaluate_basis_at_points(BasisSetManager *basis,
                                std::vector<std::vector<double>> &basis_values,
                                std::vector<std::vector<std::vector<double>>> &basis_gradients);
};

*/

}    // namespace LAMMPS_NS

#endif
#endif
