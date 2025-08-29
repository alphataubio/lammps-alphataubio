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
#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include <xc.h>
#include <libint2.hpp>
#include "fmt/format.h"

namespace LAMMPS_NS {

class PairDFT : public Pair {
 public:
  PairDFT(class LAMMPS *);
  ~PairDFT() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void allocate();

 protected:
  // Core parameters
  double cut_global;
  std::string functional_name;
  std::string basis_file;
  
  // LibXC functional handles
  xc_func_type *xc_func_x;
  xc_func_type *xc_func_c;
  xc_func_type *xc_func_xc;
  
  int xc_functional_x;
  int xc_functional_c;
  int xc_functional_xc;
  
  bool use_combined_xc;
  bool is_hybrid;
  bool is_meta_gga;
  bool is_range_separated;
  
  double hybrid_coeff;
  double range_separation_param;
  
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
  double xc_energy;
  
  // Grid parameters
  int grid_size;
  std::string grid_type;
  
  // Grid data
  std::vector<std::vector<double>> grid_points;
  std::vector<double> grid_weights;
  std::vector<double> becke_weights;
  
  // Matrices
  Eigen::MatrixXd overlap_matrix;
  Eigen::MatrixXd kinetic_matrix;
  Eigen::MatrixXd nuclear_matrix;
  Eigen::MatrixXd coulomb_matrix;
  Eigen::MatrixXd exchange_matrix;
  Eigen::MatrixXd fock_matrix;
  Eigen::MatrixXd density_matrix;
  Eigen::MatrixXd density_matrix_prev;
  Eigen::MatrixXd mo_coefficients;
  Eigen::VectorXd mo_energies;
  
  // Basis set data
  int n_basis_functions;
  int n_shells;
  int n_occupied_orbitals;
  int n_electrons;
  
  std::vector<int> shell_to_atom;
  std::vector<int> angular_momentum;
  std::vector<std::vector<double>> exponents;
  std::vector<std::vector<double>> coefficients;
  std::vector<std::vector<double>> normalized_coefficients;
  //std::vector<std::vector<double>> atom_positions;
  
  // LibInt2 data
  std::unique_ptr<libint2::BasisSet> libint_basis;
  std::vector<std::unique_ptr<libint2::Engine>> engines;
  
  // =============================================
  // Basis set methods (pair_dft_basis.hpp)
  // =============================================
  void initialize_basis_set();
  void load_basis_from_json(const std::string &filename);
  void normalize_basis_functions();
  double compute_normalization(int l, double exponent);
  void set_atom_positions(const std::vector<std::vector<double>> &positions);
  int get_angular_momentum(int shell) const;
  std::vector<double> get_exponents(int shell) const;
  std::vector<double> get_coefficients(int shell) const;
  void evaluate_basis_at_points(const std::vector<std::vector<double>> &points,
                                std::vector<std::vector<double>> &basis_values,
                                std::vector<std::vector<std::vector<double>>> &basis_gradients);
  
  // =============================================
  // LibXC methods (pair_dft_libxc.hpp)
  // =============================================
  void parse_functional_name(const char *name);
  void initialize_libxc_functional();
  void cleanup_libxc();
  void evaluate_xc_functional(const std::vector<double> &rho,
                              const std::vector<double> &sigma,
                              const std::vector<double> &lapl,
                              const std::vector<double> &tau,
                              std::vector<double> &exc,
                              std::vector<double> &vrho,
                              std::vector<double> &vsigma,
                              std::vector<double> &vlapl,
                              std::vector<double> &vtau);
  
  // =============================================
  // LibInt2 methods (pair_dft_libint2.hpp)
  // =============================================
  void initialize_libint();
  void cleanup_libint();
  void compute_overlap_integrals();
  void compute_kinetic_integrals();
  void compute_nuclear_integrals();
  void compute_eri_with_density(const Eigen::MatrixXd &D,
                                Eigen::MatrixXd &J,
                                Eigen::MatrixXd &K);
  void compute_one_electron_integrals();
  void compute_overlap_gradient(std::vector<Eigen::MatrixXd> &dS);
  void compute_nuclear_gradient(std::vector<Eigen::MatrixXd> &dV);
  
  // =============================================
  // SCF methods (pair_dft_scf.hpp)
  // =============================================
  void perform_scf();
  void build_fock_matrix();
  void solve_roothaan_hall();
  void update_density_matrix();
  void mix_density_matrices(double mixing_param);
  double compute_density_change();
  void compute_energy();
  double compute_nuclear_repulsion_energy();
  void initialize_density_guess();
  
  // Grid integration
  void generate_integration_grid();
  void generate_lebedev_grid(int n_points, 
                             std::vector<std::vector<double>> &points,
                             std::vector<double> &weights);
  void generate_radial_grid(int n_points, double Z,
                           std::vector<double> &r, std::vector<double> &w);
  void compute_becke_weights(const std::vector<std::vector<double>> &atoms);
  void integrate_xc_on_grid(double &exc_energy, Eigen::MatrixXd &vxc_matrix);
  
  // Force methods
  void compute_hellmann_feynman_forces();
  void compute_pulay_forces();
  
  // Output methods
  void print_scf_header();
  void print_scf_iteration();
  void print_scf_summary();
};

}    // namespace LAMMPS_NS

#endif
#endif
