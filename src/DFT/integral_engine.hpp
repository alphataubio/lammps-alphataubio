/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS DFT Package - Integral Engine
   
   Handles computation of molecular integrals using libint2
------------------------------------------------------------------------- */

#ifndef LMP_DFT_INTEGRAL_ENGINE_H
#define LMP_DFT_INTEGRAL_ENGINE_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <libint2.hpp>
#include "basis_manager.hpp"

namespace LAMMPS_NS {

class IntegralEngine {
 public:
  IntegralEngine(BasisSetManager *basis);
  ~IntegralEngine();
  
  // One-electron integrals
  void compute_overlap(Eigen::MatrixXd &S);
  void compute_kinetic(Eigen::MatrixXd &T);
  void compute_nuclear(Eigen::MatrixXd &V, 
                       const std::vector<double> &charges,
                       const std::vector<std::vector<double>> &positions);
  
  // Two-electron integrals
  void compute_eri(std::vector<double> &eri_tensor);
  void compute_eri_with_density(const Eigen::MatrixXd &D, 
                                Eigen::MatrixXd &J, 
                                Eigen::MatrixXd &K);
  
  // Schwarz screening for integral prescreening
  void compute_schwarz_ints();
  double get_schwarz_bound(int i, int j) const { 
    return schwarz_matrix(i, j); 
  }
  
  // Gradient integrals for forces
  void compute_overlap_gradient(std::vector<Eigen::MatrixXd> &dS);
  void compute_kinetic_gradient(std::vector<Eigen::MatrixXd> &dT);
  void compute_nuclear_gradient(std::vector<Eigen::MatrixXd> &dV);
  void compute_eri_gradient(std::vector<Eigen::MatrixXd> &dERI);
  
  // Range-separated integrals for hybrid functionals
  void compute_erf_eri(double omega, std::vector<double> &eri_lr);
  void compute_erfc_eri(double omega, std::vector<double> &eri_sr);
  
  // Multipole integrals
  void compute_dipole(std::vector<Eigen::MatrixXd> &mu);
  void compute_quadrupole(std::vector<Eigen::MatrixXd> &Q);
  
  // Integral derivatives for CPSCF
  void compute_overlap_hessian(std::vector<std::vector<Eigen::MatrixXd>> &d2S);
  void compute_eri_derivatives(const Eigen::MatrixXd &D,
                               std::vector<Eigen::MatrixXd> &dJ,
                               std::vector<Eigen::MatrixXd> &dK);
  
  // Utility functions
  void set_precision(double thresh) { precision_threshold = thresh; }
  double get_precision() const { return precision_threshold; }
  void enable_screening(bool screen) { use_screening = screen; }
  
  // Memory estimate
  size_t estimate_memory_usage() const;
  
 private:
  BasisSetManager *basis_set;
  std::unique_ptr<libint2::BasisSet> libint_basis;
  
  // Libint2 engines for different integral types
  std::vector<libint2::Engine> overlap_engines;
  std::vector<libint2::Engine> kinetic_engines;
  std::vector<libint2::Engine> nuclear_engines;
  std::vector<libint2::Engine> eri_engines;
  std::vector<libint2::Engine> eri_grad_engines;
  
  // Integral screening
  Eigen::MatrixXd schwarz_matrix;
  double precision_threshold;
  bool use_screening;
  
  // Basis set information cache
  int n_basis;
  int n_shells;
  std::vector<size_t> shell2bf;
  std::vector<size_t> bf2shell;
  
  // Initialize libint2
  void initialize_libint();
  void cleanup_libint();
  void setup_engines();
  
  // Helper functions for integral computation
  void compute_1body_ints(libint2::Engine &engine,
                          const libint2::BasisSet &bs,
                          Eigen::MatrixXd &result);
  
  void compute_2body_ints(libint2::Engine &engine,
                          const libint2::BasisSet &bs,
                          std::vector<double> &result);
  
  void compute_2body_fock(libint2::Engine &engine,
                          const libint2::BasisSet &bs,
                          const Eigen::MatrixXd &D,
                          Eigen::MatrixXd &F);
  
  // Shell pair data for efficient computation
  struct ShellPair {
    int shell1, shell2;
    double screen_value;
    bool significant;
  };
  std::vector<ShellPair> significant_shell_pairs;
  
  void identify_significant_shell_pairs();
  
  // Integral transformation
  void transform_integrals(const Eigen::MatrixXd &C,
                           const std::vector<double> &ao_ints,
                           std::vector<double> &mo_ints);
  
  // Density fitting / resolution of identity
  bool use_density_fitting;
  std::unique_ptr<libint2::BasisSet> aux_basis;
  void setup_density_fitting(const std::string &aux_basis_name);
  void compute_3center_ints(std::vector<double> &three_center);
  void compute_2center_ints(Eigen::MatrixXd &two_center);
};

// Helper class for managing integral batches
class IntegralBatch {
 public:
  IntegralBatch(int max_nprim, int max_l);
  ~IntegralBatch();
  
  void allocate(int n_shells);
  void compute(libint2::Engine &engine, 
               const libint2::Shell &s1,
               const libint2::Shell &s2);
  void compute(libint2::Engine &engine,
               const libint2::Shell &s1,
               const libint2::Shell &s2,
               const libint2::Shell &s3,
               const libint2::Shell &s4);
  
  const double* data() const { return buffer.data(); }
  size_t size() const { return buffer.size(); }
  
 private:
  std::vector<double> buffer;
  int max_nprim;
  int max_l;
  size_t buffer_size;
};

// Specialized integral engines for specific methods
class DFTIntegralEngine : public IntegralEngine {
 public:
  DFTIntegralEngine(BasisSetManager *basis);
  
  // DFT-specific integrals
  void compute_xc_matrix(const std::vector<double> &vxc_grid,
                        const std::vector<double> &weights,
                        const std::vector<std::vector<double>> &basis_values,
                        Eigen::MatrixXd &Vxc);
  
  void compute_grid_density(const Eigen::MatrixXd &D,
                           const std::vector<std::vector<double>> &basis_values,
                           std::vector<double> &rho);
  
  void compute_grid_gradient(const Eigen::MatrixXd &D,
                            const std::vector<std::vector<double>> &basis_values,
                            const std::vector<std::vector<std::vector<double>>> &basis_grads,
                            std::vector<std::vector<double>> &grad_rho);
  
  void compute_grid_laplacian(const Eigen::MatrixXd &D,
                              const std::vector<std::vector<double>> &basis_values,
                              const std::vector<std::vector<std::vector<double>>> &basis_grads,
                              const std::vector<std::vector<std::vector<std::vector<double>>>> &basis_hess,
                              std::vector<double> &lapl_rho);
  
  void compute_grid_tau(const Eigen::MatrixXd &D,
                       const std::vector<std::vector<std::vector<double>>> &basis_grads,
                       std::vector<double> &tau);
};

// Helper functions for special integrals
namespace SpecialIntegrals {
  
  // Ewald summation for periodic systems
  void compute_ewald_integrals(const libint2::BasisSet &basis,
                               const std::vector<double> &cell_vectors,
                               double alpha,
                               Eigen::MatrixXd &S_ewald,
                               Eigen::MatrixXd &T_ewald,
                               Eigen::MatrixXd &V_ewald);
  
  // Pseudopotential integrals
  void compute_pseudopotential_integrals(const libint2::BasisSet &basis,
                                         const std::vector<int> &atomic_numbers,
                                         const std::vector<std::vector<double>> &positions,
                                         Eigen::MatrixXd &V_pp);
  
  // Effective core potential integrals
  void compute_ecp_integrals(const libint2::BasisSet &basis,
                             const std::vector<int> &atomic_numbers,
                             const std::vector<std::vector<double>> &positions,
                             Eigen::MatrixXd &V_ecp);
  
  // Relativistic integrals (scalar relativistic)
  void compute_scalar_relativistic_integrals(const libint2::BasisSet &basis,
                                             Eigen::MatrixXd &T_rel,
                                             Eigen::MatrixXd &V_rel);
  
  // Spin-orbit coupling integrals
  void compute_spin_orbit_integrals(const libint2::BasisSet &basis,
                                    const std::vector<int> &atomic_numbers,
                                    const std::vector<std::vector<double>> &positions,
                                    std::vector<Eigen::MatrixXcd> &H_so);
}

}  // namespace LAMMPS_NS

#endif
