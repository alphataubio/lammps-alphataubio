/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS DFT Package - Basis Set Manager
   
   Handles loading and managing Gaussian basis sets from JSON files
   compatible with the Basis Set Exchange format
------------------------------------------------------------------------- */

#ifndef LMP_DFT_BASIS_MANAGER_H
#define LMP_DFT_BASIS_MANAGER_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>

namespace LAMMPS_NS {

using json = nlohmann_lmp::json;

// Structure to hold basis function information
struct BasisShell {
  int atom_index;           // Which atom this shell belongs to
  int angular_momentum;     // L value (0=s, 1=p, 2=d, 3=f, etc.)
  std::vector<double> exponents;
  std::vector<double> coefficients;
  std::vector<double> normalized_coefficients;
  double center[3];         // Position of the atom
};

// Structure for a single Gaussian primitive
struct GaussianPrimitive {
  double exponent;
  double coefficient;
  double norm_const;        // Normalization constant
};

// Structure for a contracted Gaussian function
struct ContractedGaussian {
  int l, m, n;              // Cartesian angular momentum components
  std::vector<GaussianPrimitive> primitives;
  double center[3];
  
  double evaluate(double x, double y, double z) const;
  void evaluate_gradient(double x, double y, double z, double grad[3]) const;
};

class BasisSetManager {
 public:
  BasisSetManager();
  ~BasisSetManager();
  
  // Load basis set from JSON file (BSE format)
  void load_from_json(const std::string &filename);
  
  // Load basis set from internal library
  void load_from_library(const std::string &basis_name, 
                         const std::vector<int> &atomic_numbers);
  
  // Access basis information
  int get_n_basis() const { return n_basis_functions; }
  int get_n_shells() const { return n_shells; }
  int get_n_primitives() const { return n_primitives; }
  
  // Get shell information
  const BasisShell& get_shell(int idx) const { return shells[idx]; }
  std::vector<BasisShell> get_shells() const { return shells; }
  
  // Get basis functions for a specific atom
  std::vector<int> get_atom_shell_indices(int atom_idx) const;
  int get_atom_n_basis(int atom_idx) const;
  
  // Evaluate basis functions at a point
  void evaluate_basis(double x, double y, double z,
                      std::vector<double> &values) const;
  void evaluate_basis_gradient(double x, double y, double z,
                               std::vector<double> &values,
                               std::vector<std::vector<double>> &gradients) const;
  
  // Get Cartesian Gaussian functions
  std::vector<ContractedGaussian> get_cartesian_gaussians() const;
  
  // Normalization
  void normalize_basis_functions();
  static double compute_normalization(int l, double exponent);
  static double compute_overlap_integral(const GaussianPrimitive &a,
                                        const GaussianPrimitive &b,
                                        int la, int lb);
  
  // Angular momentum utilities
  static int n_cartesian(int l);  // Number of Cartesian functions for angular momentum l
  static int n_spherical(int l);  // Number of spherical functions for angular momentum l
  static void get_cartesian_components(int l, std::vector<std::vector<int>> &components);
  
  // Element information
  void set_element_info(int atomic_number, const std::string &symbol);
  std::string get_element_symbol(int atomic_number) const;
  
 private:
  int n_basis_functions;
  int n_shells;
  int n_primitives;
  int n_atoms;
  
  std::vector<BasisShell> shells;
  std::vector<int> shell_to_atom;
  std::vector<int> shell_to_basis_function;
  std::vector<ContractedGaussian> contracted_gaussians;
  
  // Element information
  std::map<int, std::string> element_symbols;
  std::map<std::string, int> element_numbers;
  
  // JSON parsing helpers
  void parse_bse_json(const json &j);
  void parse_element_basis(const json &element_data, int atomic_number);
  void parse_electron_shell(const json &shell_data, int atom_idx, 
                           double center[3]);
  
  // Convert between Cartesian and spherical harmonics
  void convert_to_cartesian();
  
  // Basis set library (common basis sets)
  void load_sto3g(const std::vector<int> &atomic_numbers);
  void load_321g(const std::vector<int> &atomic_numbers);
  void load_631g(const std::vector<int> &atomic_numbers);
  void load_631gss(const std::vector<int> &atomic_numbers);
  void load_def2svp(const std::vector<int> &atomic_numbers);
  void load_def2tzvp(const std::vector<int> &atomic_numbers);
  void load_ccpvdz(const std::vector<int> &atomic_numbers);
  void load_ccpvtz(const std::vector<int> &atomic_numbers);
  void load_augccpvdz(const std::vector<int> &atomic_numbers);
  void load_augccpvtz(const std::vector<int> &atomic_numbers);
  
  // Initialize standard element data
  void initialize_elements();
};

// Helper functions for Gaussian integrals
namespace GaussianIntegrals {
  
  // Boys function for electron-nuclear integrals
  double boys_function(int n, double x);
  
  // Overlap integral between two Gaussian primitives
  double overlap(double alpha_a, const double *Ra, int la, int ma, int na,
                double alpha_b, const double *Rb, int lb, int mb, int nb);
  
  // Kinetic energy integral
  double kinetic(double alpha_a, const double *Ra, int la, int ma, int na,
                double alpha_b, const double *Rb, int lb, int mb, int nb);
  
  // Nuclear attraction integral
  double nuclear(double alpha_a, const double *Ra, int la, int ma, int na,
                double alpha_b, const double *Rb, int lb, int mb, int nb,
                const double *Rc, double Z);
  
  // Electron repulsion integral (two-electron)
  double eri(double alpha_a, const double *Ra, int la, int ma, int na,
            double alpha_b, const double *Rb, int lb, int mb, int nb,
            double alpha_c, const double *Rc, int lc, int mc, int nc,
            double alpha_d, const double *Rd, int ld, int md, int nd);
  
  // Helper functions
  double gaussian_product_center(double alpha_a, double xa,
                                 double alpha_b, double xb);
  double gaussian_product_exponent(double alpha_a, double alpha_b);
  double gaussian_product_prefactor(double alpha_a, const double *Ra,
                                    double alpha_b, const double *Rb);
  
  // Hermite integrals
  double hermite_integral(int i, int j, int t, double Qx, 
                         double alpha_a, double alpha_b);
  
  // Factorial and double factorial
  int factorial(int n);
  int double_factorial(int n);
  
  // Binomial coefficient
  int binomial(int n, int k);
}

}  // namespace LAMMPS_NS

#endif
