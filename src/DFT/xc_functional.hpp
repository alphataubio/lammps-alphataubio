/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS DFT Package - XC Functional Wrapper
   
   Provides interface to LibXC functionals with special implementations
   for complex functionals like wB97M-V
------------------------------------------------------------------------- */

#ifndef LMP_DFT_XC_FUNCTIONAL_H
#define LMP_DFT_XC_FUNCTIONAL_H

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <xc.h>

namespace LAMMPS_NS {

// XC functional families
enum class FunctionalFamily {
  LDA,
  GGA,
  MGGA,
  HYB_LDA,
  HYB_GGA,
  HYB_MGGA,
  RANGE_SEPARATED
};

// Special functional types
enum class SpecialFunctional {
  NONE,
  WB97MV,    // wB97M-V
  WB97XV,    // wB97X-V
  B97MV,     // B97M-V
  SCANL,     // SCAN-L
  R2SCAN,    // r2SCAN
  R2SCANL    // r2SCAN-L
};

class XCFunctional {
 public:
  XCFunctional(const std::string &name);
  ~XCFunctional();
  
  // Main evaluation interface
  void evaluate(const std::vector<double> &rho,
               const std::vector<double> &sigma,  // |grad rho|^2
               const std::vector<double> &lapl,   // laplacian
               const std::vector<double> &tau,    // kinetic density
               std::vector<double> &exc,
               std::vector<double> &vrho,
               std::vector<double> &vsigma,
               std::vector<double> &vlapl,
               std::vector<double> &vtau);
  
  // Simplified interfaces for different functional types
  void evaluate_lda(const std::vector<double> &rho,
                   std::vector<double> &exc,
                   std::vector<double> &vrho);
  
  void evaluate_gga(const std::vector<double> &rho,
                   const std::vector<double> &sigma,
                   std::vector<double> &exc,
                   std::vector<double> &vrho,
                   std::vector<double> &vsigma);
  
  void evaluate_mgga(const std::vector<double> &rho,
                    const std::vector<double> &sigma,
                    const std::vector<double> &lapl,
                    const std::vector<double> &tau,
                    std::vector<double> &exc,
                    std::vector<double> &vrho,
                    std::vector<double> &vsigma,
                    std::vector<double> &vlapl,
                    std::vector<double> &vtau);
  
  // Spin-polarized evaluation
  void evaluate_spin(const std::vector<double> &rho_a,
                    const std::vector<double> &rho_b,
                    const std::vector<double> &sigma_aa,
                    const std::vector<double> &sigma_ab,
                    const std::vector<double> &sigma_bb,
                    std::vector<double> &exc,
                    std::vector<double> &vrho_a,
                    std::vector<double> &vrho_b,
                    std::vector<double> &vsigma_aa,
                    std::vector<double> &vsigma_ab,
                    std::vector<double> &vsigma_bb);
  
  // Query functional properties
  bool is_lda() const { return family == FunctionalFamily::LDA; }
  bool is_gga() const { return family == FunctionalFamily::GGA || 
                               family == FunctionalFamily::HYB_GGA; }
  bool is_meta() const { return family == FunctionalFamily::MGGA || 
                                family == FunctionalFamily::HYB_MGGA; }
  bool is_hybrid() const { return is_hybrid_functional; }
  bool is_range_separated() const { return is_rs_functional; }
  
  double get_exact_exchange_fraction() const { return exx_fraction; }
  double get_range_separation_omega() const { return omega; }
  double get_short_range_fraction() const { return sr_exx_fraction; }
  double get_long_range_fraction() const { return lr_exx_fraction; }
  
  // Special functional implementations
  void evaluate_wb97mv(const std::vector<double> &rho,
                      const std::vector<double> &sigma,
                      const std::vector<double> &tau,
                      std::vector<double> &exc,
                      std::vector<double> &vrho,
                      std::vector<double> &vsigma,
                      std::vector<double> &vtau);
  
  void evaluate_wb97xv(const std::vector<double> &rho,
                      const std::vector<double> &sigma,
                      const std::vector<double> &tau,
                      std::vector<double> &exc,
                      std::vector<double> &vrho,
                      std::vector<double> &vsigma,
                      std::vector<double> &vtau);
  
  void evaluate_b97mv(const std::vector<double> &rho,
                     const std::vector<double> &sigma,
                     const std::vector<double> &tau,
                     std::vector<double> &exc,
                     std::vector<double> &vrho,
                     std::vector<double> &vsigma,
                     std::vector<double> &vtau);
  
  // VV10 non-local correlation
  void compute_vv10_nlc(const std::vector<double> &rho,
                       const std::vector<double> &sigma,
                       double b_vv10, double C_vv10,
                       std::vector<double> &enlc,
                       std::vector<double> &vnlc_rho,
                       std::vector<double> &vnlc_sigma);
  
  // Get functional information
  std::string get_name() const { return functional_name; }
  std::string get_reference() const { return reference; }
  std::string get_doi() const { return doi; }
  
  // Set custom parameters (for development/testing)
  void set_parameter(const std::string &param_name, double value);
  double get_parameter(const std::string &param_name) const;
  
 private:
  // LibXC functional pointers
  xc_func_type *func_x;
  xc_func_type *func_c;
  xc_func_type *func_xc;
  
  // Functional properties
  FunctionalFamily family;
  SpecialFunctional special_type;
  bool is_gga_functional;
  bool is_meta_functional;
  bool is_hybrid_functional;
  bool is_rs_functional;
  bool needs_laplacian;
  
  // Hybrid/range-separated parameters
  double exx_fraction;      // Exact exchange mixing
  double sr_exx_fraction;   // Short-range exact exchange
  double lr_exx_fraction;   // Long-range exact exchange
  double omega;             // Range-separation parameter
  
  // Functional name and metadata
  std::string functional_name;
  std::string reference;
  std::string doi;
  
  // Custom parameters for special functionals
  std::map<std::string, double> custom_params;
  
  // Initialization
  void initialize_functional(const std::string &name);
  void parse_functional_string(const std::string &name);
  void setup_libxc_functional(int func_id);
  void identify_special_functional(const std::string &name);
  
  // Helper functions for B97-style functionals
  double b97_enhancement_factor(double s, double a0, double a1, double a2);
  void compute_b97_exchange(const std::vector<double> &rho,
                           const std::vector<double> &sigma,
                           const std::vector<double> &tau,
                           const double *params,
                           std::vector<double> &ex,
                           std::vector<double> &vx_rho,
                           std::vector<double> &vx_sigma,
                           std::vector<double> &vx_tau);
  
  void compute_b97_correlation(const std::vector<double> &rho,
                              const std::vector<double> &sigma,
                              const std::vector<double> &tau,
                              const double *params_ss,
                              const double *params_os,
                              std::vector<double> &ec,
                              std::vector<double> &vc_rho,
                              std::vector<double> &vc_sigma,
                              std::vector<double> &vc_tau);
  
  // Parameters for specific functionals
  struct WB97MVParams {
    // Exchange parameters
    double a0_x = 0.85;
    double a1_x = 1.007;
    double a2_x = 0.259;
    
    // Same-spin correlation
    double a0_c_ss = 1.0;
    double a1_c_ss = -3.382;
    double a2_c_ss = -0.892;
    
    // Opposite-spin correlation
    double a0_c_os = 1.0;
    double a1_c_os = -1.855;
    double a2_c_os = 0.653;
    
    // Range separation
    double omega = 0.3;
    double c_x_lr = 1.0;
    double c_x_sr = 0.15;
    
    // VV10 NLC
    double b_vv10 = 6.0;
    double C_vv10 = 0.01;
  } wb97mv_params;
  
  struct WB97XVParams {
    // Exchange parameters
    double a0_x = 0.804;
    double a1_x = 0.725;
    double a2_x = 0.212;
    
    // Correlation parameters
    double a0_c_ss = 1.0;
    double a1_c_ss = -2.544;
    double a2_c_ss = -0.133;
    
    double a0_c_os = 1.0;
    double a1_c_os = -4.528;
    double a2_c_os = 2.321;
    
    // Range separation
    double omega = 0.3;
    double c_x_lr = 1.0;
    double c_x_sr = 0.167;
    
    // VV10 NLC
    double b_vv10 = 6.0;
    double C_vv10 = 0.01;
  } wb97xv_params;
  
  struct B97MVParams {
    // Exchange parameters
    double a0_x = 1.0;
    double a1_x = 0.933;
    double a2_x = -0.209;
    
    // Correlation parameters
    double a0_c_ss = 1.0;
    double a1_c_ss = -2.270;
    double a2_c_ss = 0.827;
    
    double a0_c_os = 1.0;
    double a1_c_os = -5.572;
    double a2_c_os = 6.378;
    
    // VV10 NLC
    double b_vv10 = 11.0;
    double C_vv10 = 0.01;
  } b97mv_params;
};

// Factory function for creating XC functionals
std::unique_ptr<XCFunctional> create_xc_functional(const std::string &name);

// Registry of available functionals
class XCFunctionalRegistry {
 public:
  static XCFunctionalRegistry& instance();
  
  void register_functional(const std::string &name, 
                          const std::string &description,
                          FunctionalFamily family,
                          bool is_hybrid = false,
                          double exx_fraction = 0.0);
  
  bool is_available(const std::string &name) const;
  std::vector<std::string> get_available_functionals() const;
  std::string get_description(const std::string &name) const;
  
 private:
  XCFunctionalRegistry();
  
  struct FunctionalInfo {
    std::string description;
    FunctionalFamily family;
    bool is_hybrid;
    double exx_fraction;
  };
  
  std::map<std::string, FunctionalInfo> functionals;
  
  void initialize_registry();
};

// Helper functions for XC evaluation
namespace XCHelpers {
  
  // Compute reduced density gradient
  void compute_reduced_gradient(const std::vector<double> &rho,
                               const std::vector<double> &sigma,
                               std::vector<double> &s);
  
  // Compute dimensionless density gradient for meta-GGA
  void compute_alpha(const std::vector<double> &rho,
                    const std::vector<double> &sigma,
                    const std::vector<double> &tau,
                    std::vector<double> &alpha);
  
  // Compute z = tau_W / tau for meta-GGA
  void compute_z(const std::vector<double> &rho,
                const std::vector<double> &sigma,
                const std::vector<double> &tau,
                std::vector<double> &z);
  
  // LDA exchange and correlation
  double lda_exchange(double rho);
  double lda_correlation_vwn(double rho);
  double lda_correlation_pw92(double rho);
  
  // PBE enhancement factor
  double pbe_enhancement(double s, double kappa, double mu);
  
  // B88 exchange enhancement
  double b88_enhancement(double s, double beta);
  
  // LYP correlation
  double lyp_correlation(double rho_a, double rho_b, 
                        double sigma_aa, double sigma_ab, double sigma_bb);
}

}  // namespace LAMMPS_NS

#endif
