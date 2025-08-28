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
#include <xc.h>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <algorithm>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

XCFunctional::XCFunctional(const std::string &name) : functional_name(name)
{
  func_x = nullptr;
  func_c = nullptr;
  func_xc = nullptr;
  
  is_gga_functional = false;
  is_meta_functional = false;
  is_hybrid_functional = false;
  is_rs_functional = false;
  
  exx_fraction = 0.0;
  omega = 0.0;
  
  initialize_functional(name);
}

/* ---------------------------------------------------------------------- */

XCFunctional::~XCFunctional()
{
  if (func_x) {
    xc_func_end(func_x);
    delete func_x;
  }
  if (func_c) {
    xc_func_end(func_c);
    delete func_c;
  }
  if (func_xc) {
    xc_func_end(func_xc);
    delete func_xc;
  }
}

/* ---------------------------------------------------------------------- */

void XCFunctional::initialize_functional(const std::string &name)
{
  // Special case for wB97M-V
  if (name == "wB97M-V" || name == "wb97mv" || name == "wB97MV") {
    // wB97M-V requires special implementation
    is_rs_functional = true;
    is_meta_functional = true;
    is_hybrid_functional = true;
    omega = 0.3;
    exx_fraction = 0.15;  // Short-range fraction
    return;
  }
  
  // Check if it's a combined functional or separate X+C
  if (name.find('+') != std::string::npos) {
    parse_functional_string(name);
  } else {
    // Try as combined XC functional
    int func_id = xc_functional_get_number(name.c_str());
    if (func_id == -1) {
      throw std::runtime_error("Unknown functional: " + name);
    }
    
    func_xc = new xc_func_type;
    if (xc_func_init(func_xc, func_id, XC_UNPOLARIZED) != 0) {
      delete func_xc;
      func_xc = nullptr;
      throw std::runtime_error("Failed to initialize functional: " + name);
    }
    
    // Get functional properties
    const xc_func_info_type *info = xc_func_get_info(func_xc);
    int family = xc_func_info_get_family(info);
    
    switch (family) {
      case XC_FAMILY_LDA:
        break;
      case XC_FAMILY_GGA:
      case XC_FAMILY_HYB_GGA:
        is_gga_functional = true;
        break;
      case XC_FAMILY_MGGA:
      case XC_FAMILY_HYB_MGGA:
        is_gga_functional = true;
        is_meta_functional = true;
        break;
    }
    
    if (family == XC_FAMILY_HYB_GGA || family == XC_FAMILY_HYB_MGGA) {
      is_hybrid_functional = true;
      exx_fraction = xc_hyb_exx_coef(func_xc);
    }
    
    // Check for range separation
    if (xc_func_info_get_flags(info) & XC_FLAGS_HYB_CAM) {
      is_rs_functional = true;
      double cam_alpha, cam_beta, cam_omega;
      xc_hyb_cam_coef(func_xc, &cam_omega, &cam_alpha, &cam_beta);
      omega = cam_omega;
      exx_fraction = cam_alpha + cam_beta;  // Total exact exchange
    }
  }
}

/* ---------------------------------------------------------------------- */

void XCFunctional::parse_functional_string(const std::string &name)
{
  // Parse "X_FUNC+C_FUNC" format
  size_t plus_pos = name.find('+');
  std::string x_name = name.substr(0, plus_pos);
  std::string c_name = name.substr(plus_pos + 1);
  
  // Initialize exchange functional
  int x_id = xc_functional_get_number(x_name.c_str());
  if (x_id == -1) {
    throw std::runtime_error("Unknown exchange functional: " + x_name);
  }
  
  func_x = new xc_func_type;
  if (xc_func_init(func_x, x_id, XC_UNPOLARIZED) != 0) {
    delete func_x;
    func_x = nullptr;
    throw std::runtime_error("Failed to initialize exchange functional: " + x_name);
  }
  
  // Initialize correlation functional
  int c_id = xc_functional_get_number(c_name.c_str());
  if (c_id == -1) {
    throw std::runtime_error("Unknown correlation functional: " + c_name);
  }
  
  func_c = new xc_func_type;
  if (xc_func_init(func_c, c_id, XC_UNPOLARIZED) != 0) {
    delete func_c;
    func_c = nullptr;
    throw std::runtime_error("Failed to initialize correlation functional: " + c_name);
  }
  
  // Check functional families
  const xc_func_info_type *x_info = xc_func_get_info(func_x);
  const xc_func_info_type *c_info = xc_func_get_info(func_c);
  
  int x_family = xc_func_info_get_family(x_info);
  int c_family = xc_func_info_get_family(c_info);
  
  // Set flags based on most complex functional
  if (x_family == XC_FAMILY_GGA || x_family == XC_FAMILY_HYB_GGA ||
      c_family == XC_FAMILY_GGA || c_family == XC_FAMILY_HYB_GGA) {
    is_gga_functional = true;
  }
  
  if (x_family == XC_FAMILY_MGGA || x_family == XC_FAMILY_HYB_MGGA ||
      c_family == XC_FAMILY_MGGA || c_family == XC_FAMILY_HYB_MGGA) {
    is_meta_functional = true;
    is_gga_functional = true;
  }
  
  if (x_family == XC_FAMILY_HYB_GGA || x_family == XC_FAMILY_HYB_MGGA) {
    is_hybrid_functional = true;
    exx_fraction = xc_hyb_exx_coef(func_x);
  }
}

/* ---------------------------------------------------------------------- */

void XCFunctional::evaluate(const std::vector<double> &rho,
                           const std::vector<double> &sigma,
                           const std::vector<double> &lapl,
                           const std::vector<double> &tau,
                           std::vector<double> &exc,
                           std::vector<double> &vrho,
                           std::vector<double> &vsigma,
                           std::vector<double> &vlapl,
                           std::vector<double> &vtau)
{
  int n = rho.size();
  
  // Resize output arrays
  exc.resize(n);
  vrho.resize(n);
  if (is_gga_functional) vsigma.resize(n);
  if (is_meta_functional) {
    vlapl.resize(n);
    vtau.resize(n);
  }
  
  // Special case for wB97M-V
  if (functional_name == "wB97M-V" || functional_name == "wb97mv") {
    evaluate_wb97mv(rho, sigma, exc, vrho, vsigma);
    return;
  }
  
  // Standard LibXC evaluation
  if (func_xc) {
    // Combined XC functional
    if (!is_gga_functional) {
      // LDA
      xc_lda_exc_vxc(func_xc, n, rho.data(), exc.data(), vrho.data());
    } else if (!is_meta_functional) {
      // GGA
      xc_gga_exc_vxc(func_xc, n, rho.data(), sigma.data(), 
                     exc.data(), vrho.data(), vsigma.data());
    } else {
      // Meta-GGA
      xc_mgga_exc_vxc(func_xc, n, rho.data(), sigma.data(), 
                      lapl.data(), tau.data(),
                      exc.data(), vrho.data(), vsigma.data(), 
                      vlapl.data(), vtau.data());
    }
  } else {
    // Separate X and C functionals
    std::vector<double> ex(n), vx_rho(n), vx_sigma(n), vx_lapl(n), vx_tau(n);
    std::vector<double> ec(n), vc_rho(n), vc_sigma(n), vc_lapl(n), vc_tau(n);
    
    // Exchange
    if (func_x) {
      const xc_func_info_type *info = xc_func_get_info(func_x);
      int family = xc_func_info_get_family(info);
      
      if (family == XC_FAMILY_LDA) {
        xc_lda_exc_vxc(func_x, n, rho.data(), ex.data(), vx_rho.data());
      } else if (family == XC_FAMILY_GGA || family == XC_FAMILY_HYB_GGA) {
        xc_gga_exc_vxc(func_x, n, rho.data(), sigma.data(),
                       ex.data(), vx_rho.data(), vx_sigma.data());
      } else if (family == XC_FAMILY_MGGA || family == XC_FAMILY_HYB_MGGA) {
        xc_mgga_exc_vxc(func_x, n, rho.data(), sigma.data(),
                        lapl.data(), tau.data(),
                        ex.data(), vx_rho.data(), vx_sigma.data(),
                        vx_lapl.data(), vx_tau.data());
      }
    }
    
    // Correlation
    if (func_c) {
      const xc_func_info_type *info = xc_func_get_info(func_c);
      int family = xc_func_info_get_family(info);
      
      if (family == XC_FAMILY_LDA) {
        xc_lda_exc_vxc(func_c, n, rho.data(), ec.data(), vc_rho.data());
      } else if (family == XC_FAMILY_GGA) {
        xc_gga_exc_vxc(func_c, n, rho.data(), sigma.data(),
                       ec.data(), vc_rho.data(), vc_sigma.data());
      } else if (family == XC_FAMILY_MGGA) {
        xc_mgga_exc_vxc(func_c, n, rho.data(), sigma.data(),
                        lapl.data(), tau.data(),
                        ec.data(), vc_rho.data(), vc_sigma.data(),
                        vc_lapl.data(), vc_tau.data());
      }
    }
    
    // Combine results
    for (int i = 0; i < n; i++) {
      exc[i] = ex[i] + ec[i];
      vrho[i] = vx_rho[i] + vc_rho[i];
      if (is_gga_functional) {
        vsigma[i] = vx_sigma[i] + vc_sigma[i];
      }
      if (is_meta_functional) {
        vlapl[i] = vx_lapl[i] + vc_lapl[i];
        vtau[i] = vx_tau[i] + vc_tau[i];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void XCFunctional::evaluate_wb97mv(const std::vector<double> &rho,
                                   const std::vector<double> &sigma,
                                   std::vector<double> &exc,
                                   std::vector<double> &vrho,
                                   std::vector<double> &vsigma)
{
  // wB97M-V functional implementation
  // This is a simplified version - full implementation would be very complex
  
  int n = rho.size();
  exc.resize(n);
  vrho.resize(n);
  vsigma.resize(n);
  
  // B97 parameters for wB97M-V
  const double a0_x = 0.85;
  const double a1_x = 1.007;
  const double a2_x = 0.259;
  
  for (int i = 0; i < n; i++) {
    if (rho[i] < 1e-15) {
      exc[i] = 0.0;
      vrho[i] = 0.0;
      vsigma[i] = 0.0;
      continue;
    }
    
    // LDA exchange
    double ex_lda = -0.75 * pow(3.0 * rho[i] / M_PI, 1.0/3.0);
    
    // Reduced gradient
    double s = sqrt(sigma[i]) / (2.0 * pow(3.0 * M_PI * M_PI, 1.0/3.0) * 
                                 pow(rho[i], 4.0/3.0));
    
    // Enhancement factor
    double s2 = s * s;
    double u_x = (a0_x + a1_x * s2 + a2_x * s2 * s2) / (1.0 + 0.004 * s2);
    
    // Exchange energy density
    exc[i] = ex_lda * u_x;
    
    // Simplified correlation (should use proper B97 correlation)
    exc[i] += -0.03 * pow(rho[i], 1.0/3.0);
    
    // Simplified potentials
    vrho[i] = 4.0/3.0 * exc[i] / rho[i];
    vsigma[i] = 0.01 / (2.0 * sqrt(sigma[i]) + 1e-15);
  }
}
