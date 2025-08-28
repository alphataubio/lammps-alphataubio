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

// LibXC functional IDs for common functionals
#ifndef XC_LDA_X
#define XC_LDA_X 1
#endif
#ifndef XC_LDA_C_PW
#define XC_LDA_C_PW 12
#endif
#ifndef XC_GGA_X_PBE
#define XC_GGA_X_PBE 101
#endif
#ifndef XC_GGA_C_PBE
#define XC_GGA_C_PBE 130
#endif
#ifndef XC_GGA_X_B88
#define XC_GGA_X_B88 106
#endif
#ifndef XC_GGA_C_LYP
#define XC_GGA_C_LYP 131
#endif
#ifndef XC_HYB_GGA_XC_B3LYP
#define XC_HYB_GGA_XC_B3LYP 402
#endif

/* ---------------------------------------------------------------------- 

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

*/

/* ----------------------------------------------------------------------
   parse functional name to LibXC ID
------------------------------------------------------------------------- */

void PairDFT::parse_functional_name(const char *name)
{

  utils::logmesg(lmp, "*** parse_functional_name({})\n", name);

  // Special handling for wB97M-V
  // HYB_MGGA_XC_WB97M_V (id=531): wB97M-V exchange-correlation functional
  // N. Mardirossian and M. Head-Gordon., J. Chem. Phys. 144, 214110 (2016) (doi: 10.1063/1.4952647)

  if (strcasecmp(name, "wB97M-V") == 0 || strcasecmp(name, "wb97mv") == 0) {
    // wB97M-V is not directly in LibXC, we need special implementation
    is_range_separated = true;
    is_meta_gga = true;
    is_hybrid = true;
    range_separation_param = 0.3;
    hybrid_coeff = 0.15;  // Short-range exact exchange
    use_combined_xc = false;
    
    // We'll use a special implementation for wB97M-V
    xc_functional_xc = -999;  // Special flag for wB97M-V
    
    if (comm->me == 0) {
      utils::logmesg(lmp, "Using special implementation for wB97M-V functional\n");
    }
    return;
  }
  
  // Direct mapping of common functionals to LibXC IDs
  std::string func_str(name);
  
  // Handle common functional names with direct ID mapping
  if (func_str == "PBE") {
    // PBE = PBE exchange + PBE correlation
    use_combined_xc = false;
    xc_functional_x = XC_GGA_X_PBE;  // ID 101
    xc_functional_c = XC_GGA_C_PBE;  // ID 130
    
    xc_func_x = new xc_func_type;
    xc_func_c = new xc_func_type;
    
    if (xc_func_init(xc_func_x, xc_functional_x, XC_UNPOLARIZED) != 0) {
      error->all(FLERR, "Failed to initialize PBE exchange functional");
    }
    if (xc_func_init(xc_func_c, xc_functional_c, XC_UNPOLARIZED) != 0) {
      error->all(FLERR, "Failed to initialize PBE correlation functional");
    }
    
    if (comm->me == 0) {
      utils::logmesg(lmp, "DFT: Using PBE functional (GGA)\n");
    }
    return;
  }
  else if (func_str == "B3LYP") {
    // B3LYP is a combined hybrid functional
    use_combined_xc = true;
    xc_functional_xc = XC_HYB_GGA_XC_B3LYP;  // ID 402
    
    xc_func_xc = new xc_func_type;
    if (xc_func_init(xc_func_xc, xc_functional_xc, XC_UNPOLARIZED) != 0) {
      error->all(FLERR, "Failed to initialize B3LYP functional");
    }
    
    // Check if it's hybrid
    is_hybrid = true;
    hybrid_coeff = xc_hyb_exx_coef(xc_func_xc);
    
    if (comm->me == 0) {
      utils::logmesg(lmp, fmt::format("DFT: Using B3LYP hybrid functional ({}% HF exchange)\n", 
                                      hybrid_coeff * 100));
    }
    return;
  }
  else if (func_str == "PBE0" || func_str == "PBEH") {
    // PBE0/PBEh is a hybrid functional
    use_combined_xc = true;
    xc_functional_xc = 406;  // XC_HYB_GGA_XC_PBEH
    
    xc_func_xc = new xc_func_type;
    if (xc_func_init(xc_func_xc, xc_functional_xc, XC_UNPOLARIZED) != 0) {
      error->all(FLERR, "Failed to initialize PBE0 functional");
    }
    
    // Check if it's hybrid
    is_hybrid = true;
    hybrid_coeff = xc_hyb_exx_coef(xc_func_xc);
    
    if (comm->me == 0) {
      utils::logmesg(lmp, fmt::format("DFT: Using PBE0 hybrid functional ({}% HF exchange)\n", 
                                      hybrid_coeff * 100));
    }
    return;
  }
  
  // If not a common name, try to parse as LibXC format
  // Check if it contains a "+" for separate functionals
  if (strchr(name, '+') != nullptr) {
    // Separate X and C functionals
    use_combined_xc = false;
    char *name_copy = strdup(name);
    char *x_func = strtok(name_copy, "+");
    char *c_func = strtok(nullptr, "+");
    
    if (x_func) {
      xc_functional_x = xc_functional_get_number(x_func);
      if (xc_functional_x == -1) {
        char errmsg[256];
        snprintf(errmsg, 256, "Unknown exchange functional: %s", x_func);
        error->all(FLERR, errmsg);
      }
      xc_func_x = new xc_func_type;
      if (xc_func_init(xc_func_x, xc_functional_x, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize exchange functional");
      }
    }
    
    if (c_func) {
      xc_functional_c = xc_functional_get_number(c_func);
      if (xc_functional_c == -1) {
        char errmsg[256];
        snprintf(errmsg, 256, "Unknown correlation functional: %s", c_func);
        error->all(FLERR, errmsg);
      }
      xc_func_c = new xc_func_type;
      if (xc_func_init(xc_func_c, xc_functional_c, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize correlation functional");
      }
    }
    
    free(name_copy);
  } else {
    // Try as a combined XC functional
    use_combined_xc = true;
    xc_functional_xc = xc_functional_get_number(name);
    if (xc_functional_xc == -1) {
      char errmsg[512];
      snprintf(errmsg, 512, "Unknown functional: %s\n"
               "Supported common names: PBE, LDA, B3LYP, BLYP, BP86, PBE0, TPSS, SCAN, wB97M-V\n"
               "Or use LibXC format: HYB_GGA_XC_B3LYP, GGA_X_PBE+GGA_C_PBE, etc.", name);
      error->all(FLERR, errmsg);
    }
    xc_func_xc = new xc_func_type;
    if (xc_func_init(xc_func_xc, xc_functional_xc, XC_UNPOLARIZED) != 0) {
      error->all(FLERR, "Failed to initialize XC functional");
    }
    
    // Check functional family
    const xc_func_info_type *info = xc_func_get_info(xc_func_xc);
    int family = xc_func_info_get_family(info);
    
    if (family == XC_FAMILY_HYB_GGA || family == XC_FAMILY_HYB_MGGA) {
      is_hybrid = true;
      hybrid_coeff = xc_hyb_exx_coef(xc_func_xc);
    }
    if (family == XC_FAMILY_MGGA || family == XC_FAMILY_HYB_MGGA) {
      is_meta_gga = true;
    }
  }
}


/* ---------------------------------------------------------------------- */

void PairDFT::initialize_functional(const std::string &name)
{
  
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

void PairDFT::parse_functional_string(const std::string &name)
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

void PairDFT::evaluate(const std::vector<double> &rho,
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
      if (is_gga_functional) vsigma[i] = vx_sigma[i] + vc_sigma[i];
      if (is_meta_functional) {
        vlapl[i] = vx_lapl[i] + vc_lapl[i];
        vtau[i] = vx_tau[i] + vc_tau[i];
      }
    }
  }
}
