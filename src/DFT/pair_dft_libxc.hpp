/* ----------------------------------------------------------------------
   LibXC methods for PairDFT
------------------------------------------------------------------------- */

#include <xc_funcs.h>

#include <cstring>
#include <stdexcept>
#include <algorithm>

/* ----------------------------------------------------------------------
   Parse functional name and initialize LibXC
------------------------------------------------------------------------- */

void PairDFT::parse_functional_name(const char *name)
{
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("Parsing functional: {}\n", name));
  }
  
  // Initialize to nullptr
  xc_func_x = nullptr;
  xc_func_c = nullptr;
  xc_func_xc = nullptr;
  
  xc_functional_x = -1;
  xc_functional_c = -1;
  xc_functional_xc = -1;
  
  use_combined_xc = false;
  is_hybrid = false;
  is_meta_gga = false;
  is_range_separated = false;
  
  hybrid_coeff = 0.0;
  range_separation_param = 0.0;
  
  // Direct mapping of common functionals to LibXC IDs
  std::string func_str(name);
  
  if (func_str == "PBE") {
    // PBE = PBE exchange + PBE correlation
    use_combined_xc = false;
    xc_functional_x = XC_GGA_X_PBE;  // ID 101
    xc_functional_c = XC_GGA_C_PBE;  // ID 130
  }
  else if (func_str == "LDA") {
    // LDA = Slater exchange + PW correlation
    use_combined_xc = false;
    xc_functional_x = XC_LDA_X;     // ID 1
    xc_functional_c = XC_LDA_C_PW;  // ID 12
  }
  else if (func_str == "B3LYP") {
    // B3LYP is a combined hybrid functional
    use_combined_xc = true;
    xc_functional_xc = XC_HYB_GGA_XC_B3LYP;  // ID 402
    is_hybrid = true;
  }
  else if (func_str == "PBE0" || func_str == "PBEH") {
    // PBE0/PBEh is a hybrid functional
    use_combined_xc = true;
    xc_functional_xc = XC_HYB_GGA_XC_PBEH;  // ID 406
    is_hybrid = true;
  }
  else if (func_str == "BLYP") {
    // BLYP = B88 exchange + LYP correlation  
    use_combined_xc = false;
    xc_functional_x = XC_GGA_X_B88;  // ID 106
    xc_functional_c = XC_GGA_C_LYP;  // ID 131
  }
  else if (func_str == "BP86") {
    // BP86 = B88 exchange + P86 correlation
    use_combined_xc = false;
    xc_functional_x = XC_GGA_X_B88;   // ID 106
    xc_functional_c = XC_GGA_C_P86;   // ID 132
  }
  else if (func_str == "TPSS") {
    // TPSS = TPSS exchange + TPSS correlation
    use_combined_xc = false;
    xc_functional_x = XC_MGGA_X_TPSS;  // ID 202
    xc_functional_c = XC_MGGA_C_TPSS;  // ID 231
    is_meta_gga = true;
  }
  else if (func_str == "SCAN") {
    // SCAN = SCAN exchange + SCAN correlation
    use_combined_xc = false;
    xc_functional_x = XC_MGGA_X_SCAN;  // ID 263
    xc_functional_c = XC_MGGA_C_SCAN;  // ID 267
    is_meta_gga = true;
  }
  else {
    // Try to parse as LibXC format
    if (strchr(name, '+') != nullptr) {
      // Separate X and C functionals
      use_combined_xc = false;
      char *name_copy = strdup(name);
      char *x_func = strtok(name_copy, "+");
      char *c_func = strtok(nullptr, "+");
      
      if (x_func) {
        xc_functional_x = xc_functional_get_number(x_func);
        if (xc_functional_x == -1) {
          error->all(FLERR, fmt::format("Unknown exchange functional: {}", x_func));
        }
      }
      
      if (c_func) {
        xc_functional_c = xc_functional_get_number(c_func);
        if (xc_functional_c == -1) {
          error->all(FLERR, fmt::format("Unknown correlation functional: {}", c_func));
        }
      }
      
      free(name_copy);
    } else {
      // Try as a combined XC functional
      use_combined_xc = true;
      xc_functional_xc = xc_functional_get_number(name);
      if (xc_functional_xc == -1) {
        error->all(FLERR, fmt::format("Unknown functional: {}\n"
                 "Supported: PBE, LDA, B3LYP, BLYP, BP86, PBE0, TPSS, SCAN\n"
                 "Or use LibXC format: HYB_GGA_XC_B3LYP, GGA_X_PBE+GGA_C_PBE", name));
      }
    }
  }
  
  // Initialize the functionals
  initialize_libxc_functional();
}

/* ----------------------------------------------------------------------
   Initialize LibXC functional(s)
------------------------------------------------------------------------- */

void PairDFT::initialize_libxc_functional()
{
  if (use_combined_xc) {
    // Combined XC functional
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
    
    if (comm->me == 0) {
      const char *func_name = xc_func_info_get_name(info);
      utils::logmesg(lmp, fmt::format("Using functional: {}\n", func_name));
      if (is_hybrid) {
        utils::logmesg(lmp, fmt::format("  Hybrid with {}% exact exchange\n", hybrid_coeff * 100));
      }
    }
  } else {
    // Separate exchange and correlation
    if (xc_functional_x >= 0) {
      xc_func_x = new xc_func_type;
      if (xc_func_init(xc_func_x, xc_functional_x, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize exchange functional");
      }
      
      // Check if exchange is hybrid
      const xc_func_info_type *info = xc_func_get_info(xc_func_x);
      int family = xc_func_info_get_family(info);
      
      if (family == XC_FAMILY_HYB_GGA || family == XC_FAMILY_HYB_MGGA) {
        is_hybrid = true;
        hybrid_coeff = xc_hyb_exx_coef(xc_func_x);
      }
      if (family == XC_FAMILY_MGGA || family == XC_FAMILY_HYB_MGGA) {
        is_meta_gga = true;
      }
    }
    
    if (xc_functional_c >= 0) {
      xc_func_c = new xc_func_type;
      if (xc_func_init(xc_func_c, xc_functional_c, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize correlation functional");
      }
      
      const xc_func_info_type *info = xc_func_get_info(xc_func_c);
      int family = xc_func_info_get_family(info);
      
      if (family == XC_FAMILY_MGGA) {
        is_meta_gga = true;
      }
    }
    
    if (comm->me == 0) {
      if (xc_func_x) {
        const xc_func_info_type *info = xc_func_get_info(xc_func_x);
        utils::logmesg(lmp, fmt::format("Exchange: {}\n", xc_func_info_get_name(info)));
      }
      if (xc_func_c) {
        const xc_func_info_type *info = xc_func_get_info(xc_func_c);
        utils::logmesg(lmp, fmt::format("Correlation: {}\n", xc_func_info_get_name(info)));
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Cleanup LibXC functionals
------------------------------------------------------------------------- */

void PairDFT::cleanup_libxc()
{
  if (xc_func_x) {
    xc_func_end(xc_func_x);
    delete xc_func_x;
    xc_func_x = nullptr;
  }
  if (xc_func_c) {
    xc_func_end(xc_func_c);
    delete xc_func_c;
    xc_func_c = nullptr;
  }
  if (xc_func_xc) {
    xc_func_end(xc_func_xc);
    delete xc_func_xc;
    xc_func_xc = nullptr;
  }
}

/* ----------------------------------------------------------------------
   Evaluate XC functional at grid points
------------------------------------------------------------------------- */

void PairDFT::evaluate_xc_functional(const std::vector<double> &rho,
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
  if (!sigma.empty()) vsigma.resize(n);
  if (!lapl.empty()) vlapl.resize(n);
  if (!tau.empty()) vtau.resize(n);
  
  // Standard LibXC evaluation
  if (xc_func_xc) {
    // Combined XC functional
    if (!is_meta_gga) {
      if (sigma.empty()) {
        // LDA
        xc_lda_exc_vxc(xc_func_xc, n, rho.data(), exc.data(), vrho.data());
      } else {
        // GGA
        xc_gga_exc_vxc(xc_func_xc, n, rho.data(), sigma.data(), 
                       exc.data(), vrho.data(), vsigma.data());
      }
    } else {
      // Meta-GGA
      xc_mgga_exc_vxc(xc_func_xc, n, rho.data(), sigma.data(), 
                      lapl.data(), tau.data(),
                      exc.data(), vrho.data(), vsigma.data(), 
                      vlapl.data(), vtau.data());
    }
  } else {
    // Separate X and C functionals
    std::vector<double> ex(n), vx_rho(n), vx_sigma(n), vx_lapl(n), vx_tau(n);
    std::vector<double> ec(n), vc_rho(n), vc_sigma(n), vc_lapl(n), vc_tau(n);
    
    // Exchange
    if (xc_func_x) {
      const xc_func_info_type *info = xc_func_get_info(xc_func_x);
      int family = xc_func_info_get_family(info);
      
      if (family == XC_FAMILY_LDA) {
        xc_lda_exc_vxc(xc_func_x, n, rho.data(), ex.data(), vx_rho.data());
      } else if (family == XC_FAMILY_GGA || family == XC_FAMILY_HYB_GGA) {
        xc_gga_exc_vxc(xc_func_x, n, rho.data(), sigma.data(),
                       ex.data(), vx_rho.data(), vx_sigma.data());
      } else if (family == XC_FAMILY_MGGA || family == XC_FAMILY_HYB_MGGA) {
        xc_mgga_exc_vxc(xc_func_x, n, rho.data(), sigma.data(),
                        lapl.data(), tau.data(),
                        ex.data(), vx_rho.data(), vx_sigma.data(),
                        vx_lapl.data(), vx_tau.data());
      }
    }
    
    // Correlation
    if (xc_func_c) {
      const xc_func_info_type *info = xc_func_get_info(xc_func_c);
      int family = xc_func_info_get_family(info);
      
      if (family == XC_FAMILY_LDA) {
        xc_lda_exc_vxc(xc_func_c, n, rho.data(), ec.data(), vc_rho.data());
      } else if (family == XC_FAMILY_GGA) {
        xc_gga_exc_vxc(xc_func_c, n, rho.data(), sigma.data(),
                       ec.data(), vc_rho.data(), vc_sigma.data());
      } else if (family == XC_FAMILY_MGGA) {
        xc_mgga_exc_vxc(xc_func_c, n, rho.data(), sigma.data(),
                        lapl.data(), tau.data(),
                        ec.data(), vc_rho.data(), vc_sigma.data(),
                        vc_lapl.data(), vc_tau.data());
      }
    }
    
    // Combine results
    for (int i = 0; i < n; i++) {
      exc[i] = ex[i] + ec[i];
      vrho[i] = vx_rho[i] + vc_rho[i];
      if (!vsigma.empty()) {
        vsigma[i] = vx_sigma[i] + vc_sigma[i];
      }
      if (!vlapl.empty()) {
        vlapl[i] = vx_lapl[i] + vc_lapl[i];
      }
      if (!vtau.empty()) {
        vtau[i] = vx_tau[i] + vc_tau[i];
      }
    }
  }
}
