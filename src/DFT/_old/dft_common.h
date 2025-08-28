/* ----------------------------------------------------------------------
   Common includes and definitions for DFT package
------------------------------------------------------------------------- */

#ifndef LMP_DFT_COMMON_H
#define LMP_DFT_COMMON_H

// Standard library
#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

// Eigen for linear algebra
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// LibXC for exchange-correlation functionals
#include <xc.h>

// Libint2 for integrals
#include <libint2.hpp>
#include <libint2/basis.h>
#include <libint2/shell.h>
#include <libint2/engine.h>

// JSON for basis set files
#include <nlohmann/json.hpp>

// LAMMPS headers
#include "error.h"
#include "memory.h"
#include "comm.h"

#endif
