# DFT package using NWChemEx SCF (which includes GauXC)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build options
set(BUILD_TESTING OFF CACHE BOOL "")
set(BUILD_PYBIND11_PYBINDINGS OFF CACHE BOOL "")

# Fetch SCF which will handle all its dependencies including GauXC
FetchContent_Declare(
  scf
  GIT_REPOSITORY https://github.com/NWChemEx/SCF.git
  GIT_TAG main
)

# SCF will pull in all needed dependencies:
# - GauXC (grid integration)
# - Chemist (molecules and basis sets)
# - SimDE (simulation framework)
# - PluginPlay (module system)
# - ExchCXX (XC functionals)
# - IntegratorXX (quadratures)
# - Libint2 (integrals)

message(STATUS "DFT: Building SCF and dependencies...")
FetchContent_MakeAvailable(scf)

# DFT sources
target_sources(lammps PRIVATE
  ${LAMMPS_SOURCE_DIR}/DFT/pair_dft.cpp
  ${LAMMPS_SOURCE_DIR}/DFT/pair_dft_scf.cpp
)

# Link SCF
target_link_libraries(lammps PRIVATE nwx::scf)

# Compile definitions
target_compile_definitions(lammps PRIVATE -DLAMMPS_DFT)

message(STATUS "DFT package enabled with NWChemEx SCF")
