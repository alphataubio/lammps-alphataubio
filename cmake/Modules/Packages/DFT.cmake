# DFT package using NWChemEx SCF module
# This properly integrates SCF with all its dependencies

# First, we need to get the NWX CMake infrastructure
include(FetchContent)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Get NWX CMake modules
FetchContent_Declare(
  nwx_cmake
  GIT_REPOSITORY https://github.com/NWChemEx/NWXCMake
  GIT_TAG master
)
FetchContent_MakeAvailable(nwx_cmake)
list(APPEND CMAKE_MODULE_PATH "${nwx_cmake_SOURCE_DIR}/cmake")

# Include NWX versions first
include(nwx_versions)

# Get CMakePPLang (dependency of CMaize)
FetchContent_Declare(
  cmakepp_lang
  GIT_REPOSITORY https://github.com/CMakePP/CMakePPLang
  GIT_TAG v1.0.4
)
FetchContent_MakeAvailable(cmakepp_lang)

# Add CMakePPLang to the module path
list(APPEND CMAKE_MODULE_PATH "${cmakepp_lang_SOURCE_DIR}/cmake")

# Get CMaize build system
FetchContent_Declare(
  cmaize
  GIT_REPOSITORY https://github.com/CMakePP/CMaize
  GIT_TAG ${NWX_CMAIZE_VERSION}
)
FetchContent_MakeAvailable(cmaize)

# Add CMaize to the module path - need both paths for proper includes
list(APPEND CMAKE_MODULE_PATH "${cmaize_SOURCE_DIR}/cmake")

# Include CMaize
include(cmaize/cmaize)

# Build options
set(BUILD_TESTING OFF CACHE BOOL "")
set(BUILD_PYBIND11_PYBINDINGS OFF CACHE BOOL "")
set(BUILD_TAMM_SCF OFF CACHE BOOL "")  # We don't need TAMM for basic SCF
set(ENABLE_SIGMA OFF CACHE BOOL "")    # Disable SIGMA for now

# Set NWX module directory (where plugins will be installed)
set(NWX_MODULE_DIR "${CMAKE_BINARY_DIR}/nwx_modules" CACHE PATH "")

# Now we can use the NWX CMake functions to get dependencies
# Get SimDE (required by SCF)
include(get_simde)

# Get GauXC (required by SCF for XC integration)  
include(get_gauxc)

# Build SCF using the NWX infrastructure
cmaize_find_or_build_dependency(
  scf
  URL github.com/NWChemEx/SCF
  BUILD_TARGET scf
  FIND_TARGET nwx::scf
  CMAKE_ARGS BUILD_TESTING=OFF
             BUILD_PYBIND11_PYBINDINGS=OFF
             BUILD_TAMM_SCF=OFF
             ENABLE_SIGMA=OFF
)

# DFT sources for LAMMPS
target_sources(lammps PRIVATE
  ${LAMMPS_SOURCE_DIR}/DFT/pair_dft.cpp
)

# Link SCF and its dependencies
target_link_libraries(lammps PRIVATE scf)

# Compile definitions
target_compile_definitions(lammps PRIVATE -DLAMMPS_DFT)

message(STATUS "DFT package enabled with NWChemEx SCF module")
