# DFT package bootstrap for LAMMPS using FetchContent
# - Finds system installs first
# - Otherwise uses FetchContent to download and build dependencies
# - Installs Eigen3 locally to provide Eigen3Config.cmake for libint2
# - Builds LibXC (6.2.2) and libint2 (2.7.2) via FetchContent
# - Avoids cache collisions by using DFT_* variables

include(FetchContent)

# Set FetchContent to be quiet by default
set(FETCHCONTENT_QUIET ON)

# Set CMake policies to avoid warnings with external dependencies
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)  # option() honors normal variables
endif()
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)  # Keep FindBoost module for older dependencies
endif()
if(POLICY CMP0169)
  cmake_policy(SET CMP0169 OLD)  # Allow deprecated FetchContent_Populate if needed
endif()
if(POLICY CMP0148)
  cmake_policy(SET CMP0148 OLD)  # Keep FindPythonInterp for older dependencies
endif()

# ---------------- System checks ----------------
find_package(PkgConfig QUIET)
set(LIBXC_FOUND FALSE)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(LIBXC QUIET libxc>=5.0.0)
endif()
if(NOT LIBXC_FOUND)
  find_path(LIBXC_INCLUDE_DIR NAMES xc.h PATH_SUFFIXES libxc)
  find_library(LIBXC_LIBRARY NAMES xc libxc)
  if(LIBXC_INCLUDE_DIR AND LIBXC_LIBRARY)
    set(LIBXC_FOUND TRUE)
    set(LIBXC_INCLUDE_DIRS ${LIBXC_INCLUDE_DIR})
    set(LIBXC_LIBRARIES ${LIBXC_LIBRARY})
  endif()
endif()

set(LIBINT2_FOUND FALSE)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(LIBINT2 QUIET libint2)
endif()
if(NOT LIBINT2_FOUND)
  find_path(LIBINT2_INCLUDE_DIR NAMES libint2.hpp)
  find_library(LIBINT2_LIBRARY NAMES int2 libint2 libint)
  if(LIBINT2_INCLUDE_DIR AND LIBINT2_LIBRARY)
    set(LIBINT2_FOUND TRUE)
    set(LIBINT2_INCLUDE_DIRS ${LIBINT2_INCLUDE_DIR})
    set(LIBINT2_LIBRARIES ${LIBINT2_LIBRARY})
  endif()
endif()

find_package(Eigen3 QUIET)

# ---------------- Eigen3 vendoring via FetchContent ----------------
if(NOT Eigen3_FOUND)
  message(STATUS "[DFT] Fetching Eigen3...")
  
  # Configure Eigen3 build options before fetching
  set(BUILD_TESTING OFF CACHE BOOL "Disable Eigen3 testing" FORCE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "Disable Eigen3 documentation" FORCE)
  set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "Disable Eigen3 pkgconfig" FORCE)
  
  FetchContent_Declare(
    dft_eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH MD5=4c527a9171d71a72a9d4186e65bea559
    DOWNLOAD_EXTRACT_TIMESTAMP ON
  )
  
  # Use modern FetchContent approach
  FetchContent_MakeAvailable(dft_eigen3)
  
  # Set up Eigen3 paths for libint2 compatibility
  set(EIGEN3_INCLUDE_DIR "${dft_eigen3_SOURCE_DIR}" CACHE PATH "Eigen3 include directory" FORCE)
  set(EIGEN3_ROOT_DIR "${dft_eigen3_SOURCE_DIR}" CACHE PATH "Eigen3 root directory" FORCE)
  set(Eigen3_FOUND TRUE CACHE BOOL "Eigen3 found flag" FORCE)
  
  # Add the Eigen3 alias target if it doesn't exist
  if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${dft_eigen3_SOURCE_DIR}"
    )
  endif()
  
  message(STATUS "[DFT] Eigen3 configured from source: ${dft_eigen3_SOURCE_DIR}")
  
  # Add include directory to lammps target
  target_include_directories(lammps PRIVATE "${dft_eigen3_SOURCE_DIR}")
endif()

# ---------------- Minimal Boost Preprocessor for libint2 ----------------
# Download only the Boost Preprocessor headers that libint2 needs
message(STATUS "[DFT] Setting up minimal Boost Preprocessor for libint2...")

# Check if we already have Boost
find_package(Boost 1.56 QUIET)

if(NOT Boost_FOUND)
  # Download minimal Boost - just preprocessor and its dependencies
  # This is much smaller than the full Boost distribution  
  FetchContent_Declare(
    dft_boost_preprocessor
    URL https://github.com/boostorg/preprocessor/archive/refs/tags/boost-1.82.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP ON
  )
  
  FetchContent_Declare(
    dft_boost_config
    URL https://github.com/boostorg/config/archive/refs/tags/boost-1.82.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP ON
  )
  
  # Download both
  FetchContent_MakeAvailable(dft_boost_preprocessor dft_boost_config)
  
  # Create the proper boost include structure
  set(DFT_BOOST_ROOT "${CMAKE_BINARY_DIR}/_boost_minimal")
  file(MAKE_DIRECTORY "${DFT_BOOST_ROOT}/boost")
  
  # Copy headers maintaining proper structure
  # First copy config headers
  file(GLOB_RECURSE CONFIG_HEADERS "${dft_boost_config_SOURCE_DIR}/include/boost/*.hpp")
  foreach(header ${CONFIG_HEADERS})
    file(RELATIVE_PATH rel_path "${dft_boost_config_SOURCE_DIR}/include" "${header}")
    configure_file("${header}" "${DFT_BOOST_ROOT}/${rel_path}" COPYONLY)
  endforeach()
  
  # Then copy preprocessor headers 
  file(GLOB_RECURSE PP_HEADERS "${dft_boost_preprocessor_SOURCE_DIR}/include/boost/*.hpp")
  foreach(header ${PP_HEADERS})
    file(RELATIVE_PATH rel_path "${dft_boost_preprocessor_SOURCE_DIR}/include" "${header}")
    configure_file("${header}" "${DFT_BOOST_ROOT}/${rel_path}" COPYONLY)
  endforeach()
  
  # Create boost/version.hpp for compatibility
  file(WRITE "${DFT_BOOST_ROOT}/boost/version.hpp"
"//  Boost version.hpp configuration header
#ifndef BOOST_VERSION_HPP
#define BOOST_VERSION_HPP
#define BOOST_VERSION 108200
#define BOOST_LIB_VERSION \"1_82\"
#define BOOST_VERSION_MAJOR 1
#define BOOST_VERSION_MINOR 82
#define BOOST_VERSION_PATCH 0
#endif
")
  
  # Also create a main preprocessor.hpp that includes all sub-headers
  file(WRITE "${DFT_BOOST_ROOT}/boost/preprocessor.hpp"
"// Main preprocessor header
#ifndef BOOST_PREPROCESSOR_HPP
#define BOOST_PREPROCESSOR_HPP
#include <boost/preprocessor/arithmetic.hpp>
#include <boost/preprocessor/array.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/comparison.hpp>
#include <boost/preprocessor/config/config.hpp>
#include <boost/preprocessor/control.hpp>
#include <boost/preprocessor/debug.hpp>
#include <boost/preprocessor/facilities.hpp>
#include <boost/preprocessor/iteration.hpp>
#include <boost/preprocessor/library.hpp>
#include <boost/preprocessor/list.hpp>
#include <boost/preprocessor/logical.hpp>
#include <boost/preprocessor/punctuation.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/selection.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/slot.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/tuple.hpp>
#include <boost/preprocessor/variadic.hpp>
#endif
")
  
  # Set up Boost variables for find_package and libint2
  set(BOOST_ROOT "${DFT_BOOST_ROOT}" CACHE PATH "Boost root directory" FORCE)
  set(BOOST_INCLUDEDIR "${DFT_BOOST_ROOT}" CACHE PATH "Boost include directory" FORCE) 
  set(Boost_INCLUDE_DIR "${DFT_BOOST_ROOT}" CACHE PATH "Boost include directory" FORCE)
  set(Boost_INCLUDE_DIRS "${DFT_BOOST_ROOT}" CACHE PATH "Boost include directories" FORCE)
  set(Boost_FOUND TRUE CACHE BOOL "Boost found flag" FORCE)
  set(Boost_VERSION "1.82.0" CACHE STRING "Boost version" FORCE)
  set(Boost_VERSION_STRING "1.82.0" CACHE STRING "Boost version string" FORCE)
  set(Boost_VERSION_MACRO "108200" CACHE STRING "Boost version macro" FORCE)
  set(Boost_MAJOR_VERSION "1" CACHE STRING "Boost major version" FORCE)
  set(Boost_MINOR_VERSION "82" CACHE STRING "Boost minor version" FORCE)
  set(Boost_SUBMINOR_VERSION "0" CACHE STRING "Boost subminor version" FORCE)
  
  # Add to CMAKE_PREFIX_PATH and include directories
  list(APPEND CMAKE_PREFIX_PATH "${DFT_BOOST_ROOT}")
  set(CMAKE_INCLUDE_PATH "${DFT_BOOST_ROOT};${CMAKE_INCLUDE_PATH}" CACHE STRING "" FORCE)
  
  message(STATUS "[DFT] Minimal Boost Preprocessor configured: ${DFT_BOOST_ROOT}")
else()
  message(STATUS "[DFT] Using system Boost: ${Boost_INCLUDE_DIRS}")
  set(DFT_BOOST_ROOT "${Boost_INCLUDE_DIRS}" CACHE PATH "Boost root directory" FORCE)
endif()

# ---------------- libint2 vendoring via FetchContent ----------------
if(NOT LIBINT2_FOUND AND NOT Libint2_FOUND)
  message(STATUS "[DFT] Fetching libint2...")
  
  FetchContent_Declare(
    dft_libint2
    URL https://github.com/evaleev/libint/releases/download/v2.7.2/libint-2.7.2.tgz
    DOWNLOAD_EXTRACT_TIMESTAMP ON
  )
  
  # Set libint2 build options
  set(REQUIRE_EIGEN TRUE CACHE BOOL "" FORCE)
  set(LIBINT2_REQUIRE_EIGEN TRUE CACHE BOOL "" FORCE)
  set(LIBINT2_BUILD_SHARED_AND_STATIC_LIBS OFF CACHE BOOL "" FORCE)
  set(ENABLE_GENERIC ON CACHE BOOL "" FORCE)
  set(LIBINT2_REALTYPE "double" CACHE STRING "" FORCE)
  
  # Tell libint2 where to find dependencies
  set(EIGEN3_INCLUDE_DIR "${EIGEN3_INCLUDE_DIR}" CACHE PATH "" FORCE)
  set(BOOST_ROOT "${DFT_BOOST_ROOT}" CACHE PATH "" FORCE)
  set(BOOST_INCLUDEDIR "${DFT_BOOST_ROOT}" CACHE PATH "" FORCE)
  set(Boost_INCLUDE_DIR "${DFT_BOOST_ROOT}" CACHE PATH "" FORCE)
  set(Boost_INCLUDE_DIRS "${DFT_BOOST_ROOT}" CACHE PATH "" FORCE)
  
  # Configure and build libint2
  FetchContent_MakeAvailable(dft_libint2)
  
  target_link_libraries(lammps PRIVATE libint2)
  target_include_directories(lammps PRIVATE 
    ${dft_libint2_SOURCE_DIR}/include 
    ${dft_libint2_BINARY_DIR}/include
    ${DFT_BOOST_ROOT}  # Add our Boost headers to the include path!
  )
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBINT2)
  
  message(STATUS "[DFT] libint2 configured")
endif()

# ---------------- LibXC vendoring via FetchContent ----------------
if(NOT LIBXC_FOUND)
  message(STATUS "[DFT] Fetching LibXC...")
  
  # Set LibXC build options before fetching
  set(DISABLE_KXC ON CACHE BOOL "Disable LibXC KXC functionals" FORCE)
  set(DISABLE_LXC ON CACHE BOOL "Disable LibXC LXC functionals" FORCE)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static LibXC" FORCE)
  
  FetchContent_Declare(
    dft_libxc
    URL https://gitlab.com/libxc/libxc/-/archive/6.2.2/libxc-6.2.2.tar.gz
    URL_HASH MD5=d96888788fda9864d17271049a6fa06d
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    PATCH_COMMAND ${CMAKE_COMMAND} -E echo "Patching LibXC CMakeLists.txt..." &&
                  sed -i.bak "s/cmake_minimum_required(VERSION [0-9.]*/cmake_minimum_required(VERSION 3.5/" <SOURCE_DIR>/CMakeLists.txt || true
  )
  
  # Use modern FetchContent approach
  FetchContent_MakeAvailable(dft_libxc)
  
  target_link_libraries(lammps PRIVATE xc)
  target_include_directories(lammps PRIVATE 
    ${dft_libxc_SOURCE_DIR}/src 
    ${dft_libxc_BINARY_DIR}
  )
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
  
  message(STATUS "[DFT] LibXC configured")
else()
  target_include_directories(lammps PRIVATE ${LIBXC_INCLUDE_DIRS})
  target_link_libraries(lammps PRIVATE ${LIBXC_LIBRARIES})
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
  message(STATUS "[DFT] Using system LibXC")
endif()

# ---------------- Optional BLAS/LAPACK ----------------
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
if(BLAS_FOUND AND LAPACK_FOUND)
  target_link_libraries(lammps PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  message(STATUS "[DFT] BLAS/LAPACK found")
endif()

# ---------------- Optional OpenMP/MPI ----------------
if(BUILD_OMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(lammps PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(lammps PRIVATE -DDFT_USE_OPENMP)
    message(STATUS "[DFT] OpenMP support enabled")
  endif()
endif()

if(BUILD_MPI)
  find_package(MPI QUIET)
  if(MPI_CXX_FOUND)
    target_include_directories(lammps PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(lammps PRIVATE ${MPI_CXX_LIBRARIES})
    target_compile_definitions(lammps PRIVATE -DDFT_USE_MPI)
    message(STATUS "[DFT] MPI support enabled")
  endif()
endif()

message(STATUS "================ DFT dependency summary ================")
message(STATUS "  LibXC_FOUND   = ${LIBXC_FOUND}")
message(STATUS "  Libint2_FOUND = ${LIBINT2_FOUND}")
message(STATUS "  Eigen3_FOUND  = ${Eigen3_FOUND}")
message(STATUS "  Boost_FOUND   = ${Boost_FOUND}")
if(BLAS_FOUND AND LAPACK_FOUND)
  message(STATUS "  BLAS/LAPACK   = Found")
else()
  message(STATUS "  BLAS/LAPACK   = Not found")
endif()
if(BUILD_OMP AND OpenMP_CXX_FOUND)
  message(STATUS "  OpenMP        = Enabled")
endif()
if(BUILD_MPI AND MPI_CXX_FOUND)
  message(STATUS "  MPI           = Enabled")
endif()
message(STATUS "========================================================")
