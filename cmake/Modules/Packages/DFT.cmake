# DFT package for LAMMPS
# Requires: LibXC, libint2, Eigen3, nlohmann_json

# Find required packages
find_package(PkgConfig QUIET)

# ===== LibXC Configuration =====
set(LIBXC_FOUND FALSE)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(LIBXC QUIET libxc>=5.0.0)
endif()

if(NOT LIBXC_FOUND)
  # Try to find LibXC manually
  find_path(LIBXC_INCLUDE_DIR
    NAMES xc.h
    PATHS /usr/include /usr/local/include /opt/local/include
    PATH_SUFFIXES libxc
  )
  
  find_library(LIBXC_LIBRARY
    NAMES xc libxc
    PATHS /usr/lib /usr/local/lib /opt/local/lib /usr/lib64 /usr/local/lib64
  )
  
  if(LIBXC_INCLUDE_DIR AND LIBXC_LIBRARY)
    set(LIBXC_FOUND TRUE)
    set(LIBXC_INCLUDE_DIRS ${LIBXC_INCLUDE_DIR})
    set(LIBXC_LIBRARIES ${LIBXC_LIBRARY})
  endif()
endif()

# Download and build LibXC if not found
if(NOT LIBXC_FOUND)
  message(STATUS "LibXC not found in system. Will download and build LibXC automatically.")
  
  if(POLICY CMP0135)
    cmake_policy(SET CMP0135 OLD)
  endif()
  
  set(LIBXC_URL "https://gitlab.com/libxc/libxc/-/archive/6.2.2/libxc-6.2.2.tar.gz" CACHE STRING "URL for LibXC library sources")
  set(LIBXC_MD5 "6c866168040a0c49879c8963ec4a2fb7" CACHE STRING "MD5 checksum of LibXC library tarball")
  mark_as_advanced(LIBXC_URL)
  mark_as_advanced(LIBXC_MD5)
  
  GetFallbackURL(LIBXC_URL LIBXC_FALLBACK)
  
  if(LOCAL_DFT_LIBXC)
    set(libxc_source_dir "${LOCAL_DFT_LIBXC}")
    message(STATUS "Using local LibXC directory: ${libxc_source_dir}")
  else()
    if(EXISTS ${CMAKE_BINARY_DIR}/libxc.tar.gz)
      file(MD5 ${CMAKE_BINARY_DIR}/libxc.tar.gz DL_MD5)
    endif()
    
    if(NOT "${DL_MD5}" STREQUAL "${LIBXC_MD5}")
      message(STATUS "Downloading ${LIBXC_URL}")
      file(DOWNLOAD ${LIBXC_URL} ${CMAKE_BINARY_DIR}/libxc.tar.gz STATUS DL_STATUS SHOW_PROGRESS)
      file(MD5 ${CMAKE_BINARY_DIR}/libxc.tar.gz DL_MD5)
      if((NOT DL_STATUS EQUAL 0) OR (NOT "${DL_MD5}" STREQUAL "${LIBXC_MD5}"))
        if(LIBXC_FALLBACK)
          message(WARNING "Download from primary URL ${LIBXC_URL} failed\nTrying fallback URL ${LIBXC_FALLBACK}")
          file(DOWNLOAD ${LIBXC_FALLBACK} ${CMAKE_BINARY_DIR}/libxc.tar.gz EXPECTED_HASH MD5=${LIBXC_MD5} SHOW_PROGRESS)
        else()
          message(FATAL_ERROR "Failed to download LibXC from ${LIBXC_URL}")
        endif()
      endif()
    else()
      message(STATUS "Using already downloaded LibXC archive ${CMAKE_BINARY_DIR}/libxc.tar.gz")
    endif()
    
    message(STATUS "Extracting LibXC sources...")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E remove_directory libxc-*
      COMMAND ${CMAKE_COMMAND} -E tar xzf libxc.tar.gz
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    
    get_newest_file(${CMAKE_BINARY_DIR}/libxc-* libxc_source_dir)
  endif()
  
  # Build LibXC with CMake
  if(EXISTS ${libxc_source_dir}/CMakeLists.txt)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static LibXC library" FORCE)
    set(ENABLE_FORTRAN OFF CACHE BOOL "Disable Fortran interface" FORCE)
    set(ENABLE_CUDA OFF CACHE BOOL "Disable CUDA support" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "Disable LibXC tests" FORCE)
    
    add_subdirectory(${libxc_source_dir} build-libxc EXCLUDE_FROM_ALL)
    
    target_link_libraries(lammps PRIVATE xc)
    get_target_property(LIBXC_INCLUDE_DIRS xc INTERFACE_INCLUDE_DIRECTORIES)
    if(LIBXC_INCLUDE_DIRS)
      target_include_directories(lammps PRIVATE ${LIBXC_INCLUDE_DIRS})
    endif()
    target_include_directories(lammps PRIVATE ${CMAKE_BINARY_DIR}/build-libxc)
    target_include_directories(lammps PRIVATE ${libxc_source_dir}/src)
    
    target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
    
    message(STATUS "LibXC will be built automatically from ${libxc_source_dir}")
  else()
    message(FATAL_ERROR "LibXC source directory does not contain CMakeLists.txt")
  endif()
else()
  # Use system LibXC
  target_include_directories(lammps PRIVATE ${LIBXC_INCLUDE_DIRS})
  target_link_libraries(lammps PRIVATE ${LIBXC_LIBRARIES})
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
  
  message(STATUS "Found LibXC:")
  message(STATUS "  Include directories: ${LIBXC_INCLUDE_DIRS}")
  message(STATUS "  Libraries: ${LIBXC_LIBRARIES}")
endif()

# ===== libint2 Configuration =====
set(LIBINT2_FOUND FALSE)

# First try to find libint2 using its CMake config
find_package(Libint2 QUIET CONFIG)

if(NOT Libint2_FOUND)
  # Try pkg-config
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(LIBINT2 QUIET libint2>=2.6.0)
  endif()
  
  if(NOT LIBINT2_FOUND)
    # Manual search
    find_path(LIBINT2_INCLUDE_DIR
      NAMES libint2.hpp
      PATHS /usr/include /usr/local/include /opt/local/include
      PATH_SUFFIXES libint2
    )
    
    find_library(LIBINT2_LIBRARY
      NAMES int2 libint2
      PATHS /usr/lib /usr/local/lib /opt/local/lib /usr/lib64 /usr/local/lib64
    )
    
    if(LIBINT2_INCLUDE_DIR AND LIBINT2_LIBRARY)
      set(LIBINT2_FOUND TRUE)
      set(LIBINT2_INCLUDE_DIRS ${LIBINT2_INCLUDE_DIR})
      set(LIBINT2_LIBRARIES ${LIBINT2_LIBRARY})
    endif()
  endif()
endif()

# Download and build libint2 if not found
if(NOT LIBINT2_FOUND AND NOT Libint2_FOUND)
  message(STATUS "libint2 not found in system. Will download and build libint2 automatically.")
  
  set(LIBINT2_URL "https://github.com/evaleev/libint/releases/download/v2.7.2/libint-2.7.2.tgz" CACHE STRING "URL for libint2 library sources")
  set(LIBINT2_MD5 "37f9e2a0f4e0c22e0e3e6f2b5e2f5b8f" CACHE STRING "MD5 checksum of libint2 library tarball")
  mark_as_advanced(LIBINT2_URL)
  mark_as_advanced(LIBINT2_MD5)
  
  GetFallbackURL(LIBINT2_URL LIBINT2_FALLBACK)
  
  if(LOCAL_DFT_LIBINT2)
    set(libint2_source_dir "${LOCAL_DFT_LIBINT2}")
    message(STATUS "Using local libint2 directory: ${libint2_source_dir}")
  else()
    if(EXISTS ${CMAKE_BINARY_DIR}/libint2.tar.gz)
      file(MD5 ${CMAKE_BINARY_DIR}/libint2.tar.gz DL_MD5)
    endif()
    
    if(NOT "${DL_MD5}" STREQUAL "${LIBINT2_MD5}")
      message(STATUS "Downloading ${LIBINT2_URL}")
      file(DOWNLOAD ${LIBINT2_URL} ${CMAKE_BINARY_DIR}/libint2.tar.gz STATUS DL_STATUS SHOW_PROGRESS)
      file(MD5 ${CMAKE_BINARY_DIR}/libint2.tar.gz DL_MD5)
      if(NOT DL_STATUS EQUAL 0)
        if(LIBINT2_FALLBACK)
          message(WARNING "Download from primary URL ${LIBINT2_URL} failed\nTrying fallback URL ${LIBINT2_FALLBACK}")
          file(DOWNLOAD ${LIBINT2_FALLBACK} ${CMAKE_BINARY_DIR}/libint2.tar.gz SHOW_PROGRESS)
        else()
          message(FATAL_ERROR "Failed to download libint2 from ${LIBINT2_URL}")
        endif()
      endif()
    else()
      message(STATUS "Using already downloaded libint2 archive ${CMAKE_BINARY_DIR}/libint2.tar.gz")
    endif()
    
    message(STATUS "Extracting libint2 sources...")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E remove_directory libint-*
      COMMAND ${CMAKE_COMMAND} -E tar xzf libint2.tar.gz
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    
    get_newest_file(${CMAKE_BINARY_DIR}/libint-* libint2_source_dir)
  endif()
  
  # Build libint2 with CMake
  if(EXISTS ${libint2_source_dir}/CMakeLists.txt)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libint2 library" FORCE)
    set(LIBINT2_BUILD_SHARED OFF CACHE BOOL "Build static libint2" FORCE)
    set(ENABLE_FORTRAN OFF CACHE BOOL "Disable Fortran interface" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "Disable libint2 tests" FORCE)
    
    # Configure libint2 for standard integrals
    set(LIBINT2_SHGAUSS_ORDERING "standard" CACHE STRING "Use standard shell ordering" FORCE)
    set(WITH_MAX_AM 4 CACHE STRING "Maximum angular momentum" FORCE)
    set(WITH_OPT_AM 3 CACHE STRING "Optimized angular momentum" FORCE)
    
    add_subdirectory(${libint2_source_dir} build-libint2 EXCLUDE_FROM_ALL)
    
    target_link_libraries(lammps PRIVATE libint2)
    target_include_directories(lammps PRIVATE ${libint2_source_dir}/include)
    target_include_directories(lammps PRIVATE ${CMAKE_BINARY_DIR}/build-libint2/include)
    
    target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBINT2)
    
    message(STATUS "libint2 will be built automatically from ${libint2_source_dir}")
  else()
    message(WARNING "libint2 source directory does not contain CMakeLists.txt. DFT integrals will be limited.")
    set(LIBINT2_FOUND FALSE)
  endif()
else()
  if(Libint2_FOUND)
    # Use CMake config
    target_link_libraries(lammps PRIVATE Libint2::libint2)
    target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBINT2)
    message(STATUS "Found libint2 via CMake config")
  else()
    # Use system libint2
    target_include_directories(lammps PRIVATE ${LIBINT2_INCLUDE_DIRS})
    target_link_libraries(lammps PRIVATE ${LIBINT2_LIBRARIES})
    target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBINT2)
    
    message(STATUS "Found libint2:")
    message(STATUS "  Include directories: ${LIBINT2_INCLUDE_DIRS}")
    message(STATUS "  Libraries: ${LIBINT2_LIBRARIES}")
  endif()
endif()

# ===== Eigen3 Configuration =====
find_package(Eigen3 3.3 QUIET NO_MODULE)

if(NOT Eigen3_FOUND)
  # Try to find Eigen3 manually
  find_path(EIGEN3_INCLUDE_DIR
    NAMES Eigen/Core
    PATHS /usr/include /usr/local/include /opt/local/include
    PATH_SUFFIXES eigen3
  )
  
  if(EIGEN3_INCLUDE_DIR)
    set(Eigen3_FOUND TRUE)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
    target_include_directories(lammps PRIVATE ${EIGEN3_INCLUDE_DIR})
  else()
    # Download Eigen3 headers
    message(STATUS "Eigen3 not found. Will download Eigen3 headers.")
    
    set(EIGEN3_URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz" CACHE STRING "URL for Eigen3 headers")
    set(EIGEN3_MD5 "4c527a9171d71a72a9d4186e65bea559" CACHE STRING "MD5 checksum of Eigen3 tarball")
    
    if(EXISTS ${CMAKE_BINARY_DIR}/eigen3.tar.gz)
      file(MD5 ${CMAKE_BINARY_DIR}/eigen3.tar.gz DL_MD5)
    endif()
    
    if(NOT "${DL_MD5}" STREQUAL "${EIGEN3_MD5}")
      message(STATUS "Downloading ${EIGEN3_URL}")
      file(DOWNLOAD ${EIGEN3_URL} ${CMAKE_BINARY_DIR}/eigen3.tar.gz EXPECTED_HASH MD5=${EIGEN3_MD5} SHOW_PROGRESS)
    endif()
    
    message(STATUS "Extracting Eigen3 headers...")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E remove_directory eigen-*
      COMMAND ${CMAKE_COMMAND} -E tar xzf eigen3.tar.gz
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    
    get_newest_file(${CMAKE_BINARY_DIR}/eigen-* eigen3_source_dir)
    target_include_directories(lammps PRIVATE ${eigen3_source_dir})
    message(STATUS "Using Eigen3 from ${eigen3_source_dir}")
  endif()
else()
  target_link_libraries(lammps PRIVATE Eigen3::Eigen)
  message(STATUS "Found Eigen3 via CMake config")
endif()

# ===== Optional: BLAS/LAPACK for better performance =====
find_package(BLAS QUIET)
find_package(LAPACK QUIET)

if(BLAS_FOUND AND LAPACK_FOUND)
  target_link_libraries(lammps PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  target_compile_definitions(lammps PRIVATE -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACK)
  message(STATUS "Using BLAS/LAPACK for linear algebra")
endif()

# ===== Optional: OpenMP for parallelization =====
if(BUILD_OMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(lammps PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(lammps PRIVATE -DDFT_USE_OPENMP)
    message(STATUS "DFT package will use OpenMP for parallelization")
  endif()
endif()

# ===== Optional: MPI for parallel DFT =====
if(BUILD_MPI)
  target_compile_definitions(lammps PRIVATE -DDFT_USE_MPI)
  message(STATUS "DFT package will use MPI for parallel calculations")
endif()

# ===== Configure DFT package compilation flags =====
target_compile_definitions(lammps PRIVATE -DLAMMPS_DFT)

# Enable C++14 or higher for modern features
target_compile_features(lammps PRIVATE cxx_std_14)

# Add warning flags for development
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(lammps PRIVATE 
    $<$<CONFIG:Debug>:-Wall -Wextra -Wpedantic>
  )
endif()

# Summary message
message(STATUS "================================================")
message(STATUS "DFT Package Configuration Summary:")
message(STATUS "  LibXC: ${LIBXC_FOUND}")
if(LIBINT2_FOUND OR Libint2_FOUND)
  message(STATUS "  libint2: Found")
else()
  message(STATUS "  libint2: Not found (limited integral support)")
endif()
message(STATUS "  Eigen3: ${Eigen3_FOUND}")
if(BLAS_FOUND AND LAPACK_FOUND)
  message(STATUS "  BLAS/LAPACK: Enabled")
else()
  message(STATUS "  BLAS/LAPACK: Disabled")
endif()
if(BUILD_OMP AND OpenMP_CXX_FOUND)
  message(STATUS "  OpenMP: Enabled")
else()
  message(STATUS "  OpenMP: Disabled")
endif()
if(BUILD_MPI)
  message(STATUS "  MPI: Enabled")
else()
  message(STATUS "  MPI: Disabled")
endif()
message(STATUS "================================================")
