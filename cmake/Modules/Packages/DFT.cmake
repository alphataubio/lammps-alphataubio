# LibXC library support for DFT package
find_package(PkgConfig QUIET)

# First try to find an existing libxc installation
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

# If libxc is not found, download and build it
if(NOT LIBXC_FOUND)
  message(STATUS "LibXC not found in system. Will download and build LibXC automatically.")
  
  # Set policy to silence warnings about timestamps of downloaded files
  if(POLICY CMP0135)
    cmake_policy(SET CMP0135 OLD)
  endif()
  
  # LibXC download configuration
  # Using LibXC 6.2.2 which has good CMake support
  set(LIBXC_URL "https://gitlab.com/libxc/libxc/-/archive/6.2.2/libxc-6.2.2.tar.gz" CACHE STRING "URL for LibXC library sources")
  set(LIBXC_MD5 "6c866168040a0c49879c8963ec4a2fb7" CACHE STRING "MD5 checksum of LibXC library tarball")
  mark_as_advanced(LIBXC_URL)
  mark_as_advanced(LIBXC_MD5)
  
  # Get fallback URL for backup download location
  GetFallbackURL(LIBXC_URL LIBXC_FALLBACK)
  
  # Option for using a local LibXC directory (useful for development)
  # Use LOCAL_DFT to be consistent with other packages like LOCAL_ML-PACE
  if(LOCAL_DFT)
    set(libxc_source_dir "${LOCAL_DFT}")
    message(STATUS "Using local LibXC directory: ${libxc_source_dir}")
  else()
    # Download LibXC if not already present or if checksum doesn't match
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
    
    # Extract LibXC sources
    message(STATUS "Extracting LibXC sources...")
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E remove_directory libxc-*
      COMMAND ${CMAKE_COMMAND} -E tar xzf libxc.tar.gz
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
    
    # Find the extracted directory
    get_newest_file(${CMAKE_BINARY_DIR}/libxc-* libxc_source_dir)
  endif()
  
  # LibXC 6.x includes CMake support, so we can use add_subdirectory
  # Check if CMakeLists.txt exists in the source directory
  if(EXISTS ${libxc_source_dir}/CMakeLists.txt)
    # Configure LibXC build options
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static LibXC library" FORCE)
    set(ENABLE_FORTRAN OFF CACHE BOOL "Disable Fortran interface" FORCE)
    set(ENABLE_CUDA OFF CACHE BOOL "Disable CUDA support" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "Disable LibXC tests" FORCE)
    
    # Add LibXC as a subdirectory
    add_subdirectory(${libxc_source_dir} build-libxc EXCLUDE_FROM_ALL)
    
    # Link LAMMPS with the built LibXC
    target_link_libraries(lammps PRIVATE xc)
    
    # Get the include directory from the LibXC target
    get_target_property(LIBXC_INCLUDE_DIRS xc INTERFACE_INCLUDE_DIRECTORIES)
    if(LIBXC_INCLUDE_DIRS)
      target_include_directories(lammps PRIVATE ${LIBXC_INCLUDE_DIRS})
    endif()
    
    # Also add the source directory for generated headers
    target_include_directories(lammps PRIVATE ${CMAKE_BINARY_DIR}/build-libxc)
    target_include_directories(lammps PRIVATE ${libxc_source_dir}/src)
    
    target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
    
    message(STATUS "LibXC will be built automatically from ${libxc_source_dir}")
  else()
    # Fallback for older LibXC versions without CMake support
    # Use ExternalProject to handle autotools build
    include(ExternalProject)
    
    set(LIBXC_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/libxc_install)
    
    # Configure command for autotools
    set(LIBXC_CONFIGURE_CMD ${libxc_source_dir}/configure 
        --prefix=${LIBXC_INSTALL_PREFIX}
        --disable-shared
        --enable-static
        --disable-fortran
        CC=${CMAKE_C_COMPILER}
        CXX=${CMAKE_CXX_COMPILER}
    )
    
    ExternalProject_Add(
      libxc_build
      SOURCE_DIR ${libxc_source_dir}
      CONFIGURE_COMMAND ${LIBXC_CONFIGURE_CMD}
      BUILD_COMMAND make -j
      INSTALL_COMMAND make install
      BUILD_IN_SOURCE 1
      BUILD_BYPRODUCTS ${LIBXC_INSTALL_PREFIX}/lib/libxc.a
    )
    
    # Create imported target for LibXC
    add_library(libxc_imported STATIC IMPORTED GLOBAL)
    set_target_properties(libxc_imported PROPERTIES
      IMPORTED_LOCATION ${LIBXC_INSTALL_PREFIX}/lib/libxc.a
    )
    
    # Make sure LibXC is built before LAMMPS
    add_dependencies(lammps libxc_build)
    
    # Link LAMMPS with the built LibXC
    target_link_libraries(lammps PRIVATE libxc_imported)
    target_include_directories(lammps PRIVATE ${LIBXC_INSTALL_PREFIX}/include)
    target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
    
    message(STATUS "LibXC (autotools version) will be built automatically from ${libxc_source_dir}")
    message(STATUS "LibXC will be installed to ${LIBXC_INSTALL_PREFIX}")
  endif()
  
else()
  # Use system LibXC
  target_include_directories(lammps PRIVATE ${LIBXC_INCLUDE_DIRS})
  target_link_libraries(lammps PRIVATE ${LIBXC_LIBRARIES})
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
  
  message(STATUS "Found LibXC:")
  message(STATUS "  Include directories: ${LIBXC_INCLUDE_DIRS}")
  message(STATUS "  Libraries: ${LIBXC_LIBRARIES}")
  
  # Check LibXC version if available
  if(LIBXC_VERSION)
    message(STATUS "  Version: ${LIBXC_VERSION}")
    if(LIBXC_VERSION VERSION_LESS "5.0.0")
      message(WARNING "LibXC version ${LIBXC_VERSION} is older than recommended (5.0.0). Some features may not be available.")
    endif()
  endif()
endif()

# Optional: Enable OpenMP for LibXC if available
if(BUILD_OMP)
  target_compile_definitions(lammps PRIVATE -DLIBXC_OPENMP)
endif()

# Optional: Enable MPI for LibXC if available
if(BUILD_MPI)
  target_compile_definitions(lammps PRIVATE -DLIBXC_MPI)
endif()
