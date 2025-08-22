# Find LibXC library for DFT package
find_package(PkgConfig QUIET)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(LIBXC libxc>=5.0.0)
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

if(LIBXC_FOUND)
  target_include_directories(lammps PRIVATE ${LIBXC_INCLUDE_DIRS})
  target_link_libraries(lammps PRIVATE ${LIBXC_LIBRARIES})
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
  message(STATUS "LibXC found:")
  message(STATUS "  Include directories: ${LIBXC_INCLUDE_DIRS}")
  message(STATUS "  Libraries: ${LIBXC_LIBRARIES}")
else()
  message(FATAL_ERROR "LibXC library not found. Please install libxc (version >= 5.0.0) or specify its location.")
endif()

# Optional: Support for different LibXC versions
if(LIBXC_VERSION)
  message(STATUS "  Version: ${LIBXC_VERSION}")
  if(LIBXC_VERSION VERSION_LESS "5.0.0")
    message(WARNING "LibXC version ${LIBXC_VERSION} is older than recommended (5.0.0). Some features may not be available.")
  endif()
endif()
