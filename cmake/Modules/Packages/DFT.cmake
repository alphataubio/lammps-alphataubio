# DFT package configuration for LAMMPS
# Uses system-installed libint2, libxc, and eigen3

# Find required packages
find_package(PkgConfig REQUIRED)

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Find LibXC via pkg-config
pkg_check_modules(LIBXC REQUIRED IMPORTED_TARGET libxc)

# Find libint2 via pkg-config
pkg_check_modules(LIBINT2 REQUIRED IMPORTED_TARGET libint2)

# Add includes and libraries to lammps target
target_link_libraries(lammps PRIVATE
  Eigen3::Eigen
  PkgConfig::LIBXC
  PkgConfig::LIBINT2
)

# Add compile definitions
target_compile_definitions(lammps PRIVATE 
  -DLAMMPS_LIBINT2
  -DLAMMPS_LIBXC
)

# Optional BLAS/LAPACK
find_package(BLAS)
find_package(LAPACK)
if(BLAS_FOUND AND LAPACK_FOUND)
  target_link_libraries(lammps PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  message(STATUS "[DFT] BLAS/LAPACK found")
endif()

# Optional MPI support
if(BUILD_MPI)
  find_package(MPI)
  if(MPI_CXX_FOUND)
    target_include_directories(lammps PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(lammps PRIVATE ${MPI_CXX_LIBRARIES})
    target_compile_definitions(lammps PRIVATE -DDFT_USE_MPI)
    message(STATUS "[DFT] MPI support enabled")
  endif()
endif()

# Optional OpenMP support
if(BUILD_OMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(lammps PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(lammps PRIVATE -DDFT_USE_OPENMP)
    message(STATUS "[DFT] OpenMP support enabled")
  endif()
endif()

message(STATUS "================ DFT package configuration ================")
message(STATUS "  Eigen3:      ${EIGEN3_INCLUDE_DIR}")
message(STATUS "  LibXC:       ${LIBXC_INCLUDE_DIRS}")
message(STATUS "  libint2:     ${LIBINT2_INCLUDE_DIRS}")
message(STATUS "  Libraries:   ${LIBXC_LIBRARIES} ${LIBINT2_LIBRARIES}")
message(STATUS "===========================================================")
