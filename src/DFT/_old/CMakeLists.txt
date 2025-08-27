# Automatically include all source files in this directory
target_sources(lammps PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/pair_dft.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/basis_set_manager.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/integral_engine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/density_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/xc_functional.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/grid_integrator.cpp
)
