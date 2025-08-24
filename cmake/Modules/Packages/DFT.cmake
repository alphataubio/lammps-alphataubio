# DFT package bootstrap for LAMMPS
# - Finds system installs first
# - Otherwise downloads sources to deterministic tarball names
# - Extracts into staging dirs (never confuses .tar.gz with folders)
# - Installs Eigen3 locally to provide Eigen3Config.cmake for libint2
# - Builds LibXC (6.2.2) and libint2 (2.7.2) via add_subdirectory
# - Avoids cache collisions by using DFT_* variables

function(download_with_fallback url outpath expected_md5)
  message(STATUS "[download] ${url}")
  file(DOWNLOAD "${url}" "${outpath}" STATUS _st SHOW_PROGRESS TLS_VERIFY ON)
  list(GET _st 0 _code)
  if(NOT _code EQUAL 0)
    find_program(CURL_EXE curl)
    if(CURL_EXE)
      execute_process(COMMAND "${CURL_EXE}" -L --fail --retry 3 -o "${outpath}" "${url}"
                      RESULT_VARIABLE _curl_rc)
      if(NOT _curl_rc EQUAL 0)
        message(WARNING "[download] curl failed rc=${_curl_rc}")
      endif()
    endif()
    if(NOT EXISTS "${outpath}")
      find_program(WGET_EXE wget)
      if(WGET_EXE)
        execute_process(COMMAND "${WGET_EXE}" -O "${outpath}" "${url}"
                        RESULT_VARIABLE _wget_rc)
        if(NOT _wget_rc EQUAL 0)
          message(WARNING "[download] wget failed rc=${_wget_rc}")
        endif()
      endif()
    endif()
  endif()
  if(NOT EXISTS "${outpath}")
    message(FATAL_ERROR "[download] Failed to download ${url}")
  endif()
  if(NOT "${expected_md5}" STREQUAL "")
    file(MD5 "${outpath}" _md5)
    if(NOT "${_md5}" STREQUAL "${expected_md5}")
      message(FATAL_ERROR "MD5 mismatch for ${outpath}: got ${_md5}, expected ${expected_md5}")
    endif()
  endif()
endfunction()

function(download_and_extract name url tarball stage_dir expected_md5 OUT_VAR)
  if(EXISTS "${tarball}")
    file(MD5 "${tarball}" _dl_md5)
  endif()
  if("${expected_md5}" STREQUAL "" OR NOT "${_dl_md5}" STREQUAL "${expected_md5}")
    download_with_fallback("${url}" "${tarball}" "${expected_md5}")
  else()
    message(STATUS "[${name}] Using cached tarball ${tarball}")
  endif()
  file(REMOVE_RECURSE "${stage_dir}")
  file(MAKE_DIRECTORY "${stage_dir}")
  message(STATUS "[${name}] Extracting into ${stage_dir}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf "${tarball}"
                  WORKING_DIRECTORY "${stage_dir}" RESULT_VARIABLE _xrc)
  if(NOT _xrc EQUAL 0)
    message(FATAL_ERROR "[${name}] Extraction failed rc=${_xrc}")
  endif()
  file(GLOB _dirs "${stage_dir}/*")
  foreach(d IN LISTS _dirs)
    if(IS_DIRECTORY "${d}")
      set(${OUT_VAR} "${d}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  message(FATAL_ERROR "[${name}] No extracted directory found")
endfunction()

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


# ---------------- Eigen3 vendoring ----------------
if(NOT Eigen3_FOUND)
  set(DFT_EIGEN_URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz")
  set(DFT_EIGEN_MD5 "4c527a9171d71a72a9d4186e65bea559")
  set(_eigen_tar   "${CMAKE_BINARY_DIR}/eigen-3.4.0.tar.gz")
  set(_eigen_stage "${CMAKE_BINARY_DIR}/_eigen_src")
  download_and_extract("Eigen3" "${DFT_EIGEN_URL}" "${_eigen_tar}" "${_eigen_stage}" "${DFT_EIGEN_MD5}" _e_srcdir)
  
  set(EIGEN_PREFIX "${CMAKE_BINARY_DIR}/_eigen_install")
  file(REMOVE_RECURSE "${EIGEN_PREFIX}")  # Clean previous install
  
  execute_process(
    COMMAND ${CMAKE_COMMAND} -S "${_e_srcdir}" -B "${CMAKE_BINARY_DIR}/build-eigen"
      -DCMAKE_INSTALL_PREFIX="${EIGEN_PREFIX}"
      -DBUILD_TESTING=OFF
      -DEIGEN_BUILD_DOC=OFF
    RESULT_VARIABLE _rc
  )
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "Failed to configure Eigen3")
  endif()
  
  execute_process(
    COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}/build-eigen" --target install
    RESULT_VARIABLE _rc
  )
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "Failed to install Eigen3")
  endif()
  
  # Set these BEFORE including libint2
  set(Eigen3_DIR "${EIGEN_PREFIX}/share/eigen3/cmake" CACHE PATH "" FORCE)
  set(EIGEN3_INCLUDE_DIR "${EIGEN_PREFIX}/include/eigen3" CACHE PATH "" FORCE)
  set(EIGEN3_ROOT_DIR "${EIGEN_PREFIX}" CACHE PATH "" FORCE)
  
  # Create an Eigen3::Eigen target if it doesn't exist
  if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}"
    )
  endif()
  
  # Mark Eigen3 as found
  set(Eigen3_FOUND TRUE CACHE BOOL "" FORCE)
  
  target_include_directories(lammps PRIVATE "${EIGEN3_INCLUDE_DIR}")
endif()

# ---------------- libint2 vendoring ----------------
if(NOT LIBINT2_FOUND AND NOT Libint2_FOUND)
  set(DFT_LIBINT2_URL "https://github.com/evaleev/libint/releases/download/v2.7.2/libint-2.7.2.tgz")
  set(DFT_LIBINT2_MD5 "")
  set(_libint_tar   "${CMAKE_BINARY_DIR}/libint-2.7.2.tgz")
  set(_libint_stage "${CMAKE_BINARY_DIR}/_libint2_src")
  download_and_extract("libint2" "${DFT_LIBINT2_URL}" "${_libint_tar}" "${_libint_stage}" "${DFT_LIBINT2_MD5}" libint2_source_dir)
  
  # CRITICAL: Set these BEFORE add_subdirectory to prevent libint2 from finding the wrong Eigen
  set(EIGEN3_INCLUDE_DIR "${EIGEN_PREFIX}/include/eigen3" CACHE PATH "" FORCE)
  set(EIGEN3_ROOT_DIR "${EIGEN_PREFIX}" CACHE PATH "" FORCE)
  set(Eigen3_DIR "${EIGEN_PREFIX}/share/eigen3/cmake" CACHE PATH "" FORCE)
  set(Eigen3_FOUND TRUE CACHE BOOL "" FORCE)
  
  # Disable libint2's own Eigen detection
  set(REQUIRE_EIGEN FALSE CACHE BOOL "" FORCE)
  set(LIBINT2_REQUIRE_EIGEN FALSE CACHE BOOL "" FORCE)
  
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)
  set(CMAKE_POLICY_VERSION "3.5...${CMAKE_VERSION}")
  
  add_subdirectory(${libint2_source_dir} ${CMAKE_BINARY_DIR}/build-libint2)
  
  target_link_libraries(lammps PRIVATE libint2)
  target_include_directories(lammps PRIVATE ${libint2_source_dir}/include ${CMAKE_BINARY_DIR}/build-libint2/include)
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBINT2)
endif()

# ---------------- LibXC vendoring ----------------
if(NOT LIBXC_FOUND)
  set(DFT_LIBXC_URL "https://gitlab.com/libxc/libxc/-/archive/6.2.2/libxc-6.2.2.tar.gz")
  set(DFT_LIBXC_MD5 "d96888788fda9864d17271049a6fa06d")
  set(_libxc_tar   "${CMAKE_BINARY_DIR}/libxc-6.2.2.tar.gz")
  set(_libxc_stage "${CMAKE_BINARY_DIR}/_libxc_src")
  download_and_extract("LibXC" "${DFT_LIBXC_URL}" "${_libxc_tar}" "${_libxc_stage}" "${DFT_LIBXC_MD5}" libxc_source_dir)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)
  set(CMAKE_POLICY_VERSION "3.5...${CMAKE_VERSION}")
  add_subdirectory(${libxc_source_dir} ${CMAKE_BINARY_DIR}/build-libxc)
  target_link_libraries(lammps PRIVATE xc)
  target_include_directories(lammps PRIVATE ${libxc_source_dir}/src ${CMAKE_BINARY_DIR}/build-libxc)
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
else()
  target_include_directories(lammps PRIVATE ${LIBXC_INCLUDE_DIRS})
  target_link_libraries(lammps PRIVATE ${LIBXC_LIBRARIES})
  target_compile_definitions(lammps PRIVATE -DLAMMPS_LIBXC)
endif()

# ---------------- Optional BLAS/LAPACK ----------------
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
if(BLAS_FOUND AND LAPACK_FOUND)
  target_link_libraries(lammps PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# ---------------- Optional OpenMP/MPI ----------------
if(BUILD_OMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(lammps PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(lammps PRIVATE -DDFT_USE_OPENMP)
  endif()
endif()
if(BUILD_MPI)
  find_package(MPI QUIET)
  if(MPI_CXX_FOUND)
    target_include_directories(lammps PRIVATE ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(lammps PRIVATE ${MPI_CXX_LIBRARIES})
    target_compile_definitions(lammps PRIVATE -DDFT_USE_MPI)
  endif()
endif()

message(STATUS "================ DFT dependency summary ================")
message(STATUS "  LibXC_FOUND   = ${LIBXC_FOUND}")
message(STATUS "  Libint2_FOUND = ${LIBINT2_FOUND}")
message(STATUS "  Eigen3_FOUND  = ${Eigen3_FOUND}")
message(STATUS "========================================================")

