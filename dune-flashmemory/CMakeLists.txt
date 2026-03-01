# CMakeLists.txt for Parallel Flash Memory High-K Simulator
# Build: dunecontrol --only=dune-flashhighk all
# Run:   mpiexec -np 8 ./dune-flashhighk/build-cmake/src/dune_highk_parallel

cmake_minimum_required(VERSION 3.13)
project(dune-flashhighk CXX)

# Find DUNE modules
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")
include(DuneMacros)

dune_project()

find_package(dune-grid REQUIRED)
find_package(dune-alugrid REQUIRED)
find_package(dune-istl REQUIRED)
find_package(dune-geometry REQUIRED)

# Find PETSc
find_package(PkgConfig)
pkg_check_modules(PETSC PETSc)
if(NOT PETSC_FOUND)
  # Fallback: manual PETSc path
  set(PETSC_DIR "$ENV{PETSC_DIR}" CACHE PATH "PETSc directory")
  set(PETSC_ARCH "$ENV{PETSC_ARCH}" CACHE STRING "PETSc architecture")
  set(PETSC_INCLUDE_DIRS "${PETSC_DIR}/include" "${PETSC_DIR}/${PETSC_ARCH}/include")
  set(PETSC_LIBRARY_DIRS "${PETSC_DIR}/${PETSC_ARCH}/lib")
  set(PETSC_LIBRARIES "petsc")
endif()

# Find SLEPc
set(SLEPC_DIR "$ENV{SLEPC_DIR}" CACHE PATH "SLEPc directory")
if(SLEPC_DIR)
  set(SLEPC_INCLUDE_DIRS "${SLEPC_DIR}/include" "${SLEPC_DIR}/${PETSC_ARCH}/include")
  set(SLEPC_LIBRARY_DIRS "${SLEPC_DIR}/${PETSC_ARCH}/lib")
  set(SLEPC_LIBRARIES "slepc")
endif()

# Find MPI
find_package(MPI REQUIRED)

add_subdirectory(src)

dune_project()
