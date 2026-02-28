cmake_minimum_required(VERSION 3.13)
project(dune-p3d CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)
SET ( CMAKE_CXX_FLAGS "/opt/local/petsc-3.15.2/arch-linux-cxx-debug/lib/libpetsc.so.3.15 /opt/local/slepc-3.15.2/arch-linux-cxx-debug/lib/libslepc.so.3.15")

# start a dune project with information from dune.module
dune_project()
find_package(MPI)
find_package(PkgConfig)
pkg_search_module(PETSC REQUIRED PETSc)
pkg_search_module(SLEPC REQUIRED slepc)
find_package(PETSc)
find_package(slepc)

dune_enable_all_packages()

add_subdirectory(src)
add_subdirectory(dune)
add_subdirectory(doc)
add_subdirectory(cmake/modules)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
