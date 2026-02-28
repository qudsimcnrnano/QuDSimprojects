if(NOT dune-p1_FOUND)
# Whether this module is installed or not
set(dune-p1_INSTALLED OFF)

# Settings specific to the module

# Package initialization
# Set prefix to source dir
set(PACKAGE_PREFIX_DIR /home/athira/dune-projects/dune-p1)
macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

#report other information
set_and_check(dune-p1_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-p1_INCLUDE_DIRS "/home/athira/dune-projects/dune-p1")
set(dune-p1_CXX_FLAGS "-std=c++17 /opt/local/petsc-3.15.2/arch-linux-cxx-debug/lib/libpetsc.so.3.15 /opt/local/slepc-3.15.2/arch-linux-cxx-debug/lib/libslepc.so.3.15")
set(dune-p1_CXX_FLAGS_DEBUG "-g")
set(dune-p1_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-p1_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-p1_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-p1_DEPENDS "dune-common;dune-geometry;dune-grid;dune-istl;dune-localfunctions;dune-alugrid;dune-grid-howto;dune-gmsh4")
set(dune-p1_SUGGESTS "")
set(dune-p1_MODULE_PATH "/home/athira/dune-projects/dune-p1/cmake/modules")
set(dune-p1_LIBRARIES "")
set(dune-p1_HASPYTHON 0)
set(dune-p1_PYTHONREQUIRES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-p1_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-p1-targets.cmake")
endif()

endif()
