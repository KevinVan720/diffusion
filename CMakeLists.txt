cmake_minimum_required (VERSION 2.6)
project (diffusion Fortran)

# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
  set(CMAKE_INSTALL_PREFIX_ALT "/usr/lib64"
      CACHE PATH "Install another path prefix"
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX} ${CMAKE_INSTALL_PREFIX_ALT}")

find_library(HDF5_FORTRAN_LIB hdf5_fortran)
if("${HDF5_FORTRAN_LIB}" STREQUAL "HDF5_FORTRAN_LIB-NOTFOUND")
  message(STATUS "HDF5 FORTRAN library not found")
endif()

find_path(HDF5_INCLUDE_DIRS hdf5.mod)
if("${HDF5_INCLUDE_DIRS}" STREQUAL "HDF5_INCLUDE_DIRS-NOTFOUND")
  message(STATUS "HDF5 FORTRAN include not found")
  set(HDF5_INCLUDE_DIRS "/usr/lib64/gfortran/modules/")
endif()

include_directories(${HDF5_INCLUDE_DIRS})

# set default architecture flags to optimize fot the furrent system,
# but don't overide user settings
if(NOT CMAKE_Fortran_FLAGS MATCHES "-march|-mtune")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fno-align-commons")

add_subdirectory (src)

# install dependencies
install (FILES gamma-table.dat parameters_df.dat dNg_over_dt_cD6.dat dNg_over_dt_bD6.dat DESTINATION ${PROJECT_NAME})

