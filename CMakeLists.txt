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
#  set(CMAKE_INSTALL_PREFIX_ALT "/usr/lib64"
#      CACHE PATH "Install another path prefix"
#      FORCE)
endif()



message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
#message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX} ${CMAKE_INSTALL_PREFIX_ALT}")

find_package(HDF5 COMPONENTS Fortran)
if(HDF5_FOUND)
  message(STATUS "hdf5 found, using package")
  message(STATUS ${HDF5_INCLUDE_DIRS})
  message(STATUS ${HDF5_Fortran_LIBRARIES})
  message(STATUS ${HDF5_LIBRARIES})
  include_directories(${HDF5_INCLUDE_DIRS})
else()
  message(STATUS "hdf5 not found")
endif()


# set default architecture flags to optimize fot the furrent system,
# but don't overide user settings
if(NOT CMAKE_Fortran_FLAGS MATCHES "-march|-mtune")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fno-align-commons")

add_subdirectory (src)

# install dependencies
#install (FILES gamma-table_charm.dat  gamma-table_bottom.dat parameters_df.dat dNg_over_dt_cD6.dat dNg_over_dt_bD6.dat DESTINATION ${PROJECT_NAME})
#install (FILES parameters_df.dat dNg_over_dt_cD6.dat dNg_over_dt_bD6.dat DESTINATION share)
install (FILES dNg_over_dt_cD6.dat dNg_over_dt_bD6.dat DESTINATION share)
install (FILES parameters_df.dat DESTINATION results)
