# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.18)

if (CMAKE_VERSION VERSION_LESS "3.1")
  if (NOT (${CMAKE_VERSION} VERSION_LESS "2.8.12"))
      add_compile_options("-std=c++11")
  else ()
      add_definitions("-std=c++11")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if ("${CMAKE_DataType}" STREQUAL "double")
  set(USE_QUADMATH FALSE)
  add_definitions(-D DATATYPE_DOUBLE)
endif ()
if ("${CMAKE_DataType}" STREQUAL "float128")
  set(USE_QUADMATH TRUE)
  add_definitions(-D DATATYPE_FLOAT128)
endif ()
if ("${CMAKE_DataType}" STREQUAL "custom")
  set(USE_QUADMATH FALSE)
  add_definitions(-D DATATYPE_CUSTOM)
endif ()

# Project's name
project(claw1dArena)
set (claw1dArena_VERSION_MAJOR 1)
set (claw1dArena_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
"${PROJECT_SOURCE_DIR}/Claw1dArena_config.h.in"
"${PROJECT_BINARY_DIR}/Claw1dArena_config.h"
)
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

add_executable(testrec src/testReconstructions)
target_compile_definitions(testrec PUBLIC DEBUG_REC)

add_executable(lintra src/claw1d)
target_compile_definitions(lintra PUBLIC LINTRA)
add_executable(burgers src/claw1d)
target_compile_definitions(burgers PUBLIC BURGERS)
#add_executable(gasdynamics src/claw1d)
#target_compile_definitions(gasdynamics PUBLIC GASDIN)
add_executable(swe src/claw1d)
target_compile_definitions(swe PUBLIC SWE)


configure_file(src/python/graphit.py python/graphit.py COPYONLY)
configure_file(src/python/graphitGas.py python/graphitGas.py COPYONLY)
configure_file(src/python/graphitSWE.py python/graphitSWE.py COPYONLY)

find_package( Boost REQUIRED COMPONENTS program_options)
include_directories( ${Boost_INCLUDE_DIRS} )
set(EXTRA_LIBS ${Boost_LIBRARIES})

find_package(QUADMATH)
list(APPEND EXTRA_LIBS ${QUADMATH_LIBRARIES})

target_link_libraries(testrec ${EXTRA_LIBS})
target_link_libraries(lintra ${EXTRA_LIBS})
target_link_libraries(burgers ${EXTRA_LIBS})
#target_link_libraries(gasdynamics ${EXTRA_LIBS})
target_link_libraries(swe ${EXTRA_LIBS})

#LAPACKE is used for characteristic reconstructions in gasdynamics
set(LAPACKE_STATIC 0 CACHE BOOL "Require LAPACKE static libraries")
#mark_as_advanced(LAPACKE_STATIC)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
#target_link_libraries(gasdynamics ${LAPACKE_LIBRARIES} ${LAPACK_LIBRARIES})

#utils library
add_subdirectory("${PROJECT_SOURCE_DIR}/src/utils")
target_link_libraries(testrec utils)
target_link_libraries(lintra utils)
target_link_libraries(burgers utils)
#target_link_libraries(gasdynamics utils)
target_link_libraries(swe utils)

#dof library
add_subdirectory("${PROJECT_SOURCE_DIR}/src/dof")
target_link_libraries(testrec dof)
target_link_libraries(lintra dof)
target_link_libraries(burgers dof)
#target_link_libraries(gasdynamics dof)
target_link_libraries(swe dof)

#grid librrary
add_subdirectory("${PROJECT_SOURCE_DIR}/src/grid")
target_link_libraries(testrec grid)
target_link_libraries(lintra grid)
target_link_libraries(burgers grid)
#target_link_libraries(gasdynamics grid)
target_link_libraries(swe grid)

#reconstructions libraries
FILE(GLOB REC_SOURCES ${PROJECT_SOURCE_DIR}/src/reconstructions/*)
add_library(reconstructions       ${REC_SOURCES})
add_library(reconstructions_debug ${REC_SOURCES})
target_compile_definitions(reconstructions_debug PRIVATE DEBUG_REC)
target_link_libraries(testrec reconstructions_debug)
target_link_libraries(lintra reconstructions)
target_link_libraries(burgers reconstructions)
#target_link_libraries(gasdynamics reconstructions)
target_link_libraries(swe reconstructions)

#claws library
add_subdirectory("${PROJECT_SOURCE_DIR}/src/claws")
target_link_libraries(testrec claws)
target_link_libraries(lintra claws)
target_link_libraries(burgers claws)
#target_link_libraries(gasdynamics claws)
target_link_libraries(swe claws)

#numfluxes library
add_subdirectory("${PROJECT_SOURCE_DIR}/src/numfluxes")
target_link_libraries(lintra numfluxes)
target_link_libraries(burgers numfluxes)
#target_link_libraries(gasdynamics numfluxes)
target_link_libraries(swe numfluxes)

#timeintegration library
add_subdirectory("${PROJECT_SOURCE_DIR}/src/timeintegration")
target_link_libraries(lintra timeintegration)
target_link_libraries(burgers timeintegration)
#target_link_libraries(gasdynamics timeintegration)
target_link_libraries(swe timeintegration)

if (BUILD_DOC)
  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  else()
    set(DOXYFILE_IN doc/Doxyfile.in)
    set(DOXYFILE doc/Doxyfile)
    configure_file(${DOXYFILE_IN} ${DOXYFILE} )
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${DOXYFILE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()
endif()

#FILE(GLOB CLAWS_SOURCES ${PROJECT_SOURCE_DIR}/src/claws/)
FILE(GLOB UTILS_SOURCES ${PROJECT_SOURCE_DIR}/src/utils/*)
FILE(GLOB DOF_SOURCES ${PROJECT_SOURCE_DIR}/src/dof/*)
FILE(GLOB GRID_SOURCES ${PROJECT_SOURCE_DIR}/src/grid/*)
FILE(GLOB NUMFLUX_SOURCES ${PROJECT_SOURCE_DIR}/src/numfluxes/*)
FILE(GLOB TIME_SOURCES ${PROJECT_SOURCE_DIR}/src/timeintegration/*)
set(ALL_SOURCES
#  ${CLAWS_SOURCES}
  ${UTILS_SOURCES} ${DOF_SOURCES} ${GRID_SOURCES} ${REC_SOURCES} ${NUMFLUX_SOURCES} ${TIME_SOURCES})
add_library(shared_claw1darena SHARED ${ALL_SOURCES})

add_executable(test src/claw1d_test.cpp)
target_link_libraries(test    ${EXTRA_LIBS}  utils dof grid reconstructions claws numfluxes timeintegration )

#------------------------------------------------------------------------------------#
# Variables visible to the user
#------------------------------------------------------------------------------------#
set(VERBOSE_MAKE 0 CACHE BOOL "Set appropriate compiler and cmake flags to enable verbose output from compilation")

set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
	set(CMAKE_VERBOSE_MAKEFILE 1)
endif()

# Set build type variable
if(NOT CMAKE_BUILD_TYPE)
	# Set default build type to Debug
	set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
		"Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel GNUProfiling ScalascaProfiling."
		FORCE)

	# Set the possible values of build type for the GUI
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
	"MinSizeRel" "RelWithDebInfo" "GNUProfiling" "ScalascaProfiling")
endif()

#build documentation?
set(BUILD_DOC 0 CACHE BOOL "Build documentation?")

#datatype to use
set(CMAKE_DataType "double" CACHE STRING "datatype for variables")
set_property(CACHE CMAKE_DataType PROPERTY STRINGS double float128 custom)
