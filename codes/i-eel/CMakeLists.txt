#
# This file is a part of the i-eel project
#
#  Include and configure this version of feedforward.
#
cmake_minimum_required(VERSION 3.1)
cmake_policy(SET CMP0057 NEW)
project(feedforward)

######## Feedforward version ########
set(IEEL_VERSION_MAJOR "0")
set(IEEL_VERSION_MINOR "0")
set(IEEL_VERSION_PATCH "0")
set(IEEL_VERSION "${IEEL_VERSION_MAJOR}.${IEEL_VERSION_MINOR}.${IEEL_VERSION_PATCH}")

message(STATUS "")
message(STATUS "configuring i-eel version " ${IEEL_VERSION})

######## Command line options ########
OPTION(BUILD_SHARED_LIBS "Builds shared libraries if ON, static libraries if OFF" ON)
OPTION(SET_RPATH "set the rpath of the executable to the found libraries" ON)
OPTION(WITH_HDF5CXX "enable legacy file format using hdf5 cxx" OFF)

# Set compiler flags that are build-type related
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

######## Basic settings ########

# Ask for C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message(STATUS "Using GNU compiler options")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(STATUS "Using Apple compiler options")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL Intel)
    message(STATUS "Using Intel compiler options")
endif ()
set(COMPILER_VERSION "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# rpath settings
# From: http://www.itk.org/Wiki/CMake_RPATH_handling
if (SET_RPATH)
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    if ("${isSystemDir}" STREQUAL "-1")
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif ("${isSystemDir}" STREQUAL "-1")
else ()
    message(STATUS "not setting any rpath")
    set(CMAKE_SKIP_RPATH TRUE)
endif ()

######## Library finding and system introspection ########

message(STATUS "")
message(STATUS "Starting system introspection")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})

if (BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so)
else ()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .a)
endif ()

# Basic system checks
include(CheckLibraryExists)
include(CheckFunctionExists)
include(CheckIncludeFiles)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compiler_helper.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/target_helper.cmake)


# Find all the dependencies

# NetCDF
find_package(NetCDF REQUIRED)
message(STATUS "NetCDF include dir: ${NETCDF_INCLUDE_DIRS}")
set(CMAKE_REQUIRED_INCLUDES ${NETCDF_INCLUDE_DIRS})
CHECK_INCLUDE_FILES(netcdf.h HAVE_NETCDF_H)
if (NOT HAVE_NETCDF_H)
    message(FATAL ERROR "could not include netcdf.h")
endif ()

# Armadillo
find_package(ARMA REQUIRED)
set(CMAKE_REQUIRED_INCLUDES ${ARMA_INCLUDE_DIR})

######## Create config.h from config.h.in ########
configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/ieel/config.h")

# Set up bundled projects
add_subdirectory(basics)
add_subdirectory(src)
add_subdirectory(programs)

######## Summary ########
message("\n")
message("   ###############################################")
message("   ###########  Configuration summary  ###########")
message("   ###############################################")
message("   Compiler:                       ${COMPILER_VERSION}")
message("   Build type:                     ${CMAKE_BUILD_TYPE}")
message("   Install prefix:                 ${CMAKE_INSTALL_PREFIX}")

if (BUILD_SHARED_LIBS)
    message("   Building shared library:        yes")
    message("   Linking programs to:            shared library")
else ()
    message("   Building shared library:        no")
    message("   Linking programs to:            static library")
endif ()

message("\n   Libraries")
message("   netcdf (default file format):   enabled")
message("\n")

