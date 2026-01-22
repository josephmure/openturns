# FindRanger.cmake
# -----------------
#
# Locate the Ranger C++ random forest library
#
# This module defines:
#
# ::
#
#   RANGER_FOUND - System has Ranger
#   RANGER_INCLUDE_DIRS - The Ranger include directories
#   RANGER_LIBRARIES - The libraries needed to use Ranger
#   RANGER_VERSION - The version of Ranger found
#
# and the following imported target:
#
# ::
#
#   Ranger::Ranger - The Ranger library
#
# You can set the following variables to help guide the search:
#
# ::
#
#   RANGER_ROOT_DIR - Root directory of Ranger installation
#   RANGER_INCLUDE_DIR - Directory containing Ranger headers
#   RANGER_LIBRARY_DIR - Directory containing Ranger library

#=============================================================================
# Copyright 2005-2024 Airbus-EDF-IMACS-ONERA-Phimeca
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# First, check for user-specified paths
if(RANGER_ROOT_DIR)
  set(_RANGER_SEARCH_OPTS NO_DEFAULT_PATH)
else()
  set(_RANGER_SEARCH_OPTS)
endif()

# Find the main Ranger header file
find_path(RANGER_INCLUDE_DIR
  NAMES Forest.h
  HINTS
    ${RANGER_ROOT_DIR}
    $ENV{RANGER_ROOT_DIR}
    ${RANGER_ROOT_DIR}/src
    $ENV{RANGER_ROOT_DIR}/src
  PATH_SUFFIXES
    include
    include/ranger
    ranger
    src
  ${_RANGER_SEARCH_OPTS}
  DOC "Directory containing Ranger header files"
)

# Find the Ranger library
find_library(RANGER_LIBRARY
  NAMES ranger libranger
  HINTS
    ${RANGER_ROOT_DIR}
    $ENV{RANGER_ROOT_DIR}
    ${RANGER_LIBRARY_DIR}
    $ENV{RANGER_LIBRARY_DIR}
  PATH_SUFFIXES
    lib
    lib64
    build
    cpp_version/build
  ${_RANGER_SEARCH_OPTS}
  DOC "Ranger library"
)

# Extract version information if possible
if(RANGER_INCLUDE_DIR)
  # Try to find version from parent directory if it's a git repo
  get_filename_component(_ranger_parent "${RANGER_INCLUDE_DIR}" DIRECTORY)
  if(EXISTS "${_ranger_parent}/.git")
    find_package(Git QUIET)
    if(Git_FOUND)
      execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
        WORKING_DIRECTORY ${_ranger_parent}
        OUTPUT_VARIABLE RANGER_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
    endif()
  endif()
  
  # Fallback to a default version
  if(NOT RANGER_VERSION)
    set(RANGER_VERSION "0.16.0")
  endif()
endif()

# Verify that we found all required headers
set(_RANGER_REQUIRED_HEADERS
  Forest.h
  ForestClassification.h
  ForestRegression.h
  ForestSurvival.h
  Data.h
  Tree.h
)

set(_RANGER_HEADERS_FOUND TRUE)
if(RANGER_INCLUDE_DIR)
  foreach(_header ${_RANGER_REQUIRED_HEADERS})
    if(NOT EXISTS "${RANGER_INCLUDE_DIR}/${_header}")
      set(_RANGER_HEADERS_FOUND FALSE)
      if(Ranger_FIND_REQUIRED)
        message(STATUS "Missing required Ranger header: ${_header}")
      endif()
    endif()
  endforeach()
endif()

# Handle the QUIETLY and REQUIRED arguments and set RANGER_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Ranger
  FOUND_VAR RANGER_FOUND
  REQUIRED_VARS RANGER_LIBRARY RANGER_INCLUDE_DIR _RANGER_HEADERS_FOUND
  VERSION_VAR RANGER_VERSION
  FAIL_MESSAGE "Could not find Ranger C++ library. Set RANGER_ROOT_DIR to specify the installation directory, or build Ranger and set RANGER_LIBRARY_DIR to the build directory."
)

if(RANGER_FOUND)
  set(RANGER_INCLUDE_DIRS ${RANGER_INCLUDE_DIR})
  set(RANGER_LIBRARIES ${RANGER_LIBRARY})
  
  # Create imported target
  if(NOT TARGET Ranger::Ranger)
    add_library(Ranger::Ranger UNKNOWN IMPORTED)
    set_target_properties(Ranger::Ranger PROPERTIES
      IMPORTED_LOCATION "${RANGER_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${RANGER_INCLUDE_DIRS}"
    )
  endif()
  
  # Mark as advanced
  mark_as_advanced(RANGER_INCLUDE_DIR RANGER_LIBRARY)
  
  if(NOT Ranger_FIND_QUIETLY)
    message(STATUS "Found Ranger: ${RANGER_LIBRARY} (version ${RANGER_VERSION})")
  endif()
endif()
