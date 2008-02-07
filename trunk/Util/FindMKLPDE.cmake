# - Find MKLPDE library
# This module finds an installed fortran library that implements the MKLPDE
# linear-algebra interface (see http://www.netlib.org/MKLPDE/).
#
# The approach follows that taken for the autoconf macro file, acx_MKLPDE.m4
# (distributed at http://ac-archive.sourceforge.net/ac-archive/acx_MKLPDE.html).
#
# This module sets the following variables:
#  MKLPDE_FOUND - set to true if a library implementing the MKLPDEinterface
#    is found
#  MKLPDE_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  MKLPDE_LIBRARIES - uncached list of libraries (using full path name) to 
#    link against to use MKLPDE
 
include(CheckFortranFunctionExists)
set(MKLPDE_FOUND FALSE)

macro(Check_MKLPDE_Libraries LIBRARIES _prefix _name _flags _list _BLAS)
# This macro checks for the existence of the combination of fortran libraries
# given by _list.  If the combination is found, this macro checks (using the 
# Check_Fortran_Function_Exists macro) whether can link against that library
# combination using the name of a routine given by _name using the linker
# flags given by _flags.  If the combination of libraries is found and passes
# the link test, LIBRARIES is set to the list of complete library paths that
# have been found.  Otherwise, LIBRARIES is set to FALSE.
 
# N.B. _prefix is the prefix applied to the names of all cached variables that
# are generated internally and marked advanced by this macro.

set(_libraries_work TRUE)
set(${LIBRARIES})
set(_combined_name)
foreach(_library ${_list})
  set(_combined_name ${_combined_name}_${_library})

  if(_libraries_work)
IF (WIN32)
    find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS ENV LIB 
    )
ENDIF (WIN32)

  if(APPLE)
    find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV DYLD_LIBRARY_PATH
    )
    else(APPLE)
        find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV LD_LIBRARY_PATH
    )
    endif(APPLE)

    mark_as_advanced(${_prefix}_${_library}_LIBRARY)
    set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
    set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
  endif(_libraries_work)
endforeach(_library ${_list})

if(_libraries_work)
  # Test this combination of libraries.
  set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_BLAS})
  #message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
  check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
  set(CMAKE_REQUIRED_LIBRARIES)
  mark_as_advanced(${_prefix}${_combined_name}_WORKS)
  set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endif(_libraries_work)

if(NOT _libraries_work)
  set(${LIBRARIES} FALSE)
endif(NOT _libraries_work)

endmacro(Check_MKLPDE_Libraries)


set(MKLPDE_LINKER_FLAGS)
set(MKLPDE_LIBRARIES)

if(MKLPDE_FIND_QUIETLY OR NOT MKLPDE_FIND_REQUIRED)
  find_package(BLAS)
else(MKLPDE_FIND_QUIETLY OR NOT MKLPDE_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
endif(MKLPDE_FIND_QUIETLY OR NOT MKLPDE_FIND_REQUIRED)

if(BLAS_FOUND)
  set(MKLPDE_LINKER_FLAGS ${BLAS_LINKER_FLAGS})

#intel MKLPDE
if ( NOT MKLPDE_LIBRARIES )
    check_MKLPDE_libraries(
    MKLPDE_LIBRARIES
    MKLPDE
    cheev
    ""
    "mkl_lapack95;mkl_lapack;mkl"
    "${BLAS_LIBRARIES}"
    )
  endif ( NOT MKLPDE_LIBRARIES )

else(BLAS_FOUND)
  message(STATUS "MKLPDErequires BLAS")
endif(BLAS_FOUND)

if(MKLPDE_LIBRARIES)
  set(MKLPDE_FOUND TRUE)
else(MKLPDE_LIBRARIES)
  set(MKLPDE_FOUND FALSE)
endif(MKLPDE_LIBRARIES)

if(NOT MKLPDE_FIND_QUIETLY)
  if(MKLPDE_FOUND)
    message(STATUS "A library with MKLPDEAPI found.")
  else(MKLPDE_FOUND)
    if(MKLPDE_FIND_REQUIRED)
      message(FATAL_ERROR 
      "A required library with MKLPDEAPI not found. Please specify library location."
      )
    else(MKLPDE_FIND_REQUIRED)
      message(STATUS
      "A library with MKLPDEAPI not found. Please specify library location."
      )
    endif(MKLPDE_FIND_REQUIRED)
  endif(MKLPDE_FOUND)
endif(NOT MKLPDE_FIND_QUIETLY)
