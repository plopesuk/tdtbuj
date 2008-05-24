# - Find MKLSOLVER library
# This module finds an installed fortran library that implements the MKLSOLVER
# linear-algebra interface (see http://www.netlib.org/MKLSOLVER/).
#
# The approach follows that taken for the autoconf macro file, acx_MKLSOLVER.m4
# (distributed at http://ac-archive.sourceforge.net/ac-archive/acx_MKLSOLVER.html).
#
# This module sets the following variables:
#  MKLSOLVER_FOUND - set to true if a library implementing the MKLSOLVERinterface
#    is found
#  MKLSOLVER_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  MKLSOLVER_LIBRARIES - uncached list of libraries (using full path name) to 
#    link against to use MKLSOLVER
 
include(CheckFortranFunctionExists)
set(MKLSOLVER_FOUND FALSE)

macro(Check_MKLSOLVER_Libraries LIBRARIES _prefix _name _flags _list _BLAS)
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

endmacro(Check_MKLSOLVER_Libraries)


set(MKLSOLVER_LINKER_FLAGS)
set(MKLSOLVER_LIBRARIES)

if(MKLSOLVER_FIND_QUIETLY OR NOT MKLSOLVER_FIND_REQUIRED)
  find_package(BLAS)
  find_package(LAPACK)
else(MKLSOLVER_FIND_QUIETLY OR NOT MKLSOLVER_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif(MKLSOLVER_FIND_QUIETLY OR NOT MKLSOLVER_FIND_REQUIRED)

if(BLAS_FOUND)
  set(MKLSOLVER_LINKER_FLAGS ${BLAS_LINKER_FLAGS})

#intel MKLSOLVER

IF (WIN32)
if ( NOT MKLSOLVER_LIBRARIES )
    check_MKLSOLVER_libraries(
    MKLSOLVER_LIBRARIES
    MKLSOLVER
    djacobi_solve
    ""
    "mkl_solver"
    "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}"
    )
  endif (NOT MKLSOLVER_LIBRARIES)
ELSE (WIN32)
if ( NOT MKLSOLVER_LIBRARIES )
    check_MKLSOLVER_libraries(
    MKLSOLVER_LIBRARIES
    MKLSOLVER
    djacobi_solve
    ""
    "mkl_solver_lp64"
    "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}"
    )
  endif (NOT MKLSOLVER_LIBRARIES)
if ( NOT MKLSOLVER_LIBRARIES )
    check_MKLSOLVER_libraries(
    MKLSOLVER_LIBRARIES
    MKLSOLVER
    djacobi_solve
    ""
    "mkl_solver"
    "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}"
    )
  endif (NOT MKLSOLVER_LIBRARIES)



ENDIF (WIN32)
else(BLAS_FOUND)
  message(STATUS "MKLSOLVE Rrequires BLAS")
endif(BLAS_FOUND)

if(MKLSOLVER_LIBRARIES)
  set(MKLSOLVER_FOUND TRUE)
else(MKLSOLVER_LIBRARIES)
  set(MKLSOLVER_FOUND FALSE)
endif(MKLSOLVER_LIBRARIES)

if(NOT MKLSOLVER_FIND_QUIETLY)
  if(MKLSOLVER_FOUND)
    message(STATUS "A library with MKLSOLVER API found.")
  else(MKLSOLVER_FOUND)
    if(MKLSOLVER_FIND_REQUIRED)
      message(FATAL_ERROR 
      "A required library with MKLSOLVER API not found. Please specify library location."
      )
    else(MKLSOLVER_FIND_REQUIRED)
      message(STATUS
      "A library with MKLSOLVER API not found. Please specify library location."
      )
    endif(MKLSOLVER_FIND_REQUIRED)
  endif(MKLSOLVER_FOUND)
endif(NOT MKLSOLVER_FIND_QUIETLY)
