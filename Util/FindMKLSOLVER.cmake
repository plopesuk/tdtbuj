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
#    link against to use MKLSOLVER, if BLA_F95 is TRUE f90/f95 interfaces would be added
#  BLA_STATIC  if set on this determines what kind of linkage we do (static)
#  BLA_F95     if set on tries to find the f95 interfaces for BLAS/LAPACK
#  BLA_VENDOR  if set checks only the specified vendor, if not set checks
#     Intel generic
#     List of vendors (BLA_VENDOR) valid in this module
##  Intel10_64lp (mkl mkl_solver_lp64), Intel (generic intel, default if not set)

include(CheckFortranFunctionExists)
set(MKLSOLVER_FOUND FALSE)

macro(Check_MKLSOLVER_Libraries LIBRARIES _prefix _name _flags _list _lapack _threads)
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
   if(BLA_STATIC)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.dll")
    endif(BLA_STATIC)
    find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS ENV LIB
    )
ENDIF (WIN32)

  if(APPLE)
    if(BLA_STATIC)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so;.dylib")
    endif(BLA_STATIC)
    find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV DYLD_LIBRARY_PATH
    )
    else(APPLE)
    if(BLA_STATIC)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
    endif(BLA_STATIC)
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
  set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_lapack} ${_threads})
  #message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
  check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
  set(CMAKE_REQUIRED_LIBRARIES)
  mark_as_advanced(${_prefix}${_combined_name}_WORKS)
  set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endif(_libraries_work)

if(_libraries_work)
  set(${LIBRARIES} ${${LIBRARIES}} ${_lapack})
else(_libraries_work)
  set(${LIBRARIES} FALSE)
endif(_libraries_work)

endmacro(Check_MKLSOLVER_Libraries)


set(MKLSOLVER_LINKER_FLAGS)
set(MKLSOLVER_LIBRARIES)

if(MKLSOLVER_FIND_QUIETLY OR NOT MKLSOLVER_FIND_REQUIRED)
  find_package(LAPACK)
  find_package(Threads)
else(MKLSOLVER_FIND_QUIETLY OR NOT MKLSOLVER_FIND_REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(Threads REQUIRED)
endif(MKLSOLVER_FIND_QUIETLY OR NOT MKLSOLVER_FIND_REQUIRED)

if(LAPACK_FOUND)
  set(MKLSOLVER_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
#intel MKLSOLVER
  if ($ENV{BLA_VENDOR} MATCHES ".+")
    set(BLA_VENDOR $ENV{BLA_VENDOR})
  else ($ENV{BLA_VENDOR} MATCHES ".+")
    if(NOT BLA_VENDOR)
      set(BLA_VENDOR "All")
    endif(NOT BLA_VENDOR)
  endif ($ENV{BLA_VENDOR} MATCHES ".+")
  if (BLA_VENDOR STREQUAL "Intel10_64lp")
    if(BLA_F95)
      if ( NOT MKLSOLVER_LIBRARIES )
        check_MKLSOLVER_libraries(
        MKLSOLVER_LIBRARIES
        MKLSOLVER
        djacobi_solve
        ""
        "mkl_solver_lp64"
        "${LAPACK_LIBRARIES}"
        "${CMAKE_THREAD_LIBS_INIT}"
        )
      endif (NOT MKLSOLVER_LIBRARIES)
    else(BLA_F95)
      if ( NOT MKLSOLVER_LIBRARIES )
        check_MKLSOLVER_libraries(
        MKLSOLVER_LIBRARIES
        MKLSOLVER
        djacobi_solve
        ""
        "mkl_solver_lp64"
        "${LAPACK_LIBRARIES}"
        "${CMAKE_THREAD_LIBS_INIT}"
        )
      endif (NOT MKLSOLVER_LIBRARIES)
    endif(BLA_F95)
  else(BLA_VENDOR STREQUAL "Intel10_64lp")
    if(BLA_F95)
      if ( NOT MKLSOLVER_LIBRARIES )
        check_MKLSOLVER_libraries(
        MKLSOLVER_LIBRARIES
        MKLSOLVER
        djacobi_solve
        ""
        "mkl_solver"
        "${LAPACK_LIBRARIES}"
        "${CMAKE_THREAD_LIBS_INIT}"
        )
      endif (NOT MKLSOLVER_LIBRARIES)
    else(BLA_F95)
      if ( NOT MKLSOLVER_LIBRARIES )
        check_MKLSOLVER_libraries(
        MKLSOLVER_LIBRARIES
        MKLSOLVER
        djacobi_solve
        ""
        "mkl_solver"
        "${LAPACK_LIBRARIES}"
        "${CMAKE_THREAD_LIBS_INIT}"
        )
      endif (NOT MKLSOLVER_LIBRARIES)
    endif(BLA_F95)
  endif(BLA_VENDOR STREQUAL "Intel10_64lp")
else(LAPACK_FOUND)
  message(STATUS "MKLSOLVE requires LAPACK")
endif(LAPACK_FOUND)

if(MKLSOLVER_LIBRARIES)
  set(MKLSOLVER_FOUND TRUE)
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

