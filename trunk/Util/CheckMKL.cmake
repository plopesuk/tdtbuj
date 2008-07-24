# - Find Intel mkl library
# This module finds an installed intel mathematical kernel library


get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if(NOT _LANGUAGES_ MATCHES Fortran)
  if(LAPACK_FIND_REQUIRED)
    message(FATAL_ERROR
      "FindLAPACK is Fortran-only so Fortran must be enabled.")
  else(LAPACK_FIND_REQUIRED)
    message(STATUS "Looking for LAPACK... - NOT found (Fortran not enabled)")
    return()
  endif(LAPACK_FIND_REQUIRED)
endif(NOT _LANGUAGES_ MATCHES Fortran)

include(CheckFortranFunctionExists)

macro(Check_MKL_Libraries LIBRARIES _prefix _name _flags _other _clusterLayer _threadingLayer _computationalLayer _guideLayer _threads _aux)
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
set(_list "${_other}")
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

if(_libraries_work AND BLA_STATIC)
  set(${LIBRARIES} ${${LIBRARIES}} ";-Wl,--start-group")
endif(_libraries_work AND BLA_STATIC)
set(_list "${_clusterLayer};${_threadingLayer};${_computationalLayer}")
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
if(_libraries_work AND BLA_STATIC)
  set(${LIBRARIES} ${${LIBRARIES}} ";-Wl,--end-group")
endif(_libraries_work AND BLA_STATIC)

foreach(_library ${_guideLayer})
  set(_combined_name ${_combined_name}_${_library})
  if(_libraries_work)
    IF (WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll,.lib")
      find_library(${_prefix}_${_library}_LIBRARY
      NAMES ${_library}
      PATHS ENV LIB
      )
    ENDIF (WIN32)
    if(APPLE)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.dylib;.a")
      find_library(${_prefix}_${_library}_LIBRARY
      NAMES ${_library}
      PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV DYLD_LIBRARY_PATH
      )
    else(APPLE)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
      find_library(${_prefix}_${_library}_LIBRARY
      NAMES ${_library}
      PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV LD_LIBRARY_PATH
      )
    endif(APPLE)
    mark_as_advanced(${_prefix}_${_library}_LIBRARY)
    set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
    set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
  endif(_libraries_work)
endforeach(_library ${_guideLayer})

if(_libraries_work)
  # Test this combination of libraries.
  set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_threads} ${_aux})
  #message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
  check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
  mark_as_advanced(${_prefix}${_combined_name}_WORKS)
  set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endif(_libraries_work)

 if(_libraries_work)
     set(${LIBRARIES} ${${LIBRARIES}})
 else(_libraries_work)
    set(${LIBRARIES} FALSE)
 endif(_libraries_work)
 set(CMAKE_REQUIRED_LIBRARIES)
endmacro(Check_MKL_Libraries)
