PROJECT(tdtbuj Fortran)

set(CMAKE_BUILD_TYPE "Debug")

set(VERSION_MAJOR "0")
set(VERSION_MINOR "9")
set(VERSION_PATCH "1")


find_program(SVN "svn")
if(SVN)
  execute_process(COMMAND ${SVN} info ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE SVN_OUT)
  string(REGEX REPLACE "\n" " " svn_int "${SVN_OUT}")
  string(REGEX REPLACE ".*Revision: ([^ ]+).*" "\\1" VERSION_PATCH "${svn_int}" )
  set(VERSION_PATCH "svn_${VERSION_PATCH}")
endif(SVN)

set(AUTHOR "Alin M Elena")
set(AUTHOR_DETAILS "alin.elena@qub.ac.uk")
set(DESCRITION "Time dependent tight binding +UJ for finite systems")

cmake_minimum_required(VERSION 2.6.0)
if(COMMAND cmake_policy)
  cmake_policy(VERSION 2.6)
endif(COMMAND cmake_policy)

##########################
###source files
set(src_dir source)
set(target_name tdtbuj)
##########################
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Util")

##########################
## external libs
FIND_PACKAGE(BLAS)
MESSAGE("Using this version of blas95 ${BLAS95_LIBRARIES} with flags
${BLAS_LINKER_FLAGS}")

FIND_PACKAGE(LAPACK)
MESSAGE("Using this version of lapack95 ${LAPACK95_LIBRARIES} with flags
${LAPACK_LINKER_FLAGS}")

find_package(MKLSOLVER)
message("Using this version of solver: ${MKLSOLVER_LIBRARIES} with flags  ${MKLSOLVER_LINKER_FLAGS} ")
set(LIBS ${MKLSOLVER_LINKER_FLAGS} ${LAPACK95_LIBRARIES} ${MKLSOLVER_LIBRARIES} ${BLAS95_LIBRARIES} )
##########################


######################
###documentation
set(doc_dir doc)
find_package(Doxygen)
######################

##set the output folder for libs and bins
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
###################


######
## adds user libs

add_library(parser
${src_dir}/Constants.f90
${src_dir}/Parser.f90
${src_dir}/Types.f90
${src_dir}/Useful.f90
)

add_library(linalg
${src_dir}/Constants.f90
${src_dir}/Types.f90
${src_dir}/Useful.f90
${src_dir}/LinearAlgebra.F90
)

####################
#~ set_source_files_properties(${src_dir}/SCF2.f90
#~ 	PROPERTIES COMPILE_FLAGS /Od)


####################
##### adds the main target
add_executable(${target_name}
${src_dir}/Constants.f90
${src_dir}/Types.f90
${src_dir}/Useful.f90
${src_dir}/Gutenberg.f90
${src_dir}/ReadData.f90
${src_dir}/TailFunctions.f90
${src_dir}/TightBinding.f90
${src_dir}/Electrostatics.f90
${src_dir}/SlaterKoster.f90
${src_dir}/DensityMatrix.f90
${src_dir}/Hamiltonian.f90
${src_dir}/Mixing.f90
${src_dir}/SCF.f90
${src_dir}/Simplex.f90
${src_dir}/SA.f90
${src_dir}/SimplexSA.f90
${src_dir}/Fit.f90
${src_dir}/Dynamics.f90
${src_dir}/LBFGS.f90
${src_dir}/BFGS.f90
${src_dir}/DriverRoutines.f90
${src_dir}/Testing.f90
${src_dir}/main.f90
)
####################

#################
## sets the linking
link_directories (${LIBRARY_OUTPUT_PATH})
target_link_libraries(${target_name} parser linalg ${LIBS})
#################

if(DOXYGEN)
  configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.cmake Doxyfile)
  add_custom_target(docs
    ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile)
endif(DOXYGEN)

####################
####install target
####################

if(NOT ${PREFIX})
  MESSAGE("Install folder ${PREFIX}")
  set(CMAKE_INSTALL_PREFIX "${PREFIX}")
endif(NOT ${PREFIX})

install(TARGETS ${target_name} parser linalg 
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib/${target_name}
	    ARCHIVE DESTINATION lib/${target_name}
	)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/examples
        DESTINATION share
        REGEX ".svn" EXCLUDE)
if(DOXYGEN)
  install(CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} docs)")
  install(DIRECTORY ${CMAKE_BINARY_DIR}/doc/html/ DESTINATION doc)
endif(DOXYGEN)
#######################
## uninstall target
#######################

###create the uninstall script
set(unfile ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
file(WRITE ${unfile} "IF(NOT EXISTS \"install_manifest.txt\")\n")
file(APPEND ${unfile} "MESSAGE(\"FATAL_ERROR Cannot find \\\"install manifest\\\": install_manifest.txt\")\n")
file(APPEND ${unfile} "ENDIF(NOT EXISTS \"install_manifest.txt\")\n")
file(APPEND ${unfile} "FILE(READ \"install_manifest.txt\" files)\n")
file(APPEND ${unfile} "STRING(REGEX REPLACE \"\\n\" \";\" files \"\${files}\")\n")
file(APPEND ${unfile} "FOREACH(file \${files})\n")
file(APPEND ${unfile} "  MESSAGE(STATUS \"Uninstalling \\\"\${file}\\\"\")\n")
file(APPEND ${unfile} "  IF(EXISTS \"\${file}\")\n")
file(APPEND ${unfile} "    EXEC_PROGRAM(\n")
file(APPEND ${unfile} "      \"\${CMAKE_COMMAND}\" ARGS \"-E remove \\\"\${file}\\\"\"\n")
file(APPEND ${unfile} "    OUTPUT_VARIABLE rm_out\n")
file(APPEND ${unfile} "      RETURN_VALUE rm_retval\n")
file(APPEND ${unfile} "      )\n")
file(APPEND ${unfile} "    IF(\"\${rm_retval}\" STREQUAL 0\)\n")
file(APPEND ${unfile} "    ELSE(\"\${rm_retval}\" STREQUAL 0\)\n")
file(APPEND ${unfile} "    MESSAGE(FATAL_ERROR \"Problem when removing \\\"\${file}\\\"\")\n")
file(APPEND ${unfile} "    ENDIF(\"\${rm_retval}\" STREQUAL 0)\n")
file(APPEND ${unfile} "  ELSE(EXISTS \"\${file}\")\n")
file(APPEND ${unfile} "    MESSAGE(STATUS \"File \\\"\${file}\\\" does not exist.\")\n")
file(APPEND ${unfile} "  ENDIF(EXISTS \"\${file}\")\n")
file(APPEND ${unfile} "ENDFOREACH(file)\n")
########################

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
###########################

###########################
##Packaging
###########################

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${DESCRIPTION})
SET(CPACK_PACKAGE_EXECUTABLES ${target_name})
SET(CPACK_PACKAGE_VENDOR ${AUTHOR})
SET(CPACK_PACKAGE_CONTACT ${AUTHOR_DETAILS})
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${target_name}-${VERSION_MAJOR}.${VERSION_MINOR}")
SET(CPACK_GENERATOR "TGZ;TBZ2;ZIP")
SET(CPACK_CMAKE_GENERATOR ${CMAKE_GENERATOR})
SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2;ZIP")
SET(CPACK_SOURCE_IGNORE_FILES "/thesis/;/\\\\.svn/")

INCLUDE(CPack)
