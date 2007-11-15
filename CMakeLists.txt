PROJECT(tdtbuj Fortran)

##########################
###source files

set(src_dir source)
set(target_name tdtbuj.x)
##########################

######
## external libs
FIND_PACKAGE(BLAS)
MESSAGE("Using this version of blas ${BLAS_LIBRARIES} with flags
${BLAS_LINKER_FLAGS}")

find_package(LAPACK)
message("Using this version of lapack: ${LAPACK_LIBRARIES} with flags  ${LAPACK_LINKER_FLAGS} ")


set(LIBS ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
########



####
###documentation
set(doc_dir doc)

find_package(Doxygen)

find_file(LATEX
NAMES pdflatex
PATHS ENV PATH
)
MESSAGE("Using this version of pdflatex ${LATEX}")

find_file(MKIDX
NAMES makeindex
PATHS ENV PATH
)

MESSAGE("Using this version of makeindex ${MKIDX}")
####


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
${src_dir}/LinearAlgebra.f90
)

####################

#################
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
${src_dir}/DriverRoutines.f90
${src_dir}/main.f90
)
####################

############
## sets the linking
link_directories (${LIBRARY_OUTPUT_PATH})
target_link_libraries(${target_name} parser linalg ${LIBS})
#################


add_custom_target(docs
${DOXYGEN_EXECUTABLE}
)


######
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

#######################
## uninstall target
##################

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
