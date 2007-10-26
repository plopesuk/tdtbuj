!****d* Module/constants
! NAME
! constants
! SYNOPSIS
! use constants
! DESCRIPTION
! This is the module that defines all the constants used in the program
! AUTHOR
! Alin M.Elena (Queen's University Belfast)
! CREATION DATE
! 14th of January, 2006
!****** 


module constants
  implicit none
  private
  
! define precision for reals
  integer, parameter, public :: pr=kind(1.0d0)
! number of character per line
  integer, parameter, public :: ml=255 
! number of character per word
  integer, parameter, public :: mw=40 
! verbosity levels
  integer, parameter,public :: low_verbos=5
  integer, parameter,public :: medium_verbos=15
  integer, parameter,public :: high_verbos=25
! debug levels
  integer, parameter,public :: low_debug=5
  integer, parameter,public :: medium_debug=15
  integer, parameter,public :: high_debug=25

!!!!!!!!!!



! RunType
 integer, parameter, public :: run_sp = 1
 integer, parameter, public :: run_bo = 2
 integer, parameter, public :: run_ehrenfest = 3
 integer, parameter, public :: run_fit = 4
 integer, parameter, public :: run_force_test = 5
 integer, parameter, public :: run_force_testx = 6
 integer, parameter, public :: run_force_testy = 7
 integer, parameter, public :: run_force_testz = 8
 integer, parameter, public :: run_ehrenfest_damped = 9 

! SCFType
 integer, parameter, public :: scf_tbuj = 1

!SCFMixType
integer, parameter, public :: scfmix_broyden = 1
integer, parameter, public :: scfmix_pulay = 2

!electrostatics
integer, parameter, public :: electrostatics_point=1
integer, parameter, public :: electrostatics_multipoles=2

!spins
  integer, parameter, public :: spin_down=0
  integer, parameter, public :: spin_up=1
! bond type 
  integer, parameter, public :: bond_gsp=1
  integer, parameter, public :: bond_harrison=2

! smearing methods
  integer, parameter, public :: sm_fd=1
  integer, parameter, public :: sm_mp=2
  integer, parameter, public :: sm_cmu=3
! units
  integer,parameter, public :: units_ev=1
  integer,parameter, public :: units_au=2
  integer,parameter, public :: units_si=3
! mathematics
  real(pr), parameter, public :: pi = 3.14159265358979323846264338327950_pr
  
  
! physics


end module constants
