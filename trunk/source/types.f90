!> \brief defines different data types used around the program
!> \details 
!> \author Alin M Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \remarks



module types
  use constants
  implicit none
  private
!> I/O data type
  type, public :: io_type

  character(len=mw) :: inp_file !<  name of the input file from which the data are read
!>  name of the error input file (parsing errors, before opening specific output files)  
character(len=mw) :: inp_err
!>  name of the file in which the output is put 
  character(len=mw)  :: output_file
!>  how much information to be printed in the output file 
  integer :: verbosity
!>   name of the file in which the debug data are written
  character(len=mw)  :: debug_file
!>  the debug level
  integer :: debug
!>  output on screen?
  logical :: stdout
!>  unit number for input file
  integer :: uinp
!>  unit number for output file  
  integer :: uout
!>  unit number for debug file
  integer :: udeb=-1
!>  unit number for error file
  integer :: uerr
!> is it first time when is read?
  logical :: first_time = .false.
  end type io_type

  type ,private :: time_type
  real(pr) :: start
  real(pr) :: end
  real(pr) :: int
  end type time_type 

  type,public :: general_type
  real(pr) :: electronic_temperature
  real(pr) :: electronic_mu
  real(pr) :: ionic_temperature
  real(pr) :: netcharge
  real(pr) :: deltat
  integer  :: nsteps	
  integer  :: runtype
  logical  :: scf
  integer  :: scf_type
  integer  :: maxscf
  real(pr) :: scfmix
  integer  :: scfmixtype
  real(pr) :: scftol
  integer  :: scfmixn
  logical  :: velocity_scaling
  real(pr) :: dm_occupation_tolerance
  integer  :: max_orbitals_per_atom 
  real(pr) :: h_element_threshold
  real(pr) :: hessian_displacement
  logical  :: spin 
  logical  :: collinear
  real(pr) :: sdu
  integer  :: Euler_steps
  integer  :: electrostatics
!   precompute the q's and v's for tb+u multipole case if is set to false	
  logical  :: comp_elec 
  integer  :: units
  integer  :: bond
  logical  :: embedding
!   for the force test
  integer :: f_steps
  real(pr) :: f_start
  real(pr) :: f_dx
! for smearing methods Methfessel&Paxton
  integer :: mp_N
  integer :: smethod
  real(pr) :: gamma
! create excitation
  integer :: hole
  integer :: excite
  integer :: hole_spin
  integer :: excite_spin
  character(len=20) :: job_name
  integer :: ran3_seed
  logical :: write_ani
  logical :: write_ene
  logical :: write_density
  logical :: read_density
  logical :: read_velocity
  logical :: first_time = .false.
!!!!!!
  type(time_type):: time
  end type general_type

! 	type, public :: fit_p
! 	integer :: neps,feval,ns,ipr
! 	integer ::iter,nt
! 	real(pr) :: fit_ftol
! 	logical :: restart_fit
! 	character(15) :: fit_type
! 	real(pr):: temp
! 	real(pr) :: step
! 	real(pr) ::step_ad
! 	real(pr) :: rt
! 	end type fit_p


  type, public :: atomic_type              
  logical               :: created 
  integer               :: natoms
  integer               :: nmoving
  integer,  allocatable :: id(:)
  integer,  allocatable :: sp(:)
  integer,  allocatable :: moving(:)
  integer,  allocatable :: isscf(:)
  real(pr), allocatable :: x(:),y(:),z(:)
  real(pr), allocatable :: dx(:),dy(:),dz(:) 
  real(pr), allocatable :: chrg(:),chrg0(:)
  real(pr), allocatable :: fx(:), fy(:), fz(:)
  real(pr), allocatable :: vx(:), vy(:), vz(:)
  real(pr), allocatable :: xo(:),yo(:),zo(:)
  real(pr), allocatable :: fxo(:), fyo(:), fzo(:)
  real(pr), allocatable :: vxo(:), vyo(:), vzo(:)
  real(pr), allocatable :: bias(:)
  integer,  allocatable  :: orbs(:,:)
  real(pr) :: tdipx,tdipy,tdipz
  integer :: nacceptor
  integer :: ndonor
  integer :: nspacer
  integer, allocatable :: acceptor(:)
  integer, allocatable :: donor(:)
  integer, allocatable :: spacer(:)

  end type atomic_type
! 
  type, public :: species_type
  logical :: created
  integer :: nspecies
  integer, allocatable :: id(:)
  real(pr), allocatable :: mass(:)
  integer,  allocatable :: z(:)
  real(pr), allocatable ::zval(:)
  logical,  allocatable :: isqm(:)
  real(pr), allocatable :: ulocal(:)
  real(pr), allocatable :: jlocal(:)
  real(pr), allocatable :: uinter(:)
  integer,  allocatable :: norbs(:)	
  end type species_type


  type, public :: orbital_type
  integer :: atom,sp,n,l,m
  logical :: spin
  real(pr) :: occup
  end type orbital_type

  type, public :: basis_type
  integer :: norbital
  integer  :: max_orbitals_per_atom 
  type(orbital_type), allocatable :: orbitals(:)
  end type basis_type

  type,public :: delta_type
  integer :: sp,l
  real(pr),allocatable :: d(:,:,:)
  end type delta_type


! variables now



end module types
