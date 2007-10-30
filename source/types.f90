!> \brief defines different data types used around the program
!> \author Alin M Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \remarks

module types
  use constants
  implicit none
  private
!> \brief I/O data type
!> \details see read_data::read_io
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

!> \brief time keeping type 
  type ,public :: time_type
  real(pr) :: start !< start time
  real(pr) :: end !< end time
  real(pr) :: int !< intermediate time
  end type time_type 

!> \brief see read_data::read_general module for the meaning
  type,public :: general_type
    real(pr) :: bias
!   real(pr) :: electronic_temperature
!   real(pr) :: electronic_mu
!   real(pr) :: ionic_temperature
!   real(pr) :: netcharge
!   real(pr) :: deltat
!   integer  :: nsteps	
!   integer  :: runtype
   logical  :: scf
!   integer  :: scf_type
!   integer  :: maxscf
!   real(pr) :: scfmix
!   integer  :: scfmixtype
!   real(pr) :: scftol
!   integer  :: scfmixn
!   logical  :: velocity_scaling
!   real(pr) :: dm_occupation_tolerance
   integer  :: max_orbitals_per_atom
!   real(pr) :: h_element_threshold
!   real(pr) :: hessian_displacement
   logical  :: spin
!   logical  :: collinear
!   real(pr) :: sdu
!   integer  :: Euler_steps
!   integer  :: electrostatics
! !   precompute the q's and v's for tb+u multipole case if is set to false	
!   logical  :: comp_elec 
!   integer  :: units
!   integer  :: bond
!   logical  :: embedding
! !   for the force test
!   integer :: f_steps
!   real(pr) :: f_start
!   real(pr) :: f_dx
! ! for smearing methods Methfessel&Paxton
!   integer :: mp_N
!   integer :: smethod
!   real(pr) :: gamma
! ! create excitation
!   integer :: hole
!   integer :: excite
!   integer :: hole_spin
!   integer :: excite_spin
  character(len=mw) :: job_name !< name of the job
  integer :: ran3_seed !< seed to initialize the random number generator
!   logical :: write_ani
!   logical :: write_ene
!   logical :: write_density
!   logical :: read_density
   logical :: read_velocity !< read velocity block?
   logical :: first_time = .false. !< is it first time when is read 
!   logical :: scc
!   logical :: SymRefRho
!   logical :: bIsSCFConverged
!!!!!!
  type(time_type):: time
  
  end type general_type
  
!> \brief fitting data type
!> \todo complete documentation
  type, public :: fit_p
    integer :: neps,feval,ns,ipr
    integer ::iter,nt
    real(pr) :: fit_ftol
    logical :: restart_fit
    character(15) :: fit_type
    real(pr):: temp
    real(pr) :: step
    real(pr) ::step_ad
    real(pr) :: rt
    integer :: iNoParams !< no of parameters in optimization
  end type fit_p

!> \brief data type for atoms properties
  type, public :: atomic_type
    logical               :: created !< was it read?
    integer               :: natoms !< no of atoms
    integer,  allocatable :: id(:) !< id of the atoms
    integer,  allocatable :: sp(:) !< specie of the atoms 
    logical,  allocatable :: ismoving(:) !< is it moving?
    logical,  allocatable :: isscf(:) !< is it SCF?
    integer, allocatable :: scf(:) !< list of SCCF atoms
    integer, allocatable :: moving(:) !< list of moving atoms
    real(pr), allocatable :: x(:),y(:),z(:) !< cartesian coordinates
    real(pr), allocatable :: dx(:),dy(:),dz(:) !< cartesian dipole moments
    real(pr), allocatable :: chrg(:),chrg0(:) !< excess charge, initial excess charge 
    real(pr), allocatable :: fx(:), fy(:), fz(:) !< cartesian forces
    real(pr), allocatable :: vx(:), vy(:), vz(:) !< cartesian velocities
    real(pr), allocatable :: xo(:),yo(:),zo(:) !< old cartesian coordinates
    real(pr), allocatable :: fxo(:), fyo(:), fzo(:) !< old cartesian forces
    real(pr), allocatable :: vxo(:), vyo(:), vzo(:) !< old cartesian velocities
    real(pr), allocatable :: bias(:) !< bias of the atoms
integer(pr), allocatable :: norbs(:) !< no of orbitals of atom
    integer,  allocatable  :: orbs(:,:) !< orbitals, first index atom, second orbital
    real(pr) :: tdipx,tdipy,tdipz !< total dipole moment, cartesian components
    integer :: nacceptor !< no of atoms in acceptor
    integer :: ndonor !< no of atoms in donor
    integer :: nspacer !< no of atoms in spacer
    integer, allocatable :: acceptor(:) !< atoms in acceptor
    integer, allocatable :: donor(:) !< atoms in donor
    integer, allocatable :: spacer(:) !< atoms in spacer
    integer :: nscf !< no of scf atoms
    integer :: nmoving !< no of moving atoms
  end type atomic_type


!> \brief species type 
  type, public :: species_type
    logical :: created = .false. !< is it created?
    integer :: nspecies !< no of species
    integer, allocatable :: id(:) !< if of specie
    real(pr), allocatable :: mass(:) !< atomic mass of specie
    integer,  allocatable :: z(:) !< Z of specie
    real(pr), allocatable :: zval(:) !< no of valence electrons for each specie
    real(pr), allocatable :: ulocal(:) !< local U (Hubbard U) 
    real(pr), allocatable :: jlocal(:) !< local J
    real(pr), allocatable :: uinter(:) !< Uinter (used for the screaning of electrostatics)
    integer,  allocatable :: norbs(:) !< no of orbitals	
  end type species_type

!> \brief atomic orbital data type
  type, private :: orbital_type
    integer :: atom,sp,n,l,m !< atom to which belongs and specie, quantum numbers n,l,m
    logical :: spin !< spin up or down up=.true. down=.false.
    real(pr) :: occup !< occupancy of the orbital
  end type orbital_type

!> \brief atomic basis data type
  type, private :: basis_type
    integer :: norbitals !< no of orbitals  
    type(orbital_type), allocatable :: orbitals(:) !< orbitals
  end type basis_type
!> \brief to be added
  type,private :: delta_type
    integer :: sp,l
    real(pr),allocatable :: d(:,:,:)
  end type delta_type

!> \brief all atomic data type 
  type, public :: atomicx_type
    type(atomic_type) :: atoms !< atmoic properties
    type(species_type) :: species !< specie properties
    type(basis_type) :: basis !< system basis
    type(delta_type),allocatable :: delta(:) !< delta params for easch specie
    type(orbital_type), allocatable :: species_basis(:,:) !< species basis set
  end type atomicx_type


end module types
