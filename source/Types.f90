!> \brief defines different data types used around the program
!> \author Alin M Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \remarks

module m_Types
  use m_Constants
  implicit none
  private
!> \brief I/O data type
!> \details see m_ReadData::read_io
  type, public :: ioType
  character(len=k_mw) :: inpFile !<  name of the input file from which the data is read
!>  name of the error input file (parsing errors, before opening specific output files)
  character(len=k_mw) :: inpErr
!>  name of the file in which the output is put
  character(len=k_mw)  :: outputFile
!>  how much information to be Printed in the output file
  integer :: verbosity
!>   name of the file in which the debug data is written
  character(len=k_mw)  :: debugFile
!>   name of the file in which the animation data is written
  character(len=k_mw)  :: aniFile
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
!>  unit number for animation file
  integer :: uani
!> is it first time when is read?
  logical :: firstTime = .false.
  end type ioType

!> \brief time keeping type
  type ,public :: timeType
  real(k_pr) :: start !< start time
  real(k_pr) :: end !< end time
  real(k_pr) :: int !< interk_mediate time
  end type timeType

!> \brief fitting data type
!> \internal point to more explanations
  type, public :: fitType
    integer :: neps,feval,ns
    integer ::iter,nt
    real(k_pr) :: fitTol
    logical :: restartFit
    integer :: fitMethod
    real(k_pr):: temp
    real(k_pr) :: step
    real(k_pr) ::stepAd
    real(k_pr) :: rt
    integer :: iNoParams !< no of parameters in optimization
  end type fitType

!> \brief see m_ReadData::ReadGeneral module for the meaning
  type,public :: generalType
    real(k_pr) :: bias
    integer :: BiasRampSteps
    real(k_pr) :: electronicTemperature
    real(k_pr) :: electronicMu
    real(k_pr) :: ionicTemperature
    real(k_pr) :: netcharge
    real(k_pr) :: deltat
    real(k_pr) :: forceTolerance
    integer  :: nsteps
    integer  :: runType
    logical  :: scf
    integer  :: scfType
    integer  :: maxscf
    real(k_pr) :: scfMix
    real(k_pr) :: scftol
    integer  :: scfMixn
    logical  :: scaleVelocities !< if true scales the velocities
    real(k_pr) :: dmOccupationTolerance
    integer  :: maxOrbitalsPerAtom
    real(k_pr) :: hElementThreshold
!   real(k_pr) :: hessian_displacement
    logical  :: spin
    logical  :: collinear
!   real(k_pr) :: sdu
    integer  :: EulerSteps
    integer  :: electrostatics
    logical  :: compElec !<precompute the q's and v's for tb+u multipole case if is set to false
    integer  :: units
    integer  :: bond
    logical  :: embedding !< Has the model embedding?
! !   for the force test
    integer :: fsteps   !< no of steps to sample the force testing step
    integer :: fatom  !< atom on which we test
    real(k_pr) :: fstart !< the starting position of the atom
    real(k_pr) :: fdx !< step for numerical derivative used in testing
! ! for smearing methods Methfessel&Paxton
    integer :: mpN
    real(k_pr) :: mpW
    integer :: smearMethod
    real(k_pr) :: gamma
! ! create excitation
    logical :: lIsExcited
    integer :: holeState
    integer :: exciteState
    integer :: holeSpin
    integer :: exciteSpin
    character(len=k_mw) :: jobName !< name of the job
    real(k_pr) :: ranseed !< seed to Initialize the random number generator
    logical :: writeAnimation !< if true will write the coordinates in a file
    integer :: AnimationSteps !< writes only steps which are multiple of it
    logical :: ReadVelocity !< read velocity block?
    logical :: firstTime = .false. !< is it first time when is read
    logical :: screened
    logical :: SymRefRho
    logical :: lIsSCFConverged
    integer :: maxIt !< no of iterations to use in bisection method
    real(k_pr) :: qTolerance !< charge tolerance used to find Fermi level
    type(timeType):: time
    type(fitType)  :: fit
    real(k_pr) :: CurrSimTime=0.0_k_pr !< keeps the time in a td simulation
    real(k_pr) :: epsG,gtol
    real(k_pr) :: epsX,xtol
    real(k_pr) :: epsF,ftol
    integer :: maxFEval
    integer :: HessianM
    integer :: geomAlg !< the geometry optimisation algorithm
    integer :: wdensity !< what kind of density shloud be used for Eherenfest Damped
    logical :: hasElectricField !< do we apply an external electric field?
    real(k_pr) :: E(1:3)  !< the components of the electric field
    integer :: etd !< selects the type of function that gives time dependent
                   !< part for the electric field
    real(k_pr) ::freq !< frequency of the electric field
    real(k_pr) :: phi0 !< initial phase of the electric field
    real(k_pr) :: t0 !< time t0 at which we apply the electric field
    real(k_pr) :: sigma !< the width of the gaussian
  end type generalType

!> data type that keeps the neighbours list
type, private :: neighbourType
  integer :: n=0 !< no of neighbours
  integer, allocatable :: a(:) !< the list of neighbours
  logical :: created=.false. !< was the space allocated?
end type neighbourType

!> \brief data type for atoms properties
  type, public :: atomicType
    logical               :: created !< was it read?
    integer               :: natoms !< no of atoms
    integer,  allocatable :: id(:) !< id of the atoms
    integer,  allocatable :: sp(:) !< specie of the atoms
    logical,  allocatable :: ismoving(:) !< is it moving?
    logical,  allocatable :: isscf(:) !< is it SCF?
    integer, allocatable :: scf(:) !< list of SCCF atoms
    integer, allocatable :: moving(:) !< list of moving atoms
    real(k_pr), allocatable :: x(:),y(:),z(:) !< cartesian coordinates
    real(k_pr), allocatable :: dx(:),dy(:),dz(:) !< cartesian dipole moments
    real(k_pr), allocatable :: chrg(:),chrg0(:) !< excess charge, initial excess charge
    real(k_pr), allocatable :: fx(:), fy(:), fz(:) !< cartesian forces
    real(k_pr), allocatable :: vx(:), vy(:), vz(:) !< cartesian velocities
    real(k_pr), allocatable :: xo(:),yo(:),zo(:) !< old cartesian coordinates
    real(k_pr), allocatable :: fxo(:), fyo(:), fzo(:) !< old cartesian forces
    real(k_pr), allocatable :: vxo(:), vyo(:), vzo(:) !< old cartesian velocities
    real(k_pr), allocatable :: bias(:) !< bias of the atoms
    integer(k_pr), allocatable :: norbs(:) !< no of orbitals of atom
    integer,  allocatable  :: orbs(:,:) !< orbitals, first index atom, second orbital
    type(neighbourType),allocatable :: neighbours(:) !< keeps a list ofneighbours
    real(k_pr) :: tdipx,tdipy,tdipz !< total dipole moment, cartesian components
    integer :: nacceptor !< no of atoms in acceptor
    integer :: ndonor !< no of atoms in donor
    integer :: nspacer !< no of atoms in spacer
    integer, allocatable :: acceptor(:) !< atoms in acceptor
    integer, allocatable :: donor(:) !< atoms in donor
    integer, allocatable :: spacer(:) !< atoms in spacer
    integer :: nscf !< no of scf atoms
    integer :: nmoving !< no of moving atoms
    real(k_pr),allocatable :: MagMom(:) !< local magnetic momment
    real(k_pr) :: MagneticMoment !< total magnetic momment
    integer :: ncurrent !< no of atoms on which to compute the current
    integer, allocatable :: current(:) !< the atoms on which to compute the current
    integer :: nCurrentOnBonds=0 !< number of the currents on bonds that will be computed
    integer, allocatable :: currentOnBonds(:) !< keeps the orbitals on which we compute the current
  end type atomicType

!> \brief species type
  type, public :: speciesType
    logical :: created = .false. !< is it created?
    integer :: nspecies !< no of species
    integer, allocatable :: id(:) !< if of specie
    real(k_pr), allocatable :: mass(:) !< atomic mass of specie
    integer,  allocatable :: z(:) !< Z of specie
    real(k_pr), allocatable :: zval(:) !< no of valence electrons for each specie
    real(k_pr), allocatable :: ulocal(:,:) !< local U (Hubbard U) (first index deals with the species counter second with the orbital one)
    real(k_pr), allocatable :: jlocal(:,:) !< local J (first index deals with the species counter second with the orbital one)
    real(k_pr), allocatable :: uinter(:) !< Uinter (used for the screaning of k_electrostatics)
    integer,  allocatable :: norbs(:) !< no of orbitals
  end type speciesType

!> \brief atomic orbital data type
  type, public :: orbitalType
    integer :: atom,sp,n,l,m !< atom to which belongs and specie, quantum numbers n,l,m
    logical :: spin !< spin up or down up=.true. down=.false.
    real(k_pr) :: occup !< occupancy of the orbital
  end type orbitalType

!> \brief atomic basis data type
  type, private :: basisType
    integer :: norbitals !< no of orbitals
    type(orbitalType), allocatable :: orbitals(:) !< orbitals
  end type basisType

!> \brief to be added
  type,private :: deltaType
    integer :: sp,l
    real(k_pr),allocatable :: d(:,:,:)
  end type deltaType

!> \brief all atomic data type
  type, public :: atomicxType
    type(atomicType) :: atoms !< atomic properties
    type(speciesType) :: species !< specie properties
    type(basisType) :: basis !< system basis
    type(orbitalType), allocatable :: speciesBasis(:,:) !< species basis set
  end type atomicxType

  type, public :: qvs
    integer :: dim
    real(k_pr), allocatable :: a(:)
  end type qvs

!> \brief tail parameters
  type, public :: tailType
    real(k_pr) :: rIn
    real(k_pr) :: a, b, c, rOut
  end type tailType

!> \brief hopping parameters
  type, public :: hopMatrixType
    integer :: l1,l2,ll ! nn -->l1, mm-->l2,  ll-->m
    real(k_pr) :: a1,a2,a3,a4 !< embedding parameters
    real(k_pr) :: phi0,r0,rc,r1,rcut,n,nc,d0,dc,d1,dcut,m,mc
    real(k_pr), allocatable  :: eps(:)
    real(k_pr), allocatable :: a(:,:,:)
    type(tailType)  :: repTail
    type(tailType),allocatable:: hmnTail(:,:,:)
  end type hopMatrixType

!> \brief data type for tight-binding model
  type, public :: modelType
    type(deltaType), allocatable :: delta(:) !< delta params for each specie
    type(hopMatrixType), allocatable :: hopping(:,:) !< hopping parameters
  end type modelType

!> square matrix type it may keep both dense and sparse matrices
  type, public :: matrixType
    logical :: isSparse !< is it sparse?
    logical :: created !< is it created?
    integer :: dim !< dik_mesion of the matrix
    integer :: nonZero !< no of non zero elements in a sparse matrix
    integer, allocatable :: indx(:), jndx(:) !< the indeces of the non zero elements for sparse
    complex(kind=k_pr), allocatable :: a(:,:) !< keeps the data
  end type matrixType
 type :: buffersType
    real(k_pr), allocatable :: densityin(:)
    real(k_pr), allocatable :: densityout(:)
    real(k_pr), allocatable :: densitynext(:)
    real(k_pr), allocatable:: dins(:,:)
    real(k_pr), allocatable:: douts(:,:)
    real(k_pr), allocatable::  res(:,:)
    real(k_pr),allocatable :: tmpA(:)
    type(matrixType) :: tmpB, h
    real(k_pr),allocatable  :: f(:) ,g(:)
    complex(k_pr),allocatable :: a(:,:)
    integer, allocatable :: pos1(:), pos2(:),itmp(:)
    real(k_pr), allocatable :: nstart(:)
 end type buffersType

!> defines the type for the the data needed to solve the problem
  type,public :: solutionType
    type(matrixType) :: h !< total hamiltonian
    type(matrixType) :: hin !< hin hamiltonian without scf "condiments" added
    type(matrixType) :: h2 !< h2 hamiltonian the scf "condiments"
    type(matrixType) :: hdown !< spin down hamiltonian
    type(matrixType) :: hup !< spin up hamiltonian
    type(matrixType) :: eigenvecs !< eigenvectors of h
    type(matrixType) :: forceOp !< force operator matrix
    real(k_pr), allocatable :: eigenvals(:) !< eigenvalues of h
    real(k_pr), allocatable :: n0(:) !< n0 initial density matrix in vector representation
    real(k_pr), allocatable :: potential(:) !< electrostatic potential on each atom
    real(k_pr), allocatable :: density(:) !< keeps the \f$ \delta q \f$
    real(k_pr), allocatable :: field(:,:) !< electrostatic field
    real(k_pr), allocatable :: gcoeff(:,:,:) !< precomputed Gaunt Coefficients
    real(k_pr), allocatable ::  rgc(:,:,:) !< precomputed real Gaunt coefficients
    real(k_pr), pointer,dimension(:) :: fact(:)=> NULL() !< precomputed factorial values
    real(k_pr) :: seed(97) !< seeds for the random number generator \see rmarin
    type(matrixType) :: rho !< density matrix
    real(k_pr) :: electronicEntropy !< electronic entropy of the system
    real(k_pr) :: totalEnergy !< total energy of the system
    type(buffersType) :: buff
    real(k_pr),pointer, dimension(:) :: hermite => NULL() !< hermite polynomials used for Methfessel&Paxton smearing method
    real(k_pr),allocatable :: Distances(:,:) !< euclidean distances matrix
    real(k_pr), allocatable :: CurrentMatrix(:,:) !< bond currents matrix by orbitals
    real(k_pr), allocatable :: CurrentMatrix2(:,:) !< bond currents matrix by atoms, the diagonal term contains the total bond current, only specified atoms entries will be computed and populated
    type(MatrixType) :: rhodot,deltaRho,rhonew,rhoold,rho0 !< used for Ehrenfest dynamic
    type(qvs), allocatable :: delq(:), vs(:)
  end type solutionType
end module m_Types
