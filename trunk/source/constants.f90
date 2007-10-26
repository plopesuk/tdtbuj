!> \brief defines constants that are needed in the code.                  
!> \author Alin M Elena
!> \date January 2007


module constants
  implicit none
  private
  

  integer, parameter, public :: pr=kind(1.0d0) !< define precision for reals
  integer, parameter, public :: ml=255 !< number of character per line
  integer, parameter, public :: mw=40 !< number of character per word 

  integer, parameter,public :: low_verbos=5
  integer, parameter,public :: medium_verbos=15
  integer, parameter,public :: high_verbos=25
!> verbosity levels

  integer, parameter,public :: low_debug=5
  integer, parameter,public :: medium_debug=15
  integer, parameter,public :: high_debug=25
!> debug levels
!!!!!!!!!!


 integer, parameter, public :: run_sp = 1
 integer, parameter, public :: run_bo = 2
 integer, parameter, public :: run_ehrenfest = 3
 integer, parameter, public :: run_fit = 4
 integer, parameter, public :: run_force_test = 5
 integer, parameter, public :: run_force_testx = 6
 integer, parameter, public :: run_force_testy = 7
 integer, parameter, public :: run_force_testz = 8
 integer, parameter, public :: run_ehrenfest_damped = 9 
!> RunType

 integer, parameter, public :: scf_tbuj = 1
!> SCFType

integer, parameter, public :: scfmix_broyden = 1
integer, parameter, public :: scfmix_pulay = 2
!>SCFMixType

integer, parameter, public :: electrostatics_point=1
integer, parameter, public :: electrostatics_multipoles=2
!> electrostatics

  integer, parameter, public :: spin_down=0
  integer, parameter, public :: spin_up=1
!>spins 
  integer, parameter, public :: bond_gsp=1
  integer, parameter, public :: bond_harrison=2
!> bond type

  integer, parameter, public :: sm_fd=1
  integer, parameter, public :: sm_mp=2
  integer, parameter, public :: sm_cmu=3
!> smearing methods

   integer,parameter, public :: units_ev=1
   integer,parameter, public :: units_au=2
   integer,parameter, public :: units_si=3
   integer,parameter, public :: units_ry=4
!> units
   
 
 ! mathematics
  real(pr), parameter, public ::   pi = 3.14159265358979323846264338327950_pr !< \f$\pi\f$
  real(pr),parameter, public :: infinity = huge(1.0_pr) !< \f$ \infinity \f$  
  
! physics
  integer,parameter,public  :: max_orbitals_per_atom = 25
   real(pr),  public :: amu_to_internal !< atomic mass units to internal mass units
   real(pr),  public :: boltzmann_k  !< Boltzmann constant \f$ k_B \f$
   real(pr), parameter,  public ::  ev_to_hartree = 0.03674932600434263_pr !< electron volt to Hartree conversion factor
   real(pr), parameter, public ::   hatree_to_ev = 27.211383411_pr !< hartree to electron volt conversion factor
   real(pr),  public :: hbar !< Plank's constant \f$\hbar\f$  
   real(pr), parameter, public ::   bohr_to_ang = 0.5291772083_pr !< bohr to \f$\AA\f$ conversion factor
   real(pr), parameter, public ::    efs_to_microamps = 160.2176462_pr
   real(pr),  public :: epsilon0 !< \f$ \epsilon_0 \f$ vacuum permiability
   real(pr),  public :: e2 !< square of the elctronic charge
   real(pr),  public :: me !< mass of the electron
   real(pr), public :: charge2SI !< charge to SI conversion factor
   real(pr), public :: length2SI !< length to SI conversion factor
   real(pr), parameter,public :: debye2SI=3.335640952D-30 !< debye to SI (Cm) conversion factor


   public :: initialize_constants

contains 

!> \brief Initializes the constants according to the system of units chosen
!> \author Alin M Elena
!> \date 2007
!> \param units selects the system of units

   subroutine initialize_constants(units)
      integer,intent(in) :: units     
      select case(units)
      case(units_ev)
         amu_to_internal = 103.6426867_pr
         boltzmann_k = 8.61734215D-5 ! in eV / K
         hbar = 0.65821188926_pr ! in eV fs
         epsilon0 = 5.526349954D-3
         e2=1.0_pr
         me=1.0_pr ! 
         charge2SI= 1.60217653D-19 !Cm
         length2SI=1.0D-10 !m
      case(units_au)
         amu_to_internal =  1836.15267247_pr
         boltzmann_k = 3.1668154D-6 !Eh/K
         hbar = 1.0_pr 
         epsilon0 = 1/(4.0_pr*pi)
         e2=1.0_pr
         me=1.0_pr
         charge2SI= 1.60217653D-19 !Cm
         length2SI=0.5291772108D-10 !m
      case(units_si)
         amu_to_internal = 1.66053886D-27 !kg
         boltzmann_k = 1.3806505D-23 !J/K 
         hbar = 1.05457168D-34 !J s 
         epsilon0 = 8.854187817D-12 !F/m
         e2=1.60217653D-19**2 !C^2
         me= 9.1093826D-31! kg
         charge2SI= 1.0_pr !Cm
         length2SI=1.0_pr !m
        case(units_ry)
        !!!!!!!!!! to be checked
        
         amu_to_internal =  1836.15267247_pr
         boltzmann_k = 6.3361921D-6 !Rydberg/K
         hbar = 1.0_pr 
         epsilon0 = 1.0_pr/(4.0_pr*pi)
         e2=2.0_pr
         me=0.5_pr
         ! it should be devided by sqrt(2.0) but because we do want to see the charges in meaningful 
         ! numbers we do not do it                        
         charge2SI= 1.60217653D-19!/sqrt(2.0_pr) 
         length2SI=0.5291772108D-10 !m

         
      end select         

   end subroutine initialize_constants




end module constants
