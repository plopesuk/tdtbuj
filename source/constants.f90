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
!> SCFMixType

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
  real(pr), parameter, public ::   pi = 3.14159265358979323846264338327950_pr !< \f$ \pi \f$
  real(pr),parameter, public :: infinity = huge(1.0_pr) !< infinity
  
! physics
  integer,parameter,public  :: max_orbitals_per_atom = 25
   real(pr),  public :: amu_to_internal !< atomic mass units to internal mass units
   real(pr),  public :: boltzmann_k  !< Boltzmann constant \f$ k_B \f$
   real(pr), parameter,  public ::  ev_to_hartree = 0.03674932600434263_pr !< electron volt to Hartree conversion factor
   real(pr), parameter, public ::   hatree_to_ev = 27.211383411_pr !< hartree to electron volt conversion factor
   real(pr),  public :: hbar !< Plank's constant \f$ \hbar \f$  
   real(pr), parameter, public ::   bohr_to_ang = 0.5291772083_pr !< bohr to  Angstroms conversion factor
   real(pr), parameter, public ::    efs_to_microamps = 160.2176462_pr
   real(pr),  public :: epsilon0 !< \f$ \epsilon_0 \f$ vacuum permiability
   real(pr),  public :: e2 !< square of the elctronic charge
   real(pr),  public :: me !< mass of the electron
   real(pr), public :: charge2SI !< charge to SI conversion factor
   real(pr), public :: length2SI !< length to SI conversion factor
  real(pr), parameter,public :: debye2SI=3.335640952D-30 !< debye to SI (Cm) conversion factor
  integer, parameter :: nz=110 !< no of atomic elements that have info associated with


  public :: initialize_constants
  public :: symbol
  public :: el_name
  public :: el_config
  public :: covalent_radius
  public :: period
  public :: vdw_radius
  public :: weight
  public :: group


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

!> \brief returns the symbol for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:19:04
!> \param iz integer Z of the atom
  character(len=2) function Symbol(iz)  
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'symbol'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=2) :: dummy(1:nz)

    data dummy /"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",&
      "P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",&
      "Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru",&
      "Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr",&
      "Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W",&
      "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",&
      "Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",&
      "Rf","Db","Sg","Bh","Hs","Mt","Ds"/

    if ((iz<1).or.(iz>nz)) then
      write(6,*) ' symbol: out of range iz =',iz
      symbol = ' '
    else
      symbol = dummy(iz)
    endif
  end function symbol

!> \brief returns the name for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:18:26
!> \param iz integer Z of the atom
  character(len=mw) function el_name(iz)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'el_name'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=mw) :: dummy(1:nz)

    data dummy /"Hydrogen","Helium","Lithium","Beryllium","Boron","Carbon",&
      "Nitrogen","Oxygen","Fluorine","Neon","Sodium","Magnesium","Aluminium",&
      "Silicon","Phosphorus","Sulphur","Chlorine","Argon","Potassium","Calcium",&
      "Scandium","Titanium","Vanadium","Chromium","Manganese","Iron","Cobalt",&
      "Nickel","Copper","Zinc","Gallium","Germanium","Arsenic","Selenium",&
      "Bromine","Krypton","Rubidium","Strontium","Yttrium","Zirconium","Niobium",&
      "Molybdenum","Technetium","Ruthenium","Rhodium","Palladium","Silver",&
      "Cadmium","Indium","Tin","Antimony","Tellurium","Iodine","Xenon","Caesium",&
      "Barium","Lanthanum","Cerium","Praseodymium","Neodymium","Promethium",&
      "Samarium","Europium","Gadolinium","Terbium","Dysprosium","Holmium",&
      "Erbium","Thulium","Ytterbium","Lutetium","Hafnium","Tantalum","Tungsten",&
      "Rhenium","Osmium","Iridium","Platinum","Gold","Mercury","Thallium",&
      "Lead","Bismuth","Polonium","Astatine","Radon","Francium","Radium",&
      "Actinium","Thorium","Protactinium","Uranium","Neptunium","Plutonium",&
      "Americium","Curium","Berkelium","Californium","Einsteinium","Fermium",&
      "Mendelevium","Nobelium","Lawrencium","Rutherfordium","Dubnium",&
      "Seaborgium","Bohrium","Hassium","Meitnerium","Darmstadtium"/

    if ((iz<1).or.(iz>nz)) then
      write(6,*) ' element name: out of range iz =',iz
      el_name = ' '
    else
      el_name = dummy(iz)
    endif
  end function el_name

!> \brief returns the period from periodic for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:17:39
!> \param iz integer Z of the atom

  character(len=4) function period(iz)  
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'period'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=4) :: dummy(1:nz)

    data dummy /"1","1","2","2","2","2","2","2","2","2",&
      "3","3","3","3","3","3","3","3",&
      "4","4","4","4","4","4","4","4","4","4","4","4","4","4","4","4","4","4",&
      "5","5","5","5","5","5","5","5","5","5","5","5","5","5","5"," 5"," 5"," 5",&
      "6","6","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan",&
      "6Lan","6Lan","6Lan","6Lan","6","6","6","6","6","6","6","6","6","6",&
      "6","6","6","6","6","6","7","7","7Act","7Act","7Act","7Act","7Act","7Act",&
      "7Act","7Act","7Act","7Act","7Act","7Act","7Act","7Act","7","7","7","7",&
      "7","7","7","7"/

    if ((iz<1).or.(iz>nz)) then
      write(6,*) 'period: out of range iz =',iz
      period= " "
    else
      period= dummy(iz)
    endif
  end function period

!> \brief returns the group from periodic table for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:17:12
!> \param iz integer Z of the atom

  integer function group(iz)  
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'group'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    integer :: dummy(1:nz)

    data dummy /1,18,1,2,13,14,15,16,17,18,1,2,13,14,15,16,17,18,1,2,3,4,5,6,7,8,&
      9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,&
      1,2,19,19,19,19,19,19,19,19,19,19,19,19,19,19,3,4,5,6,7,8,9,10,11,12,13,&
      14,15,16,17,18,1,2,20,20,20,20,20,20,20,20,20,20,20,20,20,20,3,4,5,6,7,8,&
      9,10/

    if ((iz<1).or.(iz>nz)) then
      write(6,*) 'group: out of range iz =',iz
      group= 0
    else
      group= dummy(iz)
    endif
  end function group

!> \brief returns electronic configuration of valence band for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:15:58
!> \param iz integer Z of the atom
!> \todo complete the electronic configuration
  character(len=mw) function el_config(iz)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'el_config'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=mw) :: dummy(1:nz)

      data dummy /"1s1","1s2","2s1","2s2","2s2 2p1","2s2 2p2","2s2 2p3","2s2 2p4","2s2 2p5",&
      "2s2 2p6","3s1","3s2","3s2 3p1","3s2 3p2","3s2 3p3",&
      "p","p","p","s","s","d","d","d","d","d","d","d","d","d","d","p","p",&
      "p","p","p","p","s","s","d","d","d","d","d","d","d","d","d","d","p",&
      "p","p","p","p","p","s","s","f","f","f","f","f","f","f","f","f","f",&
      "f","f","f","f","d","d","d","d","d","d","d","d","d","d","p","p","p",&
      "p","p","p","s","s","f","f","f","f","f","f","f","f","f","f","f","f",&
      "f","f","d","d","d","d","d","d","d","d"/

    if ((iz<1).or.(iz>nz)) then
      write(6,*) 'electronic configuration: out of range iz =',iz
      el_config = " "
    else
      el_config = dummy(iz)
    endif
  end function el_config

!> \brief returns atomic weight for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:15:38
!> \param iz integer Z of the atom

  real(pr) function weight(iz)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'weight'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    real(pr) :: dummy(1:nz)

    data dummy /1.008,4.003,6.941,9.012,10.811,12.011,14.007,15.999,18.998,20.18,&
      22.991,24.305,26.982,28.086,30.974,32.066,35.453,39.948,39.098,40.078,&
      44.956,47.867,50.942,51.996,54.938,55.845,58.933,58.693,63.546,65.39,&
      69.723,72.61,74.922,78.96,79.904,83.8,85.468,87.62,88.906,91.224,92.906,&
      95.94,98,101.07,102.906,106.42,107.868,112.411,114.818,118.71,121.76,&
      127.6,126.904,131.29,132.905,137.327,138.906,140.116,140.908,144.24,&
      145,150.36,151.964,157.25,158.925,162.5,164.93,167.26,168.934,173.04,&
      174.967,178.49,180.948,183.84,186.207,190.23,192.217,195.078,196.967,&
      200.59,204.383,207.2,208.98,210,210,222,223,226,227,232.038,231.036,&
      238.029,237,244,243,247,247,251,252,257,258,259,262,261,262,266,264,&
      269,268,271/
! values are in atomic units of mass    
    if ((iz<1).or.(iz>nz)) then
      write(6,*) 'weight: out of range iz =',iz
      weight = 1.0_pr
    else
      weight= real(dummy(iz),pr)
    endif
  end function weight

!> \brief returns covalent radius for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:15:07
!> \param iz integer Z of the atom

  real(pr) function covalent_radius(iz)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'covalent_radius'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    real(pr) :: dummy(1:nz)

    data dummy /0.23,1.5,0.68,0.35,0.83,0.68,0.68,0.68,0.64,1.5,0.97,1.1,1.35,&
      1.2,1.05,1.02,0.99,1.51,1.33,0.99,1.44,1.47,1.33,1.35,1.35,1.34,1.33,&
      1.5,1.52,1.45,1.22,1.17,1.21,1.22,1.21,1.5,1.47,1.12,1.78,1.56,1.48,&
      1.47,1.35,1.4,1.45,1.5,1.59,1.69,1.63,1.46,1.46,1.47,1.4,1.5,1.67,1.34,&
      1.87,1.83,1.82,1.81,1.8,1.8,1.99,1.79,1.76,1.75,1.74,1.73,1.72,1.94,1.72,&
      1.57,1.43,1.37,1.35,1.37,1.32,1.5,1.5,1.7,1.55,1.54,1.54,1.68,1.21,1.5,&
      1.5,1.9,1.88,1.79,1.61,1.58,1.55,1.53,1.51,0.99,1.54,1.83,1.5,1.5,1.5,&
      1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/
! values are in Angstroms   
    if ((iz<1).or.(iz>nz)) then
      write(6,*) 'covalent_radius: out of range iz =',iz
      covalent_radius = 2.0_pr
    else
      covalent_radius= real(dummy(iz),pr)
    endif
  end function covalent_radius  

!> \brief returns the van der Waals radius for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:14:09
!> \param iz integer Z of the atom

  real(pr) function vdw_radius(iz)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'vdw_radius'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz

    real(pr) :: dummy(1:nz)
    data dummy /1.2,1.4,1.82,2,2,1.7,1.55,1.52,1.47,1.54,2.27,1.73,2,2.1,1.8,1.8,&
      1.75,1.88,2.75,2,2,2,2,2,2,2,2,1.63,1.4,1.39,1.87,2,1.85,1.9,1.85,2.02,&
      2,2,2,2,2,2,2,2,2,1.63,1.72,1.58,1.93,2.17,2,2.06,1.98,2.16,2,2,2,2,2,2,&
      2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1.72,1.66,1.55,1.96,2.02,2,2,2,2,2,2,&
      2,2,2,1.86,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
! values are in Angstroms   
    if ((iz<1).or.(iz>nz)) then
      write(6,*) 'vdw_radius: out of range iz =',iz
      vdw_radius = 2.0_pr
    else
      vdw_radius= real(dummy(iz),pr)
    endif
  end function vdw_radius 

end module constants
