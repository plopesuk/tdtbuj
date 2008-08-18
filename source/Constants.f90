!> \brief defines constants that are needed in the code.
!> \author Alin M Elena
!> \date January 2007
module m_Constants
  implicit none
  private

  integer, parameter, public :: k_pr=kind(1.0d0) !< define precision for reals
  integer, parameter, public :: k_ml=255 !< number of character per line
  integer, parameter, public :: k_mw=40 !< number of character per word

  integer, parameter,public :: k_lowVerbos=5
  integer, parameter,public :: k_mediumVerbos=15
  integer, parameter,public :: k_highVerbos=25
!> verbosity levels

  integer, parameter,public :: k_lowDebug=5
  integer, parameter,public :: k_mediumDebug=15
  integer, parameter,public :: k_highDebug=25
!> debug levels
  integer, parameter, public :: k_runSp = 1
  integer, parameter, public :: k_runBO = 2
  integer, parameter, public :: k_runEhrenfest = 3
  integer, parameter, public :: k_runFit = 4
  integer, parameter, public :: k_runForceTest = 5
  integer, parameter, public :: k_runForceTestx = 6
  integer, parameter, public :: k_runForceTesty = 7
  integer, parameter, public :: k_runForceTestz = 8
  integer, parameter, public :: k_runEhrenfestDamped = 9
  integer, parameter, public :: k_runFragments = 10
  integer, parameter, public :: k_runGeometryOptimisation = 11
  integer, parameter, public :: k_runtestTails = 12
  integer, parameter, public :: k_runSpecial = 99
!> RunType

  integer, parameter, public :: k_scfTbuj = 1
  integer, parameter, public :: k_scfTbujo = 2
  integer, parameter, public :: k_scfTbu = 3
  integer, parameter, public :: k_scfTbuo = 4
!> SCFType

  integer, parameter, public :: k_electrostaticsPoint=1
  integer, parameter, public :: k_electrostaticsMultipoles=2
!> electrostatics

  integer, parameter, public :: k_spinDown=0
  integer, parameter, public :: k_spinUp=1
!>spins
  integer, parameter, public :: k_bondGSP=1
  integer, parameter, public :: k_bondHarrison=2
!> bond type
  integer, parameter, public :: k_smFD=1
  integer, parameter, public :: k_smMP=2
  integer, parameter, public :: k_smCS=4
  integer, parameter, public :: k_smCMU=3
!> smearing methods
  integer,parameter, public :: k_unitsEV=1
  integer,parameter, public :: k_unitsAU=2
  integer,parameter, public :: k_unitsSI=3
  integer,parameter, public :: k_unitsARU=4
!> units
  integer, parameter, public :: k_simplex=1
  integer, parameter, public :: k_SA=2
  integer, parameter, public :: k_simplexSA=3
  integer, parameter, public :: k_TrustRegion=4
!> fitting methods

  integer, parameter, public :: k_constE=1
  integer, parameter, public :: k_trigonometricE=2
  integer, parameter, public :: k_gaussianE=3
  integer, parameter, public :: k_customE=4
!> time dependent electric field methods



  integer, parameter, public :: k_wrSCF=1
  integer, parameter, public :: k_wrnSCF=2
  integer, parameter, public :: k_wrTailored=3
!> what type of density to be used in EhrenfestDamped
  integer, parameter, public :: k_lbfgs=0
  integer, parameter, public :: k_bfgs=1
!> geometry optimisation algorithms
 ! mathematics
  real(k_pr), parameter, public :: k_pi = 3.14159265358979323846264338327950_k_pr !< \f$ \pi \f$
  real(k_pr), parameter, public :: k_infinity = huge(1.0_k_pr) !< infinity
  real(k_pr), parameter, public :: k_zero=0.0_k_pr !< _zero
  real(k_pr), parameter, public :: k_one=1.0_k_pr !< _one
  complex(k_pr), parameter, public :: k_cOne=(1.0_k_pr,0.0_k_pr) !< complex one
! physics
   real(k_pr),  public :: k_amuToInternal !< atomic mass units to internal mass units
   real(k_pr),  public :: k_kb  !< Boltzmann constant \f$ k_B \f$
   real(k_pr), parameter,  public ::  k_evToHartree = 0.03674932600434263_k_pr !< electron volt to Hartree conversion factor
  real(k_pr), parameter, public ::   k_hartreeToEv = 27.211383411_k_pr !< hartree to electron volt conversion factor
  real(k_pr),  public :: k_hbar !< Plank's constant \f$ \hbar \f$
  real(k_pr), parameter, public ::   k_bohrToAng = 0.5291772083_k_pr !< bohr to  Angstroms conversion factor
  real(k_pr), parameter, public ::    k_efsToMicroamps = 160.2176462_k_pr
  real(k_pr),  public :: k_epsilon0 !< \f$ \epsilon_0 \f$ vacuum permiability
  real(k_pr),  public :: k_e !< elctronic charge
  real(k_pr),  public :: k_e2 !< square of the elctronic charge
  real(k_pr),  public :: k_me !< mass of the electron
  real(k_pr), public :: k_charge2SI !< charge to SI conversion factor
  real(k_pr), public :: k_length2SI !< length to SI conversion factor
  real(k_pr), public :: k_time2SI !< time to SI conversion factor
  real(k_pr), parameter,public :: k_debye2SI=3.335640952D-30 !< debye to SI (Cm) conversion factor
  real(k_pr), parameter,public :: k_chargeSI=1.60217653D-19 !< charge in SI (C)
  integer, parameter :: k_nz=110 !< no of atomic elements that have info associated with

  ! error messages codes for Trust region algorithm
  integer, parameter, public :: TR_SUCCESS = 1501
  integer, parameter, public :: TR_INVALID_OPTION = 1502
  integer, parameter, public :: TR_OUT_OF_MEMORY = 1503

  public :: InitializeConstants
  public :: Symbol
  public :: ElName
  public :: ElConfig
  public :: CovalentRadius
  public :: Period
  public :: VdWRadius
  public :: weight
  public :: group


contains

!> \brief Initializes the m_Constants according to the system of units chosen
!> \author Alin M Elena
!> \date 2007
!> \param units selects the system of units

   subroutine InitializeConstants(units)
      integer,intent(in) :: units
      select case(units)
      case(k_unitsEV)
        k_amuToInternal = 103.6426867_k_pr
        k_kb = 8.61734215D-5 ! in eV / K
        k_hbar = 0.65821188926_k_pr ! in eV fs
        k_epsilon0 = 5.526349954D-3
        k_e=1.0_k_pr
        k_e2=k_e*k_e
        k_me=1.0_k_pr !
        k_charge2SI= 1.60217653D-19 !Cm
        k_length2SI=1.0D-10 !m
        k_length2SI=0.048397933_k_pr !fs
      case(k_unitsAU)
        k_amuToInternal =  1836.15267247_k_pr
        k_kb = 3.1668154D-6 !Eh/K
        k_hbar = 1.0_k_pr
        k_epsilon0 = 1/(4.0_k_pr*k_pi)
        k_e=1.0_k_pr
        k_e2=k_e*k_e
        k_me=1.0_k_pr
        k_charge2SI= 1.60217653D-19 !Cm
        k_length2SI=0.5291772108D-10 !m
        k_time2SI=2.418884326505D-2 !fs
      case(k_unitsSI)
        k_amuToInternal = 1.66053886D-27 !kg
        k_kb= 1.3806505D-23 !J/K
        k_hbar = 1.05457168D-34 !J s
        k_epsilon0 = 8.854187817D-12 !F/m
        k_e=1.60217653D-19 ! C
        k_e2=k_e*k_e
        k_me= 9.1093826D-31! kg
        k_charge2SI= 1.0_k_pr !Cm
        k_length2SI=1.0_k_pr !m
        k_time2SI=1.0D-15 !fs
      case(k_unitsARU)
        k_amuToInternal =  1836.15267247_k_pr
        k_kb = 6.3361921D-6 !Rydberg/K
        k_hbar = 1.0_k_pr
        k_epsilon0 = 1.0_k_pr/(4.0_k_pr*k_pi)
        k_e=sqrt(2.0_k_pr)
        k_e2=2.0_k_pr
        k_me=0.5_k_pr
        k_charge2SI= 1.60217653D-19/sqrt(2.0_k_pr)
        k_length2SI=0.5291772108D-10 !m
        k_time2SI=4.83776865301D-2!fs
      end select

   end subroutine InitializeConstants

!> \brief returns the symbol for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:19:04
!> \param iz integer Z of the atom
  character(len=2) function Symbol(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'Symbol'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=2) :: dummy(1:k_nz)

    data dummy /"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",&
      "P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",&
      "Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru",&
      "Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr",&
      "Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W",&
      "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",&
      "Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",&
      "Rf","Db","Sg","Bh","Hs","Mt","Ds"/

    if ((iz<1).or.(iz>k_nz)) then
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
  character(len=k_mw) function ElName(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'ElName'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=k_mw) :: dummy(1:k_nz)

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

    if ((iz<1).or.(iz>k_nz)) then
      write(6,*) ' element name: out of range iz =',iz
      ElName = ' '
    else
      ElName = dummy(iz)
    endif
  end function ElName

!> \brief returns the Period from Periodic for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:17:39
!> \param iz integer Z of the atom

  character(len=4) function Period(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'Period'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=4) :: dummy(1:k_nz)

    data dummy /"1","1","2","2","2","2","2","2","2","2",&
      "3","3","3","3","3","3","3","3",&
      "4","4","4","4","4","4","4","4","4","4","4","4","4","4","4","4","4","4",&
      "5","5","5","5","5","5","5","5","5","5","5","5","5","5","5"," 5"," 5"," 5",&
      "6","6","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan","6Lan",&
      "6Lan","6Lan","6Lan","6Lan","6","6","6","6","6","6","6","6","6","6",&
      "6","6","6","6","6","6","7","7","7Act","7Act","7Act","7Act","7Act","7Act",&
      "7Act","7Act","7Act","7Act","7Act","7Act","7Act","7Act","7","7","7","7",&
      "7","7","7","7"/

    if ((iz<1).or.(iz>k_nz)) then
      write(6,*) 'Period: out of range iz =',iz
      Period= " "
    else
      Period= dummy(iz)
    endif
  end function Period

!> \brief returns the group from Periodic table for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:17:12
!> \param iz integer Z of the atom

  integer function group(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'group'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    integer :: dummy(1:k_nz)

    data dummy /1,18,1,2,13,14,15,16,17,18,1,2,13,14,15,16,17,18,1,2,3,4,5,6,7,8,&
      9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,&
      1,2,19,19,19,19,19,19,19,19,19,19,19,19,19,19,3,4,5,6,7,8,9,10,11,12,13,&
      14,15,16,17,18,1,2,20,20,20,20,20,20,20,20,20,20,20,20,20,20,3,4,5,6,7,8,&
      9,10/

    if ((iz<1).or.(iz>k_nz)) then
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
!> \internal complete the electronic configuration
  character(len=k_mw) function ElConfig(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'ElConfig'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    character(len=k_mw) :: dummy(1:k_nz)

      data dummy /"1s1","1s2","2s1","2s2","2s2 2p1","2s2 2p2","2s2 2p3","2s2 2p4","2s2 2p5",&
      "2s2 2p6","3s1","3s2","3s2 3p1","3s2 3p2","3s2 3p3",&
      "p","p","p","s","s","d","d","d","d","d","d","d","d","d","d","p","p",&
      "p","p","p","p","s","s","d","d","d","d","d","d","d","d","d","d","p",&
      "p","p","p","p","p","s","s","f","f","f","f","f","f","f","f","f","f",&
      "f","f","f","f","d","d","d","d","d","d","d","d","d","d","p","p","p",&
      "p","p","p","s","s","f","f","f","f","f","f","f","f","f","f","f","f",&
      "f","f","d","d","d","d","d","d","d","d"/

    if ((iz<1).or.(iz>k_nz)) then
      write(6,*) 'electronic configuration: out of range iz =',iz
      ElConfig = " "
    else
      ElConfig = dummy(iz)
    endif
  end function ElConfig

!> \brief returns atomic weight for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:15:38
!> \param iz integer Z of the atom

  real(k_pr) function weight(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'weight'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    real(k_pr) :: dummy(1:k_nz)

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
    if ((iz<1).or.(iz>k_nz)) then
      write(6,*) 'weight: out of range iz =',iz
      weight = 1.0_k_pr
    else
      weight= real(dummy(iz),k_pr)
    endif
  end function weight

!> \brief returns covalent radius for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:15:07
!> \param iz integer Z of the atom

  real(k_pr) function CovalentRadius(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'CovalentRadius'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz
    real(k_pr) :: dummy(1:k_nz)

    data dummy /0.23,1.5,0.68,0.35,0.83,0.68,0.68,0.68,0.64,1.5,0.97,1.1,1.35,&
      1.2,1.05,1.02,0.99,1.51,1.33,0.99,1.44,1.47,1.33,1.35,1.35,1.34,1.33,&
      1.5,1.52,1.45,1.22,1.17,1.21,1.22,1.21,1.5,1.47,1.12,1.78,1.56,1.48,&
      1.47,1.35,1.4,1.45,1.5,1.59,1.69,1.63,1.46,1.46,1.47,1.4,1.5,1.67,1.34,&
      1.87,1.83,1.82,1.81,1.8,1.8,1.99,1.79,1.76,1.75,1.74,1.73,1.72,1.94,1.72,&
      1.57,1.43,1.37,1.35,1.37,1.32,1.5,1.5,1.7,1.55,1.54,1.54,1.68,1.21,1.5,&
      1.5,1.9,1.88,1.79,1.61,1.58,1.55,1.53,1.51,0.99,1.54,1.83,1.5,1.5,1.5,&
      1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/
! values are in Angstroms
    if ((iz<1).or.(iz>k_nz)) then
      write(6,*) 'CovalentRadius: out of range iz =',iz
      CovalentRadius = 2.0_k_pr
    else
      CovalentRadius= real(dummy(iz),k_pr)
    endif
  end function CovalentRadius

!> \brief returns the van der Waals radius for element iz
!> \author Alin M Elena
!> \date 29/10/07, 23:14:09
!> \param iz integer Z of the atom

  real(k_pr) function VdWRadius(iz)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'VdWRadius'
    !--subroutine parameters--------------------------!
    integer, intent(in) :: iz

    real(k_pr) :: dummy(1:k_nz)
    data dummy /1.2,1.4,1.82,2,2,1.7,1.55,1.52,1.47,1.54,2.27,1.73,2,2.1,1.8,1.8,&
      1.75,1.88,2.75,2,2,2,2,2,2,2,2,1.63,1.4,1.39,1.87,2,1.85,1.9,1.85,2.02,&
      2,2,2,2,2,2,2,2,2,1.63,1.72,1.58,1.93,2.17,2,2.06,1.98,2.16,2,2,2,2,2,2,&
      2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1.72,1.66,1.55,1.96,2.02,2,2,2,2,2,2,&
      2,2,2,1.86,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
! values are in Angstroms
    if ((iz<1).or.(iz>k_nz)) then
      write(6,*) 'VdWRadius: out of range iz =',iz
      VdWRadius = 2.0_k_pr
    else
      VdWRadius= real(dummy(iz),k_pr)
    endif
  end function VdWRadius

end module m_Constants
