!> \brief deals with reading the input file(s)
!> \details it saves the gathered info in the right variables
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \remarks
!
module m_ReadData
  use m_Constants
  use m_Useful
  use m_Parser
  use m_Types
  use m_Gutenberg
  use ISO_FORTRAN_ENV
  use m_LinearAlgebra
  implicit none
  private
  public :: Initialize
  public :: CloseIoGeneral
  public :: CleanMemory
!
contains
!
!
!> \brief reads from ioLoc different fields of genLoc
!> \author Alin M Elena
!> \date ~2007
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to run
!> \param atomic type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine Initialize (ioLoc, genLoc, atomic, tbMod)
    character (len=*), parameter :: myname = 'Initialize'
    type (ioType), intent (inout) :: ioLoc
    type (generalType), intent (inout) :: genLoc
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tbMod
    real (k_pr) :: aux
    integer :: errno = - 1
!
    ioLoc%uerr = GetUnit ()
    ioLoc%inpErr = trim (ioLoc%inpFile) // ".err"
    open (unit=ioLoc%uerr, file=trim(ioLoc%inpErr), status="replace", action="write", iostat=errno)
    if (errno /= 0) then
      write (*,*) "I can not create error file ", trim (ioLoc%inpFile) // ".err"
      stop
    end if
    ioLoc%uinp = GetUnit ()
    open (unit=ioLoc%uinp, file=ioLoc%inpFile, status='old', action='read', iostat=errno)
!
    if (errno /= 0) then
      write (ioLoc%uerr,*) "I can not open file ", ioLoc%inpFile
      stop
    end if
!
! parse file and in the same time
! makes a minimal check of corectness
    call cpu_time (genLoc%time%Int)
!
    call ParseFile (ioLoc)
    call cpu_time (genLoc%time%end)
    aux = genLoc%time%end - genLoc%time%Int
!
    call cpu_time (genLoc%time%Int)
!
    call ReadIo (ioLoc)
!
    if (ioLoc%debug >= k_mediumDebug) then
      call cpu_time (genLoc%time%end)
      write (ioLoc%udeb, '(a,f16.6,a)') "Parsing time ", aux, " seconds"
      write (ioLoc%udeb, '(a,f16.6,a)') "Setting ioInfo ", genLoc%time%end - genLoc%time%Int, " seconds"
    end if
!
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 1)
    call ReadGeneral (ioLoc, genLoc)
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 2)
    call InitializeConstants (genLoc%units)
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 1)
    call ReadSpecies (ioLoc, genLoc, atomic%species)
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 2)
!
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 1)
    call ReadAtoms (ioLoc, genLoc, atomic)
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 2)
!
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 1)
    call ReadBasis (ioLoc, genLoc, atomic)
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 2)
!
    if ((genLoc%scfType == k_scftbuo) .or. (genLoc%scfType == k_scftbujo)) then
      if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 1)
      call ReadSpecies (ioLoc, genLoc, atomic%species, atomic%speciesBasis)
      if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 2)
    end if
!
!
    call PrintSpecies (genLoc, ioLoc, atomic%species)
    call PrintAtoms (ioLoc, genLoc, atomic)
    call PrintBasis (ioLoc, genLoc, atomic)
!
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 1)
    call ReadTBModel (ioLoc, genLoc, atomic, tbMod)
    if (genLoc%electrostatics == k_electrostaticsMultipoles) then
      call ReadDelta (ioLoc, genLoc, atomic, tbMod)
      call PrintDelta (ioLoc, genLoc, atomic, tbMod)
    end if
    if (ioLoc%debug >= k_mediumDebug) call timing (ioLoc, genLoc%time, 2)
!
    call EndParse
  end subroutine Initialize
!
!> \brief reads the names for I/O files and the output levels
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param ioLoc type(ioType)
!> \details
!> A short description of the input file ioInfo part
!
  subroutine ReadIo (ioLoc)
    character (len=*), parameter :: myname = "ReadIo"
    type (ioType), intent (inout) :: ioLoc
    integer :: errno
!comm_io DebugFile & string & input file name.dbg & Name of the debug file \\
!
! read the name of output file plus debug (y/n) and output level
    ioLoc%debugFile = GetString (ioLoc, "DebugFile", trim(ioLoc%inpFile)//".dbg", .false.)
    ioLoc%udeb = GetUnit ()
    open (unit=ioLoc%udeb, file=trim(ioLoc%debugFile), status="replace", action="write", iostat=errno)
    if (errno /= 0) call error ("I can not create file "//trim(ioLoc%debugFile), myname, .true., ioLoc)
!
    ioLoc%aniFile = GetString (ioLoc, "aniFile", trim(ioLoc%inpFile)//".xyz", .false.)
    ioLoc%uani = GetUnit ()
    open (unit=ioLoc%uani, file=trim(ioLoc%aniFile), status="replace", action="write", iostat=errno)
    if (errno /= 0) call error ("I can not create file "//trim(ioLoc%aniFile), myname, .true., ioLoc)
!
    ioLoc%outputFile = GetString (ioLoc, "OutputFile", trim(ioLoc%inpFile)//".out")
    ioLoc%stdout = GetLogical (ioLoc, "OnScreen", .false.)
! prepares the output according to the settings from input file(s)
    if (ioLoc%stdout) then
      ioLoc%uout = 6
    else
      ioLoc%uout = GetUnit ()
      open (unit=ioLoc%uout, file=trim(ioLoc%outputFile), status="replace", action="write", iostat=errno)
      if (errno /= 0) call error ("I can not create file "//trim(ioLoc%outputFile), myname, .true., ioLoc)
    end if
!
! levels
    ioLoc%verbosity = GetInteger (ioLoc, "OutputLevel", 5)
    ioLoc%debug = GetInteger (ioLoc, "DebugLevel", 5)
    ioLoc%firstTime = .not. ioLoc%firstTime
!
    if (ioLoc%verbosity >= k_highDebug) then
      write (ioLoc%udeb, '(a,l1)') "Read ioInfo for the first time: ", ioLoc%firstTime
    end if
!
    write (ioLoc%uout, '(a,a,a,i0)') "Input file: ", trim (ioLoc%inpFile), " on unit: ", ioLoc%uinp
    write (ioLoc%uout, '(a,a,a,i0)') "Output file(OutputFile): ", trim (ioLoc%outputFile), " on unit: ", ioLoc%uout
    write (ioLoc%uout, '(a,i0)') "Output Level(OutputLevel): ", ioLoc%verbosity
    write (ioLoc%uout, '(a,a,a,i0)') "Debug file(DebugFile): ", trim (ioLoc%debugFile), " on unit: ", ioLoc%udeb
    write (ioLoc%uout, '(a,i0)') "Debug Level(DebugLevel): ", ioLoc%debug
    write (ioLoc%uout, '(a,a,a,i0)') "Error file: ", trim (ioLoc%inpErr), " on unit: ", ioLoc%uerr
    write (ioLoc%uout, '(a,a,a,i0)') "Animation file: ", trim (ioLoc%aniFile), " on unit: ", ioLoc%uani
    write (ioLoc%uout, '(a,l1)') "Standard Output(OnScreen): ", ioLoc%stdout
  end subroutine ReadIo
!
!
!> \page ioVars I/O Variables
!> \latexonly
!> \begin{tabular}{|c||p{0.2\textwidth}|c||p{0.25\textwidth}|}
!>\hline
!> \textbf{VarName} & \textbf{Values} & \textbf{Default Value} & \textbf{Description} \\\
!>\hline \hline
!> DebugFile & string & input file name.dbg & Name of the debug file \\\
!> \hline
!> OutputFile & string & input file name.out & Name of the output file \\\
!>\hline
!> OnScreen & logical & .false. & on/off Printing on the screen, on means that nothing will be written in the output file\\\
!>\hline
!> DebugLevel & 5,15,25 & 5 & debug information level (low, medium,high) \\\
!>\hline
!> OutputLevel & 5,15,25 & 5 & output information level (low, medium,high) \\\
!> \hline
!> AniFile & string & animation file name.xyz & Name of the animation file \\\
!>\hline
!> \end{tabular}
!> \endlatexonly
!> \htmlonly
!> <TABLE CELLPADDING=3 BORDER="1">
!> <TR><TH ALIGN="CENTER"><B>VarName</B></TH>
!> <TH ALIGN="LEFT" VALIGN="TOP" WIDTH=100><B>Values</B></TH>
!> <TH ALIGN="CENTER"><B>Default Value</B></TH>
!> <TH ALIGN="LEFT" VALIGN="TOP" WIDTH=125><B>Description</B></TH>
!> </TR>
!> <TR><TD ALIGN="CENTER">DebugFile</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>string</TD>
!> <TD ALIGN="CENTER">input file name.dbg</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Name of the debug file</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">OutputFile</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>string</TD>
!> <TD ALIGN="CENTER">input file name.out</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Name of the output file</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">AniFile</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>string</TD>
!> <TD ALIGN="CENTER">animation file name.xyz</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Name of the animation file</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">OnScreen</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off Printing on the screen, on
!>      means that nothing will be written in the output file</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">DebugLevel</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>5,15,25</TD>
!> <TD ALIGN="CENTER">5</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>debug information level (low, medium,high)</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">OutputLevel</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>5,15,25</TD>
!> <TD ALIGN="CENTER">5</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>output information level (low, medium,high)</TD>
!> </TR>
!> </TABLE>
!> \endhtmlonly
!
!
!
!
!> \brief reads the parameters necessary for running
!> \author Alin M Elena
!> \date 21st of January, 2007
!> \param  ioLoc type(ioType) I/O details (see m_Types::ioType)
!> \param genLoc type(generalType), keeps all the general information about the parameters of the program
!> \details
!> A short description of the input file general variables part
!>
  subroutine ReadGeneral (ioLoc, genLoc)
    character (len=*), parameter :: name = "ReadGeneral"
    type (ioType), intent (inout) :: ioLoc
    type (generalType), intent (inout) :: genLoc
!
    character (len=k_mw) :: saux
    integer :: nt, errno
!
    if ( .not. genLoc%firstTime) genLoc%firstTime = .not. genLoc%firstTime
!
    genLoc%jobName = GetString (ioLoc, "JobName", "no name")
    write (ioLoc%uout, '(a,a)') "Job Name(JobName): ", genLoc%jobName
!
    genLoc%ranseed = GetReal (ioLoc, "RanSeed", 0.5_k_pr)
    write (ioLoc%uout, '(a,f0.16)') "Random number generator seed(RanSeed): ", genLoc%ranseed
!
    genLoc%ReadVelocity = GetLogical (ioLoc, "ReadVel", .false.)
    write (ioLoc%uout, '(a,l1)') "Read Velocity(ReadVel): ", genLoc%ReadVelocity
!
    genLoc%bias = GetReal (ioLoc, "VBias", 0.0_k_pr)
    write (ioLoc%uout, '(a,f0.8)') "Bias(VBias): ", genLoc%bias
!
    genLoc%maxOrbitalsPerAtom = GetInteger (ioLoc, "MaxOrbsPerAtom", 8)
    write (ioLoc%uout, '(a,i0)') "Maximum Orbitals Per Atom (MaxOrbsPerAtom): ", genLoc%maxOrbitalsPerAtom
!
    genLoc%spin = GetLogical (ioLoc, "Spin", .false.)
    write (ioLoc%uout, '(a,l1)') "spin polarisation(Spin): ", genLoc%spin
    saux = GetString (ioLoc, "Units", "AU")
    write (ioLoc%uout, '(a,a)') "System of Units (Units): ", trim (saux)
    if (cstr(trim(saux), 'AU')) then
      genLoc%units = k_unitsAU
    else if (cstr(trim(saux), 'eVA')) then
      genLoc%units = k_unitsEV
    else if (cstr(trim(saux), 'SI')) then
      genLoc%units = k_unitsSI
    else if (cstr(trim(saux), 'ARU')) then
      genLoc%units = k_unitsARU
    else
      call error ("The requested system of units is not implemented", name, .true., ioLoc)
    end if
!
    saux = GetString (ioLoc, "BondType", "Harrison")
    write (ioLoc%uout, '(a,a)') "bond type (BondType): ", trim (saux)
    if (cstr(trim(saux), 'Harrison')) then
      genLoc%bond = k_bondHarrison
    else if (cstr(trim(saux), 'GSP')) then
      genLoc%bond = k_bondGSP
    else
      call error ("The requested bond type is not implemented", name, .true., ioLoc)
    end if
!
    genLoc%embedding = GetLogical (ioLoc, "Embedding", .true.)
    write (ioLoc%uout, '(a,l1)') "has embedding(Embedding): ", genLoc%embedding
!
!comm_gen Electrostatics & PointCharges, Multipoles & PointCharges & method use to compute electrostatic interaction\\
    saux = GetString (ioLoc, "Electrostatics", "PointCharges")
    write (ioLoc%uout, '(a,a)') "Electrostatics type (Electrostatics): ", trim (saux)
    if (cstr(trim(saux), 'PointCharges')) then
      genLoc%electrostatics = k_electrostaticsPoint
    else if (cstr(trim(saux), 'Multipoles')) then
      genLoc%electrostatics = k_electrostaticsMultipoles
    else
      call error ("The requested Electrostatics type is not implemented", name, .true., ioLoc)
    end if
!
!comm_gen PrecomputeMultipoles & logical & .true. & on/off precompute multipoles\\
    genLoc%compElec = .not. GetLogical (ioLoc, "PrecomputeMultipoles", .true.)
    write (ioLoc%uout, '(a,l1)') "precompute multipoles(PrecomputeMultipoles): ", genLoc%compElec
!
    genLoc%scf = GetLogical (ioLoc, "SCF", .false.)
    write (ioLoc%uout, '(a,l1)') "Is SCF?(SCF): ", genLoc%scf
!comm_gen SCFType & TB+UJ & TB+UJ & SCF method\\ 
    saux = GetString (ioLoc, "SCFType", "TB+UJ")
    write (ioLoc%uout, '(a,a)') "SCF type(SCFType): ", trim (saux)
    if (cstr(trim(saux), 'TB+UJ')) then
      genLoc%scfType = k_scfTbuj
    else if (cstr(trim(saux), 'TB+U')) then
      genLoc%scfType = k_scfTbu
    else if (cstr(trim(saux), 'TB+UJO')) then
      genLoc%scfType = k_scftbujo
    else if (cstr(trim(saux), 'TB+UO')) then
      genLoc%scfType = k_scftbuo
    else
      call error ("The requested SCFType is not implemented", name, .true., ioLoc)
    end if
    genLoc%maxscf = GetInteger (ioLoc, "SCFSteps", 100)
    write (ioLoc%uout, '(a,i0)') "SCF Maximum number of steps(SCFSteps): ", genLoc%maxscf
    genLoc%scfMix = GetReal (ioLoc, "SCFMix", 0.85_k_pr)
    write (ioLoc%uout, '(a,f0.8)') "SCF Mixing parameter(SCFMix): ", genLoc%scfMix
    genLoc%scftol = GetReal (ioLoc, "SCFTol", 1.0e-8_k_pr)
    write (ioLoc%uout, '(a,ES12.4)') "SCF convergence tolerance(SCFTol): ", genLoc%scftol
    genLoc%scfMixn = GetInteger (ioLoc, "SCFMixN", 4)
    write (ioLoc%uout, '(a,i0)') "number of iterations to mix(SCFMixN): ", genLoc%scfMixn
    saux = GetString (ioLoc, "RunType", "SinglePoint")
    write (ioLoc%uout, '(a,a)') "Calculation type(RunType): ", trim (saux)
    if (cstr(trim(saux), 'SinglePoint')) then
      genLoc%runType = k_runSp
    else if (cstr(trim(saux), 'bodynamics')) then
      genLoc%runType = k_runBO
    else if (cstr(trim(saux), 'ehrenfest')) then
      genLoc%runType = k_runEhrenfest
    else if (cstr(trim(saux), 'fit')) then
      genLoc%runType = k_runFit
    else if (cstr(trim(saux), 'forcetest')) then
      genLoc%runType = k_runForceTest
    else if (cstr(trim(saux), 'forcetestx')) then
      genLoc%runType = k_runForceTestx
    else if (cstr(trim(saux), 'forcetesty')) then
      genLoc%runType = k_runForceTesty
    else if (cstr(trim(saux), 'forcetestz')) then
      genLoc%runType = k_runForceTestz
    else if (cstr(trim(saux), 'ehrenfestdamped')) then
      genLoc%runType = k_runEhrenfestDamped
    else if (cstr(trim(saux), 'GeometryOptimization')) then
      genLoc%runType = k_runGeometryOptimisation
    else if (cstr(trim(saux), 'TestTails')) then
      genLoc%runType = k_runtestTails
    else if (cstr(trim(saux), 'Special')) then
      genLoc%runType = k_runSpecial
    else
      call error ("The requested RunType is not implemented", name, .true., ioLoc)
    end if
!
    genLoc%hElementThreshold = GetReal (ioLoc, "HElTres", 1.0e-10_k_pr)
    write (ioLoc%uout, '(a,ES12.4)') "hamiltionian element thresold(HElThres): ", genLoc%hElementThreshold
!
    genLoc%collinear = GetLogical (ioLoc, "CollinearSpins", .false.)
    write (ioLoc%uout, '(a,l1)') "collinear spins(CollinearSpins): ", genLoc%collinear
!
    genLoc%maxIt = GetInteger (ioLoc, "MaxIt", 500)
    write (ioLoc%uout, '(a,i0)') "Maximum number of iterations to find Fermi level(MaxIt): ", genLoc%maxIt
!
    genLoc%qTolerance = GetReal (ioLoc, "ChargeTol", 1.0e-10_k_pr)
    write (ioLoc%uout, '(a,ES12.4)') "Charge Tolerance used to find Fermi level(ChargeTol): ", genLoc%qTolerance
!
    genLoc%netcharge = GetReal (ioLoc, "NetCharge", 0.0_k_pr)
    write (ioLoc%uout, '(a,ES12.4)') "Net charge(NetCharge): ", genLoc%netcharge
    saux = GetString (ioLoc, "SmearingMethod", "FD")
    write (ioLoc%uout, '(a,a)') "Smearing Methos used to find Fermi level (SmearingMethod): ", trim (saux)
    if (cstr(trim(saux), 'FD')) then
      genLoc%smearMethod = k_smFD
    else if (cstr(trim(saux), 'MP')) then
      genLoc%smearMethod = k_smMP
    else if (cstr(trim(saux), 'CS')) then
      genLoc%smearMethod = k_smCS
    else if (cstr(trim(saux), 'CMU')) then
      genLoc%smearMethod = k_smCMU
    else
      call error ("The requested SmearingMethod is not implemented", name, .true., ioLoc)
    end if
!
    select case (genLoc%smearMethod)
    case (k_smFD)
      genLoc%electronicTemperature = GetReal (ioLoc, "ElectronicTemperature", 300.0_k_pr)
      write (ioLoc%uout, '(a,ES12.4)') "Electronic Temperature(ElectronicTemperature): ", genLoc%electronicTemperature
    case (k_smMP)
      genLoc%MPW = GetReal (ioLoc, "ElectronicW", 0.005_k_pr)
      write (ioLoc%uout, '(a,ES12.4)') "Electronic W(ElectronicW): ", genLoc%MPW
      genLoc%mpN = GetInteger (ioLoc, "MPN", 2)
      write (ioLoc%uout, '(a,i0)') "Order of Hermite polynomials used find Fermi level(MPN): ", genLoc%mpN
    case (k_smCMU)
      genLoc%electronicMu = GetReal (ioLoc, "ElectronicMu", 0.0_k_pr)
      write (ioLoc%uout, '(a,ES12.4)') "Chemical Potential(ElectronicMu): ", genLoc%electronicMu
    case (k_smCS)
      genLoc%MPW = GetReal (ioLoc, "ElectronicW", 0.005_k_pr)
      write (ioLoc%uout, '(a,ES12.4)') "Electronic W(ElectronicW): ", genLoc%MPW
    end select
    genLoc%dmOccupationTolerance = GetReal (ioLoc, "DMOccTol", 1.0e-10_k_pr)
    write (ioLoc%uout, '(a,ES12.4)') "density matrix occupation tolerance(DMOccTol): ", genLoc%dmOccupationTolerance !
!
    genLoc%SymRefRho = GetLogical (ioLoc, "SymRefRho", .false.)
    write (ioLoc%uout, '(a,l1)') "Symetric Reference Density Matrix (SymRefRho): ", genLoc%SymRefRho
!
    genLoc%screened = GetLogical (ioLoc, "Screened", .false.)
    write (ioLoc%uout, '(a,l1)') "Screened Coulomb interation (Screened): ", genLoc%screened
!
    genLoc%writeAnimation = GetLogical (ioLoc, "WriteAnimation", .true.)
    write (ioLoc%uout, '(a,l1)') "Write animation(WriteAnimation): ", genLoc%writeAnimation
!
    if (genLoc%writeAnimation) then
      genLoc%AnimationSteps = GetInteger (ioLoc, "AnimationSteps", 100)
      write (ioLoc%uout, '(a,i0)') "Which steps from the animations are written (AnimationSteps): ", genLoc%AnimationSteps
    end if
!
!
    genLoc%deltat = GetReal (ioLoc, "DeltaT", 0.001_k_pr)
    write (ioLoc%uout, '(a,ES12.4)') "Time step(DeltaT): ", genLoc%deltat
!
    genLoc%nsteps = GetInteger (ioLoc, "Nsteps", 100)
    write (ioLoc%uout, '(a,i0)') "Number of steps(Nsteps): ", genLoc%nsteps
!
    genLoc%scaleVelocities = GetLogical (ioLoc, "VelScale", .false.)
    write (ioLoc%uout, '(a,l1)') "Scale  Velocities?(VelocitiesScale): ", genLoc%scaleVelocities
!
    genLoc%ionicTemperature = GetReal (ioLoc, "IonicTemperature", 300.0_k_pr)
    write (ioLoc%uout, '(a,ES12.4)') "Ionic Temperature(IonicTemperature): ", genLoc%ionicTemperature
!
    genLoc%BiasRampSteps = GetInteger (ioLoc, "BiasRampSteps",-1)
    write (ioLoc%uout, '(a,i0)') "No of Bias Ramp Steps (BiasRampSteps): ", genLoc%BiasRampSteps
!
    if (genLoc%runType == k_runEhrenfestDamped) then
      genLoc%eulerSteps = GetInteger (ioLoc, "EulerSteps", 100)
      write (ioLoc%uout, '(a,i0)') "Euler steps (EulerSteps): ", genLoc%eulerSteps
      genLoc%gamma = GetReal (ioLoc, "Gamma", 0.5_k_pr)
      write (ioLoc%uout, '(a,ES12.4)') "damping factor in Ehrenfest Equation(Gamma): ", genLoc%gamma
      saux = GetString (ioLoc, "WhatDensity", "SCF")
      write (ioLoc%uout, '(a,a)') "What Density (WhatDensity) = ", trim (saux)
      if (cstr(trim(saux), "SCF")) then
        genLoc%wdensity = k_wrSCF
      else if (cstr(trim(saux), "nonSCF")) then
        genLoc%wdensity = k_wrnSCF
      else if (cstr(trim(saux), "Tailored")) then
        genLoc%wdensity = k_wrTailored
      else
        call error ("The requested density matrix method is not implemented", name, .true., ioLoc)
      end if
    end if
    if (genLoc%runType == k_runFit) then
      saux = GetString (ioLoc, "FitMethod", "Simplex")
      write (ioLoc%uout, '(a,a)') "Fit method = ", trim (saux)
      if (cstr(trim(saux), "simplex")) then
        genLoc%fit%fitMethod = k_simplex
      else if (cstr(trim(saux), "SA")) then
        genLoc%fit%fitMethod = k_sa
      else if (cstr(trim(saux), "simplexSA")) then
        genLoc%fit%fitMethod = k_simplexSA
      else if (cstr(trim(saux), "TrustRegion")) then
        genLoc%fit%fitMethod = k_TrustRegion
      else
        call error ("The requested fitting method is not implemented", name, .true., ioLoc)
      end if
      genLoc%fit%iNoParams = GetInteger (ioLoc, "NoFitParams", 0)
      write (ioLoc%uout, '(a,i0)') "No of parameters to fit = ", genLoc%fit%iNoParams
      if (genLoc%fit%iNoParams <= 0) then
        call error ("NoFitParams has to be positive", name, .true., ioLoc)
      end if
!
      genLoc%fit%fitTol = GetReal (ioLoc, "FitTol", 1.0e-4_k_pr)
      write (ioLoc%uout, '(a,ES12.4)') "Fit Tolerance = ", genLoc%fit%fitTol
!
      genLoc%fit%restartFit = GetLogical (ioLoc, "RestartFit", .false.)
      write (ioLoc%uout, '(a,l1)') "Restart Fit = ", genLoc%fit%restartFit
!
      genLoc%fit%neps = GetInteger (ioLoc, "FitNeps", 4)
      write (ioLoc%uout, '(a,i8)') "Number of terms used to compute termination criteria= ", genLoc%fit%neps
      select case (genLoc%fit%fitMethod)
      case (k_simplex)
!
        genLoc%fit%iter = GetInteger (ioLoc, "SimplexMaxIter", 10000)
        write (ioLoc%uout, '(a,i8)') "MAximum number of iterations for simplex = ", genLoc%fit%iter
      case (k_sa)
!
        genLoc%fit%feval = GetInteger (ioLoc, "SAMaxFeval", 100000)
        write (ioLoc%uout, '(a,i8)') "Maximum number of cost function evaluations = ", genLoc%fit%feval
!
        genLoc%fit%ns = GetInteger (ioLoc, "SACycles", 20)
        write (ioLoc%uout, '(a,i8)') "Maximum number of cycles before step adjustment  = ", genLoc%fit%ns
!
        genLoc%fit%nt = GetInteger (ioLoc, "SAnBefRed", 5)
        write (ioLoc%uout, '(a,i8)') "Number of iterations before temperature reduction = ", genLoc%fit%nt
!
        genLoc%fit%temp = GetReal (ioLoc, "SATemp", 100.0_k_pr)
        write (ioLoc%uout, '(a,f8.3)') "SA temperature = ", genLoc%fit%temp
!
        genLoc%fit%step = GetReal (ioLoc, "SAstep", 1.0_k_pr)
        write (ioLoc%uout, '(a,f8.3)') "SA initial step = ", genLoc%fit%step
!
        genLoc%fit%stepAd = GetReal (ioLoc, "SAStepAdj", 2.0_k_pr)
        write (ioLoc%uout, '(a,f8.3)') "SA step adjustment = ", genLoc%fit%stepAd
!
        genLoc%fit%rt = GetReal (ioLoc, "SATempRed", 0.5_k_pr)
        write (ioLoc%uout, '(a,f8.3)') "SA temperature reduction factor = ", genLoc%fit%rt
      case (k_simplexSA)
        genLoc%fit%iter = GetInteger (ioLoc, "SimplexMaxIter", 10000)
        write (ioLoc%uout, '(a,i8)') "MAximum number of iterations for simplex = ", genLoc%fit%iter
!
        genLoc%fit%rt = GetReal (ioLoc, "SATempRed", 0.5_k_pr)
        write (ioLoc%uout, '(a,f8.3)') "SA temperature reduction factor = ", genLoc%fit%rt
!
        genLoc%fit%temp = GetReal (ioLoc, "SATemp", 100.0_k_pr)
        write (ioLoc%uout, '(a,f8.3)') "SA temperature = ", genLoc%fit%temp
      case (k_TrustRegion)
        genLoc%fit%step = GetReal (ioLoc, "TRstep", 100.0_k_pr)
        write (ioLoc%uout, '(a,f8.3)') "TR initial step = ", genLoc%fit%step
        genLoc%fit%iter = GetInteger (ioLoc, "TrMaxIter", 10000)
        write (ioLoc%uout, '(a,i8)') "MAximum number of iterations for TrustRegion = ", genLoc%fit%iter
      end select
    end if
!
    if (genLoc%runType == k_runForceTest) then
      genLoc%fdx = GetReal (ioLoc, "DerivStep", 0.001_k_pr)
      write (ioLoc%uout, '(a,f8.3)') "Step for numerical derivative (DerivStep) = ", genLoc%fdx
    end if
!
    if ((genLoc%runType == k_runForceTestx) .or. (genLoc%runType == k_runForceTesty) .or. (genLoc%runType == k_runForceTestz)) then
      genLoc%fdx = GetReal (ioLoc, "DerivStep", 0.001_k_pr)
      write (ioLoc%uout, '(a,f8.3)') "Step for numerical derivative (DerivStep) = ", genLoc%fdx
      genLoc%fstart = GetReal (ioLoc, "ForceStart", 0.0_k_pr)
      write (ioLoc%uout, '(a,f8.3)') "Start for numerical derivative (ForceStart) = ", genLoc%fstart
      genLoc%fsteps = GetInteger (ioLoc, "ForceSteps", 100)
      write (ioLoc%uout, '(a,i0)') "No of steps for sampling (ForceSteps) = ", genLoc%fsteps
      genLoc%fatom = GetInteger (ioLoc, "ForceOnAtom", 1)
      write (ioLoc%uout, '(a,i0)') "Compute force on atom (ForceOnAtom) = ", genLoc%fatom
    end if
!
!
! !comm_gen Excite & logical & .false. & on/off create an excited state \\
    genLoc%lIsExcited = GetLogical (ioLoc, "Excite", .false.)
    write (ioLoc%uout, '(a,l1)') "Create excitation(Excite): ", genLoc%lIsExcited
!
    if (genLoc%lIsExcited) then
! !comm_gen HoleState & integer & 0 & no of level where we create the hole\\  
      genLoc%holeState = GetInteger (ioLoc, "HoleState", 0)
      write (ioLoc%uout, '(a,i0)') "Hole level(HoleState): ", genLoc%holeState
!
! !comm_gen Excite & integer & 0 & no of level to create excitation\\  
      genLoc%exciteState = GetInteger (ioLoc, "ExciteState", 0)
      write (ioLoc%uout, '(a,i0)') "Excitation level(ExciteState): ", genLoc%exciteState
! !comm_gen HoleSpin & D,U & D & spin of the hole\\  
      saux = GetString (ioLoc, "HoleSpin", "D")
      write (ioLoc%uout, '(a,a)') "Spin of the Hole (HoleSpin): ", trim (saux)
      if (cstr(trim(saux), 'D')) then
        genLoc%holeSpin = k_spinDown
      else if (cstr(trim(saux), 'U')) then
        genLoc%holeSpin = k_spinUp
      else
        call error ("The requested spin type is not implemented", name, .true., ioLoc)
      end if
!
!comm_gen ExciteSpin & D,U & D & spin of the excitation\\  
      saux = GetString (ioLoc, "ExciteSpin", "D")
      write (ioLoc%uout, '(a,a)') "Spin of the excitation (ExciteSpin): ", trim (saux)
      if (cstr(trim(saux), 'D')) then
        genLoc%exciteSpin = k_spinDown
      else if (cstr(trim(saux), 'U')) then
        genLoc%exciteSpin = k_spinUp
      else
        call error ("The requested spin type is not implemented", name, .true., ioLoc)
      end if
    end if
    if (genLoc%runType == k_runGeometryOptimisation) then
      saux = GetString (ioLoc, "GeomOptAlg", "LBFGS")
      write (ioLoc%uout, '(a,a)') "Geometry Optimization algorithm (GeomOptAlg): ", trim (saux)
      if (cstr(trim(saux), "LBFGS")) then
        genLoc%geomAlg = k_lbfgs
      else if (cstr(trim(saux), "BFGS")) then
        genLoc%geomAlg = k_bfgs
      else
        call error ("The requested fitting method is not implemented", name, .true., ioLoc)
      end if
!
      select case (genLoc%geomAlg)
      case (k_lbfgs)
        genLoc%fTol = GetReal (ioLoc, "EnergyTolerance", 1.0e-4_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Energy Tolerance (during geometry optimisation)(EnergyTolerance): ", genLoc%fTol
        genLoc%gTol = GetReal (ioLoc, "ForceTolerance", 1.0e-10_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Force Tolerance (during geometry optimisation)(ForceTolerance): ", genLoc%gTol
        genLoc%xTol = GetReal (ioLoc, "XTolerance", 1.0e-10_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') " Coordinates Tolerance (during geometry optimisation)(XTolerance): ", genLoc%xTol
        genLoc%epsf = GetReal (ioLoc, "EpsilonE", 1.0e-10_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Epsilon Energy Tolerance (during geometry optimisation)(EpsilonE): ", genLoc%epsf
        genLoc%epsg = GetReal (ioLoc, "EpsilonG", 1.0e-10_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Epsilon Gradient Tolerance (during geometry optimisation)(EpsilonG): ", genLoc%epsg
        genLoc%epsX = GetReal (ioLoc, "EpsilonX", 1.0e-10_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Epsilon Coordinates Tolerance (during geometry optimisation)(EpsilonX): ", genLoc%epsX
        genLoc%maxFEval = GetInteger (ioLoc, "MaxFEval", 20)
        write (ioLoc%uout, '(a,i8)') "Maximum number of evaluations for energy during each step of geometry optimisation (MaxFEval)&
       &", genLoc%maxFEval
        genLoc%nsteps = GetInteger (ioLoc, "GeometryNSteps", 100)
        write (ioLoc%uout, '(a,i0)') "Number of steps for geometry optimisation (GeometryNSteps): ", genLoc%nsteps
        genLoc%HessianM = GetInteger (ioLoc, "HessianM", 7)
        write (ioLoc%uout, '(a,i8)') "Number of corrections used (3<= m <=7) by bfgs to compute hessian geometry optimisation (Hess&
       &ianM) ", genLoc%HessianM
      case (k_bfgs)
        genLoc%fTol = GetReal (ioLoc, "EnergyTolerance", 1.0e-4_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Energy Tolerance (during geometry optimisation)(EnergyTolerance): ", genLoc%fTol
        genLoc%gTol = GetReal (ioLoc, "ForceTolerance", 1.0e-10_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Force Tolerance (during geometry optimisation)(ForceTolerance): ", genLoc%gTol
        genLoc%xTol = GetReal (ioLoc, "XTolerance", 1.0e-10_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') " Coordinates Tolerance (during geometry optimisation)(XTolerance): ", genLoc%xTol
        genLoc%nsteps = GetInteger (ioLoc, "GeometryNSteps", 100)
        write (ioLoc%uout, '(a,i0)') "Number of steps for geometry optimisation (GeometryNSteps): ", genLoc%nsteps
      end select
    end if
!
    genLoc%hasElectricField = GetLogical (ioLoc, "hasElectricField", .false.)
    write (ioLoc%uout, '(a,l1)') "Do we apply an external electric field(hasElectricField): ", genLoc%hasElectricField
!
    if (genLoc%hasElectricField) then
      if (GetBlock(ioLoc, "ElectricField", nt)) then
        read (nt, fmt=*, iostat=errno) genLoc%E(1), genLoc%E(2), genLoc%E(3)
        if (errno /= 0) then
          call error ("block ElectricField is not in the right format", name, .true., ioLoc)
        end if
      else
        call error ("block ElectricField is missing (you try to apply an external electric field but you do not specify it)", name, &
       & .true., ioLoc)
      end if
      write (ioLoc%uout, '(a,3f16.8)') "External Electric Field(block ElectricField): ", genLoc%E
      saux = GetString (ioLoc, "EtimeDependent", "Constant")
      write (ioLoc%uout, '(a,a)') "Time dependent part of the electric field (EtimeDependent): ", trim (saux)
      if (cstr(trim(saux), "Constant")) then
        genLoc%etd = k_constE
      else if (cstr(trim(saux), "Trigonometric")) then
        genLoc%etd = k_trigonometricE
      else if (cstr(trim(saux), "Gaussian")) then
        genLoc%etd = k_gaussianE
      else if (cstr(trim(saux), "Custom")) then
        genLoc%etd = k_customE
      else
        genLoc%etd = k_constE
      end if
      select case (genLoc%etd)
      case (k_constE, k_customE)
      case (k_trigonometricE)
        genLoc%freq = GetReal (ioLoc, "FreqE", 1.0_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Frequency of the field)(freqE): ", genLoc%freq
        genLoc%phi0 = GetReal (ioLoc, "Phi0E", 0.0_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "Initial phase of the field)(Phi0E): ", genLoc%phi0
      case (k_gaussianE)
        genLoc%t0 = GetReal (ioLoc, "t0E", 0.0_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "at what time we apply the fiield)(t0E): ", genLoc%t0
        genLoc%sigma = GetReal (ioLoc, "sigmaE", 0.001_k_pr)
        write (ioLoc%uout, '(a,ES12.4)') "the width of the gaussian for the fiield)(sigmaE): ", genLoc%sigma
      end select
    else
      genLoc%E = 0.0_k_pr
      genLoc%etd = k_constE
    end if
!
!
  end subroutine ReadGeneral
!
!> \page control
!> \latexonly
!> \begin{longtable}{|c||p{0.2\textwidth}|c||p{0.25\textwidth}|}
!> \hline
!> \textbf{VarName} & \textbf{Values} & \textbf{Default Value} & \textbf{Description} \\\
!> \hline \hline
!> JobName & string & no name & a name for the job \\\
!> \hline
!> RanSeed & real(0,1)& 0.5 & A seed for the random number generator \\\
!> \hline
!> ReadVel & logical & .false. & on/off reading velocity block \\\
!> \hline
!> SCF & logical & .false. & on/off self consistent field method \\\
!> \hline
!> VBias & real & 0.0 & bias factor\\\
!> \hline
!> MaxOrbsPerAtom & integer & 8 & maximum number of orbitals per atom\\\
!> \hline
!> Units & AU, EVA, SI, ARU & AU & system of units: atomic units(AU), electronVolt-Angstrom(eVA),
!> International(SI) Rydeberg Atomic Units (ARU)\\\
!> \hline
!> BondType & Harrison, GSP &Harrison& bond type \\\
!> \hline
!> Embedding & logical & .true. & on/off embedding method\\\
!> \hline
!> Electrostatics & PointCharges, Multipoles & PointCharges & method use to compute electrostatic interaction\\\
!> \hline
!> PrecomputeMultipoles & logical & .true. & on/off precompute multipoles\\\
!> \hline
!> SCFType & TB+UJ, TB+UJO,TB+U, TB+UO & TB+UJ & SCF method TB+UJ is for average U and J (spin calculations only) TD+U average U for
!> \hline
!> SCFSteps & integer & 100 & maximum number of steps used for scf\\\
!> \hline
!> SCFMix & real & 0.85 & mixing parameter\\\
!> \hline
!> SCFTol & real & 1e-8 & convergence tolerance\\\
!> \hline
!> SCFMixN & integer & 4 & number of iterations to mix\\\
!> \hline
!> RunType & SinglePoint, BODynamics, Ehrenfest, EhrenfestDamped, Fit, GeometryOptimization, ForceTest, ForceTestX, ForceTestY, Forc
!> \hline
!> HElThres & real & 1e-10 & hamiltionian element thresold. Any element smaller that the thresold is made zero.\\\
!> \hline
!> CollinearSpins & logical & .false. & on/off collinear spins\\\
!> \hline
!> MaxIt & integer & 500 & maximum number of iterations used to find Fermi level\\\
!> \hline
!> ChargeTol & real & 1e-10 & charge tolerance used to find Fermi level\\\
!> \hline
!> NetCharge & real & 0.0 & Net charge on the system\\\
!> \hline
!> ElectronicTemperature & real & 300.0 & Electronic temperature, used to compute occupation numbers
!> if you choose Fermi-Dirac method\\\
!> \hline
!> ElectronicW & real & 0.05 & Electronic W, used to compute occupation numbers if you choose Methfessel-Paxton method\\\
!>  \hline
!> MPN & integer & 2 & the order of Hermite polynomials used to find Fermi level by Methfessel-Paxton method \\\
!> \hline
!> ElectronicMu & real & 0.0 & chemical potential, used to compute occupation numbers if you choose constant $\mu$ method\\\
!>  \hline
!> DMOccTol & real & 1e-10 & density matrix occupation tolerance\\\
!> \hline
!> SymRefRho & logical & .false. & on/off symmetric reference density matrix\\\
!> \hline
!> Screened & logical & .false. & on/off bare Coulomb or screened electrostatic interaction \\\
!> \hline
!> WriteAnimation & logical & .true. & on/off writing animation file \\\
!> \hline
!> DeltaT & real & 0.001 & time step used to evolve equtions of motion\\\
!> \hline
!> Nsteps & integer & 100 & number of steps used to evolve equtions of motion\\\
!> \hline
!> VelScale & logical & .false. & on/off scaling velocities\\\
!> \hline
!> IonicTemperature & real & 300.0 & Ionic temperature\\\
!> \hline
!> BiasRampSteps & integer & 0 & for how many steps to apply bias\\\
!> \hline
!> EulerSteps & integer & 100 & after each EulerSteps apply an Euler integration of equations of motions\\\
!>  \hline
!> Gamma & real & 0.3 & dumping factor for Ehrenfest equation\\\
!> \hline
!> FitType & Simplex,SA, SimplexSA & Simplex & fitting method (Simplex, Simplex-Simulated Annealing, Simulated Annealing)\\\
!> \hline
!> NoFitParams & integer & 0 & numbers of parameters to fit\\\
!> \hline
!> FitTol & real & 1.0e-4 & the tolerance accepted in fit\\\
!> \hline
!> RestartFit & logical & .false. & on/off restarting the fit from a previous point\\\
!> \hline
!> FitNeps & integer & 4 & number of terms used to compute the termination criteria\\\
!> \hline
!> SimplexMaxIter & integer & 10000 & maximum number of iterations in the simplex method\\\
!> \hline
!> SAMaxFeval & integer & 100000 & maximum number of cost function evaluations in SA method\\\
!> \hline
!> SACyscles & integer & 20 & maximum number of cycles in SA method\\\
!> \hline
!> SAnBefRed & integer & 5 & number of iterations before temperature reduction in SA method\\\
!> \hline
!> SATemp & real & 100.0 & iinitial temperature in SA method\\\
!> \hline
!> SAStep & real & 1.0 & SA initial step\\\
!> \hline
!> SAStepAdj & real & 2.0 & SA step adjustment\\\
!> \hline
!> SATempRed & real & 0.5 & SA temperature reduction factor\\\
!> \hline
!> SimplexMaxIter & integer & 10000 & maximum number of iterations in the simplex method\\\
!> \hline
!> SATempRed & real & 0.5 & SA temperature reduction factor\\\
!> \hline
!> SATemp & real & 100.0 & initial temperature in SA method\\\
!> \hline
!> ForceTolerance & real & 10E-10 & Force tolerance for geometry optimisation\\\
!> \hline
!> GeomOptAlg & LBFGS, BFGS & LBFGS &Geometry Optimisation Algorithm Linear BFGS or BFGS (Broyden-Fletcher-Goldfarb-Shanno))\\\
!> \hline
!> DerivStep & real & 10E-3 & step for numerical derivatives\\\
!> \hline
!> \end{longtable}
!> \endlatexonly
!> \htmlonly
!><TABLE CELLPADDING=3 BORDER="1">
!><TR><TH ALIGN="CENTER"><B>VarName</B></TH>
!><TH ALIGN="LEFT" VALIGN="TOP" WIDTH=100><B>Values</B></TH>
!><TH ALIGN="CENTER"><B>Default Value</B></TH>
!><TH ALIGN="LEFT" VALIGN="TOP" WIDTH=125><B>Description</B></TH>
!></TR>
!><TR><TD ALIGN="CENTER">JobName</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>string</TD>
!><TD ALIGN="CENTER">no name</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>a name for the job</TD>
!></TR>
!><TR><TD ALIGN="CENTER">RanSeed</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real(0,1)</TD>
!><TD ALIGN="CENTER">0.5</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>A seed for the random number generator</TD>
!></TR>
!><TR><TD ALIGN="CENTER">ReadVel</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.false.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off reading velocity block</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SCF</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.false.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off self consistent field method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">VBias</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>bias factor</TD>
!></TR>
!><TR><TD ALIGN="CENTER">MaxOrbsPerAtom</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">8</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of orbitals per atom</TD>
!></TR>
!><TR><TD ALIGN="CENTER">Units</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>AU, eVA, SI, ARU</TD>
!><TD ALIGN="CENTER">AU</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>system of units: atomic units(AU), electronVolt-Angstrom(eVA),
!> International(SI), Rydberg Atomic Units(ARU)</TD>
!></TR>
!><TR><TD ALIGN="CENTER">BondType</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>Harrison, GSP</TD>
!><TD ALIGN="CENTER">Harrison</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>bond type</TD>
!></TR>
!><TR><TD ALIGN="CENTER">Embedding</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.true.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off embedding method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">Electrostatics</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>PointCharges, Multipoles</TD>
!><TD ALIGN="CENTER">PointCharges</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>method use to compute electrostatic interaction</TD>
!></TR>
!><TR><TD ALIGN="CENTER">PrecomputeMultipoles</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.true.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off precompute multipoles</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SCFType</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>TB+UJ, TB+UJO,TB+U, TB+UO</TD>
!><TD ALIGN="CENTER">TB+UJ</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>SCF method TB+UJ is for average U and J (spin calculations only) TD+U average U for non-sp
!></TR>
!><TR><TD ALIGN="CENTER">SCFSteps</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">100</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of steps used for scf</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SCFMix</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.85</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>mixing parameter</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SCFTol</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">1e-8</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>convergence tolerance</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SCFMixN</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">4</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>number of iterations to mix</TD>
!></TR>
!><TR><TD ALIGN="CENTER">RunType</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>SinglePoint, BODynamics, Ehrenfest, EhrenfestDamped, Fit, GeometryOptimization, ForceTest,
!><TD ALIGN="CENTER">SinglePoint</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Type of calculation</TD>
!></TR>
!><TR><TD ALIGN="CENTER">HElThres</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">1e-10</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>hamiltionian element thresold. Any element smaller that the thresold is made zero.</TD>
!></TR>
!><TR><TD ALIGN="CENTER">CollinearSpins</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.false.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off collinear spins</TD>
!></TR>
!><TR><TD ALIGN="CENTER">MaxIt</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">500</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of iterations used to find Fermi level</TD>
!></TR>
!><TR><TD ALIGN="CENTER">ChargeTol</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">1e-10</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>charge tolerance used to find Fermi level</TD>
!></TR>
!><TR><TD ALIGN="CENTER">NetCharge</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Net charge on the system</TD>
!></TR>
!><TR><TD ALIGN="CENTER">ElectronicTemperature</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">300.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Electronic temperature, used to compute occupation numbers if you
!> choose Fermi-Dirac method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">ElectronicW</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.05</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Electronic W, used to compute occupation numbers if you choose
!> Methfessel-Paxton method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">MPN</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">2</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>the order of Hermite polynomials used to find Fermi level by Methfessel-Paxton method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">ElectronicMu</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>chemical potential, used to compute occupation numbers if you choose constant  method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">DMOccTol</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">1e-10</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>density matrix occupation tolerance</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SymRefRho</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.false.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off symmetric reference density matrix</TD>
!></TR>
!><TR><TD ALIGN="CENTER">Screened</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.false.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off bare Coulomb or screened electrostatic interaction</TD>
!></TR>
!><TR><TD ALIGN="CENTER">WriteAnimation</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.true.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off writing animation file</TD>
!></TR>
!><TR><TD ALIGN="CENTER">DeltaT</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.001</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>time step used to evolve equtions of motion</TD>
!></TR>
!><TR><TD ALIGN="CENTER">Nsteps</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">100</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>number of steps used to evolve equtions of motion</TD>
!></TR>
!><TR><TD ALIGN="CENTER">VelScale</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.false.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off scaling velocities</TD>
!></TR>
!><TR><TD ALIGN="CENTER">IonicTemperature</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">300.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Ionic temperature</TD>
!></TR>
!><TR><TD ALIGN="CENTER">BiasRampSteps</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>for how many steps to apply bias</TD>
!></TR>
!><TR><TD ALIGN="CENTER">EulerSteps</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">100</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>after each EulerSteps apply an Euler integration of equations of motions</TD>
!></TR>
!><TR><TD ALIGN="CENTER">Gamma</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.3</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>dumping factor for Ehrenfest equation</TD>
!></TR>
!><TR><TD ALIGN="CENTER">FitType</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>Simplex,SA, SimplexSA</TD>
!><TD ALIGN="CENTER">Simplex</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>fitting method (Simplex, Simplex-Simulated Annealing, Simulated Annealing)</TD>
!></TR>
!><TR><TD ALIGN="CENTER">NoFitParams</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>numbers of parameters to fit</TD>
!></TR>
!><TR><TD ALIGN="CENTER">FitTol</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">1.0e-4</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>the tolerance accepted in fit</TD>
!></TR>
!><TR><TD ALIGN="CENTER">RestartFit</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!><TD ALIGN="CENTER">.false.</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off restarting the fit from a previous point</TD>
!></TR>
!><TR><TD ALIGN="CENTER">FitNeps</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">4</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>number of terms used to compute the termination criteria</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SimplexMaxIter</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">10000</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of iterations in the simplex method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SAMaxFeval</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">100000</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of cost function evaluations in SA method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SACyscles</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">20</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of cycles in SA method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SAnBefRed</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">5</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>number of iterations before temperature reduction in SA method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SATemp</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">100.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>iinitial temperature in SA method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SAStep</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">1.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>SA initial step</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SAStepAdj</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">2.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>SA step adjustment</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SATempRed</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.5</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>SA temperature reduction factor</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SimplexMaxIter</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!><TD ALIGN="CENTER">10000</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of iterations in the simplex method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SATempRed</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">0.5</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>SA temperature reduction factor</TD>
!></TR>
!><TR><TD ALIGN="CENTER">SATemp</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">100.0</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>initial temperature in SA method</TD>
!></TR>
!><TR><TD ALIGN="CENTER">ForceTolerance</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">10E-10</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>force tolerance for geometry optimisation</TD>
!></TR>
!><TR><TD ALIGN="CENTER">GeomOptAlg</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>LBFGS, BFGS</TD>
!><TD ALIGN="CENTER">LBFGS</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Geometry Optimisation Algorithm Linear BFGS or BFGS (Broyden-Fletcher-Goldfarb-Shanno))</T
!></TR>
!><TR><TD ALIGN="CENTER">DerivStep</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!><TD ALIGN="CENTER">10E-3</TD>
!><TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>step used for numercal derivatives</TD>
!></TR>
!></TABLE>
!> \endhtmlonly
!
!> \brief closes all the I/O units
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param ioLoc type(ioType)
!
  subroutine CloseIoGeneral (ioLoc)
    character (len=*), parameter :: myname = 'CloseIoGeneral'
    type (ioType), intent (inout) :: ioLoc
    integer :: errno
!
    close (unit=ioLoc%udeb, iostat=errno)
    if (errno /= 0) call error ("I can not close file "//trim(ioLoc%debugFile), myname, .false., ioLoc)
    close (unit=ioLoc%uerr, iostat=errno)
    if (errno /= 0) call error ("I can not close file "//trim(ioLoc%inpErr), myname, .false., ioLoc)
    close (unit=ioLoc%uani, iostat=errno)
    if (errno /= 0) call error ("I can not close file "//trim(ioLoc%aniFile), myname, .false., ioLoc)
!
    if ( .not. ioLoc%stdout) then
      close (unit=ioLoc%uout, iostat=errno)
      if (errno /= 0) call error ("I can not close file "//trim(ioLoc%outputFile), myname, .false., ioLoc)
    end if
  end subroutine CloseIoGeneral
!
!
!> \brief times different events
!> \details it should be called in pairs with stage = 1 to Initialize
!> and stage = 2 to write down the info
!> \author Alin M Elena
!> \date 29th of October 2007
!> \param io type(ioType) i/o units
!> \param times type(timeType) time keepers
!> \param stage selects the stage
!
  subroutine timing (io, times, stage)
    character (len=*), parameter :: myname = "timing"
    type (ioType), intent (in) :: io
    type (timeType), intent (inout) :: times
    integer, intent (in) :: stage
!
    if (stage == 1) then
      call cpu_time (times%end)
    else if (stage == 2) then
      call cpu_time (times%end)
      write (io%udeb, '(a,f16.6,a)') "parsing and setting ", times%end - times%Int, " seconds"
    end if
!
  end subroutine timing
!
!
!> \brief reads and allocates the info about atoms
!> \details this is all that should be done in this module.
!> the rest should happen in a specialized module
!> \author Alin M Elena
!> \date 29/10/07, 17:33:50
!> \param io type(ioType) i/o units
!> \param general type(generalType) general data
!> \param atomix type(atomicxType) contains info about atoms
!
  subroutine ReadAtoms (io, general, atomix)
    character (len=*), parameter :: sMyName = "ReadAtoms"
    type (ioType), intent (inout) :: io
    type (generalType), intent (in) :: general
    type (atomicxType), intent (inout) :: atomix
! arguments of the subroutine    
    integer :: nt, errno, i, k, atomId
    character (len=k_ml) :: saux
    integer, allocatable :: ids (:)
    real (k_pr) :: vx, vy, vz
!
    atomix%atoms%natoms = GetInteger (io, "NumberOfAtoms",-1)
    if (atomix%atoms%natoms <= 0) then
      call error (sMyName, "NumberOfAtoms token is missing or is negative!", .true., io)
    end if
    atomix%atoms%created = .true.
    allocate (atomix%atoms%id(1:atomix%atoms%natoms))
    allocate (atomix%atoms%sp(1:atomix%atoms%natoms))
    allocate (atomix%atoms%bias(1:atomix%atoms%natoms))
    allocate (atomix%atoms%isscf(1:atomix%atoms%natoms))
    allocate (atomix%atoms%ismoving(1:atomix%atoms%natoms))
    allocate (atomix%atoms%x(1:atomix%atoms%natoms))
    allocate (atomix%atoms%y(1:atomix%atoms%natoms))
    allocate (atomix%atoms%z(1:atomix%atoms%natoms))
    allocate (atomix%atoms%vx(1:atomix%atoms%natoms))
    allocate (atomix%atoms%vy(1:atomix%atoms%natoms))
    allocate (atomix%atoms%vz(1:atomix%atoms%natoms))
    allocate (atomix%atoms%fx(1:atomix%atoms%natoms))
    allocate (atomix%atoms%fy(1:atomix%atoms%natoms))
    allocate (atomix%atoms%fz(1:atomix%atoms%natoms))
    allocate (atomix%atoms%dx(1:atomix%atoms%natoms))
    allocate (atomix%atoms%dy(1:atomix%atoms%natoms))
    allocate (atomix%atoms%dz(1:atomix%atoms%natoms))
    allocate (atomix%atoms%xo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%yo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%zo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%vxo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%vyo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%vzo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%fxo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%fyo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%fzo(1:atomix%atoms%natoms))
    allocate (atomix%atoms%chrg(1:atomix%atoms%natoms))
    allocate (atomix%atoms%chrg0(1:atomix%atoms%natoms))
    allocate (atomix%atoms%norbs(1:atomix%atoms%natoms))
    allocate (atomix%atoms%MagMom(1:atomix%atoms%natoms))
    allocate (atomix%atoms%neighbours(1:atomix%atoms%natoms))
!
    atomix%atoms%id (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%sp (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%bias (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%x (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%y (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%z (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%vx (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%vy (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%vz (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%vxo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%vyo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%vzo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%xo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%yo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%zo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%dx (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%dy (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%dz (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%fx (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%fy (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%fz (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%fxo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%fyo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%fzo (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%chrg (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%chrg0 (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%ismoving (1:atomix%atoms%natoms) = .true.
    atomix%atoms%isscf (1:atomix%atoms%natoms) = .true.
    atomix%atoms%MagMom (1:atomix%atoms%natoms) = 0.0_k_pr
    atomix%atoms%norbs (1:atomix%atoms%natoms) = 0
    atomix%atoms%nmoving = 0
    atomix%atoms%nscf = 0
!
!
    if (GetBlock(io, "AtomsData", nt)) then
      do i = 1, atomix%atoms%natoms
        read (nt, fmt=*, iostat=errno) atomix%atoms%sp(i), atomix%atoms%x(i), atomix%atoms%y(i), atomix%atoms%z(i), &
       & atomix%atoms%bias(i), atomix%atoms%isscf(i), atomix%atoms%ismoving(i)
        if (errno /= 0) then
          call error ("block AtomsData is not in the right format", sMyName, .true., io)
        end if
        atomix%atoms%id (i) = i
        write (io%udeb, '(i6,4f16.6,1x,l1,1x,l1)') atomix%atoms%sp(i), atomix%atoms%x(i), atomix%atoms%y(i), atomix%atoms%z(i), &
       & atomix%atoms%bias(i), atomix%atoms%isscf(i), atomix%atoms%ismoving(i)
        if ( .not. isInList(atomix%atoms%sp(i), atomix%species%id)) then
          write (saux,*) "block AtomsData contains undefined specie: ", atomix%atoms%sp(i)
          call error (trim(saux), sMyName, .true., io)
        end if
        atomix%atoms%bias (i) = atomix%atoms%bias(i) * general%bias
        if (atomix%atoms%isscf(i)) atomix%atoms%nscf = atomix%atoms%nscf + 1
        if (atomix%atoms%ismoving(i)) atomix%atoms%nmoving = atomix%atoms%nmoving + 1
      end do
    else
      call error (sMyName, "AtomsData block is missing!", .true., io)
    end if
!
    if (atomix%atoms%nscf /= 0) then
      k = 0
      allocate (atomix%atoms%scf(1:atomix%atoms%nscf))
      do i = 1, atomix%atoms%natoms
        if (atomix%atoms%isscf(i)) then
          k = k + 1
          atomix%atoms%scf (k) = atomix%atoms%id(i)
        end if
      end do
    end if
!
    if (atomix%atoms%nmoving /= 0) then
      k = 0
      allocate (atomix%atoms%moving(1:atomix%atoms%nmoving))
      do i = 1, atomix%atoms%natoms
        if (atomix%atoms%ismoving(i)) then
          k = k + 1
          atomix%atoms%moving (k) = atomix%atoms%id(i)
        end if
      end do
    end if
    call ComputeEuclideanMatrix (atomix%atoms, io)
!
    if (general%ReadVelocity) then
      if (GetBlock(io, "VelocitiesData", nt)) then
        errno = 0
        k = 0
        allocate (ids(1:atomix%atoms%natoms))
        ids = 0
        do
          read (nt, fmt=*, iostat=errno) atomId, vx, vy, vz
          if ((errno /= 0) .and. (errno /= IOSTAT_END)) then
            call error ("block VelocitiesData is not in the right format", sMyName, .true., io)
          end if
          if (errno == IOSTAT_END) then
            exit
          end if
          write (io%udeb, '(i6,1x,f0.8,1x,f0.8,1x,f0.8)') atomId, vx, vy, vz
          k = k + 1
          if (k > atomix%atoms%natoms) then
            write (saux, '(a,i0,a)') "block VelocitiesData contains more than ", atomix%atoms%natoms, " entries"
            call error (trim(saux), sMyName, .true., io)
          end if
          if (isInList(atomId, ids)) then
            write (saux,*) "block VelocitiesData contains duplicated atom ", atomId
            call error (trim(saux), sMyName, .true., io)
          end if
          ids (k) = atomId
          if ( .not. isInList(atomId, atomix%atoms%id)) then
            write (saux,*) "block VelocitiesData contains undefined atom: ", atomix%atoms%id(i)
            call error (trim(saux), sMyName, .true., io)
          end if
          atomix%atoms%vx (atomId) = vx
          atomix%atoms%vy (atomId) = vy
          atomix%atoms%vz (atomId) = vz
        end do
        deallocate (ids)
      else
        call error ("Block VelocitiesData is missing even if ReadVel=T!", sMyName, .true., io)
      end if
    end if
!
!! donor acceptor spacer atoms
    atomix%atoms%nacceptor = GetInteger (io, "NAcceptor", 0)
    atomix%atoms%ndonor = GetInteger (io, "Ndonor", 0)
!
    if (atomix%atoms%ndonor /= 0) then
      allocate (atomix%atoms%donor(1:atomix%atoms%ndonor))
      atomix%atoms%donor = 0
      if (GetBlock(io, 'DonorAtoms', nt)) then
        read (nt,*, iostat=errno) (atomix%atoms%donor(i), i=1, atomix%atoms%ndonor)
        if (errno /= 0) then
          call error ("Unexpexted end of DonorAtoms block", sMyName, .true., io)
        end if
        do i = 1, atomix%atoms%ndonor
          if ( .not. isInList(atomix%atoms%donor(i), atomix%atoms%id)) then
            write (saux, '(a,i0)') "Invalid atom specified for donor", atomix%atoms%donor(i)
            call error (trim(saux), sMyName, .true., io)
          end if
        end do
      else
        call error ("No atoms for Donor found", sMyName, .true., io)
      end if
    end if
!
    if (atomix%atoms%nacceptor /= 0) then
      allocate (atomix%atoms%acceptor(1:atomix%atoms%nacceptor))
      atomix%atoms%acceptor = 0
      if (GetBlock(io, 'AcceptorAtoms', nt)) then
        read (nt,*, iostat=errno) (atomix%atoms%acceptor(i), i=1, atomix%atoms%nacceptor)
        if (errno /= 0) then
          call error ("Unexpexted end of AcceptorAtoms block", sMyName, .true., io)
        end if
        do i = 1, atomix%atoms%nacceptor
          if ( .not. isInList(atomix%atoms%acceptor(i), atomix%atoms%id)) then
            write (saux, '(a,i0)') "Invalid atom specified for acceptor", atomix%atoms%acceptor(i)
            call error (trim(saux), sMyName, .true., io)
          end if
        end do
      else
        call error ("No atoms for Acceptor found", sMyName, .true., io)
      end if
    end if
    atomix%atoms%nspacer = atomix%atoms%natoms - atomix%atoms%ndonor - atomix%atoms%nacceptor
    if (atomix%atoms%nspacer /= 0) then
      allocate (atomix%atoms%spacer(1:atomix%atoms%nspacer))
      atomix%atoms%spacer = 0
      k = 0
      do i = 1, atomix%atoms%natoms
        if (( .not. isInList(atomix%atoms%id(i), atomix%atoms%donor)) .and. ( .not. isInList(atomix%atoms%id(i), &
       & atomix%atoms%acceptor))) then
          k = k + 1
          atomix%atoms%spacer (k) = atomix%atoms%id(i)
        end if
      end do
    end if
    atomix%atoms%ncurrent = GetInteger (io, "NCurrent", 0)
!
    if (atomix%atoms%ncurrent /= 0) then
      allocate (atomix%atoms%current(1:atomix%atoms%ncurrent))
      if (atomix%atoms%ncurrent == atomix%atoms%natoms) then
        atomix%atoms%current = atomix%atoms%id
      else
        atomix%atoms%current = 0
        if (GetBlock(io, 'CurrentAtoms', nt)) then
          read (nt,*, iostat=errno) (atomix%atoms%current(i), i=1, atomix%atoms%ncurrent)
          if (errno /= 0) then
            call error ("Unexpected end of CurrentAtoms block", sMyName, .true., io)
          end if
          do i = 1, atomix%atoms%ncurrent
            if ( .not. isInList(atomix%atoms%current(i), atomix%atoms%id)) then
              write (saux, '(a,i0)') "Invalid atom specified for current calculations", atomix%atoms%current(i)
              call error (trim(saux), sMyName, .true., io)
            end if
          end do
        else
          call error ("No atoms for CurrentAtoms found", sMyName, .true., io)
        end if
      end if
    end if
    atomix%atoms%ncurrentOnBonds = GetInteger (io, "NCurrentOnBonds", 0)
    if (atomix%atoms%ncurrent > 0) then
      allocate (atomix%atoms%currentonBonds(1:2*atomix%atoms%ncurrentOnBonds))
      atomix%atoms%currentonBonds = 0
      if (GetBlock(io, 'CurrentOnBonds', nt)) then
        read (nt,*, iostat=errno) (atomix%atoms%currentonBonds(2*i-1), atomix%atoms%currentonBonds(2*i), i=1, &
       & atomix%atoms%ncurrentOnBonds)
        if (errno /= 0) then
          call error ("Unexpected end of CurrentOnBonds block", sMyName, .true., io)
        end if
      else
        call error ("No atoms for CurrentAtoms found", sMyName, .true., io)
      end if
    end if
  end subroutine ReadAtoms
!
!> \page atoms Atoms data  
!> The strucure of a typical block will be  
!> \verbatim  
!> NumberOfAtoms 1  
!> block AtomsData  
!> 1 0.0 0.0 0.0 0.0 T F  
!> endblock AtomsData  
!> \endverbatim  
!> NumberOfAtoms specifies the number of Atoms to be read from the block AtomsData\n  
!> each atom has a line with 7 entities\n  
!> <1> is the id of the specie of the atom (integer)\n  
!> <2-4> carteasian coordinates (reals)\n  
!> <5> local bias (real) gets multiplied internal with the values speciefied by \c Bias \n  
!> <6> should the atom be treated scf in a scf calculation (logical)\n  
!> <7> should the atom move in a molecular dynamics simulation (logical)\n  
!> The velocities block  
!> \verbatim  
!> block VelocitiesData  
!> 1 0.0 0.0 3.0  
!> endblock VelocitiesData  
!> \endverbatim  
!> VelocitiesData block allows to specify initial velocties for atoms. Not all the atoms have to be present.  
!> The missing ones have the velocities Initialized with k_zero. Each atom has a line\n  
!> <1> the id of the atom (integer) \n  
!> <2-4> the carteasian components of the velocity (reals) \n  
!> AcceptorAtoms and DonorAtoms blocks specify which atoms belong to the acceptor and which to the donor.  
!> The rest of the atoms are considered to be part of the spacer. Each block waits for a line with a  
!> list of atoms in between the \c block \c endblock. \c Ndonor and \c NAcceptor specify how many atoms to be  
!> expected in each list  
!>\verbatim  
!> NDonor 1  
!> block DonorAtoms  
!> 3  
!> endblock DonorAtoms  
!>  
!> NAcceptor 2  
!> block AcceptorAtoms  
!> 2 1  
!> endblock AcceptorAtoms  
!>\endverbatim  
!> - if you read the VelocitiesData block not all the exceptions are caught.  
!> check the output to be sure that it read the proper stuff. (fortran will automatically convert an integer to real  
!> so will generate a valid entry when is not the case)  
!> - atoms ids are set according to their position in the AtomsData block  
!
!> \brief reads and allocates the info about species of atoms
!> \details this is all that should be done in this module.
!> the rest should happen in a specialized module
!> \author Alin M Elena
!> \date 29/10/07, 17:49:04
!> \param io type(ioType) i/o units
!> \param general type(generalType) general data
!> \param specs type(speciesType) contains info about atoms
!> \param specBasis type(orbitalType) optional contains the information about the basis set for each specie
  subroutine ReadSpecies (io, general, specs, specBasis)
    character (len=*), parameter :: sMyName = "ReadSpecies"
    type (ioType), intent (inout) :: io
    type (generalType), intent (in) :: general
    type (speciesType), intent (inout) :: specs
    type (orbitalType), intent (inout), optional :: specBasis (:, :)
!internal variables    
    integer :: nt, i, j, errno, sp, maxL, shift
    character (len=k_ml) :: saux
!
    specs%nspecies = GetInteger (io, "NumberOfSpecies",-1)
    if (specs%nspecies <= 0) then
      call error (sMyName, "NumberOfSpecies token is missing or is negative! ---", .true., io)
    end if
    specs%created = .true.
    if ( .not. general%scf) then
      allocate (specs%id(1:specs%nspecies))
      allocate (specs%mass(1:specs%nspecies))
      allocate (specs%z(1:specs%nspecies))
      allocate (specs%zval(1:specs%nspecies))
      allocate (specs%ulocal(1:specs%nspecies, 1:1))
      allocate (specs%jlocal(1:specs%nspecies, 1:1))
      allocate (specs%uinter(1:specs%nspecies))
      allocate (specs%norbs(1:specs%nspecies))
      specs%norbs = 0
      specs%zval = 0
      if (GetBlock(io, "SpeciesData", nt)) then
        do i = 1, specs%nspecies
          read (nt, fmt=*, iostat=errno) sp, specs%z(i), specs%mass(i), specs%ulocal(i, 1), specs%jlocal(i, 1), specs%uinter(i)
          specs%id (i) = i
          if (sp /= i) then
            write (saux, '(a,i0,a,i0,a,i0)') "the id of the specie differs from the one expected by convention found ", sp, " expec&
           &ted ", i, " will make it ", i
            call error (trim(saux), sMyName, .true., io)
          end if
          if (errno /= 0) then
            call error ("block SpeciesData is not in the right format", sMyName, .true., io)
          end if
          write (io%udeb, '(2i6,4f16.6)') specs%id(i), specs%z(i), specs%mass(i), specs%ulocal(i, 1), specs%jlocal(i, 1), &
         & specs%uinter(i)
          specs%mass (i) = specs%mass(i) * k_amuToInternal
        end do
        if (isUnique(specs%id(1:specs%nspecies))) then
          call error ("id must be unique for each specie", sMyName, .true., io)
        end if
      else
        call error ("Specie block is missing!!!", sMyName, .true., io)
      end if
    else
      select case (general%scfType)
      case (k_scfTbuj)
        if ( .not. general%spin) then
          call error (sMyName, "This model suppose to be spin polarised try SCFType = TB+U ---", .true., io)
        end if
        allocate (specs%id(1:specs%nspecies))
        allocate (specs%mass(1:specs%nspecies))
        allocate (specs%z(1:specs%nspecies))
        allocate (specs%zval(1:specs%nspecies))
        allocate (specs%ulocal(1:specs%nspecies, 1:1))
        allocate (specs%jlocal(1:specs%nspecies, 1:1))
        allocate (specs%uinter(1:specs%nspecies))
        allocate (specs%norbs(1:specs%nspecies))
        specs%norbs = 0
        specs%zval = 0
        if (GetBlock(io, "SpeciesData", nt)) then
          do i = 1, specs%nspecies
            read (nt, fmt=*, iostat=errno) sp, specs%z(i), specs%mass(i), specs%ulocal(i, 1), specs%jlocal(i, 1), specs%uinter(i)
            specs%id (i) = i
            if (sp /= i) then
              write (saux, '(a,i0,a,i0,a,i0)') "the id of the specie differs from the one expected by convention found ", sp, " exp&
             &ected ", i, " will make it ", i
              call error (trim(saux), sMyName, .true., io)
            end if
            if (errno /= 0) then
              call error ("block SpeciesData is not in the right format", sMyName, .true., io)
            end if
            write (io%udeb, '(2i6,4f16.6)') specs%id(i), specs%z(i), specs%mass(i), specs%ulocal(i, 1), specs%jlocal(i, 1), &
           & specs%uinter(i)
            specs%mass (i) = specs%mass(i) * k_amuToInternal
          end do
          if (isUnique(specs%id(1:specs%nspecies))) then
            call error ("id must be unique for each specie", sMyName, .true., io)
          end if
        else
          call error ("Specie block is missing!!!", sMyName, .true., io)
        end if
      case (k_scfTbu)
        if (general%spin) then
          call error (sMyName, "This model suppose to be spinless try SCFType = TB+UJ ---", .true., io)
        end if
        allocate (specs%id(1:specs%nspecies))
        allocate (specs%mass(1:specs%nspecies))
        allocate (specs%z(1:specs%nspecies))
        allocate (specs%zval(1:specs%nspecies))
        allocate (specs%ulocal(1:specs%nspecies, 1:1))
        allocate (specs%jlocal(1:specs%nspecies, 1:1))
        allocate (specs%uinter(1:specs%nspecies))
        allocate (specs%norbs(1:specs%nspecies))
        specs%norbs = 0
        specs%zval = 0
        specs%jlocal = 0.0_k_pr
        if (GetBlock(io, "SpeciesData", nt)) then
          do i = 1, specs%nspecies
            read (nt, fmt=*, iostat=errno) sp, specs%z(i), specs%mass(i), specs%ulocal(i, 1), specs%uinter(i)
            specs%id (i) = i
            if (sp /= i) then
              write (saux, '(a,i0,a,i0,a,i0)') "the id of the specie differs from the one expected by convention found ", sp, " exp&
             &ected ", i, " will make it ", i
              call error (trim(saux), sMyName, .true., io)
            end if
            if (errno /= 0) then
              call error ("block SpeciesData is not in the right format", sMyName, .true., io)
            end if
            write (io%udeb, '(2i6,3f16.6)') specs%id(i), specs%z(i), specs%mass(i), specs%ulocal(i, 1), specs%uinter(i)
            specs%mass (i) = specs%mass(i) * k_amuToInternal
          end do
          if (isUnique(specs%id(1:specs%nspecies))) then
            call error ("id must be unique for each specie", sMyName, .true., io)
          end if
        else
          call error ("Specie block is missing!!!", sMyName, .true., io)
        end if
      case (k_scftbuo, k_scftbujo)
        if ( .not. present(specBasis)) then
          allocate (specs%id(1:specs%nspecies))
          allocate (specs%mass(1:specs%nspecies))
          allocate (specs%z(1:specs%nspecies))
          allocate (specs%zval(1:specs%nspecies))
          allocate (specs%uinter(1:specs%nspecies))
          allocate (specs%norbs(1:specs%nspecies))
          allocate (specs%ulocal(1:specs%nspecies, 1:1))
          allocate (specs%jlocal(1:specs%nspecies, 1:1))
          specs%jlocal = 0.0_k_pr
          specs%ulocal = 0.0_k_pr
          specs%norbs = 0
          specs%zval = 0
!
          if (GetBlock(io, "SpeciesData", nt)) then
            do i = 1, specs%nspecies
              read (nt, fmt=*, iostat=errno) sp, specs%z(i), specs%mass(i), specs%uinter(i)
              specs%id (i) = i
              if (sp /= i) then
                write (saux, '(a,i0,a,i0,a,i0)') "the id of the specie differs from the one expected by convention found ", sp, " e&
               &xpected ", i, " will make it ", i
                call error (trim(saux), sMyName, .true., io)
              end if
              if (errno /= 0) then
                call error ("block SpeciesData is not in the right format", sMyName, .true., io)
              end if
              write (io%udeb, '(2i6,2f16.6)') specs%id(i), specs%z(i), specs%mass(i), specs%uinter(i)
              specs%mass (i) = specs%mass(i) * k_amuToInternal
            end do
            if (isUnique(specs%id(1:specs%nspecies))) then
              call error ("id must be unique for each specie", sMyName, .true., io)
            end if
          else
            call error ("Specie block is missing!!!", sMyName, .true., io)
          end if
        else
          maxL = Lmax (specBasis, specs) + 1
          deallocate (specs%ulocal)
          deallocate (specs%jlocal)
          if (general%scfType == k_scftbuo) then
            allocate (specs%ulocal(1:specs%nspecies, 0:maxL))
            allocate (specs%jlocal(1:specs%nspecies, 0:maxL))
            specs%jlocal = 0.0_k_pr
            specs%ulocal = 0.0_k_pr
            if (general%spin) then
              call error (sMyName, "This model suppose to be spinless try SCFType = TB+UJ ---", .true., io)
            end if
            if (GetBlock(io, "HubbardU", nt)) then
              do i = 1, specs%nspecies
                read (nt, fmt=*, iostat=errno) sp
                if (sp /= i) then
                  write (saux, '(a,i0,a,i0)') "the id of the specie differs from the one expected by convention found ", sp, " expe&
                 &cted ", i
                  call error (trim(saux), sMyName, .true., io)
                end if
                if (errno /= 0) then
                  call error ("block HubbardU is not in the right format", sMyName, .true., io)
                end if
                write (io%udeb, '(i0)') sp
                specs%ulocal (i, 0) = getLMax (i, specBasis, specs) + 1
                do j = 1, ceiling (specs%ulocal(i, 0))
                  read (nt, fmt=*, iostat=errno) specs%ulocal(i, j)
                  if (errno /= 0) then
                    call error ("block SpeciesData is not in the right format", sMyName, .true., io)
                  end if
                  write (io%udeb, '(f16.6)') specs%ulocal(i, j)
                end do
              end do
            else
              call error ("HubbardU block is missing!!!", sMyName, .true., io)
            end if
          end if
          if (general%scfType == k_scftbujo) then
            allocate (specs%ulocal(1:specs%nspecies, 0:2*maxL))
            allocate (specs%jlocal(1:specs%nspecies, 0:2*maxL))
            specs%jlocal = 0.0_k_pr
            specs%ulocal = 0.0_k_pr
            if ( .not. general%spin) then
              call error (sMyName, "This model suppose to be spin polarised try SCFType = TB+UJO ---", .true., io)
            end if
            if (GetBlock(io, "HubbardUJ", nt)) then
              do i = 1, specs%nspecies
                read (nt, fmt=*, iostat=errno) sp
                if (sp /= i) then
                  write (saux, '(a,i0,a,i0)') "the id of the specie differs from the one expected by convention found ", sp, " expe&
                 &cted ", i
                  call error (trim(saux), sMyName, .true., io)
                end if
                if (errno /= 0) then
                  call error ("block HubbardUJ is not in the right format", sMyName, .true., io)
                end if
                write (io%udeb, '(i0)') sp
                specs%ulocal (i, 0) = getLMax (i, specBasis, specs) + 1
                specs%jlocal (i, 0) = getLMax (i, specBasis, specs) + 1
                shift = specs%ulocal (i, 0)
                do j = 1, ceiling (specs%ulocal(i, 0))
                  read (nt, fmt=*, iostat=errno) specs%ulocal(i, j), specs%ulocal(i, j+shift), specs%jlocal(i, j), specs%jlocal(i, &
                 & j+shift)
                  if (errno /= 0) then
                    call error ("block SpeciesData is not in the right format", sMyName, .true., io)
                  end if
                  write (io%udeb, '(4f16.6)') specs%ulocal(i, j), specs%ulocal(i, j+shift), specs%jlocal(i, j), specs%jlocal(i, &
                 & j+shift)
                end do
              end do
            else
              call error ("HubbardUJ block is missing!!!", sMyName, .true., io)
            end if
          end if
        end if
      end select
    end if
  end subroutine ReadSpecies
!
!> \page species "Atomic Species Data"  
!> The structure of the SpeciesData block  
!> \verbatim  
!> NumberOfSpecies 2  
!> block SpeciesData  
!> 1 2 0.0 1.0 2.0 3.0  
!> 3 3 4.0 2.0 3.0 4.0  
!> endblock SpeciesData  
!> \endverbatim  
!> It expects to read a number NumberofSpecies lines in beetween the \c block \c endblock  
!> The structue of a line is\n  
!> <1> id of the specie (integer). It has to be unique. \n  
!> <2> Z the atomic number (integer) \n  
!> <3> atomic mass (in a.m.u.) it gets converted internally to the chosen system of units \n  
!> <4-5> U J (reals) U,J for the SCF calculation \n  
!> <6> the screening factor for tha electrostatic interaction if the option is chosen. \n  
!> \f[ V_I=\frac{e^2}{4\pi\epsilon_0}\cfrac{q_J}{\sqrt{r_{IJ}^2+\left ( \cfrac{1}{U_I}+\cfrac{1}{U_J}\right )^2}} \f]  
!
!
!> \brief reads and Initializes the basis set
!> \author Alin M Elena
!> \date 30/10/07, 16:43:38
!> \param io type(ioType) i/o units
!> \param gen type(generalType) general data
!> \param atomix type(atomicxType) contains info about atoms
  subroutine ReadBasis (io, gen, atomix)
    character (len=*), parameter :: sMyName = "ReadBasis"
    type (ioType), intent (inout) :: io
    type (generalType), intent (in) :: gen
    type (atomicxType), intent (inout) :: atomix
!
    integer :: i, j, nt, errno, sp, norbitals
    character (len=k_ml) :: saux
    integer :: n, l, m, k
    real (k_pr) :: qu, qd
    integer, allocatable :: tmpid (:)
!
    allocate (atomix%speciesBasis(1:atomix%species%nspecies, 1:gen%maxOrbitalsPerAtom))
    do i = 1, atomix%species%nspecies
      do j = 1, gen%maxOrbitalsPerAtom
        atomix%speciesBasis(i, j)%sp = 0
        atomix%speciesBasis(i, j)%atom = 0
        atomix%speciesBasis(i, j)%spin = .false.
        atomix%speciesBasis(i, j)%n = 0
        atomix%speciesBasis(i, j)%l = 0
        atomix%speciesBasis(i, j)%m = 0
        atomix%speciesBasis(i, j)%occup = 0.0_k_pr
      end do
    end do
    allocate (tmpid(1:atomix%species%nspecies))
    tmpid = 0
    if (GetBlock(io, "Basis", nt)) then
      do i = 1, atomix%species%nspecies
        read (nt, fmt=*, iostat=errno) sp, norbitals
        if (errno /= 0) then
          write (saux, '(a,i0)') "error reading specie no ", i
          call error (trim(saux), sMyName, .true., io)
        end if
        if ( .not. isInList(sp, atomix%species%id)) then
          write (saux, '(a,i0)') "undefined specie ", sp
          call error (trim(saux), sMyName, .true., io)
        end if
        if (isInList(sp, tmpid)) then
          write (saux, '(a,i0)') "specie basis already specified ", sp
          call error (trim(saux), sMyName, .true., io)
        end if
        tmpid (i) = sp
        if (gen%spin) then
          atomix%species%norbs (i) = 2 * norbitals
        else
          atomix%species%norbs (i) = norbitals
        end if
        if (atomix%species%norbs(i) > gen%maxOrbitalsPerAtom) then
          write (saux, '(a,i0)') "Maximum number of orbitals per atom exceeded for specie ", sp
          call error (trim(saux), sMyName, .true., io)
        end if
        do j = 1, norbitals
          if (gen%spin) then
            read (nt, fmt=*, iostat=errno) n, l, m, qd, qu
          else
            read (nt, fmt=*, iostat=errno) n, l, m, qd
          end if
          if (errno /= 0) then
            write (saux, '(a,i0,a,i0)') "error reading orbital ", j, " specie ", sp
            call error (trim(saux), sMyName, .true., io)
          end if
          if (gen%spin) then
! spin down
            atomix%speciesBasis(i, j)%sp = sp
            atomix%speciesBasis(i, j)%atom = 0
            atomix%speciesBasis(i, j)%spin = .false.
            atomix%speciesBasis(i, j)%n = n
            atomix%speciesBasis(i, j)%l = l
            atomix%speciesBasis(i, j)%m = m
            atomix%speciesBasis(i, j)%occup = qd
! spin up
            atomix%speciesBasis(i, j+norbitals)%sp = sp
            atomix%speciesBasis(i, j+norbitals)%atom = 0
            atomix%speciesBasis(i, j+norbitals)%spin = .true.
            atomix%speciesBasis(i, j+norbitals)%n = atomix%speciesBasis(i, j)%n
            atomix%speciesBasis(i, j+norbitals)%l = atomix%speciesBasis(i, j)%l
            atomix%speciesBasis(i, j+norbitals)%m = atomix%speciesBasis(i, j)%m
            atomix%speciesBasis(i, j+norbitals)%occup = qu
            atomix%species%zval (i) = atomix%species%zval(i) + atomix%speciesBasis(i, j+norbitals)%occup + atomix%speciesBasis(i, &
           & j)%occup
          else
            atomix%speciesBasis(i, j)%sp = sp
            atomix%speciesBasis(i, j)%atom = 0
            atomix%speciesBasis(i, j)%spin = .false.
            atomix%speciesBasis(i, j)%n = n
            atomix%speciesBasis(i, j)%l = l
            atomix%speciesBasis(i, j)%m = m
            atomix%speciesBasis(i, j)%occup = qd
            atomix%species%zval (i) = atomix%species%zval(i) + atomix%speciesBasis(i, j)%occup
          end if
        end do
      end do
!
! allocate the basis for the atoms
      k = maxval (atomix%species%norbs)
      allocate (atomix%atoms%orbs(1:atomix%atoms%natoms, 1:k))
      do i = 1, atomix%atoms%natoms
        do j = 1, k
          atomix%atoms%orbs (i, j) = 0
        end do
      end do
      atomix%basis%norbitals = 0
      do i = 1, atomix%atoms%natoms
        atomix%basis%norbitals = atomix%basis%norbitals + atomix%species%norbs(atomix%atoms%sp(i))
      end do
      allocate (atomix%basis%orbitals(1:atomix%basis%norbitals))
!     ! build the basis
      k = 0
      if (gen%spin) then
! spin down        
        do i = 1, atomix%atoms%natoms
          sp = atomix%atoms%sp(i)
          atomix%atoms%norbs (i) = atomix%species%norbs(sp)
          do j = 1, atomix%species%norbs(sp)
            if ( .not. atomix%speciesBasis(sp, j)%spin) then
              k = k + 1
              atomix%basis%orbitals (k) = atomix%speciesBasis(sp, j)
              atomix%basis%orbitals(k)%atom = i
!                 ! this array maps which orbitals belong to atom i
              atomix%atoms%orbs (i, j) = k
            end if
          end do
        end do
!spin up
        do i = 1, atomix%atoms%natoms
          sp = atomix%atoms%sp(i)
          atomix%atoms%norbs (i) = atomix%species%norbs(sp)
          do j = 1, atomix%species%norbs(sp)
            if (atomix%speciesBasis(sp, j)%spin) then
              k = k + 1
              atomix%basis%orbitals (k) = atomix%speciesBasis(sp, j)
              atomix%basis%orbitals(k)%atom = i
!                 ! this array maps which orbitals belong to atom i
              atomix%atoms%orbs (i, j) = k
            end if
          end do
        end do
      else
        do i = 1, atomix%atoms%natoms
          sp = atomix%atoms%sp(i)
          atomix%atoms%norbs (i) = atomix%species%norbs(sp)
          do j = 1, atomix%species%norbs(sp)
            k = k + 1
            atomix%basis%orbitals (k) = atomix%speciesBasis(sp, j)
            atomix%basis%orbitals(k)%atom = i
!                 ! this array maps which orbitals belong to atom i
            atomix%atoms%orbs (i, j) = k
          end do
        end do
      end if
    else
      call error ("Basis block not found", sMyName, .true., io)
    end if
    deallocate (tmpid)
  end subroutine ReadBasis
!
!> \page basis Basis Set
!> The \c basis block
!> \verbatim
!> block Basis
!> 1 4
!> 0 0  0 1.0 1.0
!> 0 1  1 0.2 0.2
!> 0 1 -1 0.3 0.3
!> 0 1  0 0.4 0.4
!> 2 1
!> 0 0  0 0.5 0.5
!> endblock Basis
!> \endverbatim
!> The structure for each specie is:\n
!> A line with two integers is expected. First is the id of the specie, second the number of orbitals for that specie.
!> A number of lines equal with the number of orbitals specified previously is expected.
!> the structure of each line is \n
!> <1-3> n,l,m (integers) n principal quantum number, l angular quantum number m magnetic quantum number
!> n is not used at anything for the moment
!> <4-5> occupation numbers for the orbital. In the case of a spin polarised calculation (\c Spin T) both o them are read
!> for (\c Spin F) only the first one is read. They are used to compute \c Zval
!
!> \brief selects the initialization routine for the tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:07:56
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine ReadTBModel (io, gen, atomix, tbMod)
    character (len=*), parameter :: sMyName = "ReadTBModel"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomix
    type (modelType), intent (inout) :: tbMod
!
    integer :: i, j, k
!
    allocate (tbMod%hopping(1:atomix%species%nspecies, 1:atomix%species%nspecies))
!
    do i = 1, atomix%species%nspecies
      do j = 1, atomix%species%nspecies
        tbMod%hopping(i, j)%l1 = atomix%speciesBasis(i, 1)%l
        do k = 2, atomix%species%norbs(i)
          if (tbMod%hopping(i, j)%l1 < atomix%speciesBasis(i, k)%l) then
            tbMod%hopping(i, j)%l1 = atomix%speciesBasis(i, k)%l
          end if
        end do
!
        tbMod%hopping(i, j)%l2 = atomix%speciesBasis(j, 1)%l
        do k = 2, atomix%species%norbs(j)
          if (tbMod%hopping(i, j)%l2 < atomix%speciesBasis(j, k)%l) then
            tbMod%hopping(i, j)%l2 = atomix%speciesBasis(j, k)%l
          end if
        end do
!
        tbMod%hopping(i, j)%ll = Min (tbMod%hopping(i, j)%l1, tbMod%hopping(i, j)%l2)
        allocate (tbMod%hopping(i, j)%a(0:tbMod%hopping(i, j)%l1, 0:tbMod%hopping(i, j)%l2, 0:tbMod%hopping(i, j)%ll))
        allocate (tbMod%hopping(i, j)%eps(0:atomix%species%norbs(i)-1))
!
      end do
    end do
!
    do i = 1, atomix%species%nspecies
      do j = 1, atomix%species%nspecies
        tbMod%hopping(i, j)%eps = 0.0_k_pr
        tbMod%hopping(i, j)%a1 = 0.0_k_pr
        tbMod%hopping(i, j)%a2 = 0.0_k_pr
        tbMod%hopping(i, j)%a3 = 0.0_k_pr
        tbMod%hopping(i, j)%a4 = 0.0_k_pr
        tbMod%hopping(i, j)%phi0 = 0.0_k_pr
        tbMod%hopping(i, j)%r0 = 0.0_k_pr
        tbMod%hopping(i, j)%rc = 0.0_k_pr
        tbMod%hopping(i, j)%r1 = 0.0_k_pr
        tbMod%hopping(i, j)%rcut = 0.0_k_pr
        tbMod%hopping(i, j)%n = 0.0_k_pr
        tbMod%hopping(i, j)%nc = 0.0_k_pr
        tbMod%hopping(i, j)%d0 = 0.0_k_pr
        tbMod%hopping(i, j)%dc = 0.0_k_pr
        tbMod%hopping(i, j)%d1 = 0.0_k_pr
        tbMod%hopping(i, j)%dcut = 0.0_k_pr
        tbMod%hopping(i, j)%m = 0.0_k_pr
        tbMod%hopping(i, j)%mc = 0.0_k_pr
        tbMod%hopping(i, j)%a = 0.0_k_pr
      end do
    end do
    select case (gen%bond)
    case (k_bondGSP)
      call ReadTbGSP (io, gen, atomix, tbMod)
      call PrintTbGSP (io, gen, atomix, tbMod)
    case (k_bondHarrison)
      call ReadTbHarrison (io, gen, atomix, tbMod)
      call PrintTbHarrison (io, gen, atomix, tbMod)
    end select
!
  end subroutine ReadTBModel
!
!
!> \brief selects the initialization routine for the Goodwin-Skinner-Petifor tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:15:15
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
!> \remarks Using tight binding model as described in: \n
!> A.P. Horsfield,P.D. Godwin,D.G. Pettifor,A.P. Sutton \n
!> Computational materials synthesis. I. A tight-binding schek_me for hydrocarbons - Phys. Rev. B 54, 22, 15 773
!
  subroutine ReadTbGSP (io, gen, atomix, tbMod)
    character (len=*), parameter :: sMyName = "ReadTbGSP"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomix
    type (modelType), intent (inout) :: tbMod
    integer :: i, j, i1, k, k1, k2, p
    character (len=k_mw) :: readVar
!
    do i = 1, atomix%species%nspecies
      do j = 1, atomix%species%nspecies
        if (i == j) then
          if (gen%spin) then
            i1 = 0
            p = 0
            do while (i1 <= atomix%species%norbs(i)/ 2-1)
              readVar = "eps"
              readVar = trim (ccvar(i, i, readVar))
              tbMod%hopping(i, j)%eps(i1:i1+2*p) = GetReal (io, trim(ccvar(i1, i1, readVar)), 0.0_k_pr)
              tbMod%hopping(i, j)%eps(i1+atomix%species%norbs(i)/2:i1+atomix%species%norbs(i)/2+2*p) = tbMod%hopping(i, j)%eps(i1)
              i1 = i1 + 2 * p + 1
              p = p + 1
            end do
          else
            i1 = 0
            p = 0
            do while (i1 <= atomix%species%norbs(i)- 1)
              readVar = "eps"
              readVar = trim (ccvar(i, i, readVar))
              tbMod%hopping(i, j)%eps(i1:i1+2*p) = GetReal (io, trim(ccvar(i1, i1, readVar)), 0.0_k_pr)
              i1 = i1 + 2 * p + 1
              p = p + 1
            end do
          end if
          if (gen%embedding) then
            readVar = "a1"
            tbMod%hopping(i, j)%a1 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
            readVar = "a2"
            tbMod%hopping(i, j)%a2 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
            readVar = "a3"
            tbMod%hopping(i, j)%a3 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
            readVar = "a4"
            tbMod%hopping(i, j)%a4 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
          end if
        end if
        readVar = "phi0"
        tbMod%hopping(i, j)%phi0 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "r0"
        tbMod%hopping(i, j)%r0 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "rc"
        tbMod%hopping(i, j)%rc = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "r1"
        tbMod%hopping(i, j)%r1 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "rcut"
        tbMod%hopping(i, j)%rcut = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "n"
        tbMod%hopping(i, j)%n = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "nc"
        tbMod%hopping(i, j)%nc = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "d0"
        tbMod%hopping(i, j)%d0 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "dc"
        tbMod%hopping(i, j)%dc = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "d1"
        tbMod%hopping(i, j)%d1 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "dcut"
        tbMod%hopping(i, j)%dcut = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "m"
        tbMod%hopping(i, j)%m = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "mc"
        tbMod%hopping(i, j)%mc = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        do k = 0, tbMod%hopping(i, j)%l1
          do k1 = k, tbMod%hopping(i, j)%l2
            do k2 = 0, k
              tbMod%hopping(i, j)%a(k, k1, k2) = GetReal (io, trim(ccnlm(i, j, k, k1, k2)), 0.0_k_pr)
            end do
          end do
        end do
      end do
    end do
!
    do i = 1, atomix%species%nspecies
      do j = 1, atomix%species%nspecies
        do k = 0, tbMod%hopping(i, j)%l1
          do k1 = 0, tbMod%hopping(i, j)%l2
            do k2 = 0, Min (k, k1)
              if (k > k1) then
                tbMod%hopping(i, j)%a(k, k1, k2) = (-1.0_k_pr) ** ((k-k1+Abs(k-k1))/2.0_k_pr) * tbMod%hopping(j, i)%a(k1, k, k2)
              end if
            end do
          end do
        end do
      end do
    end do
!
  end subroutine ReadTbGSP
!
!> \brief selects the initialization routine for the Harrison tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:16:15
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine ReadTbHarrison (io, gen, atomix, tbMod)
    character (len=*), parameter :: sMyName = "ReadTbHarrison"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomix
    type (modelType), intent (inout) :: tbMod
!
    integer i, j, i1, k, k1, k2, p
    character (len=k_mw) :: readVar
!
    do i = 1, atomix%species%nspecies
      do j = 1, atomix%species%nspecies
        if (i == j) then
          if (gen%spin) then
            i1 = 0
            p = 0
            do while (i1 <= atomix%species%norbs(i)/ 2-1)
              readVar = "eps"
              readVar = trim (ccvar(i, i, readVar))
              tbMod%hopping(i, j)%eps(i1:i1+2*p) = GetReal (io, trim(ccvar(i1, i1, readVar)), 0.0_k_pr)
              tbMod%hopping(i, j)%eps(i1+atomix%species%norbs(i)/2:i1+atomix%species%norbs(i)/2+2*p) = tbMod%hopping(i, j)%eps(i1)
              i1 = i1 + 2 * p + 1
              p = p + 1
            end do
          else
            i1 = 0
            p = 0
            do while (i1 <= atomix%species%norbs(i)/ 2-1)
              readVar = "eps"
              readVar = trim (ccvar(i, i, readVar))
              tbMod%hopping(i, j)%eps(i1:i1+2*p) = GetReal (io, trim(ccvar(i1, i1, readVar)), 0.0_k_pr)
              i1 = i1 + 2 * p + 1
              p = p + 1
            end do
          end if
          if (gen%embedding) then
            readVar = "a1"
            tbMod%hopping(i, j)%a1 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
            readVar = "a2"
            tbMod%hopping(i, j)%a2 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
            readVar = "a3"
            tbMod%hopping(i, j)%a3 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
            readVar = "a4"
            tbMod%hopping(i, j)%a4 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
          end if
        end if
!
        readVar = "phi0"
        tbMod%hopping(i, j)%phi0 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "r1"
        tbMod%hopping(i, j)%r1 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "rcut"
        tbMod%hopping(i, j)%rcut = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "n"
        tbMod%hopping(i, j)%n = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "d1"
        tbMod%hopping(i, j)%d1 = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "dcut"
        tbMod%hopping(i, j)%dcut = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
        readVar = "m"
        tbMod%hopping(i, j)%m = GetReal (io, trim(ccvar(i, j, readVar)), 0.0_k_pr)
!
        do k = 0, tbMod%hopping(i, j)%l1
          do k1 = k, tbMod%hopping(i, j)%l2
            do k2 = 0, k
              tbMod%hopping(i, j)%a(k, k1, k2) = GetReal (io, trim(ccnlm(i, j, k, k1, k2)), 0.0_k_pr)
            end do
          end do
        end do
      end do
    end do
!
!
    do i = 1, atomix%species%nspecies
      do j = 1, atomix%species%nspecies
        do k = 0, tbMod%hopping(i, j)%l1
          do k1 = 0, tbMod%hopping(i, j)%l2
            do k2 = 0, Min (k, k1)
              if (k > k1) then
                tbMod%hopping(i, j)%a(k, k1, k2) = (-1.0_k_pr) ** ((k-k1+Abs(k-k1))/2.0_k_pr) * tbMod%hopping(j, i)%a(k1, k, k2)
              end if
            end do
          end do
        end do
      end do
    end do
!
  end subroutine ReadTbHarrison
!
!> \brief reads the delta block
!> \details the block is read only if multipoles electrostatics is asked for
!> \author Alin M Elena
!> \date 31/10/07, 17:52:20
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine ReadDelta (io, gen, atomix, tbMod)
    character (len=*), parameter :: sMyName = "ReadDelta"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomix
    type (modelType), intent (inout) :: tbMod
    integer :: nt, sp, j, i, k, l, errno
    character (len=k_ml) :: saux
    integer, allocatable :: tmpid (:)
!
    if (GetBlock(io, "DeltaPol", nt)) then
      allocate (tbMod%delta(1:atomix%species%nspecies))
      allocate (tmpid(1:atomix%species%nspecies))
      tmpid (1:atomix%species%nspecies) = 0
      do i = 1, atomix%species%nspecies
        read (nt, fmt=*, iostat=errno) sp
        if (errno /= 0) then
          call error ("Specie indicator missing, please fix it!(reading delta)", sMyName, .true., io)
        end if
        if ( .not. isInList(sp, atomix%species%id)) then
          write (saux, "(a,i0)") "undefined specie: ", sp
          call error (trim(saux), sMyName, .true., io)
        end if
        if (isInList(sp, tmpid)) then
          write (saux, "(a,i0)") "specie has been already read: ", sp
          call error (trim(saux), sMyName, .true., io)
        end if
        tmpid (i) = sp
        tbMod%delta(i)%sp = sp
        if (gen%spin) then
          j = maxval (atomix%speciesBasis(tbMod%delta(i)%sp, 1:atomix%species%norbs(tbMod%delta(i)%sp)/2)%l)
        else
          j = maxval (atomix%speciesBasis(tbMod%delta(i)%sp, 1:atomix%species%norbs(tbMod%delta(i)%sp))%l)
        end if
        tbMod%delta(i)%l = j
        allocate (tbMod%delta(i)%d(0:2*j, 0:j, 0:j))
        tbMod%delta(i)%d = 0.0_k_pr
        do j = 0, tbMod%delta(i)%l
          do k = j, tbMod%delta(i)%l
            do l = Abs (j-k), j + k
              if (Mod(k+j+l, 2) == 0) then
                if ((k == j) .and. (l == 0)) then
                  tbMod%delta(i)%d(0, k, j) = 1.0_k_pr
                else
                  read (nt,*, iostat=errno) tbMod%delta(i)%d(l, j, k)
                  if (errno /= 0) then
                    write (saux, '(a,i0,a,3i0)') "Error reading delta for specie ", sp, " indices: ", l, j, k
                    call error (trim(saux), sMyName, .true., io)
                  end if
                  tbMod%delta(i)%d(l, k, j) = tbMod%delta(i)%d(l, j, k)
                end if
              end if
            end do
          end do
        end do
      end do
      deallocate (tmpid)
    else
      call error ("DeltaPol block is missing", sMyName, .true., io)
    end if
  end subroutine ReadDelta
!
!> \page delta DeltaPol Block  
!> Note that: delta (l,l1,l2) with |l1-l2|<=l<=|l1+l2| \n  
!>  we read delta(l,l1,l2) \n  
!> the rules are that we read only l1<l2 if l1=l2 then we read only if l/=0 \n  
!> for each specie we have \n  
!>  we start with a line specifying the specie then first line is delta(1,0,1) second delta(2,1,1) (for a sp specie) \n  
!> if one specie has no element to be read it is mandatory to be present  
!> \verbatim  
!>      block DeltaPol  
!>      1  
!>      1.0  
!>      0.0  
!>      2  
!>      endblock DeltaPol  
!> \endverbatim  
!
!> \brief deallocates the memory
!> \author Alin M Elena
!> \date 29/10/07, 17:38:16
!> \param atomic type(atomicxType) contains atomic data
!> \param general type(generalType) general data that controls the flow of the program
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param io type(ioType) contains all the info about I/O files
  subroutine CleanMemory (io, atomic, general, tbMod, Sol)
    character (len=*), parameter :: sMyName = "CleanMemory"
    type (ioType), intent (in) :: io
    type (atomicxType), intent (inout) :: atomic
    type (generalType), intent (inout) :: general
    type (modelType), intent (inout) :: tbMod
    type (solutionType), intent (inout) :: Sol
    integer :: i, j, ierr
!
    deallocate (atomic%speciesBasis)
    deallocate (atomic%species%id)
    deallocate (atomic%species%mass)
    deallocate (atomic%species%z)
    deallocate (atomic%species%zval)
    deallocate (atomic%species%ulocal)
    deallocate (atomic%species%jlocal)
    deallocate (atomic%species%uinter)
    deallocate (atomic%species%norbs)
!
    deallocate (atomic%atoms%id)
    deallocate (atomic%atoms%sp)
    deallocate (atomic%atoms%bias)
    deallocate (atomic%atoms%isscf)
    deallocate (atomic%atoms%ismoving)
    deallocate (atomic%atoms%x)
    deallocate (atomic%atoms%y)
    deallocate (atomic%atoms%z)
    deallocate (atomic%atoms%vx)
    deallocate (atomic%atoms%vy)
    deallocate (atomic%atoms%vz)
    deallocate (atomic%atoms%fx)
    deallocate (atomic%atoms%fy)
    deallocate (atomic%atoms%fz)
    deallocate (atomic%atoms%dx)
    deallocate (atomic%atoms%dy)
    deallocate (atomic%atoms%dz)
    deallocate (atomic%atoms%xo)
    deallocate (atomic%atoms%yo)
    deallocate (atomic%atoms%zo)
    deallocate (atomic%atoms%vxo)
    deallocate (atomic%atoms%vyo)
    deallocate (atomic%atoms%vzo)
    deallocate (atomic%atoms%fxo)
    deallocate (atomic%atoms%fyo)
    deallocate (atomic%atoms%fzo)
    deallocate (atomic%atoms%chrg)
    deallocate (atomic%atoms%chrg0)
    deallocate (atomic%atoms%norbs)
    if (allocated(atomic%atoms%spacer)) deallocate (atomic%atoms%spacer)
    if (allocated(atomic%atoms%donor)) deallocate (atomic%atoms%donor)
    if (allocated(atomic%atoms%current)) deallocate (atomic%atoms%current)
    if (allocated(atomic%atoms%acceptor)) deallocate (atomic%atoms%acceptor)
    if (allocated(atomic%atoms%scf)) deallocate (atomic%atoms%scf)
    if (allocated(atomic%atoms%moving)) deallocate (atomic%atoms%moving)
    deallocate (atomic%basis%orbitals)
    deallocate (atomic%atoms%orbs)
    deallocate (atomic%atoms%MagMom)
    if (allocated(atomic%atoms%currentonBonds)) deallocate (atomic%atoms%currentonBonds)
    do i = 1, atomic%atoms%natoms
      if (atomic%atoms%neighbours(i)%created) then
        deallocate (atomic%atoms%neighbours(i)%a)
      end if
    end do
    deallocate (atomic%atoms%neighbours)
    do i = 1, atomic%species%nspecies
      do j = 1, atomic%species%nspecies
        deallocate (tbMod%hopping(i, j)%a)
        deallocate (tbMod%hopping(i, j)%eps)
        deallocate (tbMod%hopping(i, j)%hmnTail)
      end do
    end do
    deallocate (tbMod%hopping)
    if (general%electrostatics == k_electrostaticsMultipoles) then
      do i = 1, atomic%species%nspecies
        deallocate (tbMod%delta(i)%d)
      end do
      deallocate (tbMod%delta)
    end if
!
    call DestroyMatrix (Sol%h, io)
    call DestroyMatrix (Sol%forceOp, io)
    call DestroyMatrix (Sol%eigenvecs, io)
    deallocate (Sol%eigenvals)
!
    if (general%scf) then
      call DestroyMatrix (Sol%hin, io)
      call DestroyMatrix (Sol%h2, io)
      deallocate (Sol%potential)
      deallocate (Sol%field)
    end if
    if (general%spin) then
      call DestroyMatrix (Sol%hdown, io)
      call DestroyMatrix (Sol%hup, io)
    end if
    call DestroyMatrix (Sol%rho, io)
    if ((general%scf) .and. (general%electrostatics == k_electrostaticsMultipoles)) then
      deallocate (Sol%gcoeff)
      deallocate (Sol%rgc)
    end if
!
    deallocate (Sol%n0)
    deallocate (Sol%density)
    deallocate (Sol%buff%dins)
    deallocate (Sol%buff%douts)
    deallocate (Sol%buff%res)
    deallocate (Sol%buff%densityin)
    deallocate (Sol%buff%densityout)
    deallocate (Sol%buff%densitynext)
    deallocate (Sol%buff%tmpA)
    call DestroyMatrix (Sol%buff%tmpB, io)
    deallocate (Sol%buff%pos1)
    deallocate (Sol%buff%pos2)
    deallocate (Sol%buff%f)
    deallocate (Sol%buff%g)
    deallocate (Sol%buff%a)
    deallocate (Sol%buff%nstart)
    deallocate (Sol%fact)
    deallocate (Sol%Distances)
    deallocate (Sol%buff%itmp)
    deallocate (Sol%CurrentMatrix)
    deallocate (Sol%CurrentMatrix2)
    call DestroyMatrix (Sol%buff%h, io)
    if (general%smearMethod == k_smMP) then
      deallocate (Sol%hermite)
    end if
    if ((general%runType == k_runEhrenfestDamped) .or. (general%runType == k_runEhrenfest)) then
      call DestroyMatrix (Sol%rhodot, io)
      call DestroyMatrix (Sol%deltaRho, io)
      call DestroyMatrix (Sol%rhoold, io)
      call DestroyMatrix (Sol%rhonew, io)
      call DestroyMatrix (Sol%rho0, io)
    end if
    if (( .not. general%compElec) .and. (general%electrostatics == k_electrostaticsMultipoles)) then
      do i = 1, atomic%atoms%natoms
        deallocate (Sol%delq(i)%a)
        deallocate (Sol%vs(i)%a)
      end do
    end if
    call MKL_FreeBuffers ()
  end subroutine CleanMemory
end module m_ReadData
