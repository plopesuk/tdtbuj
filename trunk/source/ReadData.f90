!> \brief deals with reading the input file(s)
!> \details it saves the gathered info in the right variables
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \remarks

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

contains


!> \brief reads from ioLoc different fields of genLoc
!> \author Alin M Elena
!> \date ~2007
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to run
!> \param atomic type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters

  subroutine Initialize(ioLoc,genLoc,atomic,tbMod)
    character(len=*), parameter :: myname = 'Initialize'
    type(ioType), intent(inout) :: ioLoc
    type(generalType), intent(inout) :: genLoc
    type(atomicxType), intent(inout) :: atomic
    type(modelType), intent(inout) :: tbMod
    real(k_pr) :: aux
    integer :: errno=-1,nt,i

    ioLoc%uerr=GetUnit()
    ioLoc%inpErr=trim(ioLoc%inpFile)//".err"
    open(unit=ioLoc%uerr,file=trim(ioLoc%inpErr),status="replace",&
      action="write",iostat=errno)
    if (errno /= 0) then
      write(*,*)"I can not create error file ", trim(ioLoc%inpFile)//".err"
      stop
    endif 
    ioLoc%uinp=GetUnit()
    open(unit=ioLoc%uinp,file=ioLoc%inpFile,status='old',&
      action='read',iostat=errno)

    if (errno /= 0) then
      write(ioLoc%uerr,*)"I can not open file ", ioLoc%inpFile
      stop
    endif

! parse file and in the same time 
! makes a minimal check of corectness
    call cpu_time(genLoc%time%int)

    call ParseFile(ioLoc)
    call cpu_time(genLoc%time%end)
    aux=genLoc%time%end-genLoc%time%int

    call cpu_time(genLoc%time%int)

    call ReadIo(ioLoc)

    if (ioLoc%debug>=k_mediumDebug) then
      call cpu_time(genLoc%time%end)
      write(ioLoc%udeb,'(a,f16.6,a)')"Parsing time "&
        ,aux," seconds"
      write(ioLoc%udeb,'(a,f16.6,a)')"Setting ioInfo "&
        ,genLoc%time%end-genLoc%time%int," seconds"
    endif

    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,1)
    call ReadGeneral(ioLoc,genLoc)
    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,2)
    call InitializeConstants(genLoc%units)
    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,1)
    call ReadSpecies(ioLoc,genLoc,atomic%species)
    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,2)

    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,1)
    call ReadAtoms(ioLoc,genLoc,atomic)
    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,2)

    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,1)
    call ReadBasis(ioLoc,genLoc,atomic)
    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,2)

    call PrintSpecies(ioLoc,atomic%species)
    call PrintAtoms(ioLoc,genLoc,atomic)
    call PrintBasis(ioLoc,genLoc,atomic)

    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,1)
    call ReadTBModel(ioLoc,genLoc,atomic,tbMod)
    if (genLoc%electrostatics==k_electrostaticsMultipoles) then
      call ReadDelta(ioLoc,genLoc,atomic,tbMod)
      call PrintDelta(ioLoc,genLoc,atomic,tbMod)
    endif
    if (ioLoc%debug>=k_mediumDebug) call timing(ioLoc,genLoc%time,2)

    call EndParse
  end subroutine Initialize

!> \brief reads the names for I/O files and the output levels
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param ioLoc type(ioType)
!> \details
!> A short description of the input file ioInfo part
!> \latexonly
!>
!> \begin{tabular}{|c||p{0.2\textwidth}|c||p{0.25\textwidth}|}
!> \hline
!> \textbf{VarName} & \textbf{Values} & \textbf{Default Value} & \textbf{Description} \\\
!> \hline \hline
!> DebugFile & string & input file name.dbg & Name of the debug file \\ 
!> \hline
!> OutputFile & string & input file name.out & Name of the output file \\ 
!> \hline
!> OnScreen & logical & .false. & on/off Printing on the screen, on means that nothing will be written in the output file \\ 
!> \hline
!> DebugLevel & 5,15,25 & 5 & debug information level (low, medium,high) \\
!> \hline
!> OutputLevel & 5,15,25 & 5 & output information level (low, medium,high) \\
!> \hline
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

subroutine ReadIo(ioLoc)
character(len=*),parameter :: myname="ReadIo"
  type(ioType), intent(inout) :: ioLoc
  integer :: errno


!comm_io DebugFile & string & input file name.dbg & Name of the debug file \\

! read the name of output file plus debug (y/n) and output level
  ioLoc%debugFile=GetString(ioLoc,"DebugFile",&
    trim(ioLoc%inpFile)//".dbg",.false.)
  ioLoc%udeb=GetUnit()
  open(unit=ioLoc%udeb,file=trim(ioLoc%debugFile),status="replace",&
    action="write",iostat=errno)
  if (errno /= 0) call error("I can not create file "//&
    trim(ioLoc%debugFile),myname,.true.,ioLoc)

!comm_io OutputFile & string & input file name.out & Name of the output file \\
!comm_io OnScreen & logical & .false. & on/off Printing on the screen, on means that nothing will be written in the output file \\
  ioLoc%outputFile=GetString(ioLoc,"OutputFile",&
    trim(ioLoc%inpFile)//".out")
  ioLoc%stdout=GetLogical(ioLoc,"OnScreen",.false.)
! prepares the output according to the settings from input file(s)
  if (ioLoc%stdout) then
    ioLoc%uout=6
  else
    ioLoc%uout=GetUnit()
    open(unit=ioLoc%uout,file=trim(ioLoc%outputFile),status="replace",&
      action="write",iostat=errno)
    if (errno /= 0) call error("I can not create file "//trim(ioLoc%outputFile)&
      ,myname,.true.,ioLoc)
  endif

! levels
  ioLoc%verbosity=GetInteger(ioLoc,"OutputLevel",5)
  ioLoc%debug=GetInteger(ioLoc,"DebugLevel",5)
  ioLoc%firstTime=.not. ioLoc%firstTime

  if (ioLoc%verbosity>=k_highDebug) then
    write( ioLoc%udeb,'(a,l1)')"Read ioInfo for the first time: "&
      ,ioLoc%firstTime
  endif

  write( ioLoc%uout,'(a,a)')"Input file: "&
    ,trim(ioLoc%inpFile)
  write( ioLoc%uout,'(a,a)')"Output file(OutputFile): "&
    ,trim(ioLoc%outputFile)
  write( ioLoc%uout,'(a,i0)')"Output Level(OutputLevel): "&
    ,ioLoc%verbosity
  write( ioLoc%uout,'(a,a)')"Debug file(DebugFile): "&
    ,trim(ioLoc%debugFile)
  write( ioLoc%uout,'(a,i0)')"Debug Level(DebugLevel): "&
    ,ioLoc%debug
  write( ioLoc%uout,'(a,a)')"Error file: "&
    ,trim(ioLoc%inpErr)
  write( ioLoc%uout,'(a,l1)')"Standard Output(OnScreen): "&
    ,ioLoc%stdout


end subroutine ReadIo


!> \brief reads the parameters necessary for running
!> \author Alin M Elena
!> \date 21st of January, 2007 
!> \param  ioLoc type(ioType) I/O details (see m_Types::ioType)
!> \param genLoc type(generalType), keeps all the general information about the parameters of the program
!> \details
!> A short description of the input file general variables part
!> \latexonly
!> \begin{longtable}{|c||p{0.2\textwidth}|c||p{0.25\textwidth}|}
!> \hline
!> \textbf{VarName} & \textbf{Values} & \textbf{Default Value} & \textbf{Description} \\\
!> \hline \hline
!> JobName & string & no name & a name for the job \\ 
!> \hline
!> RanSeed & real(0,1)& 0.5 & A seed for the random number generator \\
!> \hline
!> ReadVel & logical & .false. & on/off reading velocity block \\
!> \hline
!> SCF & logical & .false. & on/off self consistent field method \\
!> \hline
!> \end{longtable}
!> \endlatexonly
!> \htmlonly
!> <TABLE CELLPADDING=3 BORDER="1">
!> <TR><TH ALIGN="CENTER"><B>VarName</B></TH>
!> <TH ALIGN="LEFT" VALIGN="TOP" WIDTH=100><B>Values</B></TH>
!> <TH ALIGN="CENTER"><B>Default Value</B></TH>
!> <TH ALIGN="LEFT" VALIGN="TOP" WIDTH=125><B>Description</B></TH>
!> </TR>
!> <TR><TD ALIGN="CENTER">JobName</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>string</TD>
!> <TD ALIGN="CENTER">no name</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>a name for the job</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">RanSeed</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real (0,1) </TD>
!> <TD ALIGN="CENTER">0.5</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>A seed for the random number generator</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">ReadVel</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off reading velocity block</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Spin</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off spin polarisation</TD>
!> </TR>
!> </TABLE>
!> \endhtmlonly

subroutine ReadGeneral(ioLoc,genLoc)
  character(len=*),parameter :: name="ReadGeneral"
  type(ioType), intent(inout) :: ioLoc
  type(generalType), intent(inout) :: genLoc

  character(len=k_mw) :: saux

    if (.not. genLoc%firstTime) genLoc%firstTime=.not. genLoc%firstTime

    genLoc%jobName=GetString(ioLoc,"JobName","no name")
    write( ioLoc%uout,'(a,a)')"Job Name(JobName): "&
      ,genLoc%jobName

    genLoc%ranseed=GetReal(ioLoc,"RanSeed",0.5_k_pr)
    write( ioLoc%uout,'(a,f0.16)')"Random number generator seed(RanSeed): "&
      ,genLoc%ranseed

    genLoc%ReadVelocity=GetLogical(ioLoc,"ReadVel",.false.)
    write( ioLoc%uout,'(a,l1)')"Read Velocity(ReadVel): "&
     ,genLoc%ReadVelocity

!comm_gen VBias & real & 0.0 & bias factor\\
    genLoc%bias=GetReal(ioLoc,"VBias",0.0_k_pr)
    write( ioLoc%uout,'(a,f0.8)')"Bias(VBias): "&
     ,genLoc%bias


!comm_gen MaxOrbsPerAtom & integer & 8 & maximum number of orbitals per atom\\
    genLoc%maxOrbitalsPerAtom=GetInteger(ioLoc,"MaxOrbsPerAtom",8)
    write( ioLoc%uout,'(a,i0)')"Maximum Orbitals Per Atom (MaxOrbsPerAtom): "&
      ,genLoc%maxOrbitalsPerAtom

    genLoc%spin=GetLogical(ioLoc,"Spin",.false.)
    write( ioLoc%uout,'(a,l1)')"spin polarisation(Spin): "&
      ,genLoc%spin

!comm_gen Units & AU EV SI & AU & system of units: atomic units(AU), electronVolt-Angstrom(eVA), International(SI)\\
    saux=GetString(ioLoc,"Units","AU")
    write( ioLoc%uout,'(a,a)')"System of Units (Units): "&
      ,trim(saux)
    if (cstr(trim(saux),'AU')) then
      genLoc%units = k_unitsAU
    elseif (cstr(trim(saux),'eVA')) then
      genLoc%units = k_unitsEV
    elseif (cstr(trim(saux),'SI')) then
      genLoc%units = k_unitsSI
    elseif (cstr(trim(saux),'RY')) then
      genLoc%units = k_unitsRY
    else
      call error("The requested system of units is not implemented",name,.true.,ioLoc)
    endif

 !comm_gen BondType & Harrison, GSP &Harrison& bond type \\
! \hline
    saux=GetString(ioLoc,"BondType","Harrison")
    write( ioLoc%uout,'(a,a)')"bond type (BondType): "&
      ,trim(saux)
    if (cstr(trim(saux),'Harrison')) then
      genLoc%bond = k_bondHarrison
    elseif (cstr(trim(saux),'GSP')) then
      genLoc%bond = k_bondGSP
    else
      call error("The requested bond type is not implemented",name,.true.,ioLoc)
    endif

   !comm_gen Embedding & logical & .true. & on/off embedding method\\
    genLoc%embedding=GetLogical(ioLoc,"Embedding",.true.)
    write( ioLoc%uout,'(a,l1)')"has embedding(Embedding): "&
     ,genLoc%embedding

!comm_gen Electrostatics & PointCharges, Multipoles & PointCharges & method use to compute electrostatic interaction\\
    saux=GetString(ioLoc,"Electrostatics","PointCharges")
    write( ioLoc%uout,'(a,a)')"Electrostatics type (Electrostatics): "&
     ,trim(saux)
    if (cstr(trim(saux),'PointCharges')) then
      genLoc%electrostatics = k_electrostaticsPoint
    elseif (cstr(trim(saux),'Multipoles')) then
      genLoc%electrostatics = k_electrostaticsMultipoles
    else
      call error("The requested Electrostatics type is not implemented",name,.true.,ioLoc)
    endif

!comm_gen PrecomputeMultipoles & logical & .true. & on/off precompute multipoles\\
    genLoc%compElec=GetLogical(ioLoc,"PrecomputeMultipoles",.true.)
    write( ioLoc%uout,'(a,l1)')"precompute multipoles(PrecomputeMultipoles): "&
      ,genLoc%compElec

    genLoc%scf=GetLogical(ioLoc,"SCF",.false.)
    write( ioLoc%uout,'(a,l1)')"Is SCF?(SCF): "&
      ,genLoc%scf
 !comm_gen SCFType & TB+UJ & TB+UJ & SCF method\\
    saux=GetString(ioLoc,"SCFType","TB+UJ")
     write( ioLoc%uout,'(a,a)')"SCF type(SCFType): "&
       ,trim(saux)
    if (cstr(trim(saux),'TB+UJ')) then
       genLoc%scfType = k_scfTbuj
    else
       call error("The requested SCFType is not implemented",name,.true.,ioLoc)
    endif

 !comm_gen SCFSteps & integer & 100 & maximum number of steps used for scf\\
    genLoc%maxscf=GetInteger(ioLoc,"SCFSteps",100)
    write( ioLoc%uout,'(a,i0)')"SCF Maximum number of steps(SCFSteps): "&
      ,genLoc%maxscf
  !comm_gen SCFMix & real & 0.85 & mixing parameter\\
    genLoc%scfMix=GetReal(ioLoc,"SCFMix",0.85_k_pr)
    write( ioLoc%uout,'(a,f0.8)')"SCF Mixing parameter(SCFMix): "&
      ,genLoc%scfMix

!comm_gen SCFTol & real & 1e-8 & convergence tolerance\\
    genLoc%scftol=GetReal(ioLoc,"SCFTol",1.0e-8_k_pr)
    write( ioLoc%uout,'(a,ES12.4)')"SCF convergence tolerance(SCFTol): "&
      ,genLoc%scftol

  !comm_gen SCFMixN & integer & 4 & number of iterations to mix\\
    genLoc%scfMixn=GetInteger(ioLoc,"SCFMixN",4)
    write( ioLoc%uout,'(a,i0)')"number of iterations to mix(SCFMixN): "&
      ,genLoc%scfMixn

!comm_gen RunType & SinglePoint, BODynamics, Ehrenfest, EhrenfestDamped, Fit, ForceTest, ForceTestX,
! ForceTestY, ForceTestZ & SinglePoint & Type of calculation\\
    saux=GetString(ioLoc,"RunType","SinglePoint")
    write( ioLoc%uout,'(a,a)')"Calculation type(RunType): "&
      ,trim(saux)
    if (cstr(trim(saux),'SinglePoint')) then
      genLoc%runType = k_runSp
    elseif (cstr(trim(saux),'bodynamics')) then
      genLoc%runType = k_runBO
    elseif (cstr(trim(saux),'ehrenfest')) then
      genLoc%runType = k_runEhrenfest
    elseif (cstr(trim(saux),'fit')) then
      genLoc%runType = k_runFit
    elseif (cstr(trim(saux),'forcetest')) then
      genLoc%runType = k_runForceTest
    elseif (cstr(trim(saux),'forcetestx')) then
      genLoc%runType = k_runForceTestx
    elseif (cstr(trim(saux),'forcetesty')) then
      genLoc%runType = k_runForceTesty
    elseif (cstr(trim(saux),'forcetestz')) then
      genLoc%runType = k_runForceTestz
    elseif (cstr(trim(saux),'ehrenfestdumped')) then
      genLoc%runType = k_runEhrenfestDumped
    elseif (cstr(trim(saux),'GeometryOptimization')) then
      genLoc%runType = k_runGeomBFGS
    elseif (cstr(trim(saux),'Special')) then
      genLoc%runType = k_runSpecial
    else
      call error("The requested RunType is not implemented",name,.true.,ioLoc)
    endif

! comm_gen HElThres & real & 1e-10 & hamiltionian element thresold. Any element smaller that the thresold is made zero.\\
    genLoc%hElementThreshold=GetReal(ioLoc,"HElTres",1.0e-10_k_pr)
    write( ioLoc%uout,'(a,ES12.4)')"hamiltionian element thresold(HElThres): "&
      ,genLoc%hElementThreshold

!comm_gen CollinearSpins & logical & .false. & on/off collinear spins\\
      genLoc%collinear=GetLogical(ioLoc,"CollinearSpins",.false.)
      write( ioLoc%uout,'(a,l1)')"collinear spins(CollinearSpins): "&
      ,genLoc%collinear
!comm_gen MaxIt & integer & 500 & maximum number of iterations used to find Fermi level\\
      genLoc%maxIt=GetInteger(ioLoc,"MaxIt",500)
      write( ioLoc%uout,'(a,i0)')"Maximum number of iterations to find Fermi level(MaxIt): "&
      ,genLoc%maxit
! comm_gen ChargeTol & real & 1e-10 & charge tolerance used to find Fermi level\\
    genLoc%qTolerance=GetReal(ioLoc,"ChargeTol",1.0e-10_k_pr)
    write( ioLoc%uout,'(a,ES12.4)')"Charge Tolerance used to find Fermi level(ChargeTol): "&
      ,genLoc%qTolerance
!comm_gen NetCharge & real & 0.0 & Net charge on the system\\
      genLoc%netcharge=GetReal(ioLoc,"NetCharge",0.0_k_pr)
      write( ioLoc%uout,'(a,ES12.4)')"Net charge(NetCharge): "&
      ,genLoc%netcharge
    saux=GetString(ioLoc,"SmearingMethod","FD")
    write( ioLoc%uout,'(a,a)')"Smearing Methos used to find Fermi level (SmearingMethod): "&
      ,trim(saux)
    if (cstr(trim(saux),'FD')) then
      genLoc%smearMethod = k_smFD
    elseif (cstr(trim(saux),'MP')) then
      genLoc%smearMethod = k_smMP
    elseif (cstr(trim(saux),'CMU')) then
      genLoc%smearMethod = k_smCMU
    else
      call error("The requested SmearingMethod is not implemented",name,.true.,ioLoc)
    endif

    select case(genLoc%smearMethod)
      case(k_smFD)
!comm_gen ElectronicTemperature & real & 300.0 & Electronic temperature, used to compute occupation numbers if you choose Fermi-Dirac method\\
        genLoc%electronicTemperature=GetReal(ioLoc,"ElectronicTemperature",300.0_k_pr)
        write( ioLoc%uout,'(a,ES12.4)')"Electronic Temperature(ElectronicTemperature): "&
          ,genLoc%electronicTemperature
      case(k_smMP)
!comm_gen ElectronicW & real & 0.05 & Electronic W, used to compute occupation numbers if you choose Methfessel-Paxton method\\
        genLoc%MPW=GetReal(ioLoc,"ElectronicW",0.05_k_pr)
        write( ioLoc%uout,'(a,ES12.4)')"Electronic W(ElectronicW): "&
          ,genLoc%MPW
!comm_gen MPN & integer & 2 & the order of Hermite polynomials used to find Fermi level by Methfessel-Paxton method \\
      genLoc%mpN=GetInteger(ioLoc,"MPN",500)
      write( ioLoc%uout,'(a,i0)')"Order of Germite polynomials used find Fermi level(MPN): "&
      ,genLoc%maxit
      case(k_smCMU)
!comm_gen ElectronicMu & real & 0.0 & chemical potential, used to compute occupation numbers if you choose constant $\mu$ method\\
        genLoc%electronicMu=GetReal(ioLoc,"ElectronicMu",0.0_k_pr)
        write( ioLoc%uout,'(a,ES12.4)')"Chemical Potential(ElectronicMu): "&
          ,genLoc%electronicMu
    end select
! !comm_gen DMOccTol & real & 1e-10 & density matrix occupation tolerance\\
   genLoc%dmOccupationTolerance=GetReal(ioLoc,"DMOccTol",1.0e-10_k_pr)
   write( ioLoc%uout,'(a,ES12.4)')"density matrix occupation tolerance(DMOccTol): "&
     ,genLoc%dmOccupationTolerance!
!comm_gen SymRefRho & logical & .false. & on/off symmetric reference density matrix\\
  genLoc%SymRefRho=GetLogical(ioLoc,"SymRefRho",.false.)
  write( ioLoc%uout,'(a,l1)')"Symetric Reference Density Matrix (SymRefRho): "&
   ,genLoc%SymRefRho
!comm_gen Screened & logical & .false. & on/off bare Coulomb or screened electrostatic interaction \\
  genLoc%screened=GetLogical(ioLoc,"Screened",.false.)
  write( ioLoc%uout,'(a,l1)')"Screened Coulomb interation (Screened): "&
    ,genLoc%screened
!comm_gen WriteAnimation & logical & .true. & on/off writing animation file \\
  genLoc%writeAnimation=GetLogical(ioLoc,"WriteAnimation",.true.)
  write( ioLoc%uout,'(a,l1)')"Write animation(WriteAnimation): "&
    ,genLoc%writeAnimation
!comm_gen DeltaT & real & 0.001 & time step used to evolve equtions of motion\\
  genLoc%deltat=GetReal(ioLoc,"DeltaT",0.001_k_pr)
  write( ioLoc%uout,'(a,ES12.4)')"Time step(DeltaT): "&
    ,genLoc%deltat
 
!comm_gen Nsteps & integer & 100 & number of steps used to evolve equtions of motion\\
  genLoc%nsteps=GetInteger(ioLoc,"Nsteps",100)
  write( ioLoc%uout,'(a,i0)')"Number of steps(Nsteps): "&
    ,genLoc%nsteps

!comm_gen VelScale & logical & .false. & on/off scaling velocities\\
  genLoc%scaleVelocities=GetLogical(ioLoc,"VelScale",.false.)
  write( ioLoc%uout,'(a,l1)')"Scale  Velocities?(VelocitiesScale): "&
     ,genLoc%scaleVelocities

!comm_gen IonicTemperature & real & 300.0 & Ionic temperature\\
  genLoc%ionicTemperature=GetReal(ioLoc,"IonicTemperature",300.0_k_pr)
  write( ioLoc%uout,'(a,g)')"Ionic Temperature(IonicTemperature): "&
     ,genLoc%ionicTemperature
 !comm_gen BiasRampSteps & integer & 0 & for how many steps to apply bias\\
  genLoc%BiasRampSteps=GetInteger(ioLoc,"BiasRampSteps",100)
  write( ioLoc%uout,'(a,i0)')"No of Bias Ramp Steps (BiasRampSteps): "&
    ,genLoc%BiasRampSteps

 !comm_gen EulerSteps & integer & 100 & after each EulerSteps apply an Euler integration of equations of motions\\
   genLoc%eulerSteps=GetInteger(ioLoc,"EulerSteps",100)
   write( ioLoc%uout,'(a,i0)')"Euler steps (EulerSteps): "&
     ,genLoc%EulerSteps
 !comm_gen Gamma & real & 0.3 & dumping factor for Ehrenfest equation\\
   genLoc%gamma=GetReal(ioLoc,"Gamma",0.5_k_pr)
   write( ioLoc%uout,'(a,g)')"damping factor in Ehrenfest Equation(Gamma): "&
     ,genLoc%gamma

! !comm_gen WriteEne & logical & .true. & on/off writing energy file \\
!   genLoc%write_ene=GetLogical(ioLoc,"WriteEne",.true.)
!   write( ioLoc%uout,'(a,l)')"Write Energy(WriteEne): "&
!     ,genLoc%write_ene
! 
! !comm_gen WriteDen & logical & .false. & on/off writing density file \\
!   genLoc%write_density=GetLogical(ioLoc,"WriteDen",.false.)
!   write( ioLoc%uout,'(a,l)')"Write density(WriteDen): "&
!     ,genLoc%write_density
! !comm_gen ReadDen & logical & .false. & on/off reading density file \\
!   genLoc%read_density=GetLogical(ioLoc,"ReadDen",.false.)
!   write( ioLoc%uout,'(a,l)')"Read Density(ReadDen): "&
!     ,genLoc%read_density
! 
!

! !comm_gen SpinDU & real & 0.0 & spin down spin up difference\\
!   genLoc%sdu=GetReal(ioLoc,"SpinDU",0.0_k_pr)
!   write( ioLoc%uout,'(a,g)')"spin down spin up difference(SpinDU): "&
!     ,genLoc%sdu
! 
! !comm_gen FSteps & integer & 100 & Number of steps used to test the force/energy consitency\\
!   genLoc%f_steps=GetInteger(ioLoc,"FSteps",100)
!   write( ioLoc%uout,'(a,i0)')"force steps (FSteps): "&
!     ,genLoc%f_steps
! 
! !comm_gen FStart & real & 0.0 & position at which the force/energy consistency test starts\\
!   genLoc%f_start=GetReal(ioLoc,"FStart",0.0_k_pr)
!   write( ioLoc%uout,'(a,g)')"initial position for force(FStart): "&
!     ,genLoc%f_start
! 
! !comm_gen Fdx & real & 0.001 & space step for the force/energy consistency test starts\\
!   genLoc%f_dx=GetReal(ioLoc,"FStart",0.001_k_pr)
!   write( ioLoc%uout,'(a,g)')"space step for force(Fdx): "&
!     ,genLoc%f_dx
! 
! !comm_gen Gamma & real & 0.3 & dumping factor for Ehrenfest equation\\
!   genLoc%gamma=GetReal(ioLoc,"Gamma",0.5_k_pr)
!   write( ioLoc%uout,'(a,g)')"damping factor in Ehrenfest Equation(Gamma): "&
!     ,genLoc%gamma
! 
 
! !comm_gen Hole & integer & 0 & no of level to create hole\\
!   genLoc%hole=GetInteger(ioLoc,"Hole",0)
!   write( ioLoc%uout,'(a,i0)')"Hole level(Hole): "&
!     ,genLoc%hole
! 
! !comm_gen Excite & integer & 0 & no of level to create excitation\\
!   genLoc%excite=GetInteger(ioLoc,"Excite",0)
!   write( ioLoc%uout,'(a,i0)')"Excitation level(Excite): "&
!     ,genLoc%excite
! !comm_gen HoleSpin & D,U & D & spin of the hole\\
!   saux=GetString(ioLoc,"HoleSpin","D")
!   write( ioLoc%uout,'(a,a)')"Spin of the Hole (HoleSpin): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'D')) then
!     genLoc%hole_spin = spin_down
!   elseif (cstr(trim(saux),'U')) then
!     genLoc%hole_spin = spin_up
!   else  
!     call error("The requested spin type is not implemented",name,.true.,ioLoc)
!   endif
! 
! 
! !comm_gen ExciteSpin & D,U & D & spin of the excitation\\
!   saux=GetString(ioLoc,"ExciteSpin","D")
!   write( ioLoc%uout,'(a,a)')"Spin of the excitation (ExciteSpin): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'D')) then
!     genLoc%excite_spin = spin_down
!   elseif (cstr(trim(saux),'U')) then
!     genLoc%excite_spin = spin_up
!   else  
!     call error("The requested spin type is not implemented",name,.true.,ioLoc)
!   endif
! 
! 

end subroutine ReadGeneral


!> \brief closes all the I/O units
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param ioLoc type(ioType)

subroutine CloseIoGeneral(ioLoc)
  character(len=*), parameter :: myname = 'CloseIoGeneral'
  type(ioType), intent(inout) :: ioLoc
  integer :: errno

  if (.not.ioLoc%stdout) then
    close(unit=ioLoc%uout,iostat=errno)
    if (errno /= 0) call error("I can not close file "//trim(ioLoc%outputFile)&
      ,myname,.false.,ioLoc)
  endif
  close(unit=ioLoc%udeb,iostat=errno)
  if (errno /= 0) call error("I can not close file "//&
    trim(ioLoc%debugFile),myname,.false.,ioLoc)
  close(unit=ioLoc%uerr,iostat=errno)
  if (errno /= 0) call error("I can not close file "//&
    trim(ioLoc%inpErr),myname,.false.,ioLoc)
end subroutine CloseIoGeneral


!> \brief times different events
!> \details it should be called in pairs with stage = 1 to Initialize 
!> and stage = 2 to write down the info
!> \author Alin M Elena
!> \date 29th of October 2007
!> \param io type(ioType) i/o units
!> \param times type(timeType) time keepers
!> \param stage selects the stage

  subroutine timing(io,times,stage)
    character(len=*),parameter :: myname = "timing"
    type(ioType), intent(in) :: io
    type(timeType), intent(inout) :: times
    integer, intent(in) :: stage

    if (stage==1) then
      call cpu_time(times%end)
    elseif (stage==2) then
      call cpu_time(times%end)
      write(io%udeb,'(a,f16.6,a)')"parsing and setting "&
        ,times%end-times%int," seconds"
    endif

  end subroutine timing


!> \brief reads and allocates the info about atoms
!> \details this is all that should be done in this module.
!> the rest should happen in a specialized module
!> \author Alin M Elena
!> \date 29/10/07, 17:33:50
!> \param io type(ioType) i/o units
!> \param general type(generalType) general data
!> \param atomix type(atomicxType) contains info about atoms
!> \remarks
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
!> - atoms ids are set according to theri position in the AtomsData block

  subroutine ReadAtoms(io,general,atomix)
    character(len=*), parameter :: sMyName="ReadAtoms"
    type(ioType), intent(inout) :: io
    type(generalType), intent(in)  :: general
    type(atomicxType), intent(inout) :: atomix
    ! arguments of the subroutine
    integer :: nt, errno,i,k, atom_id
    character(len=k_ml) :: saux
    integer, allocatable :: ids(:)
    real(k_pr) :: vx,vy,vz

    atomix%atoms%natoms = GetInteger(io,"NumberOfAtoms",-1)
    if (atomix%atoms%natoms<=0) then
      call error(sMyName,"NumberOfAtoms token is missing or is negative!",.true.,io)
    endif
    atomix%atoms%created=.true.
    allocate(atomix%atoms%id(1:atomix%atoms%natoms))
    allocate(atomix%atoms%sp(1:atomix%atoms%natoms))
    allocate(atomix%atoms%bias(1:atomix%atoms%natoms))
    allocate(atomix%atoms%isscf(1:atomix%atoms%natoms))
    allocate(atomix%atoms%ismoving(1:atomix%atoms%natoms))
    allocate(atomix%atoms%x(1:atomix%atoms%natoms))
    allocate(atomix%atoms%y(1:atomix%atoms%natoms))
    allocate(atomix%atoms%z(1:atomix%atoms%natoms))
    allocate(atomix%atoms%vx(1:atomix%atoms%natoms))
    allocate(atomix%atoms%vy(1:atomix%atoms%natoms))
    allocate(atomix%atoms%vz(1:atomix%atoms%natoms))
    allocate(atomix%atoms%fx(1:atomix%atoms%natoms))
    allocate(atomix%atoms%fy(1:atomix%atoms%natoms))
    allocate(atomix%atoms%fz(1:atomix%atoms%natoms))
    allocate(atomix%atoms%dx(1:atomix%atoms%natoms))
    allocate(atomix%atoms%dy(1:atomix%atoms%natoms))
    allocate(atomix%atoms%dz(1:atomix%atoms%natoms))
    allocate(atomix%atoms%xo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%yo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%zo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%vxo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%vyo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%vzo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%fxo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%fyo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%fzo(1:atomix%atoms%natoms))
    allocate(atomix%atoms%chrg(1:atomix%atoms%natoms))
    allocate(atomix%atoms%chrg0(1:atomix%atoms%natoms))
    allocate(atomix%atoms%norbs(1:atomix%atoms%natoms))
    allocate(atomix%atoms%MagMom(1:atomix%atoms%natoms))

    atomix%atoms%vx=0.0_k_pr
    atomix%atoms%vy=0.0_k_pr
    atomix%atoms%vz=0.0_k_pr
    atomix%atoms%vxo=0.0_k_pr
    atomix%atoms%vyo=0.0_k_pr
    atomix%atoms%vzo=0.0_k_pr
    atomix%atoms%xo=0.0_k_pr
    atomix%atoms%yo=0.0_k_pr
    atomix%atoms%zo=0.0_k_pr
    atomix%atoms%dx=0.0_k_pr
    atomix%atoms%dy=0.0_k_pr
    atomix%atoms%dz=0.0_k_pr
    atomix%atoms%fx=0.0_k_pr
    atomix%atoms%fy=0.0_k_pr
    atomix%atoms%fz=0.0_k_pr
    atomix%atoms%fxo=0.0_k_pr
    atomix%atoms%fyo=0.0_k_pr
    atomix%atoms%fzo=0.0_k_pr
    atomix%atoms%chrg=0.0_k_pr
    atomix%atoms%chrg0=0.0_k_pr
    atomix%atoms%ismoving=.true.
    atomix%atoms%isscf=.true.
    atomix%atoms%nmoving=0
    atomix%atoms%nscf=0
    atomix%atoms%norbs=0
    atomix%atoms%MagMom=0.0_k_pr

    if (GetBlock(io,"AtomsData",nt)) then
      do i = 1, atomix%atoms%natoms
        read(nt,fmt=*,iostat=errno)&
          atomix%atoms%sp(i),atomix%atoms%x(i),atomix%atoms%y(i),atomix%atoms%z(i),&
          atomix%atoms%bias(i),atomix%atoms%isscf(i),atomix%atoms%ismoving(i)
        if (errno/=0) then
          call error("block AtomsData is not in the right format",sMyName,.true.,io)
        endif
          atomix%atoms%id(i)=i
          write(io%udeb,'(i6,4f16.6,1x,l1,1x,l1)')&
            atomix%atoms%sp(i),atomix%atoms%x(i),atomix%atoms%y(i),atomix%atoms%z(i),&
            atomix%atoms%bias(i),atomix%atoms%isscf(i),atomix%atoms%ismoving(i)
          if(.not.isInList(atomix%atoms%sp(i),atomix%species%id)) then
            write(saux,*)"block AtomsData contains undefined specie: ", atomix%atoms%sp(i)
            call error(trim(saux),sMyName,.true.,io)
          endif
          atomix%atoms%bias(i)=atomix%atoms%bias(i)*general%bias
          if (atomix%atoms%isscf(i)) atomix%atoms%nscf=atomix%atoms%nscf+1
          if (atomix%atoms%ismoving(i)) atomix%atoms%nmoving=atomix%atoms%nmoving+1
      enddo
    else
      call error(sMyName,"AtomsData block is missing!",.true.,io)
    endif

    if (atomix%atoms%nscf/=0) then
      k=0
      allocate(atomix%atoms%scf(1:atomix%atoms%nscf))
        do i=1,atomix%atoms%natoms
          if (atomix%atoms%isscf(i)) then
            k=k+1
            atomix%atoms%scf(k)=atomix%atoms%id(i)
          endif
        enddo
    endif

    if (atomix%atoms%nmoving/=0) then
      k=0
      allocate(atomix%atoms%moving(1:atomix%atoms%nmoving))
        do i=1,atomix%atoms%natoms
          if (atomix%atoms%ismoving(i)) then
            k=k+1
            atomix%atoms%moving(k)=atomix%atoms%id(i)
          endif
        enddo
    endif

    if (general%ReadVelocity) then
      if (GetBlock(io,"VelocitiesData",nt)) then
        errno=0
        k=0
        allocate(ids(1:atomix%atoms%natoms))
        ids=0
        do 
          read(nt,fmt=*,iostat=errno)&
            atom_id,vx,vy,vz
          if ((errno/=0).and.(errno/=IOSTAT_END)) then
            call error("block VelocitiesData is not in the right format",sMyName,.true.,io)
          endif
          if (errno==IOSTAT_END) then
            exit
          endif
          write(io%udeb,'(i6,1x,f0.8,1x,f0.8,1x,f0.8)')&
              atom_id,vx,vy,vz
          k=k+1
          if (k>atomix%atoms%natoms) then
            write(saux,'(a,i0,a)')"block VelocitiesData contains more than ", atomix%atoms%natoms, " entries"
            call error(trim(saux),sMyName,.true.,io)
          endif
          if(isInList(atom_id,ids)) then
            write(saux,*)"block VelocitiesData contains duplicated atom ", atom_id
            call error(trim(saux),sMyName,.true.,io)
          endif
          ids(k)=atom_id
          if(.not.isInList(atom_id,atomix%atoms%id)) then
            write(saux,*)"block VelocitiesData contains undefined atom: ", atomix%atoms%id(i)
            call error(trim(saux),sMyName,.true.,io)
          endif
          atomix%atoms%vx(atom_id)=vx
          atomix%atoms%vy(atom_id)=vy
          atomix%atoms%vz(atom_id)=vz
        enddo
        deallocate(ids)
      else
        call error("Block VelocitiesData is missing even if ReadVel=T!",sMyName,.true.,io)
      endif
    endif
    
!! donor acceptor spacer atoms
    atomix%atoms%nacceptor=GetInteger(io,"NAcceptor",0)
    atomix%atoms%ndonor=GetInteger(io,"Ndonor",0)

    if (atomix%atoms%ndonor/=0) then
      allocate(atomix%atoms%donor(1:atomix%atoms%ndonor))
      atomix%atoms%donor=0
      if (GetBlock(io,'DonorAtoms',nt)) then
        read(nt,*,iostat=errno)(atomix%atoms%donor(i),i=1,atomix%atoms%ndonor)
        if (errno/=0) then
          call error("Unexpexted end of DonorAtoms block",sMyName,.true.,io)
        endif
        do i=1,atomix%atoms%ndonor
          if (.not.isInList(atomix%atoms%donor(i),atomix%atoms%id)) then
            write(saux,'(a,i0)')"Invalid atom specified for donor",atomix%atoms%donor(i)
            call error(trim(saux),sMyName,.true.,io)
          endif
        enddo
      else
        call error("No atoms for Donor found",sMyName,.true.,io)
      endif
    endif

    if (atomix%atoms%nacceptor/=0) then
      allocate(atomix%atoms%acceptor(1:atomix%atoms%nacceptor))
      atomix%atoms%acceptor=0
      if (GetBlock(io,'AcceptorAtoms',nt)) then
        read(nt,*,iostat=errno)(atomix%atoms%acceptor(i),i=1,atomix%atoms%nacceptor)
        if (errno/=0) then
          call error("Unexpexted end of AcceptorAtoms block",sMyName,.true.,io)
        endif
        do i=1,atomix%atoms%nacceptor
          if (.not.isInList(atomix%atoms%acceptor(i),atomix%atoms%id)) then
            write(saux,'(a,i0)')"Invalid atom specified for acceptor",atomix%atoms%acceptor(i)
            call error(trim(saux),sMyName,.true.,io)
          endif
        enddo
      else
        call error("No atoms for Acceptor found",sMyName,.true.,io)
      endif
    endif
    atomix%atoms%nspacer=atomix%atoms%natoms-atomix%atoms%ndonor-atomix%atoms%nacceptor
    if (atomix%atoms%nspacer/=0) then
      allocate(atomix%atoms%spacer(1:atomix%atoms%nspacer))
      atomix%atoms%spacer=0
      k=0
      do i=1,atomix%atoms%natoms
        if((.not.isInList(atomix%atoms%id(i),atomix%atoms%donor)).and. &
            (.not.isInList(atomix%atoms%id(i),atomix%atoms%acceptor))) then
          k=k+1
          atomix%atoms%spacer(k)=atomix%atoms%id(i)
        endif
      enddo
    endif
  end subroutine ReadAtoms



!> \brief reads and allocates the info about species of atoms 
!> \details this is all that should be done in this module.
!> the rest should happen in a specialized module
!> \author Alin M Elena
!> \date 29/10/07, 17:49:04
!> \param io type(ioType) i/o units
!> \param general type(generalType) general data
!> \param specs type(speciesType) contains info about atoms
!> \remarks
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

  subroutine ReadSpecies(io,general,specs)
    character(len=*), parameter :: sMyName="ReadSpecies"
    type(ioType), intent(inout) :: io
    type(generalType), intent(in)  :: general
    type(speciesType), intent(inout) :: specs
    !internal variables
    integer :: nt,i,errno,sp
    character(len=k_ml) :: saux

    specs%nspecies=GetInteger(io,"NumberOfSpecies",-1)
    if (specs%nspecies <= 0) then
      call error(sMyName,"NumberOfSpecies token is missing or is negative! ---",.true.,io)
    endif
    specs%created = .true.

    allocate(specs%id(1:specs%nspecies))
    allocate(specs%mass(1:specs%nspecies))
    allocate(specs%z(1:specs%nspecies))
    allocate(specs%zval(1:specs%nspecies))
    allocate(specs%ulocal(1:specs%nspecies))
    allocate(specs%jlocal(1:specs%nspecies))
    allocate(specs%uinter(1:specs%nspecies))
    allocate(specs%norbs(1:specs%nspecies))
    specs%norbs=0
    specs%zval=0
    if (GetBlock(io,"SpeciesData",nt)) then
      do i = 1, specs%nspecies
        read(nt,fmt=*,iostat=errno) sp,specs%z(i),specs%mass(i),specs%ulocal(i),specs%jlocal(i),specs%uinter(i)
        specs%id(i)=i
        if (sp/=i) then
          write(saux,'(a,i0,a,i0,a,i0)')"the id of the specie differs from the one expected by convention found ",&
            sp," expected ",i," will make it ",i
          call error(trim(saux),sMyName,.false.,io)
        endif
        if (errno/=0) then
          call error("block SpeciesData is not in the right format",sMyName,.true.,io)
        endif
        write(io%udeb,'(2i6,4f16.6)') specs%id(i),specs%z(i),specs%mass(i),specs%ulocal(i),specs%jlocal(i),specs%uinter(i)
      enddo
      if (isUnique(specs%id(1:specs%nspecies))) then
        call error("id must be unique for each specie",sMyName,.true.,io)
      endif
    else
      call error("Specie block is missing!!!",sMyName,.true.,io)
    endif
    
  end subroutine ReadSpecies



!> \brief reads and Initializes the basis set
!> \author Alin M Elena
!> \date 30/10/07, 16:43:38
!> \param io type(ioType) i/o units
!> \param gen type(generalType) general data
!> \param atomix type(atomicxType) contains info about atoms
!> \remarks
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
  subroutine ReadBasis(io,gen,atomix)
    character(len=*), parameter :: sMyName="ReadBasis"
    type(ioType), intent(inout) :: io
    type(generalType), intent(in)  :: gen
    type(atomicxType), intent(inout) :: atomix

    integer :: i,j,nt,errno,sp, norbitals
    character(len=k_ml) :: saux
    integer :: n,l,m,k
    real(k_pr) :: qu,qd
    integer, allocatable :: tmpid(:)

    allocate(atomix%speciesBasis(1:atomix%species%nspecies,1:gen%maxOrbitalsPerAtom))
    allocate(tmpid(1:atomix%species%nspecies))
    tmpid=0
    if (GetBlock(io,"Basis",nt)) then
      do i = 1,atomix%species%nspecies
        read(nt,fmt=*,iostat=errno) sp, norbitals
        if(errno/=0) then
          write(saux,'(a,i0)')"error reading specie no ",i
          call error(trim(saux),sMyName,.true.,io)
        endif
        if(.not.isInList(sp,atomix%species%id)) then
          write(saux,'(a,i0)')"undefined specie ",sp
          call error(trim(saux),sMyName,.true.,io)
        endif
        if(isInList(sp,tmpid)) then
          write(saux,'(a,i0)')"specie basis already specified ",sp
          call error(trim(saux),sMyName,.true.,io)
        endif
        tmpid(i)=sp
        if (gen%spin) then
          atomix%species%norbs(i)=2*norbitals
        else
          atomix%species%norbs(i)=norbitals
        endif
        if (atomix%species%norbs(i) > gen%maxOrbitalsPerAtom) then
          write(saux,'(a,i0)')"Maximum number of orbitals per atom exceeded for specie ",sp
          call error(trim(saux),sMyName,.true.,io)
        endif
        do j=1,norbitals
          read(nt,fmt=*,iostat=errno) n,l,m,qd,qu
          if (errno/=0) then
            write(saux,'(a,i0,a,i0)')"error reading orbital ",j, " specie ",sp
            call error(trim(saux),sMyName,.true.,io)
          endif
          if (gen%spin) then
! spin down
            atomix%speciesBasis(i,j)%sp = sp
            atomix%speciesBasis(i,j)%atom = 0
            atomix%speciesBasis(i,j)%spin=.false.
            atomix%speciesBasis(i,j)%n=n
            atomix%speciesBasis(i,j)%l=l
            atomix%speciesBasis(i,j)%m=m
            atomix%speciesBasis(i,j)%occup=qd
! spin up   
            atomix%speciesBasis(i,j+norbitals)%sp = sp
            atomix%speciesBasis(i,j+norbitals)%atom = 0
            atomix%speciesBasis(i,j+norbitals)%spin=.true.
            atomix%speciesBasis(i,j+norbitals)%n=atomix%speciesBasis(i,j)%n
            atomix%speciesBasis(i,j+norbitals)%l=atomix%speciesBasis(i,j)%l
            atomix%speciesBasis(i,j+norbitals)%m=atomix%speciesBasis(i,j)%m
            atomix%speciesBasis(i,j+norbitals)%occup=qu
            atomix%species%zval(i)=atomix%species%zval(i)+atomix%speciesBasis(i,j+norbitals)%occup&
                +atomix%speciesBasis(i,j)%occup
          else
            atomix%speciesBasis(i,j)%sp = sp
            atomix%speciesBasis(i,j)%atom = 0
            atomix%speciesBasis(i,j)%spin=.false.
            atomix%speciesBasis(i,j)%n=n
            atomix%speciesBasis(i,j)%l=l
            atomix%speciesBasis(i,j)%m=m
            atomix%speciesBasis(i,j)%occup=qd
            atomix%species%zval(i)=atomix%species%zval(i)+atomix%speciesBasis(i,j)%occup
          endif
        enddo
      enddo

! allocate the basis for the atoms
    k=maxval(atomix%species%norbs)
    allocate(atomix%atoms%orbs(1:atomix%atoms%natoms,1:k))
    atomix%basis%norbitals = 0
    do i=1,atomix%atoms%natoms
      atomix%basis%norbitals = atomix%basis%norbitals + atomix%species%norbs(atomix%atoms%sp(i))
    enddo
    allocate(atomix%basis%orbitals(1:atomix%basis%norbitals))
!     ! build the basis
     k = 0
     if (gen%spin) then
        ! spin down
        do i=1,atomix%atoms%natoms
          sp=atomix%atoms%sp(i)
          atomix%atoms%norbs(i)=atomix%species%norbs(sp)
          do j=1,atomix%species%norbs(sp)
            if (.not.atomix%speciesBasis(sp,j)%spin) then
              k = k + 1
              atomix%basis%orbitals(k) = atomix%speciesBasis(sp,j)
              atomix%basis%orbitals(k)%atom = i
!                 ! this array maps which orbitals belong to atom i
              atomix%atoms%orbs(i,j) = k
            endif
          enddo
        enddo
!spin up
        do i=1,atomix%atoms%natoms
          sp=atomix%atoms%sp(i)
          atomix%atoms%norbs(i)=atomix%species%norbs(sp)
          do j=1,atomix%species%norbs(sp)
            if (atomix%speciesBasis(sp,j)%spin) then
              k = k + 1
              atomix%basis%orbitals(k) = atomix%speciesBasis(sp,j)
              atomix%basis%orbitals(k)%atom = i
!                 ! this array maps which orbitals belong to atom i
              atomix%atoms%orbs(i,j) = k
            endif
          enddo
        enddo
      else
        do i=1,atomix%atoms%natoms
          sp=atomix%atoms%sp(i)
          atomix%atoms%norbs(i)=atomix%species%norbs(sp)
          do j=1,atomix%species%norbs(sp)
              k = k + 1
              atomix%basis%orbitals(k) = atomix%speciesBasis(sp,j)
              atomix%basis%orbitals(k)%atom = i
!                 ! this array maps which orbitals belong to atom i
              atomix%atoms%orbs(i,j) = k
          enddo
        enddo
      endif
    else
      call error("Basis block not found",sMyName,.true.,io)
    endif
    deallocate(tmpid)
end subroutine ReadBasis



!> \brief selects the initialization routine for the tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:07:56
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine ReadTBModel(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="ReadTBModel"
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomix
    type(modelType), intent(inout) :: tbMod

    integer :: i,j,k,k1,k2

    allocate(tbMod%hopping(1:atomix%species%nspecies,1:atomix%species%nspecies))

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        tbMod%hopping(i,j)%l1=atomix%speciesBasis(i,1)%l
        do k=2,atomix%species%norbs(i)
          if (tbMod%hopping(i,j)%l1 < atomix%speciesBasis(i,k)%l) then
            tbMod%hopping(i,j)%l1 = atomix%speciesBasis(i,k)%l
          endif
        enddo

        tbMod%hopping(i,j)%l2=atomix%speciesBasis(j,1)%l
        do k=2,atomix%species%norbs(j)
          if (tbMod%hopping(i,j)%l2<atomix%speciesBasis(j,k)%l) then
            tbMod%hopping(i,j)%l2=atomix%speciesBasis(j,k)%l
          endif
        enddo

        tbMod%hopping(i,j)%ll=min(tbMod%hopping(i,j)%l1,tbMod%hopping(i,j)%l2)
        allocate(tbMod%hopping(i,j)%a(0:tbMod%hopping(i,j)%l1,0:tbMod%hopping(i,j)%l2,0:tbMod%hopping(i,j)%ll))
        allocate(tbMod%hopping(i,j)%eps(0:atomix%species%norbs(i)-1))

      enddo
    enddo

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        tbMod%hopping(i,j)%eps=0.0_k_pr
        tbMod%hopping(i,j)%a1=0.0_k_pr
        tbMod%hopping(i,j)%a2=0.0_k_pr
        tbMod%hopping(i,j)%a3=0.0_k_pr
        tbMod%hopping(i,j)%a4=0.0_k_pr
        tbMod%hopping(i,j)%phi0=0.0_k_pr
        tbMod%hopping(i,j)%r0=0.0_k_pr
        tbMod%hopping(i,j)%rc=0.0_k_pr
        tbMod%hopping(i,j)%r1=0.0_k_pr
        tbMod%hopping(i,j)%rcut=0.0_k_pr
        tbMod%hopping(i,j)%n=0.0_k_pr
        tbMod%hopping(i,j)%nc=0.0_k_pr
        tbMod%hopping(i,j)%d0=0.0_k_pr
        tbMod%hopping(i,j)%dc=0.0_k_pr
        tbMod%hopping(i,j)%d1=0.0_k_pr
        tbMod%hopping(i,j)%dcut=0.0_k_pr
        tbMod%hopping(i,j)%m=0.0_k_pr
        tbMod%hopping(i,j)%mc=0.0_k_pr
        tbMod%hopping(i,j)%a=0.0_k_pr
      enddo
    enddo
    select case (gen%bond)
      case (k_bondGSP)
        call ReadTbGSP(io,gen,atomix,tbMod)
        call PrintTbGSP(io,gen,atomix,tbMod)
      case (k_bondHarrison)
        call ReadTbHarrison(io,gen,atomix,tbMod)
        call PrintTbHarrison(io,gen,atomix,tbMod)
    end select  

  end subroutine ReadTBModel


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

  subroutine ReadTbGSP(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="ReadTbGSP"
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomix
    type(modelType), intent(inout) :: tbMod
    integer :: i,j,i1,k,k1,k2,p
    character(len=k_mw) :: read_var

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        if (i == j) then
          if (gen%spin) then
            i1=0
            p=0
            do while (i1<=atomix%species%norbs(i)/2-1)
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=GetReal(io,trim(ccvar(i1,i1,read_var)),0.0_k_pr)
              tbMod%hopping(i,j)%eps(i1+atomix%species%norbs(i)/2:i1+atomix%species%norbs(i)/2+2*p)=tbMod%hopping(i,j)%eps(i1)
              i1=i1+2*p+1
              p=p+1
            enddo
          else
            i1=0
            p=0
            do while (i1<=atomix%species%norbs(i)-1)
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=GetReal(io,trim(ccvar(i1,i1,read_var)),0.0_k_pr)
              i1=i1+2*p+1
              p=p+1
            enddo
          endif
          if (gen%embedding) then
            read_var="a1"
            tbMod%hopping(i,j)%a1=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
            read_var="a2"
            tbMod%hopping(i,j)%a2=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
            read_var="a3"
            tbMod%hopping(i,j)%a3=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
            read_var="a4"
            tbMod%hopping(i,j)%a4=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
          endif
        endif
        read_var="phi0"
        tbMod%hopping(i,j)%phi0=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="r0"
        tbMod%hopping(i,j)%r0=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="rc"
        tbMod%hopping(i,j)%rc=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="r1"
        tbMod%hopping(i,j)%r1=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="rcut"
        tbMod%hopping(i,j)%rcut=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="n"
        tbMod%hopping(i,j)%n=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="nc"
        tbMod%hopping(i,j)%nc=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="d0"
        tbMod%hopping(i,j)%d0=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="dc"
        tbMod%hopping(i,j)%dc=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="d1"
        tbMod%hopping(i,j)%d1=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="dcut"
        tbMod%hopping(i,j)%dcut=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="m"
        tbMod%hopping(i,j)%m=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="mc"
        tbMod%hopping(i,j)%mc=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        do k=0,tbMod%hopping(i,j)%l1
          do k1=k,tbMod%hopping(i,j)%l2
            do k2=0,k
              tbMod%hopping(i,j)%a(k,k1,k2)=GetReal(io,trim(ccnlm(i,j,k,k1,k2)),0.0_k_pr)
            enddo
          enddo
        enddo
      enddo
    enddo

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
              if (k>k1) then
                tbMod%hopping(i,j)%a(k,k1,k2)=(-1.0_k_pr)**((k-k1+abs(k-k1))/2.0_k_pr)*tbMod%hopping(j,i)%a(k1,k,k2)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine ReadTbGSP

!> \brief selects the initialization routine for the Harrison tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:16:15
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine ReadTbHarrison(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="ReadTbHarrison"
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomix
    type(modelType), intent(inout) :: tbMod

    integer i,j,i1,k,k1,k2,p
    character(len=k_mw) :: read_var
    
    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        if (i==j) then
          if (gen%spin) then
            i1=0
            p=0
            do while (i1<=atomix%species%norbs(i)/2-1)
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=GetReal(io,trim(ccvar(i1,i1,read_var)),0.0_k_pr)
              tbMod%hopping(i,j)%eps(i1+atomix%species%norbs(i)/2:i1+atomix%species%norbs(i)/2+2*p)=tbMod%hopping(i,j)%eps(i1)
              i1=i1+2*p+1
              p=p+1
            enddo
          else
            i1=0
            p=0
            do while (i1<=atomix%species%norbs(i)/2-1)
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=GetReal(io,trim(ccvar(i1,i1,read_var)),0.0_k_pr)
              i1=i1+2*p+1
              p=p+1
            enddo
          endif
          if (gen%embedding) then
            read_var="a1"
            tbMod%hopping(i,j)%a1=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
            read_var="a2"
            tbMod%hopping(i,j)%a2=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
            read_var="a3"
            tbMod%hopping(i,j)%a3=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
            read_var="a4"
            tbMod%hopping(i,j)%a4=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
          endif
        endif

        read_var="phi0"
        tbMod%hopping(i,j)%phi0=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="r1"
        tbMod%hopping(i,j)%r1=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="rcut"
        tbMod%hopping(i,j)%rcut=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="n"
        tbMod%hopping(i,j)%n=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="d1"
        tbMod%hopping(i,j)%d1=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="dcut"
        tbMod%hopping(i,j)%dcut=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)
        read_var="m"
        tbMod%hopping(i,j)%m=GetReal(io,trim(ccvar(i,j,read_var)),0.0_k_pr)

        do k=0,tbMod%hopping(i,j)%l1
          do k1=k,tbMod%hopping(i,j)%l2
            do k2=0,k
              tbMod%hopping(i,j)%a(k,k1,k2)=GetReal(io,trim(ccnlm(i,j,k,k1,k2)),0.0_k_pr)
            enddo
          enddo
        enddo
      enddo
    enddo


    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
              if (k>k1) then
                tbMod%hopping(i,j)%a(k,k1,k2)=(-1.0_k_pr)**((k-k1+abs(k-k1))/2.0_k_pr)*tbMod%hopping(j,i)%a(k1,k,k2)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine ReadTbHarrison

!> \brief reads the delta block
!> \details the block is read only if multipoles k_electrostatics is asked for
!> \author Alin M Elena
!> \date 31/10/07, 17:52:20
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
!> \remarks Note that: delta (l,l1,l2) with |l1-l2|<=l<=|l1+l2| \n
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
  subroutine ReadDelta(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="ReadDelta"
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomix
    type(modelType), intent(inout) :: tbMod
    integer :: nt,sp,j,i,k,l,errno
    character(len=k_ml) :: saux
    integer, allocatable :: tmpid(:)
        
    if (GetBlock(io,"DeltaPol",nt)) then
      allocate(tbMod%delta(1:atomix%species%nspecies))
      allocate(tmpid(1:atomix%species%nspecies))
      tmpid=0
      do i = 1,atomix%species%nspecies
        read(nt,fmt=*,iostat=errno) sp
        if (errno/=0) then
          call error("Specie indicator missing, please fix it!(reading delta)",sMyName,.true.,io)
        endif
        if (.not. isInList(sp,atomix%species%id)) then
          write(saux,"(a,i0)")"undefined specie: ",sp
          call error(trim(saux),sMyName,.true.,io)
        endif
        if (isInList(sp,tmpid)) then
          write(saux,"(a,i0)")"specie has been already read: ",sp
          call error(trim(saux),sMyName,.true.,io)
        endif
        tmpid(i)=sp
        tbMod%delta(i)%sp=sp
        j=maxval(atomix%speciesBasis(tbMod%delta(i)%sp,1:atomix%species%norbs(tbMod%delta(i)%sp)/2)%l)
        tbMod%delta(i)%l=j
        allocate(tbMod%delta(i)%d(0:2*j,0:j,0:j))
        tbMod%delta(i)%d=0.0_k_pr
        do j=0,tbMod%delta(i)%l
          do k=j,tbMod%delta(i)%l
            do l=abs(j-k),j+k
              if (mod(k+j+l,2)==0) then
                if ((k==j).and.(l==0)) then
                  tbMod%delta(i)%d(0,k,j)=1.0_k_pr
                else
                  read(nt,*,iostat=errno) tbMod%delta(i)%d(l,j,k)
                  if (errno/=0) then
                    write(saux,'(a,i0,a,3i0)')"Error reading delta for specie ",sp," indices: ",l,j,k
                    call error(trim(saux),sMyName,.true.,io)
                  endif
                  tbMod%delta(i)%d(l,k,j)=tbMod%delta(i)%d(l,j,k)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
      else
        call error("DeltaPol block is missing",sMyName,.true.,io)
      endif
  end subroutine ReadDelta




!> \brief deallocates the memory
!> \author Alin M Elena
!> \date 29/10/07, 17:38:16
!> \param atomic type(atomicxType) contains atomic data
!> \param general type(generalType) general data that controls the flow of the program
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param io type(ioType) contains all the info about I/O files


        
  subroutine CleanMemory(io,atomic,general,tbMod,Sol)
    character(len=*), parameter :: sMyName="CleanMemory"
    type(ioType), intent(in) :: io
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: general
    type(modelType),intent(inout) :: tbMod
    type(solutionType), intent(inout) :: sol

    integer :: i,j

    deallocate(atomic%species%id)
    deallocate(atomic%species%mass)
    deallocate(atomic%species%z)
    deallocate(atomic%species%zval)
    deallocate(atomic%species%ulocal)
    deallocate(atomic%species%jlocal)
    deallocate(atomic%species%uinter)
    deallocate(atomic%species%norbs)

    deallocate(atomic%atoms%id)
    deallocate(atomic%atoms%sp)
    deallocate(atomic%atoms%bias)
    deallocate(atomic%atoms%isscf)
    deallocate(atomic%atoms%ismoving)
    deallocate(atomic%atoms%x)
    deallocate(atomic%atoms%y)
    deallocate(atomic%atoms%z)
    deallocate(atomic%atoms%vx)
    deallocate(atomic%atoms%vy)
    deallocate(atomic%atoms%vz)
    deallocate(atomic%atoms%fx)
    deallocate(atomic%atoms%fy)
    deallocate(atomic%atoms%fz)
    deallocate(atomic%atoms%dx)
    deallocate(atomic%atoms%dy)
    deallocate(atomic%atoms%dz)
    deallocate(atomic%atoms%xo)
    deallocate(atomic%atoms%yo)
    deallocate(atomic%atoms%zo)
    deallocate(atomic%atoms%vxo)
    deallocate(atomic%atoms%vyo)
    deallocate(atomic%atoms%vzo)
    deallocate(atomic%atoms%fxo)
    deallocate(atomic%atoms%fyo)
    deallocate(atomic%atoms%fzo)
    deallocate(atomic%atoms%chrg)
    deallocate(atomic%atoms%chrg0)
    deallocate(atomic%atoms%norbs)
    if (allocated(atomic%atoms%spacer)) deallocate(atomic%atoms%spacer)
    if (allocated(atomic%atoms%donor)) deallocate(atomic%atoms%donor)
    if (allocated(atomic%atoms%acceptor)) deallocate(atomic%atoms%acceptor)
    if (allocated(atomic%atoms%scf)) deallocate(atomic%atoms%scf)
    if (allocated(atomic%atoms%moving)) deallocate(atomic%atoms%moving)
    deallocate(atomic%speciesBasis)
    deallocate(atomic%basis%orbitals)
    deallocate(atomic%atoms%orbs)

    do i=1,atomic%species%nspecies
      do j=1,atomic%species%nspecies
        deallocate(tbMod%hopping(i,j)%a)
        deallocate(tbMod%hopping(i,j)%eps)
        deallocate(tbMod%hopping(i,j)%hmnTail)
      enddo
    enddo
    deallocate(tbMod%hopping)
    if(general%electrostatics==k_electrostaticsMultipoles) then
      do i=1,atomic%species%nspecies
        deallocate(tbMod%delta(i)%d)
      enddo
      deallocate(tbMod%delta)
    endif

    call DestroyMatrix(sol%h,io)
    call DestroyMatrix(sol%eigenvecs,io)
    deallocate(sol%eigenvals)

    if (general%scf) then
      call DestroyMatrix(sol%hin,io)
      call DestroyMatrix(sol%h2,io)
      deallocate(sol%potential)
      deallocate(sol%field)
    endif
    if (general%spin) then
      call DestroyMatrix(sol%hdown,io)
      call DestroyMatrix(sol%hup,io)
    endif
    call DestroyMatrix(sol%rho,io)
    if ((general%scf).and.(general%electrostatics==k_electrostaticsMultipoles)) then
      deallocate(sol%gcoeff)
      deallocate(sol%rgc)
    endif

    deallocate(sol%n0)
    deallocate(sol%density)

  end subroutine CleanMemory
end module m_ReadData