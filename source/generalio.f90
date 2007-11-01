!> \brief deals with reading the input file(s)
!> \details it saves the gathered info in the right variables
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \remarks

module read_data
  use constants
  use useful
  use parser
  use types
  use ISO_FORTRAN_ENV


  implicit none
  private
  public :: initialize
  public :: close_io_general
  public :: clean_memory

contains


!> \brief reads from io_loc different fields of gen_loc
!> \author Alin M Elena
!> \date ~2007
!> \param io_loc type(io_type) contains all the info about I/O files
!> \param gen_loc type(general_type) contains the info needed by the program to run
!> \param atomic type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters

  subroutine initialize(io_loc,gen_loc,atomic,tbMod)
    character(len=*), parameter :: myname = 'initialize'
    type(io_type), intent(inout) :: io_loc
    type(general_type), intent(inout) :: gen_loc
    type(atomicx_type), intent(inout) :: atomic
    type(model_type), intent(inout) :: tbMod
    real(pr) :: aux
    integer :: errno=-1,nt,i

    io_loc%uerr=get_unit()
    io_loc%inp_err=trim(io_loc%inp_file)//".err"
    open(unit=io_loc%uerr,file=trim(io_loc%inp_err),status="replace",&
      action="write",iostat=errno)
    if (errno /= 0) then
      write(*,*)"I can not create error file ", trim(io_loc%inp_file)//".err"
      stop
    endif 
    io_loc%uinp=get_unit()
    open(unit=io_loc%uinp,file=io_loc%inp_file,status='old',&
      action='read',iostat=errno)

    if (errno /= 0) then
      write(io_loc%uerr,*)"I can not open file ", io_loc%inp_file
      stop
    endif


! parse file and in the same time 
! makes a minimal check of corectness
    call cpu_time(gen_loc%time%int)

    call parse_file(io_loc)
    call cpu_time(gen_loc%time%end)
    aux=gen_loc%time%end-gen_loc%time%int

    call cpu_time(gen_loc%time%int)

    call read_io(io_loc)

    if (io_loc%debug>=medium_debug) then
      call cpu_time(gen_loc%time%end)
      write(io_loc%udeb,'(a,f16.6,a)')"Parsing time "&
        ,aux," seconds"
      write(io_loc%udeb,'(a,f16.6,a)')"Setting io_info "&
        ,gen_loc%time%end-gen_loc%time%int," seconds"
    endif

    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,1)
    call read_general(io_loc,gen_loc)
    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,2)

    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,1)
    call read_species(io_loc,gen_loc,atomic%species)
    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,2)

    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,1)
    call read_atoms(io_loc,gen_loc,atomic)
    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,2)

    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,1)
    call read_basis(io_loc,gen_loc,atomic)
    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,2)

    call print_species(io_loc,atomic%species)
    call print_atoms(io_loc,gen_loc,atomic)
    call print_basis(io_loc,gen_loc,atomic)

    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,1)
    call read_tb_model(io_loc,gen_loc,atomic,tbMod)
    if (gen_loc%electrostatics==electrostatics_multipoles) then
      call read_delta(io_loc,gen_loc,atomic,tbMod)
    endif
    if (io_loc%debug>=medium_debug) call timing(io_loc,gen_loc%time,2)

    call end_parse
  end subroutine initialize

!> \brief reads the names for I/O files and the output levels
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param io_loc type(io_type) 
!> \details
!> A short description of the input file io_info part
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
!> OnScreen & logical & .false. & on/off printing on the screen, on means that nothing will be written in the output file \\ 
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
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off printing on the screen, on
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


subroutine read_io(io_loc)
character(len=*),parameter :: myname="read_io"
  type(io_type), intent(inout) :: io_loc
  integer :: errno


!comm_io DebugFile & string & input file name.dbg & Name of the debug file \\

! read the name of output file plus debug (y/n) and output level
  io_loc%debug_file=get_string(io_loc,"DebugFile",&
    trim(io_loc%inp_file)//".dbg",.false.)
  io_loc%udeb=get_unit()
  open(unit=io_loc%udeb,file=trim(io_loc%debug_file),status="replace",&
    action="write",iostat=errno)
  if (errno /= 0) call error("I can not create file "//&
    trim(io_loc%debug_file),myname,.true.,io_loc)

!comm_io OutputFile & string & input file name.out & Name of the output file \\
!comm_io OnScreen & logical & .false. & on/off printing on the screen, on means that nothing will be written in the output file \\
  io_loc%output_file=get_string(io_loc,"OutputFile",&
    trim(io_loc%inp_file)//".out")
  io_loc%stdout=get_logical(io_loc,"OnScreen",.false.)
! prepares the output according to the settings from input file(s)
  if (io_loc%stdout) then
    io_loc%uout=6
  else
    io_loc%uout=get_unit()
    open(unit=io_loc%uout,file=trim(io_loc%output_file),status="replace",&
      action="write",iostat=errno)
    if (errno /= 0) call error("I can not create file "//trim(io_loc%output_file)&
      ,myname,.true.,io_loc)
  endif

! levels
  io_loc%verbosity=get_integer(io_loc,"OutputLevel",5)
  io_loc%debug=get_integer(io_loc,"DebugLevel",5)
  io_loc%first_time=.not. io_loc%first_time

  if (io_loc%verbosity>=high_debug) then
    write( io_loc%udeb,'(a,l1)')"Read io_info for the first time: "&
      ,io_loc%first_time
  endif

  write( io_loc%uout,'(a,a)')"Input file: "&
    ,trim(io_loc%inp_file)
  write( io_loc%uout,'(a,a)')"Output file(OutputFile): "&
    ,trim(io_loc%output_file)
  write( io_loc%uout,'(a,i0)')"Output Level(OutputLevel): "&
    ,io_loc%verbosity  
  write( io_loc%uout,'(a,a)')"Debug file(DebugFile): "&
    ,trim(io_loc%debug_file)
  write( io_loc%uout,'(a,i0)')"Debug Level(DebugLevel): "&
    ,io_loc%debug
  write( io_loc%uout,'(a,a)')"Error file: "&
    ,trim(io_loc%inp_err)
  write( io_loc%uout,'(a,l1)')"Standard Output(OnScreen): "&
    ,io_loc%stdout


end subroutine read_io


!> \brief reads the parameters necessary for running
!> \author Alin M Elena
!> \date 21st of January, 2007 
!> \param  io_loc type(io_type) I/O details (see types::io_type)
!> \param gen_loc type(general_type), keeps all the general information about the parameters of the program
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

subroutine read_general(io_loc,gen_loc)
  character(len=*),parameter :: name="read_general"
  type(io_type), intent(inout) :: io_loc
  type(general_type), intent(inout) :: gen_loc

  character(len=mw) :: saux

    if (.not. gen_loc%first_time) gen_loc%first_time=.not. gen_loc%first_time

    gen_loc%job_name=get_string(io_loc,"JobName","no name")
    write( io_loc%uout,'(a,a)')"Job Name(JobName): "&
      ,gen_loc%job_name

    gen_loc%ranseed=get_real(io_loc,"RanSeed",0.5_pr)
    write( io_loc%uout,'(a,f0.16)')"Random number generator seed(RanSeed): "&
      ,gen_loc%ranseed

    gen_loc%read_velocity=get_logical(io_loc,"ReadVel",.false.)
    write( io_loc%uout,'(a,l1)')"Read Velocity(ReadVel): "&
     ,gen_loc%read_velocity

!comm_gen VBias & real & 0.0 & bias factor\\
    gen_loc%bias=get_real(io_loc,"VBias",0.0_pr)
    write( io_loc%uout,'(a,f0.8)')"Bias(VBias): "&
     ,gen_loc%bias


!comm_gen MaxOrbsPerAtom & integer & 8 & maximum number of orbitals per atom\\
    gen_loc%max_orbitals_per_atom=get_integer(io_loc,"MaxOrbsPerAtom",8)
    write( io_loc%uout,'(a,i0)')"Maximum Orbitals Per Atom (MaxOrbsPerAtom): "&
      ,gen_loc%max_orbitals_per_atom

    gen_loc%spin=get_logical(io_loc,"Spin",.false.)
    write( io_loc%uout,'(a,l1)')"spin polarisation(Spin): "&
      ,gen_loc%spin

    gen_loc%scf=get_logical(io_loc,"SCF",.false.)
    write( io_loc%uout,'(a,l1)')"Is SCF?(SCF): "&
      ,gen_loc%scf

!comm_gen Units & AU EV SI & AU & system of units: atomic units(AU), electronVolt-Angstrom(eVA), International(SI)\\
    saux=get_string(io_loc,"Units","AU")
    write( io_loc%uout,'(a,a)')"System of Units (Units): "&
      ,trim(saux)
    if (cstr(trim(saux),'AU')) then
      gen_loc%units = units_au
    elseif (cstr(trim(saux),'eVA')) then
      gen_loc%units = units_ev
    elseif (cstr(trim(saux),'SI')) then
      gen_loc%units = units_si
    else
      call error("The requested system of units is not implemented",name,.true.,io_loc)
    endif

 !comm_gen BondType & Harrison, GSP &Harrison& bond type \\
! \hline
    saux=get_string(io_loc,"BondType","Harrison")
    write( io_loc%uout,'(a,a)')"bond type (BondType): "&
      ,trim(saux)
    if (cstr(trim(saux),'Harrison')) then
      gen_loc%bond = bond_harrison
    elseif (cstr(trim(saux),'GSP')) then
      gen_loc%bond = bond_gsp
    else
      call error("The requested bond type is not implemented",name,.true.,io_loc)
    endif

   !comm_gen Embedding & logical & .true. & on/off embedding method\\
    gen_loc%embedding=get_logical(io_loc,"Embedding",.true.)
    write( io_loc%uout,'(a,l1)')"has embedding(Embedding): "&
     ,gen_loc%embedding

!comm_gen Electrostatics & PointCharges, Multipoles & PointCharges & method use to compute electrostatic interaction\\
    saux=get_string(io_loc,"Electrostatics","PointCharges")
    write( io_loc%uout,'(a,a)')"Electrostatics type (Electrostatics): "&
     ,trim(saux)
    if (cstr(trim(saux),'PointCharges')) then
      gen_loc%electrostatics = electrostatics_point
    elseif (cstr(trim(saux),'Multipoles')) then
      gen_loc%electrostatics = electrostatics_multipoles
    else
      call error("The requested Electrostatics type is not implemented",name,.true.,io_loc)
    endif
! 
!comm_gen PrecomputeMultipoles & logical & .true. & on/off precompute multipoles\\
    gen_loc%comp_elec=get_logical(io_loc,"PrecomputeMultipoles",.true.)
    write( io_loc%uout,'(a,l1)')"precompute multipoles(PrecomputeMultipoles): "&
      ,gen_loc%comp_elec
! 

! !comm_gen WriteAni & logical & .false. & on/off writing animation file \\
!   gen_loc%write_ani=get_logical(io_loc,"WriteAni",.false.)
!   write( io_loc%uout,'(a,l)')"Write animation(WriteAni): "&
!     ,gen_loc%write_ani
! 
! !comm_gen WriteEne & logical & .true. & on/off writing energy file \\
!   gen_loc%write_ene=get_logical(io_loc,"WriteEne",.true.)
!   write( io_loc%uout,'(a,l)')"Write Energy(WriteEne): "&
!     ,gen_loc%write_ene
! 
! !comm_gen WriteDen & logical & .false. & on/off writing density file \\
!   gen_loc%write_density=get_logical(io_loc,"WriteDen",.false.)
!   write( io_loc%uout,'(a,l)')"Write density(WriteDen): "&
!     ,gen_loc%write_density
! !comm_gen ReadDen & logical & .false. & on/off reading density file \\
!   gen_loc%read_density=get_logical(io_loc,"ReadDen",.false.)
!   write( io_loc%uout,'(a,l)')"Read Density(ReadDen): "&
!     ,gen_loc%read_density
! 
!
! !comm_gen IonicTemperature & real & 300.0 & Ionic temperature\\
!   gen_loc%ionic_temperature=get_real(io_loc,"IonicTemperature",300.0_pr)
!   write( io_loc%uout,'(a,g)')"Ionic Temperature(IonicTemperature): "&
!     ,gen_loc%ionic_temperature
! 
! !comm_gen NetCharge & real & 0.0 & Net charge on the system\\
!   gen_loc%netcharge=get_real(io_loc,"NetCharge",0.0_pr)
!   write( io_loc%uout,'(a,g)')"Net charge(NetCharge): "&
!     ,gen_loc%netcharge
! 
! !comm_gen ElectronicTemperature & real & 300.0 & Electronic temperature, used to compute
!       occupation numbers if you choose Fermi-Dirac or Methfessel-Paxton method\\
!   gen_loc%electronic_temperature=get_real(io_loc,"ElectronicTemperature",300.0_pr)
!   write( io_loc%uout,'(a,g)')"Electronic Temperature(ElectronicTemperature): "&
!     ,gen_loc%electronic_temperature
! 
! !comm_gen ElectronicMu & real & 0.0 & chemical potential, used to compute occupation numbers if you choose constant $\mu$ method\\
!   gen_loc%electronic_mu=get_real(io_loc,"ElectronicMu",0.0_pr)
!   write( io_loc%uout,'(a,g)')"Chemical Potential(ElectronicMu): "&
!     ,gen_loc%electronic_mu
! 
! !comm_gen DeltaT & real & 0.001 & time step used to evolve equtions of motion\\
!   gen_loc%deltat=get_real(io_loc,"DeltaT",0.001_pr)
!   write( io_loc%uout,'(a,g)')"Time step(DeltaT): "&
!     ,gen_loc%deltat
! 
! !comm_gen Nsteps & integer & 100 & number of steps used to evolve equtions of motion\\
!   gen_loc%nsteps=get_integer(io_loc,"Nsteps",100)
!   write( io_loc%uout,'(a,i0)')"Number of steps(Nsteps): "&
!     ,gen_loc%nsteps
! 
! !comm_gen RunType & SinglePoint, BODynamics, Ehrenfest, EhrenfestDamped, Fit, ForceTest, ForceTestX,
!! ForceTestY, ForceTestZ & SinglePoint & Type of calculation\\
!   saux=get_string(io_loc,"RunType","SinglePoint")
!   write( io_loc%uout,'(a,a)')"Calculation type(RunType): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'SinglePoint')) then
!     gen_loc%runtype = run_sp
!   elseif (cstr(trim(saux),'bodynamics')) then
!     gen_loc%runtype = run_bo
!   elseif (cstr(trim(saux),'ehrenfest')) then
!     gen_loc%runtype = run_ehrenfest
!   elseif (cstr(trim(saux),'fit')) then
!     gen_loc%runtype = run_fit
!   elseif (cstr(trim(saux),'forcetest')) then
!     gen_loc%runtype = run_force_test
!   elseif (cstr(trim(saux),'forcetestx')) then
!     gen_loc%runtype = run_force_testx
!   elseif (cstr(trim(saux),'forcetesty')) then
!     gen_loc%runtype = run_force_testy
!   elseif (cstr(trim(saux),'forcetestz')) then
!     gen_loc%runtype = run_force_testz
!   elseif (cstr(trim(saux),'ehrenfestdamped')) then
!     gen_loc%runtype = run_ehrenfest_damped
!   else 
!     call error("The requested RunType is not implemented",name,.true.,io_loc)
!   endif
! 

! 
! !comm_gen SCFType & TB+UJ & TB+UJ & SCF method\\
!   saux=get_string(io_loc,"SCFType","TB+UJ")
!   write( io_loc%uout,'(a,a)')"SCF type(SCFType): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'TB+UJ')) then
!     gen_loc%scf_type = scf_tbuj
!   else 
!     call error("The requested SCFType is not implemented",name,.true.,io_loc)
!   endif
! 
! !comm_gen SCFSteps & integer & 100 & maximum number of steps used for scf\\
!   gen_loc%maxscf=get_integer(io_loc,"SCFSteps",100)
!   write( io_loc%uout,'(a,i0)')"SCF Maximum number of steps(SCFSteps): "&
!     ,gen_loc%maxscf
! !comm_gen SCFMix & real & 0.85 & mixing parameter\\
!   gen_loc%scfmix=get_real(io_loc,"SCFMix",0.85_pr)
!   write( io_loc%uout,'(a,g)')"SCF Mixing parameter(SCFMix): "&
!     ,gen_loc%scfmix
! 
! !comm_gen SCFMixType & Broyden, Pulay & Broyden & SCF mixing method\\
!   saux=get_string(io_loc,"SCFMixType","Broyden")
!   write( io_loc%uout,'(a,a)')"SCF Mixing type(SCFMixType): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'Broyden')) then
!     gen_loc%scfmixtype = scfmix_broyden
!   elseif (cstr(trim(saux),'Pulay')) then
!     gen_loc%scfmixtype = scfmix_pulay
!   else  
!     call error("The requested SCFMixType is not implemented",name,.true.,io_loc)
!   endif
! 
! !comm_gen SCFTol & real & 1e-8 & convergence tolerance\\
!   gen_loc%scftol=get_real(io_loc,"SCFTol",1.0e-8_pr)
!   write( io_loc%uout,'(a,g)')"SCF convergence tolerance(SCFTol): "&
!     ,gen_loc%scftol
! 
! !comm_gen SCFMixN & integer & 4 & number of iterations to mix\\
!   gen_loc%scfmixn=get_integer(io_loc,"SCFMixN",4)
!   write( io_loc%uout,'(a,i0)')"number of iterations to mix(SCFMixN): "&
!     ,gen_loc%scfmixn
! 
! !comm_gen VelScale & logical & .false. & on/off scaling velocities\\
!   gen_loc%velocity_scaling=get_logical(io_loc,"VelScale",.false.)
!   write( io_loc%uout,'(a,l)')"scaling  Velocities(VelScale): "&
!     ,gen_loc%velocity_scaling
! !comm_gen DMOccTol & real & 1e-10 & density matrix occupation tolerance\\
!   gen_loc%dm_occupation_tolerance=get_real(io_loc,"DMOccTol",1.0e-10_pr)
!   write( io_loc%uout,'(a,g)')"density matrix occupation tolerance(DMOccTol): "&
!     ,gen_loc%dm_occupation_tolerance
! !comm_gen HElThres & real & 1e-10 & hamiltionian element thresold. Any element smaller that the thresold is made zero.\\
!   gen_loc%h_element_threshold=get_real(io_loc,"DMOccTol",1.0e-10_pr)
!   write( io_loc%uout,'(a,g)')"hamiltionian element thresold(HElThres): "&
!     ,gen_loc%h_element_threshold
! 


! 
! !comm_gen CollinearSpins & logical & .false. & on/off collinear spins\\
!   gen_loc%collinear=get_logical(io_loc,"CollinearSpins",.false.)
!   write( io_loc%uout,'(a,l)')"collinear spins(CollinearSpins): "&
!     ,gen_loc%collinear
! 
! !comm_gen SpinDU & real & 0.0 & spin down spin up difference\\
!   gen_loc%sdu=get_real(io_loc,"SpinDU",0.0_pr)
!   write( io_loc%uout,'(a,g)')"spin down spin up difference(SpinDU): "&
!     ,gen_loc%sdu
! 
! !comm_gen EulerSteps & integer & 100 & after each EulerSteps apply an Euler integration of equations of motions\\
!   gen_loc%euler_Steps=get_integer(io_loc,"EulerSteps",100)
!   write( io_loc%uout,'(a,i0)')"Euler steps (EulerSteps): "&
!     ,gen_loc%euler_steps
! 

! 
! !comm_gen FSteps & integer & 100 & Number of steps used to test the force/energy consitency\\
!   gen_loc%f_steps=get_integer(io_loc,"FSteps",100)
!   write( io_loc%uout,'(a,i0)')"force steps (FSteps): "&
!     ,gen_loc%f_steps
! 
! !comm_gen FStart & real & 0.0 & position at which the force/energy consistency test starts\\
!   gen_loc%f_start=get_real(io_loc,"FStart",0.0_pr)
!   write( io_loc%uout,'(a,g)')"initial position for force(FStart): "&
!     ,gen_loc%f_start
! 
! !comm_gen Fdx & real & 0.001 & space step for the force/energy consistency test starts\\
!   gen_loc%f_dx=get_real(io_loc,"FStart",0.001_pr)
!   write( io_loc%uout,'(a,g)')"space step for force(Fdx): "&
!     ,gen_loc%f_dx
! 
! !comm_gen Gamma & real & 0.3 & damping factor for Ehrenfest equation\\
!   gen_loc%gamma=get_real(io_loc,"Gamma",0.5_pr)
!   write( io_loc%uout,'(a,g)')"damping factor in Ehrenfest Equation(Gamma): "&
!     ,gen_loc%gamma
! 
! !comm_gen MPN & integer & 2 & the order used for Methfessel-Paxton method\\
!   gen_loc%mp_N=get_integer(io_loc,"MPN",2)
!   write( io_loc%uout,'(a,i0)')"Order for Methfessel-Paxton method(MPN): "&
!     ,gen_loc%mp_N
! 
! !comm_gen Hole & integer & 0 & no of level to create hole\\
!   gen_loc%hole=get_integer(io_loc,"Hole",0)
!   write( io_loc%uout,'(a,i0)')"Hole level(Hole): "&
!     ,gen_loc%hole
! 
! !comm_gen Excite & integer & 0 & no of level to create excitation\\
!   gen_loc%excite=get_integer(io_loc,"Excite",0)
!   write( io_loc%uout,'(a,i0)')"Excitation level(Excite): "&
!     ,gen_loc%excite
! !comm_gen HoleSpin & D,U & D & spin of the hole\\
!   saux=get_string(io_loc,"HoleSpin","D")
!   write( io_loc%uout,'(a,a)')"Spin of the Hole (HoleSpin): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'D')) then
!     gen_loc%hole_spin = spin_down
!   elseif (cstr(trim(saux),'U')) then
!     gen_loc%hole_spin = spin_up
!   else  
!     call error("The requested spin type is not implemented",name,.true.,io_loc)
!   endif
! 
! 
! !comm_gen ExciteSpin & D,U & D & spin of the excitation\\
!   saux=get_string(io_loc,"ExciteSpin","D")
!   write( io_loc%uout,'(a,a)')"Spin of the excitation (ExciteSpin): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'D')) then
!     gen_loc%excite_spin = spin_down
!   elseif (cstr(trim(saux),'U')) then
!     gen_loc%excite_spin = spin_up
!   else  
!     call error("The requested spin type is not implemented",name,.true.,io_loc)
!   endif
! 
! 
! !comm_gen Units & AU EV SI & AU & system of units atomic units, electronVolt-Angstrom, International\\
!   saux=get_string(io_loc,"Units","AU")
!   write( io_loc%uout,'(a,a)')"System of Units (Units): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'Au')) then
!     gen_loc%units = units_au
!   elseif (cstr(trim(saux),'EV')) then
!     gen_loc%units = units_ev
!   elseif (cstr(trim(saux),'SI')) then
!     gen_loc%units = units_si
!   else  
!     call error("The requested system of units is not implemented",name,.true.,io_loc)
!   endif
! 
! !comm_gen SmearingMethod & FD,MP,CMU & FD & smearing method (Femi-Dirac, Methfessel-Paxton, constant $\mu$)\\
!   saux=get_string(io_loc,"SmearingMethod","FD")
!   write( io_loc%uout,'(a,a)')"Smearing Method (SmearingMethod): "&
!     ,trim(saux)
!   if (cstr(trim(saux),'FD')) then
!     gen_loc%smethod = sm_fd
!   elseif (cstr(trim(saux),'MP')) then
!     gen_loc%smethod = sm_mp
!   elseif (cstr(trim(saux),'CMU')) then
!     gen_loc%smethod = sm_cmu
!   else  
!     call error("The requested smearing method is not implemented",name,.true.,io_loc)
!   endif!

end subroutine read_general


!> \brief closes all the I/O units
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param io_loc type(io_type)

subroutine close_io_general(io_loc)
  character(len=*), parameter :: myname = 'close_io_general'
  type(io_type), intent(inout) :: io_loc
  integer :: errno

  if (.not.io_loc%stdout) then
    close(unit=io_loc%uout,iostat=errno)
    if (errno /= 0) call error("I can not close file "//trim(io_loc%output_file)&
      ,myname,.false.,io_loc)
  endif
  close(unit=io_loc%udeb,iostat=errno)
  if (errno /= 0) call error("I can not close file "//&
    trim(io_loc%debug_file),myname,.false.,io_loc)
  close(unit=io_loc%uerr,iostat=errno)
  if (errno /= 0) call error("I can not close file "//&
    trim(io_loc%inp_err),myname,.false.,io_loc)
end subroutine close_io_general


!> \brief times different events
!> \details it should be called in pairs with stage = 1 to initialize 
!> and stage = 2 to write down the info
!> \author Alin M Elena
!> \date 29th of October 2007
!> \param io type(io_type) i/o units
!> \param times type(time_type) time keepers
!> \param stage selects the stage

  subroutine timing(io,times,stage)
    character(len=*),parameter :: myname = "timing"
    type(io_type), intent(in) :: io
    type(time_type), intent(inout) :: times
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
!> \param io type(io_type) i/o units
!> \param general type(general_type) general data
!> \param atomix type(atomicx_type) contains info about atoms
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
!> The missing ones have the velocities initialized with zero. Each atom has a line\n
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

  subroutine read_atoms(io,general,atomix)
    character(len=*), parameter :: sMyName="read_atoms"
    type(io_type), intent(inout) :: io
    type(general_type), intent(in)  :: general
    type(atomicx_type), intent(inout) :: atomix
    ! arguments of the subroutine
    integer :: nt, errno,i,k, atom_id
    character(len=ml) :: saux
    integer, allocatable :: ids(:)
    real(pr) :: vx,vy,vz

    atomix%atoms%natoms = get_integer(io,"NumberOfAtoms",-1)
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

    atomix%atoms%vx=0.0_pr
    atomix%atoms%vy=0.0_pr
    atomix%atoms%vz=0.0_pr
    atomix%atoms%vxo=0.0_pr
    atomix%atoms%vyo=0.0_pr
    atomix%atoms%vzo=0.0_pr
    atomix%atoms%xo=0.0_pr
    atomix%atoms%yo=0.0_pr
    atomix%atoms%zo=0.0_pr
    atomix%atoms%dx=0.0_pr
    atomix%atoms%dy=0.0_pr
    atomix%atoms%dz=0.0_pr
    atomix%atoms%fx=0.0_pr
    atomix%atoms%fy=0.0_pr
    atomix%atoms%fz=0.0_pr
    atomix%atoms%fxo=0.0_pr
    atomix%atoms%fyo=0.0_pr
    atomix%atoms%fzo=0.0_pr
    atomix%atoms%chrg=0.0_pr
    atomix%atoms%chrg0=0.0_pr
    atomix%atoms%ismoving=.true.
    atomix%atoms%isscf=.true.
    atomix%atoms%nmoving=0
    atomix%atoms%nscf=0
    atomix%atoms%norbs=0

    if (get_block(io,"AtomsData",nt)) then
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

    if (general%read_velocity) then
      if (get_block(io,"VelocitiesData",nt)) then
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
    atomix%atoms%nacceptor=get_integer(io,"NAcceptor",0)
    atomix%atoms%ndonor=get_integer(io,"Ndonor",0)

    if (atomix%atoms%ndonor/=0) then
      allocate(atomix%atoms%donor(1:atomix%atoms%ndonor))
      atomix%atoms%donor=0
      if (get_block(io,'DonorAtoms',nt)) then
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
      if (get_block(io,'AcceptorAtoms',nt)) then
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
  end subroutine read_atoms



!> \brief prints coordinates, velocities and type of atoms
!> \author Alin M Elena
!> \date 30/10/07, 13:22:04
!> \param io type(io_type) i/o units
!> \param general type(general_type) general data
!> \param atomix type(atomicx_type) contains info about atoms

  subroutine print_atoms(io,general,atomix)
    character(len=*), parameter :: sMyName="print_atoms"
    type(io_type), intent(inout) :: io
    type(general_type), intent(in)  :: general
    type(atomicx_type), intent(in) :: atomix
    integer :: i

    write(io%uout,'(a)')  "=StartAtomsData==============================================================================="
    write(io%uout,'(a)')  "      Id  Sp El               X               Y               Z            Bias IsSCF IsMoving"
    write(io%uout,'(a)')"----------------------------------------------------------------------------------------------"
    do i=1,atomix%atoms%natoms
          write(io%uout,'(i8,i4,1x,a2,4f16.8,5x,l1,8x,l1)')atomix%atoms%id(i),atomix%atoms%sp(i),symbol(get_Z(atomix,i)),&
              atomix%atoms%x(i),atomix%atoms%y(i),atomix%atoms%z(i),atomix%atoms%bias(i),&
              atomix%atoms%isscf(i),atomix%atoms%ismoving(i)
    enddo
    write(io%uout,'(a)')   "=EndAtomsData================================================================================="
    if (general%read_velocity) then
      write(io%uout,'(a)')  "=StartVelocitiesData==========================================="
      write(io%uout,'(a)')  "      Id  Sp El              VX              VY              VZ"
      write(io%uout,'(a)')  "---------------------------------------------------------------"
      do i=1,atomix%atoms%natoms
            write(io%uout,'(i8,i4,1x,a2,3f16.8,5x)')atomix%atoms%id(i),atomix%atoms%sp(i),symbol(get_Z(atomix,i)),&
              atomix%atoms%vx(i),atomix%atoms%vy(i),atomix%atoms%vz(i)
      enddo
      write(io%uout,'(a)')  "=EndVelocitiesData============================================="
    endif
    write(io%uout,'(a)')"=AtomicLists=========================================="
    if(atomix%atoms%nscf/=0) then
      write(io%uout,'(a)')"SCFAtoms:"
      do i=1,atomix%atoms%nscf
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%scf(i),"(",symbol(get_Z(atomix,atomix%atoms%scf(i))),") "
      enddo
      write(io%uout,*)
    endif
    if(atomix%atoms%nmoving/=0) then
      write(io%uout,'(a)')"Moving Atoms:"
      do i=1,atomix%atoms%nmoving
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%moving(i),"(",symbol(get_Z(atomix,atomix%atoms%moving(i))),") "
      enddo
      write(io%uout,*)
    endif
    write(io%uout,'(a,i0)')"Atoms in Donor Group(NDonor): "&
        ,atomix%atoms%ndonor
    if(atomix%atoms%ndonor/=0) then
      write(io%uout,'(a)')"Donor Atoms:"
      do i=1,atomix%atoms%ndonor
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%donor(i),"(",symbol(get_Z(atomix,atomix%atoms%donor(i))),") "
      enddo
      write(io%uout,*)
    endif
    write( io%uout,'(a,i0)')"Atoms in Acceptor Group(NAcceptor): "&
        ,atomix%atoms%nacceptor
    if(atomix%atoms%nacceptor/=0) then
      write(io%uout,'(a)')"Acceptor Atoms:"
      do i=1,atomix%atoms%nacceptor
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%acceptor(i),"(",symbol(get_Z(atomix,atomix%atoms%acceptor(i))),") "
      enddo
      write(io%uout,*)
    endif
    write(io%uout,'(a,i0)')"Atoms in Spacer Group(NSpacer): "&
        ,atomix%atoms%nspacer
    if(atomix%atoms%nspacer/=0) then
      write(io%uout,'(a)')"Spacer Atoms:"
      do i=1,atomix%atoms%nspacer
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%spacer(i),"(",symbol(get_Z(atomix,atomix%atoms%spacer(i))),") "
      enddo
      write(io%uout,*)
    endif
    write(io%uout,'(a)')"=EndLists============================================="
  end subroutine print_atoms

!> \brief returns the atomic no Z
!> \author Alin M Elena
!> \date 30th of October 2007
!> \param atomic type(atomicx_type) data about atoms
!> \param atom integer the id of an atom
  integer function get_Z(atomic,atom)
    character(len=*), parameter :: sMyName="get_Z"
    type(atomicx_type), intent(in) :: atomic
    integer, intent(in) :: atom

    get_Z=atomic%species%z(get_specie(atomic,atom))
  end function get_Z

!> \brief returns the id of the specie of an atom
!> \author Alin M Elena
!> \date 30th of October 2007
!> \param atomic type(atomicx_type) data about atoms
!> \param atom integer the id of an atom
  integer function get_specie(atomic,atom)
    character(len=*), parameter :: sMyName="get_specie"
    type(atomicx_type), intent(in) :: atomic
    integer, intent(in) :: atom

    get_specie=atomic%atoms%sp(atom)
  end function get_specie

!> \brief reads and allocates the info about species of atoms 
!> \details this is all that should be done in this module.
!> the rest should happen in a specialized module
!> \author Alin M Elena
!> \date 29/10/07, 17:49:04
!> \param io type(io_type) i/o units
!> \param general type(general_type) general data
!> \param specs type(species_type) contains info about atoms
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

  subroutine read_species(io,general,specs)
    character(len=*), parameter :: sMyName="read_species"
    type(io_type), intent(inout) :: io
    type(general_type), intent(in)  :: general
    type(species_type), intent(inout) :: specs
    !internal variables
    integer :: nt,i,errno,sp
    character(len=ml) :: saux

    specs%nspecies=get_integer(io,"NumberOfSpecies",-1)
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
    if (get_block(io,"SpeciesData",nt)) then
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
    
  end subroutine read_species


!> \brief prints the info about species
!> \details at high verbosity some intristing info is printed about each specie
!> \author Alin M Elena
!> \date 29/10/07, 19:59:04
!> \param io type(io_type) i/o units
!> \param specs type(species_type) species data

subroutine print_species(io,specs)
  character(len=*), parameter :: sMyName="print_species"
  type(io_type), intent(in) :: io
  type(species_type), intent(inout) :: specs
  integer :: i
  write(io%uout,'(a)')"=SpeciesData============================================================================"
  write(io%uout,'(a)')"  Id  Z el    Zval            Mass          Ulocal          Jlocal          Uinter Norbs"
  write(io%uout,'(a)')"----------------------------------------------------------------------------------------"
  do i=1,specs%nspecies
    write(io%uout,'(i4,i3,1x,a,f8.4,4f16.8,i5)')specs%id(i),specs%z(i),symbol(specs%z(i)),specs%zval(i),specs%mass(i),&
            specs%ulocal(i),specs%jlocal(i),specs%uinter(i),specs%norbs(i)
    specs%mass(i)=specs%mass(i)*amu_to_internal
  enddo
  write(io%uout,'(a)')"=EndSpeciesData========================================================================="

  if (io%verbosity>=high_verbos) then
    do i=1, specs%nspecies
      write(io%uout,'(a)')"========================================================"
      write(io%uout,'(a,i0,a)')"Specie: ",specs%id(i)," detailed info"
      write(io%uout,'(a,a,a,a)')"name: ",trim(el_name(specs%z(i)))," symbol: ",trim(symbol(specs%z(i)))
      write(io%uout,'(a,i0,a,f0.3,a,f0.3,a,f0.3,a)')"Z=",specs%z(i)," Zvalence=",specs%zval(i),&
        " mass=",specs%mass(i)," (",weight(specs%z(i)),") a.u."
      write(io%uout,'(a,a,a,i0,a,a)')"period: ",trim(period(specs%z(i)))," group: ",group(specs%z(i)),&
        " electronic config:",trim(el_config(specs%z(i)))
      write(io%uout,'(a,f0.3,a,f0.3,a)')"van der Waals radius=",vdw_radius(specs%z(i)),&
        "A covalent radius=",covalent_radius(specs%z(i)),"A"
    enddo
    write(io%uout,'(a)')"========================================================"
  endif
end subroutine print_species



!> \brief reads and initializes the basis set
!> \author Alin M Elena
!> \date 30/10/07, 16:43:38
!> \param io type(io_type) i/o units
!> \param gen type(general_type) general data
!> \param atomix type(atomicx_type) contains info about atoms
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
  subroutine read_basis(io,gen,atomix)
    character(len=*), parameter :: sMyName="read_basis"
    type(io_type), intent(inout) :: io
    type(general_type), intent(in)  :: gen
    type(atomicx_type), intent(inout) :: atomix

    integer :: i,j,nt,errno,sp, norbitals
    character(len=ml) :: saux
    integer :: n,l,m,k
    real(pr) :: qu,qd
    integer, allocatable :: tmpid(:)

    allocate(atomix%species_basis(1:atomix%species%nspecies,1:gen%max_orbitals_per_atom))
    allocate(tmpid(1:atomix%species%nspecies))
    tmpid=0
    if (get_block(io,"Basis",nt)) then
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
        if (atomix%species%norbs(i) > gen%max_orbitals_per_atom) then
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
            atomix%species_basis(i,j)%sp = sp
            atomix%species_basis(i,j)%atom = 0
            atomix%species_basis(i,j)%spin=.false.
            atomix%species_basis(i,j)%n=n
            atomix%species_basis(i,j)%l=l
            atomix%species_basis(i,j)%m=m
            atomix%species_basis(i,j)%occup=qd
! spin up   
            atomix%species_basis(i,j+norbitals)%sp = sp
            atomix%species_basis(i,j+norbitals)%atom = 0
            atomix%species_basis(i,j+norbitals)%spin=.true.
            atomix%species_basis(i,j+norbitals)%n=atomix%species_basis(i,j)%n
            atomix%species_basis(i,j+norbitals)%l=atomix%species_basis(i,j)%l
            atomix%species_basis(i,j+norbitals)%m=atomix%species_basis(i,j)%m
            atomix%species_basis(i,j+norbitals)%occup=qu
            atomix%species%zval(i)=atomix%species%zval(i)+atomix%species_basis(i,j+norbitals)%occup&
                +atomix%species_basis(i,j)%occup
          else
            atomix%species_basis(i,j)%sp = sp
            atomix%species_basis(i,j)%atom = 0
            atomix%species_basis(i,j)%spin=.false.
            atomix%species_basis(i,j)%n=n
            atomix%species_basis(i,j)%l=l
            atomix%species_basis(i,j)%m=m
            atomix%species_basis(i,j)%occup=qd
            atomix%species%zval(i)=atomix%species%zval(i)+atomix%species_basis(i,j)%occup
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
            if (.not.atomix%species_basis(sp,j)%spin) then
              k = k + 1
              atomix%basis%orbitals(k) = atomix%species_basis(sp,j)
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
            if (atomix%species_basis(sp,j)%spin) then
              k = k + 1
              atomix%basis%orbitals(k) = atomix%species_basis(sp,j)
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
              atomix%basis%orbitals(k) = atomix%species_basis(sp,j)
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
end subroutine read_basis

!> \brief prints the basis set information
!> \author Alin M Elena
!> \date 30/10/07, 17:56:02
!> \param io type(io_type) i/o units
!> \param gen type(general_type) general data
!> \param atomix type(atomicx_type) contains info about atoms
  subroutine print_basis(io,gen,atomix)
    character(len=*), parameter :: sMyName="print_basis"
    type(io_type), intent(inout) :: io
    type(general_type), intent(in)  :: gen
    type(atomicx_type), intent(inout) :: atomix
    integer :: i,j,sp
    character(len=1) :: spin=""
    character(len=2) :: el=""
    write(io%uout,'(a)')"==BasisSetInfo==================================="
    write(io%uout,'(a)')"==Listed by species==================="
    write(io%uout,'(a)')"  #Orb  Sp El  N  L  M Spin Occupation"
    do i=1,atomix%species%nspecies
      sp=atomix%species%id(i)
      do j=1,atomix%species%norbs(i)
        if (gen%spin) then
          if (atomix%species_basis(i,j)%spin) then
            spin="U"
          else
            spin="D"
          endif
        endif
        write(io%uout,'(i6,i4,1x,a,3i3,4x,a,1x,f0.8)')j,atomix%species_basis(i,j)%sp,&
              symbol(atomix%species%z(atomix%species_basis(i,j)%sp)),&
            atomix%species_basis(i,j)%n,atomix%species_basis(i,j)%l,atomix%species_basis(i,j)%m,&
                spin,atomix%species_basis(i,j)%occup
      enddo
    enddo
    write(io%uout,'(a)')"==End Listed by species==============="
    write(io%uout,'(a)')"==Listed by atoms============================"
    write(io%uout,'(a)')"  #Orb  Atom  Sp  El  N  L  M Spin Occupation"
    spin=""
    el="el"
    do i=1,atomix%basis%norbitals
      if (atomix%basis%orbitals(i)%spin) then
        spin="U"
      else
        spin="D"
      endif
      el=symbol(get_Z(atomix,atomix%basis%orbitals(i)%atom))
      write(io%uout,'(2i6,i4,a4,3i3,1x,a4,1x,f0.8)')i,atomix%basis%orbitals(i)%atom,atomix%basis%orbitals(i)%sp,el,&
            atomix%basis%orbitals(i)%n,atomix%basis%orbitals(i)%l,atomix%basis%orbitals(i)%m,spin,&
            atomix%basis%orbitals(i)%occup
    enddo
    write(io%uout,'(a)')"==End Listed by atoms========================"
    write(io%uout,'(a)')"==EndBasisSetInfo================================"
  end subroutine print_basis


!> \brief selects the initialization routine for the tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:07:56
!> \param io type(io_type) contains all the info about I/O files
!> \param gen type(general_type) contains the info needed by the program to run
!> \param atomix type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters
  subroutine read_tb_model(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="read_tb_model"
    type(io_type), intent(inout) :: io
    type(general_type), intent(inout) :: gen
    type(atomicx_type), intent(inout) :: atomix
    type(model_type), intent(inout) :: tbMod

    integer :: i,j,k,k1,k2

    allocate(tbMod%hopping(1:atomix%species%nspecies,1:atomix%species%nspecies))

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        tbMod%hopping(i,j)%l1=atomix%species_basis(i,1)%l
        do k=2,atomix%species%norbs(i)
          if (tbMod%hopping(i,j)%l1 < atomix%species_basis(i,k)%l) then
            tbMod%hopping(i,j)%l1 = atomix%species_basis(i,k)%l
          endif
        enddo

        tbMod%hopping(i,j)%l2=atomix%species_basis(j,1)%l
        do k=2,atomix%species%norbs(j)
          if (tbMod%hopping(i,j)%l2<atomix%species_basis(j,k)%l) then
            tbMod%hopping(i,j)%l2=atomix%species_basis(j,k)%l
          endif
        enddo

        tbMod%hopping(i,j)%ll=min(tbMod%hopping(i,j)%l1,tbMod%hopping(i,j)%l2)
        allocate(tbMod%hopping(i,j)%a(0:tbMod%hopping(i,j)%l1,0:tbMod%hopping(i,j)%l2,0:tbMod%hopping(i,j)%ll))
        allocate(tbMod%hopping(i,j)%eps(0:atomix%species%norbs(i)-1))

      enddo
    enddo

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        tbMod%hopping(i,j)%eps=0.0_pr
        tbMod%hopping(i,j)%a1=0.0_pr
        tbMod%hopping(i,j)%a2=0.0_pr
        tbMod%hopping(i,j)%a3=0.0_pr
        tbMod%hopping(i,j)%a4=0.0_pr
        tbMod%hopping(i,j)%phi0=0.0_pr
        tbMod%hopping(i,j)%r0=0.0_pr
        tbMod%hopping(i,j)%rc=0.0_pr
        tbMod%hopping(i,j)%r1=0.0_pr
        tbMod%hopping(i,j)%rcut=0.0_pr
        tbMod%hopping(i,j)%n=0.0_pr
        tbMod%hopping(i,j)%nc=0.0_pr
        tbMod%hopping(i,j)%d0=0.0_pr
        tbMod%hopping(i,j)%dc=0.0_pr
        tbMod%hopping(i,j)%d1=0.0_pr
        tbMod%hopping(i,j)%dcut=0.0_pr
        tbMod%hopping(i,j)%m=0.0_pr
        tbMod%hopping(i,j)%mc=0.0_pr
        tbMod%hopping(i,j)%a=0.0_pr
      enddo
    enddo
    select case (gen%bond)
      case (bond_gsp)
        call read_tb_gsp(io,gen,atomix,tbMod)
        call print_tb_gsp(io,gen,atomix,tbMod)
      case (bond_harrison)
        call read_tb_harrison(io,gen,atomix,tbMod)
        call print_tb_harrison(io,gen,atomix,tbMod)
    end select  

  end subroutine read_tb_model


!> \brief selects the initialization routine for the Goodwin-Skinner-Petifor tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:15:15
!> \param io type(io_type) contains all the info about I/O files
!> \param gen type(general_type) contains the info needed by the program to run
!> \param atomix type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters
!> \remarks Using tight binding model as described in: \n
!> A.P. Horsfield,P.D. Godwin,D.G. Pettifor,A.P. Sutton \n
!> Computational materials synthesis. I. A tight-binding scheme for hydrocarbons - Phys. Rev. B 54, 22, 15 773

  subroutine read_tb_gsp(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="read_tb_gsp"
    type(io_type), intent(inout) :: io
    type(general_type), intent(inout) :: gen
    type(atomicx_type), intent(inout) :: atomix
    type(model_type), intent(inout) :: tbMod
    integer :: i,j,i1,k,k1,k2,p
    character(len=mw) :: read_var

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        if (i == j) then
          if (gen%spin) then
            i1=0
            p=0
            do while (i1<=atomix%species%norbs(i)/2-1)
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=get_real(io,trim(ccvar(i1,i1,read_var)),0.0_pr)
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
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=get_real(io,trim(ccvar(i1,i1,read_var)),0.0_pr)
              i1=i1+2*p+1
              p=p+1
            enddo
          endif
          if (gen%embedding) then
            read_var="a1"
            tbMod%hopping(i,j)%a1=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
            read_var="a2"
            tbMod%hopping(i,j)%a2=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
            read_var="a3"
            tbMod%hopping(i,j)%a3=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
            read_var="a4"
            tbMod%hopping(i,j)%a4=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
          endif
        endif
        read_var="phi0"
        tbMod%hopping(i,j)%phi0=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="r0"
        tbMod%hopping(i,j)%r0=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="rc"
        tbMod%hopping(i,j)%rc=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="r1"
        tbMod%hopping(i,j)%r1=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="rcut"
        tbMod%hopping(i,j)%rcut=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="n"
        tbMod%hopping(i,j)%n=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="nc"
        tbMod%hopping(i,j)%nc=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="d0"
        tbMod%hopping(i,j)%d0=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="dc"
        tbMod%hopping(i,j)%dc=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="d1"
        tbMod%hopping(i,j)%d1=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="dcut"
        tbMod%hopping(i,j)%dcut=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="m"
        tbMod%hopping(i,j)%m=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="mc"
        tbMod%hopping(i,j)%mc=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        do k=0,tbMod%hopping(i,j)%l1
          do k1=k,tbMod%hopping(i,j)%l2
            do k2=0,k
              tbMod%hopping(i,j)%a(k,k1,k2)=get_real(io,trim(ccnlm(i,j,k,k1,k2)),0.0_pr)
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
                tbMod%hopping(i,j)%a(k,k1,k2)=(-1.0_pr)**((k-k1+abs(k-k1))/2.0_pr)*tbMod%hopping(j,i)%a(k1,k,k2)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine read_tb_gsp



!> \brief prints the parameters of the gsp model that will be used 
!> \author Alin M Elena
!> \date 31/10/07, 16:55:21
!> \param io type(io_type) contains all the info about I/O files
!> \param gen type(general_type) contains the info needed by the program to run
!> \param atomix type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters
  subroutine print_tb_gsp(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="print_tb_gsp"
    type(io_type), intent(inout) :: io
    type(general_type), intent(inout) :: gen
    type(atomicx_type), intent(inout) :: atomix
    type(model_type), intent(inout) :: tbMod
    integer :: i,j,i1,k,k1,k2
    character(len=mw) :: read_var

    write(io%uout,'(a)')"==GSP Model======================================="
    write(io%uout,*) &
        "Using tight binding model as described in:"
    write(io%uout,*) &
        "A.P. Horsfield,P.D. Godwin,D.G. Pettifor,A.P. Sutton"
    write(io%uout,*) &
        "Computational materials synthesis. I. A tight-binding scheme"
    write(io%uout,*) &
        "for hydrocarbons - Phys. Rev. B 54, 22, 15 773"
    write(io%uout,'(a)')"===================================================="

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        write(io%uout,'(a,2i2)')"# Species",i,j
        if (i.eq.j) then
          do i1=0,atomix%species%norbs(i)-1
            read_var="eps"
            read_var=trim(ccvar(i,i,read_var))
            write(io%uout,*) trim(ccvar(i1,i1,read_var)),tbMod%hopping(i,j)%eps(i1)
          enddo
          if (gen%embedding) then
            read_var="a1"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a1
            read_var="a2"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a2
            read_var="a3"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a3
            read_var="a4"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a4
          endif
        endif
        read_var="phi0"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%phi0
        read_var="r0"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%r0
        read_var="rc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%rc
        read_var="r1"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%r1
        read_var="rcut"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%rcut
        read_var="n"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%n
        read_var="nc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),  tbMod%hopping(i,j)%nc
        read_var="d0"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%d0
        read_var="dc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%dc
        read_var="d1"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%d1
        read_var="dcut"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%dcut
        read_var="m"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),  tbMod%hopping(i,j)%m
        read_var="mc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%mc
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
              write(io%uout,*) &
              trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2)
            enddo
          enddo
        enddo
      enddo
    enddo
    write(io%uout,*) &
      "========================Finish========================="


end subroutine print_tb_gsp

 

!> \brief selects the initialization routine for the Harrison tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 15:16:15
!> \param io type(io_type) contains all the info about I/O files
!> \param gen type(general_type) contains the info needed by the program to run
!> \param atomix type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters

  subroutine read_tb_harrison(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="read_tb_harrison"
    type(io_type), intent(inout) :: io
    type(general_type), intent(inout) :: gen
    type(atomicx_type), intent(inout) :: atomix
    type(model_type), intent(inout) :: tbMod

    integer i,j,i1,k,k1,k2,p
    character(len=mw) :: read_var
    
    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        if (i==j) then
          if (gen%spin) then
            i1=0
            p=0
            do while (i1<=atomix%species%norbs(i)/2-1)
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=get_real(io,trim(ccvar(i1,i1,read_var)),0.0_pr)
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
              tbMod%hopping(i,j)%eps(i1:i1+2*p)=get_real(io,trim(ccvar(i1,i1,read_var)),0.0_pr)
              i1=i1+2*p+1
              p=p+1
            enddo
          endif
          if (gen%embedding) then
            read_var="a1"
            tbMod%hopping(i,j)%a1=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
            read_var="a2"
            tbMod%hopping(i,j)%a2=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
            read_var="a3"
            tbMod%hopping(i,j)%a3=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
            read_var="a4"
            tbMod%hopping(i,j)%a4=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
          endif
        endif

        read_var="phi0"
        tbMod%hopping(i,j)%phi0=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="r1"
        tbMod%hopping(i,j)%r1=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="rcut"
        tbMod%hopping(i,j)%rcut=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="n"
        tbMod%hopping(i,j)%n=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="d1"
        tbMod%hopping(i,j)%d1=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="dcut"
        tbMod%hopping(i,j)%dcut=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)
        read_var="m"
        tbMod%hopping(i,j)%m=get_real(io,trim(ccvar(i,j,read_var)),0.0_pr)

        do k=0,tbMod%hopping(i,j)%l1
          do k1=k,tbMod%hopping(i,j)%l2
            do k2=0,k
              tbMod%hopping(i,j)%a(k,k1,k2)=get_real(io,trim(ccnlm(i,j,k,k1,k2)),0.0_pr)
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
                tbMod%hopping(i,j)%a(k,k1,k2)=(-1.0_pr)**((k-k1+abs(k-k1))/2.0_pr)*tbMod%hopping(j,i)%a(k1,k,k2)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine read_tb_harrison

!> \brief prints the initialization routine for the Harrison tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 23:20:15
!> \param io type(io_type) contains all the info about I/O files
!> \param gen type(general_type) contains the info needed by the program to run
!> \param atomix type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters

  subroutine print_tb_harrison(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="print_tb_harrison"
    type(io_type), intent(inout) :: io
    type(general_type), intent(inout) :: gen
    type(atomicx_type), intent(inout) :: atomix
    type(model_type), intent(inout) :: tbMod

    integer i,j,i1,k,k1,k2,p
    character(len=mw) :: read_var

    write(io%uout,*)"==Harrison Parameters================================================"
    write(io%uout,*) &
        "Using Harrison type tight binding model "
      do i=1,atomix%species%nspecies
        do j=1,atomix%species%nspecies
          write(io%uout,'(a,2i2)')"#",i,j
          if (i==j) then
            do i1=0,atomix%species%norbs(i)-1
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              write(io%uout,*) trim(ccvar(i1,i1,read_var)),tbMod%hopping(i,j)%eps(i1)
            enddo
            if (gen%embedding) then
              read_var="a1"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a1
              read_var="a2"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a2
              read_var="a3"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a3
              read_var="a4"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a4
            endif
          endif
          read_var="phi0"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%phi0

          read_var="r1"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%r1
          read_var="rcut"
          write(io%uout,*) &
              trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%rcut
          read_var="n"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%n
          read_var="d1"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%d1
          read_var="dcut"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%dcut
          read_var="m"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%m
          do k=0,tbMod%hopping(i,j)%l1
            do k1=0,tbMod%hopping(i,j)%l2
              do k2=0,min(k,k1)
                write(io%uout,*) &
                trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2)
              enddo
            enddo
          enddo
        enddo
      enddo

      write(io%uout,*)"==End Harrison Parameters============================================"
  end subroutine print_tb_harrison
!> \brief reads the delta block
!> \details the block is read only if multipoles electrostatics is asked for
!> \author Alin M Elena
!> \date 31/10/07, 17:52:20
!> \param io type(io_type) contains all the info about I/O files
!> \param gen type(general_type) contains the info needed by the program to run
!> \param atomix type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters
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
  subroutine read_delta(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="read_delta"
    type(io_type), intent(inout) :: io
    type(general_type), intent(inout) :: gen
    type(atomicx_type), intent(inout) :: atomix
    type(model_type), intent(inout) :: tbMod
    integer :: nt,sp,j,i,k,l,errno
    character(len=ml) :: saux
    integer, allocatable :: tmpid(:)
        
    if (get_block(io,"DeltaPol",nt)) then
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
        j=maxval(atomix%species_basis(tbMod%delta(i)%sp,1:atomix%species%norbs(tbMod%delta(i)%sp)/2)%l)
        tbMod%delta(i)%l=j
        allocate(tbMod%delta(i)%d(0:2*j,0:j,0:j))
        tbMod%delta(i)%d=0.0_pr
        do j=0,tbMod%delta(i)%l
          do k=j,tbMod%delta(i)%l
            do l=abs(j-k),j+k
              if (mod(k+j+l,2)==0) then
                if ((k==j).and.(l==0)) then
                  tbMod%delta(i)%d(0,k,j)=1.0_pr
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
      write(io%uout,*) &
        '==DeltaPol Block====================================================='
      write(io%uout,*) "Note that: delta (l,l1,l2) with |l1-l2|<=l<=|l1+l2|"
      do i = 1,atomix%species%nspecies
        write(io%uout,'(a,i2)') "Species ",tbMod%delta(i)%sp
        do j=0,tbMod%delta(i)%l
          do k=j,tbMod%delta(i)%l
            if (k==j) then
              write(io%uout,'(100(a,3i3,a,f0.8))')(' delta(',l,j,k,') = ', tbMod%delta(i)%d(l,j,k),l=abs(j-k),abs(j+k))
            else
              write(io%uout,'(100(a,3i3,a,3i3,a,f14.8))')(' delta(',l,j,k,') = delta(',l,k,j,') = ', &
                    tbMod%delta(i)%d(l,j,k),l=abs(j-k),abs(j+k))
            endif
          enddo
        enddo
      enddo
      write(io%uout,*) &
        '==End DeltaPol Block================================================='
      else
        call error("DeltaPol block is missing",sMyName,.true.,io)
      endif
  end subroutine read_delta

!> \brief deallocates the memory
!> \author Alin M Elena
!> \date 29/10/07, 17:38:16
!> \param atomic type(atomicx_type) contains atomic data
!> \param general type(general_type) general data that controls the flow of the program
!> \param tbMod type(model_type) contains iformation about the tight binding model parameters

  subroutine clean_memory(atomic,general,tbMod)
    character(len=*), parameter :: sMyName="clean_memory"
    type(atomicx_type), intent(inout) :: atomic
    type(general_type), intent(in) :: general
    type(model_type),intent(in) :: tbMod

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
    deallocate(atomic%species_basis)
    deallocate(atomic%basis%orbitals)
    deallocate(atomic%atoms%orbs)

    do i=1,atomic%species%nspecies
      do j=1,atomic%species%nspecies
        deallocate(tbMod%hopping(i,j)%a)
        deallocate(tbMod%hopping(i,j)%eps)
      enddo
    enddo
    deallocate(tbMod%hopping)
    if(general%electrostatics==electrostatics_multipoles) then
      do i=1,atomic%species%nspecies
        deallocate(tbMod%delta(i)%d)
      enddo
      deallocate(tbMod%delta)
    endif
  
    !internal variables
  end subroutine clean_memory
end module read_data