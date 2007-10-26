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


  implicit none
  private
  public :: initialize,close_io_general

contains


!> \brief reads from io_loc different fields of gen_loc
!> \author Alin M Elena
!> \date ~2007
!> \param io_loc type(io_type) contains all the info about I/O files
!> \param gen_loc type(general_type) contains the info needed by the program to run

  subroutine initialize(io_loc,gen_loc)
    character(len=*), parameter :: myname = 'initialize'
    type(io_type), intent(inout) :: io_loc
    type(general_type), intent(inout) :: gen_loc		
    real(pr) :: aux
    integer :: errno,nt,i
    real(pr) :: a1,a2,a3


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

call read_general(io_loc,gen_loc)

if (get_block(io_loc,"coords",nt)) then
  write(io_loc%udeb,'(a)')"block coords"
  do i=1,2
    read(nt,fmt=*,iostat=errno)a1,a2,a3
    if (errno/=0) call error("block coords is not in the right format"&
      ,myname,.true.,io_loc)
    write(io_loc%udeb,'(3f16.6)')a1,a2,a3
  enddo
  write(io_loc%udeb,'(a)')"endblock coords"
endif
a1=get_real(io_loc,"gf1",0.0_pr)

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
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off printing on the screen, on means that nothing will be written in the output file</TD>
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
!comm_io DebugLevel & 5,15,25 & 5 & debug information level (low, medium,high) \\
!comm_io OutputLevel & 5,15,25 & 5 & output information level (low, medium,high) \\
  io_loc%verbosity=get_integer(io_loc,"OutputLevel",5)   
  io_loc%debug=get_integer(io_loc,"DebugLevel",5)
  io_loc%first_time=.not. io_loc%first_time

  if (io_loc%verbosity>=high_debug) then
    write( io_loc%udeb,'(a,l)')"Read io_info for the first time: "&
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
  write( io_loc%uout,'(a,l)')"Standard Output(OnScreen): "&
    ,io_loc%stdout


end subroutine read_io


      
      
!> \brief reads the parameters necessary for running
!> \details
!> \author Alin M Elena
!> \date 21st of January, 2007 
!> \param  io_loc type(io_type) I/O details (see types::io_type)
!> \param gen_loc type(general_type), keeps all the general information about the parameters of the program
!> \remarks
!> A short description of the input file general variables part
!> \latexonly
!>  \begin{longtable}{|c||p{0.2\textwidth}|c||p{0.25\textwidth}|}
!> \hline
!> \textbf{VarName} & \textbf{Values} & \textbf{Default Value} & \textbf{Description} \\\
!> \hline \hline
!> JobName & string & no name & a name for the job \\ 
!> \hline
!> RanSeed & integer & 12345 & A seed for the random number generator \\ 
!> \hline
!> WriteAni & logical & .false. & on/off writing animation file \\ 
!> \hline
!> WriteEne & logical & .true. & on/off writing energy file \\ 
!> \hline
!> WriteDen & logical & .false. & on/off writing density file \\ 
!> \hline
!> ReadDen & logical & .false. & on/off reading density file \\ 
!> \hline
!> ReadVel & logical & .false. & on/off reading velocity block \\ 
!> \hline
!> IonicTemperature & real & 300.0 & Ionic temperature\\ 
!> \hline
!> NetCharge & real & 0.0 & Net charge on the system\\ 
!> \hline
!> ElectronicTemperature & real & 300.0 & Electronic temperature, used to compute occupation numbers if you choose Fermi-Dirac or Methfessel-Paxton method\\ 
!> \hline
!> ElectronicMu & real & 0.0 & chemical potential, used to compute occupation numbers if you choose constant $\mu$ method\\ 
!> \hline
!> DeltaT & real & 0.001 & time step used to evolve equtions of motion\\ 
!> \hline
!> Nsteps & integer & 100 & number of steps used to evolve equtions of motion\\ 
!> \hline
!> RunType & SinglePoint, BODynamics, Ehrenfest, EhrenfestDamped, Fit, ForceTest, ForceTestX, ForceTestY, ForceTestZ & SinglePoint & Type of calculation\\ 
!> \hline
!> SCF & logical & .false. & on/off self consistent field method \\ 
!> \hline
!> SCFType & TB+UJ & TB+UJ & SCF method\\ 
!> \hline
!> SCFSteps & integer & 100 & maximum number of steps used for scf\\ 
!> \hline
!> SCFMix & real & 0.85 & mixing parameter\\ 
!> \hline
!> SCFMixType & Broyden, Pulay & Broyden & SCF mixing method\\ 
!> \hline
!> SCFTol & real & 1e-8 & convergence tolerance\\ 
!> \hline
!> SCFMixN & integer & 4 & number of iterations to mix\\ 
!> \hline
!> VelScale & logical & .false. & on/off scaling velocities\\ 
!> \hline
!> DMOccTol & real & 1e-10 & density matrix occupation tolerance\\ 
!> \hline
!> HElThres & real & 1e-10 & hamiltionian element thresold. Any element smaller that the thresold is made zero.\\ 
!> \hline
!> Spin & logical & .false. & on/off spin polarisation\\ 
!> \hline
!> CollinearSpins & logical & .false. & on/off collinear spins\\ 
!> \hline
!> SpinDU & real & 0.0 & spin down spin up difference\\ 
!> \hline
!> EulerSteps & integer & 100 & after each EulerSteps apply an Euler integration of equations of motions\\ 
!> \hline
!> Electrostatics & PointCharges, Multipoles & PointCharges & method use to compute electrostatic interaction\\ 
!> \hline
!> PrecomputeMultipoles & logical & .true. & on/off precompute multipoles\\ 
!> \hline
!> Embedding & logical & .true. & on/off embedding method\\ 
!> \hline
!> FSteps & integer & 100 & Number of steps used to test the force/energy consitency\\ 
!> \hline
!> FStart & real & 0.0 & position at which the force/energy consistency test starts\\ 
!> \hline
!> Fdx & real & 0.001 & space step for the force/energy consistency test starts\\ 
!> \hline
!> Gamma & real & 0.3 & damping factor for Ehrenfest equation\\ 
!> \hline
!> MPN & integer & 2 & the order used for Methfessel-Paxton method\\ 
!> \hline
!> Hole & integer & 0 & no of level to create hole\\ 
!> \hline
!> Excite & integer & 0 & no of level to create excitation\\ 
!> \hline
!> HoleSpin & D,U & D & spin of the hole\\ 
!> \hline
!> ExciteSpin & D,U & D & spin of the excitation\\ 
!> \hline
!> Units & AU EV SI & AU & system of units atomic units, electronVolt-Angstrom, International\\ 
!> \hline
!> SmearingMethod & FD,MP,CMU & FD & smearing method (Femi-Dirac, Methfessel-Paxton, constant $\mu$)\\ 
!> \hline
!> BondType & Harrison, GSP &Harrison& bond type\\ 
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
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer </TD>
!> <TD ALIGN="CENTER">12345</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>A seed for the random number generator</TD>
!> </TR>
!> <TR> <TD ALIGN="CENTER">WriteAni </TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical </TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off writing animation file </TD>
!> </TR>
!> <TR> <TD ALIGN="CENTER">WriteEne </TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.true. </TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off writing energy file </TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">WriteDen </TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false. </TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off writing density file</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">ReadDen</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off reading density file</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">ReadVel</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off reading velocity block</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">IonicTemperature</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">300.0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Ionic temperature</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">NetCharge</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Net charge on the system</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">ElectronicTemperature</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">300.0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Electronic temperature, used to compute occupation numbers if you choose Fermi-Dirac or Methfessel-Paxton method</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">ElectronicMu</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>chemical potential, used to compute occupation numbers if you choose constant  method</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">DeltaT</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.001</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>time step used to evolve equtions of motion</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Nsteps</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">100</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>number of steps used to evolve equtions of motion</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">RunType</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>SinglePoint, BODynamics, Ehrenfest, EhrenfestDamped, Fit, ForceTest, ForceTestX, ForceTestY, ForceTestZ</TD>
!> <TD ALIGN="CENTER">SinglePoint</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Type of calculation</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SCF</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off self consistent field method</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SCFType</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>TB+UJ</TD>
!> <TD ALIGN="CENTER">TB+UJ</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>SCF method</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SCFSteps</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">100</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>maximum number of steps used for scf</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SCFMix</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.85</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>mixing parameter</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SCFMixType</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>Broyden, Pulay</TD>
!> <TD ALIGN="CENTER">Broyden</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>SCF mixing method</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SCFTol</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">1e-8</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>convergence tolerance</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SCFMixN</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">4</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>number of iterations to mix</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">VelScale</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off scaling velocities</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">DMOccTol</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">1e-10</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>density matrix occupation tolerance</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">HElThres</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">1e-10</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>hamiltionian element thresold. Any element smaller that the thresold is made zero.</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Spin</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off spin polarisation</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">CollinearSpins</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.false.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off collinear spins</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SpinDU</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>spin down spin up difference</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">EulerSteps</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">100</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>after each EulerSteps apply an Euler integration of equations of motions</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Electrostatics</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>PointCharges, Multipoles</TD>
!> <TD ALIGN="CENTER">PointCharges</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>method use to compute electrostatic interaction</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">PrecomputeMultipoles</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.true.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off precompute multipoles</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Embedding</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>logical</TD>
!> <TD ALIGN="CENTER">.true.</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>on/off embedding method</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">FSteps</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">100</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>Number of steps used to test the force/energy consitency</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">FStart</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>position at which the force/energy consistency test starts</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Fdx</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.001</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>space step for the force/energy consistency test starts</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Gamma</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>real</TD>
!> <TD ALIGN="CENTER">0.3</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>damping factor for Ehrenfest equation</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">MPN</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">2</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>the order used for Methfessel-Paxton method</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Hole</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>no of level to create hole</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Excite</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>integer</TD>
!> <TD ALIGN="CENTER">0</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>no of level to create excitation</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">HoleSpin</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>D,U</TD>
!> <TD ALIGN="CENTER">D</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>spin of the hole</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">ExciteSpin</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>D,U</TD>
!> <TD ALIGN="CENTER">D</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>spin of the excitation</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">Units</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>AU EV SI</TD>
!> <TD ALIGN="CENTER">AU</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>system of units atomic units, electronVolt-Angstrom, International</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">SmearingMethod</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>FD,MP,CMU</TD>
!> <TD ALIGN="CENTER">FD</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>smearing method (Femi-Dirac, Methfessel-Paxton, constant )</TD>
!> </TR>
!> <TR><TD ALIGN="CENTER">BondType</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=100>Harrison, GSP</TD>
!> <TD ALIGN="CENTER">Harrison</TD>
!> <TD ALIGN="LEFT" VALIGN="TOP" WIDTH=125>bond type</TD>
!> </TR>
!> </TABLE>




!> \endhtmlonly

subroutine read_general(io_loc,gen_loc)
character(len=*),parameter :: name="read_general"
  type(io_type), intent(inout) :: io_loc
  type(general_type), intent(inout) :: gen_loc

character(len=mw) :: saux

      if (.not. gen_loc%first_time) gen_loc%first_time=.not. gen_loc%first_time

!comm_gen JobName & string & no name & a name for the job \\
  gen_loc%job_name=get_string(io_loc,"JobName","no name")
  write( io_loc%uout,'(a,a)')"Job Name(JobName): "&
    ,gen_loc%job_name

!comm_gen RanSeed & integer & 12345 & A seed for the random number generator \\
  gen_loc%ran3_seed=get_integer(io_loc,"RanSeed",12345)
  write( io_loc%uout,'(a,i0)')"Random number generator seed(RanSeed): "&
    ,gen_loc%ran3_seed

!comm_gen WriteAni & logical & .false. & on/off writing animation file \\
  gen_loc%write_ani=get_logical(io_loc,"WriteAni",.false.)
  write( io_loc%uout,'(a,l)')"Write animation(WriteAni): "&
    ,gen_loc%write_ani

!comm_gen WriteEne & logical & .true. & on/off writing energy file \\
  gen_loc%write_ene=get_logical(io_loc,"WriteEne",.true.)
  write( io_loc%uout,'(a,l)')"Write Energy(WriteEne): "&
    ,gen_loc%write_ene

!comm_gen WriteDen & logical & .false. & on/off writing density file \\
  gen_loc%write_density=get_logical(io_loc,"WriteDen",.false.)
  write( io_loc%uout,'(a,l)')"Write density(WriteDen): "&
    ,gen_loc%write_density
!comm_gen ReadDen & logical & .false. & on/off reading density file \\
  gen_loc%read_density=get_logical(io_loc,"ReadDen",.false.)
  write( io_loc%uout,'(a,l)')"Read Density(ReadDen): "&
    ,gen_loc%read_density

!comm_gen ReadVel & logical & .false. & on/off reading velocity block \\
  gen_loc%read_velocity=get_logical(io_loc,"ReadVel",.false.)
  write( io_loc%uout,'(a,l)')"Read Velocity(ReadVel): "&
    ,gen_loc%read_velocity

!comm_gen IonicTemperature & real & 300.0 & Ionic temperature\\
  gen_loc%ionic_temperature=get_real(io_loc,"IonicTemperature",300.0_pr)
  write( io_loc%uout,'(a,g)')"Ionic Temperature(IonicTemperature): "&
    ,gen_loc%ionic_temperature

!comm_gen NetCharge & real & 0.0 & Net charge on the system\\
  gen_loc%netcharge=get_real(io_loc,"NetCharge",0.0_pr)
  write( io_loc%uout,'(a,g)')"Net charge(NetCharge): "&
    ,gen_loc%netcharge

!comm_gen ElectronicTemperature & real & 300.0 & Electronic temperature, used to compute occupation numbers if you choose Fermi-Dirac or Methfessel-Paxton method\\
  gen_loc%electronic_temperature=get_real(io_loc,"ElectronicTemperature",300.0_pr)
  write( io_loc%uout,'(a,g)')"Electronic Temperature(ElectronicTemperature): "&
    ,gen_loc%electronic_temperature

!comm_gen ElectronicMu & real & 0.0 & chemical potential, used to compute occupation numbers if you choose constant $\mu$ method\\
  gen_loc%electronic_mu=get_real(io_loc,"ElectronicMu",0.0_pr)
  write( io_loc%uout,'(a,g)')"Chemical Potential(ElectronicMu): "&
    ,gen_loc%electronic_mu

!comm_gen DeltaT & real & 0.001 & time step used to evolve equtions of motion\\
  gen_loc%deltat=get_real(io_loc,"DeltaT",0.001_pr)
  write( io_loc%uout,'(a,g)')"Time step(DeltaT): "&
    ,gen_loc%deltat

!comm_gen Nsteps & integer & 100 & number of steps used to evolve equtions of motion\\
  gen_loc%nsteps=get_integer(io_loc,"Nsteps",100)
  write( io_loc%uout,'(a,i0)')"Number of steps(Nsteps): "&
    ,gen_loc%nsteps

!comm_gen RunType & SinglePoint, BODynamics, Ehrenfest, EhrenfestDamped, Fit, ForceTest, ForceTestX, ForceTestY, ForceTestZ & SinglePoint & Type of calculation\\
  saux=get_string(io_loc,"RunType","SinglePoint")
  write( io_loc%uout,'(a,a)')"Calculation type(RunType): "&
    ,trim(saux)
  if (cstr(trim(saux),'SinglePoint')) then
    gen_loc%runtype = run_sp
  elseif (cstr(trim(saux),'bodynamics')) then
    gen_loc%runtype = run_bo
  elseif (cstr(trim(saux),'ehrenfest')) then
    gen_loc%runtype = run_ehrenfest
  elseif (cstr(trim(saux),'fit')) then
    gen_loc%runtype = run_fit
  elseif (cstr(trim(saux),'forcetest')) then
    gen_loc%runtype = run_force_test
  elseif (cstr(trim(saux),'forcetestx')) then
    gen_loc%runtype = run_force_testx
  elseif (cstr(trim(saux),'forcetesty')) then
    gen_loc%runtype = run_force_testy
  elseif (cstr(trim(saux),'forcetestz')) then
    gen_loc%runtype = run_force_testz
  elseif (cstr(trim(saux),'ehrenfestdamped')) then
    gen_loc%runtype = run_ehrenfest_damped
  else 
    call error("The requested RunType is not implemented",name,.true.,io_loc)
  endif

!comm_gen SCF & logical & .false. & on/off self consistent field method \\
  gen_loc%scf=get_logical(io_loc,"SCF",.false.)
  write( io_loc%uout,'(a,l)')"SCF?(SCF): "&
    ,gen_loc%scf

!comm_gen SCFType & TB+UJ & TB+UJ & SCF method\\
  saux=get_string(io_loc,"SCFType","TB+UJ")
  write( io_loc%uout,'(a,a)')"SCF type(SCFType): "&
    ,trim(saux)
  if (cstr(trim(saux),'TB+UJ')) then
    gen_loc%scf_type = scf_tbuj
  else 
    call error("The requested SCFType is not implemented",name,.true.,io_loc)
  endif

!comm_gen SCFSteps & integer & 100 & maximum number of steps used for scf\\
  gen_loc%maxscf=get_integer(io_loc,"SCFSteps",100)
  write( io_loc%uout,'(a,i0)')"SCF Maximum number of steps(SCFSteps): "&
    ,gen_loc%maxscf
!comm_gen SCFMix & real & 0.85 & mixing parameter\\
  gen_loc%scfmix=get_real(io_loc,"SCFMix",0.85_pr)
  write( io_loc%uout,'(a,g)')"SCF Mixing parameter(SCFMix): "&
    ,gen_loc%scfmix

!comm_gen SCFMixType & Broyden, Pulay & Broyden & SCF mixing method\\
  saux=get_string(io_loc,"SCFMixType","Broyden")
  write( io_loc%uout,'(a,a)')"SCF Mixing type(SCFMixType): "&
    ,trim(saux)
  if (cstr(trim(saux),'Broyden')) then
    gen_loc%scfmixtype = scfmix_broyden
  elseif (cstr(trim(saux),'Pulay')) then
    gen_loc%scfmixtype = scfmix_pulay
  else  
    call error("The requested SCFMixType is not implemented",name,.true.,io_loc)
  endif

!comm_gen SCFTol & real & 1e-8 & convergence tolerance\\
  gen_loc%scftol=get_real(io_loc,"SCFTol",1.0e-8_pr)
  write( io_loc%uout,'(a,g)')"SCF convergence tolerance(SCFTol): "&
    ,gen_loc%scftol

!comm_gen SCFMixN & integer & 4 & number of iterations to mix\\
  gen_loc%scfmixn=get_integer(io_loc,"SCFMixN",4)
  write( io_loc%uout,'(a,i0)')"number of iterations to mix(SCFMixN): "&
    ,gen_loc%scfmixn

!comm_gen VelScale & logical & .false. & on/off scaling velocities\\
  gen_loc%velocity_scaling=get_logical(io_loc,"VelScale",.false.)
  write( io_loc%uout,'(a,l)')"scaling  Velocities(VelScale): "&
    ,gen_loc%velocity_scaling
!comm_gen DMOccTol & real & 1e-10 & density matrix occupation tolerance\\
  gen_loc%dm_occupation_tolerance=get_real(io_loc,"DMOccTol",1.0e-10_pr)
  write( io_loc%uout,'(a,g)')"density matrix occupation tolerance(DMOccTol): "&
    ,gen_loc%dm_occupation_tolerance
!comm_gen HElThres & real & 1e-10 & hamiltionian element thresold. Any element smaller that the thresold is made zero.\\
  gen_loc%h_element_threshold=get_real(io_loc,"DMOccTol",1.0e-10_pr)
  write( io_loc%uout,'(a,g)')"hamiltionian element thresold(HElThres): "&
    ,gen_loc%h_element_threshold

!comm_gen Spin & logical & .false. & on/off spin polarisation\\
  gen_loc%spin=get_logical(io_loc,"Spin",.false.)
  write( io_loc%uout,'(a,l)')"spin polarisation(Spin): "&
    ,gen_loc%spin

!comm_gen CollinearSpins & logical & .false. & on/off collinear spins\\
  gen_loc%collinear=get_logical(io_loc,"CollinearSpins",.false.)
  write( io_loc%uout,'(a,l)')"collinear spins(CollinearSpins): "&
    ,gen_loc%collinear

!comm_gen SpinDU & real & 0.0 & spin down spin up difference\\
  gen_loc%sdu=get_real(io_loc,"SpinDU",0.0_pr)
  write( io_loc%uout,'(a,g)')"spin down spin up difference(SpinDU): "&
    ,gen_loc%sdu

!comm_gen EulerSteps & integer & 100 & after each EulerSteps apply an Euler integration of equations of motions\\
  gen_loc%euler_Steps=get_integer(io_loc,"EulerSteps",100)
  write( io_loc%uout,'(a,i0)')"Euler steps (EulerSteps): "&
    ,gen_loc%euler_steps

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

!comm_gen PrecomputeMultipoles & logical & .true. & on/off precompute multipoles\\
  gen_loc%comp_elec=get_logical(io_loc,"PrecomputeMultipoles",.true.)
  write( io_loc%uout,'(a,l)')"precompute multipoles(PrecomputeMultipoles): "&
    ,gen_loc%comp_elec

!comm_gen Embedding & logical & .true. & on/off embedding method\\
  gen_loc%embedding=get_logical(io_loc,"Embedding",.true.)
  write( io_loc%uout,'(a,l)')"has embedding(Embedding): "&
    ,gen_loc%embedding

!comm_gen FSteps & integer & 100 & Number of steps used to test the force/energy consitency\\
  gen_loc%f_steps=get_integer(io_loc,"FSteps",100)
  write( io_loc%uout,'(a,i0)')"force steps (FSteps): "&
    ,gen_loc%f_steps

!comm_gen FStart & real & 0.0 & position at which the force/energy consistency test starts\\
  gen_loc%f_start=get_real(io_loc,"FStart",0.0_pr)
  write( io_loc%uout,'(a,g)')"initial position for force(FStart): "&
    ,gen_loc%f_start

!comm_gen Fdx & real & 0.001 & space step for the force/energy consistency test starts\\
  gen_loc%f_dx=get_real(io_loc,"FStart",0.001_pr)
  write( io_loc%uout,'(a,g)')"space step for force(Fdx): "&
    ,gen_loc%f_dx

!comm_gen Gamma & real & 0.3 & damping factor for Ehrenfest equation\\
  gen_loc%gamma=get_real(io_loc,"Gamma",0.5_pr)
  write( io_loc%uout,'(a,g)')"damping factor in Ehrenfest Equation(Gamma): "&
    ,gen_loc%gamma

!comm_gen MPN & integer & 2 & the order used for Methfessel-Paxton method\\
  gen_loc%mp_N=get_integer(io_loc,"MPN",2)
  write( io_loc%uout,'(a,i0)')"Order for Methfessel-Paxton method(MPN): "&
    ,gen_loc%mp_N

!comm_gen Hole & integer & 0 & no of level to create hole\\
  gen_loc%hole=get_integer(io_loc,"Hole",0)
  write( io_loc%uout,'(a,i0)')"Hole level(Hole): "&
    ,gen_loc%hole

!comm_gen Excite & integer & 0 & no of level to create excitation\\
  gen_loc%excite=get_integer(io_loc,"Excite",0)
  write( io_loc%uout,'(a,i0)')"Excitation level(Excite): "&
    ,gen_loc%excite
!comm_gen HoleSpin & D,U & D & spin of the hole\\
  saux=get_string(io_loc,"HoleSpin","D")
  write( io_loc%uout,'(a,a)')"Spin of the Hole (HoleSpin): "&
    ,trim(saux)
  if (cstr(trim(saux),'D')) then
    gen_loc%hole_spin = spin_down
  elseif (cstr(trim(saux),'U')) then
    gen_loc%hole_spin = spin_up
  else  
    call error("The requested spin type is not implemented",name,.true.,io_loc)
  endif


!comm_gen ExciteSpin & D,U & D & spin of the excitation\\
  saux=get_string(io_loc,"ExciteSpin","D")
  write( io_loc%uout,'(a,a)')"Spin of the excitation (ExciteSpin): "&
    ,trim(saux)
  if (cstr(trim(saux),'D')) then
    gen_loc%excite_spin = spin_down
  elseif (cstr(trim(saux),'U')) then
    gen_loc%excite_spin = spin_up
  else  
    call error("The requested spin type is not implemented",name,.true.,io_loc)
  endif


!! to be moved in read_atoms
! ! ! !comm_atom NAcceptor & integer & 0 & number of atoms in Acceptor group\\
! ! ! gen_loc%nacceptor=get_integer(io_loc,"NAcceptor",0)
! ! !  write( io_loc%uout,'(a,i0)')"Atoms in Acceptor Group(NAcceptor): "&
! ! !         ,gen_loc%naccpetor
! ! ! 
! ! ! !comm_atom NDonor & integer & 0 & number of atoms in Donor group\\
! ! ! gen_loc%nacceptor=get_integer(io_loc,"NAcceptor",0)
! ! !  write( io_loc%uout,'(a,i0)')"Atoms in Acceptor Group(NAcceptor): "&
! ! !         ,gen_loc%naccpetor


!comm_gen Units & AU EV SI & AU & system of units atomic units, electronVolt-Angstrom, International\\
  saux=get_string(io_loc,"Units","AU")
  write( io_loc%uout,'(a,a)')"System of Units (Units): "&
    ,trim(saux)
  if (cstr(trim(saux),'Au')) then
    gen_loc%units = units_au
  elseif (cstr(trim(saux),'EV')) then
    gen_loc%units = units_ev
  elseif (cstr(trim(saux),'SI')) then
    gen_loc%units = units_si
  else  
    call error("The requested system of units is not implemented",name,.true.,io_loc)
  endif

!comm_gen SmearingMethod & FD,MP,CMU & FD & smearing method (Femi-Dirac, Methfessel-Paxton, constant $\mu$)\\
  saux=get_string(io_loc,"SmearingMethod","FD")
  write( io_loc%uout,'(a,a)')"Smearing Method (SmearingMethod): "&
    ,trim(saux)
  if (cstr(trim(saux),'FD')) then
    gen_loc%smethod = sm_fd
  elseif (cstr(trim(saux),'MP')) then
    gen_loc%smethod = sm_mp
  elseif (cstr(trim(saux),'CMU')) then
    gen_loc%smethod = sm_cmu
  else  
    call error("The requested smearing method is not implemented",name,.true.,io_loc)
  endif

!comm_gen BondType & Harrison, GSP &Harrison& bond type\\
  saux=get_string(io_loc,"BondType","Harrison")
  write( io_loc%uout,'(a,a)')"bond type (BondType): "&
    ,trim(saux)
  if (cstr(trim(saux),'Harrison')) then
    gen_loc%bond = bond_harrison
  elseif (cstr(trim(saux),'MP')) then
    gen_loc%bond = bond_gsp
  else  
    call error("The requested bond type is not implemented",name,.true.,io_loc)
  endif


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

end module read_data
