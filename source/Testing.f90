module m_Testing
  use m_Constants
  use m_Types
  use m_Useful
  use m_DriverRoutines
  implicit none
  private


  public :: ForceTest
  public :: forceTestx
  public :: forceTesty
  public :: forceTestz

contains
!> \brief tests force on x
!> \details  moves atom on the x axis from start in steps steps by dx each step
!>   in order to test that the energy and force are consistent
!>    generates two outputs:
!>    fort.777 data file with has 3 columns (x,force,-total energy)
!>   can be visualised with xmgrace (xmgrace -nxy fort.777)
!>   writes in the spciefied file the maximum error done in (%)
!> author Alin M. Elena (Belfast)
!> \date 22nd of June, 2006 generalized
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine forceTestx(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'forceTestx'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb

    integer :: i,atom
    real(k_pr) :: dx,tot
    real(k_pr), allocatable :: measured(:),analytic(:)

    if((gen%fatom <=0) .or.(gen%fatom > atomic%atoms%natoms)) then
      call error(myname,"can not calculate force on non existing atom",.true.,io)
    endif

    allocate(measured(1:gen%fsteps),analytic(1:gen%fsteps))
    dx=gen%fdx
    atom = gen%fatom
    do i=1,gen%fsteps
      atomic%atoms%x(atom)=gen%fstart+(i-1)*dx
      call SinglePoint(io,gen,atomic,tb,sol)
      tot=sol%totalEnergy
      write(777,*)atomic%atoms%x(atom),atomic%atoms%fx(atom),-tot
      measured(i)=-tot
      analytic(i)=atomic%atoms%fx(atom)
    end do
     write(io%uout,'(a,f12.8,a)')&
       'Maximum error',norm(measured,analytic,dx),'%'
     deallocate(measured,analytic)
   end subroutine forceTestx

!> \brief tests force on y
!> \details  moves atom on the y axis from start in steps steps by dx each step
!>   in order to test that the energy and force are consistent
!>    generates two outputs:
!>    fort.777 data file with has 3 columns (x,force,-total energy)
!>   can be visualised with xmgrace (xmgrace -nxy fort.777)
!>   writes in the spciefied file the maximum error done in (%)
!> author Alin M. Elena (Belfast)
!> \date 22nd of June, 2006 generalized
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine forceTesty(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'forceTesty'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb

    integer :: i,atom
    real(k_pr) :: dx,tot
    real(k_pr), allocatable :: measured(:),analytic(:)

    if((gen%fatom <=0) .or.(gen%fatom > atomic%atoms%natoms)) then
      call error(myname,"can not calculate force on non existing atom",.true.,io)
    endif

    allocate(measured(1:gen%fsteps),analytic(1:gen%fsteps))
    dx=gen%fdx
    atom = gen%fatom
    do i=1,gen%fsteps
      atomic%atoms%y(atom)=gen%fstart+(i-1)*dx
      call SinglePoint(io,gen,atomic,tb,sol)
      tot=sol%totalEnergy
      write(777,*)atomic%atoms%y(atom),atomic%atoms%fy(atom),-tot
      measured(i)=-tot
      analytic(i)=atomic%atoms%fy(atom)
    end do
     write(io%uout,'(a,f12.8,a)')&
       'Maximum error',norm(measured,analytic,dx),'%'
     deallocate(measured,analytic)
   end subroutine forceTesty

!> \brief tests force on z
!> \details  moves atom on the z axis from start in steps steps by dx each step
!>   in order to test that the energy and force are consistent
!>    generates two outputs:
!>    fort.777 data file with has 3 columns (x,force,-total energy)
!>   can be visualised with xmgrace (xmgrace -nxy fort.777)
!>   writes in the spciefied file the maximum error done in (%)
!> author Alin M. Elena (Belfast)
!> \date 22nd of June, 2006 generalized
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine forceTestz(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'forceTestx'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb

    integer :: i,atom
    real(k_pr) :: dx,tot
    real(k_pr), allocatable :: measured(:),analytic(:)

    if((gen%fatom <=0) .or.(gen%fatom > atomic%atoms%natoms)) then
      call error(myname,"can not calculate force on non existing atom",.true.,io)
    endif

    allocate(measured(1:gen%fsteps),analytic(1:gen%fsteps))
    dx=gen%fdx
    atom = gen%fatom
    do i=1,gen%fsteps
      atomic%atoms%z(atom)=gen%fstart+(i-1)*dx
      call SinglePoint(io,gen,atomic,tb,sol)
      tot=sol%totalEnergy
      write(777,*)atomic%atoms%z(atom),atomic%atoms%fz(atom),-tot
      measured(i)=-tot
      analytic(i)=atomic%atoms%fz(atom)
    end do
     write(io%uout,'(a,f12.8,a)')&
       'Maximum error',norm(measured,analytic,dx),'%'
     deallocate(measured,analytic)
   end subroutine forceTestz

!> \brief  computes the analytic and the numeric force on atoms
!> \details  for the numeric force shifts the x,y,z coordinates by a small
!>  quantity dx. In the end display a useful summary.
!> \author Alin M. Elena (Belfast)
!> \date 1st of April, 2006
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine ForceTest(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'ForceTest'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    integer       :: i
    real(k_pr)      :: fa_x, fa_y, fa_z, fn_x, fn_y, fn_z
    real(k_pr)      :: eenergy, renergy, scfenergy, minusts
    real(k_pr)      :: ep, em, dx,totnx,totny,totnz
    !-------------------------------------------------!

    totnx=0.0_k_pr
    totny=0.0_k_pr
    totnz=0.0_k_pr
    dx = gen%fdx
    do i=1,atomic%atoms%natoms
       ! Analytic forces
      call SinglePoint(io,gen,atomic,tb,sol)
      fa_x = atomic%atoms%fx(i)
      fa_y = atomic%atoms%fy(i)
      fa_z = atomic%atoms%fz(i)
       ! Numeric x force
      atomic%atoms%x(i) = atomic%atoms%x(i) + dx
      call SinglePoint(io,gen,atomic,tb,sol)
      ep = sol%totalEnergy
      atomic%atoms%x(i) = atomic%atoms%x(i) - 2.0_k_pr*dx
      call SinglePoint(io,gen,atomic,tb,sol)
      em = sol%totalEnergy
      atomic%atoms%x(i) = atomic%atoms%x(i) + dx
      fn_x = -(ep - em) / (dx + dx)
      totnx=totnx+fn_x
       ! Numeric y force
      atomic%atoms%y(i) = atomic%atoms%y(i) + dx
      call SinglePoint(io,gen,atomic,tb,sol)
      ep = sol%totalEnergy
      atomic%atoms%y(i) = atomic%atoms%y(i) - 2.0_k_pr*dx
      call SinglePoint(io,gen,atomic,tb,sol)
      em = sol%totalEnergy
      atomic%atoms%y(i) = atomic%atoms%y(i) + dx
      fn_y = -(ep - em) / (dx + dx)
      totny=totny+fn_y
       ! Numeric z force
      atomic%atoms%z(i) = atomic%atoms%z(i) + dx
      call SinglePoint(io,gen,atomic,tb,sol)
      ep = sol%totalEnergy
      atomic%atoms%z(i) = atomic%atoms%z(i) - 2.0_k_pr*dx
      call SinglePoint(io,gen,atomic,tb,sol)
      em = sol%totalEnergy
      atomic%atoms%z(i) = atomic%atoms%z(i) + dx
      fn_z = -(ep - em) / (dx + dx)
      totnz=totnz+fn_z

      write(*,'(a,i0)') "Atom: ",i
      write (*,'(a,F16.8,2X,F16.8,2X,F16.8)') 'Analytic: ',fa_x, fa_y, fa_z
      write (*,'(a,F16.8,2X,F16.8,2X,F16.8)') ' Numeric: ',fn_x, fn_y, fn_z
      write (*,'(a,F16.8,2X,F16.8,2X,F16.8)') 'An.-Num.: ',fa_x-fn_x, fa_y-fn_y, fa_z-fn_z
      write (*,*)
    end do
    fa_x=sum(atomic%atoms%fx(:))
    fa_y=sum(atomic%atoms%fy(:))
    fa_z=sum(atomic%atoms%fz(:))
    write(*,*)"Total forces:           x            y             z"
    write(*,'(a,F16.8,2X,F16.8,2X,F16.8)')"Analytic: ",fa_x,fa_y,fa_z
    write(*,'(a,F16.8,2X,F16.8,2X,F16.8)')" Numeric: ",totnx,totny,totnz
    write(*,'(a,F16.8,2X,F16.8,2X,F16.8)')"An.-Num.: ",fa_x-totnx,fa_y-totny,fa_z-totnz
  end subroutine forceTest

end module m_Testing