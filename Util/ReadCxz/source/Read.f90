module m_Read
use m_Useful
implicit none
private

public :: InitRead
public :: ReadFrames

contains

  subroutine InitRead(file,cxz,unitout)
  character(len=*), parameter :: myname="InitRead"
  character(len=k_mw), intent(inout) :: file
  type(cxzType),intent(inout) :: cxz
  integer, intent(inout) :: unitout

  integer :: errno,n,k,frames,z,i,j
  character(len=k_ml) :: line
  integer :: unitinp
  logical :: status

    unitinp=GetUnit()
    open(unit=unitinp,file=trim(file),status="old",action="read",iostat=errno)

    if (errno /= 0 ) then
      call error("Error reading file "//trim(file),myname,.true.,unitout)
    endif
    ! determine no of atoms and no of frames
    read(unitinp,*,iostat=errno) n
    if (errno /= 0 ) then
      call error("Error reading no of atoms from file "//trim(file),myname,.true.,unitout)
    endif
    cxz%natoms=n
    status=GetLine(unitinp,line)
    do i=1,n
      status=GetLine(unitinp,line)
      if (.not. status) then
        call error("Error reading atoms from file "//trim(file),myname,.true.,unitout)
      endif
    enddo
    k=n+2
    read(unitinp,*,iostat=errno) n
    if (errno /= 0 ) then
      call error("Error reading no of currents from file "//trim(file),myname,.true.,unitout)
    endif
    k=k+1
    do i=1,n
      k=k+1
      read(unitinp,*,iostat=errno) z
      if (errno /= 0 ) then
        call error("Error reading no of currents from file "//trim(file),myname,.true.,unitout)
      endif
      do j=1,z
        status=GetLine(unitinp,line)
        if (.not. status) then
          call error("Error reading atoms from file "//trim(file),myname,.true.,unitout)
        endif
        k=k+1
      enddo
    enddo
    z=k
    do while (GetLine(unitinp,line))
      k =k +1
    enddo
    if (mod(k,z) == 0) then
      cxz%nframes = k/z
    else
      call error("Non integer number of frames in file "//trim(file),myname,.true.,unitout)
    endif
    close(unitinp)
  end subroutine InitRead


  subroutine ReadFrames(file,cxz,unitout)
    character(len=*), parameter :: myname="ReadFrames"
    character(len=k_mw), intent(inout) :: file
    type(cxzType),intent(inout) :: cxz
    integer, intent(inout) :: unitout

    integer :: errno,err,i,n,j,k,zz
    character(len=k_ml) :: line
    character(len=k_mw) :: saux,a,b,c
    real(k_pr) :: x,y,z,curr1,curr2,q
    integer :: unitinp
    logical :: status

    unitinp=GetUnit()
    open(unit=unitinp,file=trim(file),status="old",action="read",iostat=errno)
    do i=1,cxz%nframes
      read(unitinp,*,iostat=errno) n
      write(saux,'(i0)')i
      if (errno /= 0 ) then
        call error("Error reading first line of frame "//trim(saux),myname,.true.,unitout)
      endif
      read(unitinp,fmt='(a)',iostat=errno) a
      if (errno /= 0 ) then
        call error("Error reading label of frame "//trim(saux),myname,.true.,unitout)
      endif
      cxz%frames(i)%frameLabel = trim(a)
      read(a,*,iostat=errno)b,c,q
      if (errno /= 0 ) then
        call error("could non determine the timestamp of frame "//trim(saux),myname,.false.,unitout)
      endif
      cxz%frames(i)%timestamp = q
      err=0
      do j=1,cxz%natoms
         read(unitinp,*,iostat=errno) cxz%frames(i)%element(j),cxz%frames(i)%x(j),cxz%frames(i)%y(j),cxz%frames(i)%z(j),&
            cxz%frames(i)%curr1(j),cxz%frames(i)%curr2(j)
         err=err+abs(errno)
      enddo
      if (err /= 0 ) then
        call error("error readind atoms data from frame "//trim(saux),myname,.true.,unitout)
      endif
      read(unitinp,*,iostat=errno) n
      if (errno /= 0 ) then
        call error("Error reading no of currents from file "//trim(file),myname,.true.,unitout)
      endif
      do k=1,n
        read(unitinp,*,iostat=errno) zz
        if (errno /= 0 ) then
          call error("Error reading no of currents from file "//trim(file),myname,.true.,unitout)
        endif
        do j=1,zz
          status=GetLine(unitinp,line)
          if (.not. status) then
            call error("Error reading atoms from file "//trim(file),myname,.true.,unitout)
          endif
        enddo
      enddo
    enddo
    close(unitinp)
  end subroutine ReadFrames

!> \brief logical function reads a line from a unit
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param uno integer, unit number from where to read the line
!> \param line character(len=*), the line that was read
!> \return .true. if successfull in reading the line, .false. otherwise

  logical function GetLine(uno,line)
    character(len=*), parameter :: myname = 'GetLine'
    character(len=k_ml), intent(out) :: line
    integer, intent(in) ::uno
    integer :: errno

    inquire(unit=uno,iostat=errno)
    GetLine=.false.
    if (errno/=0) then
      write(*,*)"Unexpected error opening the input file(s)"
      stop
    endif
    read(uno,fmt='(a)',iostat=errno)line
    if (errno==0) GetLine=.true.
  end function GetLine

end module m_Read