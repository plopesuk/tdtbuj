module m_Read
  use m_Useful
  implicit none
  private
!
  public :: InitRead
  public :: ReadAndConvertFrames
!
contains
!
  subroutine InitRead (file, cxz, unitout)
    character (len=*), parameter :: myname = "InitRead"
    character (len=k_mw), intent (inout) :: file
    type (cxzType), intent (inout) :: cxz
    integer, intent (inout) :: unitout
!
    integer :: errno, n, k, frames, z, i, j
    character (len=k_ml) :: line
    integer :: unitinp
    logical :: status
!
    unitinp = GetUnit ()
    open (unit=unitinp, file=trim(file), status="old", action="read", iostat=errno)
!
    if (errno /= 0) then
      call error ("Error reading file "//trim(file), myname, .true., unitout)
    end if
! determine no of atoms and no of frames
    read (unitinp,*, iostat=errno) n
    if (errno /= 0) then
      call error ("Error reading no of atoms from file "//trim(file), myname, .true., unitout)
    end if
    cxz%natoms = n
    status = GetLine (unitinp, line)
    do i = 1, n
      status = GetLine (unitinp, line)
      if ( .not. status) then
        call error ("Error reading atoms from file "//trim(file), myname, .true., unitout)
      end if
    end do
    k = n + 2
    read (unitinp,*, iostat=errno) n
    if (errno /= 0) then
      call error ("Error reading no of currents from file "//trim(file), myname, .true., unitout)
    end if
    k = k + 1
    do i = 1, n
      k = k + 1
      read (unitinp,*, iostat=errno) z
      if (errno /= 0) then
        call error ("Error reading no of currents from file "//trim(file), myname, .true., unitout)
      end if
      do j = 1, z
        status = GetLine (unitinp, line)
        if ( .not. status) then
          call error ("Error reading atoms from file "//trim(file), myname, .true., unitout)
        end if
        k = k + 1
      end do
    end do
    z = k
    do while (GetLine(unitinp, line))
      k = k + 1
    end do
    if (Mod(k, z) == 0) then
      cxz%nframes = k / z
    else
      call error ("Non integer number of frames in file "//trim(file), myname, .true., unitout)
    end if
    close (unitinp)
  end subroutine InitRead
!
!
  subroutine ReadAndConvertFrames (file, cxz, unitout)
    character (len=*), parameter :: myname = "ReadAndConvertFrames"
    character (len=k_mw), intent (inout) :: file
    type (cxzType), intent (inout) :: cxz
    integer, intent (inout) :: unitout
!
    integer :: errno, err, i, n, j, k, zz,at1,at2
    character (len=k_ml) :: line
    character (len=k_mw) :: saux, a, b, c,element
    real (k_pr) :: x, y, z, curr1, q,x1min,x1max,x2min,x2max
    integer :: unitinp
    logical :: status
!
    unitinp = GetUnit ()
    open (unit=unitinp, file=trim(file), status="old", action="read", iostat=errno)
    do i = 1, cxz%nframes
      read (unitinp,*, iostat=errno) n
      write (saux, '(i0)') i
      if (errno /= 0) then
        call error ("Error reading first line of frame "//trim(saux), myname, .true., unitout)
      end if
      write(unitout,*) n
      read (unitinp, fmt='(a)', iostat=errno) a
      if (errno /= 0) then
        call error ("Error reading label of frame "//trim(saux), myname, .true., unitout)
      end if
      write(unitout,*) trim(a)
      if (errno /= 0) then
        call error ("could non determine the timestamp of frame "//trim(saux), myname, .false., unitout)
      end if
      err = 0
      do j = 1, cxz%natoms
        read (unitinp,*, iostat=errno) element, x, y, z,curr1
        if (errno /= 0) then
          call error ("error readind atoms data from frame "//trim(saux), myname, .true., unitout)
        end if
        if ((i==1).and.(j==1)) then
          x1min=k_conv*curr1
          x1max=k_conv*curr1
        else
          if (x1min>k_conv*curr1) then
            x1min=k_conv*curr1
          endif
          if (x1max<k_conv*curr1) then
            x1max=k_conv*curr1
          endif
        endif
        write(unitout,'(a,1x,4g)')trim(element),x,y,z,k_conv*curr1
      end do

      read (unitinp,*, iostat=errno) n
      if (errno /= 0) then
        call error ("Error reading no of currents from file "//trim(file), myname, .true., unitout)
      end if
      write(unitout,*)n
      do k = 1, n
        read (unitinp,*, iostat=errno) zz
        if (errno /= 0) then
          call error ("Error reading no of currents from file "//trim(file), myname, .true., unitout)
        end if
        write(unitout,*)zz
        do j = 1, zz
          read (unitinp,*, iostat=errno)at1,at2,curr1
          if (errno /=0) then
            call error ("Error reading atoms from file "//trim(file), myname, .true., unitout)
          end if
          write(unitout,'(i0,x,i0,g)')at1,at2,k_conv*curr1
          if ((i==1).and.(j==1)) then
            x2min=k_conv*curr1
            x2max=k_conv*curr1
          else
            if (x2min>k_conv*curr1) then
              x2min=k_conv*curr1
            endif
            if (x2max<k_conv*curr1) then
              x2max=k_conv*curr1
            endif
          endif
        end do
      end do
    end do
    write(*,'(a,2g)')"Min/Max for current on atoms: ",x1min,x1max
    write(*,'(a,2g)')"Min/Max for bond currents: ",x2min,x2max
    close (unitinp)
  end subroutine ReadAndConvertFrames
!
!> \brief logical function reads a line from a unit
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param uno integer, unit number from where to read the line
!> \param line character(len=*), the line that was read
!> \return .true. if successfull in reading the line, .false. otherwise
!
  logical function GetLine (uno, line)
    character (len=*), parameter :: myname = 'GetLine'
    character (len=k_ml), intent (out) :: line
    integer, intent (in) :: uno
    integer :: errno
!
    inquire (unit=uno, iostat=errno)
    GetLine = .false.
    if (errno /= 0) then
      write (*,*) "Unexpected error opening the input file(s)"
      stop
    end if
    read (uno, fmt='(a)', iostat=errno) line
    if (errno == 0) GetLine = .true.
  end function GetLine
!
end module m_Read
