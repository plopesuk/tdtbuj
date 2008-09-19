module m_Read
  use m_Useful
  implicit none
  private
!
  public :: InitRead
  public :: ReadFrames
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
    character (len=k_ml) :: line, saux
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
    do i = 1, cxz%nring
      if ((cxz%ring(i) > n) .or. (cxz%ring(i) < 1)) then
        write (saux, '(a,i0,a,a)') "Atom ", cxz%ring(i), " is not contained by the file ", trim (file)
        call error (trim(saux), myname, .true., unitout)
      end if
    end do
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
  subroutine ReadFrames (file, cxz, unitout)
    character (len=*), parameter :: myname = "ReadFrames"
    character (len=k_mw), intent (inout) :: file
    type (cxzType), intent (inout) :: cxz
    integer, intent (inout) :: unitout
!
    integer :: errno, err, i, n, j, k, zz
    character (len=k_ml) :: line
    character (len=k_mw) :: saux, a, b, c, element
    real (k_pr) :: x, y, z, curr1, curr2, q
    integer :: unitinp, pj, pj2, at1, at2
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
      read (unitinp, fmt='(a)', iostat=errno) a
      if (errno /= 0) then
        call error ("Error reading label of frame "//trim(saux), myname, .true., unitout)
      end if
      cxz%frames(i)%frameLabel = trim (a)
      read (a,*, iostat=errno) b, c, q
      if (errno /= 0) then
        call error ("could non determine the timestamp of frame "//trim(saux), myname, .false., unitout)
      end if
      cxz%frames(i)%timestamp = q
      err = 0
      do j = 1, cxz%natoms
        read (unitinp,*, iostat=errno) element, x, y, z, curr1, curr2
        err = err + Abs (errno)
        if (err /= 0) then
          call error ("error reading atoms data from frame "//trim(saux), myname, .true., unitout)
        end if
        if (isInList(j, cxz%ring(1:cxz%nring), pj)) then
          cxz%frames(i)%element(pj) = trim (element)
          cxz%frames(i)%x(pj) = x
          cxz%frames(i)%y(pj) = y
          cxz%frames(i)%z(pj) = z
          cxz%frames(i)%curr1(pj) = curr1
          cxz%frames(i)%curr2(pj) = curr2
        end if
      end do
      read (unitinp,*, iostat=errno) n
      if (errno /= 0) then
        call error ("Error reading no of currents from file "//trim(file), myname, .true., unitout)
      end if
      do k = 1, n
        read (unitinp,*, iostat=errno) zz
        if (errno /= 0) then
          call error ("Error reading no of currents from file "//trim(file), myname, .true., unitout)
        end if
        do j = 1, zz
          read (unitinp,*, iostat=errno) at1, at2, curr1
          do pj = 1, cxz%nring
            if ((at1 == cxz%ring(pj)) .and. (at2 == cxz%ring(pj+1))) then
              cxz%frames(i)%bondCurrent(pj) = curr1
            end if
          end do
          if (errno /= 0) then
            call error ("Error reading atoms from file "//trim(file), myname, .true., unitout)
          end if
        end do
      end do
      cxz%frames(i)%element(cxz%nring+1) = cxz%frames(i)%element(1)
      cxz%frames(i)%x(cxz%nring+1) = cxz%frames(i)%x(1)
      cxz%frames(i)%y(cxz%nring+1) = cxz%frames(i)%y(1)
      cxz%frames(i)%z(cxz%nring+1) = cxz%frames(i)%z(1)
      cxz%frames(i)%curr1(cxz%nring+1) = cxz%frames(i)%curr1(1)
      cxz%frames(i)%curr2(cxz%nring+1) = cxz%frames(i)%curr2(1)
    end do
    close (unitinp)
  end subroutine ReadFrames
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
