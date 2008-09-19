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
  subroutine InitRead (file, mxz, unitout)
    character (len=*), parameter :: myname = "InitRead"
    character (len=k_mw), intent (inout) :: file
    type (mxzType), intent (inout) :: mxz
    integer, intent (inout) :: unitout
!
    integer :: errno, n, k, frames
    character (len=k_ml) :: line
    integer :: unitinp
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
    mxz%natoms = n
    k = 1
    do while (GetLine(unitinp, line))
      k = k + 1
    end do
    if (Mod(k, n+2) == 0) then
      mxz%nframes = k / (n+2)
    else
      call error ("Non integer number of frames in file "//trim(file), myname, .true., unitout)
    end if
    close (unitinp)
  end subroutine InitRead
!
!
  subroutine ReadFrames (file, mxz, unitout)
    character (len=*), parameter :: myname = "ReadFrames"
    character (len=k_mw), intent (inout) :: file
    type (mxzType), intent (inout) :: mxz
    integer, intent (inout) :: unitout
!
    integer :: errno, err, i, n, j
    character (len=k_ml) :: line
    character (len=k_mw) :: saux, a, b, c
    real (k_pr) :: x, y, z, dx, dy, dz, q, tdx, tdy, tdz, charge
    integer :: unitinp
!
    unitinp = GetUnit ()
    open (unit=unitinp, file=trim(file), status="old", action="read", iostat=errno)
    do i = 1, mxz%nframes
      read (unitinp,*, iostat=errno) n
      write (saux, '(i0)') i
      if (errno /= 0) then
        call error ("Error reading first line of frame "//trim(saux), myname, .true., unitout)
      end if
      read (unitinp, fmt='(a)', iostat=errno) a
      if (errno /= 0) then
        call error ("Error reading label of frame "//trim(saux), myname, .true., unitout)
      end if
      mxz%frames(i)%frameLabel = trim (a)
      read (a,*, iostat=errno) b, c, q
      if (errno /= 0) then
        call error ("could non determine the timestamp of frame "//trim(saux), myname, .false., unitout)
      end if
      mxz%frames(i)%timestamp = q
      err = 0
      do j = 1, mxz%natoms
        read (unitinp,*, iostat=errno) mxz%frames(i)%element(j), mxz%frames(i)%x(j), mxz%frames(i)%y(j), mxz%frames(i)%z(j), &
       & mxz%frames(i)%q(j), mxz%frames(i)%dx(j), mxz%frames(i)%dy(j), mxz%frames(i)%dz(j)
        err = err + Abs (errno)
      end do
      if (err /= 0) then
        call error ("error readind atoms data from frame "//trim(saux), myname, .true., unitout)
      end if
!! now compute the dipole it may seem like a waste of loops to do it here and not in previous loop but thing what mwy happen it     
!! reading fails    
      tdx = 0.0_k_pr
      tdy = 0.0_k_pr
      tdz = 0.0_k_pr
      charge = 0.0_k_pr
!$omp  parallel  private(j)      
!$omp do schedule(static) reduction(+ : tdx, tdy,tdz,charge)      
      do j = 1, mxz%natoms
        tdx = mxz%frames(i)%dx(j) + tdx
        tdy = mxz%frames(i)%dy(j) + tdy
        tdz = mxz%frames(i)%dz(j) + tdz
        charge = charge + mxz%frames(i)%q(j)
      end do
!$omp end do      
!$omp end parallel      
      mxz%frames(i)%tdx = tdx
      mxz%frames(i)%tdy = tdy
      mxz%frames(i)%tdz = tdz
      mxz%frames(i)%charge = charge
      mxz%frames(i)%dipole = Sqrt (mxz%frames(i)%tdx**2+mxz%frames(i)%tdy**2+mxz%frames(i)%tdz**2)
      if (Abs(mxz%frames(i)%dipole) > tiny(1.0_k_pr)) then
        mxz%frames(i)%tdx = mxz%frames(i)%tdx / mxz%frames(i)%dipole
        mxz%frames(i)%tdy = mxz%frames(i)%tdy / mxz%frames(i)%dipole
        mxz%frames(i)%tdz = mxz%frames(i)%tdz / mxz%frames(i)%dipole
      end if
      mxz%frames(i)%projDipole = mxz%frames(i)%dipole * &
     & (mxz%frames(1)%tdx*mxz%frames(i)%tdx+mxz%frames(1)%tdy*mxz%frames(i)%tdy+mxz%frames(1)%tdz*mxz%frames(i)%tdz)
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
