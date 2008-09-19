program ReadCxz
  use m_Useful
  use m_Read
  use m_Gutenberg
  implicit none
  character (len=*), parameter :: myname = "Readcxz"
  integer :: narguments, outputunit, errno, err, i, unitCurrent
  character (len=k_mw) :: inpFile, arg, fileCurrent
  type (cxzType) :: cxz
!
!
!
  fileCurrent = "currents.out"
  outputunit = 6
  narguments = iargc ()
!
  if (narguments == 1) then
    call getarg (1, arg)
    inpFile = arg
  else
    call error ("wrong number of arguments!!", myname, .true., outputunit)
  end if
!
  call InitRead (inpFile, cxz, outputunit)
!! allocate and initialize it memory
  allocate (cxz%frames(1:cxz%nframes), stat=errno)
  err = Abs (errno)
  do i = 1, cxz%nframes
    allocate (cxz%frames(i)%x(1:cxz%natoms), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%y(1:cxz%natoms), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%z(1:cxz%natoms), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%element(1:cxz%natoms), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%curr1(1:cxz%natoms), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%curr2(1:cxz%natoms), stat=errno)
    err = err + Abs (errno)
  end do
  if (errno > 0) then
    call error ("Error allocating the memory ", myname, .true., outputunit)
  end if
  do i = 1, cxz%nframes
    cxz%frames(i)%x = 0.0_k_pr
    cxz%frames(i)%y = 0.0_k_pr
    cxz%frames(i)%z = 0.0_k_pr
    cxz%frames(i)%curr1 = 0.0_k_pr
    cxz%frames(i)%curr2 = 0.0_k_pr
    cxz%frames(i)%element = ""
  end do
!
!
  call ReadFrames (inpFile, cxz, outputunit)
  unitCurrent = GetUnit ()
  open (unit=unitCurrent, file=trim(fileCurrent), status="replace", action="write")
  call PrintCurrents (cxz, unitCurrent)
  close (unitCurrent)
!
  do i = 1, cxz%nframes
    deallocate (cxz%frames(i)%x)
    deallocate (cxz%frames(i)%y)
    deallocate (cxz%frames(i)%z)
    deallocate (cxz%frames(i)%curr1)
    deallocate (cxz%frames(i)%curr2)
    deallocate (cxz%frames(i)%element)
  end do
  deallocate (cxz%frames)
!
end program ReadCxz
