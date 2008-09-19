program BfrCxz
  use m_Useful
  use m_Read
  use m_Gutenberg
  use m_MagneticField
  implicit none
  character (len=*), parameter :: myname = "Readcxz"
  integer :: narguments, outputunit, errno, err, i, unitCurrent
  character (len=k_mw) :: inpFile, arg, fileCurrent, fileBfield
  type (cxzType) :: cxz
  integer :: at
!
  fileCurrent = "currents.out"
  fileBfield = "bfield.out"
  outputunit = 6
  narguments = iargc ()
!
  if (narguments >= 2) then
    call getarg (1, arg)
    inpFile = arg
    call getarg (2, arg)
    read (arg,*) at
    cxz%nring = at
  else
    call error ("to few arguments!! BFromCxz filename <no of atoms in ring> <the list of atoms>", myname, .true., outputunit)
  end if
  if (cxz%nring <= 0) then
    call error ("your ring should have a positive number of atoms!!", myname, .true., outputunit)
  else
! it is a ring so the last element would be the same as the first one    
    allocate (cxz%ring(1:cxz%nring+1))
  end if
!
!
  if (narguments >= 2+cxz%nring) then
    do i = 1, cxz%nring
      call getarg (2+i, arg)
      read (arg,*) cxz%ring(i)
    end do
    cxz%ring (cxz%nring+1) = cxz%ring(1)
  else
    call error ("to few arguments!! BFromCxz filename <no of atoms in ring> <the list of atoms>", myname, .true., outputunit)
  end if
!
  call InitRead (inpFile, cxz, outputunit)
!! allocate and initialize it memory
  allocate (cxz%frames(1:cxz%nframes), stat=errno)
  err = Abs (errno)
  do i = 1, cxz%nframes
    allocate (cxz%frames(i)%x(1:cxz%nring+1), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%y(1:cxz%nring+1), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%z(1:cxz%nring+1), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%element(1:cxz%nring+1), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%curr1(1:cxz%nring+1), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%curr2(1:cxz%nring+1), stat=errno)
    err = err + Abs (errno)
    allocate (cxz%frames(i)%bondCurrent(1:cxz%nring+1), stat=errno)
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
    cxz%frames(i)%bondCurrent = 0.0_k_pr
    cxz%frames(i)%element = ""
  end do
!
!
  call ReadFrames (inpFile, cxz, outputunit)
  unitCurrent = GetUnit ()
  open (unit=unitCurrent, file=trim(fileCurrent), status="replace", action="write")
  call PrintCurrents (cxz, unitCurrent)
  close (unitCurrent)
  call computeB (cxz)
!
  unitCurrent = GetUnit ()
  open (unit=unitCurrent, file=trim(fileBfield), status="replace", action="write")
  call PrintMagneticField (cxz, unitCurrent)
  close (unitCurrent)
  do i = 1, cxz%nframes
    deallocate (cxz%frames(i)%x)
    deallocate (cxz%frames(i)%y)
    deallocate (cxz%frames(i)%z)
    deallocate (cxz%frames(i)%curr1)
    deallocate (cxz%frames(i)%curr2)
    deallocate (cxz%frames(i)%bondCurrent)
    deallocate (cxz%frames(i)%element)
  end do
  deallocate (cxz%ring)
  deallocate (cxz%frames)
!
end program BfrCxz
