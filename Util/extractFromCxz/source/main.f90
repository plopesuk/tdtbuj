program ReadCxz
  use m_Useful
  use m_Read
  use m_Gutenberg
  implicit none
  character (len=*), parameter :: myname = "Readcxz"
  integer :: narguments, outputunit, errno, err, i, unitCurrent
  character (len=k_mw) :: inpFile, arg, fileCurrent
  type (cxzType) :: cxz
  integer :: at
!
!
!
  outputunit = 6
  narguments = iargc ()
!
  if (narguments == 3) then
    call getarg (1, arg)
    inpFile = arg
    call getarg (2, arg)
    read (arg,*) at
    cxz%atom1 = at
    call getarg (3, arg)
    read (arg,*) at
    cxz%atom2 = at
    if (cxz%atom1 == cxz%atom2) then
      call error ("the atoms are the same!!", myname, .true., outputunit)
    end if
  else
    call error ("wrong number of arguments!!--- usage extractCxz filename at1 at2", myname, .true., outputunit)
  end if
!
  write (fileCurrent, '(a,i0,a1,i0,a)') "currents", cxz%atom1, "_", cxz%atom2, ".out"
!
  call InitRead (inpFile, cxz, outputunit)
!! allocate and initialize it memory
  allocate (cxz%frames(1:cxz%nframes), stat=errno)
  err = Abs (errno)
  if (errno > 0) then
    call error ("Error allocating the memory ", myname, .true., outputunit)
  end if
  call ReadFrames (inpFile, cxz, outputunit)
  unitCurrent = GetUnit ()
  open (unit=unitCurrent, file=trim(fileCurrent), status="replace", action="write")
  call PrintCurrents (cxz, unitCurrent)
  close (unitCurrent)
!
  deallocate (cxz%frames)
!
end program ReadCxz
