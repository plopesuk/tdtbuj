program changeUnitsCxz
  use m_Useful
  use m_Read
  implicit none
  character (len=*), parameter :: myname = "changeUnitsCxz"
  integer :: narguments, outputunit, errno, err, i, unitCurrent
  character (len=k_mw) :: inpFile, arg, fileCurrent
  type (cxzType) :: cxz
!
!
!
  narguments = iargc ()
!
  if (narguments == 1) then
    call getarg (1, arg)
    inpFile = arg
  else
    call error ("wrong number of arguments!!", myname, .true., outputunit)
  end if
  fileCurrent = "microAmps"//trim(inpfile)
  outputunit = 6
  open (unit=outputunit, file=trim(fileCurrent), status="unknown", action="write", iostat=errno)
!
  call InitRead (inpFile, cxz, outputunit)
!
!
  call ReadAndConvertFrames (inpFile, cxz, outputunit)
!
  close (unit=outputunit)
end program changeUnitsCxz
