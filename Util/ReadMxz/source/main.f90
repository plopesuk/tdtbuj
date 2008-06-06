program ReadMxz
use m_Useful
use m_Read
use m_Gutenberg
implicit none
character(len=*), parameter :: myname="ReadMxz"
integer :: narguments,outputunit,errno,err,i,unitcharge,unitdipole
character(len=k_mw) :: inpFile,arg,fileCharge,fileDipole
type(mxzType) :: mxz


  fileCharge="charges.out"
  fileDipole="dipoles.out"
  outputunit=6
  narguments=iargc()

  if (narguments==1) then
    call getarg(1,arg)
    inpFile=arg
  else
    call error("wrong number of arguments!!",myname,.true.,outputunit)
  endif

  call InitRead(inpfile,mxz,outputunit)
!! allocate and initialize it memory
  allocate(mxz%frames(1:mxz%nframes),stat=errno)
  err=abs(errno)
  do i=1, mxz%nframes
    allocate(mxz%frames(i)%x(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
    allocate(mxz%frames(i)%y(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
    allocate(mxz%frames(i)%z(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
    allocate(mxz%frames(i)%dx(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
    allocate(mxz%frames(i)%dy(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
    allocate(mxz%frames(i)%dz(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
    allocate(mxz%frames(i)%q(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
    allocate(mxz%frames(i)%element(1:mxz%natoms),stat=errno)
    err=err+abs(errno)
  enddo
  if (errno > 0 ) then
    call error("Error allocating the memory ",myname,.true.,outputunit)
  endif
  do i=1, mxz%nframes
    mxz%frames(i)%x=0.0_k_pr
    mxz%frames(i)%y=0.0_k_pr
    mxz%frames(i)%z=0.0_k_pr
    mxz%frames(i)%dx=0.0_k_pr
    mxz%frames(i)%dy=0.0_k_pr
    mxz%frames(i)%dz=0.0_k_pr
    mxz%frames(i)%q=0.0_k_pr
    mxz%frames(i)%element=""
  enddo



  call ReadFrames(inpfile,mxz,outputunit)
  unitcharge=GetUnit()
  open(unit=unitcharge,file=trim(fileCharge),status="replace",action="write")
  call PrintCharges(mxz,unitcharge)
  close(unitcharge)
  unitdipole=GetUnit()
  open(unit=unitdipole,file=trim(fileDipole),status="replace",action="write")
  call PrintDipoles(mxz,unitdipole)
  close(unitdipole)

  do i=1, mxz%nframes
    deallocate(mxz%frames(i)%x)
    deallocate(mxz%frames(i)%y)
    deallocate(mxz%frames(i)%z)
    deallocate(mxz%frames(i)%dx)
    deallocate(mxz%frames(i)%dy)
    deallocate(mxz%frames(i)%dz)
    deallocate(mxz%frames(i)%q)
    deallocate(mxz%frames(i)%element)
  enddo
  deallocate(mxz%frames)

end program ReadMxz