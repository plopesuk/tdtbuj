module m_Gutenberg
use m_Useful
implicit none
private

public :: PrintCharges
public :: PrintDipoles

contains

  subroutine PrintCharges(mxz,unitout)
    character(len=*), parameter :: myname="PrintCharges"
    type(mxzType),intent(inout) :: mxz
    integer, intent(inout) :: unitout

    character(len=k_mw) :: fmtContainer,saux
    integer :: i

    write(unitout,'(a1,a16,a16)',advance="no")"#","Time ","Total Charge "
    do i=1,mxz%natoms
     write(unitout,'(a4,i8,a1,a2,a1)',advance="no")" Col",i+2,"-",trim(mxz%frames(1)%element(i)),";"
    enddo
    write(unitout,*)
    write(unitout,*)"#to get the atom number substract 2 from col number"
    write(fmtContainer,'("(f16.8,1x,f16.8,",i0 ,"f16.8)")') mxz%natoms
    do i=1,mxz%nframes
      write(unitout,trim(fmtContainer))mxz%frames(i)%timestamp,mxz%frames(i)%charge,mxz%frames(i)%q
    enddo
  end subroutine PrintCharges

  subroutine PrintDipoles(mxz,unitout)
    character(len=*), parameter :: myname="PrintDipoles"
    type(mxzType),intent(inout) :: mxz
    integer, intent(inout) :: unitout

    character(len=k_ml) :: fmtContainer,saux
    integer :: i

    write(unitout,'(a1,5a16)',advance="no")"#","Time ","Total Dipole ","DX ", "DY ", "DZ "
    do i=1,mxz%natoms
      write(unitout,'(a4,i8,a1,a2,a1)',advance="no")"XCol",i+5,"-",trim(mxz%frames(1)%element(i)),";"
    enddo
    do i=1,mxz%natoms
      write(unitout,'(a4,i8,a1,a2,a1)',advance="no")"YCol",i+5+mxz%natoms,"-",trim(mxz%frames(1)%element(i)),";"
    enddo
    do i=1,mxz%natoms
      write(unitout,'(a4,i8,a1,a2,a1)',advance="no")"ZCol",i+5+2*mxz%natoms,"-",trim(mxz%frames(1)%element(i)),";"
    enddo
    write(unitout,*)
    write(unitout,*)"#to get the atom number for X component substract 5 from last number after XCol"
    write(unitout,'(a,i0,a)')"#to get the atom number for Y component substract ",5+mxz%natoms," from last number after YCol"
    write(unitout,'(a,i0,a)')"#to get the atom number for Z component substract ",5+2*mxz%natoms," from last number after ZCol"
    write(fmtContainer,'("(f16.8,1x,4f16.8,",i0,"f16.8,",i0,"f16.8,",i0,"f16.8)")') mxz%natoms,mxz%natoms,mxz%natoms
    do i=1,mxz%nframes
      write(unitout,trim(fmtContainer))mxz%frames(i)%timestamp,mxz%frames(i)%dipole,mxz%frames(i)%tdx,mxz%frames(i)%tdy,mxz%frames(i)%tdz,&
        mxz%frames(i)%dx,mxz%frames(i)%dy,mxz%frames(i)%dz
    enddo
  end subroutine PrintDipoles


end module m_Gutenberg