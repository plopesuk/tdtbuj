module m_MagneticField
use m_Useful
implicit none
private

public :: computeB

contains

subroutine computeB(cxz)
  character(len=k_mw), parameter :: myname="computeB"
  type(cxzType), intent(inout) :: cxz

  real(k_pr) :: rc(1:3),r,sina,cosa,rot(1:3,1:3),rm(1:3),rcp(1:3),invrot(1:3,1:3),a(1:3),b(1:3),mag(1:3),xa,gB(1:3),length
  integer :: i,j,k,m

  do j=1,cxz%nframes
  ! compute the centre of the ring
    rc=0.0_k_pr
    do i=1,cxz%nring
      rc(1)=rc(1)+cxz%frames(j)%x(i)
      rc(2)=rc(2)+cxz%frames(j)%y(i)
      rc(3)=rc(3)+cxz%frames(j)%z(i)
    enddo
    rc=rc/real(cxz%nring,k_pr)
    gb=0.0_k_pr
    do i=1,cxz%nring
     ! determine the rotation angle for each bond
     ! we start by computing the middle of the bond
     a(1)=cxz%frames(j)%x(i)
     a(2)=cxz%frames(j)%y(i)
     a(3)=cxz%frames(j)%z(i)
     b(1)=cxz%frames(j)%x(i+1)
     b(2)=cxz%frames(j)%y(i+1)
     b(3)=cxz%frames(j)%z(i+1)
     length=norm(b-a)
     sina=(a(1)-b(1))/length
     cosa=(b(2)-a(2))/length
     call Buildrotz(rot,cosa,-sina)
     ! matrix rot S'=rotS
     call Buildrotz(invrot,cosa,sina)
     ! matrix invrot S=rotS'
     rcp=matmul(rot,rc)
     a=matmul(rot,a)
     b=matmul(rot,b)
     xa=a(2)
     r=a(1)
     call localB(mag,rcp,xa,length,r,cxz%frames(j)%bondCurrent(i)*k_mju0/(4.0_k_pr*k_pi))
     gB=gB+matmul(invrot,mag)
    enddo
    cxz%frames(j)%b=gb*k_B2tesla

  enddo
  end subroutine computeB

!> \brief builds the rotation matrix around z
!> \author Alin M Elena
!> \date 18/08/2008, 10:41
!> \param rot the rotation matrix
!> \param ca real cos(alpha)
!> \param sa real sin(alpha)
!> \remarks \f[
!> R_z(\alpha)=\left(
!> \begin{array}{ccc}
!> \cos \alpha & -\sin \alpha & 0 \\\\
!> \sin \alpha & \cos \alpha & 0 \\\\
!> 0 & 0 & 1 \\\\
!> \end{array}
!> \right)
!> \f]
  subroutine Buildrotz(rot,ca,sa)
    character(len=*), parameter :: myname="Buildrotz"
    real(k_pr), intent(inout) :: rot(:,:)
    real(k_pr), intent(in) :: ca,sa

    rot=0.0_k_pr
    rot(1,1)=ca
    rot(1,2)=-sa
    rot(1,3)=0.0_k_pr
    rot(2,1)=sa
    rot(2,2)=ca
    rot(2,3)=0.0_k_pr
    rot(3,1)=0.0_k_pr
    rot(3,2)=0.0_k_pr
    rot(3,3)=1.0_k_pr

  end subroutine Buildrotz

  subroutine localB(b,r,a,length,xr,pre)
  character(len=*), parameter :: myname="localB"
  real(k_pr), intent(inout) :: b(:)
  real(k_pr), intent(in) :: r(:),a,xr,pre,length
  real(k_pr):: d,h

    d=sqrt((r(1)-xr)*(r(1)-xr)+r(3)*r(3))
    h=((a+length-r(2))/sqrt(d*d+(a+length-r(2))*(a+length-r(2)))-(a-r(2))/sqrt(d*d+(a-r(2))*(a-r(2))))/(d*d)
    b(2) = 0.0_k_pr
    b(1) = pre*h*r(3)
    b(3) =-pre*h*(r(1)-xr)

  end subroutine localB




end module m_MagneticField