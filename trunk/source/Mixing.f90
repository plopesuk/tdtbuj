!> \brief does the Pulay mixing in the scf method
!> \author Alin M Elena
!> \date 08/11/07, 23:55:07
!> \remarks
!> http://dx.doi.org/10.1002/jcc.540030413
module m_Mixing
  use m_Constants
  use m_Types
  implicit none
  private

  public :: InitMix
  public :: MixDensity

contains
!> \brief shifts the vectors needed in mixing
!> \details the oldest from the colection will be discarded
!> \author Alin M Elena
!> \date 08/11/07, 15:11:47
!> \param ci,co,d matrices each line keeps a previous vecotr used in mixing
!> (ci input vectors, co output vectors, d residue vectors)
!> \param din,dout arrays new vectors that will be added
!> \param p integer number of vectors in collection
  subroutine InitMix(ci,co,d,din,dout,p)
    character(len=*), parameter :: myname="InitMix"
    real(k_pr), intent(inout) :: ci(:,:),co(:,:),d(:,:)
    real(k_pr),intent(in) :: din(:),dout(:)
    integer, intent(in) :: p
    integer :: i

! first time there is nothing to shift    
    do i=p,2,-1
      ci(:,i)=ci(:,i-1)
      co(:,i)=co(:,i-1)
      d(:,i)=d(:,i-1)
    enddo
    ci(:,1)=din(:)
    co(:,1)=dout(:)
    d(:,1)=co(:,1)-ci(:,1)
  end subroutine InitMix

!> \brief Pulay mixing
!> \author Alin M Elena
!> \date 08/11/07, 15:20:09
!> \param ci, co, d matrices each line keeps a previous vecotr used in mixing
!> (ci input vectors, co output vectors, d residue vectors)
!> \param dnext array the new vector generated
!> \param residue the residue vector
!> \param dmax real the maximum difference between the new vector and the old one.
!> smaller it is closer to the scf solution one is
!> \param alpha real the mixing factor
!> \param n integer the dimension of the vector to be mixed
!> \param p integer how deep to go back in "history"
!> \param nit integer current iteration
!> \param ierr integer error indicator 0 means success other numbers indicate failure
  subroutine  MixDensity(ci,co,d,dnext,residue,dmax,alpha,n,p,nit,ierr)
    character(len=*), parameter :: myname="MixDensity"
    real(k_pr), intent(inout) :: ci(:,:),co(:,:),d(:,:)
    real(k_pr),intent(out) :: dnext(:)
    integer, intent(in) :: n,p,nit
    integer, intent(out) :: ierr
    real(k_pr), intent(out) :: residue,dmax
    real(k_pr), intent(in)::alpha
    real(k_pr), allocatable ::c(:,:), beta(:),di(:),dou(:)
    integer :: i,j,nmix,info
    integer, allocatable :: wk(:)

    allocate(di(1:n))
    nmix=min(p,nit)
    info=0

    if (nmix==1) then
      dnext=alpha*co(:,1) + (1-alpha)*ci(:,1)
    else
      allocate(c(nmix+1,nmix+1))
      allocate(beta(nmix+1))
      allocate(wk(nmix+1))
      allocate(dou(1:n))
      do i=1,nmix
        do j=i+1,nmix
          c(i,j) = sum( d(:,nmix-i+1)*d(:,nmix-j+1) )
          c(j,i) = c(i,j)
        enddo
        c(i,i) = sum( d(:,nmix-i+1)*d(:,nmix-i+1) )
      enddo
      beta=0.0_k_pr
      c(1:nmix,nmix+1)=1.0_k_pr
      c(nmix+1,1:nmix)=1.0_k_pr
      c(nmix+1,1+nmix)=0.0_k_pr
      beta(1:nmix)=0.0_k_pr
      beta(nmix+1)=1.0_k_pr
      wk=0
      call dgesv(nmix+1,1,c,nmix+1,wk,beta,nmix+1,info)
      if (info==0) then
        di=0.0_k_pr
        dou=0.0_k_pr
        do i=1,nmix
          di=di+beta(i)*ci(:,nmix-i+1)
          dou=dou+beta(i)*co(:,nmix-i+1)
        enddo
        dnext=alpha*dou + (1-alpha)*di
      endif
      deallocate(c)
      deallocate(beta)
      deallocate(wk)
      deallocate(dou)
    endif
      di = co(:,nmix) - ci(:,nmix)
      dmax = abs(maxval(di))
      di = di*di
      residue = sqrt(sum(di*di))/real(n,k_pr)
      ierr=info
      deallocate(di)
   end subroutine MixDensity

end module m_Mixing