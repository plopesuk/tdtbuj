!> \brief does the Pulay mixing in the scf method
!> \author Alin M Elena
!> \date 08/11/07, 23:55:07
!> \remarks
!> http://dx.doi.org/10.1002/jcc.540030413
module m_Mixing
  use m_Constants
  use m_Types
  use m_Useful
#ifdef MKL95
  use mkl95_LAPACK, only: gesv
  use mkl95_BLAS, only: dot
#endif
  implicit none
  private
!
  public :: InitMix
  public :: MixDensity
!
contains
!> \brief shifts the vectors needed in mixing
!> \details the oldest from the colection will be discarded
!> \author Alin M Elena
!> \date 08/11/07, 15:11:47
!> \param ci,co,d matrices each line keeps a previous vecotr used in mixing
!> (ci input vectors, co output vectors, d residue vectors)
!> \param din,dout arrays new vectors that will be added
!> \param p integer number of vectors in collection
  subroutine InitMix (ci, co, d, din, dout, p)
    character (len=*), parameter :: myname = "InitMix"
    real (k_pr), intent (inout) :: ci (:, :), co (:, :), d (:, :)
    real (k_pr), intent (in) :: din (:), dout (:)
    integer, intent (in) :: p
    integer :: i
!
! first time there is nothing to shift
    do i = p, 2, - 1
      ci (:, i) = ci (:, i-1)
      co (:, i) = co (:, i-1)
      d (:, i) = d (:, i-1)
    end do
    ci (:, 1) = din (:)
    co (:, 1) = dout (:)
    d (:, 1) = co (:, 1) - ci (:, 1)
  end subroutine InitMix
!
!> \brief Pulay mixing
!> \author Alin M Elena
!> \date 08/11/07, 15:20:09
!> \param ci, co, d matrices each column keeps a previous vector used in mixing
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
!> \param io type(ioType) contains all the info about I/O files
  subroutine MixDensity (ci, co, d, dnext, residue, dmax, alpha, n, p, nit, ierr,io)
    character (len=*), parameter :: myname = "MixDensity"
    real (k_pr), intent (inout) :: ci (:, :), co (:, :), d (:, :)
    real (k_pr), intent (out) :: dnext (:)
    integer, intent (in) :: n, p, nit
    integer, intent (out) :: ierr
    real (k_pr), intent (out) :: residue, dmax
    real (k_pr), intent (in) :: alpha
    type(ioType), intent(in) :: io
    real (k_pr), allocatable :: c (:, :), beta (:), di (:), dou (:)
    integer :: i, j, nmix, info
    integer, allocatable :: wk (:)
#ifndef MKL95
    real(k_pr),external :: ddot
#endif
!
    allocate (di(1:n))
    di (1:n) = 0.0_k_pr
    nmix = Min (p, nit)
    info = 0
!
    if (nmix == 1) then
      dnext = alpha * co (:, 1) + (1-alpha) * ci (:, 1)
    else
      allocate (c(1:nmix+1, 1:nmix+1))
      allocate (beta(1:nmix+1))
      allocate (wk(1:nmix+1))
      allocate (dou(1:n))

      dou (1:n) = 0.0_k_pr
      do i = 1, nmix - 1
        do j = i + 1, nmix
#ifdef MKL95
          c (i, j) = dot (d(:, nmix-i+1), d(:, nmix-j+1))
#else
          c (i, j) = ddot(n,d(:, nmix-i+1), k_ione, d(:, nmix-j+1), k_ione )
#endif
          c (j, i) = c (i, j)
        end do
#ifdef MKL95
        c (i, i) = dot (d(:, nmix-i+1), d(:, nmix-i+1))
#else
        c (i, i) = ddot (n,d(:, nmix-i+1),k_ione, d(:, nmix-i+1),k_ione)
#endif
      end do

#ifdef MKL95
      c (nmix, nmix) = dot (d(:, 1), d(:, 1))
#else
      c (nmix, nmix) = ddot (n,d(:, 1),k_ione, d(:, 1),k_ione)
#endif
      beta (1:nmix+1) = 0.0_k_pr
!
      c (1:nmix, nmix+1) = 1.0_k_pr
      c (nmix+1, 1:nmix) = 1.0_k_pr
      c (nmix+1, 1+nmix) = 0.0_k_pr
      beta (1:nmix) = 0.0_k_pr
      beta (nmix+1) = 1.0_k_pr
      wk (1:nmix+1) = 0
#ifdef MKL95
      call gesv (c, beta, wk, info)
#else
      call dgesv(nmix+1,k_ione,c,nmix+1,wk,beta,nmix+1,info)
#endif
!
      if (info == 0) then
        do i = 1, nmix
          do j = 1, n
            di (j) = di (j) + beta (i) * ci (j, nmix-i+1)
            dou (j) = dou (j) + beta (i) * co (j, nmix-i+1)
          end do
        end do
        dnext = alpha * dou + (1-alpha) * di
      end if

      deallocate (c)
      deallocate (beta)
      deallocate (wk)
      deallocate (dou)
    end if
!
    di = co (:, nmix) - ci (:, nmix)
    dmax = Abs (maxval(di))
    di = di * di
#ifdef MKL95
    residue = Sqrt (dot(di, di)) / real (n, k_pr)
#else
    residue = Sqrt (ddot(n,di,k_ione, di,k_ione)) / real (n, k_pr)
#endif
    ierr = info
    deallocate (di)
  end subroutine MixDensity
!
end module m_Mixing
