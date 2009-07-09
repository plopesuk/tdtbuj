!>\brief This module defines the matrix type and
!> implements all linear algebra subroutines
!> \details by apropriate blas/lapack calls or special sparse
!> capable subroutines.
!> \author Alin M Elena
!> based on a  version by Cristian G. Sanchez (Belfast)
!> \remarks
!> changed by Alin M Elena: the routines do not depend on global variables
!> changed to comply with the coding style document
!
module m_LinearAlgebra
  use m_Constants
  use m_Useful
  use m_Types
#ifdef MKL95
  use mkl95_lapack, only: heevr
  use mkl95_blas, only: gemm, herk
#endif
  implicit none
  private
!
  public :: CreateMatrix
  public :: CreateSparseMatrix
  public :: ResetSparseMatrix
  public :: SpmPut
  public :: DestroyMatrix
  public :: ZeroMatrix
  public :: ZeroDiagonalMatrix
  public :: CopyMatrix
!
  public :: MatrixCeaApbB
  public :: MatrixTrace
  public :: ProductTrace
  public :: ScalarTMatrix
  public :: Commutator
  public :: DiagonalizeMatrix
  public :: aastar
!
contains
!-------------------------------------------------!
! Allocates memory and posibly Initializes the !
! matrix a !
!-------------------------------------------------!
  subroutine CreateMatrix (a, m, zeroout)
    character (len=*), parameter :: myname = 'CreateMatrix'
    type (matrixType), intent (inout) :: a
    logical, intent (in) :: zeroout
    integer, intent (in) :: m
    integer :: i, j
!
    a%isSparse = .false.
    a%dim = m
    a%nonZero = 1
    allocate (a%indx(1:1))
    allocate (a%jndx(1:1))
    allocate (a%a(m, m))
    a%indx = 0
    a%jndx = 0
    if (zeroout) then
      a%a (1:a%dim, 1:a%dim) = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
    end if
    a%created = .true.
  end subroutine CreateMatrix
!-------------------------------------------------!
! Allocates memory and and Initializes the sparse !
! matrix a !
!-------------------------------------------------!
!
  subroutine CreateSparseMatrix (a, m, zeroout)
    character (len=*), parameter :: myname = 'CreateSparseMatrix'
    type (matrixType), intent (inout) :: a
    logical, intent (in) :: zeroout
    integer, intent (in) :: m
    integer :: i, j
!
    a%isSparse = .true.
    a%dim = m
    a%nonZero = 0
    allocate (a%a(1:m, 1:m))
    allocate (a%indx(1:m*m))
    allocate (a%jndx(1:m*m))
    if (zeroout) then
      a%a (1:a%dim, 1:a%dim) = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
      a%indx (1:a%dim*a%dim) = 0
      a%jndx (1:a%dim*a%dim) = 0
    end if
    a%created = .true.
  end subroutine CreateSparseMatrix
!-------------------------------------------------!
! Resets the index arrays in the sparse matrix a !
!-------------------------------------------------!
  subroutine ResetSparseMatrix (a)
    character (len=*), parameter :: myname = 'ResetSparseMatrix'
    type (matrixType), intent (inout) :: a
    integer :: i
    do i = 1, a%nonZero
      a%a (a%indx(i), a%jndx(i)) = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
    end do
    a%indx (1:a%nonZero) = 0
    a%jndx (1:a%nonZero) = 0
    a%nonZero = 0
  end subroutine ResetSparseMatrix
!
!-------------------------------------------------!
! Puts a new element into the sparse matrix a !
!-------------------------------------------------!
  subroutine SpmPut (a, i, j, aij)
    character (len=*), parameter :: myname = 'SpmPut'
    type (matrixType), intent (inout) :: a
    complex (k_pr), intent (in) :: aij
    integer, intent (in) :: i, j
!
!     if ( .not. IsIn(a, i, j)) then
!       a%nonZero = a%nonZero + 1
!       a%indx (a%nonZero) = i
!       a%jndx (a%nonZero) = j
!     end if
    a%a (i, j) = aij
  end subroutine SpmPut
!
!> \brief checks if a certain element is non zero in a sparse matrix
!> \author Alin M Elena
!> \date 09/11/07, 20:53:19
!> \param r,c integers row, column
!> \param a type(matrixType) the matrix
  logical function IsIn (a, r, c)
    integer, intent (in) :: r, c
    type (matrixType), intent (in) :: a
    integer :: t
    IsIn = .false.
    do t = 1, a%nonZero
      if ((a%indx(t) == r) .and. (a%jndx(t) == c)) then
        IsIn = .true.
        exit
      end if
    end do
  end function IsIn
!-------------------------------------------------!
! k_zeroes matrix a !
!-------------------------------------------------!
  subroutine ZeroMatrix (a, io)
    character (len=*), parameter :: myname = 'ZeroMatrix'
    type (matrixType), intent (inout) :: a
    type (ioType), intent (in) :: io
    integer :: i, j
!
    if ( .not. a%created) then
      call error ("Cannot fill an unexistant matrix", myname, .true., io)
    end if
    a%a = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
!
  end subroutine ZeroMatrix
!
!
  subroutine ZeroDiagonalMatrix (a, io)
    character (len=*), parameter :: myname = 'ZeroDiagonalMatrix'
    type (matrixType), intent (inout) :: a
    type (ioType), intent (in) :: io
    integer :: i
!
    if ( .not. a%created) then
      call error ("Cannot zero the diagonal of an unexistant matrix", myname, .true., io)
    end if
    do i = 1, a%dim
      a%a (i, i) = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
    end do

!
  end subroutine ZeroDiagonalMatrix
!
!-------------------------------------------------!
! Diagonalize a (hermitian) matrix !
!-------------------------------------------------!
  subroutine DiagonalizeMatrix (a, c, lambda, io)
    character (len=*), parameter :: myname = 'DiagonalizeMatrix'
    type (matrixType), intent (inout) :: a, c
    real (kind=k_pr), intent (inout) :: lambda (:)
    type (matrixType) :: h
    type (ioType), intent (in) :: io
    integer :: info
#ifndef MKL95
    real(k_pr) :: rtmp(1)
    real(k_pr), allocatable :: rwork(:)
    integer :: lwork,liwork,lrwork,ierror
    character(len=1) :: ta,rnge
    real(k_pr) :: vl,vu,abstol
    integer :: il,iu,itmp(1),m
    integer, allocatable :: iwork(:),issupz(:)
    complex(k_pr) :: wtmp(1)
    complex(k_pr), allocatable :: work(:)
    abstol=tiny(1.0_k_pr)*1000
#endif

    if ( .not. a%created) then
      call error ("Cannot diagonalize unexistant matrix", myname, .true., io)
    end if
#ifdef MKL95
    call heevr (a%a, lambda, uplo='U', z=c%a, info=info)
#else
    allocate(issupz(2*a%dim),stat=ierror)
    if (ierror /= 0) then
      call error ("Cannot allocate space for buffer issupz", myname, .true., io)
    endif
    lwork=-1
    liwork=-1
    lrwork=-1
    call zheevr("V","A","U", a%dim, a%a, a%dim,vl,vu,il,iu,abstol,m, lambda,&
       c%a,c%dim,issupz,wtmp, lwork,rtmp,lrwork,itmp,liwork, info)
    lwork=int(wtmp(1))
    liwork=itmp(1)
    lrwork=int(rtmp(1))
    allocate(work(lwork),iwork(liwork),rwork(lrwork),stat=ierror)
    if (ierror /= 0) then
      call error ("Cannot allocate space for buffers *work", myname, .true., io)
    endif
    call zheevr("V","A","U", a%dim, a%a, a%dim,vl,vu,il,iu,abstol,m, lambda,&
       c%a,c%dim,issupz,work, lwork,rwork,lrwork,iwork,liwork, info)
     deallocate(work,issupz,iwork,rwork,stat=ierror)
     if (ierror /= 0) then
       call error ("Cannot deallocate space for buffers", myname, .true., io)
     endif
#endif
    if (info /= 0) then
      call error ("Diagonalization failed", myname, .true., io)
    end if
  end subroutine DiagonalizeMatrix
!-------------------------------------------------!
! Deallocates memory in the matrix a !
!-------------------------------------------------!
  subroutine DestroyMatrix (a, io)
    character (len=*), parameter :: myname = 'DestroyMatrix'
    type (matrixType), intent (inout) :: a
    type (ioType), intent (in) :: io
!
    if (a%created) then
      deallocate (a%a)
      a%dim = 0
      a%created = .false.
      deallocate (a%indx)
      deallocate (a%jndx)
      if (a%isSparse) then
        a%isSparse = .false.
      end if
    else
      call error ("Cannot destroy unexistant matrix", myname, .true., io)
    end if
  end subroutine DestroyMatrix
!
!-------------------------------------------------!
! calculate the trace of a matrix !
!-------------------------------------------------!
  function MatrixTrace (a, io)
    character (len=*), parameter :: myname = 'MatrixTrace'
    complex (kind=k_pr) :: MatrixTrace
    type (matrixType), intent (in) :: a
    type (ioType), intent (in) :: io
    integer :: i
!
    if ( .not. a%created) then
      call error ("Cannot operate on unexistant matrix", myname, .true., io)
    end if
    MatrixTrace = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
    do i = 1, a%dim
      MatrixTrace = MatrixTrace + a%a(i, i)
    end do

!
  end function MatrixTrace
!-------------------------------------------------!
! Multiply scalar times matrix !
!-------------------------------------------------!
  subroutine ScalarTMatrix (s, a, io)
    character (len=*), parameter :: myname = 'ScalarTMatrix'
    type (matrixType), intent (inout) :: a
    complex (kind=k_pr), intent (in) :: s
    type (ioType), intent (in) :: io
    integer :: i, j
    if ( .not. a%created) then
      call error ("Cannot operate on an unexistant matrix", myname, .true., io)
    end if
    a%a (1:a%dim, 1:a%dim) = a%a(1:a%dim, 1:a%dim) * s

!
  end subroutine ScalarTMatrix
!
!-------------------------------------------------!
! Copy Matrix a into b !
!-------------------------------------------------!
  subroutine CopyMatrix (b, a, io)
    character (len=*), parameter :: myname = 'CopyMatrix'
    type (matrixType), intent (in) :: a
    type (matrixType), intent (inout) :: b
    type (ioType), intent (in) :: io
    integer :: i, j
    if (b%created .and. (b%dim /= a%dim)) then
      call error ("Dimensions are different", myname, .true., io)
    end if
    if (a%isSparse) then
      if ( .not. b%created) then
        call CreateSparseMatrix (b, a%dim, .true.)
      end if
      if (b%isSparse == .true.) then
        b%nonZero = a%nonZero
        b%indx (1:a%nonZero) = a%indx(1:a%nonZero)
        b%jndx (1:a%nonZero) = a%jndx(1:a%nonZero)
      end if
    else
      if ( .not. b%created) then
        call CreateMatrix (b, a%dim, .false.)
      end if
    end if
    b%a (1:a%dim, 1:a%dim) = a%a(1:a%dim, 1:a%dim)

  end subroutine CopyMatrix
!
!-------------------------------------------------!
! Commutator C=AB-BA !
!-------------------------------------------------!
  subroutine Commutator (c, a, b, io)
    character (len=*), parameter :: myname = 'commutator'
    type (matrixType), intent (in) :: a, b
    type (matrixType), intent (inout) :: c
    type (ioType), intent (in) :: io
    complex (kind=k_pr) :: alpha, beta
    integer :: i, j, ii, k
!
    if (c%created .and. (c%dim /= a%dim)) then
      call error ("Dimensions are different", myname, .true., io)
    end if
    if (b%dim /= a%dim) then
      call error ("Dimensions are different", myname, .true., io)
    end if
    if ( .not. c%created) then
      call CreateMatrix (c, a%dim, .true.)
    end if
    if (a%isSparse) then

      do j = 1, b%dim
        do ii = 1, a%nonZero
          i = a%indx (ii)
          k = a%jndx (ii)
          c%a (i, j) = c%a(i, j) + a%a(i, k) * b%a(k, j)
          c%a (j, k) = c%a(j, k) - b%a(j, i) * a%a(i, k)
        end do
      end do
    else if (b%isSparse) then
      do j = 1, b%dim
        do ii = 1, b%nonZero
          i = b%indx (ii)
          k = b%jndx (ii)
          c%a (i, j) = c%a(i, j) + a%a(i, k) * b%a(k, j)
          c%a (j, k) = c%a(j, k) - b%a(j, i) * a%a(i, k)
        end do
      end do
    else
      alpha = cmplx (1.0_k_pr, 0.0_k_pr, k_pr)
      beta = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
#ifdef MKL95
      call gemm (b%a, a%a, c%a, 'N', 'N', alpha, beta)
#else
      call zgemm ('N','N',b%dim,a%dim,c%dim,alpha,b%a,b%dim, a%a,a%dim,beta,c%a,c%dim)
#endif
      alpha = cmplx (1.0_k_pr, 0.0_k_pr, k_pr)
      beta = cmplx (-1.0_k_pr, 0.0_k_pr, k_pr)
#ifdef MKL95
      call gemm (a%a, b%a, c%a, 'N', 'N', alpha, beta)
#else
      call zgemm ('N','N',a%dim,b%dim,c%dim,alpha,a%a,a%dim,b%a,b%dim,beta, c%a,c%dim)
#endif
    end if
  end subroutine Commutator
!
!-------------------------------------------------!
! Calculate the trace of a matrix Product !
!-------------------------------------------------!
  function ProductTrace (a, b, io)
    character (len=*), parameter :: myname = 'ProductTrace'
    complex (kind=k_pr) :: ProductTrace
    type (matrixType), intent (in) :: a, b
    type (ioType), intent (in) :: io
    integer :: i, j, ii, ji
    complex (kind=k_pr) :: pt
!
    if (( .not. a%created) .or. ( .not. b%created)) then
      call error ("Cannot operate on unexistant matrix", myname, .true., io)
    end if
    ProductTrace = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
    pt = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
    if (a%isSparse) then
      do i = 1, a%nonZero
        ii = a%indx (i)
        ji = a%jndx (i)
        pt = pt + a%a(ii, ji) * b%a(ji, ii)
      end do
    else if (b%isSparse) then
      do i = 1, b%nonZero
        ii = b%indx (i)
        ji = b%jndx (i)
        pt = pt + a%a(ii, ji) * b%a(ji, ii)
      end do
    else
      do i = 1, a%dim
        do j = 1, a%dim
          pt = pt + a%a(i, j) * b%a(j, i)
        end do
      end do
    end if
    ProductTrace = pt
  end function ProductTrace
!> \brief computes matrix
!> \f$ {\mathbf C}=\alpha {\mathbf A}+\beta {\mathbf B}\f$
!> \author Alin M Elena
!> \date 09/11/07, 20:20:15
!> \param a,b type(matrixType) matrices to be substracted (b from a)
!> \param c type(matrixType) the result matrix
!> \param alpha, beta complex prefators for the matrices a,b
!> \param io type(ioType) i/o units
!
  subroutine MatrixCeaApbB (c, a, b, alpha, beta, io)
    character (len=*), parameter :: myname = "MatrixCeaApbB"
    type (matrixType), intent (in) :: a, b
    type (matrixType), intent (inout) :: c
    complex (k_pr), intent (in) :: alpha, beta
    type (ioType), intent (in) :: io
    integer :: i, j, k, l, m, n
    complex (k_pr) :: al, be
!
    if (c%created .and. (c%dim /= a%dim)) then
      call error ("Dimensions are different", myname, .true., io)
    end if
    if (b%dim /= a%dim) then
      call error ("Dimensions of the input matrices are different", myname, .true., io)
    end if
    if ( .not. c%created) then
      call error ("The result matrix does not exist!", myname, .true., io)
    end if
!
    if (c%isSparse .and. b%isSparse .and. a%isSparse) then
      if (a%nonZero == b%nonZero) then
        do i = 1, a%nonZero
          k = a%indx (i)
          l = a%jndx (i)
          call SpmPut (c, k, l, alpha*a%a(k, l)+beta*b%a(k, l))
          n = b%indx (i)
          m = b%jndx (i)
          call SpmPut (c, n, m, alpha*a%a(n, m)+beta*b%a(n, m))
        end do
      else if (a%nonZero < b%nonZero) then
        do i = 1, a%nonZero
          k = a%indx (i)
          l = a%jndx (i)
          call SpmPut (c, k, l, alpha*a%a(k, l)+beta*b%a(k, l))
          n = b%indx (i)
          m = b%jndx (i)
          call SpmPut (c, n, m, alpha*a%a(n, m)+beta*b%a(n, m))
        end do
        do i = a%nonZero + 1, b%nonZero
          k = b%indx (i)
          l = b%jndx (i)
          call SpmPut (c, k, l, beta*b%a(k, l))
        end do
      else
        do i = 1, b%nonZero
          k = a%indx (i)
          l = a%jndx (i)
          call SpmPut (c, k, l, alpha*a%a(k, l)+beta*b%a(k, l))
          n = b%indx (i)
          m = b%jndx (i)
          call SpmPut (c, n, m, alpha*a%a(n, m)+beta*b%a(n, m))
        end do
        do i = b%nonZero + 1, a%nonZero
          k = a%indx (i)
          l = a%jndx (i)
          call SpmPut (c, k, l, alpha*a%a(k, l))
        end do
      end if
    else
      do i = 1, c%dim
        do j = 1, c%dim
          c%a (i, j) = alpha * a%a(i, j) + beta * b%a(i, j)
        end do
      end do
    end if
  end subroutine MatrixCeaApbB

!> \brief perform a matrix-matrix operation using Hermitian matrices.
!> \details The operation is defined as C := alpha*A*conjg(A') + beta*C
!> this does (U sqrt(rho'))(U sqrt(rho'))* , only the upper triangle is calculated
!> unnocupied vectors are not multiplied not tr
!> the rest get filled after the herk call
!> \date 19th of December 2007
!> \param a, c the matrices
!> \param alpha, beta the prefactors
!> \param n dimension of the matrix c
  subroutine aastar (a, c, alpha, beta, n)
    complex (k_pr), intent (inout) :: a (:, :), c (:, :)
    real (k_pr), intent (in) :: alpha, beta
    integer, intent (in) :: n
    integer :: i
!
#ifdef MKL95
    call herk (a, c, 'U', 'N', alpha, beta)
#else
    call zherk('U', 'N',n,n,alpha,a,n,beta,c,n)
#endif
! this fills in the lower triangle
    do i = 1, n - 1
      c (i+1:n, i) = conjg (c(i, i+1:n))
    end do

  end subroutine aastar
!
end module m_LinearAlgebra
!-------------------------------------------------!
! end of Module Linear Algebra !
!-------------------------------------------------!
