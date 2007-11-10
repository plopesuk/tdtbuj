!>\brief This module defines the matrix type and
!> implements all linear algebra subroutines
!> \details by apropriate blas/lapack calls or special sparse 
!> capable subroutines.
!> \author Cristian G. Sanchez (Belfast)
!> \remarks
!> changed by Alin M Elena: the routines do not depend on global variables
!> changed to comply with the coding style document
!> \todo the calls to blas/lapack should be replaced with the f95 form
module m_LinearAlgebra
   use m_Constants
   use m_Useful
   use m_Types

   implicit none
   private

   public :: CreateMatrix
   public :: CreateSparseMatrix
   public :: ResetSparseMatrix
   public :: SpmPut
   public :: DestroyMatrix
   public :: SymmetrizeMatrix
   public :: RandomMatrix
   public :: RandomHermitianMatrix
   public :: IdentityMatrix
   public :: ZeroMatrix

   public :: CreateVector
   public :: DestroyVector
   public :: PrintVector
   public :: RandomVector
   public :: NormalizeVector

   public :: VectorNorm
   public :: MatrixTrace
   public :: ScalarProd
   public :: ProductTrace

   public :: ScalarTVector
   public :: ScalarTMatrix

   public :: DiagonalizeMatrix
   public :: FirstEigen
   public :: MatrixTVector
   public :: MatrixTMatrix
   public :: commutator
   public :: Anticommutator

   public :: ColumnFromMatrix
   public :: RowFromMatrix

   public :: TransposeMatrix
   public :: AdjointMatrix
   public :: CopyMatrix
   public :: CopySparseMatrix
   public :: CopyVector
   public :: MatrixAdd
   public :: MatrixCeaApbB
   public :: MatrixAddInPlace
   public :: DenseSparseDense

contains

  !-------------------------------------------------!
  ! Allocates memory and posibly Initializes the !
  ! matrix a !
  !-------------------------------------------------!

   subroutine CreateMatrix(a,m,zeroout)

      implicit none
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'CreateMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      logical,intent(in) :: zeroout
      integer,intent(in) :: m
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      a%isSparse = .false.
      a%dim = m

      allocate(a%a(m,m))

      if (zeroout) a%a = cmplx(0.0_k_pr,0.0_k_pr,k_pr)

      a%created = .true.

   end subroutine CreateMatrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Allocates memory and and Initializes the sparse !
  ! matrix a !
  !-------------------------------------------------!

   subroutine CreateSparseMatrix(a,m,zeroout)

      implicit none
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'CreateSparseMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      logical,intent(in) :: zeroout
      integer,intent(in) :: m
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      a%isSparse = .true.
      a%dim = m
      a%nonZero = 0
      allocate(a%a(m,m))
      allocate(a%indx(m*m))
      allocate(a%jndx(m*m))

      if (zeroout) a%a = cmplx(0.0_k_pr,0.0_k_pr,k_pr)

      a%created = .true.

   end subroutine CreateSparseMatrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Resets the index arrays in the sparse matrix a !
  !-------------------------------------------------!

   subroutine ResetSparseMatrix(a)

      implicit none
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'ResetSparseMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
    !--internal variables ----------------------------!
      integer :: i!,j
    !-------------------------------------------------!

      do i=1,a%nonZero
         a%a(a%indx(i),a%jndx(i)) = cmplx(0.0_k_pr,0.0_k_pr,k_pr)
      end do

      a%indx(1:a%nonZero) = 0
      a%jndx(1:a%nonZero) = 0

      a%nonZero = 0

   end subroutine ResetSparseMatrix


  !-------------------------------------------------!
  ! Puts a new element into the sparse matrix a !
  !-------------------------------------------------!

   subroutine SpmPut(a,i,j,aij)
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'SpmPut'
    !--subroutine parameters -------------------------!
      type(matrixType), intent(inout) :: a
      complex(k_pr), intent(in) :: aij
      integer,intent(in) :: i,j
    !--internal variables ----------------------------!
    !-------------------------------------------------!
      if(.not. IsIn(a,i,j)) then
         a%nonZero = a%nonZero + 1
         a%indx(a%nonZero) = i
         a%jndx(a%nonZero) = j
      endif
      a%a(i,j) = aij
 end subroutine SpmPut

!> \brief checks if a certain element is non zero in a sparse matrix
!> \author Alin M Elena
!> \date 09/11/07, 20:53:19
!> \param r,c integers row, column
!> \param a type(matrixType) the matrix
  logical function IsIn(a,r,c)
    integer, intent(in) :: r,c
    type(matrixType), intent(in) :: a
    integer :: t
    IsIn=.false.
    do t=1,a%nonZero
      if ((a%indx(t)==r).and.(a%jndx(t)==c)) then
        IsIn=.true.
        exit
      endif
    end do
  end function IsIn


  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Allocates memory and posibly Initializes the !
  ! vector v !
  !-------------------------------------------------!

   subroutine CreateVector(v,m,zeroout)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'CreateVector'
    !--subroutine parameters--------------------------!
      type(vectorType),intent(inout) :: v
      logical,intent(in) :: zeroout
      integer,intent(in) :: m
    !--internal variables-----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      
      v%dim = m

      allocate(v%v(m))

      if (zeroout) v%v = cmplx(0.0_k_pr,0.0_k_pr,k_pr)

      v%created = .true.

   end subroutine CreateVector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Deallocates memory in the matrix a !
  !-------------------------------------------------!

   subroutine DestroyMatrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'DestroyMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(ioType), intent(in) :: io
    !-------------------------------------------------!

      if (a%created) then
         deallocate(a%a)
         a%dim = 0
         a%created = .false.
         if (a%isSparse) then
            a%isSparse = .false.
            deallocate(a%indx)
            deallocate(a%jndx)
         endif
      else
         call error("Cannot destroy unexistant matrix",myname,.true.,io)
      endif

   end subroutine DestroyMatrix
! 
  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Deallocates memory in the vector v !
  !-------------------------------------------------!

   subroutine DestroyVector(v,io)
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'DestroyVector'
    !--subroutine parameters -------------------------!
      type(vectorType),intent(inout) :: v
      type(ioType), intent(in) :: io
    !-------------------------------------------------!

      if (v%created) then
         deallocate(v%v)
         v%dim = 0
         v%created = .false.
      else
         call error("Cannot destroy unexistant vector",myname,.true.,io)
      endif

   end subroutine DestroyVector

  !-------------------------------------------------!
  ! Prints Vector v !
  !-------------------------------------------------!

   subroutine PrintVector(v,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'PrintVector'
    !--subroutine parameters -------------------------!
      type(vectorType),intent(in) :: v
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!

      if (v%created) then
         write(io%uout,*) &
            "Output from PrintVector:"
         write(io%uout,*) &
            "dim = ",v%dim
         do i=1,v%dim
            write(io%uout,*) &
               "v(",i,") = ",v%v(i)
         end do
      else
         call error("Cannot Print unexistant vector",myname,.true.,io)
      endif

   end subroutine PrintVector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Symmetrizes Matrix a !
  !-------------------------------------------------!

   subroutine SymmetrizeMatrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'SymmetrizeMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j
      complex(kind=k_pr) :: temp
    !-------------------------------------------------!

      if (.not.a%created) then
        call error("Cannot symmetrize unexistant matrix",myname,.true.,io)
      endif

      do i=1,a%dim
         do j=i+1,a%dim
            temp = (a%a(i,j) + a%a(j,i)) / 2.0_k_pr
            a%a(i,j) = temp
            a%a(j,i) = temp
         end do
      end do

   end subroutine SymmetrizeMatrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Matrix a with random numbers !
  !-------------------------------------------------!

   subroutine RandomMatrix(a,u,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'RandomMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(ioType), intent(in) :: io
      real(k_pr),intent(inout) :: u(:)
    !--internal variables ----------------------------!
      integer :: i,j
    !-------------------------------------------------!

      if (.not.a%created) then
        call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      do i=1,a%dim
         do j=1,a%dim
            a%a(i,j) = cmplx(ranmar(u),ranmar(u),k_pr)
         end do
      end do

   end subroutine RandomMatrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Matrix a with random numbers !
  !-------------------------------------------------!

   subroutine RandomHermitianMatrix(a,u,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'RandomHermitianMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(ioType), intent(in) :: io
      real(k_pr),intent(inout) :: u(:)
    !--internal variables ----------------------------!
      integer :: i,j
    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      do i=1,a%dim
         do j=i+1,a%dim
            a%a(i,j) = cmplx(ranmar(u),ranmar(u),k_pr)
            a%a(j,i) = conjg(a%a(i,j))
         end do
      end do

      do i=1,a%dim
         a%a(i,i) = cmplx(ranmar(u),0.0_k_pr,k_pr)
      end do

   end subroutine RandomHermitianMatrix

  !-------------------------------------------------!
  ! k_zeroes matrix a !
  !-------------------------------------------------!

   subroutine ZeroMatrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'ZeroMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      a%a = cmplx(0.0_k_pr,0.0_k_pr,k_pr)

   end subroutine ZeroMatrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Vector v with random numbers !
  !-------------------------------------------------!

   subroutine RandomVector(v,u,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'RandomVector'
    !--subroutine parameters -------------------------!
      type(vectorType),intent(inout) :: v
      type(ioType), intent(in) :: io
      real(k_pr),intent(inout) :: u(:)
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!

      if (.not.v%created) then
         call error("Cannot fill an unexistant vector",myname,.true.,io)
      endif

      do i=1,v%dim
         v%v(i) = cmplx(ranmar(u),ranmar(u),k_pr)
      end do

   end subroutine RandomVector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Normalize Vector v !
  !-------------------------------------------------!

   subroutine NormalizeVector(v,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'NormalizeVector'
    !--subroutine parameters -------------------------!
      type(vectorType),intent(inout) :: v
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      real(kind=k_pr) :: norm
      integer :: i
    !-------------------------------------------------!

      if (.not.v%created) then
         call error("Cannot normalize an unexistant vector",myname,.true.,io)
      endif

      norm = 0.0_k_pr

      do i=1,v%dim
         norm = norm + real(v%v(i)*conjg(v%v(i)))
      end do

      norm = dsqrt(norm)

      v%v = v%v/norm

   end subroutine NormalizeVector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Matrix a with an identity matrix !
  !-------------------------------------------------!

   subroutine IdentityMatrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'IdentityMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i!,j
    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      a%a = cmplx(0.0_k_pr,0.0_k_pr,k_pr)

      do i=1,a%dim
         a%a(i,i) = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
      end do

   end subroutine IdentityMatrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Diagonalize a (hermitian) matrix !
  !-------------------------------------------------!

   subroutine DiagonalizeMatrix(a,c,lambda,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'DiagonalizeMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a,c
      real(kind=k_pr),intent(inout) :: lambda(:)
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!





      integer :: info,lwork
      complex(kind=k_pr), allocatable :: work(:)
      real(kind=k_pr), allocatable :: rwork(:)

    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot diagonalize unexistant matrix",myname,.true.,io)
      endif
      call CopyMatrix(c,a,io)
      allocate(rwork(3*a%dim-2))
      lwork = 5*a%dim-1
      allocate(work(lwork))
      call zheev('V','U', &
         a%dim,c%a,a%dim, &
         lambda, &
         work,lwork,rwork,info)
      deallocate(work)
      deallocate(rwork)
        if(info/=0) then
         call error("Diagonalization failed",myname,.true.,io)
      endif
   end subroutine DiagonalizeMatrix
  !-------------------------------------------------!


   subroutine FirstEigen(a,lambda,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'FirstEigen'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      real(kind=k_pr),intent(inout) :: lambda(:)
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!





      integer :: info,lwork
      complex(kind=k_pr), allocatable :: work(:)
      real(kind=k_pr), allocatable :: rwork(:)

    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot diagonalize unexistant matrix",myname,.true.,io)
      endif

      allocate(rwork(3*a%dim-2))
      lwork = 5*a%dim-1
      allocate(work(lwork))
      call zheev('N','U', &
         a%dim,a%a,a%dim, &
         lambda, &
         work,lwork,rwork,info)
      deallocate(work)
      deallocate(rwork)
      if(info/=0) then
         call error("Diagonalization failed",myname,.true.,io)
      endif
   end subroutine FirstEigen

  !-------------------------------------------------!
  ! calculate the norm of a vector !
  !-------------------------------------------------!
   function VectorNorm(v,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'VectorNorm'
    !--subroutine parameters -------------------------!
      real(kind=k_pr) :: VectorNorm
      type(vectorType),intent(in) :: v
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!
      if (.not.v%created) then
         call error("Cannot operate on unexistant vector",myname,.true.,io)
      endif
      VectorNorm= 0.0_k_pr
      do i=1,v%dim
         VectorNorm= VectorNorm + real(v%v(i)*conjg(v%v(i)))
      end do
      VectorNorm= dsqrt(VectorNorm)
   end function VectorNorm
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! calculate the trace of a matrix !
  !-------------------------------------------------!
   function MatrixTrace(a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'MatrixTrace'
    !--subroutine parameters -------------------------!
      complex(kind=k_pr) :: MatrixTrace
      type(matrixType),intent(in) :: a
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!
      if (.not.a%created) then
         call error("Cannot operate on unexistant matrix",myname,.true.,io)
      endif
      MatrixTrace= cmplx(0.0_k_pr,0.0_k_pr,k_pr)
      do i=1,a%dim
         MatrixTrace= MatrixTrace + a%a(i,i)
      end do
   end function MatrixTrace
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! calculate the scalar Product of two vectors !
  !-------------------------------------------------!
   function ScalarProd(va,vb,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'ScalarProd'
    !--subroutine parameters -------------------------!
      complex(kind=k_pr) :: ScalarProd
      type(vectorType),intent(in) :: va,vb
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!
      if ((.not.va%created).or.(.not.vb%created)) then
         call error("Cannot operate on non existant vector",myname,.true.,io)
      endif
      if (va%dim/=vb%dim) then
         call error("Cannot operate on vectors of different dimension",myname,.true.,io)
      endif
      ScalarProd= cmplx(0.0_k_pr,0.0_k_pr,k_pr)
      do i=1,va%dim
         ScalarProd= ScalarProd + conjg(va%v(i))*vb%v(i)
      end do
   end function ScalarProd
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Multiply scalar times matrix !
  !-------------------------------------------------!
   subroutine ScalarTMatrix(s,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'ScalarTMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      complex(kind=k_pr),intent(in) :: s
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!
      if (.not.a%created) then
         call error("Cannot operate on an unexistant matrix",myname,.true.,io)
      endif
      a%a = s*a%a
   end subroutine ScalarTMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Multiply scalar times vector !
  !-------------------------------------------------!
   subroutine ScalarTVector(s,v,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'ScalarTVector'
    !--subroutine parameters -------------------------!
      type(vectorType),intent(inout) :: v
      complex(kind=k_pr),intent(in) :: s
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i
    !-------------------------------------------------!
      if (.not.v%created) then
        call error("Cannot operate on an non existant vector",myname,.true.,io)
      endif
      v%v = s*v%v
   end subroutine ScalarTVector
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Convert a certain column from a matrix into a !
  ! vector !
  !-------------------------------------------------!
   subroutine ColumnFromMatrix(m,v,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'ColumnFromMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(vectorType),intent(inout) :: v
      integer,intent(in) :: m
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    !integer :: i
    !-------------------------------------------------!
      if (.not.a%created) then
         call error("Cannot operate on an unexistant matrix",myname,.true.,io)
      endif
      if(v%created.and.(v%dim/=a%dim)) then
        call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.v%created) then
        call CreateVector(v,a%dim,.false.)
      endif
      v%v = a%a(:,m)
   end subroutine ColumnFromMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Convert a certain row from a matrix into a !
  ! vector !
  !-------------------------------------------------!
   subroutine RowFromMatrix(m,v,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'RowFromMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(inout) :: a
      type(vectorType),intent(inout) :: v
      integer,intent(in) :: m
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i
    !-------------------------------------------------!
      if (.not.a%created) then
         call error("Cannot operate on an unexistant matrix",myname,.true.,io)
      endif
      if(v%created.and.(v%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.v%created) then
         call CreateVector(v,a%dim,.false.)
      endif
      v%v = a%a(m,:)
   end subroutine RowFromMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Tranpose Matrix a, put result in b !
  !-------------------------------------------------!
   subroutine TransposeMatrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'TransposeMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a
      type(matrixType),intent(inout) :: b
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j
    ! complex(kind=k_pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call CreateMatrix(b,a%dim,.false.)
      endif
      do i=1,a%dim
         do j=1,a%dim
            b%a(i,j) = a%a(j,i)
         end do
      end do
   end subroutine TransposeMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Adjoint of Matrix a, put result in b !
  !-------------------------------------------------!
   subroutine AdjointMatrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'AdjointMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a
      type(matrixType),intent(inout) :: b
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j
    ! complex(kind=k_pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call CreateMatrix(b,a%dim,.false.)
      endif

      do i=1,a%dim
         do j=1,a%dim
            b%a(i,j) = conjg(a%a(j,i))
         end do
      end do

   end subroutine AdjointMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Copy Matrix a into b !
  !-------------------------------------------------!
   subroutine CopyMatrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'CopyMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a
      type(matrixType),intent(inout) :: b
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    ! complex(kind=k_pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call CreateMatrix(b,a%dim,.false.)
      endif
      b%a = a%a
   end subroutine CopyMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Copy Sparse Matrix a into b !
  !-------------------------------------------------!
   subroutine CopySparseMatrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'CopySparseMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a
      type(matrixType),intent(inout) :: b
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    ! complex(kind=k_pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call CreateSparseMatrix(b,a%dim,.true.)
      endif
      b%nonZero=a%nonZero
      b%indx=a%indx
      b%jndx=a%jndx
      b%a = a%a
   end subroutine CopySparseMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Copy vector a into b !
  !-------------------------------------------------!
   subroutine CopyVector(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'CopyVector'
    !--subroutine parameters -------------------------!
      type(vectorType),intent(in) :: a
      type(vectorType),intent(inout) :: b
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    ! complex(kind=k_pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call CreateVector(b,a%dim,.false.)
      endif
      b%v = a%v
   end subroutine CopyVector
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix Product C=AB !
  !-------------------------------------------------!
   subroutine MatrixTMatrix(c,a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'MatrixTMatrix'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a,b
      type(matrixType),intent(inout) :: c
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      complex(kind=k_pr) :: alpha, beta
      integer :: i,j,k,ii
    !-------------------------------------------------!
      if (c%created.and.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.c%created) then
         call CreateMatrix(c,a%dim,.true.)
      endif
      if (a%isSparse) then

         do j=1,b%dim
            do ii=1,a%nonZero
               i = a%indx(ii)
               k = a%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
            enddo
         enddo

      elseif (b%isSparse) then

         do i=1,b%dim
            do ii=1,b%nonZero
               k = b%indx(ii)
               j = b%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
            enddo
         enddo

      else
         alpha = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
         beta = cmplx(0.0_k_pr,0.0_k_pr,k_pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            a%a,a%dim, &
            b%a,b%dim, &
            beta ,c%a,c%dim)
      endif
   end subroutine MatrixTMatrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix times vector Product y=Ax !
  !-------------------------------------------------!
   subroutine MatrixTVector(y,a,x,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'MatrixTVector'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a
      type(vectorType),intent(in) :: x
      type(vectorType),intent(inout) :: y
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
      complex(kind=k_pr) :: alpha, beta
    !-------------------------------------------------!
      if (y%created.and.(y%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (x%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.y%created) then
         call CreateVector(y,a%dim,.true.)
      endif
      alpha = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
      beta = cmplx(0.0_k_pr,0.0_k_pr,k_pr)
      call zgemv('N',a%dim,a%dim, &
         alpha , &
         a%a,a%dim, &
         x%v,1, &
         beta ,y%v,1)
   end subroutine MatrixTVector
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Commutator C=AB-BA !
  !-------------------------------------------------!
   subroutine commutator(c,a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'commutator'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a,b
      type(matrixType),intent(inout) :: c
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      complex(kind=k_pr) :: alpha, beta
      integer :: i,j,ii,k
    !-------------------------------------------------!
      if (c%created.and.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.c%created) then
         call CreateMatrix(c,a%dim,.true.)
      endif
      if (a%isSparse) then

         do j=1,b%dim
            do ii=1,a%nonZero
               i = a%indx(ii)
               k = a%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) - b%a(j,i)*a%a(i,k)
            enddo
         enddo

      elseif (b%isSparse) then

         do j=1,b%dim
            do ii=1,b%nonZero
               i = b%indx(ii)
               k = b%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) - b%a(j,i)*a%a(i,k)
            enddo
         enddo

      else
       ! both matrices are dense
         alpha = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
         beta = cmplx(0.0_k_pr,0.0_k_pr,k_pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            b%a,b%dim, &
            a%a,a%dim, &
            beta ,c%a,c%dim)
         alpha = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
         beta = cmplx(-1.0_k_pr,0.0_k_pr,k_pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            a%a,a%dim, &
            b%a,b%dim, &
            beta ,c%a,c%dim)
      endif
   end subroutine commutator
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Anticommutator C=AB+BA !
  !-------------------------------------------------!
   subroutine Anticommutator(c,a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'Anticommutator'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a,b
      type(matrixType),intent(inout) :: c
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      complex(kind=k_pr) :: alpha, beta
      integer :: i,j,k,ii
    !-------------------------------------------------!
      if (c%created.and.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.c%created) then
         call CreateMatrix(c,a%dim,.true.)
      endif
      if (a%isSparse) then

         do j=1,b%dim
            do ii=1,a%nonZero
               i = a%indx(ii)
               k = a%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) + b%a(j,i)*a%a(i,k)
            enddo
         enddo

      elseif (b%isSparse) then

         do j=1,b%dim
            do ii=1,b%nonZero
               i = b%indx(ii)
               k = b%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) + b%a(j,i)*a%a(i,k)
            enddo
         enddo

      else
       ! both matrices are dense
         alpha = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
         beta = cmplx(0.0_k_pr,0.0_k_pr,k_pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            b%a,b%dim, &
            a%a,a%dim, &
            beta ,c%a,c%dim)
         alpha = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
         beta = cmplx(1.0_k_pr,0.0_k_pr,k_pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            a%a,a%dim, &
            b%a,b%dim, &
            beta ,c%a,c%dim)
      endif
   end subroutine Anticommutator
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Calculate the trace of a matrix Product !
  !-------------------------------------------------!
   function ProductTrace(a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'ProductTrace'
    !--subroutine parameters -------------------------!
      complex(kind=k_pr) :: ProductTrace
      type(matrixType),intent(in) :: a,b
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j,ii,ji
      complex(kind=k_pr) :: pt
    !-------------------------------------------------!
      if ((.not.a%created).or.(.not.b%created)) then
         call error("Cannot operate on unexistant matrix",myname,.true.,io)
      endif
      ProductTrace= cmplx(0.0_k_pr,0.0_k_pr,k_pr)
      pt = cmplx(0.0_k_pr,0.0_k_pr,k_pr)
      if (a%isSparse) then
         do i=1,a%nonZero
            ii = a%indx(i)
            ji = a%jndx(i)
            pt = pt + a%a(ii,ji)*b%a(ji,ii)
         enddo
      elseif (b%isSparse) then
         do i=1,b%nonZero
            ii = b%indx(i)
            ji = b%jndx(i)
            pt = pt + a%a(ii,ji)*b%a(ji,ii)
         enddo
      else
         do i=1,a%dim
            do j=1,a%dim
               pt = pt + a%a(i,j)*b%a(j,i)
            enddo
         end do
      endif
      ProductTrace= pt
   end function ProductTrace
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix Addition C=A+B !
  !-------------------------------------------------!
   subroutine MatrixAdd(c,a,b,alpha,beta,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'MatrixAdd'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a,b
      type(matrixType),intent(inout) :: c
      complex(k_pr), intent(in), optional :: alpha, beta
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!
      if (c%created.and.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.c%created) then
         call CreateMatrix(c,a%dim,.true.)
      endif
      if (present(alpha).and.present(beta)) then
         c%a = beta*b%a + alpha*a%a
      elseif (present(alpha).and.(.not.present(beta))) then
         c%a = b%a + alpha*a%a
      elseif (present(beta).and.(.not.present(alpha))) then
         c%a = beta*b%a + a%a
      else
         c%a = b%a + a%a
      endif
   end subroutine MatrixAdd

!> \brief computes matrix 
!> \f$ {\mathbf C}=\alpha {\mathbf A}+\beta {\mathbf B}\f$
!> \author Alin M Elena
!> \date 09/11/07, 20:20:15
!> \param a,b type(matrixType) matrices to be substracted (b from a)
!> \param c type(matrixType) the result matrix
!> \param alpha, beta comples prefators for the matrices a,b, if any is absent is assumed 1
!> \param io type(ioType) i/o units

  subroutine MatrixCeaApbB(c,a,b,alpha,beta,io)
    character(len=*), parameter :: myname = "MatrixCeaApbB"
    type(matrixType),intent(in) :: a,b
    type(matrixType),intent(inout) :: c
    complex(k_pr), intent(in):: alpha, beta
    type(ioType), intent(in) :: io
    integer :: i,j,k,l,m,n
    complex(k_pr) :: al,be

    if (c%created.and.(c%dim/=a%dim)) then
        call error("Dimensions are different",myname,.true.,io)
    endif
    if (b%dim/=a%dim) then
        call error("Dimensions are different",myname,.true.,io)
    endif
    if (.not.c%created) then
        call error("The result matrix does not exist!",myname,.true.,io)
    endif

    if (c%IsSparse .and. b%IsSparse .and. a%IsSparse) then
      if (a%nonZero==b%nonZero) then
        do i=1,a%nonZero
          k = a%indx(i)
          l = a%jndx(i)
          call SpmPut(c,k,l,alpha*a%a(k,l)+beta*b%a(k,l))
          n = b%indx(i)
          m = b%jndx(i)
          call SpmPut(c,n,m,alpha*a%a(n,m)+beta*b%a(n,m))
        enddo
      elseif(a%nonZero<b%nonZero) then
        do i=1,a%nonZero
          k = a%indx(i)
          l = a%jndx(i)
          call SpmPut(c,k,l,alpha*a%a(k,l)+beta*b%a(k,l))
          n = b%indx(i)
          m = b%jndx(i)
          call SpmPut(c,n,m,alpha*a%a(n,m)+beta*b%a(n,m))
        enddo
        do i=a%nonZero+1,b%nonZero
          k = b%indx(i)
          l = b%jndx(i)
          call SpmPut(c,k,l,beta*b%a(k,l))
        enddo
      else
        do i=1,b%nonZero
          k = a%indx(i)
          l = a%jndx(i)
          call SpmPut(c,k,l,alpha*a%a(k,l)+beta*b%a(k,l))
          n = b%indx(i)
          m = b%jndx(i)
          call SpmPut(c,n,m,alpha*a%a(n,m)+beta*b%a(n,m))
        enddo
        do i=b%nonZero+1,a%nonZero
          k = a%indx(i)
          l = a%jndx(i)
          call SpmPut(c,k,l,alpha*a%a(k,l))
        enddo
      endif
    else
      c%a = alpha*a%a + beta*b%a 
    endif
  end subroutine MatrixCeaApbB


!-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix Addition (in place) A=A+B !
  !-------------------------------------------------!
   subroutine MatrixAddInPlace(a,b,alpha,beta,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'MatrixAddInPlace'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: b
      type(matrixType),intent(inout) :: a
      complex(k_pr), intent(in), optional :: alpha, beta
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (present(alpha).and.present(beta)) then
         a%a = beta*b%a + alpha*a%a
      elseif (present(alpha).and.(.not.present(beta))) then
         a%a = b%a + alpha*a%a
      elseif (present(beta).and.(.not.present(alpha))) then
         a%a = beta*b%a + a%a
      else
         a%a = b%a + a%a
      endif
   end subroutine MatrixAddInPlace
  !-------------------------------------------------!
  ! dense * sparse * dense triple Product !
  !-------------------------------------------------!
   subroutine DenseSparseDense(p,a,b,c,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'DenseSparseDense'
    !--subroutine parameters -------------------------!
      type(matrixType),intent(in) :: a,b,c
      type(matrixType),intent(inout) :: p
      type(ioType), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j,ii,k,l
    !-------------------------------------------------!
      if (p%created.and.(p%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if ((b%dim/=a%dim).or.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.p%created) then
         call CreateMatrix(p,a%dim,.true.)
      endif

      do i=1,b%dim
         do j=1,b%dim
            p%a(i,j) = 0.0_k_pr
            do ii=1,b%nonZero
               k = b%indx(ii)
               l = b%jndx(ii)
               p%a(i,j) = p%a(i,j) + a%a(i,k)*b%a(k,l)*c%a(l,j)
            enddo
         enddo
      enddo

   end subroutine DenseSparseDense
end module m_LinearAlgebra
!-------------------------------------------------!
! end of Module Linear Algebra !
!-------------------------------------------------!s