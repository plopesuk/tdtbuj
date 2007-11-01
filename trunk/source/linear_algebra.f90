!>\brief This module defines the matrix types and
!> implements all linear algebra subroutines
!> \details by apropriate blas/lapack calls or special sparse 
!> capable subroutines.
!> \author Cristian G. Sanchez (Belfast)
!> \remarks
!> changed by Alin M Elena: the routines do not depend on global variables
!> \todo the calls to blas/lapack should be replaced with the f95 form

module linear_algebra

   use constants
   use useful
   use types

   implicit none
   private

  !-------------------------------------------------!
  ! Public types !
  !-------------------------------------------------!


  !-------------------------------------------------!
  ! Public subroutines !
  !-------------------------------------------------!

   public :: create_matrix
   public :: create_sparse_matrix
   public :: reset_sparse_matrix
   public :: spm_put
   public :: destroy_matrix
   public :: print_matrix
   public :: symmetrize_matrix
   public :: random_matrix
   public :: random_hermitian_matrix
   public :: identity_matrix
   public :: zero_matrix

   public :: create_vector
   public :: destroy_vector
   public :: print_vector
   public :: random_vector
   public :: normalize_vector

   public :: vector_norm
   public :: matrix_trace
   public :: scalar_prod
   public :: product_trace

   public :: scalar_t_vector
   public :: scalar_t_matrix

   public :: diagonalize_matrix
   public :: first_eigen
   public :: matrix_t_vector
   public :: matrix_t_matrix
   public :: commutator
   public :: anticommutator

   public :: column_from_matrix
   public :: row_from_matrix

   public :: transpose_matrix
   public :: adjoint_matrix
   public :: copy_matrix
   public :: copy_sparse_matrix
   public :: copy_vector
   public :: matrix_add
   public :: matrix_add_inplace
   public :: dense_sparse_dense

contains

  !-------------------------------------------------!
  ! Allocates memory and posibly initializes the !
  ! matrix a !
  !-------------------------------------------------!

   subroutine create_matrix(a,m,zeroout)

      implicit none
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'create_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      logical,intent(in) :: zeroout
      integer,intent(in) :: m
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      a%is_sparse = .false.
      a%dim = m

      allocate(a%a(m,m))

      if (zeroout) a%a = cmplx(0.0_pr,0.0_pr,pr)

      a%created = .true.

   end subroutine create_matrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Allocates memory and and initializes the sparse !
  ! matrix a !
  !-------------------------------------------------!

   subroutine create_sparse_matrix(a,m,zeroout)

      implicit none
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'create_sparse_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      logical,intent(in) :: zeroout
      integer,intent(in) :: m
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      a%is_sparse = .true.
      a%dim = m
      a%nonzero = 0
      allocate(a%a(m,m))
      allocate(a%indx(m*m))
      allocate(a%jndx(m*m))

      if (zeroout) a%a = cmplx(0.0_pr,0.0_pr,pr)

      a%created = .true.

   end subroutine create_sparse_matrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Resets the index arrays in the sparse matrix a !
  !-------------------------------------------------!

   subroutine reset_sparse_matrix(a)

      implicit none
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'reset_sparse_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
    !--internal variables ----------------------------!
      integer :: i!,j
    !-------------------------------------------------!

      do i=1,a%nonzero
         a%a(a%indx(i),a%jndx(i)) = cmplx(0.0_pr,0.0_pr,pr)
      end do

      a%indx(1:a%nonzero) = 0
      a%jndx(1:a%nonzero) = 0

      a%nonzero = 0

   end subroutine reset_sparse_matrix


  !-------------------------------------------------!
  ! Puts a new element into the sparse matrix a !
  !-------------------------------------------------!

   subroutine spm_put(a,i,j,aij)
    !--subroutine name -------------------------------!
      character(len=*), parameter :: myname = 'spm_put'
    !--subroutine parameters -------------------------!
      type(matrix), intent(inout) :: a
      complex(pr), intent(in) :: aij
      integer,intent(in) :: i,j
    !--internal variables ----------------------------!
    !-------------------------------------------------!
      if(.not. is_in(i,j)) then
         a%nonzero = a%nonzero + 1
         a%indx(a%nonzero) = i
         a%jndx(a%nonzero) = j
      endif
      a%a(i,j) = aij
          
contains
          
      function is_in(r,c)
         integer, intent(in) :: r,c
         logical :: is_in 
         integer :: t
          
         is_in=.false.
         do t=1,a%nonzero
            if (a%indx(t)==r) then
               if(a%jndx(t)==c) then
                  is_in=.true.
                  exit
               endif
            endif
         end do
                   
      end function is_in
                   
   end subroutine spm_put


  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Allocates memory and posibly initializes the !
  ! vector v !
  !-------------------------------------------------!

   subroutine create_vector(v,m,zeroout)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'create_vector'
    !--subroutine parameters--------------------------!
      type(vector),intent(inout) :: v
      logical,intent(in) :: zeroout
      integer,intent(in) :: m
    !--internal variables-----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      
      v%dim = m

      allocate(v%v(m))

      if (zeroout) v%v = cmplx(0.0_pr,0.0_pr,pr)

      v%created = .true.

   end subroutine create_vector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Deallocates memory in the matrix a !
  !-------------------------------------------------!

   subroutine destroy_matrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'destroy_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(io_type), intent(in) :: io
    !-------------------------------------------------!

      if (a%created) then
         deallocate(a%a)
         a%dim = 0
         a%created = .false.
         if (a%is_sparse) then
            a%is_sparse = .false.
            deallocate(a%indx)
            deallocate(a%jndx)
         endif
      else
         call error("Cannot destroy unexistant matrix",myname,.true.,io)
      endif

   end subroutine destroy_matrix
! 
  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Deallocates memory in the vector v !
  !-------------------------------------------------!

   subroutine destroy_vector(v,io)
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'destroy_vector'
    !--subroutine parameters -------------------------!
      type(vector),intent(inout) :: v
      type(io_type), intent(in) :: io
    !-------------------------------------------------!

      if (v%created) then
         deallocate(v%v)
         v%dim = 0
         v%created = .false.
      else
         call error("Cannot destroy unexistant vector",myname,.true.,io)
      endif

   end subroutine destroy_vector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Prints Matrix a !
  !-------------------------------------------------!

   subroutine print_matrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'print_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j
    !-------------------------------------------------!

      if (a%created) then
         write(io%uout,*) &
            "Output from print_matrix:"
         write(io%uout,*) &
            "dim = ",a%dim
         do i=1,a%dim
            write(io%uout,'(a,i5,a,1500f9.4)') &
               'row ',i,'=',(real(a%a(i,j)),j=1,a%dim)
         enddo
         do i=1,a%dim
            write(io%uout,'(a,i5,a,1500f9.4)') &
               'row ',i,'=',(aimag(a%a(i,j)),j=1,a%dim)
         enddo
      else
         call error("Cannot print unexistant matrix",myname,.true.,io)
      endif

   end subroutine print_matrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Prints Vector v !
  !-------------------------------------------------!

   subroutine print_vector(v,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'print_vector'
    !--subroutine parameters -------------------------!
      type(vector),intent(in) :: v
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!

      if (v%created) then
         write(io%uout,*) &
            "Output from print_vector:"
         write(io%uout,*) &
            "dim = ",v%dim
         do i=1,v%dim
            write(io%uout,*) &
               "v(",i,") = ",v%v(i)
         end do
      else
         call error("Cannot print unexistant vector",myname,.true.,io)
      endif

   end subroutine print_vector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Symmetrizes Matrix a !
  !-------------------------------------------------!

   subroutine symmetrize_matrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'symmetrize_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j
      complex(kind=pr) :: temp
    !-------------------------------------------------!

      if (.not.a%created) then
        call error("Cannot symmetrize unexistant matrix",myname,.true.,io)
      endif

      do i=1,a%dim
         do j=i+1,a%dim
            temp = (a%a(i,j) + a%a(j,i)) / 2.0_pr
            a%a(i,j) = temp
            a%a(j,i) = temp
         end do
      end do

   end subroutine symmetrize_matrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Matrix a with random numbers !
  !-------------------------------------------------!

   subroutine random_matrix(a,u,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'random_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(io_type), intent(in) :: io
      real(pr),intent(inout) :: u(:)
    !--internal variables ----------------------------!
      integer :: i,j
    !-------------------------------------------------!

      if (.not.a%created) then
        call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      do i=1,a%dim
         do j=1,a%dim
            a%a(i,j) = cmplx(ranmar(u),ranmar(u),pr)
         end do
      end do

   end subroutine random_matrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Matrix a with random numbers !
  !-------------------------------------------------!

   subroutine random_hermitian_matrix(a,u,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'random_hermitian_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(io_type), intent(in) :: io
      real(pr),intent(inout) :: u(:)
    !--internal variables ----------------------------!
      integer :: i,j
    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      do i=1,a%dim
         do j=i+1,a%dim
            a%a(i,j) = cmplx(ranmar(u),ranmar(u),pr)
            a%a(j,i) = conjg(a%a(i,j))
         end do
      end do

      do i=1,a%dim
         a%a(i,i) = cmplx(ranmar(u),0.0_pr,pr)
      end do

   end subroutine random_hermitian_matrix

  !-------------------------------------------------!
  ! zeroes matrix a !
  !-------------------------------------------------!

   subroutine zero_matrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'zero_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      a%a = cmplx(0.0_pr,0.0_pr,pr)

   end subroutine zero_matrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Vector v with random numbers !
  !-------------------------------------------------!

   subroutine random_vector(v,u,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'random_vector'
    !--subroutine parameters -------------------------!
      type(vector),intent(inout) :: v
      type(io_type), intent(in) :: io
      real(pr),intent(inout) :: u(:)
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!

      if (.not.v%created) then
         call error("Cannot fill an unexistant vector",myname,.true.,io)
      endif

      do i=1,v%dim
         v%v(i) = cmplx(ranmar(u),ranmar(u),pr)
      end do

   end subroutine random_vector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Normalize Vector v !
  !-------------------------------------------------!

   subroutine normalize_vector(v,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'normalize_vector'
    !--subroutine parameters -------------------------!
      type(vector),intent(inout) :: v
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      real(kind=pr) :: norm
      integer :: i
    !-------------------------------------------------!

      if (.not.v%created) then
         call error("Cannot normalize an unexistant vector",myname,.true.,io)
      endif

      norm = 0.0_pr

      do i=1,v%dim
         norm = norm + real(v%v(i)*conjg(v%v(i)))
      end do

      norm = dsqrt(norm)

      v%v = v%v/norm

   end subroutine normalize_vector

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Fill Matrix a with an identity matrix !
  !-------------------------------------------------!

   subroutine identity_matrix(a,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'identity_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i!,j
    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot fill an unexistant matrix",myname,.true.,io)
      endif

      a%a = cmplx(0.0_pr,0.0_pr,pr)

      do i=1,a%dim
         a%a(i,i) = cmplx(1.0_pr,0.0_pr,pr)
      end do

   end subroutine identity_matrix

  !-------------------------------------------------!

  !-------------------------------------------------!
  ! Diagonalize a (hermitian) matrix !
  !-------------------------------------------------!

   subroutine diagonalize_matrix(a,c,lambda,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'diagonalize_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a,c
      real(kind=pr),intent(inout) :: lambda(:)
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!





      integer :: info,lwork
      complex(kind=pr), allocatable :: work(:)
      real(kind=pr), allocatable :: rwork(:)

    !-------------------------------------------------!

      if (.not.a%created) then
         call error("Cannot diagonalize unexistant matrix",myname,.true.,io)
      endif
      call copy_matrix(c,a,io)
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
   end subroutine diagonalize_matrix
  !-------------------------------------------------!


   subroutine first_eigen(a,lambda,io)

      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'first_eigen'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      real(kind=pr),intent(inout) :: lambda(:)
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!





      integer :: info,lwork
      complex(kind=pr), allocatable :: work(:)
      real(kind=pr), allocatable :: rwork(:)

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
   end subroutine first_eigen

  !-------------------------------------------------!
  ! calculate the norm of a vector !
  !-------------------------------------------------!
   function vector_norm(v,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'vector_norm'
    !--subroutine parameters -------------------------!
      real(kind=pr) :: vector_norm
      type(vector),intent(in) :: v
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!
      if (.not.v%created) then
         call error("Cannot operate on unexistant vector",myname,.true.,io)
      endif
      vector_norm= 0.0_pr
      do i=1,v%dim
         vector_norm= vector_norm + real(v%v(i)*conjg(v%v(i)))
      end do
      vector_norm= dsqrt(vector_norm)
   end function vector_norm
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! calculate the trace of a matrix !
  !-------------------------------------------------!
   function matrix_trace(a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'matrix_trace'
    !--subroutine parameters -------------------------!
      complex(kind=pr) :: matrix_trace
      type(matrix),intent(in) :: a
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!
      if (.not.a%created) then
         call error("Cannot operate on unexistant matrix",myname,.true.,io)
      endif
      matrix_trace= cmplx(0.0_pr,0.0_pr,pr)
      do i=1,a%dim
         matrix_trace= matrix_trace + a%a(i,i)
      end do
   end function matrix_trace
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! calculate the scalar product of two vectors !
  !-------------------------------------------------!
   function scalar_prod(va,vb,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'scalar_prod'
    !--subroutine parameters -------------------------!
      complex(kind=pr) :: scalar_prod
      type(vector),intent(in) :: va,vb
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i
    !-------------------------------------------------!
      if ((.not.va%created).or.(.not.vb%created)) then
         call error("Cannot operate on non existant vector",myname,.true.,io)
      endif
      if (va%dim/=vb%dim) then
         call error("Cannot operate on vectors of different dimension",myname,.true.,io)
      endif
      scalar_prod= cmplx(0.0_pr,0.0_pr,pr)
      do i=1,va%dim
         scalar_prod= scalar_prod + conjg(va%v(i))*vb%v(i)
      end do
   end function scalar_prod
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Multiply scalar times matrix !
  !-------------------------------------------------!
   subroutine scalar_t_matrix(s,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'scalar_t_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      complex(kind=pr),intent(in) :: s
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    !-------------------------------------------------!
      if (.not.a%created) then
         call error("Cannot operate on an unexistant matrix",myname,.true.,io)
      endif
      a%a = s*a%a
   end subroutine scalar_t_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Multiply scalar times vector !
  !-------------------------------------------------!
   subroutine scalar_t_vector(s,v,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'scalar_t_vector'
    !--subroutine parameters -------------------------!
      type(vector),intent(inout) :: v
      complex(kind=pr),intent(in) :: s
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i
    !-------------------------------------------------!
      if (.not.v%created) then
        call error("Cannot operate on an non existant vector",myname,.true.,io)
      endif
      v%v = s*v%v
   end subroutine scalar_t_vector
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Convert a certain column from a matrix into a !
  ! vector !
  !-------------------------------------------------!
   subroutine column_from_matrix(m,v,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'column_from_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(vector),intent(inout) :: v
      integer,intent(in) :: m
      type(io_type), intent(in) :: io
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
        call create_vector(v,a%dim,.false.)
      endif
      v%v = a%a(:,m)
   end subroutine column_from_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Convert a certain row from a matrix into a !
  ! vector !
  !-------------------------------------------------!
   subroutine row_from_matrix(m,v,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'row_from_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(inout) :: a
      type(vector),intent(inout) :: v
      integer,intent(in) :: m
      type(io_type), intent(in) :: io
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
         call create_vector(v,a%dim,.false.)
      endif
      v%v = a%a(m,:)
   end subroutine row_from_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Tranpose Matrix a, put result in b !
  !-------------------------------------------------!
   subroutine transpose_matrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'transpose_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a
      type(matrix),intent(inout) :: b
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j
    ! complex(kind=pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call create_matrix(b,a%dim,.false.)
      endif
      do i=1,a%dim
         do j=1,a%dim
            b%a(i,j) = a%a(j,i)
         end do
      end do
   end subroutine transpose_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Adjoint of Matrix a, put result in b !
  !-------------------------------------------------!
   subroutine adjoint_matrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'adjoint_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a
      type(matrix),intent(inout) :: b
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j
    ! complex(kind=pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call create_matrix(b,a%dim,.false.)
      endif

      do i=1,a%dim
         do j=1,a%dim
            b%a(i,j) = conjg(a%a(j,i))
         end do
      end do

   end subroutine adjoint_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Copy Matrix a into b !
  !-------------------------------------------------!
   subroutine copy_matrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'copy_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a
      type(matrix),intent(inout) :: b
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    ! complex(kind=pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call create_matrix(b,a%dim,.false.)
      endif
      b%a = a%a
   end subroutine copy_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Copy Sparse Matrix a into b !
  !-------------------------------------------------!
   subroutine copy_sparse_matrix(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'copy_sparse_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a
      type(matrix),intent(inout) :: b
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    ! complex(kind=pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call create_sparse_matrix(b,a%dim,.true.)
      endif
      b%nonzero=a%nonzero
      b%indx=a%indx
      b%jndx=a%jndx
      b%a = a%a
   end subroutine copy_sparse_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Copy vector a into b !
  !-------------------------------------------------!
   subroutine copy_vector(b,a,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'copy_vector'
    !--subroutine parameters -------------------------!
      type(vector),intent(in) :: a
      type(vector),intent(inout) :: b
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
    ! complex(kind=pr) :: temp
    !-------------------------------------------------!
      if(b%created.and.(b%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.b%created) then
         call create_vector(b,a%dim,.false.)
      endif
      b%v = a%v
   end subroutine copy_vector
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix Product C=AB !
  !-------------------------------------------------!
   subroutine matrix_t_matrix(c,a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'matrix_t_matrix'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a,b
      type(matrix),intent(inout) :: c
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      complex(kind=pr) :: alpha, beta
      integer :: i,j,k,ii
    !-------------------------------------------------!
      if (c%created.and.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.c%created) then
         call create_matrix(c,a%dim,.true.)
      endif
      if (a%is_sparse) then

         do j=1,b%dim
            do ii=1,a%nonzero
               i = a%indx(ii)
               k = a%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
            enddo
         enddo

      elseif (b%is_sparse) then

         do i=1,b%dim
            do ii=1,b%nonzero
               k = b%indx(ii)
               j = b%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
            enddo
         enddo

      else
         alpha = cmplx(1.0_pr,0.0_pr,pr)
         beta = cmplx(0.0_pr,0.0_pr,pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            a%a,a%dim, &
            b%a,b%dim, &
            beta ,c%a,c%dim)
      endif
   end subroutine matrix_t_matrix
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix times vector Product y=Ax !
  !-------------------------------------------------!
   subroutine matrix_t_vector(y,a,x,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'matrix_t_vector'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a
      type(vector),intent(in) :: x
      type(vector),intent(inout) :: y
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
    ! integer :: i,j
      complex(kind=pr) :: alpha, beta
    !-------------------------------------------------!
      if (y%created.and.(y%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (x%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.y%created) then
         call create_vector(y,a%dim,.true.)
      endif
      alpha = cmplx(1.0_pr,0.0_pr,pr)
      beta = cmplx(0.0_pr,0.0_pr,pr)
      call zgemv('N',a%dim,a%dim, &
         alpha , &
         a%a,a%dim, &
         x%v,1, &
         beta ,y%v,1)
   end subroutine matrix_t_vector
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Commutator C=AB-BA !
  !-------------------------------------------------!
   subroutine commutator(c,a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'commutator'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a,b
      type(matrix),intent(inout) :: c
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      complex(kind=pr) :: alpha, beta
      integer :: i,j,ii,k
    !-------------------------------------------------!
      if (c%created.and.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.c%created) then
         call create_matrix(c,a%dim,.true.)
      endif
      if (a%is_sparse) then

         do j=1,b%dim
            do ii=1,a%nonzero
               i = a%indx(ii)
               k = a%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) - b%a(j,i)*a%a(i,k)
            enddo
         enddo

      elseif (b%is_sparse) then

         do j=1,b%dim
            do ii=1,b%nonzero
               i = b%indx(ii)
               k = b%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) - b%a(j,i)*a%a(i,k)
            enddo
         enddo

      else
       ! both matrices are dense
         alpha = cmplx(1.0_pr,0.0_pr,pr)
         beta = cmplx(0.0_pr,0.0_pr,pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            b%a,b%dim, &
            a%a,a%dim, &
            beta ,c%a,c%dim)
         alpha = cmplx(1.0_pr,0.0_pr,pr)
         beta = cmplx(-1.0_pr,0.0_pr,pr)
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
   subroutine anticommutator(c,a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'anticommutator'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a,b
      type(matrix),intent(inout) :: c
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      complex(kind=pr) :: alpha, beta
      integer :: i,j,k,ii
    !-------------------------------------------------!
      if (c%created.and.(c%dim/=a%dim)) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (b%dim/=a%dim) then
         call error("Dimensions are different",myname,.true.,io)
      endif
      if (.not.c%created) then
         call create_matrix(c,a%dim,.true.)
      endif
      if (a%is_sparse) then

         do j=1,b%dim
            do ii=1,a%nonzero
               i = a%indx(ii)
               k = a%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) + b%a(j,i)*a%a(i,k)
            enddo
         enddo

      elseif (b%is_sparse) then

         do j=1,b%dim
            do ii=1,b%nonzero
               i = b%indx(ii)
               k = b%jndx(ii)
               c%a(i,j) = c%a(i,j) + a%a(i,k)*b%a(k,j)
               c%a(j,k) = c%a(j,k) + b%a(j,i)*a%a(i,k)
            enddo
         enddo

      else
       ! both matrices are dense
         alpha = cmplx(1.0_pr,0.0_pr,pr)
         beta = cmplx(0.0_pr,0.0_pr,pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            b%a,b%dim, &
            a%a,a%dim, &
            beta ,c%a,c%dim)
         alpha = cmplx(1.0_pr,0.0_pr,pr)
         beta = cmplx(1.0_pr,0.0_pr,pr)
         call zgemm('N','N',a%dim,a%dim,a%dim, &
            alpha , &
            a%a,a%dim, &
            b%a,b%dim, &
            beta ,c%a,c%dim)
      endif
   end subroutine anticommutator
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Calculate the trace of a matrix product !
  !-------------------------------------------------!
   function product_trace(a,b,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'product_trace'
    !--subroutine parameters -------------------------!
      complex(kind=pr) :: product_trace
      type(matrix),intent(in) :: a,b
      type(io_type), intent(in) :: io
    !--internal variables ----------------------------!
      integer :: i,j,ii,ji
      complex(kind=pr) :: pt
    !-------------------------------------------------!
      if ((.not.a%created).or.(.not.b%created)) then
         call error("Cannot operate on unexistant matrix",myname,.true.,io)
      endif
      product_trace= cmplx(0.0_pr,0.0_pr,pr)
      pt = cmplx(0.0_pr,0.0_pr,pr)
      if (a%is_sparse) then
         do i=1,a%nonzero
            ii = a%indx(i)
            ji = a%jndx(i)
            pt = pt + a%a(ii,ji)*b%a(ji,ii)
         enddo
      elseif (b%is_sparse) then
         do i=1,b%nonzero
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
      product_trace= pt
   end function product_trace
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix Addition C=A+B !
  !-------------------------------------------------!
   subroutine matrix_add(c,a,b,alpha,beta,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'matrix_add'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a,b
      type(matrix),intent(inout) :: c
      complex(pr), intent(in), optional :: alpha, beta
      type(io_type), intent(in) :: io
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
         call create_matrix(c,a%dim,.true.)
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
   end subroutine matrix_add
  !-------------------------------------------------!
  !-------------------------------------------------!
  ! Matrix Addition (in place) A=A+B !
  !-------------------------------------------------!
   subroutine matrix_add_inplace(a,b,alpha,beta,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'matrix_add_inplace'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: b
      type(matrix),intent(inout) :: a
      complex(pr), intent(in), optional :: alpha, beta
      type(io_type), intent(in) :: io
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
   end subroutine matrix_add_inplace
  !-------------------------------------------------!
  ! dense * sparse * dense triple product !
  !-------------------------------------------------!
   subroutine dense_sparse_dense(p,a,b,c,io)
      implicit none
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'dense_sparse_dense'
    !--subroutine parameters -------------------------!
      type(matrix),intent(in) :: a,b,c
      type(matrix),intent(inout) :: p
      type(io_type), intent(in) :: io
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
         call create_matrix(p,a%dim,.true.)
      endif

      do i=1,b%dim
         do j=1,b%dim
            p%a(i,j) = 0.0_pr
            do ii=1,b%nonzero
               k = b%indx(ii)
               l = b%jndx(ii)
               p%a(i,j) = p%a(i,j) + a%a(i,k)*b%a(k,l)*c%a(l,j)
            enddo
         enddo
      enddo

   end subroutine dense_sparse_dense
end module linear_algebra
!-------------------------------------------------!
! end of Module Linear Algebra !
!-------------------------------------------------!s