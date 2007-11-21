!> \brief contains the functions and subroutines to build the density matrix
!> \author Alin M Elena
!> \date 07/11/07, 13:45:34
module m_DensityMatrix
  use m_Constants
  use m_Types
  use m_Useful
  use m_Gutenberg
  use m_LinearAlgebra, only : ZeroMatrix

  implicit none
  private

  public :: CreateDensityMatrixSpin
  public :: CreateDensityMatrixNoSpin
  public :: BuildDensity

contains

!> \brief cerates density matrix for spin polarised case
!> \author Alin M Elena
!> \date 07/11/07, 13:49:34
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param sol type(solutionType) contains information about the solution space
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
  subroutine CreateDensityMatrixSpin(gen,atomic,sol,io)
    character(len=*), parameter :: MyName="CreateDensityMatrixSpin"
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    type(solutionType), intent(inout) :: sol 
    integer :: i,j
    complex(k_pr) :: trace
    real(k_pr) :: qtotal
    character(len=k_ml) :: saux

    call ZeroMatrix(sol%rho,io)
    qtotal = 0.0_k_pr
    do i=1,atomic%atoms%natoms
      qtotal = qtotal + atomic%species%zval(atomic%atoms%sp(i))
    enddo
    qtotal = qtotal - gen%netcharge
    if (abs(qtotal-int(qtotal))>gen%qTolerance) then
      write(saux,'(a,f0.8)')"Fractional charge detected found ",qtotal
      call error(trim(saux),myname,.true.,io)
    endif

    call FindFermi(gen,sol,io,qtotal)
    call GenerateRho(gen,sol,io)
      
  end subroutine CreateDensityMatrixSpin

!> \brief finds the chemical potential
!> \author Alin M Elena
!> \date 07/11/07, 15:55:22
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param sol type(solutionType) contains information about the solution space
!> \param qtotal real the total charge
  subroutine FindFermi(gen,sol,io,qtotal)
    character(len=*), parameter :: myname = 'FindFermi'
    real(k_pr), intent(in) :: qtotal
    type(generalType), intent(inout) :: gen
    type(ioType), intent(in) :: io
    type(solutionType), intent(in) :: sol
    integer      :: i,k
    real(k_pr)              :: a,b
    real(k_pr) :: q

    if (gen%smearMethod /= k_smCMU) then
      select case(gen%smearMethod)
        case(k_smFD)
          a = sol%eigenvals(1)-gen%electronicTemperature*k_kb*log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)
          b = sol%eigenvals(sol%eigenvecs%dim)+ &
            gen%electronicTemperature*k_kb*log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)
         case(k_smMP)
            a = sol%eigenvals(1)-log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)/gen%mpW
            b = sol%eigenvals(sol%eigenvecs%dim)+log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)/gen%mpW
      end select
      do i=1,gen%maxit
          ! update mu to the centre of the interval
        gen%electronicMu = (b - a)/2.0_k_pr + a
          ! calculate total charge with current mu
        q = 0.0_k_pr
        do k=1,sol%eigenvecs%dim
          q = q + fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
        enddo

        if ((q-qtotal)>0.0_k_pr) then
             ! we have more electrons than we need
             ! set mu to be the right of the interval
          b = gen%electronicMu
        else
             ! we have less electrons
             ! set mu to be the left of the interval
          a = gen%electronicMu
        endif
          if (abs(q-qtotal)<gen%qTolerance) exit
      enddo
      if (i==(gen%maxit+1)) call error("Could not find the fermi level",myname,.true.,io)
    endif

    if (io%verbosity>=k_highVerbos) then
      write(io%uout,"(a,f16.8)") &
            "chemical potential: ",gen%electronicMu
    endif

   end subroutine FindFermi


!> \brief builds the density matrix in the spin polarised case
!> \author Alin M Elena
!> \date 07/11/07, 16:21:21
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param sol type(solutionType) contains information about the solution space
!> \internal zherk should be called in Linearalgebra, implement MP
  subroutine GenerateRho(gen,sol,io)
    character(len=*), parameter :: myname = 'GenerateRho'
    type(generalType), intent(in) :: gen
    type(ioType), intent(in) :: io
    type(solutionType), intent(inout) :: sol
    integer      :: i,j,k
    real(k_pr),allocatable  :: f(:)
    real(k_pr) :: fa,fb,entropy,tiny
    complex(k_pr),allocatable :: a(:,:)
    integer, allocatable :: pos1(:), pos2(:)

      allocate(f(1:sol%eigenvecs%dim))
      allocate(a(1:sol%eigenvecs%dim,1:sol%eigenvecs%dim))


    ! set the occupations according to mu
      do k=1,sol%rho%dim
        f(k) = fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
      enddo

    ! The density matrix is built from the diagonal representation in the
    ! basis of the eigenstates of H (rho') by transforming with the
    ! matrix of eigenvectors that diagonalize H
    ! rho = U rho' U*
    ! this loop multiplies sqrt(rho') times U
      a=sol%eigenvecs%a
      do k=1,sol%eigenvecs%dim
        a(:,k) = sqrt(f(k))*a(:,k)
      enddo

    ! this does (U sqrt(rho'))(U sqrt(rho'))* , only the upper triangle is calculated
    ! ----- ---- unnocupied vectors are not multiplied not tr
      call zherk ( 'U', 'N', sol%rho%dim, sol%rho%dim, 1.0_k_pr, a, sol%rho%dim, 0.0_k_pr, sol%rho%a, sol%rho%dim )

    ! this fills in the lower triangle
      do i=1,sol%rho%dim
        sol%rho%a(i+1:sol%rho%dim,i) = conjg(sol%rho%a(i,i+1:sol%rho%dim))
      enddo

      entropy = 0.0_k_pr
      select case(gen%smearMethod)
      case(k_smFD)
        tiny = epsilon(1.0_k_pr)
        do i=1,sol%rho%dim
          fa = max(f(i),tiny)
          fb = max(1-f(i),tiny)
          entropy = entropy + fa*log(fa)+fb*log(fb)
        enddo
        entropy = k_kb * entropy
      case(k_smMP)
!          do i=1,rho%dim
!             entropy=entropy+sn((eigenvals(i)-general%electronic_mu)/general%electronic_temperature,general%mp_N)
!          enddo
      end select

      sol%electronicEntropy=-entropy
      allocate(pos1(1:sol%rho%dim),pos2(1:sol%rho%dim))
      call order(sol%eigenvals,pos1,pos2)
      if (io%Verbosity>=k_highVerbos) then
        write(io%uout,*) &
            '--Occupation Numbers-------------------------------------------'
        do i=1,sol%rho%dim
          if (pos1(i) <=sol%rho%dim/2) then
            write(io%uout,'(2f16.8,a)') &
              sol%eigenvals(pos1(i)),f(pos1(i))," down"
          else
            write(io%uout,'(2f16.8,a)') &
              sol%eigenvals(pos1(i)),f(pos1(i))," up"
          endif
        enddo
        write(io%uout,*) &
            '------------------------------------------------------------'
      endif
    deallocate (pos1,pos2)
    deallocate(f,a)

   end subroutine GenerateRho

!> \brief sorts the eigenvalues recording the old positions
!> \author Alin M Elena
!> \date 07/11/07, 16:32:01
!> \param evals array eigenvalues
!> \param tp, spin arrays keep the indeces 
   subroutine order(evals,tp,spin)
    !--subroutine name--------------------------------!   
      character(len=*), parameter :: myname = 'order'
      real(k_pr),intent(inout) :: evals(:)
      integer,intent(out) :: spin(:),tp(:)

    !--internal variables ----------------------------!
      integer :: i,j,k,n
      n=size(evals)
      do i=1,n
         tp(i)=i
      end do
      spin=0
!sorting
      do i=1,n-1
         do j=i+1,n
            if (evals(tp(i))>evals(tp(j))) then
          ! then swap them
               k=tp(i)
               tp(i)=tp(j)
               tp(j)=k
            endif
         enddo
      enddo

      do i=1,n
         k=tp(i)
         spin(k)=i
      enddo
   end subroutine order

!> \brief creates density matrix for a non spin polarised system
!> \author Alin M Elena
!> \date 07/11/07, 13:51:34
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param sol type(solutionType) contains information about the solution space
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \internal the call too zherk should be done via LinearAlgebra
!> \internal activate MP smearing method
  subroutine CreateDensityMatrixNoSpin(gen,atomic,sol,io)
    character(len=*), parameter :: MyName="CreateDensityMatrixNoSpin"
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    type(solutionType), intent(inout) :: sol
      integer      :: i,j,k
      real(k_pr)              :: a,b
      complex(k_pr)           :: trace
      real(k_pr),allocatable  :: f(:)
      real(k_pr) :: qtotal,q
      real(k_pr) :: entropy,tiny,fa,fb
      integer :: upper_occ_state, upper_non_one
    !-------------------------------------------------!
    call ZeroMatrix(sol%rho,io)
   ! calculate number of electrons
    qtotal = 0.0_k_pr
    do i=1,atomic%atoms%natoms
      qtotal = qtotal + atomic%species%zval(atomic%atoms%sp(i))
    enddo
    qtotal = qtotal - gen%netcharge
    qtotal = 0.5_k_pr * qtotal
    allocate(f(sol%rho%dim))
    f=0.0_k_pr
    if (gen%smearMethod /= k_smCMU) then
       ! bisection to find the proper mu       !
         a = sol%eigenvals(1)
         b = sol%eigenvals(sol%rho%dim)
         do i=1,gen%maxIt
          ! update mu to the centre of the interval
            gen%electronicMu = (b - a)/2.0_k_pr + a
          ! calculate total charge with current mu
            q = 0.0_k_pr
            do k=1,sol%rho%dim
               q = q + fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
            enddo

            if ((q-qtotal)>0.0_k_pr) then
             ! we have more electrons than we need
             ! set mu to be the right of the interval
               b = gen%electronicMu
            else
             ! we have less electrons
             ! set mu to be the left of the interval
               a = gen%electronicMu
            endif
          !            write(*,'(i4,4f13.6)')i,q-qtotal,gen%electronicMu,a,b
          ! if we are happy with the charge then exit
            if (abs(q-qtotal)<gen%qTolerance) exit
         enddo
         if (i==(gen%maxIt+1)) call error("Could not find the fermi level",myname,.true.,io)
      endif

     ! set the occupations according to mu
       upper_occ_state = 1
       upper_non_one = 0
       do k=1,sol%rho%dim
          f(k) = fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
          if (f(k) > gen%dmOccupationTolerance) upper_occ_state = k
          if (abs(f(k) - 1.0_k_pr) < gen%dmOccupationTolerance) upper_non_one = k
       enddo
! 
     ! The density matrix is built from the diagonal representation in the
     ! basis of the eigenstates of H (rho') by transforming with the
     ! matrix of eigenvectors that diagonalize H
     ! rho = U rho' U*
     ! this loop multiplies sqrt(rho') times U
       do k=upper_non_one+1,upper_occ_state
          sol%eigenvecs%a(:,k) = sqrt(f(k))*sol%eigenvecs%a(:,k)
       enddo
 
     ! this does (U sqrt(rho'))(U sqrt(rho'))* , only the upper triangle is calculated
     ! unnocupied vectors are not multiplied
       sol%rho%a=0.0_k_pr

       call zherk ( 'U', 'N', sol%rho%dim, upper_occ_state, 1.0_k_pr, sol%eigenvecs%a, &
                    sol%rho%dim, 0.0_k_pr, sol%rho%a, sol%rho%dim )
 
     ! this fills in the lower triangle
       do i=1,sol%rho%dim
          sol%rho%a(i+1:sol%rho%dim,i) = conjg(sol%rho%a(i,i+1:sol%rho%dim))
       enddo
! 
!       call write_dos(eigenvals,eigenvecs,gen%electronicMu)
! 
     ! calculate the contribution to the free energy
       entropy = 0.0_k_pr
        select case(gen%smearMethod)
         case(k_smFD)
           tiny = epsilon(1.0_k_pr)
           do i=upper_non_one+1,upper_occ_state
              fa = max(f(i),tiny)
             fb = max(1-f(i),tiny)
             entropy = entropy + fa*log(fa)+fb*log(fb)
          enddo
          entropy = -k_kb * 2.0_k_pr * entropy
       case(k_smMP)
!          do i=1,sol%rho%dim
!             entropy=entropy+sn((eigenvals(i)-gen%electronicMu)/gen%electronic_temperature,gen%mp_N)
!          enddo
!          entropy=-2.0_k_pr*entropy
       end select

    ! now put it in the global variable
      sol%electronicEntropy = entropy

      if (io%verbosity>=k_highVerbos) then
        write(io%uout,*) &
            "Ocupation Numbers"
         do i=1,sol%rho%dim
            write(io%uout,'(i0,1x,f16.8,1x,f16.8)') &
               i,sol%eigenvals(i),2.0_k_pr*f(i)
         enddo
      endif
      deallocate(f)

  end subroutine CreateDensityMatrixNoSpin

!> \brief computes the charge excess \f$ \delta q\f$ array
!> \author Alin M Elena
!> \date 08/11/07, 12:04:33
!> \param sol type(solutionType) contains information about the solution space
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param first logical if present and true builds the initial guess for the density matrix
!> \param gen type(generalType) contains the info needed by the program to k_run
  subroutine BuildDensity(atomic,sol,gen,first)
    character(len=*), parameter :: myname = 'BuildDensity'
    logical, intent(in),optional :: first
    type(atomicxType), intent(in) :: atomic
    type(generalType),intent(in), optional :: gen
    type(solutionType), intent(inout) :: sol
    integer :: i,n,m,j,k
    real(k_pr), allocatable :: nstart(:)

  if (present(first)) then
    if (first) then
      allocate(nstart(1:atomic%basis%norbitals))
      do i=1,atomic%atoms%natoms
        do k=1,atomic%species%norbs(atomic%atoms%sp(i))
          j= atomic%atoms%orbs(i,k)
          nstart(j) = atomic%basis%orbitals(j)%occup-gen%netcharge/atomic%basis%norbitals
        end do
      end do
      m=(atomic%basis%norbitals-1)*atomic%basis%norbitals/2
      sol%rho%a=cmplx(0.0_k_pr,0.0_k_pr)
      do i=1,atomic%basis%norbitals
        sol%rho%a(i,i)=cmplx(nstart(i),0.0_k_pr)
      end do
    deallocate(nstart)
    endif
  else
! store the off diagonal terms
! the number of off diagonal terms
    m=(atomic%basis%norbitals-1)*atomic%basis%norbitals/2
    do i=1,atomic%basis%norbitals-1
      n=ka(i,atomic%basis%norbitals)+1
! the off diagonal
      sol%density(n:n-1+atomic%basis%norbitals-i)=sol%rho%a(i,i+1:atomic%basis%norbitals)
! the diagonal
      sol%density(m+i)=sol%rho%a(i,i)-getn0(i,i,sol)
    end do
! the last diagonal term
    sol%density(m+atomic%basis%norbitals)=sol%rho%a(atomic%basis%norbitals,atomic%basis%norbitals)&
      -Getn0(atomic%basis%norbitals,atomic%basis%norbitals,sol)
    endif
  end subroutine BuildDensity

!> \brief extracts the i,j element from the reference density matrix
!> \author Alin M Elena
!> \date 08/11/07, 12:14:09
!> \param i,j integers coordinates of the matrix element
!> \param sol type(solutionType) contains information about the solution space

  real(k_pr) function Getn0(i,j,sol)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'Getn0'
    !--subroutine parameters -------------------------!
    integer, intent(in) :: i,j
    type(solutionType), intent(in) :: sol
    if (i==j) then
      Getn0=sol%n0(j)
    else
      Getn0=0.0_k_pr
    endif
  end function Getn0
end module m_DensityMatrix