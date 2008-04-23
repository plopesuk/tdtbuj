!> \brief contains the functions and subroutines to build the density matrix
!> \author Alin M Elena
!> \date 07/11/07, 13:45:34
module m_DensityMatrix
  use m_Constants
  use m_Types
  use m_Useful
  use m_Gutenberg
  use m_LinearAlgebra, only : ZeroMatrix, aastar


  implicit none
  private

  public :: CreateDensityMatrixSpin
  public :: CreateDensityMatrixNoSpin
  public :: CreateDensityMatrixExcited
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
    integer :: i
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
    type(ioType), intent(inout) :: io
    type(solutionType), intent(inout) :: sol
    integer      :: i,k
    real(k_pr)              :: a,b
    real(k_pr) :: q

    if (gen%smearMethod /= k_smCMU) then
      select case(gen%smearMethod)
        case(k_smFD)
          a = sol%eigenvals(1)-gen%electronicTemperature*k_kb*log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)
          b = sol%eigenvals(sol%eigenvecs%dim)+ &
            gen%electronicTemperature*k_kb*log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)
         case(k_smMP,k_smCS)
            a = sol%eigenvals(1)-log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)/gen%mpW
            b = sol%eigenvals(sol%eigenvecs%dim)+log(-1.0_k_pr+1.0_k_pr/gen%qTolerance*100.0_k_pr)/gen%mpW
      end select
      do i=1,gen%maxit
          ! update mu to the centre of the interval
        gen%electronicMu = (b - a)/2.0_k_pr + a
          ! calculate total charge with current mu
        q = 0.0_k_pr

       select case(gen%smearMethod)
        case(k_smFD)
          do k=1,sol%eigenvecs%dim
            q = q + fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
          enddo
        case(k_smMP)
          do k=1,sol%eigenvecs%dim
            q = q + occupMP(gen,sol,sol%eigenvals(k))
          enddo
        case(k_smCS)
          do k=1,sol%eigenvecs%dim
            q = q + MarzariF((gen%electronicMu-sol%eigenvals(k))/gen%MPW)
          enddo
        end select
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
   end subroutine FindFermi


!> \brief builds the density matrix in the spin polarised case
!> \author Alin M Elena
!> \date 07/11/07, 16:21:21
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param sol type(solutionType) contains information about the solution space
  subroutine GenerateRho(gen,sol,io)
    character(len=*), parameter :: myname = 'GenerateRho'
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    type(solutionType), intent(inout) :: sol
    integer      :: i,k    
    real(k_pr) :: fa,fb,entropy,tiny

    ! set the occupations according to mu

    select case(gen%smearMethod)
      case(k_smFD)
        do k=1,sol%rho%dim
          sol%buff%f(k) = fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
        enddo
      case(k_smMP)
        do k=1,sol%eigenvecs%dim
          sol%buff%f(k) = occupMP(gen,sol,sol%eigenvals(k))
        enddo
      case(k_smCS)
        do k=1,sol%eigenvecs%dim
          sol%buff%f(k) =MarzariF((gen%electronicMu-sol%eigenvals(k))/gen%MPW)
        enddo
      end select

    ! The density matrix is built from the diagonal representation in the
    ! basis of the eigenstates of H (rho') by transforming with the
    ! matrix of eigenvectors that diagonalize H
    ! rho = U rho' U*
    ! this loop multiplies sqrt(rho') times U
      sol%buff%a=sol%eigenvecs%a
      do k=1,sol%eigenvecs%dim
        sol%buff%a(:,k) = sqrt(sol%buff%f(k))*sol%buff%a(:,k)
      enddo
   
    ! this does (U sqrt(rho'))(U sqrt(rho'))* , only the upper triangle is calculated
    ! ----- ---- unnocupied vectors are not multiplied not tr
      call ZeroMatrix(sol%rho,io)   
      call aastar(sol%buff%a,sol%rho%a,1.0_k_pr,0.0_k_pr,sol%rho%dim)
      entropy = 0.0_k_pr
    select case(gen%smearMethod)
      case(k_smFD)
        tiny = epsilon(1.0_k_pr)
        do i=1,sol%rho%dim
          fa = max(sol%buff%f(i),tiny)
          fb = max(1-sol%buff%f(i),tiny)
          entropy = entropy + fa*log(fa)+fb*log(fb)
        enddo
        entropy = -gen%electronicTemperature*k_kb * entropy
      case(k_smMP)
          do i=1,sol%rho%dim
            entropy=entropy+sn((sol%eigenvals(i)-gen%electronicMU)/gen%MPW,gen%mpN,sol)
       enddo
      entropy=gen%MPW*entropy
      case(k_smCS)
        do i=1,sol%rho%dim
          entropy=entropy+MarzariS((sol%eigenvals(i)-gen%electronicMU)/gen%MPW)
        enddo
      entropy=gen%MPW*entropy
    end select
      sol%electronicEntropy=-entropy
      sol%buff%pos1=0
      sol%buff%pos2=0
      call order(sol%eigenvals,sol%buff%pos1,sol%buff%pos2)
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
!> \internal activate MP smearing method
  subroutine CreateDensityMatrixNoSpin(gen,atomic,sol,io)
    character(len=*), parameter :: MyName="CreateDensityMatrixNoSpin"
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    type(solutionType), intent(inout) :: sol
      integer      :: i,k
      real(k_pr)              :: a,b      
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
    
    sol%buff%f=0.0_k_pr
    if (gen%smearMethod /= k_smCMU) then
       ! bisection to find the proper mu       !
         a = sol%eigenvals(1)
         b = sol%eigenvals(sol%rho%dim)
         do i=1,gen%maxIt
          ! update mu to the centre of the interval
            gen%electronicMu = (b - a)/2.0_k_pr + a
          ! calculate total charge with current mu
            q = 0.0_k_pr
            select case(gen%smearMethod)
              case(k_smFD)
                do k=1,sol%eigenvecs%dim
                  q = q + fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
                enddo
              case(k_smMP)
                do k=1,sol%eigenvecs%dim
                  q = q + occupMP(gen,sol,sol%eigenvals(k))
                enddo
            case(k_smCS)
                do k=1,sol%eigenvecs%dim
                  q = q + MarzariF((gen%electronicMu-sol%eigenvals(k))/gen%mpW)
                enddo
            end select
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
          select case(gen%smearMethod)
            case(k_smFD)
              sol%buff%f(k) = fermi(gen%electronicTemperature,sol%eigenvals(k),gen%electronicMu)
            case(k_smMP)
                sol%buff%f(k) =  occupMP(gen,sol,sol%eigenvals(k))
            case(k_smCS)
              sol%buff%f(k) = MarzariF((gen%electronicMu-sol%eigenvals(k))/gen%mpW)
          end select
          if (sol%buff%f(k) > gen%dmOccupationTolerance) upper_occ_state = k
          if (abs(sol%buff%f(k) - 1.0_k_pr) < gen%dmOccupationTolerance) upper_non_one = k
       enddo
!
     ! The density matrix is built from the diagonal representation in the
     ! basis of the eigenstates of H (rho') by transforming with the
     ! matrix of eigenvectors that diagonalize H
     ! rho = U rho' U*
     ! this loop multiplies sqrt(rho') times U
       do k=upper_non_one+1,upper_occ_state
          sol%eigenvecs%a(:,k) = sqrt(sol%buff%f(k))*sol%eigenvecs%a(:,k)
       enddo

     ! this does (U sqrt(rho'))(U sqrt(rho'))* , only the upper triangle is calculated
     ! unnocupied vectors are not multiplied

      call aastar(sol%eigenvecs%a(:,1:upper_occ_state),sol%rho%a,1.0_k_pr,0.0_k_pr,sol%rho%dim)
       entropy = 0.0_k_pr
        select case(gen%smearMethod)
         case(k_smFD)
           tiny = epsilon(1.0_k_pr)
           do i=upper_non_one+1,upper_occ_state
              fa = max(sol%buff%f(i),tiny)
             fb = max(1-sol%buff%f(i),tiny)
             entropy = entropy + fa*log(fa)+fb*log(fb)
          enddo
          entropy = gen%electronicTemperature*k_kb * 2.0_k_pr * entropy
       case(k_smMP)
          do i=1,sol%rho%dim
            entropy=entropy+sn((sol%eigenvals(i)-gen%electronicMU)/gen%MPW,gen%mpN,sol)
          enddo
          entropy=-2.0_k_pr*entropy*gen%MPW
        case(k_smCS)  
          do i=1,sol%rho%dim
            entropy=entropy+MarzariS((-sol%eigenvals(i)+gen%electronicMU)/gen%MPW)
          enddo
          entropy=-2.0_k_pr*entropy*gen%MPW
       end select

    ! now put it in the global variable
      sol%electronicEntropy = entropy
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
    type(atomicxType), intent(inout) :: atomic
    type(generalType),intent(inout), optional :: gen
    type(solutionType), intent(inout) :: sol
    integer :: i,n,m,j,k    

  if (present(first)) then
    if (first) then
      sol%buff%nstart=0.0_k_pr
      do i=1,atomic%atoms%natoms
        do k=1,atomic%species%norbs(atomic%atoms%sp(i))
          j= atomic%atoms%orbs(i,k)
          sol%buff%nstart(j) = atomic%basis%orbitals(j)%occup-gen%netcharge/atomic%basis%norbitals
        end do
      end do
      m=(atomic%basis%norbitals-1)*atomic%basis%norbitals/2
      sol%rho%a=cmplx(0.0_k_pr,0.0_k_pr)
      do i=1,atomic%basis%norbitals
        sol%rho%a(i,i)=cmplx(sol%buff%nstart(i),0.0_k_pr)
      end do
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
    type(solutionType), intent(inout) :: sol
    if (i==j) then
      Getn0=sol%n0(j)
    else
      Getn0=0.0_k_pr
    endif
  end function Getn0


  subroutine CreateDensityMatrixExcited(gen,atomic,sol,io)
    character(len=*), parameter :: myname = 'CreateDensityMatrixExcited'
    type(solutionType), intent(inout)    :: sol
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    type(atomicxType), intent(inout) :: atomic    
    integer      :: i,homo,extra
    real(k_pr) :: qtotal
    !-------------------------------------------------!

    call ZeroMatrix(sol%rho,io)


    ! calculate number of electrons
    qtotal = 0.0_k_pr
    do i=1,atomic%atoms%natoms
      qtotal = qtotal + atomic%species%zval(atomic%atoms%sp(i))
    enddo
    qtotal = qtotal - gen%netCharge
    if (abs(qtotal-int(qtotal))>gen%qTolerance) then
      call error("Fractional charge detected for dm spin",myname,.false.,io)
    endif

!      eigenvals(pos1(i)) would give you the i-th value from the ordered eigenvals
!     array (ascending)
!      pos2(i) would give you the position of the i-th value of original eigenvals
!     in the ordered eigenvals
    call order(sol%eigenvals,sol%buff%pos1,sol%buff%pos2)
    call FindFermi(gen,sol,io,qtotal)
    call FindHomo(gen,sol,io,sol%buff%pos1,homo)
    call CheckExcitation(gen,sol,io,homo,sol%buff%pos2,extra)
!     qtotal=qtotal+real(extra,dp)
!     call find_fermi(eigenvecs,eigenvals,qtotal)
!     print *, qtotal, general%electronic_mu
    call GenerateRhoExcited(gen,sol,io,homo,sol%buff%pos1,sol%buff%pos2)
!       if (control_var%output_level>ol_verbose) then
!          trace = matrix_trace(rho)
!          write(control_var%output_file,*) &
!             '--Density Matrix Altered---------------------------------------'
!          write(control_var%output_file,'(a,2f8.3)') &
!             'trace = ',trace
!          do i=1,rho%dim
!             write(control_var%output_file,'(a,i5,a,1500f7.3)') &
!                'row ',i,'=',(real(rho%a(i,j)),j=1,rho%dim)
!          enddo
!          write(control_var%output_file,*)' '
!          do i=1,rho%dim
!             write(control_var%output_file,'(a,i5,a,1500f7.3)') &
!                'row ',i,'=',(aimag(rho%a(i,j)),j=1,rho%dim)
!          enddo
!          write(control_var%output_file,*) &
!             '---------------------------------------------------------------'
!       endif    
  end subroutine CreateDensityMatrixExcited

  subroutine FindHomo(gen,sol,io,pos,homoLevel)
    character(len=*), parameter :: myname = 'FindHomo'
    integer, intent(in) :: pos(:)
    integer,intent(out) :: homoLevel
    type(solutionType), intent(inout) :: sol
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    !--internal variables ----------------------------!
    integer      :: k    
    !-------------------------------------------------!

        ! set the occupations according to mu
    homoLevel=1

    select case(gen%smearMethod)
      case(k_smFD)
        do k=1,sol%rho%dim
          sol%buff%f(pos(k)) = fermi(gen%electronicTemperature,sol%eigenvals(pos(k)),gen%electronicMu)
        enddo
      case(k_smMP)
        do k=1,sol%eigenvecs%dim
          sol%buff%f(pos(k)) = occupMP(gen,sol,sol%eigenvals(pos(k)))
        enddo
        case(k_smCS)
          do k=1,sol%eigenvecs%dim
          sol%buff%f(pos(k)) = MarzariS((-sol%eigenvals(pos(k))+gen%electronicMu)/gen%mpW)
        enddo
    end select

    do k=1,sol%rho%dim
      if (sol%buff%f(pos(k))<gen%dmOccupationTolerance) then
        homoLevel=k-1
        exit
      endif
    enddo

    if (io%Verbosity>=k_HighVerbos) then
      write(io%uout,"(a,i0,a,i0)") &
          "HOMO is at: ",homoLevel, " in the output ", pos(homoLevel)
    endif
  end subroutine FindHomo


  subroutine CheckExcitation(gen,sol,io,homo,pos,extra)
    character(len=*), parameter :: myname = 'CheckExcitation'
    integer, intent(in) :: homo
    integer, intent(in) :: pos(:)
    integer,intent(out) :: extra
    type(solutionType), intent(inout) :: sol
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    integer :: hole,excite,n


    n=sol%rho%dim
    hole=gen%holeState+gen%holeSpin*n/2
    excite=gen%exciteState+gen%exciteSpin*n/2

    if (pos(hole)>homo) then
      call error("You can not create a hole at a level already empty",myname,.true.,io)
    endif

    if (pos(excite)<homo) then
      call error("You can not create an excitation at a level already occupied",myname,.true.,io)
    endif
    extra=pos(excite)-homo+1
  end subroutine CheckExcitation

  subroutine GenerateRhoExcited(gen,sol,io,homo,pos1,pos2)
    character(len=*), parameter :: myname = 'GenerateRhoExcited'
    type(solutionType), intent(inout) :: sol
    type(generalType), intent(inout) :: gen
    type(ioType), intent(inout) :: io
    integer, intent(in) :: pos1(:),pos2(:),homo
    integer      :: i,k,hole,excite,n    
    real(k_pr) :: fa,fb,entropy,tiny
    
    n=sol%rho%dim    
    ! set the occupations according to mu
    

   select case(gen%smearMethod)
    case(k_smFD)
        do k=1,n
          sol%buff%f(k) = fermi(gen%electronicTemperature,sol%eigenvals(pos1(k)),gen%electronicMu)
        enddo
      case(k_smMP)
        do k=1,sol%eigenvecs%dim
          sol%buff%f(k) =occupMP(gen,sol,sol%eigenvals(pos1(k)))
        enddo
       case(k_smCS)
        do k=1,sol%eigenvecs%dim
          sol%buff%f(k) =MarzariF((-sol%eigenvals(pos1(k))+gen%electronicMu)/gen%mpW)
        enddo         
    end select

    hole=gen%holeState + gen%holeSpin * n/2
    excite=gen%exciteState + gen%exciteSpin * n/2
    sol%buff%g=sol%buff%f
!     print *, pos2(excite),f(pos2(excite)),pos2(hole),f(pos2(hole))
    tiny=sol%buff%f(pos2(excite))
    sol%buff%f(pos2(excite))=sol%buff%f(pos2(hole))
    sol%buff%f(pos2(hole))=tiny
!     print *, pos2(excite),f(pos2(excite)),pos2(hole),f(pos2(hole))
!     do k=homo+1,pos2(excite)-1
!     f(k)=1.0_k_pr-f(k)
!     enddo
!     f(pos2(excite)+1)=1.0_k_pr-f(pos2(excite)+1)

    ! The density matrix is built from the diagonal representation in the
    ! basis of the eigenstates of H (rho') by transforming with the
    ! matrix of eigenvectors that diagonalize H
    ! rho = U rho' U*
    ! this loop multiplies sqrt(rho') times U
    sol%buff%a=sol%eigenvecs%a
    do k=1,n
      sol%buff%a(:,pos1(k)) = sqrt(sol%buff%f(k))*sol%buff%a(:,pos1(k))
    enddo

    call ZeroMatrix(sol%rho,io)
    ! this does (U sqrt(rho'))(U sqrt(rho'))* , only the upper triangle is calculated
    ! ----- ---- unnocupied vectors are not multiplied not tr
    call aastar(sol%buff%a,sol%rho%a,1.0_k_pr,0.0_k_pr,n)
    entropy = 0.0_k_pr
    select case(gen%smearMethod)
      case(k_smFD)
        tiny = epsilon(1.0_k_pr)
        do i=1,sol%rho%dim
          fa = max(sol%buff%f(i),tiny)
          fb = max(1-sol%buff%f(i),tiny)
          entropy = entropy + fa*log(fa)+fb*log(fb)
        enddo
        entropy = -gen%electronicTemperature*k_kb * entropy
      case(k_smMP)
        do i=1,sol%rho%dim
         entropy=entropy+sn((sol%eigenvals(i)-gen%electronicMU)/gen%MPW,gen%mpN,sol)
       enddo
       entropy=gen%MPW*entropy
      case(k_smCS)
       do i=1,sol%rho%dim
         entropy=entropy+MarzariS((-sol%eigenvals(i)+gen%electronicMU)/gen%MPW)
       enddo
       entropy=gen%MPW*entropy
    end select
    sol%electronicEntropy=-entropy

    if (io%verbosity>=k_HighVerbos) then
      write(io%uout,*) &
        '--Occupation Numbers Altered-----------------------------------'
      do i=1,n
        write(io%uout,'(3f16.8)') &
             sol%eigenvals(pos1(i)),sol%buff%f(i),sol%buff%g(i)
      enddo
      write(io%uout,*) &
            '------------------------------------------------------------'
    endif
   end subroutine GenerateRhoExcited
end module m_DensityMatrix
