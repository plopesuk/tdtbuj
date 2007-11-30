!> \brief deals with the general tight binding functions
!> \author Alin M Elena
!> \date 01/11/07, 11:48:23

module m_TightBinding
  use m_Constants
  use m_Types
  use m_Useful
  use m_TailFunctions
  use m_Gutenberg
  use m_LinearAlgebra
  implicit none

  private

  public :: SetSolutionSpace
  public :: rad
  public :: RadP
  public :: RadPp
  public :: rep
  public :: RepP
  public :: RepPp
  public :: onsite
  public :: embedding
  public :: EmbeddingP
  public :: EmbeddingPp
  public :: InitMagneticMoment
  public :: ComputeMagneticMoment
  public :: SetTails

contains
!> \brief Initializes the solution space
!> \details allocates memory for all the variables necessary for the calculation. It will
!> also call the routines that precompute some values (eg factorial, Gaunt coefficients)
!> \author Alin M Elena
!> \date 01/11/07, 11:49:38
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space

  subroutine SetSolutionSpace(ioLoc,genLoc,atomic,tbMod,sol)
    character(len=*), parameter :: sMyName="SetSolutionSpace"
    type(ioType), intent(inout) :: ioLoc
    type(generalType), intent(inout) :: genLoc
    type(atomicxType), intent(inout) :: atomic
    type(modelType), intent(inout) :: tbMod
    type(solutionType), intent(inout) :: sol
    integer :: l
    real(k_pr) :: cr

    call cpu_time(cr)
    if (.not.sol%h%created) then
      call CreateSparseMatrix(sol%h,atomic%basis%norbitals,.true.)
      call CreateSparseMatrix(sol%forceOp,sol%h%dim,.true.)
    else
      call ResetSparseMatrix(sol%h)
      call ResetSparseMatrix(sol%forceOp)
    endif
    if (.not.sol%eigenvecs%created) then
      allocate(sol%eigenvals(1:atomic%basis%norbitals))
      call CreateMatrix(sol%eigenvecs,atomic%basis%norbitals,.true.)
    else
      sol%eigenvals=0.0_k_pr
      call ZeroMatrix(sol%eigenvecs,ioLoc)
    endif
    if (genLoc%scf) then
      if (sol%hin%created) then
        call ResetSparseMatrix(sol%hin)
        call ResetSparseMatrix(sol%h2)
        sol%potential=0.0_k_pr
        sol%field=0.0_k_pr
      else
        call CreateSparseMatrix(sol%hin,atomic%basis%norbitals,.true.)
        call CreateSparseMatrix(sol%h2,atomic%basis%norbitals,.true.)
        allocate(sol%potential(1:atomic%atoms%natoms))
        allocate(sol%field(1:atomic%atoms%natoms,1:3))
        call ResetSparseMatrix(sol%hin)
        call ResetSparseMatrix(sol%h2)
        sol%potential=0.0_k_pr
        sol%field=0.0_k_pr
      endif

    endif

    if (genLoc%spin) then
       !spin down
      if (.not.sol%hdown%created) then
        call CreateMatrix(sol%hdown,atomic%basis%norbitals/2,.true.)
      else
        call ZeroMatrix(sol%hdown,ioLoc)
      endif
      if (.not.sol%hup%created) then
        call CreateMatrix(sol%hup,atomic%basis%norbitals/2,.true.)
      else
        call ZeroMatrix(sol%hup,ioLoc)
      endif
    endif
    if (.not.sol%rho%created) then
      call CreateMatrix(sol%rho,atomic%basis%norbitals,.true.)
    else
      call ZeroMatrix(sol%rho,ioLoc)
    endif

    l= maxval(tbMod%hopping(:,:)%l2)
    allocate(sol%fact(0:(6*(l+1)+1)))
    call initFact(6*(l+1)+1,sol%fact)

    if ((genLoc%scf).and.(genLoc%electrostatics==k_electrostaticsMultipoles)) then
      call InitGaunt(2*l+1,2*l+1,4*l+2,sol,ioLoc)
    endif
    if (.not.allocated(sol%n0)) then
      allocate(sol%n0(1:atomic%basis%norbitals))
    endif
    sol%n0=0.0_k_pr
    call Buildn0(genLoc,atomic,sol)
    l=sol%h%dim*(sol%h%dim-1)/2
    allocate(sol%density(1:l+sol%h%dim))
    call setTails(ioLoc,genLoc,atomic,tbMod,sol)
    call rmarin(int((cr-int(cr))*3132),int(genLoc%ranseed*30081),sol%seed,ioLoc)

  end subroutine SetSolutionSpace


!> \brief sets the tails parameters
!> \author Alin M Elena
!> \date 02/11/07, 16:18:21
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space

  subroutine setTails(ioLoc,genLoc,atomic,tbMod,sol)
    character(len=*), parameter :: sMyName="setTails"
    type(ioType), intent(inout) :: ioLoc
    type(generalType), intent(inout) :: genLoc
    type(atomicxType), intent(inout) :: atomic
    type(modelType), intent(inout) :: tbMod
    type(solutionType), intent(inout) :: sol

    integer ::i,j
    real(k_pr) :: f,fp,fpp
    integer :: k,k1,k2
!     real(k_pr) :: r
!     integer :: z,y
!    calculate the tail function parameters
    do i=1,atomic%species%nspecies
      do j=1,atomic%species%nspecies
          ! repulsions
        f = RepNoTail(tbMod%hopping(i,j)%d1,atomic%species%id(i),atomic%species%id(j),tbMod,genLoc)
        fp = RepPNoTail(tbMod%hopping(i,j)%d1,atomic%species%id(i),atomic%species%id(j),tbMod,genLoc)
        fpp = RepPpNoTail(tbMod%hopping(i,j)%d1,atomic%species%id(i),atomic%species%id(j),tbMod,genLoc)
        tbMod%hopping(i,j)%repTail=makeTailType(f,fp,fpp,tbMod%hopping(i,j)%d1,tbMod%hopping(i,j)%dcut,ioLoc)

          ! hoppings
        if (.not.allocated(tbMod%hopping(i,j)%hmnTail)) then
          allocate(tbMod%hopping(i,j)%hmnTail(0:tbMod%hopping(i,j)%l1,0:tbMod%hopping(i,j)%l2,0:tbMod%hopping(i,j)%ll))
        endif
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
!               write(io%uout,*) &
!               trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2)
              f = tbMod%hopping(i,j)%a(k,k1,k2)*RadNoTail(tbMod%hopping(i,j)%r1,&
                    atomic%species%id(i),atomic%species%id(j),genLoc,tbMod)
              fp = tbMod%hopping(i,j)%a(k,k1,k2)*RadPNoTail(3,tbMod%hopping(i,j)%r1,&
                    atomic%species%id(i),atomic%species%id(j),genLoc,tbMod)
              fpp = tbMod%hopping(i,j)%a(k,k1,k2)*RadPpNoTail(3,3,tbMod%hopping(i,j)%r1,&
                    atomic%species%id(i),atomic%species%id(j),genLoc,tbMod)
              tbMod%hopping(i,j)%hmnTail(k,k1,k2) = &
              makeTailType(f,fp,fpp,tbMod%hopping(i,j)%r1,tbMod%hopping(i,j)%rcut,ioLoc)
            enddo
          enddo
        enddo

      enddo
    enddo

    call PrintTail(atomic,tbMod,ioLoc)
!     z=600
!     do i=1,atomic%species%nspecies
!       do j=1,atomic%species%nspecies
!         do k=0,tbMod%hopping(i,j)%l1
!           do k1=0,tbMod%hopping(i,j)%l2
!             do k2=0,min(k,k1)
!               write(z,*)"# r ",trim(ccnlm(i,j,k,k1,k2))
!               do y=50,1000
!                 r=y*0.01_k_pr
!                 write(z,*)r,rad(r,atomic%species%id(i),atomic%species%id(j),genLoc,tbMod,k,k1,k2)
!               enddo
!               z=z+1
!             enddo
!           enddo
!         enddo
!       enddo
!     enddo
end subroutine setTails





!==Gaunt coefficients====================================================
!> \brief returns the Gaunt coefficient for l1,m1,l2,m2,l3,m3
!> \author Alin M Elena
!> \date 01/11/07, 17:15:07
!> \param l1,m1,l2,m2,l3,m3 angular momentum and quantum magnetic numbers
!> \param io type(ioType) contains all the info about I/O files
!> \param sol type(solutionType) contains information about the solution space
  real(k_pr) function cgaunt(l1,m1,l2,m2,l3,m3,sol,io)
    character(len=*),parameter :: myname="cgaunt"
    integer, intent(in) :: l1,m1,l2,m2,l3,m3
    type(ioType), intent(in) :: io
    type(solutionType), intent(in) :: sol
    real(k_pr) :: sq
    integer ::g
    character(len=k_ml) :: saux
     !     check for validity
    if (l1<0 .or. l2<0 .or. l3<0 ) then
      write(saux,'(a,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') "error in the input l-values in gaunt", &
            l1,m1, l2,m2,l3,m3
      call error(trim(saux),myname,.true.,io)
    end if
    if (abs(m1)>l1 .or. abs(m2)>l2  &
      .or. abs(m3)>l3) then
      write(saux,'(a,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') "error in the input m-values in gaunt", &
        l1,m1, l2,m2,l3,m3
      call error(trim(saux),myname,.true.,io)
    end if

    if ((abs(l1-l2)<=l3).and.(l3<=l1+l2).and.&
      (mod(l1+l2+l3,2)==0).and. (m1==m3-m2))then

      g=(l1+l2+l3)/2
      sq=real((2*l1+1)*(2*l2+1),k_pr)
      sq=sqrt(sq/(4.0_k_pr*k_pi))*sol%fact(g)/(sol%fact(g-l1)*sol%fact(g-l2)*sol%fact(g-l3))
      cgaunt=(-1)**(g-l3)*sq*cgc(l1,m1,l2,m2,l3,m3,sol)*dl(l1,l2,l3,sol)
    else
      cgaunt=0.0_k_pr
    endif

  end function cgaunt

!> \brief the "heart" of cgaunt
!> \author Alin M Elena
!> \date 01/11/07, 17:20:53
!> \param a,x,b,y,c,z
!> \param sol type(solutionType) contains information about the solution space
  real(k_pr) function cgc(a,x,b,y,c,z,sol)
    integer, intent(in) :: a,x,b,y,c,z
    type(solutionType), intent(in) :: sol
    integer :: kmin,kmax,k
    real(k_pr) :: sum,gamma

    gamma=sqrt(sol%fact(a+x)*sol%fact(b+y)*sol%fact(c+z)* &
           sol%fact(a-x)*sol%fact(b-y)*sol%fact(c-z)*real(2*c+1,k_pr))
    sum=0.0_k_pr
    kmin=max(b-c-x,a+y-c,0)
    kmax=min(b+y,a-x,a+b-c)
    do k=kmin,kmax
    sum=sum+(-1)**k/(sol%fact(k)*sol%fact(a+b-c-k)*sol%fact(a-x-k)* &
        sol%fact(b+y-k)*sol%fact(c-b+x+k)*sol%fact(c-a-y+k))
    enddo
    cgc=dl(a,b,c,sol)*gamma*sum
  end function cgc

  real(k_pr) function dl(l1,l2,l3,sol)
    integer, intent(in):: l1,l2,l3
    type(solutionType), intent(in) :: sol
    integer :: l
    l=l1+l2+l3
    dl=sqrt(sol%fact(l-2*l1)*sol%fact(l-2*l2)*sol%fact(l-2*l3)/sol%fact(l+1))
  end function dl

  integer function idx(l,m)
    integer :: l,m
    idx=l*l+l+1+m
  end function idx

  subroutine  ridx(k,l,m)
    integer :: l,m, k
    l=int(sqrt(real(k-1)))
    m=k-1-l*l-l
  end subroutine ridx

  complex(k_pr) function ulmmiu(m,miu)
    integer :: m,miu,z
    real(k_pr) :: sq2,r,i

    sq2=dsqrt(2.0_k_pr)
    z=0
    r=d(m,z)*d(miu,z)+h(miu)*(d(m,miu)+(-1)**m*d(m,-miu))/sq2
    i=h(-miu)*(-d(m,-miu)+d(m,miu)*(-1)**m)/sq2

    ulmmiu=cmplx(r,i,k_pr)
  end function ulmmiu

!> \brief Builds the Gaunt coefficients tables
!> \author Alin M Elena
!> \param l1m,l2m,l3m integers maxima ls for angular momenta
!> \param io type(ioType) contains all the info about I/O files
!> \param sol type(solutionType) contains information about the solution space
  subroutine InitGaunt(l1m,l2m,l3m,sol,io)
    character(len=*),parameter :: myname="InitGaunt"
    integer,intent(in) :: l1m,l2m,l3m
    type(ioType), intent(in) :: io
    type(solutionType), intent(inout) :: sol
    integer :: n1,n2,n3
    integer :: i,j,k,l1,l2,l3,m1,m2,m3,ii,kk,jj,p=0
    complex(k_pr) :: aux,tt

    n1=l1m*l1m+2*l1m+1
    n2=l2m*l2m+2*l2m+1
    n3=l3m*l3m+2*l3m+1
    allocate(sol%gcoeff(1:n1,1:n2,1:n3),sol%rgc(1:n1,1:n2,1:n3))

    sol%gcoeff=0.0_k_pr
    sol%rgc=0.0_k_pr

    do l1=0,l1m
      do m1=-l1,l1
        do l2=0,l2m
          do m2=-l2,l2
            do l3=0,l3m
              do m3=-l3,l3
                sol%gcoeff(idx(l1,m1),idx(l2,m2),idx(l3,m3))=&
                  cgaunt(l1,m1,l2,m2,l3,m3,sol,io)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    if (io%verbosity >= k_highVerbos) then
      write(io%uout,'(/a)')&
        "===Gaunt Coefficients and Real Gaunt coefficients=================="
      write(io%uout,*)"  Gaunt R-Gaunt  l1 m1 l2 m2 l3 m3 i  j  k"
    endif
    do i=1,n1
      call ridx(i,l1,m1)
      do j=1,n2
        call ridx(j,l2,m2)
        do k=1,n3
          call ridx(k,l3,m3)
          tt=cmplx(0.0_k_pr,0.0_k_pr,k_pr)

          do ii=-l1,l1
            do jj=-l2,l2
              do kk=-l3,l3
                if (sol%gcoeff(idx(l1,ii),idx(l2,jj),idx(l3,kk))/=k_zero) then
                  aux=ulmmiu(ii,m1)
                  tt=tt+(-1)**kk*ulmmiu(ii,m1)*ulmmiu(jj,m2)*ulmmiu(-kk,m3)*sol%gcoeff(idx(l1,ii),idx(l2,jj),idx(l3,kk))
                endif
              enddo
            enddo
          enddo
          sol%rgc(i,j,k)=real(tt,k_pr)
          if (io%verbosity >= k_highVerbos) then
            if ((sol%gcoeff(i,j,k)/=k_zero).or.(sol%rgc(i,j,k)/=k_zero)) then
              write(io%uout,'(2F8.4,6I3,3I3)')sol%gcoeff(i,j,k),sol%rgc(i,j,k),l1,m1,l2,m2,l3,m3,i,j,k
              p=p+1
              if (mod(p,40)==0) write(io%uout,*)"  Gaunt R-Gaunt  l1 m1 l2 m2 l3 m3 i  j  k"
            endif
          endif
        enddo
      enddo
    enddo
    if (io%verbosity >= k_highVerbos) then
      write(io%uout,'(/a/)')&
        "==End Gaunt Coefficients============================================="
    endif
  end subroutine InitGaunt
! =====End Gaunt Coefficients ================================================================================
! ====Radial dependence=========================================================================================
  function rad(r,sp1,sp2,gen,tb,l1,l2,m)
    ! Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'rad'
    !--subroutine parameters -------------------------!
    real(k_pr), intent(inout) :: r
    integer , intent(inout) :: sp1,sp2
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    integer, intent(in) :: l1,l2,m
    real(k_pr)             :: rad
    !-------------------------------------------------!
    if (r<=(tb%hopping(sp1,sp2)%r1)) then
      rad=tb%hopping(sp1,sp2)%a(l1,l2,m)*RadNoTail(r,sp1,sp2,gen,tb)
    elseif (r<=(tb%hopping(sp1,sp2)%rcut)) then
      rad=tailFunction(r,tb%hopping(sp1,sp2)%hmnTail(l1,l2,m))
    else
      rad=0.0_k_pr
    endif
  end function rad

  function RadNoTail(r,sp1,sp2,gen,tb)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'rad'
    !--subroutine parameters -------------------------!
    real(k_pr), intent(in) ::r
    integer , intent(in) ::sp1,sp2
    type(modelType),intent(in) :: tb
    type(generalType), intent(in) :: gen

    real(k_pr) :: RadNoTail
    !-------------------------------------------------!
    real(k_pr) ::r0,n,nc,rc,phi0
    select case(gen%bond)
    case(k_bondGSP)
      r0 = tb%hopping(sp1,sp2)%r0
      n  = tb%hopping(sp1,sp2)%n
      nc = tb%hopping(sp1,sp2)%nc
      rc = tb%hopping(sp1,sp2)%rc
      RadNoTail=(r0/r)**n*exp(n*(-(r/rc)**nc+(r0/rc)**nc))
    case (k_bondHarrison)
       n=tb%hopping(sp1,sp2)%n
      RadNoTail = k_hbar**2/k_me/r**(n)
    end select

  end function RadNoTail


  function RadP(alpha,r,sp1,sp2,gen,tb,l1,l2,m)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RadP'
    !--subroutine parameters -------------------------!
    real(k_pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2,alpha
    integer, intent(in) :: l1,l2,m
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    real(k_pr) :: RadP
    !-------------------------------------------------!


    if (r<=tb%hopping(sp1,sp2)%r1) then
      RadP=tb%hopping(sp1,sp2)%a(l1,l2,m)*RadPNoTail(alpha,r,sp1,sp2,gen,tb)
    elseif (r<=tb%hopping(sp1,sp2)%rcut) then
      if (alpha==3) then
        RadP=TailFunctionP(r,tb%hopping(sp1,sp2)%hmnTail(l1,l2,m))
      else
        RadP=0.0_k_pr
      endif
    else
      RadP=0.0_k_pr
    endif


  end function RadP


  function RadPNoTail(alpha,r,sp1,sp2,gen,tb)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RadP_notail'
    !--subroutine parameters -------------------------!
    real(k_pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2
    integer, intent(in) :: alpha
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    real(k_pr) :: RadPNoTail
    !-------------------------------------------------!
    real(k_pr) ::r0,n,nc,rc,phi0

    select case(gen%bond)
    case(k_bondGSP)
      r0 = tb%hopping(sp1,sp2)%r0
      n= tb%hopping(sp1,sp2)%n
      nc=tb%hopping(sp1,sp2)%nc
      rc=tb%hopping(sp1,sp2)%rc
      select case(alpha)
      case (3)
        RadPNoTail=-RadNoTail(r,sp1,sp2,gen,tb)*(1.0_k_pr+nc*(r/rc)**nc)*(n/r)
      case default
        RadPNoTail=0.0_k_pr
      endselect
    case (k_bondHarrison)
      select case(alpha)
      case (3)
        n    = tb%hopping(sp1,sp2)%n
        RadPNoTail = -n/r**(n+1)*k_hbar**2/k_me
      case default
        RadPNoTail=0.0_k_pr
      end select
    end select
  end function RadPNoTail

 function RadPp(alpha,beta,r,sp1,sp2,gen,tb,l1,l2,m)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RadPpNoTail'
    !--subroutine parameters -------------------------!
    real(k_pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2,alpha,beta
    type(modelType),intent(inout) :: tb
    integer, intent(in) :: l1,l2,m
    type(generalType), intent(inout) :: gen
    real(k_pr) :: RadPp
    !-------------------------------------------------!

    if (r<=tb%hopping(sp1,sp2)%r1) then
      RadPp=tb%hopping(sp1,sp2)%a(l1,l2,m)*RadPpNoTail(alpha,beta,r,sp1,sp2,gen,tb)
    elseif (r<=tb%hopping(sp1,sp2)%rcut) then
      if ((alpha==3).and.(beta==3)) then
        RadPp=TailFunctionPp(r,tb%hopping(sp1,sp2)%hmnTail(l1,l2,m))
      else
        RadPp=0.0_k_pr
      endif
    else
      RadPp=0.0_k_pr
    endif



  end function RadPp

  function RadPpNoTail(alpha,beta,r,sp1,sp2,gen,tb)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RadPpNoTail'
    !--subroutine parameters -------------------------!
    real(k_pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2
    integer, intent(in) ::  alpha,beta
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    real(k_pr) :: RadPpNoTail
    !-------------------------------------------------!
    real(k_pr) ::r0,n,nc,rc,phi0

    select case(gen%bond)
    case(k_bondGSP)
      r0 = tb%hopping(sp1,sp2)%r0
      n= tb%hopping(sp1,sp2)%n
      nc=tb%hopping(sp1,sp2)%nc
      rc=tb%hopping(sp1,sp2)%rc
      if ((alpha==3).and.(beta==3)) then
        RadPpNoTail=RadNoTail(r,sp1,sp2,gen,tb) *((-n/r-n*nc*(r/rc)**nc/r)**2-(-n/(r*r)+nc*(nc-1.0_k_pr)*n/(r*r)*(r/rc)**nc))
      else
        RadPpNoTail=0.0_k_pr
      endif
    case(k_bondHarrison)
      if ((alpha==3).and.(beta==3)) then
        n    = tb%hopping(sp1,sp2)%n
        RadPpNoTail = n*(n+1.0_k_pr)/r**(n+2)*k_hbar**2/k_me
      else
        RadPpNoTail = 0.0_k_pr
      endif
    end select
  end function RadPpNoTail
! ============end radial dependence==============================================================================
! ====repulsive term======================================================================================
  function rep(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'rep'
    !--subroutine parameters -------------------------!
    real(k_pr)             :: rep
    real(k_pr), intent(in) :: r
    integer , intent(in) :: sp1,sp2
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    !--internal variables ----------------------------!
    !    real(k_pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!


    if (r<=tb%hopping(sp1,sp2)%d1) then
      rep=RepNoTail(r,sp1,sp2,tb,gen)
    elseif (r<=tb%hopping(sp1,sp2)%dcut) then
      rep = tailFunction(r,tb%hopping(sp1,sp2)%repTail)
    else
      rep = 0.0_k_pr
    endif

  end function rep


  function RepNoTail(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RepNoTail'
    !--subroutine parameters -------------------------!
    real(k_pr)             :: RepNoTail
    real(k_pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    !--internal variables ----------------------------!
    real(k_pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!
    select case(gen%bond)
    case(k_bondGSP)
      phi0 = tb%hopping(sp1,sp2)%phi0
      d0   = tb%hopping(sp1,sp2)%d0
      m    = tb%hopping(sp1,sp2)%m
      mc   = tb%hopping(sp1,sp2)%mc
      dc   = tb%hopping(sp1,sp2)%dc
      RepNoTail = phi0*(d0/r)**m*exp(m*(-(r/dc)**mc+(d0/dc)**mc))

    case (k_bondHarrison)
      phi0   = tb%hopping(sp1,sp2)%phi0
      m    = tb%hopping(sp1,sp2)%m
      RepNoTail = phi0/r**m
    end select
  end function RepNoTail

  function RepP(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RepP'
    !--subroutine parameters -------------------------!
    real(k_pr)             :: RepP
    real(k_pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    !--internal variables ----------------------------!
!    real(k_pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!

    if (r<=tb%hopping(sp1,sp2)%d1) then
      RepP=RepPNoTail(r,sp1,sp2,tb,gen)
    elseif (r<=tb%hopping(sp1,sp2)%dcut) then
      RepP = TailFunctionP(r,tb%hopping(sp1,sp2)%repTail)
    else
      RepP = 0.0_k_pr
    endif

  end function RepP

  function RepPNoTail(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RepPNoTail'
    !--subroutine parameters -------------------------!
    real(k_pr)             :: RepPNoTail
    real(k_pr), intent(in) :: r
    integer ,intent(in)  :: sp1,sp2
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    !--internal variables ----------------------------!
    real(k_pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!
    select case(gen%bond)
    case(k_bondGSP)
      phi0 = tb%hopping(sp1,sp2)%phi0
      d0   = tb%hopping(sp1,sp2)%d0
      m    = tb%hopping(sp1,sp2)%m
      mc   = tb%hopping(sp1,sp2)%mc
      dc   = tb%hopping(sp1,sp2)%dc

      RepPNoTail = RepNoTail(r,sp1,sp2,tb,gen)*(-m/r-m*mc*(r/dc)**mc/r)
    case (k_bondHarrison)

      phi0   = tb%hopping(sp1,sp2)%phi0
      m    = tb%hopping(sp1,sp2)%m

      RepPNoTail = -m*phi0/r**(m+1)
    end select
  end function RepPNoTail


  function RepPp(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RepPp'
    !--subroutine parameters -------------------------!
    real(k_pr)             :: RepPp
    real(k_pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    !--internal variables ----------------------------!
!    real(k_pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!

    if (r<=tb%hopping(sp1,sp2)%d1) then
      RepPp=RepPpNoTail(r,sp1,sp2,tb,gen)
    elseif (r<=tb%hopping(sp1,sp2)%dcut) then
      RepPp = TailFunctionPp(r,tb%hopping(sp1,sp2)%repTail)
    else
      RepPp = 0.0_k_pr
    endif

  end function RepPp


  function RepPpNoTail(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'RepPpNoTail'
    !--subroutine parameters -------------------------!
    real(k_pr)             :: RepPpNoTail
    real(k_pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(modelType),intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    !--internal variables ----------------------------!
    real(k_pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!
    select case(gen%bond)
    case(k_bondGSP)
      phi0 = tb%hopping(sp1,sp2)%phi0
      d0   = tb%hopping(sp1,sp2)%d0
      m    = tb%hopping(sp1,sp2)%m
      mc   = tb%hopping(sp1,sp2)%mc
      dc   = tb%hopping(sp1,sp2)%dc

      RepPpNoTail=RepNoTail(r,sp1,sp2,tb,gen)*((-m/r- m*mc*(r/dc)**mc/r)**2-(-m/(r*r)+mc*(mc-1.0_k_pr)*m/(r*r)*(r/dc)**mc))
    case (k_bondHarrison)

      phi0   = tb%hopping(sp1,sp2)%phi0
      m    = tb%hopping(sp1,sp2)%m

      RepPpNoTail = m*(m+1.0_k_pr)*phi0/r**(m+2)
    end select
  end function RepPpNoTail
! ========end  repulsive term ============================================================================
!======== onsite term======================
  function onsite(orba,orbb,tb)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = "onsite"
    !--subroutine parameters -------------------------!
    real(k_pr) :: onsite
    type(orbitalType) , intent(in) :: orba, orbb
    type(modelType), intent(in) :: tb
    !--internal variables ----------------------------!

    if ((orba%l==orbb%l).and.(orba%m==orbb%m)) then
      onsite=tb%hopping(orba%sp,orba%sp)%eps(orba%l)
    else
      onsite=0.0_k_pr
    endif
  end function onsite
!=====end onsite term ======================
!========embedding term=========================================
  function embedding(x,sp,tb)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = "embedding"
    !--subroutine parameters -------------------------!
    real(k_pr) :: embedding
    real(k_pr) , intent(in):: x
    integer, intent(in) :: sp
    type(modelType), intent(in) :: tb
    !--internal variables ----------------------------!
    real(k_pr) a1,a2,a3,a4

    a1=tb%hopping(sp,sp)%a1
    a2=tb%hopping(sp,sp)%a2
    a3=tb%hopping(sp,sp)%a3
    a4=tb%hopping(sp,sp)%a4

    embedding= a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x

  end function embedding

  function EmbeddingP(x,sp,tb)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = "embedding"
    !--subroutine parameters -------------------------!
    real(k_pr) :: EmbeddingP
    real(k_pr) , intent(in):: x
    integer, intent(in) :: sp
    type(modelType), intent(in) :: tb
    !--internal variables ----------------------------!
    real(k_pr) a1,a2,a3,a4

    a1=tb%hopping(sp,sp)%a1
    a2=tb%hopping(sp,sp)%a2
    a3=tb%hopping(sp,sp)%a3
    a4=tb%hopping(sp,sp)%a4
    EmbeddingP = a1+2.0_k_pr*a2*x+3.0_k_pr*a3*x*x+4.0_k_pr*a4*x*x*x

  end function EmbeddingP

  function EmbeddingPp(x,sp,tb)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'EmbeddingPp'
    !--subroutine parameters -------------------------!
    real(k_pr) :: EmbeddingPp
    real(k_pr) , intent(in):: x
    integer, intent(in) :: sp
    type(modelType), intent(in) :: tb
    !--internal variables ----------------------------!
    real(k_pr) a2,a3,a4


    a2=tb%hopping(sp,sp)%a2
    a3=tb%hopping(sp,sp)%a3
    a4=tb%hopping(sp,sp)%a4

    EmbeddingPp= 2.0_k_pr*a2+6.0_k_pr*a3*x+12.0_k_pr*a4*x**2

  end function EmbeddingPp

!========end embedding term ====================================
!> \brief initializes the magnetic moment
!> \details it just computes the difference between the spin up and spin down charges
!> \author Alin M Elena
!> \date 08/11/07, 11:39:49
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
  subroutine InitMagneticMoment(atomic)
  character(len=*), parameter :: myname="InitMagneticMoment"
    type(atomicxType), intent(inout) :: atomic


    real(k_pr) :: tmom,lmom
    integer :: i,k
    tmom=0.0_k_pr

    do i=1,atomic%atoms%natoms
    lmom=0.0_k_pr
      do k=1,atomic%atoms%norbs(i)/2
        lmom=lmom-(atomic%basis%orbitals(atomic%atoms%orbs(i,k))%occup-&
        atomic%basis%orbitals(atomic%atoms%orbs(i,k+atomic%atoms%norbs(i)/2))%occup)
    enddo
    atomic%atoms%MagMom(i)=lmom
    tmom=tmom+atomic%atoms%MagMom(i)
    enddo
    atomic%atoms%MagneticMoment=tmom
  end subroutine InitMagneticMoment

!> \brief computes the n0
!> \details we suppose that we start with a diagonal "density matrix"
!> \author Alin M Elena
!> \date 08/11/07, 10:41:28
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine Buildn0(gen,atomic,sol)
    character(len=*), parameter :: myname = 'Buildn0'
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    integer :: i,k,j,m

    m=atomic%basis%norbitals/2
    do i=1,atomic%atoms%natoms
      do k=1,atomic%species%norbs(atomic%atoms%sp(i))/2
        j= atomic%atoms%orbs(i,k)
        if (gen%SymRefRho) then
          sol%n0(j) = (atomic%basis%orbitals(j)%occup+ &
          atomic%basis%orbitals(j+m)%occup)/2.0_k_pr&
            -gen%netCharge/atomic%basis%norbitals
          sol%n0(j+m)=sol%n0(j)
        else
          sol%n0(j) = atomic%basis%orbitals(j)%occup-gen%netcharge/atomic%basis%norbitals
          sol%n0(j+m)=atomic%basis%orbitals(j+m)%occup-gen%netcharge/atomic%basis%norbitals
        endif
      end do
    end do
  end subroutine Buildn0

!> \brief computes the magnetic pmoment for each atom
!> \details internally all the "hard" work is done by LocalMoment
!> \author Alin M Elena
!> \date 08/11/07, 11:41:42
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param io type(ioType) contains all the info about I/O files
 subroutine ComputeMagneticMoment(gen,atomic,sol,io)
   character(len=*), parameter :: myname = 'ComputeMagneticMoment'
    type(generalType), intent(in) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(ioType), intent(inout) :: io
    integer :: from,to,i
    real(k_pr) :: lmom,tmom
    tmom=0.0_k_pr
    if (gen%spin) then
    do i=1,atomic%atoms%natoms
      atomic%atoms%MagMom(i)=LocalMoment(atomic%atoms%id(i),atomic,.false.,io,sol)
      tmom=tmom+atomic%atoms%MagMom(i)
    enddo
    atomic%atoms%MagneticMoment=tmom
    else
      atomic%atoms%MagMom=0.0_k_pr
      atomic%atoms%MagneticMoment=0.0_k_pr
    endif
  end subroutine ComputeMagneticMoment
end module m_TightBinding