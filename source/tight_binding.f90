!> \brief deals with the general tight binding functions
!> \author Alin M Elena
!> \date 01/11/07, 11:48:23

module TightBinding
  use constants
  use types
  use useful
  use TailFunctions
  use linear_algebra

  implicit none

  private

  public :: set_solution_space
  public :: rad
  public :: rad_p
  public :: rad_pp
  public :: rep
  public :: rep_p
  public :: rep_pp
  public :: onsite
  public :: embedding
  public :: embedding_p
  public :: embedding_pp

contains
!> \brief initializes the solution space
!> \details allocates memory for all the variables necessary for the calculation. It will
!> also call the routines that precompute some values (eg factorial, Gaunt coefficients)
!> \author Alin M Elena
!> \date 01/11/07, 11:49:38
!> \param io_loc type(io_type) contains all the info about I/O files
!> \param gen_loc type(general_type) contains the info needed by the program to run
!> \param atomic type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains information about the tight binding model parameters
!> \param sol type(solution_type) contains information about the solution space

  subroutine set_solution_space(io_loc,gen_loc,atomic,tbMod,sol)
    character(len=*), parameter :: sMyName="set_solution_space"
    type(io_type), intent(inout) :: io_loc
    type(general_type), intent(inout) :: gen_loc
    type(atomicx_type), intent(inout) :: atomic
    type(model_type), intent(inout) :: tbMod
    type(solution_type), intent(inout) :: sol
    integer :: l

    if (.not.sol%h%created) then
      call create_sparse_matrix(sol%h,atomic%basis%norbitals,.true.)
    else
      call reset_sparse_matrix(sol%h)
    endif
    if (.not.sol%eigenvecs%created) then
      allocate(sol%eigenvals(1:atomic%basis%norbitals))
      call create_matrix(sol%eigenvecs,atomic%basis%norbitals,.true.)
    else
      sol%eigenvals=0.0_pr
      call zero_matrix(sol%eigenvecs,io_loc)
    endif

    if (gen_loc%spin) then
       !spin down
      if (.not.sol%hdown%created) then
        call create_matrix(sol%hdown,atomic%basis%norbitals/2,.true.)
      else
        call zero_matrix(sol%hdown,io_loc)
      endif
      if (.not.sol%hup%created) then
        call create_matrix(sol%hup,atomic%basis%norbitals/2,.true.)
      else
        call zero_matrix(sol%hup,io_loc)
      endif
    endif
    if (.not.sol%rho%created) then
      call create_matrix(sol%rho,atomic%basis%norbitals,.true.)
    else
      call zero_matrix(sol%rho,io_loc)
    endif

    l= maxval(tbMod%hopping(:,:)%l2)
    allocate(sol%fact(0:(6*(l+1)+1)))
    call init_fact(6*(l+1)+1,sol%fact)

    if ((gen_loc%scf).and.(gen_loc%electrostatics==electrostatics_multipoles)) then
      call init_gaunt(2*l+1,2*l+1,4*l+2,sol,io_loc)
 !     call build_n0(basis_var%norbital)
  !    call build_hn0(basis_var%norbital)
    endif
    call set_tails(io_loc,gen_loc,atomic,tbMod,sol)
  end subroutine set_solution_space



!> \brief sets the tails parameters
!> \author Alin M Elena
!> \date 02/11/07, 16:18:21
!> \param io_loc type(io_type) contains all the info about I/O files
!> \param gen_loc type(general_type) contains the info needed by the program to run
!> \param atomic type(atomic_type) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(model_type) contains information about the tight binding model parameters
!> \param sol type(solution_type) contains information about the solution space

  subroutine set_tails(io_loc,gen_loc,atomic,tbMod,sol)
    character(len=*), parameter :: sMyName="set_tails"
    type(io_type), intent(inout) :: io_loc
    type(general_type), intent(inout) :: gen_loc
    type(atomicx_type), intent(inout) :: atomic
    type(model_type), intent(inout) :: tbMod
    type(solution_type), intent(inout) :: sol

    integer ::i,j
    real(pr) :: f,fp,fpp
    integer :: k,k1,k2
!     real(pr) :: r
!     integer :: z,y
  ! calculate the tail function parameters
    do i=1,atomic%species%nspecies
      do j=1,atomic%species%nspecies
          ! repulsions
        f = rep_notail(tbMod%hopping(i,j)%d1,atomic%species%id(i),atomic%species%id(j),tbMod,gen_loc)
        fp = rep_p_notail(tbMod%hopping(i,j)%d1,atomic%species%id(i),atomic%species%id(j),tbMod,gen_loc)
        fpp = rep_pp_notail(tbMod%hopping(i,j)%d1,atomic%species%id(i),atomic%species%id(j),tbMod,gen_loc)
        tbMod%hopping(i,j)%rep_tail=make_tail_parameters(f,fp,fpp,tbMod%hopping(i,j)%d1,tbMod%hopping(i,j)%dcut,io_loc)

          ! hoppings
        allocate(tbMod%hopping(i,j)%hmn_tail(0:tbMod%hopping(i,j)%l1,0:tbMod%hopping(i,j)%l2,0:tbMod%hopping(i,j)%ll))
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
!               write(io%uout,*) &
!               trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2)
              f = tbMod%hopping(i,j)%a(k,k1,k2)*rad_no_tail(tbMod%hopping(i,j)%r1,&
                    atomic%species%id(i),atomic%species%id(j),gen_loc,tbMod)
              fp = tbMod%hopping(i,j)%a(k,k1,k2)*rad_p_no_tail(3,tbMod%hopping(i,j)%r1,&
                    atomic%species%id(i),atomic%species%id(j),gen_loc,tbMod)
              fpp = tbMod%hopping(i,j)%a(k,k1,k2)*rad_pp_no_tail(3,3,tbMod%hopping(i,j)%r1,&
                    atomic%species%id(i),atomic%species%id(j),gen_loc,tbMod)
              tbMod%hopping(i,j)%hmn_tail(k,k1,k2) = &
              make_tail_parameters(f,fp,fpp,tbMod%hopping(i,j)%r1,tbMod%hopping(i,j)%rcut,io_loc)
            enddo
          enddo
        enddo

      enddo
    enddo

    call print_tail_parameters(atomic,tbMod,io_loc)
!     z=600
!     do i=1,atomic%species%nspecies
!       do j=1,atomic%species%nspecies
!         do k=0,tbMod%hopping(i,j)%l1
!           do k1=0,tbMod%hopping(i,j)%l2
!             do k2=0,min(k,k1)
!               write(z,*)"# r ",trim(ccnlm(i,j,k,k1,k2))
!               do y=50,1000
!                 r=y*0.01_pr
!                 write(z,*)r,rad(r,atomic%species%id(i),atomic%species%id(j),gen_loc,tbMod,k,k1,k2)
!               enddo
!               z=z+1
!             enddo
!           enddo
!         enddo
!       enddo
!     enddo
end subroutine set_tails




                
!==Gaunt coefficients====================================================
!> \brief returns the Gaunt coefficient for l1,m1,l2,m2,l3,m3
!> \author Alin M Elena
!> \date 01/11/07, 17:15:07
!> \param l1,m1,l2,m2,l3,m3 angular momentum and quantum magnetic numbers
!> \param io type(io_type) contains all the info about I/O files
!> \param sol type(solution_type) contains information about the solution space
  real(pr) function cgaunt(l1,m1,l2,m2,l3,m3,sol,io)
    character(len=*),parameter :: myname="cgaunt"
    integer, intent(in) :: l1,m1,l2,m2,l3,m3
    type(io_type), intent(in) :: io
    type(solution_type), intent(in) :: sol
    real(pr) :: sq
    integer ::g
    character(len=ml) :: saux
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
      sq=real((2*l1+1)*(2*l2+1),pr)
      sq=sqrt(sq/(4.0_pr*pi))*sol%fact(g)/(sol%fact(g-l1)*sol%fact(g-l2)*sol%fact(g-l3))
      cgaunt=(-1)**(g-l3)*sq*cgc(l1,m1,l2,m2,l3,m3,sol)*dl(l1,l2,l3,sol)
    else
      cgaunt=0.0_pr
    endif

  end function cgaunt

!> \brief the "heart" of cgaunt
!> \author Alin M Elena
!> \date 01/11/07, 17:20:53
!> \param a,x,b,y,c,z
!> \param sol type(solution_type) contains information about the solution space
  real(pr) function cgc(a,x,b,y,c,z,sol)
    integer, intent(in) :: a,x,b,y,c,z
    type(solution_type), intent(in) :: sol
    integer :: kmin,kmax,k
    real(pr) :: sum,gamma

    gamma=sqrt(sol%fact(a+x)*sol%fact(b+y)*sol%fact(c+z)* &
           sol%fact(a-x)*sol%fact(b-y)*sol%fact(c-z)*real(2*c+1,pr))
    sum=0.0_pr
    kmin=max(b-c-x,a+y-c,0)
    kmax=min(b+y,a-x,a+b-c)
    do k=kmin,kmax
    sum=sum+(-1)**k/(sol%fact(k)*sol%fact(a+b-c-k)*sol%fact(a-x-k)* &
        sol%fact(b+y-k)*sol%fact(c-b+x+k)*sol%fact(c-a-y+k))
    enddo
    cgc=dl(a,b,c,sol)*gamma*sum
  end function cgc

  real(pr) function dl(l1,l2,l3,sol)
    integer, intent(in):: l1,l2,l3
    type(solution_type), intent(in) :: sol
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

  complex(pr) function ulmmiu(m,miu)
    integer :: m,miu,z
    real(pr) :: sq2,r,i

    sq2=dsqrt(2.0_pr)
    z=0
    r=d(m,z)*d(miu,z)+h(miu)*(d(m,miu)+(-1)**m*d(m,-miu))/sq2
    i=h(-miu)*(-d(m,-miu)+d(m,miu)*(-1)**m)/sq2

    ulmmiu=cmplx(r,i,pr)
  end function ulmmiu

!> \brief Builds the Gaunt coefficients tables
!> \author Alin M Elena
!> \param l1m,l2m,l3m integers maxima ls for angular momenta
!> \param io type(io_type) contains all the info about I/O files
!> \param sol type(solution_type) contains information about the solution space
  subroutine init_gaunt(l1m,l2m,l3m,sol,io)
    character(len=*),parameter :: myname="init_gaunt"
    integer,intent(in) :: l1m,l2m,l3m
    type(io_type), intent(in) :: io
    type(solution_type), intent(inout) :: sol
    integer :: n1,n2,n3
    integer :: i,j,k,l1,l2,l3,m1,m2,m3,ii,kk,jj,p=0
    complex(pr) :: aux,tt

    n1=l1m*l1m+2*l1m+1
    n2=l2m*l2m+2*l2m+1
    n3=l3m*l3m+2*l3m+1
    allocate(sol%gcoeff(1:n1,1:n2,1:n3),sol%rgc(1:n1,1:n2,1:n3))

    sol%gcoeff=0.0_pr
    sol%rgc=0.0_pr
    
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

    if (io%verbosity > high_verbos) then
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
          tt=cmplx(0.0_pr,0.0_pr,pr)

          do ii=-l1,l1
            do jj=-l2,l2
              do kk=-l3,l3
                if (sol%gcoeff(idx(l1,ii),idx(l2,jj),idx(l3,kk))/=zero) then
                  aux=ulmmiu(ii,m1)
                  tt=tt+(-1)**kk*ulmmiu(ii,m1)*ulmmiu(jj,m2)*ulmmiu(-kk,m3)*sol%gcoeff(idx(l1,ii),idx(l2,jj),idx(l3,kk))
                endif
              enddo
            enddo
          enddo
          sol%rgc(i,j,k)=real(tt,pr)
          if (io%verbosity > high_verbos) then
            if ((sol%gcoeff(i,j,k)/=zero).or.(sol%rgc(i,j,k)/=zero)) then
              write(io%uout,'(2F8.4,6I3,3I3)')sol%gcoeff(i,j,k),sol%rgc(i,j,k),l1,m1,l2,m2,l3,m3,i,j,k
              p=p+1
              if (mod(p,40)==0) write(io%uout,*)"  Gaunt R-Gaunt  l1 m1 l2 m2 l3 m3 i  j  k"
            endif
          endif  
        enddo
      enddo
    enddo
    if (io%verbosity > high_verbos) then
      write(io%uout,'(/a/)')&
        "==End Gaunt Coefficients============================================="
    endif
  end subroutine init_gaunt
! =====End Gaunt Coefficients ================================================================================
! ====Radial dependence=========================================================================================
  function rad(r,sp1,sp2,gen,tb,l1,l2,m)
    ! Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rad'
    !--subroutine parameters -------------------------!
    real(pr), intent(inout) :: r
    integer , intent(inout) :: sp1,sp2
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    integer, intent(in) :: l1,l2,m
    real(pr)             :: rad
    !-------------------------------------------------!

    if (r<=(tb%hopping(sp1,sp2)%r1)) then
      rad=tb%hopping(sp1,sp2)%a(l1,l2,m)*rad_no_tail(r,sp1,sp2,gen,tb)
    elseif (r<=(tb%hopping(sp1,sp2)%rcut)) then
      rad=tail_function(r,tb%hopping(sp1,sp2)%hmn_tail(l1,l2,m))
    else
      rad=0.0_pr
    endif
  end function rad

  function rad_no_tail(r,sp1,sp2,gen,tb)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rad'
    !--subroutine parameters -------------------------!
    real(pr), intent(in) ::r
    integer , intent(in) ::sp1,sp2
    type(model_type),intent(in) :: tb
    type(general_type), intent(in) :: gen
  
    real(pr) :: rad_no_tail
    !-------------------------------------------------!
    real(pr) ::r0,n,nc,rc,phi0
    select case(gen%bond)
    case(bond_gsp)
      r0 = tb%hopping(sp1,sp2)%r0
      n  = tb%hopping(sp1,sp2)%n
      nc = tb%hopping(sp1,sp2)%nc 
      rc = tb%hopping(sp1,sp2)%rc
      rad_no_tail=(r0/r)**n*exp(n*(-(r/rc)**nc+(r0/rc)**nc))
    case (bond_harrison)
       n=tb%hopping(sp1,sp2)%n
      rad_no_tail = hbar**2/me/r**(n)
    end select

  end function rad_no_tail


  function rad_p(alpha,r,sp1,sp2,gen,tb,l1,l2,m)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rad_p'
    !--subroutine parameters -------------------------!
    real(pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2,alpha
    integer, intent(in) :: l1,l2,m
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    real(pr) :: rad_p
    !-------------------------------------------------!


    if (r<=tb%hopping(sp1,sp2)%r1) then
      rad_p=tb%hopping(sp1,sp2)%a(l1,l2,m)*rad_p_no_tail(alpha,r,sp1,sp2,gen,tb)
    elseif (r<=tb%hopping(sp1,sp2)%rcut) then
      if (alpha==3) then 
        rad_p=tail_function_p(r,tb%hopping(sp1,sp2)%hmn_tail(l1,l2,m))
      else
        rad_p=0.0_pr
      endif
    else
      rad_p=0.0_pr
    endif


  end function rad_p


  function rad_p_no_tail(alpha,r,sp1,sp2,gen,tb)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rad_p_notail'
    !--subroutine parameters -------------------------!
    real(pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2
    integer, intent(in) :: alpha
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    real(pr) :: rad_p_no_tail
    !-------------------------------------------------!
    real(pr) ::r0,n,nc,rc,phi0

    select case(gen%bond)
    case(bond_gsp)
      r0 = tb%hopping(sp1,sp2)%r0
      n= tb%hopping(sp1,sp2)%n
      nc=tb%hopping(sp1,sp2)%nc
      rc=tb%hopping(sp1,sp2)%rc
      select case(alpha)
      case (3)
        rad_p_no_tail=-rad_no_tail(r,sp1,sp2,gen,tb)*(1.0_pr+nc*(r/rc)**nc)*(n/r)
      case default
        rad_p_no_tail=0.0_pr
      endselect
    case (bond_harrison)
      select case(alpha)
      case (3)
        n    = tb%hopping(sp1,sp2)%n
        rad_p_no_tail = -n/r**(n+1)*hbar**2/me
      case default
        rad_p_no_tail=0.0_pr
      end select
    end select
  end function rad_p_no_tail

 function rad_pp(alpha,beta,r,sp1,sp2,gen,tb,l1,l2,m)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rad_pp_no_tail'
    !--subroutine parameters -------------------------!
    real(pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2,alpha,beta
    type(model_type),intent(inout) :: tb
    integer, intent(in) :: l1,l2,m
    type(general_type), intent(inout) :: gen
    real(pr) :: rad_pp
    !-------------------------------------------------!

    if (r<=tb%hopping(sp1,sp2)%r1) then
      rad_pp=tb%hopping(sp1,sp2)%a(l1,l2,m)*rad_pp_no_tail(alpha,beta,r,sp1,sp2,gen,tb)
    elseif (r<=tb%hopping(sp1,sp2)%rcut) then
      if ((alpha==3).and.(beta==3)) then      
        rad_pp=tail_function_pp(r,tb%hopping(sp1,sp2)%hmn_tail(l1,l2,m))
      else
        rad_pp=0.0_pr
      endif
    else
      rad_pp=0.0_pr
    endif



  end function rad_pp

  function rad_pp_no_tail(alpha,beta,r,sp1,sp2,gen,tb)
    !Gives you the radial dependence for elemental overlaps
    ! for specie sp1, sp2
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rad_pp_no_tail'
    !--subroutine parameters -------------------------!
    real(pr), intent(inout) ::r
    integer , intent(inout) ::sp1,sp2
    integer, intent(in) ::  alpha,beta
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    real(pr) :: rad_pp_no_tail
    !-------------------------------------------------!
    real(pr) ::r0,n,nc,rc,phi0

    select case(gen%bond)
    case(bond_gsp)
      r0 = tb%hopping(sp1,sp2)%r0
      n= tb%hopping(sp1,sp2)%n
      nc=tb%hopping(sp1,sp2)%nc 
      rc=tb%hopping(sp1,sp2)%rc
      if ((alpha==3).and.(beta==3)) then
        rad_pp_no_tail=rad_no_tail(r,sp1,sp2,gen,tb) *((-n/r-n*nc*(r/rc)**nc/r)**2-(-n/(r*r)+nc*(nc-1.0_pr)*n/(r*r)*(r/rc)**nc))
      else
        rad_pp_no_tail=0.0_pr
      endif
    case(bond_harrison)
      if ((alpha==3).and.(beta==3)) then
        n    = tb%hopping(sp1,sp2)%n
        rad_pp_no_tail = n*(n+1.0_pr)/r**(n+2)*hbar**2/me
      else
        rad_pp_no_tail = 0.0_pr
      endif
    end select
  end function rad_pp_no_tail
! ============end radial dependence==============================================================================
! ====repulsive term======================================================================================
  function rep(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rep'
    !--subroutine parameters -------------------------!
    real(pr)             :: rep
    real(pr), intent(in) :: r
    integer , intent(in) :: sp1,sp2
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    !--internal variables ----------------------------!  
    !    real(pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!


    if (r<=tb%hopping(sp1,sp2)%d1) then
      rep=rep_notail(r,sp1,sp2,tb,gen)
    elseif (r<=tb%hopping(sp1,sp2)%dcut) then
      rep = tail_function(r,tb%hopping(sp1,sp2)%rep_tail)
    else
      rep = 0.0_pr
    endif

  end function rep


  function rep_notail(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rep_notail'
    !--subroutine parameters -------------------------!
    real(pr)             :: rep_notail
    real(pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    !--internal variables ----------------------------!  
    real(pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!
    select case(gen%bond)
    case(bond_gsp)
      phi0 = tb%hopping(sp1,sp2)%phi0
      d0   = tb%hopping(sp1,sp2)%d0
      m    = tb%hopping(sp1,sp2)%m
      mc   = tb%hopping(sp1,sp2)%mc
      dc   = tb%hopping(sp1,sp2)%dc
      rep_notail = phi0*(d0/r)**m*exp(m*(-(r/dc)**mc+(d0/dc)**mc))

    case (bond_harrison)
      phi0   = tb%hopping(sp1,sp2)%phi0
      m    = tb%hopping(sp1,sp2)%m
      rep_notail = phi0/r**m
    end select
  end function rep_notail

  function rep_p(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rep_p'
    !--subroutine parameters -------------------------!
    real(pr)             :: rep_p
    real(pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    !--internal variables ----------------------------!  
!    real(pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!

    if (r<=tb%hopping(sp1,sp2)%d1) then
      rep_p=rep_p_notail(r,sp1,sp2,tb,gen)
    elseif (r<=tb%hopping(sp1,sp2)%dcut) then
      rep_p = tail_function_p(r,tb%hopping(sp1,sp2)%rep_tail)
    else
      rep_p = 0.0_pr
    endif

  end function rep_p

  function rep_p_notail(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rep_p_notail'
    !--subroutine parameters -------------------------!
    real(pr)             :: rep_p_notail
    real(pr), intent(in) :: r
    integer ,intent(in)  :: sp1,sp2
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    !--internal variables ----------------------------!  
    real(pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!
    select case(gen%bond)
    case(bond_gsp)
      phi0 = tb%hopping(sp1,sp2)%phi0
      d0   = tb%hopping(sp1,sp2)%d0
      m    = tb%hopping(sp1,sp2)%m
      mc   = tb%hopping(sp1,sp2)%mc 
      dc   = tb%hopping(sp1,sp2)%dc

      rep_p_notail = rep_notail(r,sp1,sp2,tb,gen)*(-m/r-m*mc*(r/dc)**mc/r)
    case (bond_harrison)

      phi0   = tb%hopping(sp1,sp2)%phi0
      m    = tb%hopping(sp1,sp2)%m

      rep_p_notail = -m*phi0/r**(m+1)
    end select
  end function rep_p_notail


  function rep_pp(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rep_pp'
    !--subroutine parameters -------------------------!
    real(pr)             :: rep_pp
    real(pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    !--internal variables ----------------------------!  
!    real(pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!

    if (r<=tb%hopping(sp1,sp2)%d1) then
      rep_pp=rep_pp_notail(r,sp1,sp2,tb,gen)
    elseif (r<=tb%hopping(sp1,sp2)%dcut) then
      rep_pp = tail_function_pp(r,tb%hopping(sp1,sp2)%rep_tail)
    else
      rep_pp = 0.0_pr
    endif

  end function rep_pp


  function rep_pp_notail(r,sp1,sp2,tb,gen)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rep_pp_notail'
    !--subroutine parameters -------------------------!
    real(pr)             :: rep_pp_notail
    real(pr), intent(in) :: r
    integer, intent(in)  :: sp1,sp2
    type(model_type),intent(inout) :: tb
    type(general_type), intent(inout) :: gen
    !--internal variables ----------------------------!  
    real(pr) :: phi0,d0,m,dc,mc
    !-------------------------------------------------!
    select case(gen%bond)
    case(bond_gsp)
      phi0 = tb%hopping(sp1,sp2)%phi0
      d0   = tb%hopping(sp1,sp2)%d0
      m    = tb%hopping(sp1,sp2)%m
      mc   = tb%hopping(sp1,sp2)%mc 
      dc   = tb%hopping(sp1,sp2)%dc

      rep_pp_notail=rep_notail(r,sp1,sp2,tb,gen)*((-m/r- m*mc*(r/dc)**mc/r)**2-(-m/(r*r)+mc*(mc-1.0_pr)*m/(r*r)*(r/dc)**mc))
    case (bond_harrison)

      phi0   = tb%hopping(sp1,sp2)%phi0
      m    = tb%hopping(sp1,sp2)%m

      rep_pp_notail = m*(m+1.0_pr)*phi0/r**(m+2)
    end select
  end function rep_pp_notail
! ========end  repulsive term ============================================================================
!======== onsite term======================
  function onsite(orba,orbb,tb)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = "onsite"
    !--subroutine parameters -------------------------!
    real(pr) :: onsite
    type(orbital_type) , intent(in) :: orba, orbb
    type(model_type), intent(in) :: tb
    !--internal variables ----------------------------! 
    
    if ((orba%l==orbb%l).and.(orba%m==orbb%m)) then
      onsite=tb%hopping(orba%sp,orba%sp)%eps(orba%l)
    else
      onsite=0.0_pr
    endif
  end function onsite
!=====end onsite term ======================
!========embedding term=========================================
  function embedding(x,sp,tb)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = "embedding"
    !--subroutine parameters -------------------------!
    real(pr) :: embedding
    real(pr) , intent(in):: x
    integer, intent(in) :: sp
    type(model_type), intent(in) :: tb
    !--internal variables ----------------------------! 
    real(pr) a1,a2,a3,a4

    a1=tb%hopping(sp,sp)%a1
    a2=tb%hopping(sp,sp)%a2
    a3=tb%hopping(sp,sp)%a3
    a4=tb%hopping(sp,sp)%a4

    embedding= a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x

  end function embedding

  function embedding_p(x,sp,tb)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = "embedding"
    !--subroutine parameters -------------------------!
    real(pr) :: embedding_p
    real(pr) , intent(in):: x
    integer, intent(in) :: sp
    type(model_type), intent(in) :: tb
    !--internal variables ----------------------------! 
    real(pr) a1,a2,a3,a4

    a1=tb%hopping(sp,sp)%a1
    a2=tb%hopping(sp,sp)%a2
    a3=tb%hopping(sp,sp)%a3
    a4=tb%hopping(sp,sp)%a4

    embedding_p = a1+2.0_pr*a2*x+3.0_pr*a3*x*x+4.0_pr*a4*x*x*x

  end function embedding_p

  function embedding_pp(x,sp,tb)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'embedding_pp'
    !--subroutine parameters -------------------------!
    real(pr) :: embedding_pp
    real(pr) , intent(in):: x
    integer, intent(in) :: sp
    type(model_type), intent(in) :: tb
    !--internal variables ----------------------------! 
    real(pr) a2,a3,a4


    a2=tb%hopping(sp,sp)%a2
    a3=tb%hopping(sp,sp)%a3
    a4=tb%hopping(sp,sp)%a4

    embedding_pp= 2.0_pr*a2+6.0_pr*a3*x+12.0_pr*a4*x**2

  end function embedding_pp

!========end embedding term ====================================
end module TightBinding