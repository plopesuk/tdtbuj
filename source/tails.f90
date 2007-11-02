!> \brief implements all routines and functions related with tails for hopping
!> integrals and repulsive potentials
!> \author Cristian G. Sanchez and Alin M Elena (Queen's University Belfast)
!> date 1st of November 2007


module TailFunctions
  use constants
  use types
  use useful
  implicit none
  private


  public :: make_tail_parameters
  public :: print_tail_parameters
  public :: tail_function
  public :: tail_function_p
  public :: tail_function_pp

contains


!> \brief returns  tail_parameters
!> given the cutoffs and values of the function
!> derivatives at the inner cutoff
!>\details the expression of the tail function \f[
!>          t(r)=(r-r_{cut})^3(10c+5br+3ar^2+5br_{cut}+4arr_{cut}+3ar_{cut}^2)/60
!>          \f] with \f$ a,b,c\f$ given by
!> \f[
!>a=\frac{10}{(r_1-r_{cut})^5}(12A-6Br_1+Cr_1^2+6Br_{cut}-2Cr_1r_{cut}+Cr_{cut}^2)
!>\f]
!>\f[
!> b=-\frac{4}{(r_1-r_{cut})^5}(45Ar_1-21Br_1^2+3Cr_1^3+15Ar_{cut}+12Br_1r_{cut}-4Cr_1^2r_{cut}+
!> 9Br_{cut}^2-Cr_1r_{cut}^2+2Cr_{cut}^3)
!>\f]
!>\f[
!>c=-\frac{1}{(r_1-r_{cut})^5}(-60Ar_1^2+24Br_1^3-3Cr_1^4-60Ar_1r_{cut}+12Br_1^2r_{cut}-
!>36Br_1r_{cut}^2+8Cr_1^2r_{cut}^2-4Cr_1r_{cut}^3-Cr_{cut}^4)
!>\f] with \f$A=f(r_1)\f$, \f$B=f^{\prime}(r_1)\f$, \f$C=f^{\prime\prime}(r_1)\f$ \f$f(r)\f$ is the function that we tail
!> \f$ r_1\f$, \f$ r_{cut} \f$ inner and outer cuttoff
!> \author Cristian G. Sanchez
!> \date ~2004-2005
!> \param f,fp,fpp reals the function and its derivatives at inner cuttoff
!> \param r_in,r_out reals inner and outer cuttoff
!> \param io type(io_type) i/o units
  function make_tail_parameters(f,fp,fpp,r_in,r_out,io)
  !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'make_tail_parameters'
  !--subroutine parameters -------------------------!
      type(tail_parameters) :: make_tail_parameters
      real(pr) :: f, fp, fpp
      real(pr) :: r_in, r_out
      real(pr) :: x, x2, a, b, c
      type(io_type), intent(inout) :: io
  !--internal variables ----------------------------!
      logical :: status
      character(len=ml) :: saux
  !-------------------------------------------------!

      make_tail_parameters%r_in  = r_in

      x = r_in
      x2 = r_out
      a = (10*(12*f - 6*fp*x + fpp*x**2 + 6*fp*x2 - 2*fpp*x*x2 + fpp*x2**2))/(x - x2)**5
      b = (-4*(45*f*x - 21*fp*x**2 + 3*fpp*x**3 + 15*f*x2 + 12*fp*x*x2 - 4*fpp*x**2*x2 + 9*fp*x2**2 - fpp*x*x2**2 + &
          2*fpp*x2**3))/(x - x2)**5
      c =  -((-60*f*x**2 + 24*fp*x**3 - 3*fpp*x**4 - 60*f*x*x2 + 12*fp*x**2*x2 - 36*fp*x*x2**2 + 8*fpp*x**2*x2**2 - &
          4*fpp*x*x2**3 - fpp*x2**4)/(x - x2)**5)

      status = check_tail_param(f,fp,fpp,x,x2)
      if (.not.status) then
        write(saux,'(a,f0.8,1x,f0.8)') &
      "Second derivative of tail function has solution in ri,rcut region, you should increase rcut", x,x2
        call error(trim(saux),myname,.false.,io)
      endif
      make_tail_parameters%a = a
      make_tail_parameters%b = b
      make_tail_parameters%c = c
      make_tail_parameters%r_out = x2

  end function make_tail_parameters


!> \brief prints the tail parameters
!> \author Alin M Elena
!> \date 02/11/07, 09:00:08
!> \param io type(io_type) i/o units
!> \param tail type(tail_parameters) tail parameters
!> \param atomix type(atomicx_type) atomic details

  subroutine print_tail_parameters(atomix,tbmod,io)
    character(len=*), parameter :: sMyName="print_tail_parameters"
    type(io_type), intent(in) :: io
    type(atomicx_type), intent(in) :: atomix
    type(model_type),intent(in) :: tbmod

    integer i,j
    character(len=mw) :: saux
    integer :: k,k1,k2
    write(io%uout,'(a)')"==Tail functions========================================================================================"
    write(io%uout,'(a)')"==Repulsive terms tail parameters======================================================"
        write(io%uout,'(2a4,3a20,2a10)')"Sp1","Sp2","a","b","c","r_1","r_cut"
    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        write(io%uout,'(2i4,3f20.8,2f10.4)')i,j,tbMod%hopping(i,j)%rep_tail%a,&
          tbMod%hopping(i,j)%rep_tail%b,tbMod%hopping(i,j)%rep_tail%c,&
          tbMod%hopping(i,j)%rep_tail%r_in,tbMod%hopping(i,j)%rep_tail%r_out
      enddo
    enddo
        write(io%uout,'(a)')"==End Repulsive terms tail parameters================================================"
    write(io%uout,'(a)')"==Radial terms tail parameters====================================================================="
        write(io%uout,'(a10,9x,3a20,2a10)')"hopping","a","b","c","r_1","r_cut"
    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
              write(io%uout,'(a9,f10.4,3f20.8,2f10.4)')trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2),&
                  tbMod%hopping(i,j)%hmn_tail(k,k1,k2)%a,&
                  tbMod%hopping(i,j)%hmn_tail(k,k1,k2)%b,tbMod%hopping(i,j)%hmn_tail(k,k1,k2)%c,&
                  tbMod%hopping(i,j)%hmn_tail(k,k1,k2)%r_in,tbMod%hopping(i,j)%hmn_tail(k,k1,k2)%r_out
            enddo
          enddo
        enddo
      enddo
    enddo
    write(io%uout,'(a)')"==End Radial terms tail parameters================================================================="
    write(io%uout,'(a)')"==End Tail functions===================================================================================="
  end subroutine print_tail_parameters


!> \brief checks the tail parameters,  
!> \details to see wether there is an inflexion point in the
!> tail region causing a maximum in the force
!> \author Alin M. Elena,Queen's University Belfast, UK
!> \date 11th of February, 2005
!> \param f,fp,fpp reals the function and its derivatives at inner cuttoff
!> \param ri,rc reals inner and outer cuttoff

  function check_tail_param(f,fp,fpp,ri,rc)
!--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'check_tail_param'
!--subroutine parameters -------------------------!
    logical :: check_tail_param
    real(pr) :: f, fp, fpp
    real(pr) :: ri, rc
    !--internal variables ----------------------------!
    real(pr) :: aux1,aux2,x1,x2,sw,aux3
!-------------------------------------------------!

    aux2 = 2.0_pr*(-120.0_pr*f - 60.0_pr*fp*rc - 10.0_pr*fpp*rc**2 + 60.0_pr*fp*ri +&
                    20.0_pr*fpp*rc*ri - 10.0_pr*fpp*ri**2)

    aux1 = (60.0_pr*f*rc + 36.0_pr*fp*rc**2 + 8.0_pr*fpp*rc**3 + 180.0_pr*f*ri + 48.0_pr*fp*rc*ri - &
        4.0_pr*fpp*rc**2*ri - 84.0_pr*fp*ri**2 - 16.0_pr*fpp*rc*ri**2 + 12.0_pr*fpp*ri**3)**2 - &
        4.0_pr*(-120.0_pr*f - 60.0_pr*fp*rc - 10.0_pr*fpp*rc**2 + 60.0_pr*fp*ri + &
        20.0_pr*fpp*rc*ri - 10.0_pr*fpp*ri**2)*&
        (-(fpp*rc**4) - 60.0_pr*f*rc*ri - 36.0_pr*fp*rc**2*ri - 4.0_pr*fpp*rc**3*ri - 60.0_pr*f*ri**2 +&
    12.0_pr*fp*rc*ri**2 + 8.0_pr*fpp*rc**2*ri**2 + 24.0_pr*fp*ri**3 - 3.0_pr*fpp*ri**4)

    aux3 = -60.0_pr*f*rc - 36.0_pr*fp*rc**2 - 8.0_pr*fpp*rc**3 - 180.0_pr*f*ri - 48.0_pr*fp*rc*ri +&
        4.0_pr*fpp*rc**2*ri + 84.0_pr*fp*ri**2 + 16.0_pr*fpp*rc*ri**2 - 12.0_pr*fpp*ri**3

    if (aux1<0.0_pr) then
      check_tail_param=.true.
    else if (abs(aux2)<epsilon(aux2)) then
      check_tail_param=.true.
    else
      x1 = (aux3 + sqrt(aux1))/aux2
      x2 = (aux3 - sqrt(aux1))/aux2

      if (x1>x2) then
        sw=x1
        x1=x2
        x2=sw
      end if

      if ((ri>x1).and.(rc>x2)) then
        check_tail_param=.true.
      else if (ri<x2) then
        check_tail_param=.true.
      else
        if (rc<x1) then
          check_tail_param=.true.
        else
          check_tail_param=.false.
        end if
      end if
    end if

  end function check_tail_param


!>\brief returns the value of the tail
!> function at r, given the parameters in tailp
!>\author Cristian G. Sanchez
!> \date ~2004-2005
!> \param r real the point where the function is evaluated
!> \param tailp the tail parameters
  function tail_function(r,tailp)
  !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'tail_function'
!--subroutine parameters -------------------------!
    real(pr) :: tail_function
    type(tail_parameters) :: tailp
    real(pr) :: r
!--internal variables ----------------------------!
!-------------------------------------------------!

    tail_function= ((r - tailp%r_out)**3*(10*tailp%c + &
        5*tailp%b*r + 3*tailp%a*r**2 + 5*tailp%b*tailp%r_out + &
        4*tailp%a*r*tailp%r_out + 3*tailp%a*tailp%r_out**2))/60.0_pr

  end function tail_function


!> \brief returns a the value of the tail function first derivative at r
!> \author Cristian G. Sanchez
!> \date ~2004-2005
!> \param r real the point where the function is evaluated
!> \param tailp the tail parameters
  function tail_function_p(r,tailp)
  !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'tail_function_p'
!--subroutine parameters -------------------------!
    real(pr) :: tail_function_p
    type(tail_parameters) :: tailp
    real(pr) :: r
!--internal variables ----------------------------!
!-------------------------------------------------!

    tail_function_p =  ((r - tailp%r_out)**2*(6*tailp%c + &
        4*tailp%b*r + 3*tailp%a*r**2 + 2*tailp%b*tailp%r_out &
        + 2*tailp%a*r*tailp%r_out + tailp%a*tailp%r_out**2))/12.0_pr

  end function tail_function_p

!> \brief returns a the value of the tail function second derivative at r
!> \author Cristian G. Sanchez
!> \date ~2004-2005
!> \param r real the point where the function is evaluated
!> \param tailp the tail parameters
  function tail_function_pp(r,tailp)
  !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'tail_function_pp'
!--subroutine parameters -------------------------!
    real(pr) :: tail_function_pp
    type(tail_parameters) :: tailp
    real(pr) :: r
!--internal variables ----------------------------!
!-------------------------------------------------!

    tail_function_pp = (tailp%c + tailp%b*r + tailp%a*r*r)*(r - tailp%r_out)

  end function tail_function_pp

end module TailFunctions