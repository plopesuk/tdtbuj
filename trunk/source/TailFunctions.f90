!> \brief implements all routines and functions related with tails for hopping
!> integrals and repulsive potentials
!> \author Cristian G. Sanchez and Alin M Elena (Queen's University Belfast)
!> date 1st of November 2007
module m_TailFunctions
  use m_Constants
  use m_Types
  use m_Useful
  implicit none
  private

  public :: makeTailType
  public :: tailFunction
  public :: TailFunctionP
  public :: TailFunctionPp

contains
!> \brief returns  tailType
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
!> \param rIn,rOut reals inner and outer cuttoff
!> \param io type(ioType) i/o units
  function makeTailType(f,fp,fpp,rIn,rOut,io)
  !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'makeTailType'
  !--subroutine parameters -------------------------!
      type(tailType) :: makeTailType
      real(k_pr) :: f, fp, fpp
      real(k_pr) :: rIn, rOut
      real(k_pr) :: x, x2, a, b, c
      type(ioType), intent(inout) :: io
  !--internal variables ----------------------------!
      logical :: status
      character(len=k_ml) :: saux
  !-------------------------------------------------!

      makeTailType%rIn  = rIn

      x = rIn
      x2 = rOut
      a = (10*(12*f - 6*fp*x + fpp*x**2 + 6*fp*x2 - 2*fpp*x*x2 + fpp*x2**2))/(x - x2)**5
      b = (-4*(45*f*x - 21*fp*x**2 + 3*fpp*x**3 + 15*f*x2 + 12*fp*x*x2 - 4*fpp*x**2*x2 + 9*fp*x2**2 - fpp*x*x2**2 + &
          2*fpp*x2**3))/(x - x2)**5
      c =  -((-60*f*x**2 + 24*fp*x**3 - 3*fpp*x**4 - 60*f*x*x2 + 12*fp*x**2*x2 - 36*fp*x*x2**2 + 8*fpp*x**2*x2**2 - &
          4*fpp*x*x2**3 - fpp*x2**4)/(x - x2)**5)

      status = CheckTailParam(f,fp,fpp,x,x2)
      if (.not.status) then
        write(saux,'(a,f0.8,1x,f0.8)') &
      "Second derivative of tail function has solution in ri,rcut region, you should increase rcut", x,x2
        call error(trim(saux),myname,.false.,io)
      endif
      makeTailType%a = a
      makeTailType%b = b
      makeTailType%c = c
      makeTailType%rOut = x2

  end function makeTailType

!> \brief checks the tail parameters,
!> \details to see wether there is an inflexion point in the
!> tail region causing a maximum in the force
!> \author Alin M. Elena,Queen's University Belfast, UK
!> \date 11th of February, 2005
!> \param f,fp,fpp reals the function and its derivatives at inner cuttoff
!> \param ri,rc reals inner and outer cuttoff

  function CheckTailParam(f,fp,fpp,ri,rc)
!--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'CheckTailParam'
!--subroutine parameters -------------------------!
    logical :: CheckTailParam
    real(k_pr) :: f, fp, fpp
    real(k_pr) :: ri, rc
    !--internal variables ----------------------------!
    real(k_pr) :: aux1,aux2,x1,x2,sw,aux3
!-------------------------------------------------!

    aux2 = 2.0_k_pr*(-120.0_k_pr*f - 60.0_k_pr*fp*rc - 10.0_k_pr*fpp*rc**2 + 60.0_k_pr*fp*ri +&
                    20.0_k_pr*fpp*rc*ri - 10.0_k_pr*fpp*ri**2)

    aux1 = (60.0_k_pr*f*rc + 36.0_k_pr*fp*rc**2 + 8.0_k_pr*fpp*rc**3 + 180.0_k_pr*f*ri + 48.0_k_pr*fp*rc*ri - &
        4.0_k_pr*fpp*rc**2*ri - 84.0_k_pr*fp*ri**2 - 16.0_k_pr*fpp*rc*ri**2 + 12.0_k_pr*fpp*ri**3)**2 - &
        4.0_k_pr*(-120.0_k_pr*f - 60.0_k_pr*fp*rc - 10.0_k_pr*fpp*rc**2 + 60.0_k_pr*fp*ri + &
        20.0_k_pr*fpp*rc*ri - 10.0_k_pr*fpp*ri**2)*&
        (-(fpp*rc**4) - 60.0_k_pr*f*rc*ri - 36.0_k_pr*fp*rc**2*ri - 4.0_k_pr*fpp*rc**3*ri - 60.0_k_pr*f*ri**2 +&
    12.0_k_pr*fp*rc*ri**2 + 8.0_k_pr*fpp*rc**2*ri**2 + 24.0_k_pr*fp*ri**3 - 3.0_k_pr*fpp*ri**4)

    aux3 = -60.0_k_pr*f*rc - 36.0_k_pr*fp*rc**2 - 8.0_k_pr*fpp*rc**3 - 180.0_k_pr*f*ri - 48.0_k_pr*fp*rc*ri +&
        4.0_k_pr*fpp*rc**2*ri + 84.0_k_pr*fp*ri**2 + 16.0_k_pr*fpp*rc*ri**2 - 12.0_k_pr*fpp*ri**3

    if (aux1<0.0_k_pr) then
      CheckTailParam=.true.
    else if (abs(aux2)<epsilon(aux2)) then
      CheckTailParam=.true.
    else
      x1 = (aux3 + sqrt(aux1))/aux2
      x2 = (aux3 - sqrt(aux1))/aux2

      if (x1>x2) then
        sw=x1
        x1=x2
        x2=sw
      end if

      if ((ri>x1).and.(rc>x2)) then
        CheckTailParam=.true.
      else if (ri<x2) then
        CheckTailParam=.true.
      else
        if (rc<x1) then
          CheckTailParam=.true.
        else
          CheckTailParam=.false.
        end if
      end if
    end if

  end function CheckTailParam


!>\brief returns the value of the tail
!> function at r, given the parameters in tailp
!>\author Cristian G. Sanchez
!> \date ~2004-2005
!> \param r real the point where the function is evaluated
!> \param tailp the tail parameters
  function tailFunction(r,tailp)
  !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'tailFunction'
!--subroutine parameters -------------------------!
    real(k_pr) :: tailFunction
    type(tailType) :: tailp
    real(k_pr) :: r
!--internal variables ----------------------------!
!-------------------------------------------------!

    tailFunction= ((r - tailp%rOut)**3*(10*tailp%c + &
        5*tailp%b*r + 3*tailp%a*r**2 + 5*tailp%b*tailp%rOut + &
        4*tailp%a*r*tailp%rOut + 3*tailp%a*tailp%rOut**2))/60.0_k_pr

  end function tailFunction


!> \brief returns a the value of the tail function first derivative at r
!> \author Cristian G. Sanchez
!> \date ~2004-2005
!> \param r real the point where the function is evaluated
!> \param tailp the tail parameters
  function TailFunctionP(r,tailp)
  !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'TailFunctionP'
!--subroutine parameters -------------------------!
    real(k_pr) :: TailFunctionP
    type(tailType) :: tailp
    real(k_pr) :: r
!--internal variables ----------------------------!
!-------------------------------------------------!

    TailFunctionP =  ((r - tailp%rOut)**2*(6*tailp%c + &
        4*tailp%b*r + 3*tailp%a*r**2 + 2*tailp%b*tailp%rOut &
        + 2*tailp%a*r*tailp%rOut + tailp%a*tailp%rOut**2))/12.0_k_pr

  end function TailFunctionP

!> \brief returns a the value of the tail function second derivative at r
!> \author Cristian G. Sanchez
!> \date ~2004-2005
!> \param r real the point where the function is evaluated
!> \param tailp the tail parameters
  function TailFunctionPp(r,tailp)
  !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'TailFunctionPp'
!--subroutine parameters -------------------------!
    real(k_pr) :: TailFunctionPp
    type(tailType) :: tailp
    real(k_pr) :: r
!--internal variables ----------------------------!
!-------------------------------------------------!

    TailFunctionPp = (tailp%c + tailp%b*r + tailp%a*r*r)*(r - tailp%rOut)

  end function TailFunctionPp

end module m_TailFunctions