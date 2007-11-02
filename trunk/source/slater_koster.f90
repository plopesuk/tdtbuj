!> \brief Implements Slater Koster table and the first derivatives
!> \details as in the paper "Automatic Generation of Matrix Element Derivatives for Tight Binding Models" \n
!> A. M. Elena and M. Meister, Phys. Rev. B, 72, 165107 (2005)
!> \author Alin M Elena
!> \date 25th of January, 2005

module SlaterKoster
  use constants
  use types
  use useful
  use TightBinding
  implicit none

  private

  public :: hmn
  public :: hmn_p

contains

  function hmn(r,l,m,n,orba,orbb,gen,tb,sol)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'hmn'
    !--subroutine parameters -------------------------!
    real(pr) :: hmn
    real(pr), intent(inout)        :: r,l,m,n
    type(orbital_type), intent(in) :: orba, orbb
    type(model_type), intent(inout) :: tb
    type(solution_type), intent(inout) :: sol
    type(general_type), intent(inout) :: gen

    !--internal variables ----------------------------!  
    real(pr) ::auxs,cg1,sg1,cg2,sg2
    integer :: mp,ls,ss,sb,m1,m2,l1,l2
    real(pr) :: am1,am2,d1,d2,s1,s2,t1,t2,h1

    !-------------------------------------------------!

    if (abs(abs(n)-1.0_pr)>epsilon(n)) then
      cg1 = vcg(abs(orba%m),l,m,n,sol)
      cg2 = vcg(abs(orbb%m),l,m,n,sol)
      sg1 = vsg(abs(orba%m),l,m,n,sol)
      sg2 = vsg(abs(orbb%m),l,m,n,sol)
    else
      cg1=1.0_pr
      cg2=1.0_pr
      sg1=0.0_pr
      sg2=0.0_pr
    endif

    ls=min(orba%l,orbb%l)
    ss=orba%sp
    sb=orbb%sp
    m1=abs(orba%m)
    m2=abs(orbb%m)
    l1=orba%l
    l2=orbb%l
    mp=0

    am1 = am(cg1,sg1,orba%m)
    am2 = am(cg2,sg2,orbb%m)  
    d1 = rotmmpl(m1,0,l1,n,sol)
    d2 = rotmmpl(m2,0,l2,n,sol)
    h1 = rad(r,ss,sb,gen,tb,l1,l2,mp)

    auxs = 2.0_pr*am1*am2*d1*d2*h1
    do mp=1,ls
      s1 = smmpl(orba%m,mp,l1,n,cg1,sg1,sol)
      s2 = smmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
      t1 = tmmpl(orba%m,mp,l1,n,cg1,sg1,sol)
      t2 = tmmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
      auxs = auxs + (s1*s2+t1*t2)*rad(r,ss,sb,gen,tb,l1,l2,mp)
    enddo
    hmn=auxs

  end function hmn

  function hmn_p(r,l,m,n,orba,orbb,alpha,gen,tb,sol,pole)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'sppna_hmn_p'
    !--subroutine parameters -------------------------!
    real(pr) :: hmn_p
    real(pr), intent(inout) :: r,l,m,n
    integer, intent(inout) :: alpha
    type(orbital_type) , intent(in) :: orba, orbb
    type(model_type), intent(inout) :: tb
    type(solution_type), intent(inout) :: sol
      type(general_type), intent(inout) :: gen
    integer, optional :: pole
    !--internal variables ----------------------------!  
    real(pr) :: cg1,cg2,sg1,sg2,auxs
    integer :: mp,ls,ss,sb,m1,m2,l1,l2
    real(pr) :: am1,am2,d1,d2,s1,s2,t1,t2,h1
    real(pr) :: am1a,am2a,d1a,d2a,s1a,s2a,t1a,t2a,h1a
    !-------------------------------------------------!
    ls=min(orba%l,orbb%l)
    ss=orba%sp
    sb=orbb%sp
    m1=abs(orba%m)
    m2=abs(orbb%m)
    l1=orba%l
    l2=orbb%l
    mp=0

    if (present(pole)) then
      select case(pole)
      case(1)
        cg1 = 1.0_pr
        cg2 = 1.0_pr
        sg1 = 0.0_pr
        sg2 = 0.0_pr
      case(2)
        cg1 = 0.0_pr
        cg2 = 0.0_pr
        sg1 = sin(m1*Pi/2.0_pr)
        sg2 = sin(m2*Pi/2.0_pr)
      end select
    else
      cg1 = vcg(m1,l,m,n,sol)
      cg2 = vcg(m2,l,m,n,sol)
      sg1 = vsg(m1,l,m,n,sol)
      sg2 = vsg(m2,l,m,n,sol)
    endif

    select case(alpha)
    case (1)

      am1 = am(cg1,sg1,orba%m)
      am2 = am(cg2,sg2,orbb%m)  
      d1 = rotmmpl(m1,0,l1,n,sol)
      d2 = rotmmpl(m2,0,l2,n,sol)
      h1 = rad(r,ss,sb,gen,tb,l1,l2,mp)

      am1a = am_p(alpha,cg1,sg1,orba%m)
      am2a = am_p(alpha,cg2,sg2,orbb%m)  
      auxs = 2.0_pr*am1a*am2*d1*d2*h1+2.0_pr*am1*am2a*d1*d2*h1
      do mp=1,ls
        s1 = smmpl(orba%m,mp,l1,n,cg1,sg1,sol)
        s2 = smmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
        t1 = tmmpl(orba%m,mp,l1,n,cg1,sg1,sol)
        t2 = tmmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
        s1a = smmpl_p(alpha,orba%m,mp,l1,n,cg1,sg1,sol)
        s2a = smmpl_p(alpha,orbb%m,mp,l2,n,cg2,sg2,sol)
        t1a = tmmpl_p(alpha,orba%m,mp,l1,n,cg1,sg1,sol)
        t2a = tmmpl_p(alpha,orbb%m,mp,l2,n,cg2,sg2,sol)

        auxs = auxs + (s1a*s2+s1*s2a+t1a*t2+t1*t2a)*rad(r,ss,sb,gen,tb,l1,l2,mp)
      enddo
      hmn_p = auxs

    case (2)
      am1 = am(cg1,sg1,orba%m)
      am2 = am(cg2,sg2,orbb%m)  
      d1 = rotmmpl(m1,0,l1,n,sol)
      d2 = rotmmpl(m2,0,l2,n,sol)
      h1 = rad(r,ss,sb,gen,tb,l1,l2,mp)
      d1a = rotmmpl_p(alpha,m1,0,l1,n,sol)
      d2a = rotmmpl_p(alpha,m2,0,l2,n,sol)

      auxs = 2.0_pr*am1*am2*d1a*d2*h1+2.0_pr*am1*am2*d1*d2a*h1
      do mp=1,ls
        s1 = smmpl(orba%m,mp,l1,n,cg1,sg1,sol)
        s2 = smmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
        t1 = tmmpl(orba%m,mp,l1,n,cg1,sg1,sol)
        t2 = tmmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
        s1a = smmpl_p(alpha,orba%m,mp,l1,n,cg1,sg1,sol)
        s2a = smmpl_p(alpha,orbb%m,mp,l2,n,cg2,sg2,sol)
        t1a = tmmpl_p(alpha,orba%m,mp,l1,n,cg1,sg1,sol)
        t2a = tmmpl_p(alpha,orbb%m,mp,l2,n,cg2,sg2,sol)
        auxs = auxs + (s1a*s2+s1*s2a+t1a*t2+t1*t2a)*rad(r,ss,sb,gen,tb,l1,l2,mp)
      enddo
      hmn_p = auxs

    case (3)
      am1 = am(cg1,sg1,orba%m)
      am2 = am(cg2,sg2,orbb%m)  
      d1 = rotmmpl(m1,0,l1,n,sol)
      d2 = rotmmpl(m2,0,l2,n,sol)
      h1a = rad_p(alpha,r,ss,sb,gen,tb,l1,l2,mp)
      auxs = 2.0_pr*am1*am2*d1*d2*h1a

      do mp=1,ls
        s1 = smmpl(orba%m,mp,l1,n,cg1,sg1,sol)
        s2 = smmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
        t1 = tmmpl(orba%m,mp,l1,n,cg1,sg1,sol)
        t2 = tmmpl(orbb%m,mp,l2,n,cg2,sg2,sol)
        s1a = smmpl_p(alpha,orba%m,mp,l1,n,cg1,sg1,sol)
        s2a = smmpl_p(alpha,orbb%m,mp,l2,n,cg2,sg2,sol)
        t1a = tmmpl_p(alpha,orba%m,mp,l1,n,cg1,sg1,sol)
        t2a = tmmpl_p(alpha,orbb%m,mp,l2,n,cg2,sg2,sol)
        auxs = auxs + (s1*s2+t1*t2)*rad_p(alpha,r,ss,sb,gen,tb,l1,l2,mp)
      enddo
      hmn_p = auxs
    end select
  end function hmn_p

  function vsg(mp,l,m,n,sol)
       !----------------------------------------------------!
       ! vsg gives you sin(mx) knowing sin(x) and cos(x)    !
       ! input an integer mp and n,l,m direction cosines    !
       !                                                    ! 
       ! 18th of January, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!
    integer, intent(in) :: mp
    real(pr), intent(in) :: n,l,m
    type(solution_type), intent(in) :: sol
    real(pr) :: vsg
    real(pr) :: summ!,factm
    integer ::i


    summ=0.0_pr
    do i=0,int((mp-1)/2.0_pr)
      if (mp-2*i-1>=0)  summ=summ+(-1.0_pr)**i*sol%fact(mp)/(sol%fact(2*i+1)*sol%fact(mp-2*i-1))*&
        (m/sqrt(1.0_pr-n*n))**(2*i+1)*(l/sqrt(1.0_pr-n*n))**(mp-2*i-1)
    enddo
    vsg=summ

  end function vsg


  function vcg(mp,l,m,n,sol)
       !----------------------------------------------------!
       ! vsg gives you sin(mx) knowing sin(x),cos(x)        !
       ! input an integer mp and n,l,m direction cosines    !
       !                                                    ! 
       ! 18th of January, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!

    integer, intent(in) :: mp
    real(pr), intent(in) :: n,l,m
    type(solution_type), intent(in) :: sol
    real(pr) :: vcg
    real(pr) :: summ!,factm
    integer ::i



    summ=0.0_pr
    do i=0,int(mp/2)
      summ=summ+(-1.0_pr)**i*sol%fact(mp)/(sol%fact(2*i)*sol%fact(mp-2*i))*&
          (m/sqrt(1.0_pr-n*n))**(2*i)*(l/sqrt(1.0_pr-n*n))**(mp-2*i)
    enddo
    vcg=summ
  end function vcg



  function am(cg,sg,m)
       !----------------------------------------------------!
       ! am implements the equation (TBP)                   !
       ! input an integer m, cos(mx), sin(mx)               !
       !                                                    !
       ! 4th of February, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!
    integer, intent(in) :: m
    real(pr), intent(in) :: cg,sg
    real(pr) :: am
    logical :: l=.true.

    if (m/=0) then
      am = (-1.0_pr)**abs(m)*(cg*h(m,l)+sg*h(-m,l))
    else
      am = 1.0_pr/sqrt(2.0_pr)
    endif

  end function am

  function am_p(alpha,cg,sg,m)
       !----------------------------------------------------!
       ! am_p implements the derivatives of am with respect !
       ! polar coordinates                                  ! 
       ! the same as am plus an integer alpha meaning       !
       ! variable of derivative                             !
       !                                                    !
       ! 4th of February, 2005                              !
       ! Alin M. Elena, Queen's University Belfast, UK      !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!
    integer, intent(in) :: alpha,m
    real(pr), intent(in) ::cg,sg
    real(pr) :: am_p

    
    if (alpha == 1) then
      am_p=abs(m)*bm(cg,sg,m)
    else
      am_p=0.0_pr
    endif

  end function am_p

  function bm(cg,sg,m)
       !----------------------------------------------------!
       ! bm implements the equation (TBP)                   !
       ! input an integer m, cos(mx), sin(mx)               !
       !                                                    !
       ! 4th of February, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!

    integer, intent(in) :: m
    real(pr), intent(in) :: cg,sg
    real(pr) :: bm
    logical :: l=.true.
      
    if (m/=0) then
      bm = (-1.0_pr)**abs(m)*(cg*h(-m,l)-sg*h(m,l))
    else
      bm=0.0_pr
    endif

  end function bm

  function bm_p(alpha,cg,sg,m)
       !----------------------------------------------------!
       ! bm_p implements the derivatives of am with respect !
       ! polar coordinates                                  ! 
       ! the same as am plus an integer alpha meaning       !
       ! variable of derivative                             !
       !                                                    !
       ! 4th of February, 2005                              !
       ! Alin M. Elena, Queen's University Belfast, UK      !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!

    integer, intent(in) :: alpha,m
    real(pr), intent(in) ::cg,sg
    real(pr) :: bm_p

    if (alpha == 1) then
      bm_p=-abs(m)*am(cg,sg,m)
    else
      bm_p = 0.0_pr
    endif
  end function bm_p

  function rotmmpl(m,mp,l,n,sol)
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'rotmmpl'
    !--subroutine parameters -------------------------!
       !----------------------------------------------------!
       ! rotmmpl implements the equation (TBP)              !
       ! input  integers m,mp,l and dp n = cos(beta)        !
       !                                                    !
       ! 4th of February, 2005                              !
       ! Alin M. Elena,Queen's University  Belfast, UK      !
       ! Modified by Alin M Elena                                       !
       ! Date 31st of March 2006                                              !
       !----------------------------------------------------!
    integer, intent(in) :: m,mp,l
    real(pr), intent(in) :: n
    type(solution_type), intent(in) :: sol
    real(pr) ::rotmmpl

    real(pr) ::sum
    integer ::t
    sum=0.0_pr
    do t=max(0,-m-mp),min(l-m,l-mp)
      sum=sum+((1+N)**(t+(m+mp)/2.0_pr))*((1-N)**(l-t-(m+mp)/2.0_pr))*((-1)**t)/&
        (sol%fact(l-mp-t)*sol%fact(l-m-t)*sol%fact(t)*sol%fact(t+m+mp))
    enddo
    rotmmpl=(2.0_pr**(-l))*sqrt(sol%fact(l+mp)*sol%fact(l-mp)*sol%fact(l+m)*sol%fact(l-m))*sum*(-1.0_pr)**(l-mp)

  end function rotmmpl


  function rotmmpl_p(alpha,m,mp,l,n,sol)
       !----------------------------------------------------!
       ! rotmmpl_p implements the derivative of rotmmpl     !
       ! with respect to polar coordinate alpha             !
       ! input  integers m,mp,l and dp n = cos(beta), alpha !
       !                                                    !
       ! 4th of February, 2005                              !
       ! Alin M. Elena,Queen's University  Belfast, UK      !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!

    integer, intent(in) :: m,mp,l,alpha
    real(pr), intent(in) :: n
    type(solution_type), intent(in) :: sol
    real(pr) ::rotmmpl_p

    real(pr) ::t1,t2
!       integer ::t
    select case(alpha)
    case (1)
      rotmmpl_p=0.0_pr
    case (2)
      t1=0.0_pr
      t2=0.0_pr
      if (abs(mp-1)<=l) &
        t1=sqrt(real((l+mp)*(l-mp+1),pr))*rotmmpl(m,mp-1,l,n,sol)
      if (abs(mp+1)<=l) &
        t2=-sqrt(real((l-mp)*(l+mp+1),pr))*rotmmpl(m,mp+1,l,n,sol)

      rotmmpl_p=0.5_pr*(t1+t2)
    case (3)
      rotmmpl_p=0.0_pr
    end select
  end function rotmmpl_p


  function tmmpl(m,mp,l,n,cg,sg,sol)
       !----------------------------------------------------!
       ! tmmpl implements matrix T the equation (TBP)       !
       ! inputing an integer m,mp,l and dp n=cos(b)         !
       ! sin(mx)=sg,cos(mx)=cg                              !
       !                                                    !
       ! 18th of January, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!
    integer, intent(in) :: m,mp,l
    real(pr), intent(in) :: n,cg,sg
    type(solution_type), intent(inout) :: sol
    real(pr) ::tmmpl

    integer ::absm!, absmp
    if (m /= 0) then
      absm=abs(m)
      tmmpl=bm(cg,sg,m)*(((-1)**mp)*rotmmpl(absm,mp,l,n,sol)-rotmmpl(absm,-mp,l,n,sol))
    else
      tmmpl = 0.0_pr
    endif
  end function tmmpl


  function tmmpl_p(alpha,m,mp,l,n,cg,sg,sol)
       !----------------------------------------------------!
       ! tmmpl_p implements derivative of matrix T          !
       ! equation (TBP) with respect to polar coordinates   !
       ! alpha                                              !
       ! inputing an integer m,mp,l and dp n=cos(b)         !
       ! sin(mx)=sg,cos(mx)=cg                              !
       !                                                    !
       ! 18th of January, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!

    integer, intent(in) :: m,mp,l,alpha
    real(pr), intent(in) :: n,cg,sg
    type(solution_type), intent(inout) :: sol
    real(pr) ::tmmpl_p

    integer ::absm


    absm=abs(m)
    select case (alpha)
    case (1)
      if (m==0) then
        tmmpl_p=0.0_pr
      else   
        tmmpl_p=-absm*am(cg,sg,m)*(((-1)**mp)*rotmmpl(absm,mp,l,n,sol)&
          -rotmmpl(absm,-mp,l,n,sol))
      endif
    case(2)
      if (m==0) then
        tmmpl_p=0.0_pr
      else 
        tmmpl_p=bm(cg,sg,m)*(((-1)**mp)*rotmmpl_p(alpha,absm,mp,l,n,sol)&
          -rotmmpl_p(alpha,absm,-mp,l,n,sol))
      endif
    case(3)
      tmmpl_p=0.0_pr
    end select

  end function tmmpl_p

  function smmpl(m,mp,l,n,cg,sg,sol)
       !----------------------------------------------------!
       ! smmpl implements matrix S the equation (TBP)       !
       ! inputing an integer m,mp,l and dp n=cos(b)         !
       ! sin(mx)=sg,cos(mx)=cg                              !
       !                                                    !
       ! 18th of January, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!

    integer, intent(in) :: m,mp,l
    real(pr), intent(in) :: n,cg,sg
    type(solution_type), intent(inout) :: sol
    real(pr) ::smmpl

    integer ::absm

    absm=abs(m)
    smmpl=am(cg,sg,m)*(((-1)**mp)*rotmmpl(absm,mp,l,n,sol)+rotmmpl(absm,-mp,l,n,sol))

  end function smmpl

  function smmpl_p(alpha,m,mp,l,n,cg,sg,sol)
       !----------------------------------------------------!
       ! smmpl_p implements derivative of matrix S          !
       ! equation (TBP) with respect to polar coordinates   !
       ! alpha                                              !
       ! inputing an integer m,mp,l and dp n=cos(b)         !
       ! sin(mx)=sg,cos(mx)=cg                              !
       !                                                    !
       ! 25th of January, 2005                              !
       ! Alin M. Elena,Queen's University Belfast, UK       !
       ! Modified by                                        !
       ! Date                                               !
       !----------------------------------------------------!


    integer, intent(in) :: m,mp,l,alpha
    real(pr), intent(in) :: n,sg,cg
    type(solution_type), intent(inout) :: sol
    real(pr) ::smmpl_p
    integer ::absm

    absm=abs(m)
    select case(alpha)
    case (1)
      smmpl_p=absm*bm(cg,sg,m)*(((-1)**mp)*rotmmpl(absm,mp,l,n,sol)&
        +rotmmpl(absm,-mp,l,n,sol))
    case (2)
      smmpl_p=am(cg,sg,m)*(((-1)**mp)*rotmmpl_p(alpha,absm,mp,l,n,sol)&
        +rotmmpl_p(alpha,absm,-mp,l,n,sol))
    case(3)
      smmpl_p=0.0_pr
    end select
  end function smmpl_p

end module SlaterKoster