!> \brief Implements Slater Koster table and the first derivatives
!> \details as in the paper "Automatic Generation of Matrix Element Derivatives for Tight Binding Models" \n
!> A. M. Elena and M. Meister, Phys. Rev. B, 72, 165107 (2005)
!> \author Alin M Elena
!> \date 25th of January, 2005
!
module m_SlaterKoster
  use m_Constants
  use m_Types
  use m_Useful
  use m_TightBinding
  implicit none
!
  private
!
  public :: hmn
  public :: DhmnXYZ
!
contains
!
  function hmn (r, l, m, n, orba, orbb, gen, tb, sol)
!--subroutine name--------------------------------!    
    character (len=*), parameter :: myname = 'hmn'
!--subroutine parameters -------------------------!    
    real (k_pr) :: hmn
    real (k_pr), intent (inout) :: r, l, m, n
    type (orbitalType), intent (in) :: orba, orbb
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (generalType), intent (inout) :: gen
!
!--internal variables ----------------------------!    
    real (k_pr) :: auxs, cg1, sg1, cg2, sg2
    integer :: mp, ls, ss, sb, m1, m2, l1, l2
    real (k_pr) :: am1, am2, d1, d2, s1, s2, t1, t2, h1
!
!-------------------------------------------------!    
!
    if (Abs(Abs(n)-1.0_k_pr) > epsilon(n)) then
      cg1 = vcg (Abs(orba%m), l, m, n, sol)
      cg2 = vcg (Abs(orbb%m), l, m, n, sol)
      sg1 = vsg (Abs(orba%m), l, m, n, sol)
      sg2 = vsg (Abs(orbb%m), l, m, n, sol)
    else
      cg1 = 1.0_k_pr
      cg2 = 1.0_k_pr
      sg1 = 0.0_k_pr
      sg2 = 0.0_k_pr
    end if
!
    ls = Min (orba%l, orbb%l)
    ss = orba%sp
    sb = orbb%sp
    m1 = Abs (orba%m)
    m2 = Abs (orbb%m)
    l1 = orba%l
    l2 = orbb%l
    mp = 0
!
    am1 = am (cg1, sg1, orba%m)
    am2 = am (cg2, sg2, orbb%m)
    d1 = rotmmpl (m1, 0, l1, n, sol)
    d2 = rotmmpl (m2, 0, l2, n, sol)
    h1 = rad (r, ss, sb, gen, tb, l1, l2, mp)
!
    auxs = 2.0_k_pr * am1 * am2 * d1 * d2 * h1
    do mp = 1, ls
      s1 = smmpl (orba%m, mp, l1, n, cg1, sg1, sol)
      s2 = smmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
      t1 = tmmpl (orba%m, mp, l1, n, cg1, sg1, sol)
      t2 = tmmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
      auxs = auxs + (s1*s2+t1*t2) * rad (r, ss, sb, gen, tb, l1, l2, mp)
    end do
    hmn = auxs
!
  end function hmn
!
  function HmnP (r, l, m, n, orba, orbb, alpha, gen, tb, sol, pole)
!--subroutine name--------------------------------!    
    character (len=*), parameter :: myname = 'sppna_HmnP'
!--subroutine parameters -------------------------!    
    real (k_pr) :: HmnP
    real (k_pr), intent (inout) :: r, l, m, n
    integer, intent (inout) :: alpha
    type (orbitalType), intent (in) :: orba, orbb
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (generalType), intent (inout) :: gen
    integer, optional :: pole
!--internal variables ----------------------------!    
    real (k_pr) :: cg1, cg2, sg1, sg2, auxs
    integer :: mp, ls, ss, sb, m1, m2, l1, l2
    real (k_pr) :: am1, am2, d1, d2, s1, s2, t1, t2, h1
    real (k_pr) :: am1a, am2a, d1a, d2a, s1a, s2a, t1a, t2a, h1a
!-------------------------------------------------!    
    ls = Min (orba%l, orbb%l)
    ss = orba%sp
    sb = orbb%sp
    m1 = Abs (orba%m)
    m2 = Abs (orbb%m)
    l1 = orba%l
    l2 = orbb%l
    mp = 0
!
    if (present(pole)) then
      select case (pole)
      case (1)
        cg1 = 1.0_k_pr
        cg2 = 1.0_k_pr
        sg1 = 0.0_k_pr
        sg2 = 0.0_k_pr
      case (2)
        cg1 = 0.0_k_pr
        cg2 = 0.0_k_pr
        sg1 = Sin (m1*k_pi/2.0_k_pr)
        sg2 = Sin (m2*k_pi/2.0_k_pr)
      end select
    else
      cg1 = vcg (m1, l, m, n, sol)
      cg2 = vcg (m2, l, m, n, sol)
      sg1 = vsg (m1, l, m, n, sol)
      sg2 = vsg (m2, l, m, n, sol)
    end if
!
    select case (alpha)
    case (1)
!
      am1 = am (cg1, sg1, orba%m)
      am2 = am (cg2, sg2, orbb%m)
      d1 = rotmmpl (m1, 0, l1, n, sol)
      d2 = rotmmpl (m2, 0, l2, n, sol)
      h1 = rad (r, ss, sb, gen, tb, l1, l2, mp)
!
      am1a = AmP (alpha, cg1, sg1, orba%m)
      am2a = AmP (alpha, cg2, sg2, orbb%m)
      auxs = 2.0_k_pr * am1a * am2 * d1 * d2 * h1 + 2.0_k_pr * am1 * am2a * d1 * d2 * h1
      do mp = 1, ls
        s1 = smmpl (orba%m, mp, l1, n, cg1, sg1, sol)
        s2 = smmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
        t1 = tmmpl (orba%m, mp, l1, n, cg1, sg1, sol)
        t2 = tmmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
        s1a = SmmplP (alpha, orba%m, mp, l1, n, cg1, sg1, sol)
        s2a = SmmplP (alpha, orbb%m, mp, l2, n, cg2, sg2, sol)
        t1a = TmmplP (alpha, orba%m, mp, l1, n, cg1, sg1, sol)
        t2a = TmmplP (alpha, orbb%m, mp, l2, n, cg2, sg2, sol)
!
        auxs = auxs + (s1a*s2+s1*s2a+t1a*t2+t1*t2a) * rad (r, ss, sb, gen, tb, l1, l2, mp)
      end do
      HmnP = auxs
!
    case (2)
      am1 = am (cg1, sg1, orba%m)
      am2 = am (cg2, sg2, orbb%m)
      d1 = rotmmpl (m1, 0, l1, n, sol)
      d2 = rotmmpl (m2, 0, l2, n, sol)
      h1 = rad (r, ss, sb, gen, tb, l1, l2, mp)
      d1a = RotmmplP (alpha, m1, 0, l1, n, sol)
      d2a = RotmmplP (alpha, m2, 0, l2, n, sol)
!
      auxs = 2.0_k_pr * am1 * am2 * d1a * d2 * h1 + 2.0_k_pr * am1 * am2 * d1 * d2a * h1
      do mp = 1, ls
        s1 = smmpl (orba%m, mp, l1, n, cg1, sg1, sol)
        s2 = smmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
        t1 = tmmpl (orba%m, mp, l1, n, cg1, sg1, sol)
        t2 = tmmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
        s1a = SmmplP (alpha, orba%m, mp, l1, n, cg1, sg1, sol)
        s2a = SmmplP (alpha, orbb%m, mp, l2, n, cg2, sg2, sol)
        t1a = TmmplP (alpha, orba%m, mp, l1, n, cg1, sg1, sol)
        t2a = TmmplP (alpha, orbb%m, mp, l2, n, cg2, sg2, sol)
        auxs = auxs + (s1a*s2+s1*s2a+t1a*t2+t1*t2a) * rad (r, ss, sb, gen, tb, l1, l2, mp)
      end do
      HmnP = auxs
!
    case (3)
      am1 = am (cg1, sg1, orba%m)
      am2 = am (cg2, sg2, orbb%m)
      d1 = rotmmpl (m1, 0, l1, n, sol)
      d2 = rotmmpl (m2, 0, l2, n, sol)
      h1a = RadP (alpha, r, ss, sb, gen, tb, l1, l2, mp)
      auxs = 2.0_k_pr * am1 * am2 * d1 * d2 * h1a
!
      do mp = 1, ls
        s1 = smmpl (orba%m, mp, l1, n, cg1, sg1, sol)
        s2 = smmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
        t1 = tmmpl (orba%m, mp, l1, n, cg1, sg1, sol)
        t2 = tmmpl (orbb%m, mp, l2, n, cg2, sg2, sol)
        s1a = SmmplP (alpha, orba%m, mp, l1, n, cg1, sg1, sol)
        s2a = SmmplP (alpha, orbb%m, mp, l2, n, cg2, sg2, sol)
        t1a = TmmplP (alpha, orba%m, mp, l1, n, cg1, sg1, sol)
        t2a = TmmplP (alpha, orbb%m, mp, l2, n, cg2, sg2, sol)
        auxs = auxs + (s1*s2+t1*t2) * RadP (alpha, r, ss, sb, gen, tb, l1, l2, mp)
      end do
      HmnP = auxs
    end select
  end function HmnP
!
  function vsg (mp, l, m, n, sol)
!----------------------------------------------------!       
! vsg gives you sin(mx) knowing sin(x) and cos(x)    !       
! input an integer mp and n,l,m direction cosines    !       
!                                                    !       
! 18th of January, 2005                              !       
! Alin M. Elena,Queen's University Belfast, UK       !       
! Modified by                                        !       
! Date                                               !       
!----------------------------------------------------!       
    integer, intent (in) :: mp
    real (k_pr), intent (in) :: n, l, m
    type (solutionType), intent (in) :: sol
    real (k_pr) :: vsg
    real (k_pr) :: summ !,factm
    integer :: i
!
!
    summ = 0.0_k_pr
    do i = 0, Int ((mp-1)/2.0_k_pr)
      if (mp-2*i-1 >= 0) summ = summ + (-1.0_k_pr) ** i * sol%fact(mp) / (sol%fact(2*i+1)*sol%fact(mp-2*i-1)) * &
     & (m/Sqrt(1.0_k_pr-n*n)) ** (2*i+1) * (l/Sqrt(1.0_k_pr-n*n)) ** (mp-2*i-1)
    end do
    vsg = summ
!
  end function vsg
!
!
  function vcg (mp, l, m, n, sol)
!----------------------------------------------------!       
! vsg gives you sin(mx) knowing sin(x),cos(x)        !       
! input an integer mp and n,l,m direction cosines    !       
!                                                    !       
! 18th of January, 2005                              !       
! Alin M. Elena,Queen's University Belfast, UK       !       
! Modified by                                        !       
! Date                                               !       
!----------------------------------------------------!       
!
    integer, intent (in) :: mp
    real (k_pr), intent (in) :: n, l, m
    type (solutionType), intent (in) :: sol
    real (k_pr) :: vcg
    real (k_pr) :: summ !,factm
    integer :: i
!
!
!
    summ = 0.0_k_pr
    do i = 0, Int (mp/2)
      summ = summ + (-1.0_k_pr) ** i * sol%fact(mp) / (sol%fact(2*i)*sol%fact(mp-2*i)) * (m/Sqrt(1.0_k_pr-n*n)) ** (2*i) * &
     & (l/Sqrt(1.0_k_pr-n*n)) ** (mp-2*i)
    end do
    vcg = summ
  end function vcg
!
!
!
  function am (cg, sg, m)
!----------------------------------------------------!       
! am implements the equation (TBP)                   !       
! input an integer m, cos(mx), sin(mx)               !       
!                                                    !       
! 4th of February, 2005                              !       
! Alin M. Elena,Queen's University Belfast, UK       !       
! Modified by                                        !       
! Date                                               !       
!----------------------------------------------------!       
    integer, intent (in) :: m
    real (k_pr), intent (in) :: cg, sg
    real (k_pr) :: am
    logical :: l = .true.
!
    if (m /= 0) then
      am = (-1.0_k_pr) ** Abs (m) * (cg*h(m, l)+sg*h(-m, l))
    else
      am = 1.0_k_pr / Sqrt (2.0_k_pr)
    end if
!
  end function am
!
  function AmP (alpha, cg, sg, m)
!----------------------------------------------------!       
! AmP implements the derivatives of am with respect !       
! polar coordinates                                  !       
! the same as am plus an integer alpha meaning       !       
! variable of derivative                             !       
!                                                    !       
! 4th of February, 2005                              !       
! Alin M. Elena, Queen's University Belfast, UK      !       
! Modified by                                        !       
! Date                                               !       
!----------------------------------------------------!       
    integer, intent (in) :: alpha, m
    real (k_pr), intent (in) :: cg, sg
    real (k_pr) :: AmP
!
!
    if (alpha == 1) then
      AmP = Abs (m) * bm (cg, sg, m)
    else
      AmP = 0.0_k_pr
    end if
!
  end function AmP
!
  function bm (cg, sg, m)
!----------------------------------------------------!       
! bm implements the equation (TBP)                   !       
! input an integer m, cos(mx), sin(mx)               !       
!                                                    !       
! 4th of February, 2005                              !       
! Alin M. Elena,Queen's University Belfast, UK       !       
! Modified by                                        !       
! Date                                               !       
!----------------------------------------------------!       
!
    integer, intent (in) :: m
    real (k_pr), intent (in) :: cg, sg
    real (k_pr) :: bm
    logical :: l = .true.
!
    if (m /= 0) then
      bm = (-1.0_k_pr) ** Abs (m) * (cg*h(-m, l)-sg*h(m, l))
    else
      bm = 0.0_k_pr
    end if
!
  end function bm
!
  function BmP (alpha, cg, sg, m)
!----------------------------------------------------!       
! BmP implements the derivatives of am with respect !       
! polar coordinates                                  !       
! the same as am plus an integer alpha meaning       !       
! variable of derivative                             !       
!                                                    !       
! 4th of February, 2005                              !       
! Alin M. Elena, Queen's University Belfast, UK      !       
! Modified by                                        !       
! Date                                               !       
!----------------------------------------------------!       
!
    integer, intent (in) :: alpha, m
    real (k_pr), intent (in) :: cg, sg
    real (k_pr) :: BmP
!
    if (alpha == 1) then
      BmP = - Abs (m) * am (cg, sg, m)
    else
      BmP = 0.0_k_pr
    end if
  end function BmP
!
  function rotmmpl (m, mp, l, n, sol)
!--subroutine name--------------------------------!    
    character (len=*), parameter :: myname = 'rotmmpl'
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
    integer, intent (in) :: m, mp, l
    real (k_pr), intent (in) :: n
    type (solutionType), intent (in) :: sol
    real (k_pr) :: rotmmpl
!
    real (k_pr) :: sum
    integer :: t
    sum = 0.0_k_pr
    do t = Max (0,-m-mp), Min (l-m, l-mp)
      sum = sum + ((1+n)**(t+(m+mp)/2.0_k_pr)) * ((1-n)**(l-t-(m+mp)/2.0_k_pr)) * ((-1)**t) / &
     & (sol%fact(l-mp-t)*sol%fact(l-m-t)*sol%fact(t)*sol%fact(t+m+mp))
    end do
    rotmmpl = (2.0_k_pr**(-l)) * Sqrt (sol%fact(l+mp)*sol%fact(l-mp)*sol%fact(l+m)*sol%fact(l-m)) * sum * (-1.0_k_pr) ** (l-mp)
!
  end function rotmmpl
!
!
  function RotmmplP (alpha, m, mp, l, n, sol)
!----------------------------------------------------!       
! RotmmplP implements the derivative of rotmmpl     !       
! with respect to polar coordinate alpha             !       
! input  integers m,mp,l and dp n = cos(beta), alpha !       
!                                                    !       
! 4th of February, 2005                              !       
! Alin M. Elena,Queen's University  Belfast, UK      !       
! Modified by                                        !       
! Date                                               !       
!----------------------------------------------------!       
!
    integer, intent (in) :: m, mp, l, alpha
    real (k_pr), intent (in) :: n
    type (solutionType), intent (in) :: sol
    real (k_pr) :: RotmmplP
!
    real (k_pr) :: t1, t2
!       integer ::t
    select case (alpha)
    case (1)
      RotmmplP = 0.0_k_pr
    case (2)
      t1 = 0.0_k_pr
      t2 = 0.0_k_pr
      if (Abs(mp-1) <= l) t1 = Sqrt (real((l+mp)*(l-mp+1), k_pr)) * rotmmpl (m, mp-1, l, n, sol)
      if (Abs(mp+1) <= l) t2 = - Sqrt (real((l-mp)*(l+mp+1), k_pr)) * rotmmpl (m, mp+1, l, n, sol)
!
      RotmmplP = 0.5_k_pr * (t1+t2)
    case (3)
      RotmmplP = 0.0_k_pr
    end select
  end function RotmmplP
!
!
  function tmmpl (m, mp, l, n, cg, sg, sol)
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
    integer, intent (in) :: m, mp, l
    real (k_pr), intent (in) :: n, cg, sg
    type (solutionType), intent (inout) :: sol
    real (k_pr) :: tmmpl
!
    integer :: absm !, absmp
    if (m /= 0) then
      absm = Abs (m)
      tmmpl = bm (cg, sg, m) * (((-1)**mp)*rotmmpl(absm, mp, l, n, sol)-rotmmpl(absm,-mp, l, n, sol))
    else
      tmmpl = 0.0_k_pr
    end if
  end function tmmpl
!
!
  function TmmplP (alpha, m, mp, l, n, cg, sg, sol)
!----------------------------------------------------!       
! TmmplP implements derivative of matrix T          !       
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
!
    integer, intent (in) :: m, mp, l, alpha
    real (k_pr), intent (in) :: n, cg, sg
    type (solutionType), intent (inout) :: sol
    real (k_pr) :: TmmplP
!
    integer :: absm
!
!
    absm = Abs (m)
    select case (alpha)
    case (1)
      if (m == 0) then
        TmmplP = 0.0_k_pr
      else
        TmmplP = - absm * am (cg, sg, m) * (((-1)**mp)*rotmmpl(absm, mp, l, n, sol)-rotmmpl(absm,-mp, l, n, sol))
      end if
    case (2)
      if (m == 0) then
        TmmplP = 0.0_k_pr
      else
        TmmplP = bm (cg, sg, m) * (((-1)**mp)*RotmmplP(alpha, absm, mp, l, n, sol)-RotmmplP(alpha, absm,-mp, l, n, sol))
      end if
    case (3)
      TmmplP = 0.0_k_pr
    end select
!
  end function TmmplP
!
  function smmpl (m, mp, l, n, cg, sg, sol)
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
!
    integer, intent (in) :: m, mp, l
    real (k_pr), intent (in) :: n, cg, sg
    type (solutionType), intent (inout) :: sol
    real (k_pr) :: smmpl
!
    integer :: absm
!
    absm = Abs (m)
    smmpl = am (cg, sg, m) * (((-1)**mp)*rotmmpl(absm, mp, l, n, sol)+rotmmpl(absm,-mp, l, n, sol))
!
  end function smmpl
!
  function SmmplP (alpha, m, mp, l, n, cg, sg, sol)
!----------------------------------------------------!       
! SmmplP implements derivative of matrix S          !       
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
!
!
    integer, intent (in) :: m, mp, l, alpha
    real (k_pr), intent (in) :: n, sg, cg
    type (solutionType), intent (inout) :: sol
    real (k_pr) :: SmmplP
    integer :: absm
!
    absm = Abs (m)
    select case (alpha)
    case (1)
      SmmplP = absm * bm (cg, sg, m) * (((-1)**mp)*rotmmpl(absm, mp, l, n, sol)+rotmmpl(absm,-mp, l, n, sol))
    case (2)
      SmmplP = am (cg, sg, m) * (((-1)**mp)*RotmmplP(alpha, absm, mp, l, n, sol)+RotmmplP(alpha, absm,-mp, l, n, sol))
    case (3)
      SmmplP = 0.0_k_pr
    end select
  end function SmmplP
!
!==hoppings==========================================================
  real (k_pr) function DhmnXYZ (alpha, r, l, m, n, orba, orbb, gen, tb, sol)
!--subroutine name--------------------------------!    
    character (len=*), parameter :: myname = 'DhmnXYZ'
!--subroutine parameters -------------------------!    
    real (k_pr) :: r, l, m, n
    type (orbitalType) :: orba, orbb
    integer :: alpha
    real (k_pr) :: da, db, dr
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
    type (solutionType), intent (inout) :: sol
    integer :: one, two, three
    one = 1
    two = 2
    three = 3
    if (Abs(Abs(n)-1.0_k_pr) > epsilon(n)) then
      da = HmnP (r, l, m, n, orba, orbb, one, gen, tb, sol)
      db = HmnP (r, l, m, n, orba, orbb, two, gen, tb, sol)
      dr = HmnP (r, l, m, n, orba, orbb, three, gen, tb, sol)
!
!
      DhmnXYZ = da * jx (alpha, r, l, m, n) - db * jy (alpha, r, l, m, n) + dr * jz (alpha, l, m, n)
    else
! the north-south pole case   
      select case (alpha)
      case (1)
        DhmnXYZ = sign (1.0_k_pr, n) * HmnP (r, l, m, n, orba, orbb, two, gen, tb, sol, 1) / r
      case (2)
        DhmnXYZ = sign (1.0_k_pr, n) * HmnP (r, l, m, n, orba, orbb, two, gen, tb, sol, 2) / r
      case (3)
        DhmnXYZ = sign (1.0_k_pr, n) * HmnP (r, l, m, n, orba, orbb, three, gen, tb, sol, 1)
      end select
!
    end if
!      write(*,'(a,f8.4,4f24.16)')"DhmnXYZ ",gen%currSimTime,l,m,n,DhmnXYZ
  end function DhmnXYZ
!
  real (k_pr) function jx (alpha, r, l, m, n)
    character (len=*), parameter :: myname = 'jx'
    real (k_pr) :: r, l, m, n
    integer :: alpha
!
    select case (alpha)
    case (1)
      jx = - m / (r*(1.0_k_pr-n*n))
    case (2)
      jx = l / (r*(1.0_k_pr-n*n))
    case (3)
      jx = 0.0_k_pr
    end select
  end function jx
!
  real (k_pr) function jy (alpha, r, l, m, n)
!--subroutine name--------------------------------!    
    character (len=*), parameter :: myname = 'jy'
!--subroutine parameters -------------------------!    
    real (k_pr) :: r, l, m, n
    integer :: alpha
!
    select case (alpha)
    case (1)
      jy = - n * l / (r*Sqrt(1.0_k_pr-n*n))
    case (2)
      jy = - n * m / (r*Sqrt(1.0_k_pr-n*n))
    case (3)
      jy = Sqrt (1.0_k_pr-n*n) / r
    end select
!
!
  end function jy
!
!
  real (k_pr) function jz (alpha, l, m, n)
!--subroutine name--------------------------------!    
    character (len=*), parameter :: myname = 'jz'
!--subroutine parameters -------------------------!    
    real (k_pr) :: l, m, n
    integer :: alpha
!
    select case (alpha)
    case (1)
      jz = l
    case (2)
      jz = m
    case (3)
      jz = n
    end select
!
  end function jz
!
!==end hoppings======================================================
!
end module m_SlaterKoster
