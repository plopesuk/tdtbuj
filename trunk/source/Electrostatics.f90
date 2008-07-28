!> \brief electrostatics potential and field
!> \author Alin M Elena
!> \date 09/11/07, 00:01:22
module m_Electrostatics
  use m_Constants
  use m_Types
  use m_Useful
  implicit none
  private

  public :: CalcExcessCharges
  public :: CalcDipoles
  public :: BuildPotential
  public :: Charge
  public :: ChargeOnGroup
  public :: BuildField
  public :: ChargeOnL
  public :: initQvs
  public :: qlmr
  public :: vlmr
  public :: BlplR
  public :: hiujv
  public :: fip
  public :: fact3

  interface ChargeOnL
    module procedure ChargeOnLSpin, ChargeOnLNoSpin
  end interface


contains

!> \brief calculates the excess charges for all atoms
!> \details the charges are stored in atomic%atoms%chrg(:)
!> \author Alin M Elena
!> \date 08/11/07, 09:34:30
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine CalcExcessCharges(gen,atomic,sol)
    character(len=*), parameter :: myname = 'CalcExcessCharges'
    type(solutionType),intent(inout) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: gen

    integer :: at
    real(k_pr) :: aux
    ! spin down
    do at=1,atomic%atoms%natoms
      aux=PartialTrace(atomic%atoms%id(at),atomic,sol%rho,gen%spin)
      atomic%atoms%chrg(at)=aux-atomic%species%zval(atomic%atoms%sp(at))
    enddo

  end subroutine CalcExcessCharges

!> \brief computes the dipole moment
!> \author Alin M Elena
!> \date 08/11/07, 10:11:26
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine CalcDipoles(gen,atomic,sol,tb)
    character(len=*), parameter :: myname = 'CalcDipoles'
    type(solutionType),intent(inout) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: gen
    type(modelType), intent(inout) :: tb
    integer :: i

    select case(gen%electrostatics)
      case (k_electrostaticsPoint)
        do i=1,atomic%atoms%natoms
          atomic%atoms%dx(i) = atomic%atoms%x(i)*atomic%atoms%chrg(i)
          atomic%atoms%dy(i) = atomic%atoms%y(i)*atomic%atoms%chrg(i)
          atomic%atoms%dz(i) = atomic%atoms%z(i)*atomic%atoms%chrg(i)
        enddo
      case(k_electrostaticsMultipoles)
        do i=1,atomic%atoms%natoms
          if (GetLmax(atomic%atoms%sp(i),atomic%speciesBasis,atomic%species)>0) then
            atomic%atoms%dx(i)= atomic%atoms%x(i)*atomic%atoms%chrg(i)+qlmr(i,1,1,gen,sol,atomic,tb,sol%density)*sqrt(4*k_pi/3.0_k_pr)
            atomic%atoms%dy(i)= atomic%atoms%y(i)*atomic%atoms%chrg(i)+qlmr(i,1,-1,gen,sol,atomic,tb,sol%density)*sqrt(4*k_pi/3.0_k_pr)
            atomic%atoms%dz(i)= atomic%atoms%z(i)*atomic%atoms%chrg(i)+qlmr(i,1,0,gen,sol,atomic,tb,sol%density)*sqrt(4*k_pi/3.0_k_pr)
          else
            atomic%atoms%dx(i)= atomic%atoms%x(i)*atomic%atoms%chrg(i)
            atomic%atoms%dy(i)= atomic%atoms%y(i)*atomic%atoms%chrg(i)
            atomic%atoms%dz(i)= atomic%atoms%z(i)*atomic%atoms%chrg(i)
          endif
        enddo
    end select
    atomic%atoms%tdipx=0.0_k_pr
    atomic%atoms%tdipy=0.0_k_pr
    atomic%atoms%tdipz=0.0_k_pr
    do i=1,atomic%atoms%natoms
      atomic%atoms%tdipx=atomic%atoms%tdipx+atomic%atoms%dx(i)
      atomic%atoms%tdipy=atomic%atoms%tdipy+atomic%atoms%dy(i)
      atomic%atoms%tdipz=atomic%atoms%tdipz+atomic%atoms%dz(i)
    enddo
  end subroutine CalcDipoles

!> \brief builds the electrostatics potential for all the atoms
!> \details there are two kinds of electrostatics potential screened or bare Coulomb
!> they get selected by "Screened" entry
!> \author Alin M Elena
!> \date 08/11/07, 15:41:54
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \remark The screened potential is given by Klopman-Ohno approximation
!>  for more see J Phys Chem A Vol 111, No 26, pp. 5614-5621, (2007)

  subroutine BuildPotential(gen,atomic,sol)
    character(len=*), parameter :: myname = 'BuildPotential'
    type(solutionType),intent(inout) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: gen
    integer  :: i,j,k
    real(k_pr) :: ua,ub, rij,uv
    !-------------------------------------------------!
    sol%potential = 0.0_k_pr
    if (gen%screened) then
      do k=1,atomic%atoms%nscf
        i=atomic%atoms%scf(k)
        do j=1,atomic%atoms%natoms
          if (i/=j) then
            rij = Distance(atomic%atoms,i,j)
            ua=atomic%species%uinter(atomic%atoms%sp(i))
            ub=atomic%species%uinter(atomic%atoms%sp(j))
            if (abs(ua*ub)>epsilon(k_e2)) then
              uv=(1.0_k_pr/ua+1.0_k_pr/ub)/2.0_k_pr
            else
              uv=0.0_k_pr
            endif
            sol%potential(i) = sol%potential(i) + &
              k_e2*charge(j,gen,atomic,sol)/(4.0_k_pr*k_pi*k_epsilon0*&
              sqrt(rij*rij+uv*uv))
          endif
        end do
      enddo
    else
      do k=1,atomic%atoms%nscf
        i=atomic%atoms%scf(k)
        do j=1,atomic%atoms%natoms
          if (i/=j) then
            rij = Distance(atomic%atoms,i,j)
            sol%potential(i) = sol%potential(i) + &
            k_e2*charge(j,gen,atomic,sol)/(4.0_k_pr*k_pi*k_epsilon0*rij)
          endif
        end do
      enddo
    endif
  end subroutine BuildPotential

!> \brief computes the charge on an atom
!> \author Alin M Elena
!> \date 08/11/07, 15:46:03
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param at integer the atom
  real(k_pr) function charge(at,gen,atomic,sol)
    character(len=*), parameter :: myname = 'charge'
    integer, intent(in) :: at
    type(solutionType),intent(in) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(in) :: gen
    integer :: from,to,m
    real(k_pr) :: aux

    aux=0.0_k_pr
    m=atomic%basis%norbitals*(atomic%basis%norbitals-1)/2
    if (gen%spin) then
    ! spin down
      from=m+atomic%atoms%orbs(at,1)
      to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2
      aux=sum(sol%density(from:to))
    !spin up
      from=m+atomic%atoms%orbs(at,1)+atomic%basis%norbitals/2
      to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2
      aux=aux+sum(sol%density(from:to))
    else
      from=m+atomic%atoms%orbs(at,1)
      to=-1+m+atomic%atoms%orbs(at,1)+atomic%species%norbs(atomic%atoms%sp(at))
      aux=sum(sol%density(from:to))
    endif
    charge=aux
  end function charge


  subroutine ChargeOnLSpin(at,l,atomic,gen,sol,qld,qlu)
    character(len=*), parameter :: myname = 'ChargeOnLSpin'
    integer, intent(in) :: at
    integer, intent(in) :: l
    type(solutionType),intent(in) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(in) :: gen
    real(k_pr) :: qld,qlu
    integer :: from,to,m,j,i
    real(k_pr) :: aux

     qlu=0.0_k_pr
     qld=0.0_k_pr
     j=l*l
! spin down
    m=atomic%atoms%orbs(at,1)
     from=m+j
     to=-1+from+2*l+1
     do i=from,to
       qld=qld+real(sol%rho%a(i,i),k_pr)
     enddo
   !spin up
     from=m+atomic%basis%norbitals/2+j
     to=-1+from+2*l+1
     do i=from,to
       qlu=qlu+real(sol%rho%a(i,i),k_pr)
     enddo
  end subroutine ChargeOnLSpin


 subroutine ChargeOnLNoSpin(at,l,atomic,gen,sol,ql)
    character(len=*), parameter :: myname = 'ChargeOnLNoSpin'
    integer, intent(in) :: at
    integer, intent(in) :: l
    type(solutionType),intent(in) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(in) :: gen
    real(k_pr) :: ql
    integer :: from,to,m,j,i
    real(k_pr) :: aux

     ql=0.0_k_pr
     j=l*l
     m=atomic%atoms%orbs(at,1)
     from=m+j
     to=-1+from+2*l+1
     do i=from,to
       ql=ql+real(sol%rho%a(i,i),k_pr)
     enddo
  end subroutine ChargeOnLNoSpin


!> \brief builds the electrostatics potential for all the atoms
!> \details there are two kinds of electrostatics potential screened or bare Coulomb
!> they get selected by "Screened" entry
!> \author Alin M Elena
!> \date 08/11/07, 22:58:08
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine BuildField(gen,atomic,sol)
    character(len=*), parameter :: myname = 'BuildField'
    type(solutionType),intent(inout) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(inout) :: gen
    integer :: i,j,k
    real(k_pr) :: rijx, rijy, rijz, rij
    real(k_pr) :: fact,ua,ub,uv
    !-------------------------------------------------!

    sol%field = 0.0_k_pr
    if (gen%screened) then
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        do j=1,atomic%atoms%natoms
          if (i/=j) then
            call AtomDistVec(atomic%atoms,j,i,rij,rijx,rijy,rijz)
            ua=atomic%species%uinter(atomic%atoms%sp(i))
            ub=atomic%species%uinter(atomic%atoms%sp(j))
            if (abs(ua*ub)>epsilon(k_e2)) then
              uv=(1.0_k_pr/ua+1.0_k_pr/ub)/2.0_k_pr
            else
              uv=0.0_k_pr
            endif
            fact = -charge(j,gen,atomic,sol)*k_e2/(4.0_k_pr*k_pi*k_epsilon0)/&
                  (rij*rij + uv*uv)**1.5_k_pr
            sol%field(i,1) = sol%field(i,1) + fact * rijx
            sol%field(i,2) = sol%field(i,2) + fact * rijy
            sol%field(i,3) = sol%field(i,3) + fact * rijz
          endif
        end do
      end do
    else
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        do j=1,atomic%atoms%natoms
          if (i/=j) then
            call AtomDistVec(atomic%atoms,j,i,rij,rijx,rijy,rijz)
            fact =  -charge(j,gen,atomic,sol)*k_e2/(4.0_k_pr*k_pi*k_epsilon0*rij**3)
            sol%field(i,1) = sol%field(i,1) + fact * rijx
            sol%field(i,2) = sol%field(i,2) + fact * rijy
            sol%field(i,3) = sol%field(i,3) + fact * rijz
          endif
        end do
      end do
    endif
  end subroutine BuildField

!> \brief computes the charge of a group of atoms
!> \author Alin M Elena
!> \date 10/11/07, 15:40:31
!> \param group array a list of atoms
!> \param atomic type(atomicType) contains all info about the atoms
  real(k_pr) function ChargeOnGroup(group,atomic)
    character(len=*), parameter :: myname="ChargeOnGroup"
    integer, intent(in) :: group(:)
    type(atomicType) :: atomic
    integer :: i
    real(k_pr) :: sacc

    sacc=0.0_k_pr
    do i=1,size(group)
      sacc=sacc+atomic%chrg(group(i))
    enddo
    ChargeOnGroup=sacc
  end function ChargeOnGroup

   subroutine initQvs(atomic,gen,sol,tb,density)
     character(len=*), parameter :: myname = 'initQvs'
     type(atomicxType), intent(inout) :: atomic
     type(solutionType),intent(inout) :: sol
     type(modelType), intent(inout) :: tb
     real(k_pr), intent(inout) :: density(:)
     type(generalType),intent(inout) :: gen

     integer :: i,l,m,k


     do i=1,atomic%atoms%natoms
       k=GetLmax(atomic%atoms%sp(i),atomic%speciesBasis,atomic%species)
       do l=0,2*k
          do m=-l,l
            sol%delq(i)%a(idx(l,m))=compQlmR(i,l,m,sol,atomic,tb,density)
          enddo
       enddo
     enddo

     do i=1,atomic%atoms%natoms
       k=GetLmax(atomic%atoms%sp(i),atomic%speciesBasis,atomic%species)
       do l=0,2*k+1
         do m=-l,l
           sol%vs(i)%a(idx(l,m))=compVlmR(i,l,m,gen,sol,atomic,tb,density)
         enddo
       enddo
    enddo
   end subroutine initQvs

   real(k_pr) function compQlmR(at,l,m,sol,atomic,tb,density)
     character(len=*), parameter :: myname = 'compQlmR'
     integer, intent(in) :: l,m,at
     real(k_pr), intent(in) :: density(:)
     type(atomicxType), intent(inout) :: atomic
     type(solutionType),intent(inout) :: sol
     type(modelType), intent(inout) :: tb

      real(k_pr) :: sum
      integer :: k1,k2,n,l1,l2,m1,m2

      n=atomic%basis%norbitals
      sum=0.0_k_pr
      do k1=1,atomic%species%norbs(atomic%atoms%sp(at))
         l1 = atomic%basis%orbitals(atomic%atoms%orbs(at,k1))%l
         m1 = atomic%basis%orbitals(atomic%atoms%orbs(at,k1))%m
         do k2=1,atomic%species%norbs(atomic%atoms%sp(at))
           l2 = atomic%basis%orbitals(atomic%atoms%orbs(at,k2))%l
           m2 = atomic%basis%orbitals(atomic%atoms%orbs(at,k2))%m
           sum=sum+density(aidx(atomic%atoms%orbs(at,k1),atomic%atoms%orbs(at,k2),n))*&
               tb%delta(atomic%atoms%sp(at))%d(l,l1,l2)*sol%rgc(idx(l1,m1),idx(l2,m2),idx(l,m))
         enddo
      enddo
      compQlmR=sum

   end function compQlmR

  real(k_pr) function compVlmR(i,l,m,gen,sol,atomic,tb,density)
    character(len=*), parameter :: myname = 'compVlmR'

    integer, intent(in) :: i,l,m
    real(k_pr), intent(in) :: density(:)
    type(atomicxType), intent(inout) :: atomic
    type(solutionType),intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    integer :: j,lp,mp,lm
    real(k_pr) :: sum,x,y,z
    real(k_pr), allocatable :: ir(:)
    sum=0.0_k_pr
    do j=1,atomic%atoms%natoms
      if (j/=i) then
        lm=GetLmax(atomic%atoms%sp(j),atomic%speciesBasis,atomic%species)
        allocate(ir(1:(l+2*lm+2)**2))
        x=atomic%atoms%x(j)-atomic%atoms%x(i)
        y=atomic%atoms%y(j)-atomic%atoms%y(i)
        z=atomic%atoms%z(j)-atomic%atoms%z(i)
        call solidh(x,y,z,-(l+2*lm+1),ir,(l+2*lm+2)**2)
        do lp=0,2*lm
          do mp=-lp,lp
            sum=sum+qlmR(j,lp,mp,gen,sol,atomic,tb,density)*blplR(l,m,lp,mp,i,j,ir,sol)
          enddo
        enddo
        deallocate(ir)
      endif
    enddo
    if (gen%hasElectricField) then
      if (l==0) then
        sum= sum*k_e2/(4.0_k_pr*k_pi*k_epsilon0) + (atomic%atoms%x(i)*gen%E(1)+atomic%atoms%y(i)*gen%E(2)+atomic%atoms%z(i)*gen%E(3))*k_e*2.0_k_pr*sqrt(k_pi)
      elseif(l==1) then
        sum= sum*k_e2/(4.0_k_pr*k_pi*k_epsilon0) + gen%E(mod(2+m,3)+1)*k_e*sqrt(4.0_k_pr*k_pi/3.0_k_pr)
      endif
    else
      sum=sum*k_e2/(4.0_k_pr*k_pi*k_epsilon0)
    endif
    compVlmR=sum
   end function compVlmR

   real(k_pr) function qlmR(at,l,m,gen,sol,atomic,tb,density)
     character(len=*), parameter :: myname = 'qlmR'
     integer, intent(in) :: l,m,at
     real(k_pr), intent(in) :: density(:)
     type(atomicxType), intent(inout) :: atomic
     type(solutionType),intent(inout) :: sol
     type(modelType), intent(inout) :: tb
     type(generalType), intent(inout) :: gen

     if (gen%compElec) then
! we compute them on the fly
       qlmR=compQlmR(at,l,m,sol,atomic,tb,density)
     else
! they were precomputed somewhere else
       qlmR=sol%delq(at)%a(idx(l,m))
     endif
   end function qlmR



! DESCRIPTION
!    returns the structure factor for atom at1 and at2
!    |latex \begin{equation}
!    |latex B_{ll^{\prime}} = \frac{(4\pi)^2(-1)^{l^{\prime}}}{(2l+1)!!(2l^{\prime}+1)!!}
!    |latex  \sum_{m^{\prime\prime}=-l-l\prime}^{l+l\prime}\mathcal{G}_{lml^\prime m^\prime}^{l^{\prime\prime}m^{\prime\prime}}
!    |latex  \frac{X_{l^{\prime\prime}m^{\prime\prime}}(\hat{r_{12}})}{r_{12}^{l^{\prime\prime}+1}}\frac{(2l^{\prime\prime}+1)!!}{2l^{\prime\prime}+1}
!    |latex \end{equation}
!   with
! |latex $l^{\prime\prime}=l+l^{\prime}$
! USES
! AUTHOR
! Alin M. Elena (Belfast)
! CREATION DATE
! April 2006
! HISTORY
!*****

   real(k_pr) function blplR(l,m,lp,mp,at1,at2,ir,sol)
     character(len=*), parameter :: myname = 'blplR'
     integer,intent(in) :: l,m,lp,mp,at1,at2
     integer :: p,i,lpp
     real(k_pr),intent(inout) :: ir(:)
     type(solutionType),intent(inout) :: sol
     real(k_pr) :: sum,x,y,z,r

      lpp=l+lp
      r=sqrt(real(2*lpp+1,k_pr)/(4.0_k_pr*k_pi))*fact3(lpp,sol)/&
         (2.0_k_pr*lpp+1.0_k_pr)
!     r=sqrt(x*x+y*y+z*z)**(lpp+1)
      sum=0.0_k_pr
      do p=-lpp,lpp
         sum=sum+sol%rgc(idx(l,m),idx(lp,mp),idx(lpp,p))*&
            ir(idxy(lpp,p))*r
      enddo
      blplR=sum*(4.0_k_pr*k_pi)**2*(-1)**lp/(fact3(l,sol)*fact3(lp,sol))
   end function blplR



  real(k_pr) function vlmR(at,l,m,gen,sol,atomic,tb,density)
    character(len=*), parameter :: myname = 'vlmR'
    integer, intent(in) :: l,m,at
    real(k_pr), intent(inout) :: density(:)
    type(atomicxType), intent(inout) :: atomic
    type(solutionType),intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    type(generalType), intent(inout) :: gen

    if (gen%compElec) then
! we compute them on the fly
      vlmR=compVlmR(at,l,m,gen,sol,atomic,tb,density)
    else
! they were precomputed somewhere else
      vlmR=sol%vs(at)%a(idx(l,m))
    endif
  end function vlmR

!****f*   tb_sppna/fact3()
! NAME
! fact3
! SYNOPSIS
! fact3(n)
! INPUTS
! integer n
!
! DESCRIPTION
!    returns (2n+1)!! where
!    |latex \begin{equation}
!    |latex  (2n+1)!!=1\cdot 3 \cdot  ...       \cdot (2n+1)=\frac{(2n+1)!}{n!2^n}
!    |latex \end{equation}
! USES
! AUTHOR
! Alin M. Elena (Belfast)
! CREATION DATE
! 1st of May 2006
! HISTORY
!*****
  real(k_pr) function fact3(n,sol)
    character(len=*), parameter :: myname = 'fact3'
    integer, intent(in) :: n
    type(solutionType), intent(in) :: sol

    fact3 = sol%fact(2*n+1)/(sol%fact(n)*2.0_k_pr**n)
  end function fact3

  real(k_pr) function hiujv(at,lu,mu,lv,mv,gen,sol,atomic,tb,density)
    character(len=*), parameter :: myname = 'hiujv'
    integer, intent(in) :: at,lu,mu,lv,mv
    real(k_pr), intent(inout) :: density(:)
    type(atomicxType), intent(inout) :: atomic
    type(solutionType),intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    integer :: l,m,i
    real(k_pr) :: sum

    sum=0.0_k_pr
    do l=0,2*GetLmax(atomic%atoms%sp(at),atomic%speciesBasis,atomic%species)
      do m=-l,l
        sum=sum+sol%rgc(idx(lu,mu),idx(lv,mv),idx(l,m))*&
            tb%delta(atomic%atoms%sp(at))%d(l,lu,lv)*&
            vlmR(at,l,m,gen,sol,atomic,tb,density)
      enddo
    enddo
!     hiujv=sum*k_e2/(4.0_k_pr*k_pi*k_epsilon0)
    hiujv=sum
  end function hiujv

   real(k_pr) function fip(at,li,mi,p,gen,sol,atomic,tb,density)
    character(len=*), parameter :: myname = 'fip'
    integer, intent(in) :: at,li,p,mi
    real(k_pr), intent(inout) :: density(:)
    type(atomicxType), intent(inout) :: atomic
    type(solutionType),intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    type(generalType), intent(inout) :: gen
    integer :: i,lp,mp
    real(k_pr) :: sum

    sum=0.0_k_pr
    lp=li+1
    do mp=-lp,lp
      sum=sum+sol%rgc(idx(1,p),idx(li,mi),idx(lp,mp))*&
          vlmR(at,lp,mp,gen,sol,atomic,tb,density)
    enddo
    fip=sum
   end function fip

end module m_Electrostatics
