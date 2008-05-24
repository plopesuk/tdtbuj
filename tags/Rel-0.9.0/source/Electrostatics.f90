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
    type(solutionType),intent(in) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(in) :: gen
    integer :: i
    integer :: at,n
    real(k_pr) :: aux
    integer :: from,to

    n=atomic%basis%norbitals
    if (gen%spin) then
    ! spin down
      do at=1,atomic%atoms%natoms
      aux=0.0_k_pr
        from=atomic%atoms%orbs(at,1)
        to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2
        do i=from,to
          aux=aux+sol%rho%a(i,i)
        enddo
      !spin up
        from=atomic%atoms%orbs(at,1)+atomic%basis%norbitals/2
        to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2
        do i=from,to
          aux=aux+sol%rho%a(i,i)
        enddo
        atomic%atoms%chrg(at)=aux-atomic%species%zval(atomic%atoms%sp(at))
      enddo
    else
      do at=1,atomic%atoms%natoms
        aux=0.0_k_pr
        from=atomic%atoms%orbs(at,1)
        to=-1+atomic%atoms%orbs(at,1)+atomic%species%norbs(atomic%atoms%sp(at))
        do i=from,to
          aux=aux+sol%rho%a(i,i)
        enddo
          atomic%atoms%chrg(at)=aux-atomic%species%zval(atomic%atoms%sp(at))
      enddo
    endif
  end subroutine CalcExcessCharges

!> \brief computes the dipole moment
!> \author Alin M Elena
!> \date 08/11/07, 10:11:26
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine CalcDipoles(gen,atomic,sol)
    character(len=*), parameter :: myname = 'CalcDipoles'
    type(solutionType),intent(in) :: sol
    type(atomicxType), intent(inout) :: atomic
    type(generalType), intent(in) :: gen
    integer :: i

    select case(gen%electrostatics)
      case (k_electrostaticsPoint)
        do i=1,atomic%atoms%natoms
          atomic%atoms%dx(i)= atomic%atoms%x(i)*atomic%atoms%chrg(i)
          atomic%atoms%dy(i)= atomic%atoms%y(i)*atomic%atoms%chrg(i)
          atomic%atoms%dz(i)= atomic%atoms%z(i)*atomic%atoms%chrg(i)
        enddo
      case(k_electrostaticsMultipoles)
        do i=1,atomic%atoms%natoms
            if (GetLmax(atomic%atoms%sp(i),atomic%speciesBasis,atomic%species)>0) then
  !              atomic%dx(i)= atomic%x(i)*atomic%chrg(i)+qlmr(i,1,1,density)
  !              atomic%dy(i)= atomic%y(i)*atomic%chrg(i)+qlmr(i,1,-1,density)
  !              atomic%dz(i)= atomic%z(i)*atomic%chrg(i)+qlmr(i,1,0,density)
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
!> \internal add reference
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

end module m_Electrostatics
