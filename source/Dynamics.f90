!> \brief routines that deal with the dynamics
!> \details they are mainly called by the driver routines
!> \author Alin M Elena
!> \date 10/11/07, 13:21:58
module m_Dynamics
  use m_Constants
  use m_Types
  use m_Useful
  implicit none
  private
!
  public :: KineticEnergy
  public :: ScaleVelocities
  public :: InitVelocities
!
contains
!
!> \brief computes the kinetic energy
!> \author Alin M Elena
!> \date 10/11/07, 14:34:09
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
  real (k_pr) function KineticEnergy (atomic)
    character (len=*), parameter :: myname = 'KineticEnergy'
    type (atomicxType), intent (in) :: atomic
    integer :: i, k
    KineticEnergy = 0.0_k_pr
    do k = 1, atomic%atoms%nmoving
      i = atomic%atoms%moving(k)
      KineticEnergy = KineticEnergy + 0.5_k_pr * atomic%species%mass(atomic%atoms%sp(i)) * VectorModulus (atomic%atoms%vx(i), &
     & atomic%atoms%vy(i), atomic%atoms%vz(i)) ** 2
    end do
  end function KineticEnergy
!
!> \brief scales velocities to correspond to a certain temperature
!> \author Cristian Sanchez, Alin M Elena
!> \date 10/11/07, 14:36:02
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param gen type(generalType) contains the info needed by the program to k_run
  subroutine ScaleVelocities (gen, atomic)
    character (len=*), parameter :: myname = 'ScaleVelocities'
    type (atomicxType), intent (inout) :: atomic
    type (generalType), intent (in) :: gen
    integer :: i, k
    real (k_pr) :: kenergyPresent, kenergyWanted, scaleFactor
!
    kenergyPresent = KineticEnergy (atomic)
    kenergyWanted = real (atomic%atoms%nmoving, k_pr) * 1.5_k_pr * gen%ionicTemperature * k_kb
    scaleFactor = Sqrt (kenergyWanted/kenergyPresent)
    do k = 1, atomic%atoms%nmoving
      i = atomic%atoms%moving(k)
      atomic%atoms%vx (i) = atomic%atoms%vx(i) * scaleFactor
      atomic%atoms%vy (i) = atomic%atoms%vy(i) * scaleFactor
      atomic%atoms%vz (i) = atomic%atoms%vz(i) * scaleFactor
    end do
!
  end subroutine ScaleVelocities
!
!> \brief initializes velocities according to the ionic temperature
!> \author Cristian Sanchez, Alin M Elena
!> \date 10/11/07, 15:04:53
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param gen type(generalType) contains the info needed by the program to run
!> \internal for only one particle moving you get division by zero. a better scheme
!> or a correct scheme should be found
  subroutine InitVelocities (gen, atomic, sol)
    character (len=*), parameter :: myname = 'InitVelocities'
    integer :: i, k
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    real (k_pr) :: kt, mi, kex, key, kez
    real (k_pr) :: vcmx, vcmy, vcmz, mt
!
    if ( .not. gen%ReadVelocity) then
      kt = gen%ionicTemperature * k_kb
      vcmx = 0.0_k_pr
      vcmy = 0.0_k_pr
      vcmz = 0.0_k_pr
      mt = 0.0_k_pr
      do i = 1, atomic%atoms%natoms
        mi = atomic%species%mass(atomic%atoms%sp(i))
        mt = mt + mi
      end do
      do k = 1, atomic%atoms%nmoving
        i = atomic%atoms%moving(k)
        mi = atomic%species%mass(atomic%atoms%sp(i))
! generate random vector              
        atomic%atoms%vx (i) = (ranmar(sol%seed)-0.5_k_pr) * Sqrt (kt/mi)
        atomic%atoms%vy (i) = (ranmar(sol%seed)-0.5_k_pr) * Sqrt (kt/mi)
        atomic%atoms%vz (i) = (ranmar(sol%seed)-0.5_k_pr) * Sqrt (kt/mi)
! accumulate in centre of mass velocity              
        vcmx = vcmx + atomic%atoms%vx(i) * mi
        vcmy = vcmy + atomic%atoms%vy(i) * mi
        vcmz = vcmz + atomic%atoms%vz(i) * mi
      end do
      vcmx = vcmx / mt
      vcmy = vcmy / mt
      vcmz = vcmz / mt
      do k = 1, atomic%atoms%nmoving
        i = atomic%atoms%moving(k)
! substract centre of mass velocity                
        atomic%atoms%vx (i) = atomic%atoms%vx(i) - vcmx
        atomic%atoms%vy (i) = atomic%atoms%vy(i) - vcmy
        atomic%atoms%vz (i) = atomic%atoms%vz(i) - vcmz
      end do
      kex = 0.0_k_pr
      key = 0.0_k_pr
      kez = 0.0_k_pr
      do k = 1, atomic%atoms%nmoving
        i = atomic%atoms%moving(k)
! accumulate the kinetic energy              
        mi = atomic%species%mass(atomic%atoms%sp(i))
        kex = kex + 0.5_k_pr * mi * atomic%atoms%vx(i) ** 2
        key = key + 0.5_k_pr * mi * atomic%atoms%vy(i) ** 2
        kez = kez + 0.5_k_pr * mi * atomic%atoms%vz(i) ** 2
      end do
      do k = 1, atomic%atoms%nmoving
        i = atomic%atoms%moving(k)
        atomic%atoms%vx (i) = atomic%atoms%vx(i) * Sqrt (kt/kex)
        atomic%atoms%vy (i) = atomic%atoms%vy(i) * Sqrt (kt/key)
        atomic%atoms%vz (i) = atomic%atoms%vz(i) * Sqrt (kt/kez)
      end do
    end if
  end subroutine InitVelocities
end module m_Dynamics
