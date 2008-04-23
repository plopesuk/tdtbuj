!> \brief contains the routines that "drive" the program
!> \author Alin M Elena
!> \date 02/11/07, 18:04:08

module m_DriverRoutines
  use m_Constants
  use m_Useful
  use m_Types
  use m_TightBinding
  use m_Gutenberg
  use m_LinearAlgebra
  use m_SCF
  use m_Electrostatics
  use m_Hamiltonian
  use m_Dynamics
  use m_DensityMatrix
  use m_LBFGS
  implicit none

  private

  public :: SinglePoint
  public :: BornOppenheimerDynamics
  public :: EhrenfestDynamics
  public :: EhrenfestDynamicsDamped
  public :: BFGS

  contains

!> \brief the driver for total energy and electronic structure calculations
!> \author Alin M Elena
!> \date 02/11/07, 18:03:43
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine SinglePoint(ioLoc,genLoc,atomic,tbMod,sol)
    character(len=*), parameter :: sMyName="SinglePoint"
    type(ioType), intent(inout) :: ioLoc
    type(generalType), intent(inout) :: genLoc
    type(atomicxType), intent(inout) :: atomic
    type(modelType), intent(inout) :: tbMod
    type(solutionType), intent(inout) :: sol
    real(k_pr) :: eenergy,renergy, minusts,scfE

 eenergy=0.0_k_pr
 renergy=0.0_k_pr
 scfE=0.0_k_pr

    if (genLoc%spin)  then
      call InitMagneticMoment(atomic)
      write(ioLoc%uout,"(a)")"Initial Magnetic Moment"
      call PrintMagneticMoment(atomic,sol,.false.,ioLoc)
    endif        
    call FullSCF(ioLoc,genLoc,atomic,tbMod,sol)    
    eenergy = ElectronicEnergy(genLoc,sol,ioLoc)
    renergy = RepulsiveEnergy(genLoc,atomic%atoms,tbMod)
    minusts = sol%electronicEntropy
   
    write(ioLoc%uout,'(/a)')&
         '--Single Point Run----------------------------------------------'
    call PrintAtoms(ioLoc,genLoc,atomic)
    call PrintCharges(genLoc,atomic,ioLoc)
    call PrintDipoles(atomic,ioLoc)
    call PrintMagneticMoment(atomic,sol,.false.,ioLoc)
    write(ioLoc%uout,"(/a,f16.8,/a,f16.8)")&
          'Electronic energy = ',eenergy, &
       'Repulsive energy  = ',renergy
    if (genLoc%scf) then
      scfE = ScfEnergy(genLoc,atomic,sol,ioLoc)
    else
      scfE = 0.0_k_pr
    endif
    write(ioLoc%uout,"(a,f16.8,/a,f16.8,/a,f16.8,/)") &
       'Total SCF energy  = ',scfE, &
          '-TS               = ',minusts, &
          'Total energy      = ',eenergy+renergy+scfE+minusts
    call PrintForces(atomic%atoms,ioLoc)
    write(ioLoc%uout,'(/a/)')&
          '----------------------------------------------------------------'
    sol%totalEnergy=eenergy+renergy+minusts+scfE

  end subroutine SinglePoint

!> \brief driver routine for verlet velocity Born-Oppenheimer molecular dynamics
!> \author Alin M Elena
!> \date 10/11/07, 13:18:53
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine BornOppenheimerDynamics(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'BornOppenheimerDynamics'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    integer :: aniunit, eneunit,xunit,runit
    real(k_pr) :: eenergy,renergy,kenergy,penergy,scfE,minusts
    integer  :: i,istep,k
    real(k_pr) :: dt,mi
    character(len=k_ml) :: saux
    aniunit=GetUnit()
    eneunit=GetUnit()
    xunit=GetUnit()
    runit=GetUnit()
    write(io%uout,'(/a/)')&
      '--Velocity Verlet Born-Oppenheimer Dynamics---------------------'

    if (gen%writeAnimation) then
      open(file='bo_dyn.xyz',unit=aniunit)
!          open(unit=xunit,file="bo_dyn.gcd",form="UNFORMATTED",status="unknown",action="write")
!          call write_header(xunit,gen%nsteps,gen%deltat,atomic%natoms)
!          open(unit=runit,file="bo_dyn.rho",form="UNFORMATTED",status="unknown",action="write")
!          call write_header_rho(runit,gen%nsteps,gen%deltat,n)
      call PrintXYZ(aniunit,atomic,.false.,"T = 0.0")
    endif
    open(file='bo_dyn.ENE',unit=eneunit)
    ! initialize forces and velocities
    call SinglePoint(io,gen,atomic,tb,sol)
    call InitVelocities(gen,atomic,sol)
    dt = gen%deltat
    !write out the initial vels, forces and energies
    call PrintVelocities(io,atomic)
    call PrintForces(atomic%atoms,io)
    eenergy = ElectronicEnergy(gen,sol,io)
    renergy = RepulsiveEnergy(gen,atomic%atoms,tb)
    minusts = sol%electronicEntropy
    if (gen%scf) then
      scfE = ScfEnergy(gen,atomic,sol,io)
    else
      scfE = 0.0_k_pr
    endif

    write(io%uout,'(/a,f13.6,/a,f13.6,/a,f13.6,/a,f25.18,/a,f25.18,/)')&
       'Electronic energy = ',eenergy, &
      'Repulsive energy  = ',renergy, &
      'SCF energy        = ',scfE, &
         '-TS               = ',minusts, &
         'Total energy      = ',eenergy+renergy+scfE+minusts
    sol%totalEnergy=eenergy+renergy+minusts+scfE
    ! this is the time loop
    write(eneunit,'(a1,a28,6a29)')"#","Time",  "Repulsive Energy ",  "Electronic Energy",  "SCF Energy",&
      "-TS","Kinetic Energy",  "Total Energy"
    do istep=1,gen%nsteps
      !set global time variable
      gen%CurrSimTime = gen%CurrSimTime*k_time2SI
       ! calculates positions at t+dt
!       write(io%uout,*)"timestep"
!       call PrintCoordinates(io,atomic)
!       call PrintForces(atomic%atoms,io)
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        mi = atomic%species%mass(atomic%atoms%sp(i))
        atomic%atoms%x(i) = atomic%atoms%x(i) &
          + dt * atomic%atoms%vx(i) &
          + 0.5_k_pr * dt*dt * atomic%atoms%fx(i) / mi
        atomic%atoms%y(i) = atomic%atoms%y(i) &
          + dt * atomic%atoms%vy(i) &
          + 0.5_k_pr * dt*dt * atomic%atoms%fy(i) / mi
        atomic%atoms%z(i) = atomic%atoms%z(i) &
          + dt * atomic%atoms%vz(i) &
          + 0.5_k_pr * dt*dt * atomic%atoms%fz(i) / mi
      enddo
!       call PrintVelocities(io,atomic)
!       call PrintCoordinates(io,atomic)
      ! store forces at t in fold
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        atomic%atoms%fxo(i) = atomic%atoms%fx(i)
        atomic%atoms%fyo(i) = atomic%atoms%fy(i)
        atomic%atoms%fzo(i) = atomic%atoms%fz(i)
      enddo
      ! calculate forces at t+dt
      call SinglePoint(io,gen,atomic,tb,sol)
      ! calculate velocities at t+dt
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        mi = atomic%species%mass(atomic%atoms%sp(i))
        atomic%atoms%vx(i) = atomic%atoms%vx(i) &
          + 0.5_k_pr * dt * (atomic%atoms%fx(i) + atomic%atoms%fxo(i)) / mi
        atomic%atoms%vy(i) = atomic%atoms%vy(i) &
          + 0.5_k_pr * dt * (atomic%atoms%fy(i) + atomic%atoms%fyo(i)) / mi
        atomic%atoms%vz(i) = atomic%atoms%vz(i) &
          + 0.5_k_pr * dt * (atomic%atoms%fz(i) + atomic%atoms%fzo(i)) / mi
      enddo
       ! scale velocities
      if (gen%scaleVelocities) then
        call ScaleVelocities(gen,atomic)
      endif
       ! now calculate the energies
      eenergy = ElectronicEnergy(gen,sol,io)
      renergy = RepulsiveEnergy(gen,atomic%atoms,tb)
      minusts = sol%electronicEntropy
      if (gen%scf) then
        scfE = ScfEnergy(gen,atomic,sol,io)
      else
        scfE = 0.0_k_pr
      endif
      penergy = eenergy + renergy + scfE + minusts
      kenergy = KineticEnergy(atomic)
      if (gen%writeAnimation) then
!             call BuildDensity(density)
!             call calc_charges(rho)
!             call calc_dipoles(density)
!             call write_frame(xunit)
!             call writeAnimation_frame(aniunit)
! !           call write_frame_rho(runit)
        write(saux,'(a,f0.8)')"Time = ",gen%CurrSimTime
        call PrintXYZ(aniunit,atomic,.false.,trim(saux))
      endif
      if (io%Verbosity >= k_HighVerbos) then
        call PrintCoordinates(io,atomic)
        call PrintForces(atomic%atoms,io)
      endif
      write(io%uout,'(i5,a,f13.6,a,f13.6,a,f13.6)')&
            istep,' P = ',penergy,' K = ',kenergy,' E = ',penergy+kenergy
      write(eneunit,'(7f29.18)')gen%currSimTime,renergy,eenergy,scfE,minusts,kenergy,penergy+kenergy
    enddo

    if (gen%writeAnimation) then
         close(aniunit)
!          close(xunit)
!          close(runit)
    endif
    close(eneunit)

!      call write_fdf_coords()
    write(io%uout,'(/a)')&
       'Final coordinates:'
    call PrintCoordinates(io,atomic)
    write(io%uout,'(/a/)')&
         '----------------------------------------------------------------'
   end subroutine BornOppenheimerDynamics

!> \brief driver routine for verlet velocity Ehrenfest molecular dynamics
!> \author Alin M Elena
!> \date 10/11/07, 15:22:53
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!> \internal it aborts if the calculation is not scf
  subroutine EhrenfestDynamics(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'EhrenfestDynamics'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    integer :: aniunit, eneunit, popunit, xunit,runit, accUnit, donUnit, spacUnit
    real(k_pr) :: eenergy,renergy,kenergy,penergy,scfE
    integer  :: i,istep,k
    real(k_pr) :: dt,mi
    complex(k_pr) :: ihbar,trrho,st
    type(matrixType) :: rhoold,rhodot,rhonew,rho0
    real(k_pr) ::biasFactor,bfa
    character(len=k_ml) :: saux

    aniunit=GetUnit()
    eneunit=GetUnit()
    popunit=GetUnit()
    xunit=GetUnit()
    runit=GetUnit()
    accUnit=GetUnit()
    donUnit=GetUnit()
    spacUnit=GetUnit()
    write(io%uout,'(/a/)')&
         '--Velocity Verlet Ehrenfest Dynamics----------------------------'
    if (gen%writeAnimation) then
      open(file='eh_dyn.xyz',unit=aniunit)
      open(unit=xunit,file="eh_dyn.gcd",form="UNFORMATTED",status="unknown",action="write")
      open(unit=accUnit,file="eacceptor.dat",status="unknown",action="write")
      open(unit=donUnit,file="edonor.dat",status="unknown",action="write")
      open(unit=spacUnit,file="espacer.dat",status="unknown",action="write")
!      call write_header(xunit,gen%nsteps,gen%deltat,atomic%natoms)
      open(unit=runit,file="eh_dyn.rho",form="UNFORMATTED",status="unknown",action="write")
    endif
    open(file='eh_dyn.ENE',unit=eneunit)
    open(file='eh_dyn.POP',unit=popunit)
    ! prepare the electronic subsystem,
    ! which is an eigenstate of the hamitonian
    ! with the bias
    call FullSCF(io,gen,atomic,tb,sol)
    dt = gen%deltat
    ! initialize DM storage spaces
    call CreateMatrix(rhoold,sol%h%dim,.true.)
    call CreateMatrix(rhodot,sol%h%dim,.true.)
    call CreateMatrix(rhonew,sol%h%dim,.true.)
    call CreateMatrix(rho0,sol%h%dim,.true.)

    if (gen%writeAnimation) then
      call BuildDensity(atomic,sol)
      call CalcExcessCharges(gen,atomic,sol)
      call CalcDipoles(gen,atomic,sol)
      call PrintXYZ(aniunit,atomic,.false.,"T=0.0")
      write(accUnit,*) "0.0", ChargeOnGroup(atomic%atoms%acceptor,atomic%atoms)
      write(donUnit,*) "0.0", ChargeOnGroup(atomic%atoms%donor,atomic%atoms)
      write(spacUnit,*) "0.0", ChargeOnGroup(atomic%atoms%spacer,atomic%atoms)
    endif

    ! now build a hamiltonian with no bias
    call BuildHamiltonian(io,gen,atomic,tb,sol)
    call CopyMatrix(rho0,sol%rho,io)
    if (gen%BiasRampSteps>0) then
      call AddBias(1.0_k_pr,atomic,sol)
    endif
    call CopyMatrix(sol%hin,sol%h,io)
    call BuildDensity(atomic,sol)
!          if (.not.gen%comp_elec) then
!             if (gen%electrostatics==tbu_multi) call init_qvs(density)
!          endif
    call AddH2(gen,atomic,sol,tb,io)
    ihbar = cmplx(0.0_k_pr,-1.0_k_pr/k_hbar,k_pr)
    ! go back in time one step for the DM integration
    call Commutator(rhodot,sol%h,sol%rho,io)
    st = cmplx(-dt,0.0_k_pr,k_pr)
    call ScalarTMatrix(ihbar*st,rhodot,io)
    call MatrixCeaApbB(rhoold,sol%rho,rhodot,k_cone,k_cone,io)
      ! calculate the forces from the prepared DM and the present H
    call CopyMatrix(sol%h,sol%hin,io)
    call ZeroForces(atomic)
    call RepulsiveForces(gen,atomic%atoms,tb)
    call electronicForces(atomic,gen,tb,sol,io)
    ! initialize the velocities
    call InitVelocities(gen,atomic,sol)
    ! now we are ready to start the dynamics,
    ! we have the forces, velocities, positions
    ! and rho at time t
    write(eneunit,'(a1,a24,6a25)')"#","Time",  "Repuilsive Energy ",  "Electronic Energy",  "SCF Energy",&
      "Kinetic Energy",  "Total Energy",  "No of Electrons"
    st = cmplx(2.0_k_pr*dt,0.0_k_pr,k_pr)
    do istep=1,gen%nsteps
   !set global time variable
      gen%CurrSimTime = (istep-1)*dt*k_time2SI
      call BuildDensity(atomic,sol)
!             if (.not.gen%comp_elec) then
!                if (gen%electrostatics==tbu_multi) call init_qvs(density)
!             endif
      call AddH2(gen,atomic,sol,tb,io)

      call ZeroMatrix(rhodot,io)
      call Commutator(rhodot,sol%h,sol%rho,io)
      call ScalarTMatrix(ihbar*st,rhodot,io)
      call MatrixCeaApbB(rhonew,rhoold,rhodot,k_cone,k_cone,io)
       ! at this point rho contains the rho at time=t
       ! propagate the positions
       ! calculates positions at t+dt
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        mi = atomic%species%mass(atomic%atoms%sp(i))
        atomic%atoms%x(i) = atomic%atoms%x(i) &
         + dt * atomic%atoms%vx(i) &
         + 0.5_k_pr * dt*dt * atomic%atoms%fx(i) / mi
        atomic%atoms%y(i) = atomic%atoms%y(i) &
         + dt * atomic%atoms%vy(i) &
         + 0.5_k_pr * dt*dt * atomic%atoms%fy(i) / mi
        atomic%atoms%z(i) = atomic%atoms%z(i) &
         + dt * atomic%atoms%vz(i) &
         + 0.5_k_pr * dt*dt * atomic%atoms%fz(i) / mi
      enddo
       ! store forces at t in fold
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        atomic%atoms%fxo(i) = atomic%atoms%fx(i)
        atomic%atoms%fyo(i) = atomic%atoms%fy(i)
        atomic%atoms%fzo(i) = atomic%atoms%fz(i)
      enddo
       ! shuffle the DMs, rho is now rho(t+dt)
      call CopyMatrix(rhoold,sol%rho,io)
      call CopyMatrix(sol%rho,rhonew,io)
           ! calculate forces at t+dt
      call BuildHamiltonian(io,gen,atomic,tb,sol)
       ! Ramp for the bias
      if ((gen%BiasRampSteps>0).and.((istep)<=gen%BiasRampSteps)) then
        bfa = real(istep,k_pr)/real(gen%BiasRampSteps,k_pr)
        biasFactor = -(bfa-1)**3 * (1 + 3*bfa + 6*bfa**2)
        call AddBias(biasFactor,atomic,sol)
      endif

      call BuildDensity(atomic,sol)
!             if (.not.gen%comp_elec) then
!                if (gen%electrostatics==tbu_multi) call init_qvs(density)
!             endif

      call ZeroForces(atomic)
      call RepulsiveForces(gen,atomic%atoms,tb)
      call electronicForces(atomic,gen,tb,sol,io)
      call ScfForces(gen,atomic,sol,io)
       ! calculate velocities at t+dt
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        mi = atomic%species%mass(atomic%atoms%sp(i))
        atomic%atoms%vx(i) = atomic%atoms%vx(i) &
         + 0.5_k_pr * dt * (atomic%atoms%fx(i) + atomic%atoms%fxo(i)) / mi
        atomic%atoms%vy(i) = atomic%atoms%vy(i) &
         + 0.5_k_pr * dt * (atomic%atoms%fy(i) + atomic%atoms%fyo(i)) / mi
        atomic%atoms%vz(i) = atomic%atoms%vz(i) &
         + 0.5_k_pr * dt * (atomic%atoms%fz(i) + atomic%atoms%fzo(i)) / mi

      enddo

       ! scale velocities
      if (gen%scaleVelocities) then
        call ScaleVelocities(gen,atomic)
      endif
      if (gen%spin) then
        trrho   = MatrixTrace(sol%rho,io)
      else
        trrho   = 2.0_k_pr * MatrixTrace(sol%rho,io)
      endif
      eenergy = ElectronicEnergy(gen,sol,io)
      renergy = RepulsiveEnergy(gen,atomic%atoms,tb)
      if (gen%scf) then
        scfE = ScfEnergy(gen,atomic,sol,io)
      else
        scfE = 0.0_k_pr
      endif

      penergy = eenergy + renergy + scfE
      kenergy = KineticEnergy(atomic)

      if (gen%writeAnimation) then
        call CalcExcessCharges(gen,atomic,sol)
        call CalcDipoles(gen,atomic,sol)
!               call writeAnimation_frame(aniunit)
!               call write_frame(xunit)
!               call write_frame_rho(runit,rho0)
        write(accUnit,*) gen%CurrSimTime, ChargeOnGroup(atomic%atoms%acceptor,atomic%atoms)
        write(donUnit,*) gen%CurrSimTime, ChargeOnGroup(atomic%atoms%donor,atomic%atoms)
        write(spacUnit,*) gen%CurrSimTime, ChargeOnGroup(atomic%atoms%spacer,atomic%atoms)
        write(saux,'(a,f0.8)')"Time = ", gen%CurrSimTime
        call PrintXYZ(aniunit,atomic,.false.,trim(saux))
      endif
!            call write_currents(h,rho)
!            if (mod(istep,50).eq.1) call write_rho_eigenvalues(rho)
      write(eneunit,'(7f25.18)')gen%CurrSimTime,renergy,eenergy,scfE,kenergy,penergy+kenergy,real(trrho)

    enddo !istep loop
    if (gen%writeAnimation) then
      close(aniunit)
      close(xunit)
      close(runit)
      close(accUnit)
      close(donUnit)
      close(spacUnit)
    endif
    close(eneunit)
    close(popunit)
    call DestroyMatrix(rhoold,io)
    call DestroyMatrix(rhodot,io)
    call DestroyMatrix(rhonew,io)
    call DestroyMatrix(rho0,io)

    write(io%uout,'(/a/)')&
      'End Velocity Verlet-------------------------------------------------------------'
  end subroutine EhrenfestDynamics
!> \brief driver routine for verlet velocity Ehrenfest molecular dynamics
!> damped
!> \author Alin M Elena
!> \date 10/11/07, 15:22:53
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine EhrenfestDynamicsDamped(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'EhrenfestDynamicsDamped'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    integer :: aniunit, eneunit, popunit, xunit,runit, accUnit, donUnit, spacUnit
    real(k_pr) :: eenergy,renergy,kenergy,penergy,scfE
    integer  :: i,istep,k
    real(k_pr) :: dt,mi
    complex(k_pr) :: ihbar,trrho,st
    type(matrixType) :: rhoold,rhodot,rhonew,rho0,deltaRho
    real(k_pr) ::biasFactor,bfa,gamma
    character(len=k_ml) :: saux

    aniunit=GetUnit()
    eneunit=GetUnit()
    popunit=GetUnit()
    xunit=GetUnit()
    runit=GetUnit()
    accUnit=GetUnit()
    donUnit=GetUnit()
    spacUnit=GetUnit()
    write(io%uout,'(/a/)')&
         '--Velocity Verlet Ehrenfest Dynamics Damped----------------------------'
    gamma=-gen%Gamma
    if (gen%writeAnimation) then
      open(file='eh_dyn.xyz',unit=aniunit)
      open(unit=xunit,file="eh_dyn.gcd",form="UNFORMATTED",status="unknown",action="write")
      open(unit=accUnit,file="eacceptor.dat",status="unknown",action="write")
      open(unit=donUnit,file="edonor.dat",status="unknown",action="write")
      open(unit=spacUnit,file="espacer.dat",status="unknown",action="write")
!      call write_header(xunit,gen%nsteps,gen%deltat,atomic%natoms)
      open(unit=runit,file="eh_dyn.rho",form="UNFORMATTED",status="unknown",action="write")
    endif
    open(file='eh_dyn.ENE',unit=eneunit)
    open(file='eh_dyn.POP',unit=popunit)
    ! prepare the electronic subsystem,
    ! which is an eigenstate of the hamitonian
    ! with the bias
    gen%lIsExcited=.false.
    call FullSCF(io,gen,atomic,tb,sol)
    atomic%atoms%chrg0=atomic%atoms%chrg
    dt = gen%deltat
    ! initialize DM storage spaces
    call CreateMatrix(rhoold,sol%h%dim,.true.)
    call CreateMatrix(rhodot,sol%h%dim,.true.)
    call CreateMatrix(rhonew,sol%h%dim,.true.)
    call CreateMatrix(rho0,sol%h%dim,.true.)
    call CreateMatrix(deltaRho,sol%h%dim,.true.)
    call CopyMatrix(rho0,sol%rho,io)
    call CreateDensityMatrixExcited(gen,atomic,sol,io)
    ! get the starting density matrix
!     call GetRho(sol%rho)
    if (gen%writeAnimation) then
      call BuildDensity(atomic,sol)
      call CalcExcessCharges(gen,atomic,sol)
      call CalcDipoles(gen,atomic,sol)
!      call write_frame(xunit)
!        call write_frame_rho(runit,rho0)
!            call writeAnimation_frame(aniunit)
      write(accUnit,*) "0.0", ChargeOnGroup(atomic%atoms%acceptor,atomic%atoms)
      write(donUnit,*) "0.0", ChargeOnGroup(atomic%atoms%donor,atomic%atoms)
      write(spacUnit,*) "0.0", ChargeOnGroup(atomic%atoms%spacer,atomic%atoms)
      call PrintXYZ(aniunit,atomic,.false.,"T = 0.0")
    endif
!          if (.not.gen%comp_elec) then
!             if (gen%electrostatics==tbu_multi) call init_qvs(density)
!          endif
    ! now build a hamiltonian with no bias
    call BuildHamiltonian(io,gen,atomic,tb,sol)
    call MatrixCeaApbB(deltaRho,sol%rho,rho0,k_cone,-k_cone,io)
    if (gen%BiasRampSteps>0) then
      call AddBias(1.0_k_pr,atomic,sol)
    endif
    call CopyMatrix(sol%hin,sol%h,io)
    call BuildDensity(atomic,sol)
    call AddH2(gen,atomic,sol,tb,io)

    ihbar = cmplx(0.0_k_pr,-1.0_k_pr/k_hbar,k_pr)
    ! go back in time one step for the DM integration
    call Commutator(rhodot,sol%h,sol%rho,io)
    st = cmplx(-dt,0.0_k_pr,k_pr)
    call ScalarTMatrix(ihbar*st,rhodot,io)
    call ScalarTMatrix(gamma*st,deltaRho,io)
    call MatrixCeaApbB(rhoold,sol%rho,rhodot,k_cone,k_cone,io)
    call MatrixCeaApbB(rhoold,rhoold,deltaRho,k_cone,k_cone,io)
      ! calculate the forces from the prepared DM and the present H
    call CopyMatrix(sol%h,sol%hin,io)
    call ZeroForces(atomic)
    call RepulsiveForces(gen,atomic%atoms,tb)
    call electronicForces(atomic,gen,tb,sol,io)

    ! initialize the velocities
    call InitVelocities(gen,atomic,sol)
    ! now we are ready to start the dynamics,
    ! we have the forces, velocities, positions
    ! and rho at time t
    st = cmplx(2.0_k_pr*dt,0.0_k_pr,k_pr)
    write(eneunit,'(a1,a24,6a25)')"#","Time",  "Repuilsive Energy ",  "Electronic Energy",  "SCF Energy",&
      "Kinetic Energy",  "Total Energy",  "No of Electrons"
    do istep=1,gen%nsteps
   !set global time variable
      gen%CurrSimTime = (istep-1)*dt*k_time2SI
      call BuildDensity(atomic,sol)
!             if (.not.gen%comp_elec) then
!                if (gen%electrostatics==tbu_multi) call init_qvs(density)
!             endif
      call AddH2(gen,atomic,sol,tb,io)
      call ZeroMatrix(rhodot,io)
      call Commutator(rhodot,sol%h,sol%rho,io)
      call MatrixCeaApbB(deltaRho,sol%rho,rho0,k_cone,-k_cone,io)
      if (mod(istep,gen%EulerSteps)==0) then
        !euler step
        call ScalarTMatrix(ihbar*cmplx(dt,0.0_k_pr,k_pr),rhodot,io)
        call ScalarTMatrix(cmplx(gamma*dt,0.0_k_pr,k_pr),deltaRho,io)
        call MatrixCeaApbB(rhonew,sol%rho,rhodot,k_cone,k_cone,io)
      else !verlet step
        call ScalarTMatrix(ihbar*st,rhodot,io)
        call ScalarTMatrix(gamma*st,deltaRho,io)
        call MatrixCeaApbB(rhonew,rhoold,rhodot,k_cone,k_cone,io)
      end if
      call MatrixCeaApbB(rhonew,rhonew,deltaRho,k_cone,k_cone,io)
      ! at this point rho contains the rho at time=t
       ! propagate the positions
       ! calculates positions at t+dt
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        mi = atomic%species%mass(atomic%atoms%sp(i))
        atomic%atoms%x(i) = atomic%atoms%x(i) &
         + dt * atomic%atoms%vx(i) &
         + 0.5_k_pr * dt*dt * atomic%atoms%fx(i) / mi
        atomic%atoms%y(i) = atomic%atoms%y(i) &
         + dt * atomic%atoms%vy(i) &
         + 0.5_k_pr * dt*dt * atomic%atoms%fy(i) / mi
        atomic%atoms%z(i) = atomic%atoms%z(i) &
         + dt * atomic%atoms%vz(i) &
         + 0.5_k_pr * dt*dt * atomic%atoms%fz(i) / mi
      enddo
       ! store forces at t in fold
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        atomic%atoms%fxo(i) = atomic%atoms%fx(i)
        atomic%atoms%fyo(i) = atomic%atoms%fy(i)
        atomic%atoms%fzo(i) = atomic%atoms%fz(i)
      enddo
       ! shuffle the DMs, rho is now rho(t+dt)
      call CopyMatrix(rhoold,sol%rho,io)
      call CopyMatrix(sol%rho,rhonew,io)
           ! calculate forces at t+dt
      call BuildHamiltonian(io,gen,atomic,tb,sol)
       ! Ramp for the bias
      if ((gen%BiasRampSteps>0).and.((istep)<=gen%BiasRampSteps)) then
        bfa = real(istep,k_pr)/real(gen%BiasRampSteps,k_pr)
        biasFactor = -(bfa-1)**3 * (1 + 3*bfa + 6*bfa**2)
        call AddBias(biasFactor,atomic,sol)
      endif

      call BuildDensity(atomic,sol)
!             if (.not.gen%comp_elec) then
!                if (gen%electrostatics==tbu_multi) call init_qvs(density)
!             endif
      call ZeroForces(atomic)
      call RepulsiveForces(gen,atomic%atoms,tb)
      call electronicForces(atomic,gen,tb,sol,io)
      call ScfForces(gen,atomic,sol,io)
       ! calculate velocities at t+dt
      do k=1,atomic%atoms%nmoving
        i=atomic%atoms%moving(k)
        mi = atomic%species%mass(atomic%atoms%sp(i))
        atomic%atoms%vx(i) = atomic%atoms%vx(i) &
         + 0.5_k_pr * dt * (atomic%atoms%fx(i) + atomic%atoms%fxo(i)) / mi
        atomic%atoms%vy(i) = atomic%atoms%vy(i) &
         + 0.5_k_pr * dt * (atomic%atoms%fy(i) + atomic%atoms%fyo(i)) / mi
        atomic%atoms%vz(i) = atomic%atoms%vz(i) &
         + 0.5_k_pr * dt * (atomic%atoms%fz(i) + atomic%atoms%fzo(i)) / mi
      enddo
       ! scale velocities
      if (gen%scaleVelocities) then
        call ScaleVelocities(gen,atomic)
      endif
      if (gen%spin) then
        trrho   = MatrixTrace(sol%rho,io)
      else
        trrho   = 2.0_k_pr * MatrixTrace(sol%rho,io)
      endif
      eenergy = ElectronicEnergy(gen,sol,io)
      renergy = RepulsiveEnergy(gen,atomic%atoms,tb)
      if (gen%scf) then
        scfE = ScfEnergy(gen,atomic,sol,io)
      else
        scfE = 0.0_k_pr
      endif

      penergy = eenergy + renergy + scfE
      kenergy = KineticEnergy(atomic)

      if (gen%writeAnimation) then
        call CalcExcessCharges(gen,atomic,sol)
        call CalcDipoles(gen,atomic,sol)
!               call write_frame_rho(runit,rho0)
        write(accUnit,*) gen%CurrSimTime, ChargeOnGroup(atomic%atoms%acceptor,atomic%atoms)
        write(donUnit,*) gen%CurrSimTime, ChargeOnGroup(atomic%atoms%donor,atomic%atoms)
        write(spacUnit,*) gen%CurrSimTime, ChargeOnGroup(atomic%atoms%spacer,atomic%atoms)
        write(saux,'(a,f0.8)')"Time = ", gen%CurrSimTime
        call PrintXYZ(aniunit,atomic,.false.,trim(saux))
      endif
!            if (mod(istep,50).eq.1) call write_rho_eigenvalues(rho)
      write(eneunit,'(7f25.18)')gen%CurrSimTime,renergy,eenergy,scfE,kenergy,penergy+kenergy,real(trrho)
    enddo !istep loop
    if (gen%writeAnimation) then
      close(aniunit)
      close(xunit)
      close(runit)
      close(accUnit)
      close(donUnit)
      close(spacUnit)
    endif
    close(eneunit)
    close(popunit)
    call DestroyMatrix(rhoold,io)
    call DestroyMatrix(rhodot,io)
    call DestroyMatrix(rhonew,io)
    call DestroyMatrix(rho0,io)
    call DestroyMatrix(deltaRho,io)

    write(io%uout,'(/a/)')&
      'End Ehrenfest Damped-------------------------------------------------------------'
  end subroutine EhrenfestDynamicsDamped

  subroutine GetRho(rho)
    type(matrixType), intent(inout) :: rho
    integer :: ua, ub, ud,i,j,shift,n,m,ci,cj
    real(k_pr), allocatable :: a(:,:)
    ua=GetUnit()
    ub=GetUnit()
    ud=GetUnit()
    allocate(a(1:rho%dim,1:rho%dim))
    shift=rho%dim/2
    a=0.0_k_pr
    ci=1
    cj=1
    ! read acceptor
    open(unit=ua,file="rhoacceptor.dat",status="old", action="read")
    read(ua,*)n,m
    do i=1,n
      read(ua,*)(a(i,j),j=1,m)
    enddo
    do i=ci,ci+n/2-1
      do j=cj,cj+m/2-1
        rho%a(i,j)=cmplx(a(i,j),0.0_k_pr,k_pr)
        rho%a(i+shift,j+shift)=cmplx(a(i+n/2,j+m/2),0.0_k_pr,k_pr)
      enddo
    enddo
    close(ua)
!      do i=1,n
!       do j=1,m
!          write(777,'(x,f8.4,x)',advance="no") a(i,j)
!        enddo
!        write(777,*)
!      enddo
    ci=ci+n/2
    cj=cj+m/2
    open(unit=ub,file="rhobridge.dat",status="old", action="read")
    a=0.0_k_pr
    read(ub,*)n,m
    do i=1,n
      read(ub,*)(a(i,j),j=1,m)
    enddo
!     do i=1,n
!       do j=1,m
!          write(778,'(x,f8.4,x)',advance="no") a(i,j)
!        enddo
!        write(778,*)
!      enddo
    do i=ci,ci+n/2-1
      do j=cj,cj+m/2-1
        rho%a(i,j)=cmplx(a(i-ci+1,j-cj+1),0.0_k_pr,k_pr)
        rho%a(i+shift,j+shift)=cmplx(a(i-ci+1+n/2,j-cj+1+m/2),0.0_k_pr,k_pr)
      enddo
    enddo
    close(ub)
    ci=ci+n/2
    cj=cj+m/2
    open(unit=ud,file="rhodonor.dat",status="old", action="read")
    a=0.0_k_pr
    read(ud,*)n,m
    do i=1,n
      read(ud,*)(a(i,j),j=1,m)
    enddo
    close(ud)
!     do i=1,n
!       do j=1,m
!         write(779,'(x,f8.4,x)',advance="no") a(i,j)
!       enddo
!       write(779,*)
!     enddo
    do i=ci,ci+n/2-1
      do j=cj,cj+m/2-1
        rho%a(i,j)=cmplx(a(i-ci+1,j-cj+1),0.0_k_pr,k_pr)
        rho%a(i+shift,j+shift)=cmplx(a(i-ci+1+n/2,j-cj+1+m/2),0.0_k_pr,k_pr)
      enddo
    enddo

!      do i=1,rho%dim
!       do j=1,rho%dim
!          write(*,'(x,f8.4,x)',advance="no") real(rho%a(i,j))
!        enddo
!        write(*,*)
!      enddo
    deallocate(a)
  end subroutine GetRho

!> \brief driver routine for BFGS optimization
!> \author Alin M Elena
!> \date 10/11/07, 15:22:53
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine BFGS(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'BFGS'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    integer :: aniunit
    integer  :: i
    integer  :: m,n, res
    real(k_pr) :: epsx,epsf,epsg,gtol,xtol,ftol
    real(k_pr), allocatable :: x(:)
    integer  :: info,atom, iprint(1:2),maxfev

    aniunit=GetUnit()
    write(io%uout,'(/a/)')&
         '--BFGS Optimization Run-----------------------------------------'
    if (gen%writeAnimation) then
       open(file='bfgsDyn.xyz',unit=aniunit)
    endif
    n = atomic%atoms%nmoving*3
    m=min(gen%HessianM,n)
    if (io%verbosity > k_highVerbos) then
      iprint(1) = 1
    else
      iprint(1) = -1
    endif    
    iprint(2) = 0
    epsf = gen%epsF
    epsg = gen%epsG
    epsx = gen%epsX
    xtol = gen%xtol
    gtol=gen%ftol
    info = 0
    maxfev=gen%maxFeval
    ftol=gen%ftol    
    allocate(x(1:n))
    
    do i=1,atomic%atoms%nmoving
      atom=atomic%atoms%id(atomic%atoms%moving(i))
      x(3*(i-1)+1) = atomic%atoms%x(atom)
      x(3*(i-1)+2) = atomic%atoms%y(atom)
      x(3*(i-1)+3) = atomic%atoms%z(atom)
    enddo
    call lbfgs(n,m,x,epsg,epsf,epsx,xtol,gtol,ftol,maxfev,gen%nsteps,iprint,info,&
                 UpdatePoint,gen,atomic,tb,sol,io)
    res=UpdatePoint(gen,atomic,tb,sol,io,x,ftol,x)
    deallocate(x) 
   
       write(io%uout,'(/a)')&
          'Final forces:'
       call PrintForces(atomic%atoms,io)
       write(io%uout,'(/a)')&
          'Final coordinates:'
       call PrintCoordinates(io,atomic)
       if (info<0) then
          write(io%uout,'(/a,i4)')&
             'WARNING: Optimization did not converge.',info
       endif
       if (gen%writeAnimation) then
          close(aniunit)
       endif
! !      call write_fdf_coords()
       write(io%uout,'(/a)')&
          '----------------------------------------------------------------'
   end subroutine BFGS

  integer function UpdatePoint(gen,atomic,tb,sol,io,x,f,gradient)
    character(len=*), parameter :: myname = 'UpdatePoint'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    real(k_pr), intent(inout) :: f, gradient(:)
    real(k_pr), intent(in) :: x(:)
    integer :: i,atom
    UpdatePoint=0
    do i=1,atomic%atoms%nmoving
      atom=atomic%atoms%id(atomic%atoms%moving(i))
      atomic%atoms%x(atom) = x(3*(i-1)+1)
      atomic%atoms%y(atom) = x(3*(i-1)+2)
      atomic%atoms%z(atom) = x(3*(i-1)+3)
    enddo
    call SinglePoint(io,gen,atomic,tb,sol)
    if (.not.gen%lIsSCFConverged) then
      UpdatePoint=1
      gen%scf=.false.
      call error(myname,"change to non-SCF calculation for this point",.false.,io)
      call SinglePoint(io,gen,atomic,tb,sol)
      gen%scf=.true.
    endif
    do i=1,atomic%atoms%nmoving
      atom=atomic%atoms%id(atomic%atoms%moving(i))
      gradient(3*(i-1)+1) = -atomic%atoms%fx(atom)
      gradient(3*(i-1)+2) = -atomic%atoms%fy(atom)
      gradient(3*(i-1)+3) = -atomic%atoms%fz(atom)
    enddo
    f=sol%totalEnergy
   end function UpdatePoint
end module m_DriverRoutines
