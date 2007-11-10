!> \brief contains the routines that "drive" the program
!> \author Alin M Elena
!> \date 02/11/07, 18:04:08


module m_DriverRoutines
  use m_Constants
  use m_Useful
  use m_Types
  use m_TightBinding
  use m_Gutenberg
  use m_SCF
  use m_Hamiltonian
  use m_Dynamics
  implicit none

  private

  public :: SinglePoint
  public :: BornOppenheimerDynamics

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
    real(k_pr) :: eenergy=0.0_k_pr,renergy=0.0_k_pr, minusts,scfE=0.0_k_pr

    if (genLoc%spin)  then
      call InitMagneticMoment(atomic)
      write(ioLoc%uout,"(a)")"Initial Magnetic Moment"
      call PrintMagneticMoment(atomic,sol,.false.,ioLoc)
    endif
    call FullSCF(ioLoc,genLoc,atomic,tbMod,sol)
    eenergy = ElectronicEnergy(genLoc,sol,ioLoc)
    renergy = RepulsiveEnergy(genLoc,atomic%atoms,tbMod)
    minusts = - genLoc%electronicTemperature * sol%electronicEntropy
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
    integer  :: i,istep,m,n,k
    real(k_pr) :: dt,mi,currTime

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
    endif
    open(file='bo_dyn.ENE',unit=eneunit)
    ! initialize forces and velocities
    call FullSCF(io,gen,atomic,tb,sol)
    call InitVelocities(gen,atomic,sol)
    dt = gen%deltat
    !write out the initial vels, forces and energies
    call PrintVelocities(io,atomic)
    call PrintForces(atomic%atoms,io)
    eenergy = ElectronicEnergy(gen,sol,io)
    renergy = RepulsiveEnergy(gen,atomic%atoms,tb)
    minusts = - gen%electronicTemperature * sol%electronicEntropy
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
    do istep=1,gen%nsteps
      !set global time variable
      currTime = (istep-1)*dt
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
      ! calculate forces at t+dt
      call FullSCF(io,gen,atomic,tb,sol)
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
      minusts = - gen%electronicTemperature * sol%electronicEntropy
      if (gen%scf) then
        scfE = ScfEnergy(gen,atomic,sol,io)
      else
        scfE = 0.0_k_pr
      endif
      penergy = eenergy + renergy + scfE + minusts
      kenergy = KineticEnergy(atomic)
!         if (gen%writeAnimation) then
!             call BuildDensity(density)
!             call calc_charges(rho)
!             call calc_dipoles(density)
!             call write_frame(xunit)          
!             call writeAnimation_frame(aniunit)
! !           call write_frame_rho(runit)
!          endif

      call PrintCoordinates(io,atomic)
      call PrintForces(atomic%atoms,io)
      write(io%uout,'(i5,a,f13.6,a,f13.6,a,f13.6)')&
            istep,' P = ',penergy,' K = ',kenergy,' E = ',penergy+kenergy
      write(eneunit,'(7f29.18)')currTime,renergy,eenergy,scfE,minusts,kenergy,penergy+kenergy
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

!> \brief driver routine for verlet velocity (corrected) Ehrenfest molecular dynamics
!> \author Alin M Elena
!> \date 10/11/07, 15:22:53
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine EhrenfestDynamics(io,gen,atomic,tb,sol)
    character(len=*), parameter :: myname = 'EhrenfestDynamics'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb
    integer :: aniunit, eneunit, popunit, xunit,runit, accUnit, donUnit, spacUnit
    real(k_pr) :: eenergy,renergy,kenergy,penergy,scfenergy
    integer  :: i,istep,n,m,j
    real(k_pr) :: dt,mi
    complex(k_pr) :: ihbar,trrho,st
    type(matrix) :: rhoold,rhodot,rhonew,rho0
    real(k_pr) ::biasFactor,bfa

    if (gen%scf .and. gen%scftype=k_scftbuj) then

      write(io%uout,'(/a/)')&
            '--Velocity Verlet Ehrenfest Dynamics----------------------------'
      if (gen%writeAnimation) then
            open(file='eh_dyn.xyz',unit=aniunit)
            open(unit=xunit,file="eh_dyn.gcd",form="UNFORMATTED",status="unknown",action="write")
            open(unit=accUnit,file="eacceptor.dat",status="unknown",action="write")
            open(unit=donUnit,file="edonor.dat",status="unknown",action="write")
            open(unit=spacUnit,file="espacer.dat",status="unknown",action="write")
!             call write_header(xunit,gen%nsteps,gen%deltat,atomic%natoms)
            open(unit=runit,file="eh_dyn.rho",form="UNFORMATTED",status="unknown",action="write")
!             call write_header_rho(runit,gen%nsteps,gen%deltat,n)
         endif

         open(file='eh_dyn.ENE',unit=eneunit)
         open(file='eh_dyn.POP',unit=popunit)

    ! prepare the electronic subsystem,
    ! which is an eigenstate of the hamitonian
    ! with the bias

         call full_scf

         dt = gen%deltat

    ! initialize DM storage spaces

         call CreateMatrix(rhoold,h%dim,.true.)
         call CreateMatrix(rhodot,h%dim,.true.)
         call CreateMatrix(rhonew,h%dim,.true.)
         call CreateMatrix(rho0,h%dim,.true.)

         if (gen%writeAnimation) then
            call BuildDensity(density)
            call calc_charges(rho)
            call calc_dipoles(density)
            call write_frame(xunit) 
!        call write_frame_rho(runit,rho0) 
            call writeAnimation_frame(aniunit)
            write(accUnit,*) "0.0", charge_on_group(atomic%acceptor,atomic%nacceptor)
            write(donUnit,*) "0.0", charge_on_group(atomic%donor,atomic%ndonor)
            write(spacUnit,*) "0.0", charge_on_group(atomic%spacer,atomic%nspacer)
         endif       


    ! now build a hamiltonian with no bias
         call BuildHamiltonian
!      write(io%uout,'(/a/)')&1
!        '---- Start reading initial density'

!      call init_rho("dma.rho","dms.rho","dmd.rho",rho%a,rho%dim)          

!      write(io%uout,'(/a/)')&
!        '---- End reading initial density'        
         call copy_matrix(rho0,rho)             
         if (gen%BiasRampSteps.gt.0) call add_bias(1.0_k_pr)
         call add_hn0
         call copy_sparse_matrix(hin,h)
         call BuildDensity(density)
!          if (.not.gen%comp_elec) then
!             if (gen%electrostatics==tbu_multi) call init_qvs(density)
!          endif 
         call AddH2(density)       

         ihbar = cmplx(0.0_k_pr,-1.0_k_pr/hbar,k_pr)

    ! go back in time one step for the DM integration


         call commutator(rhodot,h,rho)
         st = cmplx(-dt,0.0_k_pr,k_pr)
         call scalar_t_matrix(ihbar*st,rhodot)
         call matrix_add(rhoold,rho,rhodot)
      ! calculate the forces from the prepared DM and the present H
         call copy_sparse_matrix(h,hin)

  !      call zero_forces
  !      call repulsive_forces
  !      call electronic_forces
    ! initialize the velocities

         call init_velocities

    ! now we are ready to start the dynamics,
    ! we have the forces, velocities, positions
    ! and rho at time t

         st = cmplx(2.0_k_pr*dt,0.0_k_pr,k_pr)
         do istep=1,gen%nsteps

       !set global time variable
            currTime = (istep-1)*dt

            call BuildDensity(density)
!             if (.not.gen%comp_elec) then
!                if (gen%electrostatics==tbu_multi) call init_qvs(density)
!             endif
            call AddH2(density)

            call zero_matrix(rhodot)
            call commutator(rhodot,h,rho)
            call scalar_t_matrix(ihbar*st,rhodot)       
            call matrix_add(rhonew,rhoold,rhodot)
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
            call copy_matrix(rhoold,rho)
            call copy_matrix(rho,rhonew)
           ! calculate forces at t+dt
            call BuildHamiltonian
       ! Ramp for the bias
            if ((gen%BiasRampSteps>0).and.((istep)<=gen%BiasRampSteps)) then
               bfa = real(istep,k_pr)/real(gen%BiasRampSteps,k_pr)
               biasFactor = -(bfa-1)**3 * (1 + 3*bfa + 6*bfa**2)
               call AddBias(biasFactor)
            endif

            call BuildDensity(density)
!             if (.not.gen%comp_elec) then
!                if (gen%electrostatics==tbu_multi) call init_qvs(density)
!             endif

            call zero_forces
            call repulsive_forces
            call electronic_forces
            call scf_forces

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
            call scale_velocities

            if (gen%spin) then
               trrho   = matrix_trace(rho)
            else
               trrho   = 2.0_k_pr * matrix_trace(rho)
            endif
            eenergy = electronic_energy()
            renergy = repulsive_energy()
            if (gen%scf) then
               scfenergy = scf_energy()
            else
               scfenergy = 0.0_k_pr
            endif

            penergy = eenergy + renergy + scfenergy
            kenergy = kinetic_energy()

            if (gen%writeAnimation) then
               call calc_charges(rho)
               call calc_dipoles(density)
               call writeAnimation_frame(aniunit)
               call write_frame(xunit)
               call write_frame_rho(runit,rho0)
               write(accUnit,*) istep*dt, charge_on_group(atomic%acceptor,atomic%nacceptor)
               write(donUnit,*) istep*dt, charge_on_group(atomic%donor,atomic%ndonor)
               write(spacUnit,*) istep*dt, charge_on_group(atomic%spacer,atomic%nspacer)
            endif
            call write_currents(h,rho)
            if (mod(istep,50).eq.1) call write_rho_eigenvalues(rho)
            write(eneunit,'(7f25.18)')istep*dt,renergy,eenergy,scfenergy,kenergy,penergy+kenergy,real(trrho)

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

         write(io%uout,'(/a/)')&
            'End Velocity Verlet-------------------------------------------------------------'











           
      else
!          if (gen%scf) allocate(density(basis_var%norbital))
! 
!          write(io%uout,'(/a/)')&
!             '--Velocity Verlet Ehrenfest Dynamics----------------------------'
! 
!          if (gen%writeAnimation) then
!             open(file='eh_dyn.xyz',unit=aniunit)
!          endif
! 
!          open(file='eh_dyn.ENE',unit=eneunit)
!          open(file='eh_dyn.POP',unit=popunit)
! 
!     ! prepare the electronic subsystem,
!     ! which is an eigenstate of the hamitonian
!     ! with the bias
! 
!          call full_scf
! 
!          dt = gen%deltat
! 
!     ! initialize DM storage spaces
! 
!          call CreateMatrix(rhoold,h%dim,.true.)
!          call CreateMatrix(rhodot,h%dim,.true.)
!          call CreateMatrix(rhonew,h%dim,.true.)
! 
!     ! now build a hamiltonian with no bias
! 
!          call BuildHamiltonian
!          if (gen%BiasRampSteps.gt.0) call add_bias(1.0_k_pr)
!          if (gen%scf) then
!             call build_density(density)
!             call AddH2(density)
!          endif
! 
!          ihbar = cmplx(0.0_k_pr,-1.0_k_pr/hbar,k_pr)
! 
!     ! go back in time one step for the DM integration
! 
!          call commutator(rhodot,h,rho)
!          st = cmplx(-dt,0.0_k_pr,k_pr)
!          call scalar_t_matrix(ihbar*st,rhodot)
!          call matrix_add(rhoold,rho,rhodot)
! 
!     ! calculate the forces from the prepared DM and the present H
! 
!          call clear_h_diag
!          if (gen%BiasRampSteps.gt.0) call add_bias(1.0_k_pr)
!          call zero_forces
!          call repulsive_forces
!          call electronic_forces
! 
!     ! initialize the velocities
! 
!          call init_velocities
! 
!     ! now we are ready to start the dynamics,
!     ! we have the forces, velocities, positions
!     ! and rho at time t
! 
!          st = cmplx(2.0_k_pr*dt,0.0_k_pr,k_pr)
! 
!          do istep=1,gen%nsteps
! 
!        !set global time variable
!             gen%time = (istep-1)*dt
! 
!             if (gen%scf) then
!                call build_density(density)
!                call AddH2(density)
!             endif
! 
!             call zero_matrix(rhodot)
!             call commutator(rhodot,h,rho)
!             call scalar_t_matrix(ihbar*st,rhodot) 
!             call matrix_add(rhonew,rhoold,rhodot)
! 
!        ! at this point rho contains the rho at time=t
!        ! propagate the positions 
! 
!        ! calculates positions at t+dt
!             do i=1,atomic%natoms
!                mi = species%mass(atomic%sp(i))
!                atomic%x(i) = atomic%x(i) &
!                   + dt * atomic%vx(i) &
!                   + 0.5_k_pr * dt*dt * atomic%fx(i) / mi
!                atomic%y(i) = atomic%y(i) &
!                   + dt * atomic%vy(i) &
!                   + 0.5_k_pr * dt*dt * atomic%fy(i) / mi
!                atomic%z(i) = atomic%z(i) &
!                   + dt * atomic%vz(i) &
!                   + 0.5_k_pr * dt*dt * atomic%fz(i) / mi
!             enddo
! 
!        ! store forces at t in fold
!             do i=1,atomic%natoms
!                atomic%fxo(i) = atomic%fx(i)
!                atomic%fyo(i) = atomic%fy(i)
!                atomic%fzo(i) = atomic%fz(i)
!             enddo
! 
!        ! shuffle the DMs, rho is now rho(t+dt)
!             call copy_matrix(rhoold,rho)
!             call copy_matrix(rho,rhonew)
! 
!        ! calculate forces at t+dt
!             call BuildHamiltonian
!        ! Ramp for the bias
!             if ((gen%BiasRampSteps.gt.0).and.((istep).le.gen%BiasRampSteps)) then
!                bfa = real(istep,k_pr)/real(gen%BiasRampSteps,k_pr)
!                biasFactor = -(bfa-1)**3 * (1 + 3*bfa + 6*bfa**2)
!                call add_bias(biasFactor)
!             endif
!             call zero_forces
!             call repulsive_forces
!             call electronic_forces
!             call scf_forces
! 
!        ! calculate velocities at t+dt
!             do i=1,atomic%natoms
!                mi = species%mass(atomic%sp(i))
!                atomic%vx(i) = atomic%vx(i) &
!                   + 0.5_k_pr * dt * (atomic%fx(i) + atomic%fxo(i)) / mi
!                atomic%vy(i) = atomic%vy(i) &
!                   + 0.5_k_pr * dt * (atomic%fy(i) + atomic%fyo(i)) / mi
!                atomic%vz(i) = atomic%vz(i) &
!                   + 0.5_k_pr * dt * (atomic%fz(i) + atomic%fzo(i)) / mi
!             enddo
! 
!        ! scale velocities
!             call scale_velocities
! 
!             if (gen%spin) then
!                trrho   = matrix_trace(rho)
!             else
!                trrho   = 2.0_k_pr * matrix_trace(rho)
!             endif
!             eenergy = electronic_energy()
!             renergy = repulsive_energy()
!             if (gen%scf) then
!                scfenergy = scf_energy()
!             else
!                scfenergy = 0.0_k_pr
!             endif
! 
!             penergy = eenergy + renergy + scfenergy
! 
!             kenergy = kinetic_energy()
! 
!             if (gen%writeAnimation) call writeAnimation_frame(aniunit)
! 
!        !if (mod(istep,50).eq.1) then
!             write(popunit,'(a/a/)')' '
!             do i=1,h%dim
!                write(popunit,'(i10,f13.6)')i,real(rho%a(i,i))       
!             enddo
!        !endif
! 
!        !please remove later
!        !if (mod(istep,50).eq.1) then
!        !   write(popunit,'(a/a/)')' '
!        !   do i=1,h%dim
!        !      do j=1,h%dim
!        !         write(345,*)i,j,real(rho%a(i,j)),aimag(rho%a(i,j))
!        !      enddo
!        !   enddo
!        !endif
!        !end please remove later
! 
!             call write_currents(h,rho)
!             if (mod(istep,50).eq.1) call write_rho_eigenvalues(rho)
!             write(eneunit,'(7f25.18)')istep*dt,renergy,eenergy,scfenergy,kenergy,penergy+kenergy,real(trrho)
! 
!          enddo !istep loop
! 
!          if (gen%writeAnimation) then
!             close(aniunit)
!          endif
!          close(eneunit)
!          close(popunit)
! 
!          write(io%uout,'(/a/)')&
!             '----------------------------------------------------------------'
! 
!          if (gen%scf) deallocate(density)

      endif
   end subroutine EhrenfestDynamics
      


end module m_DriverRoutines