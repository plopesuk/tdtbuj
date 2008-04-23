!> \brief self consistent field method
!> \author Alin M Elena
!> \date 05/11/07, 10:30:06
module m_SCF
  use m_Constants
  use m_Types
  use m_Useful
  use m_Gutenberg
  use m_Hamiltonian
  use m_LinearAlgebra, only : MatrixTrace,CopyMatrix,SpmPut,MatrixCeaApbB
  use m_Electrostatics
  use m_TightBinding
  use m_DensityMatrix
  use m_Mixing
  implicit none
  private

  public :: FullScf
  public :: ScfEnergy
  public :: ScfForces
  public :: AddH2

contains

!> \brief controls the scf calculation.
!> \details practically computes a single point calculation (energy and forces)
!> if you asked for a SCF calculation the scf path is followed
!> \author Alin M Elena
!> \date 05/11/07, 10:33:06
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine FullScf(ioLoc,genLoc,atomic,tbMod,sol)
    character(len=*), parameter :: sMyName="FullScf"
    type(ioType), intent(inout) :: ioLoc
    type(generalType), intent(inout) :: genLoc
    type(atomicxType), intent(inout) :: atomic
    type(modelType), intent(inout) :: tbMod
    type(solutionType), intent(inout) :: sol
    integer :: nit,m,n
    real(k_pr)  :: residual, dmax
    real(k_pr)  :: ee,re,scfe,te
    logical   :: exists,first
    integer :: ierr   
    integer :: l,ml, i,j,nmix
    complex(k_pr) :: trace
    character(len=k_ml) :: saux
    character(len=k_mw) :: labels(1:4),labelsH2(1:4)
    dmax=0.0_k_pr
    genLoc%lIsSCFConverged=.true.
    n=atomic%basis%norbitals
    m=(n-1)*n/2 
    if (genLoc%scf) then
! delta density matrix is stored in an array as upper triangular part followed by the diagonal
      select case(genLoc%scfType)
        case (k_scfTbuj)
          if (.not.genLoc%spin) then
            call error("This model should be spin polarised!",smyname,.true.,ioLoc)
          endif
!           if((.not.genLoc%compElec).and.(genLoc%k_electrostatics==k_electrostaticsMultipoles)) then
!             call allocate_qvs
!           endif
          sol%buff%dins=0.0_k_pr
          sol%buff%douts=0.0_k_pr
          sol%buff%res=0.0_k_pr

          sol%buff%densityin=0.0_k_pr
          sol%buff%densityout=0.0_k_pr
          sol%buff%densitynext=0.0_k_pr

          call BuildHamiltonian(ioLoc,genLoc,atomic,tbMod,sol)
          call AddBias(1.0_k_pr,atomic,sol)
          call CopyMatrix(sol%hin,sol%h,ioLoc)
          call DiagHamiltonian(ioLoc,genLoc,atomic,sol)
          if (ioLoc%Verbosity >= k_highVerbos) then
            write(ioLoc%uout,'(a)') "Before entering the SCF LOOP"
!
!               if (.not.genLoc%compElec) then
!                   if (genLoc%k_electrostatics==tbu_multi)call init_qvs(densityin)
!               endif
            labels(1)="Hin: Spin DD"
            labels(2)="Hin: Spin UU"
            labels(3)="Hin: Spin DU"
            labels(4)="Hin: Spin UD"
            call PrintMatrixBlocks(sol%h,labels,ioLoc,.false.,.not.genLoc%collinear)
            labels(1)="Eigenvectors: Spin DD"
            labels(2)="Eigenvectors: Spin UU"
            labels(3)="Eigenvectors: Spin DU"
            labels(4)="Eigenvectors: Spin UD"
            call PrintMatrixBlocks(sol%eigenvecs,labels,ioLoc,.false.,.not.genLoc%collinear)
            labels(1)="Hin: Spin DD"
            labels(2)="Hin: Spin UU"
            labels(3)="Hin: Spin DU"
            labels(4)="Hin: Spin UD"
            labelsH2(1)="H2: Spin DD"
            labelsH2(2)="H2: Spin UU"
            labelsH2(3)="H2: Spin DU"
            labelsH2(4)="H2: Spin UD"
            call PrintOccupationNumbers(genLoc,sol,ioLoc)
            write(ioLoc%uout,"(a,f16.8)") "chemical potential: ",genLoc%electronicMu
            write(ioLoc%uout,'(a,f16.8)')"Entropy term: ",sol%electronicEntropy
          endif
!           if (genLoc%alter_dm) then
!             call create_dm_spin_altered(eigenvec,eigenval)
!           endif
!           ! this one in fact builds the initial guess for density matrix
!           call BuildDensity(atomic,sol,genLoc,.true.)
          call BuildDensity(atomic,sol)
          sol%buff%densityin=sol%density

!           scfe=ScfEnergy(genLoc,atomic,sol,ioLoc)
          call CalcExcessCharges(genLoc,atomic,sol)
          call CalcDipoles(genLoc,atomic,sol)
          call ComputeMagneticMoment(genLoc,atomic,sol,ioLoc)

          if (ioLoc%Verbosity >= k_highVerbos) then
!
!               if (.not.genLoc%compElec) then
!                   if (genLoc%k_electrostatics==tbu_multi)call init_qvs(densityin)
!               endif
            do i=1,atomic%atoms%natoms
!                   if (genLoc%k_electrostatics==tbu_multi) then
!                     call Print_QlmR(i,densityin)
!                     call Print_VlmR(i,densityin)
!                   endif
!                   call Print_atom_rho(i,"density in ")
!                   call Print_atom_density(i,densityin,"density in ")
!
!                   if (genLoc%k_electrostatics==tbu_multi) then
!                     do j=1,atomic%natoms
!                         if (i/=j) then
!                           call Print_irregular_solid(i,j)
!                           call Print_BllpR(i,j)
!                         endif
!                     enddo
!                   endif
                call PrintAtomChargeAnalysis(i,atomic,sol,genLoc,ioLoc)
                call PrintAtomMatrix(i,atomic,sol%hin,"Hin",ioLoc,.false.)
                call PrintAtomMatrix(i,atomic,sol%rho,"Density",ioLoc,.false.)
              enddo
              call PrintCharges(genLoc,atomic,ioLoc)
              call PrintDipoles(atomic,ioLoc)
              call PrintMagneticMoment(atomic,sol,.false.,ioLoc)
            endif
!
            if (ioLoc%Verbosity >= k_highVerbos) then
              write(ioLoc%uout,*) 'SCF'
              write(ioLoc%uout,*) &
                  'nit          energy          res        drmax      Tr[rho]           mu'
            endif
            do nit=1,genLoc%maxscf
              if (ioLoc%Verbosity >=k_highVerbos)  then
                   write(ioLoc%uout,'(a,i0)') "SCF LOOP iteration: ",nit
              endif
!               if (.not.genLoc%compElec) then
!                   if (genLoc%k_electrostatics==tbu_multi) call init_qvs(densityin)
!               endif
              call AddH2(genLoc,atomic,sol,tbMod,ioLoc)
              call DiagHamiltonian(ioLoc,genLoc,atomic,sol)
!               if (genLoc%alter_dm) then
!                   call create_dm_spin_altered(eigenvec,eigenval)
!               endif
              call BuildDensity(atomic,sol)

              sol%buff%densityout=sol%density
              call CalcExcessCharges(genLoc,atomic,sol)
              call CalcDipoles(genLoc,atomic,sol)
              call ComputeMagneticMoment(genLoc,atomic,sol,ioLoc)
              if (ioLoc%Verbosity >= k_highVerbos) then
!                 if (.not.genLoc%compElec) then
!                   if (genLoc%k_electrostatics==tbu_multi) call init_qvs(densityout)
!                 endif
                call MatrixCeaApbB(sol%h2,sol%h,sol%hin,k_cOne,-k_cOne,ioLoc)
                call PrintMatrixBlocks(sol%h,labels,ioLoc,.false.,.not.genLoc%collinear)
                call PrintMatrixBlocks(sol%h2,labelsH2,ioLoc,.false.,.not.genLoc%collinear)
                call PrintOccupationNumbers(genLoc,sol,ioLoc)
                write(ioLoc%uout,"(a,f16.8)") "chemical potential: ",genLoc%electronicMu
                write(ioLoc%uout,'(a,f16.8)')"Entropy term: ",sol%electronicEntropy
                do i=1,atomic%atoms%natoms
!                     if (genLoc%k_electrostatics==tbu_multi) then
!                       call Print_QlmR(i,densityout)
!                       call Print_VlmR(i,densityout)
!                     endif
                  call PrintAtomChargeAnalysis(i,atomic,sol,genLoc,ioLoc)
                  call PrintAtomMatrix(i,atomic,sol%h,"H",ioLoc,.false.)
                  call PrintAtomMatrix(i,atomic,sol%hin,"Hin",ioLoc,.false.)
                  call PrintAtomMatrix(i,atomic,sol%h2,"H2",ioLoc,.false.)
                  call PrintAtomMatrix(i,atomic,sol%rho,"Density",ioLoc,.false.)
                enddo
                call PrintCharges(genLoc,atomic,ioLoc)
                call PrintDipoles(atomic,ioLoc)
                call PrintMagneticMoment(atomic,sol,.true.,ioLoc)
              endif
!
              ierr=-1
              nmix=genLoc%scfMixn
              do while ((genLoc%scfMixn>=1).and.(ierr/=0))              
               call InitMix(sol%buff%dins,sol%buff%douts,sol%buff%res,sol%buff%densityin,sol%buff%densityout,genLoc%scfMixn)
               call MixDensity(sol%buff%dins,sol%buff%douts,sol%buff%res,sol%buff%densitynext,residual,dmax,genLoc%scfMix,n+m,genLoc%scfMixn,nit,ierr)              
                if (ierr/=0) then
                  genLoc%scfMixn=genLoc%scfMixn-1
                  call error("Singularity in mixing matrix, no of iterations mixed reduced by one ",smyname,.false.,ioLoc)
                endif
              enddo
              if (ierr /=0) then
                call error("Singularity in mixing matrix, iterations to mix reduced up to 2",smyname,.true.,ioLoc)
              endif
              genLoc%scfMixn=nmix
              sol%buff%densityin = sol%buff%densitynext
              sol%density = sol%buff%densitynext
              call CopyMatrix(sol%h,sol%hin,ioLoc)
!               if (ioLoc%Verbosity >=k_highVerbos) then
!                   call Print_density(densityin,"density in")
!                   call Print_density(densityout,"density out")
!                   call Print_density(densitynext,"density next")
!               endif

              if (ioLoc%Verbosity >= k_highVerbos) then
                call ZeroForces(atomic)
                call RepulsiveForces(genLoc,atomic%atoms,tbMod)
                call electronicForces(atomic,genLoc,tbMod,sol,ioLoc)
                call ScfForces(genLoc,atomic,sol,ioLoc)
                sol%density=sol%buff%densitynext
                call PrintForces(atomic%atoms,ioLoc)
                ee = ElectronicEnergy(genLoc,sol,ioLoc)
                re = RepulsiveEnergy(genLoc,atomic%atoms,tbMod)
                write(ioLoc%uout,'(a)')"Energy"
                write(ioLoc%uout,'(a,f16.8)')"Electronic: ",ee
                write(ioLoc%uout,'(a,f16.8)')" Repulsive: ",re
                scfe = ScfEnergy(genLoc,atomic,sol,ioLoc)
                sol%density=sol%buff%densitynext
                write(ioLoc%uout,'(a,f16.8)')" -TS: ",sol%electronicEntropy
                write(ioLoc%uout,'(a,f16.8)')"       SCF: ",scfe
                write(ioLoc%uout,'(a,f16.8)')"     Total: ",ee+re+scfe
              endif
              if (dmax < genLoc%scftol) exit
            end do
!
            if ( ioLoc%uout /= 6) then
              if (nit>genLoc%maxscf) then
                write(6,'(a,i0,a,ES12.4)')"Warning: it did not converge after ",genLoc%maxscf," the tolerance reached is ",dmax
                genLoc%lIsSCFConverged = .false.
              else
                write(6,'(a,i0,a,ES12.4)') "converged in ",nit, " iterations up to ",dmax
              endif
            endif
            if (nit>genLoc%maxscf) then
              write(ioLoc%uout,'(a,i0,a,ES12.4)')"Warning: it did not converge after ",nit," the tolerance reached is ",dmax
              genLoc%lIsSCFConverged = .false.
            else
               write(ioLoc%uout,'(a,i0,a,ES12.4)') "converged in ",nit, " iterations up to ",dmax
            endif
            
            call CalcExcessCharges(genLoc,atomic,sol)
            call CalcDipoles(genLoc,atomic,sol)
            call ComputeMagneticMoment(genLoc,atomic,sol,ioLoc)
!           ! calculate the forces
!             if (.not.genLoc%compElec) then
!               if (genLoc%k_electrostatics==tbu_multi) call init_qvs(densityin)
!             endif

            call ZeroForces(atomic)
            call RepulsiveForces(genLoc,atomic%atoms,tbMod)
            call ElectronicForces(atomic,genLoc,tbMod,sol,ioLoc)
            call ScfForces(genLoc,atomic,sol,ioLoc)                        
!             call Print_eigens(eigenval,eigenvec,555,.false.)
!             call write_density_matrix("rho.bin",rho%a,rho%dim)
            if (ioLoc%verbosity >= k_mediumVerbos) then
              trace = MatrixTrace(sol%rho,ioLoc)
              write(saux,'(a,"(",f0.4,1x,f0.4,"i)")') "Density matrix, Trace= ",trace
             call PrintMatrix(sol%rho,trim(saux),ioLoc)
            endif            
        end select
      else
        call BuildHamiltonian(ioLoc,genLoc,atomic,tbMod,sol)
        call AddBias(1.0_k_pr,atomic,sol)
        call DiagHamiltonian(ioLoc,genLoc,atomic,sol)
        call BuildDensity(atomic,sol)
        if (ioLoc%verbosity >= k_mediumVerbos) then
          call PrintMatrix(sol%h,"Hamiltonian Matrix:",ioLoc)
          call PrintVectorA(sol%eigenvals,"Eigenvalues",.false.,.true.,ioLoc)
          call PrintMatrix(sol%eigenvecs,"Eigenvectors",ioLoc)
          trace = MatrixTrace(sol%rho,ioLoc)
          write(saux,'(a,"(",f0.4,1x,f0.4,"i)")')"Density matrix, Trace= ",trace
          call PrintMatrix(sol%rho,trim(saux),ioLoc)
        endif
        call ZeroForces(atomic)
        call RepulsiveForces(genLoc,atomic%atoms,tbMod)
        call electronicForces(atomic,genLoc,tbMod,sol,ioLoc)
        call CalcExcessCharges(genLoc,atomic,sol)
        call CalcDipoles(genLoc,atomic,sol)
        call ComputeMagneticMoment(genLoc,atomic,sol,ioLoc)
!         call Print_eigens(eigenval,eigenvec,555,.false.)
      end if
  end subroutine FullScf


!> \brief adds the \f$ {\mathbf H}_2 \f$ to the Hamiltonian
!> \author Alin M Elena
!> \date 08/11/07, 14:02:54
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param tb type(modelType) contains information about the tight binding model parameters
  subroutine AddH2(gen,atomic,sol,tb,io)
    character(len=*), parameter :: myname = 'AddH2'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    type(modelType), intent(inout) :: tb

    integer       :: i,k,j
    integer ::l1,l2,l3,l4,m1,m2,m3,m4,o1,o2,o3,o4,sp,shift
    real(k_pr) :: rAddAcc,rTmp,hij,hijf,hijd
    integer :: m,n
    real(k_pr) :: q0,q0up,q0down,aux,udq,elec

    n=sol%h%dim
    m=n*(n-1)/2
    shift=n/2
    q0=0.0_k_pr
    q0up=0.0_k_pr
    q0down=0.0_k_pr 
    elec=k_e2/(4.0_k_pr*k_pi*k_epsilon0)
    select case(gen%scfType)
      case(k_scfTbuj)
        rTmp=0.0_k_pr
        rAddAcc=0.0_k_pr
        do k=1,atomic%atoms%nscf
          i=atomic%atoms%scf(k)
          call ScfChargeNumbers(i,q0,q0up,q0down,atomic,sol)
!! ! spin down
          udq=atomic%species%ulocal(atomic%atoms%sp(i))*q0
          rAddAcc=-atomic%species%jlocal(atomic%atoms%sp(i))*q0down! spin up
          rTmp=-atomic%species%jlocal(atomic%atoms%sp(i))*q0up
          if (io%Verbosity >= k_highVerbos) then
            write(io%uout,'(a,i5,a,i5,a,a2)')" Atom: ", i,"  specie ",atomic%atoms%sp(i),&
                " element ",symbol(atomic%species%z(atomic%atoms%sp(i)))
            write(io%uout,'(a,f16.8)')" U*dQ: ", udq
            write(io%uout,'(a,f16.8,a,f16.8)')"J*dNd: ", rAddAcc, "  J*dNu",rTmp
          endif
          do o1=1,atomic%species%norbs(atomic%atoms%sp(i))/2
            o2=atomic%atoms%orbs(i,o1)
            o3=o2+shift
!! !spin down 
            hij=sol%h%a(o2,o2)+rAddAcc+udq
	    if (abs(hij)>=gen%hElementThreshold) then
		call SpmPut(sol%h,o2,o2,cmplx(hij,0.0_k_pr,k_pr))
	    endif
!! !spin up
             hij=sol%h%a(o3,o3)+rTmp+udq             
	     if (abs(hij)>=gen%hElementThreshold) then
                call SpmPut(sol%h,o3,o3,cmplx(hij,0.0_k_pr,k_pr))            
	     endif		
          enddo
        enddo
    end select
    select case(gen%electrostatics)
      case (k_electrostaticsPoint)
        call BuildPotential(gen,atomic,sol)
        do k=1,atomic%atoms%nscf
          i=atomic%atoms%scf(k)
          if (io%Verbosity >= k_highVerbos) then
            write(io%uout,'(a,i5,a,i5,a,a2)')" Atom: ", i,"  specie ",atomic%atoms%sp(i),&
             " element ",symbol(atomic%species%z(atomic%atoms%sp(i)))
            write(io%uout,'(a,f16.8)')" Vrr': ", sol%potential(i)
          endif
          do o1=1,atomic%species%norbs(atomic%atoms%sp(i))
            hij=sol%h%a(atomic%atoms%orbs(i,o1),atomic%atoms%orbs(i,o1))+sol%potential(i)
            call SpmPut(sol%h,atomic%atoms%orbs(i,o1),atomic%atoms%orbs(i,o1),cmplx(hij,0.0_k_pr,k_pr))
          enddo
        enddo
      case(k_electrostaticsMultipoles)
!             aux=0.0_k_pr
!              write(ioLoc%uout,"(a)")"Multipoles increments"
!             do k=1,atomic%nscf
!                i=atomic%isscf(k)
!                sp=atomic%sp(i)
!        ! spin down
!                do o1=1,species%norbs(sp)/2
!                   l1=atomic%basis%orbitals(atomic%orbs(i,o1))%l
!                   m1=atomic%basis%orbitals(atomic%orbs(i,o1))%m
!                   do o2=1,species%norbs(sp)/2
!                      l2=atomic%basis%orbitals(atomic%orbs(i,o2))%l
!                      m2=atomic%basis%orbitals(atomic%orbs(i,o2))%m
!                      aux=hiujv(i,l1,m1,l2,m2,density)
!                       if (ioLoc%Verbosity >= k_highVerbos) then
!                       if (abs(aux)>epsilon(aux)) &
!                       write(ioLoc%uout,"(a,i0,x,i0,x,i0,f12.8)")"d: ",i,o1,o2,aux
!                       endif
!                      hij=h%a(atomic%orbs(i,o1),atomic%orbs(i,o2))+aux
!                      call spm_put(h,atomic%orbs(i,o1),atomic%orbs(i,o2),cmplx(hij,0.0_k_pr,dp))
! !           f%a(atomic%orbs(i,o1),atomic%orbs(i,o2))=f%a(atomic%orbs(i,o1),atomic%orbs(i,o2))+cmplx(aux,0.0_k_pr,dp)
!                   enddo
!                enddo
! !              spin up
!                do o1=1+species%norbs(sp)/2,species%norbs(sp)
!                   l1=atomic%basis%orbitals(atomic%orbs(i,o1))%l
!                   m1=atomic%basis%orbitals(atomic%orbs(i,o1))%m
!                   do o2=1+species%norbs(sp)/2,species%norbs(sp)
!                      l2=atomic%basis%orbitals(atomic%orbs(i,o2))%l
!                      m2=atomic%basis%orbitals(atomic%orbs(i,o2))%m
!                      aux=hiujv(i,l1,m1,l2,m2,density)
!                      if (ioLoc%Verbosity >= k_highVerbos) then
!                      if (abs(aux)>epsilon(aux)) &
!                       write(ioLoc%uout,"(a,i0,x,i0,x,i0,f12.8)")"u: ",i,o1,o2,aux
!                       endif
!                      hij=h%a(atomic%orbs(i,o1),atomic%orbs(i,o2))+aux
!                      call spm_put(h,atomic%orbs(i,o1),atomic%orbs(i,o2),cmplx(hij,0.0_k_pr,dp))
!                   enddo
!                enddo
!             enddo
    end select    
  end subroutine AddH2

!> \brief computes the electron numbers for an atom
!> \details the toal number and in each spin channel
!> \author Alin M Elena
!> \date 08/11/07, 15:31:01
!> \param at integer the atom
!> \param q0,q0up,q0down reals electron numbers total, spin up and spin down
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine ScfChargeNumbers(at,q0,q0up,q0down,atomic,sol)
    character(len=*), parameter :: myname="ScfChargeNumbers"
    real(k_pr), intent(inout) :: q0,q0up,q0down
    integer, intent(in) :: at
    type(solutionType), intent(in) :: sol
    type(atomicxType), intent(in) :: atomic
    integer :: from,to,m

      q0up=0.0_k_pr
      q0down=0.0_k_pr
! spin down
      
      m=atomic%basis%norbitals*(atomic%basis%norbitals-1)/2
      from=m+atomic%atoms%orbs(at,1)
      to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2
      q0down=sum(sol%density(from:to))
!spin up
      from=m+atomic%atoms%orbs(at,1)+atomic%basis%norbitals/2
      to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2
      q0up=sum(sol%density(from:to))      
! total charge
      q0=q0down+q0up
   end subroutine ScfChargeNumbers

!> \brief computes the scf energy
!> \author Alin M Elena
!> \date 08/11/07, 23:13:45
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  real(k_pr) function ScfEnergy(gen,atomic,sol,io)
    character(len=*), parameter :: myname = 'ScfEnergy'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    integer       :: i,k,mp,j
    real(k_pr)      :: scfe,scfx!,el_en
! +U variables
    integer ::l1,l2,l3,l4,m1,m2,m3,m4,o1,o2,o3,o4,sp,shift,m,n
    real(k_pr) :: v_tmp,elecEn,aux
    real(k_pr) :: q0,q0up,q0down
  !-------------------------------------------------!

    scfe = 0.0_k_pr
    elecEn=0.0_k_pr
    sol%density=0.0_k_pr
    call BuildDensity(atomic,sol)
    select case(gen%scfType)
      case (k_scfTbuj)
        scfe=0.0_k_pr
        scfx=0.0_k_pr
        do k=1,atomic%atoms%nscf
          i=atomic%atoms%scf(k)
          call ScfChargeNumbers(i,q0,q0up,q0down,atomic,sol)
          scfe=scfe+atomic%species%ulocal(atomic%atoms%sp(i))*q0*q0
          scfx=scfx-atomic%species%jlocal(atomic%atoms%sp(i))*(q0up*q0up+q0down*q0down)
        enddo
        scfe=scfe*0.5_k_pr!*k_e2/(4.0_k_pr*k_pi*k_epsilon0)
        scfx=scfx*0.5_k_pr!*k_e2/(4.0_k_pr*k_pi*k_epsilon0)
        if (io%verbosity >= k_highVerbos) then
          write(io%uout,'(a,f16.8)')" SCF energy from U: ", scfe
          write(io%uout,'(a,f16.8)')" SCF energy from J: ", scfx
        endif
        scfe=scfe+scfx
    end select

    select case(gen%electrostatics)
      case (k_electrostaticsPoint)
      call BuildPotential(gen,atomic,sol)
        do k=1,atomic%atoms%nscf
          i=atomic%atoms%scf(k)
          elecEn=elecEn+0.5_k_pr*charge(i,gen,atomic,sol)*sol%potential(i)
        enddo
      case(k_electrostaticsMultipoles)
!             do k=1,atomic%nscf
!                i=atomic%isscf(k)
!                aux=0.0_k_pr
!                do l1=0,2*get_lmax(i)
!                   do m1=-l1,l1
!                      aux=aux+qlmR(i,l1,m1,density)*vlmR(i,l1,m1,density)
!                   enddo
!                enddo
!                elecEn=elecEn+0.5_k_pr*aux*e2/(4.0_k_pr*pi*epsilon0)
!             end do
    end select
    if (io%verbosity >= k_highVerbos) then
      write(io%uout,'(a,f16.8)')"SCF Electrostatics: ", elecEn
    endif
    ScfEnergy = scfe+elecEn
  end function ScfEnergy
!> \brief computes the scf contribution to the forces
!> \author Alin M Elena
!> \date 08/11/07, 23:14:30
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine ScfForces(gen,atomic,sol,io)
    character(len=*), parameter :: myname = 'ScfForces'
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    integer       :: i,m,n,li,j,mi!,j,mpp,l,m
    real(k_pr) :: q
    integer :: lm,k
    real(k_pr) :: aux,sx,sy,sz

    select case(gen%electrostatics)
      case (k_electrostaticsPoint)
        sol%density=0.0_k_pr
        call BuildDensity(atomic,sol)
        call BuildField(gen,atomic,sol)
        do k=1,atomic%atoms%nmoving
          i=atomic%atoms%moving(k)
          q=charge(i,gen,atomic,sol)
          atomic%atoms%fx(i) = atomic%atoms%fx(i)+q*sol%field(i,1)
          atomic%atoms%fy(i) = atomic%atoms%fy(i)+q*sol%field(i,2)
          atomic%atoms%fz(i) = atomic%atoms%fz(i)+q*sol%field(i,3)
        end do
      case(k_electrostaticsMultipoles)
!         n=basis_var%norbital
!         m=n*(n-1)/2
!         allocate(density(m+n))
!         call build_u_density(density)
!         do k=1,atomic%nmoving
!           i=atomic%moving(k)
!           if (control_var%output_level > ol_verbose) &
!               write(control_var%output_file,'(a,i0)')"electrostatic forces contribution by l atom ",i
!           do li=0,2*get_lmax(i)
!               aux=(2.0_dp*li+3.0_dp)*sqrt(4.0_dp*pi/3.0_dp)*e2/(4.0_dp*pi*epsilon0)
!               sx=0.0_dp
!               sy=0.0_dp
!               sz=0.0_dp
!               do mi=-li,li
!                 sx=sx+fip(i,li,mi,1,density)*qlmR(i,li,mi,density)
!                 sy=sy+fip(i,li,mi,-1,density)*qlmR(i,li,mi,density)
!                 sz=sz+fip(i,li,mi,0,density)*qlmR(i,li,mi,density)
!               end do
!               atomic%fx(i) = atomic%fx(i)-sx*aux
!               atomic%fy(i) = atomic%fy(i)-sy*aux
!               atomic%fz(i) = atomic%fz(i)-sz*aux
!               if (control_var%output_level > ol_verbose) then
!                 write(control_var%output_file,'(i0,a,3f16.8,a)')&
!                     li," (",-sx*aux,-sy*aux,-sz*aux,")"
!               endif
!           end do
!         end do
    end select
  end subroutine ScfForces

end module m_SCF