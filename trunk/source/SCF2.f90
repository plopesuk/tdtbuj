!> \brief self consistent field method
!> \author Alin M Elena
!> \date 05/11/07, 10:30:06
module m_SCF2
  use m_Constants
  use m_Types
  use m_Useful
  use m_LinearAlgebra, only : SpmPut
  use m_Electrostatics
  
  
    
  implicit none
  private

  public :: AddH2

contains

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
          udq=atomic%species%ulocal(atomic%atoms%sp(i))*q0*elec
          rAddAcc=-atomic%species%jlocal(atomic%atoms%sp(i))*q0down*elec! spin up
          rTmp=-atomic%species%jlocal(atomic%atoms%sp(i))*q0up*elec
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
    type(solutionType), intent(inout) :: sol
    type(atomicxType), intent(inout) :: atomic
    integer :: from,to,m

      q0up=0.0_k_pr
      q0down=0.0_k_pr
      q0=0.0_k_pr
! spin down
      
      m=atomic%basis%norbitals*(atomic%basis%norbitals-1)/2
      from=m+atomic%atoms%orbs(at,1)
      to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2      
      q0down=sum(sol%density(from:to))            
      !! q0down=sum(sol%density(from:to))
!spin up
      from=m+atomic%atoms%orbs(at,1)+atomic%basis%norbitals/2
      to=-1+from+atomic%species%norbs(atomic%atoms%sp(at))/2
      q0up=sum(sol%density(from:to))      
! total charge
      q0=q0down+q0up      
   end subroutine ScfChargeNumbers



end module m_SCF2