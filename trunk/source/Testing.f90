module m_Testing
  use m_Constants
  use m_Types
  use m_Useful
  implicit none
  private

contains
!> \brief tests force on x
!> \details  moves atom 1 on the x axis from start in steps steps by dx each step
!>   in order to test that the energy and force are consistent 
!>    generates two outputs:
!>    fort.777 data file with has 3 columns (x,force,-total energy)
!>   can be visualised with xmgrace (xmgrace -nxy fort.777) 
!>   writes in the spciefied file the maximum error done in (%)
!> author Alin M. Elena (Belfast)
!> \date 22nd of June, 2006 generalized
  subroutine force_testx

    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'force_testx'
    !--subroutine parameters -------------------------!
    !--internal variables ----------------------------!   
    real(dp) :: eenergy,renergy,scfenergy,minusts

    !-------------------------------------------------!
    integer :: i
    real(dp) :: dx,tot
    real(dp), allocatable :: measured(:),analytic(:)

    allocate(measured(1:general%f_steps),analytic(1:general%f_steps))

    dx=general%f_dx
    do i=1,general%f_steps
      atomic%x(1)=general%f_start+(i-1)*dx    

      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()

      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      tot=eenergy+renergy+minusts+scfenergy
      write(777,*)atomic%x(1),atomic%fx(1),-tot
      measured(i)=-tot
      analytic(i)=atomic%fx(1)

    end do
    write(control_var%output_file,'(a,f12.8,a)')&
      'Maximum error',norm(measured,analytic,dx),'%'
    deallocate(measured,analytic)   
  end subroutine force_testx


!****s*   testing/force_testy()
! NAME
! force_testy
! SYNOPSIS
! call force_testy()
! INPUTS
! 
!
! DESCRIPTION
!    moves atom 1 on the y axis from general%f_start in general%f_steps steps by general%f_dx each step
!   in order to test that the energy and force are consistent 
!    generates two outputs:
!    fort.777 data file with has 3 columns (y,force,-total energy)
!   can be visualised with xmgrace (xmgrace -nxy fort.776) 
!   writes to control_var%output_file the maximum error done in (%)   
! USES
!   subroutine full_scf
!   functions  electronic_energy, repulsive_energy, scf_energy, norm
! AUTHOR
! Alin M. Elena (Belfast)              
! CREATION DATE
! 25th of April 2006
! HISTORY
! 22nd of June, 2006 generalized  
!******  


  subroutine force_testy
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'force_testy'
    !--subroutine parameters -------------------------!
    !--internal variables ----------------------------!   
    real(dp) :: eenergy,renergy,scfenergy,minusts

    !-------------------------------------------------!
    integer :: i
    real(dp) :: dx,tot

    real(dp), allocatable :: measured(:),analytic(:)

    allocate(measured(1:general%f_steps),analytic(1:general%f_steps))

    dx=general%f_dx
    do i=1,general%f_steps
      atomic%y(1)=general%f_start+(i-1)*dx    

      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()

      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      tot=eenergy+renergy+minusts+scfenergy
      write(776,*)atomic%y(1),atomic%fy(1),-tot
      measured(i)=-tot
      analytic(i)=atomic%fy(1)

    end do
    write(control_var%output_file,'(a,f12.8,a)')&
      'Maximum error',norm(measured,analytic,dx),'%'
    deallocate(measured,analytic)   
  end subroutine force_testy  

  !****s*   testing/force_testz()
! NAME
! force_testz
! SYNOPSIS
! call force_testz()
! INPUTS
! 
!
! DESCRIPTION
!    moves atom 1 on the z axis from general%f_start in general%f_steps steps by general%f_dx each step
!   in order to test that the energy and force are consistent 
!    generates two outputs:
!    fort.777 data file with has 3 columns (x,force,-total energy)
!   can be visualised with xmgrace (xmgrace -nxy fort.775) 
!   writes to control_var%output_file the maximum error done in (%)   
! USES
!   subroutine full_scf
!   functions  electronic_energy, repulsive_energy, scf_energy, norm
! AUTHOR
! Alin M. Elena (Belfast)              
! CREATION DATE
! 25th of April 2006
! HISTORY
! 22nd of June, 2006 generalized  
!******

  subroutine force_testz
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'force_testz'
    !--subroutine parameters -------------------------!
    !--internal variables ----------------------------!   
    real(dp) :: eenergy,renergy,scfenergy,minusts

    !-------------------------------------------------!
    integer :: i
    real(dp) :: dx,tot

    real(dp), allocatable :: measured(:),analytic(:)
    allocate(measured(1:general%f_steps),analytic(1:general%f_steps))

    dx=general%f_dx
    do i=1,general%f_steps
      atomic%z(1)=general%f_start+(i-1)*dx    

      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()

      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      tot=eenergy+renergy+minusts+scfenergy
      write(775,*)atomic%z(1),atomic%fz(1),-tot
      measured(i)=-tot
      analytic(i)=atomic%fz(1)

    end do
    write(control_var%output_file,'(a,f12.8,a)')&
      'Maximum error',norm(measured,analytic,dx),'%'     
    deallocate(measured,analytic)
  end subroutine force_testz  


!****s*   testing/force_test()
! NAME
! force_test
! SYNOPSIS
! call force_test()
! INPUTS
! 
!
! DESCRIPTION
!  computes the analytic and the numeric force on atoms
!  for the numeric force shifts the x,y,z coordinates by a small 
!  quantity dx. In the end display a useful summary.   
! USES
!   subroutine full_scf
!   functions  electronic_energy, repulsive_energy, scf_energy
! AUTHOR
! Alin M. Elena (Belfast)
! CREATION DATE
! 2005
! HISTORY
! 1st of April, 2006 Modified by Alin M. Elena (Belfast) in order to display the total force on components for analytic and numeric case
!******


  subroutine force_test
    !--subroutine name--------------------------------!   
    character(len=*), parameter :: myname = 'force_test'
    !--subroutine parameters -------------------------!
    !--internal variables ----------------------------!  
    integer       :: i
    real(dp)      :: fa_x, fa_y, fa_z, fn_x, fn_y, fn_z
    real(dp)      :: eenergy, renergy, scfenergy, minusts
    real(dp)      :: ep, em, dx,totnx,totny,totnz
    !-------------------------------------------------!



    totnx=0.0_dp
    totny=0.0_dp
    totnz=0.0_dp

    dx = general%f_dx
    do i=1,atomic%natoms
       ! Analytic forces
      call full_scf
      fa_x = atomic%fx(i)
      fa_y = atomic%fy(i)
      fa_z = atomic%fz(i)

       ! Numeric x force
      atomic%x(i) = atomic%x(i) + dx
      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()
      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      ep = renergy+eenergy+minusts+scfenergy
      atomic%x(i) = atomic%x(i) - dx

      atomic%x(i) = atomic%x(i) - dx
      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()
      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      em = renergy+eenergy+minusts+scfenergy
      atomic%x(i) = atomic%x(i) + dx

      fn_x = -(ep - em) / (dx + dx)
      totnx=totnx+fn_x
       ! Numeric y force
      atomic%y(i) = atomic%y(i) + dx
      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()
      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      ep = renergy+eenergy+minusts+scfenergy
      atomic%y(i) = atomic%y(i) - dx

      atomic%y(i) = atomic%y(i) - dx
      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()
      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      em = renergy+eenergy+minusts+scfenergy
      atomic%y(i) = atomic%y(i) + dx

      fn_y = -(ep - em) / (dx + dx)
      totny=totny+fn_y
       ! Numeric z force
      atomic%z(i) = atomic%z(i) + dx
      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()
      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      ep = renergy+eenergy+minusts+scfenergy
      atomic%z(i) = atomic%z(i) - dx

      atomic%z(i) = atomic%z(i) - dx
      call full_scf

      eenergy = electronic_energy()
      renergy = repulsive_energy()
      if (general%scf) then 
        scfenergy = scf_energy()
      else
        scfenergy = 0.0_dp
      endif
      minusts = - general%electronic_temperature * electronic_entropy
      em = renergy+eenergy+minusts+scfenergy
      atomic%z(i) = atomic%z(i) + dx

      fn_z = -(ep - em) / (dx + dx)
      totnz=totnz+fn_z
      write(*,'(a,i0)') "Atom: ",i   
      write (*,'(a,F16.8,2X,F16.8,2X,F16.8)') 'Analytic: ',fa_x, fa_y, fa_z
      write (*,'(a,F16.8,2X,F16.8,2X,F16.8)') ' Numeric: ',fn_x, fn_y, fn_z
      write (*,'(a,F16.8,2X,F16.8,2X,F16.8)') 'An.-Num.: ',fa_x-fn_x, fa_y-fn_y, fa_z-fn_z
      write (*,*)

    end do
    fa_x=sum(atomic%fx(:))
    fa_y=sum(atomic%fy(:))
    fa_z=sum(atomic%fz(:))
    write(*,*)"Total forces:           x            y             z"
    write(*,'(a,F16.8,2X,F16.8,2X,F16.8)')"Analytic: ",fa_x,fa_y,fa_z
    write(*,'(a,F16.8,2X,F16.8,2X,F16.8)')" Numeric: ",totnx,totny,totnz
    write(*,'(a,F16.8,2X,F16.8,2X,F16.8)')"An.-Num.: ",fa_x-totnx,fa_y-totny,fa_z-totnz



  end subroutine force_test

end module m_Testing