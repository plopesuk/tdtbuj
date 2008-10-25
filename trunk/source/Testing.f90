module m_Testing
  use m_Constants
  use m_Types
  use m_Useful
  use m_DriverRoutines
  use m_TightBinding
  use m_Gutenberg
  implicit none
  private
!
!
  public :: ForceTest
  public :: forceTestx
  public :: forceTesty
  public :: forceTestz
  public :: testTails
  public :: centerMolecule
!
contains
!> \brief tests force on x
!> \details  moves atom on the x axis from start in steps steps by dx each step
!>   in order to test that the energy and force are consistent
!>    generates two outputs:
!>    fort.777 data file with has 3 columns (x,force,-total energy)
!>   can be visualised with xmgrace (xmgrace -nxy fort.777)
!>   writes in the spciefied file the maximum error done in (%)
!> author Alin M. Elena (Belfast)
!> \date 22nd of June, 2006 generalized
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine forceTestx (io, gen, atomic, tb, sol)
    character (len=*), parameter :: myname = 'forceTestx'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    integer :: i, atom
    real (k_pr) :: dx, tot
    real (k_pr), allocatable :: measured (:), analytic (:)
!
    if ((gen%fatom <= 0) .or. (gen%fatom > atomic%atoms%natoms)) then
      call error (myname, "can not calculate force on non existing atom", .true., io)
    end if
!
    allocate (measured(1:gen%fsteps), analytic(1:gen%fsteps))
    dx = gen%fdx
    atom = gen%fatom
    do i = 1, gen%fsteps
      atomic%atoms%x (atom) = gen%fstart + (i-1) * dx
      call SinglePoint (io, gen, atomic, tb, sol)
      tot = sol%totalEnergy
      write (777,*) atomic%atoms%x(atom), atomic%atoms%fx(atom), - tot
      measured (i) = - tot
      analytic (i) = atomic%atoms%fx(atom)
    end do
    write (io%uout, '(a,f12.8,a)') 'Maximum error', norm (measured, analytic, dx), '%'
    deallocate (measured, analytic)
  end subroutine forceTestx
!
!> \brief tests force on y
!> \details  moves atom on the y axis from start in steps steps by dx each step
!>   in order to test that the energy and force are consistent
!>    generates two outputs:
!>    fort.777 data file with has 3 columns (x,force,-total energy)
!>   can be visualised with xmgrace (xmgrace -nxy fort.777)
!>   writes in the spciefied file the maximum error done in (%)
!> author Alin M. Elena (Belfast)
!> \date 22nd of June, 2006 generalized
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine forceTesty (io, gen, atomic, tb, sol)
    character (len=*), parameter :: myname = 'forceTesty'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    integer :: i, atom
    real (k_pr) :: dx, tot
    real (k_pr), allocatable :: measured (:), analytic (:)
!
    if ((gen%fatom <= 0) .or. (gen%fatom > atomic%atoms%natoms)) then
      call error (myname, "can not calculate force on non existing atom", .true., io)
    end if
!
    allocate (measured(1:gen%fsteps), analytic(1:gen%fsteps))
    dx = gen%fdx
    atom = gen%fatom
    do i = 1, gen%fsteps
      atomic%atoms%y (atom) = gen%fstart + (i-1) * dx
      call SinglePoint (io, gen, atomic, tb, sol)
      tot = sol%totalEnergy
      write (777,*) atomic%atoms%y(atom), atomic%atoms%fy(atom), - tot
      measured (i) = - tot
      analytic (i) = atomic%atoms%fy(atom)
    end do
    write (io%uout, '(a,f12.8,a)') 'Maximum error', norm (measured, analytic, dx), '%'
    deallocate (measured, analytic)
  end subroutine forceTesty
!
!> \brief tests force on z
!> \details  moves atom on the z axis from start in steps steps by dx each step
!>   in order to test that the energy and force are consistent
!>    generates two outputs:
!>    fort.777 data file with has 3 columns (x,force,-total energy)
!>   can be visualised with xmgrace (xmgrace -nxy fort.777)
!>   writes in the spciefied file the maximum error done in (%)
!> author Alin M. Elena (Belfast)
!> \date 22nd of June, 2006 generalized
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine forceTestz (io, gen, atomic, tb, sol)
    character (len=*), parameter :: myname = 'forceTestx'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    integer :: i, atom
    real (k_pr) :: dx, tot
    real (k_pr), allocatable :: measured (:), analytic (:)
!
    if ((gen%fatom <= 0) .or. (gen%fatom > atomic%atoms%natoms)) then
      call error (myname, "can not calculate force on non existing atom", .true., io)
    end if
!
    allocate (measured(1:gen%fsteps), analytic(1:gen%fsteps))
    dx = gen%fdx
    atom = gen%fatom
    do i = 1, gen%fsteps
      atomic%atoms%z (atom) = gen%fstart + (i-1) * dx
      call SinglePoint (io, gen, atomic, tb, sol)
      tot = sol%totalEnergy
      write (777,*) atomic%atoms%z(atom), atomic%atoms%fz(atom), - tot
      measured (i) = - tot
      analytic (i) = atomic%atoms%fz(atom)
    end do
    write (io%uout, '(a,f12.8,a)') 'Maximum error', norm (measured, analytic, dx), '%'
    deallocate (measured, analytic)
  end subroutine forceTestz
!
!> \brief  computes the analytic and the numeric force on atoms
!> \details  for the numeric force shifts the x,y,z coordinates by a small
!>  quantity dx. In the end display a useful summary.
!> \author Alin M. Elena (Belfast)
!> \date 1st of April, 2006
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine ForceTest (io, gen, atomic, tb, sol)
    character (len=*), parameter :: myname = 'ForceTest'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
    integer :: i
    real (k_pr) :: fa_x, fa_y, fa_z, fn_x, fn_y, fn_z
    real (k_pr) :: eenergy, renergy, scfenergy, minusts
    real (k_pr) :: ep, em, dx, totnx, totny, totnz
!-------------------------------------------------!
!
    totnx = 0.0_k_pr
    totny = 0.0_k_pr
    totnz = 0.0_k_pr
    dx = gen%fdx
    do i = 1, atomic%atoms%natoms
! Analytic forces
      call SinglePoint (io, gen, atomic, tb, sol)
      fa_x = atomic%atoms%fx(i)
      fa_y = atomic%atoms%fy(i)
      fa_z = atomic%atoms%fz(i)
! Numeric x force
      atomic%atoms%x (i) = atomic%atoms%x(i) + dx
      call SinglePoint (io, gen, atomic, tb, sol)
      ep = sol%totalEnergy
      atomic%atoms%x (i) = atomic%atoms%x(i) - 2.0_k_pr * dx
      call SinglePoint (io, gen, atomic, tb, sol)
      em = sol%totalEnergy
      atomic%atoms%x (i) = atomic%atoms%x(i) + dx
      fn_x = - (ep-em) / (dx+dx)
      totnx = totnx + fn_x
! Numeric y force
      atomic%atoms%y (i) = atomic%atoms%y(i) + dx
      call SinglePoint (io, gen, atomic, tb, sol)
      ep = sol%totalEnergy
      atomic%atoms%y (i) = atomic%atoms%y(i) - 2.0_k_pr * dx
      call SinglePoint (io, gen, atomic, tb, sol)
      em = sol%totalEnergy
      atomic%atoms%y (i) = atomic%atoms%y(i) + dx
      fn_y = - (ep-em) / (dx+dx)
      totny = totny + fn_y
! Numeric z force
      atomic%atoms%z (i) = atomic%atoms%z(i) + dx
      call SinglePoint (io, gen, atomic, tb, sol)
      ep = sol%totalEnergy
      atomic%atoms%z (i) = atomic%atoms%z(i) - 2.0_k_pr * dx
      call SinglePoint (io, gen, atomic, tb, sol)
      em = sol%totalEnergy
      atomic%atoms%z (i) = atomic%atoms%z(i) + dx
      fn_z = - (ep-em) / (dx+dx)
      totnz = totnz + fn_z
!
      write (*, '(a,i0)') "Atom: ", i
      write (*, '(a,F16.8,2X,F16.8,2X,F16.8)') 'Analytic: ', fa_x, fa_y, fa_z
      write (*, '(a,F16.8,2X,F16.8,2X,F16.8)') ' Numeric: ', fn_x, fn_y, fn_z
      write (*, '(a,F16.8,2X,F16.8,2X,F16.8)') 'An.-Num.: ', fa_x - fn_x, fa_y - fn_y, fa_z - fn_z
      write (*,*)
    end do
    fa_x = sum (atomic%atoms%fx(:))
    fa_y = sum (atomic%atoms%fy(:))
    fa_z = sum (atomic%atoms%fz(:))
    write (*,*) "Total forces:           x            y             z"
    write (*, '(a,F16.8,2X,F16.8,2X,F16.8)') "Analytic: ", fa_x, fa_y, fa_z
    write (*, '(a,F16.8,2X,F16.8,2X,F16.8)') " Numeric: ", totnx, totny, totnz
    write (*, '(a,F16.8,2X,F16.8,2X,F16.8)') "An.-Num.: ", fa_x - totnx, fa_y - totny, fa_z - totnz
  end subroutine ForceTest
!
!
!> \brief  prints the hopping integrals and repulsive potential radial components
!> \author Alin M. Elena (Belfast)
!> \date 3rd of June, 2008
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine testTails (io, gen, atomic, tb, sol)
    character (len=*), parameter :: myname = 'testTails'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    real (k_pr) :: r
    integer :: z, y, three, i, j, k, k1, k2
    character (len=k_mw) :: saux
!
!
    three = 3
    do i = 1, atomic%species%nspecies
      do j = 1, atomic%species%nspecies
        do k = 0, tb%hopping(i, j)%l1
          do k1 = 0, tb%hopping(i, j)%l2
            do k2 = 0, Min (k, k1)
              z = GetUnit ()
              write (saux, '(a,a)') "h_", trim (ccnlm(i, j, k, k1, k2))
              open (unit=z, file=trim(saux), status="replace", action="write")
              write (z,*) "# r ", trim (ccnlm(i, j, k, k1, k2))
              do y = 50, 1000
                r = y * 0.01_k_pr
                write (z,*) r, rad (r, atomic%species%id(i), atomic%species%id(j), gen, tb, k, k1, k2), radp (r, &
               & atomic%species%id(i), atomic%species%id(j), gen, tb, k, k1, k2)
              end do
              close (z)
            end do
          end do
        end do
      end do
    end do
    do i = 1, atomic%species%nspecies
      do j = 1, atomic%species%nspecies
        z = GetUnit ()
        write (saux, '(a,2i0)') "phi", i, j
        open (unit=z, file=trim(saux), status="replace", action="write")
        write (z, '(a,i0,x,i0)') "# r ", i, j
        do y = 50, 1000
          r = y * 0.01_k_pr
          write (z,*) r, rep (r, atomic%species%id(i), atomic%species%id(j), tb, gen), repp (r, atomic%species%id(i), &
         & atomic%species%id(j), tb, gen)
        end do
        close (z)
      end do
    end do
  end subroutine testTails
!
  subroutine centerMolecule (io, gen, atomic, tb, sol)
    character (len=*), parameter :: myname = 'centerMolecule'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    integer :: i
    real (k_pr) :: cx, cy, cz, l, ca, sa, nx, ny, nz
    real (k_pr), dimension (3) :: a, b, c, n
    cx = 0.0_k_pr
    cy = 0.0_k_pr
    cz = 0.0_k_pr
    do i = 1, atomic%atoms%natoms
      cx = cx + atomic%atoms%x(i)
      cy = cy + atomic%atoms%y(i)
      cz = cz + atomic%atoms%z(i)
    end do
    cx = cx / real (atomic%atoms%natoms, k_pr)
    cy = cy / real (atomic%atoms%natoms, k_pr)
    cz = cz / real (atomic%atoms%natoms, k_pr)
    do i = 1, atomic%atoms%natoms
      atomic%atoms%x (i) = atomic%atoms%x(i) - cx
      atomic%atoms%y (i) = atomic%atoms%y(i) - cy
      atomic%atoms%z (i) = atomic%atoms%z(i) - cz
    end do
!
    a (1) = atomic%atoms%x(2) - atomic%atoms%x(1)
    a (2) = atomic%atoms%y(2) - atomic%atoms%y(1)
    a (3) = atomic%atoms%z(2) - atomic%atoms%z(1)
    b (1) = atomic%atoms%x(3) - atomic%atoms%x(1)
    b (2) = atomic%atoms%y(3) - atomic%atoms%y(1)
    b (3) = atomic%atoms%z(3) - atomic%atoms%z(1)
    n (1) = a (2) * b (3) - a (3) * b (2)
    n (2) = - a (1) * b (3) + b (1) * a (3)
    n (3) = a (1) * b (2) - b (1) * a (2)
    l = n (1) ** 2 + n (2) ** 2 + n (3) ** 2
    n (1) = n (1) / l - atomic%atoms%x(1)
    n (2) = n (2) / l - atomic%atoms%y(1)
    n (3) = n (3) / l - atomic%atoms%z(1)
    ca = n (1) / Sqrt (n(1)**2+n(2)**2)
    sa = n (2) / Sqrt (n(1)**2+n(2)**2)
    do i = 1, atomic%atoms%natoms
      nx = atomic%atoms%x(i)
      atomic%atoms%x (i) = ca * nx - sa * atomic%atoms%y(i)
      atomic%atoms%y (i) = sa * nx + ca * atomic%atoms%y(i)
      atomic%atoms%z (i) = nz
    end do
    call PrintXYZ (998, atomic, .false., "centred geometry")
!
  end subroutine centerMolecule
!
end module m_Testing
