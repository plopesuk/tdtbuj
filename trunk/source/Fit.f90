!> \brief controls the fitting process
!> \author Alin M Elena
!> \date 14/11/07, 10:11:09
!> \remarks This module is highly non standard and using it means customizing part of the code.\n
!> There are few routines that need to me modified according to the problem that needs solved.
!> - UpdateParams
!>  sets the parameters to the values contained by the yfit array. yfit array is the only parameter that
!>  should be present other may be added if necessary
!> - SetInitialParams
!>  initializes the parameters to be fit. Everything is set in the array x which is the only parameter
!>  that should be present at all time. Other may appear if necessary
!> - UpdateCost computes the cost for our fit for a set of fit parameters given by the array yfit. The rest of parameters should not
!> - PrintFit prints the fit second part of the routine may be used to print in the desired format the fit (easy to use in the progr
!
module m_Fit
  use m_Constants
  use m_Types
  use m_Useful
  use m_SA
  use m_Simplex
  use m_SimplexSA
  use m_Useful
  use m_DriverRoutines
  use m_Gutenberg
  use m_TightBinding
  implicit none
  private
  public :: fitting
!
  integer, external :: djacobi_init, djacobi_solve, djacobi_delete, djacobi
  integer, external :: dtrnlspbc_get, dtrnlspbc_init, dtrnlspbc_solve, dtrnlspbc_delete
!
  type, public :: fitDataType
    real (k_pr), pointer :: x (:)
    real (k_pr), pointer :: exper (:)
    real (k_pr), pointer :: fit (:)
    integer :: n
  end type fitDataType
!
  type (fitDataType), public :: fitData
  real (k_pr), public, pointer :: y (:), p (:, :), best (:), bounds (:, :)
!
! private variables
  character (len=100), parameter :: cExper = "experiment.dat"
  character (len=100), parameter :: cOutFit = "fitout.dat"
  character (len=100), parameter :: cBounds = "bounds.dat"
  character (len=100), parameter :: cRestart = "restartfit.dat"
  character (len=100), parameter :: cFdfFile = "new_param.fdf"
  character (len=100), parameter :: cAtomicData = "new_AtomicData.fdf"
!
contains
  subroutine SetInitialParams (x, atomic, tb)
    character (len=*), parameter :: sMyName = "SetInitialParams"
    real (k_pr), intent (inout) :: x (:)
    type (atomicxType), intent (in) :: atomic
    type (modelType), intent (in) :: tb
!
    x (1) = tb%hopping(1, 1)%eps(0)
    x (2) = 0.1_k_pr
    x (3) = tb%hopping(2, 2)%eps(0)
    x (4) = tb%hopping(3, 3)%eps(0)
    x (5) = 0.1_k_pr
    x (6) = tb%hopping(4, 4)%eps(0)
    x (7) = 0.1_k_pr
    x (8) = atomic%species%ulocal(atomic%atoms%sp(1), 1)
    x (9) = atomic%species%jlocal(atomic%atoms%sp(1), 1)
    x (10) = atomic%species%ulocal(atomic%atoms%sp(2), 1)
    x (11) = atomic%species%jlocal(atomic%atoms%sp(2), 1)
    x (12) = atomic%species%ulocal(atomic%atoms%sp(3), 1)
    x (13) = atomic%species%jlocal(atomic%atoms%sp(3), 1)
    x (14) = atomic%species%ulocal(atomic%atoms%sp(4), 1)
    x (15) = atomic%species%jlocal(atomic%atoms%sp(4), 1)
    x (16) = tb%hopping(1, 1)%a(0, 0, 0)
    x (17) = tb%hopping(1, 1)%a(0, 1, 0)
    x (18) = tb%hopping(1, 1)%a(1, 1, 0)
    x (19) = tb%hopping(1, 1)%a(1, 1, 1)
    x (20) = tb%hopping(1, 1)%r0
    x (21) = tb%hopping(1, 1)%rc
    x (22) = tb%hopping(1, 1)%n
    x (23) = tb%hopping(1, 1)%nc
    x (24) = tb%hopping(4, 4)%a(0, 0, 0)
    x (25) = tb%hopping(4, 4)%a(0, 1, 0)
    x (26) = tb%hopping(4, 4)%a(1, 1, 0)
    x (27) = tb%hopping(4, 4)%a(1, 1, 1)
    x (28) = tb%hopping(4, 4)%r0
    x (29) = tb%hopping(4, 4)%rc
    x (30) = tb%hopping(4, 4)%n
    x (31) = tb%hopping(4, 4)%nc
    x (32) = tb%hopping(2, 2)%a(0, 0, 0)
    x (33) = tb%hopping(2, 2)%r0
    x (34) = tb%hopping(2, 2)%rc
    x (35) = tb%hopping(2, 2)%n
    x (36) = tb%hopping(2, 2)%nc
    x (37) = tb%hopping(1, 2)%a(0, 0, 0)
    x (38) = tb%hopping(1, 2)%a(1, 0, 0)
    x (39) = tb%hopping(1, 2)%r0
    x (40) = tb%hopping(1, 2)%rc
    x (41) = tb%hopping(1, 2)%n
    x (42) = tb%hopping(1, 2)%nc
    x (43) = tb%hopping(3, 2)%a(0, 0, 0)
    x (44) = tb%hopping(3, 2)%a(1, 0, 0)
    x (45) = tb%hopping(3, 2)%r0
    x (46) = tb%hopping(3, 2)%rc
    x (47) = tb%hopping(3, 2)%n
    x (48) = tb%hopping(3, 2)%nc
    x (49) = tb%hopping(1, 3)%a(0, 0, 0)
    x (50) = tb%hopping(1, 3)%a(0, 1, 0)
    x (51) = tb%hopping(1, 3)%a(1, 1, 0)
    x (52) = tb%hopping(1, 3)%a(1, 1, 1)
    x (53) = tb%hopping(1, 3)%r0
    x (54) = tb%hopping(1, 3)%rc
    x (55) = tb%hopping(1, 3)%n
    x (56) = tb%hopping(1, 3)%nc
    x (57) = tb%hopping(3, 4)%a(0, 0, 0)
    x (58) = tb%hopping(3, 4)%a(0, 1, 0)
    x (59) = tb%hopping(3, 4)%a(1, 1, 0)
    x (60) = tb%hopping(3, 4)%a(1, 1, 1)
    x (61) = tb%hopping(3, 4)%r0
    x (62) = tb%hopping(3, 4)%rc
    x (63) = tb%hopping(3, 4)%n
    x (64) = tb%hopping(3, 4)%nc
  end subroutine SetInitialParams
!
  subroutine UpdateParams (yfit, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = 'UpdateParams'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
    real (k_pr), intent (in) :: yfit (:)
!
    tb%hopping(1, 1)%eps(0) = yfit (1)
    tb%hopping(1, 1)%eps(1:3) = yfit (1) + yfit (2)
    tb%hopping(1, 1)%eps(4) = yfit (1)
    tb%hopping(1, 1)%eps(5:7) = yfit (1) + yfit (2)
!
    tb%hopping(2, 2)%eps(0:1) = yfit (3)
    tb%hopping(3, 3)%eps(0) = yfit (4)
    tb%hopping(3, 3)%eps(1:3) = yfit (4) + yfit (5)
    tb%hopping(3, 3)%eps(4) = yfit (4)
    tb%hopping(3, 3)%eps(5:7) = yfit (4) + yfit (5)
!
    tb%hopping(4, 4)%eps(0) = yfit (6)
    tb%hopping(4, 4)%eps(1:3) = yfit (6) + yfit (7)
    tb%hopping(4, 4)%eps(4) = yfit (6)
    tb%hopping(4, 4)%eps(5:7) = yfit (6) + yfit (7)
    atomic%species%ulocal (atomic%atoms%sp(1), 1) = yfit (8)
    atomic%species%jlocal (atomic%atoms%sp(1), 1) = yfit (9)
    atomic%species%ulocal (atomic%atoms%sp(2), 1) = yfit (10)
    atomic%species%jlocal (atomic%atoms%sp(2), 1) = yfit (11)
    atomic%species%ulocal (atomic%atoms%sp(3), 1) = yfit (12)
    atomic%species%jlocal (atomic%atoms%sp(3), 1) = yfit (13)
    atomic%species%ulocal (atomic%atoms%sp(4), 1) = yfit (14)
    atomic%species%jlocal (atomic%atoms%sp(4), 1) = yfit (15)
    tb%hopping(1, 1)%a(0, 0, 0) = yfit (16)
    tb%hopping(1, 1)%a(0, 1, 0) = yfit (17)
    tb%hopping(1, 1)%a(1, 0, 0) = - tb%hopping(1, 1)%a(0, 1, 0)
    tb%hopping(1, 1)%a(1, 1, 0) = yfit (18)
    tb%hopping(1, 1)%a(1, 1, 1) = yfit (19)
    tb%hopping(1, 1)%r0 = yfit (20)
    tb%hopping(1, 1)%rc = yfit (21)
    tb%hopping(1, 1)%n = yfit (22)
    tb%hopping(1, 1)%nc = yfit (23)
    tb%hopping(4, 4)%a(0, 0, 0) = yfit (24)
    tb%hopping(4, 4)%a(0, 1, 0) = yfit (25)
    tb%hopping(4, 4)%a(1, 0, 0) = - tb%hopping(4, 4)%a(0, 1, 0)
    tb%hopping(4, 4)%a(1, 1, 0) = yfit (26)
    tb%hopping(4, 4)%a(1, 1, 1) = yfit (27)
    tb%hopping(4, 4)%r0 = yfit (28)
    tb%hopping(4, 4)%rc = yfit (29)
    tb%hopping(4, 4)%n = yfit (30)
    tb%hopping(4, 4)%nc = yfit (31)
    tb%hopping(2, 2)%a(0, 0, 0) = yfit (32)
    tb%hopping(2, 2)%r0 = yfit (33)
    tb%hopping(2, 2)%rc = yfit (34)
    tb%hopping(2, 2)%n = yfit (35)
    tb%hopping(2, 2)%nc = yfit (36)
    tb%hopping(1, 2)%a(0, 0, 0) = yfit (37)
    tb%hopping(2, 1)%a(0, 0, 0) = yfit (37)
    tb%hopping(1, 2)%a(1, 0, 0) = yfit (38)
    tb%hopping(2, 1)%a(0, 1, 0) = - yfit (38)
    tb%hopping(1, 2)%r0 = yfit (39)
    tb%hopping(1, 2)%rc = yfit (40)
    tb%hopping(1, 2)%n = yfit (41)
    tb%hopping(1, 2)%nc = yfit (42)
    tb%hopping(2, 1)%r0 = yfit (39)
    tb%hopping(2, 1)%rc = yfit (40)
    tb%hopping(2, 1)%n = yfit (41)
    tb%hopping(2, 1)%nc = yfit (42)
!
    tb%hopping(2, 3)%a(0, 0, 0) = yfit (43)
    tb%hopping(3, 2)%a(0, 0, 0) = yfit (43)
    tb%hopping(3, 2)%a(1, 0, 0) = yfit (44)
    tb%hopping(2, 3)%a(0, 1, 0) = - yfit (44)
    tb%hopping(2, 3)%r0 = yfit (45)
    tb%hopping(2, 3)%rc = yfit (46)
    tb%hopping(2, 3)%n = yfit (47)
    tb%hopping(2, 3)%nc = yfit (48)
    tb%hopping(3, 2)%r0 = yfit (45)
    tb%hopping(3, 2)%rc = yfit (46)
    tb%hopping(3, 2)%n = yfit (47)
    tb%hopping(3, 2)%nc = yfit (48)
!
    tb%hopping(1, 3)%a(0, 0, 0) = yfit (49)
    tb%hopping(1, 3)%a(0, 1, 0) = yfit (50)
    tb%hopping(1, 3)%a(1, 0, 0) = - tb%hopping(1, 3)%a(0, 1, 0)
    tb%hopping(1, 3)%a(1, 1, 0) = yfit (51)
    tb%hopping(1, 3)%a(1, 1, 1) = yfit (52)
    tb%hopping(1, 3)%r0 = yfit (53)
    tb%hopping(1, 3)%rc = yfit (54)
    tb%hopping(1, 3)%n = yfit (55)
    tb%hopping(1, 3)%nc = yfit (56)
    tb%hopping(3, 1)%a(0, 0, 0) = yfit (49)
    tb%hopping(3, 1)%a(0, 1, 0) = yfit (50)
    tb%hopping(3, 1)%a(1, 0, 0) = - tb%hopping(3, 1)%a(0, 1, 0)
    tb%hopping(3, 1)%a(1, 1, 0) = yfit (51)
    tb%hopping(3, 1)%a(1, 1, 1) = yfit (52)
    tb%hopping(3, 1)%r0 = yfit (53)
    tb%hopping(3, 1)%rc = yfit (54)
    tb%hopping(3, 1)%n = yfit (55)
    tb%hopping(3, 1)%nc = yfit (56)
!
    tb%hopping(3, 4)%a(0, 0, 0) = yfit (57)
    tb%hopping(3, 4)%a(0, 1, 0) = yfit (58)
    tb%hopping(3, 4)%a(1, 0, 0) = - tb%hopping(3, 4)%a(0, 1, 0)
    tb%hopping(3, 4)%a(1, 1, 0) = yfit (59)
    tb%hopping(3, 4)%a(1, 1, 1) = yfit (60)
    tb%hopping(3, 4)%r0 = yfit (61)
    tb%hopping(3, 4)%rc = yfit (62)
    tb%hopping(3, 4)%n = yfit (63)
    tb%hopping(3, 4)%nc = yfit (64)
    tb%hopping(4, 3)%a(0, 0, 0) = yfit (57)
    tb%hopping(4, 3)%a(0, 1, 0) = yfit (58)
    tb%hopping(4, 3)%a(1, 0, 0) = - tb%hopping(4, 3)%a(0, 1, 0)
    tb%hopping(4, 3)%a(1, 1, 0) = yfit (59)
    tb%hopping(4, 3)%a(1, 1, 1) = yfit (60)
    tb%hopping(4, 3)%r0 = yfit (61)
    tb%hopping(4, 3)%rc = yfit (62)
    tb%hopping(4, 3)%n = yfit (63)
    tb%hopping(4, 3)%nc = yfit (64)
    call setTails (io, gen, atomic, tb, sol)
  end subroutine UpdateParams
!
!> \brief computes the cost function for a set o parameters
!> \author Alin M Elena
!> \date 14/11/07, 09:59:55
!> \param yfit real array the list of parameters for which we compute the cost function
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  function UpdateCost (yfit, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = 'UpdateCost'
    real (k_pr), dimension (:), intent (in) :: yfit
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    real (k_pr) :: UpdateCost
    integer :: i, ilevels
    real (k_pr) :: aux
!
    aux = 0.0_k_pr
    call UpdateParams (yfit, gen, atomic, tb, sol, io)
    call SinglePoint (io, gen, atomic, tb, sol)
    if ( .not. gen%lIsSCFConverged) then
      UpdateCost = k_infinity
    else
      do i = 1, atomic%atoms%natoms
        aux = aux + (fitData%exper(i)-atomic%atoms%chrg(i)) ** 2
      end do
      UpdateCost = aux
    end if
  end function UpdateCost
!
!> \brief prints the fit results
!> \author Alin M Elena
!> \date 14/11/07, 12:18:55
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine PrintFit (gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "PrintFit"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
!
    integer :: i, j
    real (k_pr) :: aux
    aux = UpdateCost (best, gen, atomic, tb, sol, io)
    fitData%fit = atomic%atoms%chrg
    open (1, file=trim(cOutFit), status="unknown", action="write")
    do i = 1, fitData%n
      write (1, '(i6,4g16.6)') i, fitData%x(i), fitData%exper(i), fitData%fit(i), fitData%exper(i) - fitData%fit(i)
    end do
    write (1, '(a,f16.8)') "cost:", aux
    close (1)
    open (1, file=trim(cAtomicData), status="unknown", action="write")
    write (1, '(a)') "%block AtomicData"
    do j = 1, atomic%species%nspecies
      write (1, '(i0,1x,i0,4f16.8)') j, atomic%species%z(j), atomic%species%mass(j) / k_amuToInternal, atomic%species%ulocal(j, 1), &
     & atomic%species%jlocal(j, 1), atomic%species%uinter(j)
    end do
    write (1, '(a)') "%endblock AtomicData"
!
    open (1, file=trim(cFdfFile), status="unknown", action="write")
    i = io%uout
    io%uout = 1
    call PrintTbGSP (io, gen, atomic, tb)
    io%uout = i
    close (1)
  end subroutine PrintFit
!
!
!> \brief Reads the "experimental" data against which we do the fit
!> \author Alin M Elena
!> \date 14/11/07, 09:58:53
!> \param io type(ioType) contains all the info about I/O files
  subroutine InitFit (io)
    character (len=*), parameter :: myname = 'InitFit'
    type (ioType), intent (inout) :: io
    real (k_pr) :: r1, e1, ri, ei
    integer :: n, err, i
    open (1, file=trim(cExper), status="old", action="read", iostat=err)
    if (err /= 0) then
      call error ("I can not open file "//trim(cExper), myname, .true., io)
    end if
    read (1,*, iostat=err) r1, e1
    if (err /= 0) then
      call error ("I can not read from the file "//trim(cExper), myname, .true., io)
    end if
    err = 0
    n = 0
    do while (err /=- 1)
      read (1,*, iostat=err) ri, ei
      n = n + 1
    end do
    close (1)
    fitData%n = n
    allocate (fitData%x(1:n), fitData%exper(1:n), fitData%fit(1:n))
    open (1, file=trim(cExper), status="old", action="read")
    do i = 1, n
      read (1,*) fitData%x(i), fitData%exper(i)
    end do
    close (1)
  end subroutine InitFit
!> \brief reads the bounds of each parameter
!> \author Alin M Elena
!> \date 14/11/07, 10:01:40
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to run
!> \remarks column one keeps the lower limit column two the upper one
  subroutine ReadBounds (gen, io)
    character (len=*), parameter :: myname = "ReadBounds"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    integer :: i, errno
    character (len=k_ml) :: saux
!
    allocate (bounds(1:gen%fit%iNoParams, 1:2))
    open (1, file=trim(cBounds), status="old", action="read", iostat=errno)
!
    if (errno /= 0) then
      call error ("I can not find bounds file "//trim(cBounds), myname, .true., io)
    end if
    do i = 1, gen%fit%iNoParams
      read (1,*, iostat=errno) bounds (i, 1), bounds (i, 2)
      if (errno /= 0) then
        write (saux, '(a,a,i0)'), trim (cBounds), " error in line ", i
        call error (trim(saux), myname, .true., io)
      end if
    end do
    close (1)
    call PrintVectorA (bounds(:, 1), 'lower bound', .true., .false., io)
    call PrintVectorA (bounds(:, 2), 'upper bound', .true., .false., io)
  end subroutine ReadBounds
!
!> \brief Driver routine for the fitting process
!> \author Alin M Elena
!> \date 14/11/07, 09:59:55
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine fitting (io, gen, atomic, tb, sol)
    character (len=*), parameter :: myname = "Fitting"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    integer :: i, j, kk, un !,mm
    real (k_pr) :: copt
    real (k_pr), allocatable :: tol (:)
    logical :: quit
!
    call InitFit (io)
    allocate (tol(1:gen%fit%neps))
    select case (gen%fit%fitMethod)
    case (k_simplex)
! simplex fit                 
      call InitSimplex (gen, atomic, tb, sol, io)
!print on the screen initial p and y          
      write (io%uout,*) "---------------------------------------------------"
      write (io%uout,*) "initial values of parameters on first line last value the cost"
      do i = 1, gen%fit%iNoParams + 1
        write (io%uout, '(11f16.8,f20.8)') (p(i, j), j=1, gen%fit%iNoParams), y (i)
      end do
      write (io%uout,*) "---------------------------------------------------"
! end print        
      allocate (best(1:gen%fit%iNoParams), tol(1:gen%fit%neps))
      tol (1:gen%fit%neps) = 1.0e25_k_pr
      call amoeba (p, y, gen%fit%fitTol, UpdateCost, kk, bounds, gen%fit%iter, gen, atomic, tb, sol, io)
      copt = y (1)
      call PrintVectorA (p(1, :), 'current parameters:', .true., .false., io)
      do
        call amoeba (p, y, gen%fit%fitTol, UpdateCost, kk, bounds, gen%fit%iter, gen, atomic, tb, sol, io)
        do i = gen%fit%neps, 2, - 1
          tol (i) = tol (i-1)
        end do
        tol (1) = y (1)
        call PrintVectorA (p(1, :), 'current parameters:', .true., .false., io)
        quit = .false.
        if (Abs(copt-tol(1)) < gen%fit%fitTol) then
          quit = .true.
        end if
!
        if (copt > y(1)) then
          best (:) = p (1, :)
          copt = y (1)
        end if
        write (io%uout,*) "current  ", y (1), "and best cost functions", copt
        do i = 2, gen%fit%neps
          if (Abs(tol(1)-tol(i)) > gen%fit%fitTol) then
            quit = .false.
          end if
        end do
        call PrintVectorA (best, 'best parameters so far:', .true., .false., io)
        if (quit) then
          call PrintVectorA (best, 'optimal parameters:', .true., .false., io)
          exit
        end if
        p (1, :) = best (:)
        do i = 2, gen%fit%iNoParams + 1
          do j = 1, gen%fit%iNoParams
            if (i-1 == j) then
              p (i, j) = bounds (j, 1) + Abs (bounds(j, 2)-bounds(j, 1)) * ranmar (sol%seed)
            else
              p (i, j) = p (1, j)
            end if
          end do
        end do
        do i = 1, gen%fit%iNoParams + 1
          y (i) = UpdateCost (p(i, :), gen, atomic, tb, sol, io)
        end do
      end do
      call UpdateParams (best, gen, atomic, tb, sol, io)
      call PrintFit (gen, atomic, tb, sol, io)
      open (2, file=trim(cRestart), status='unknown', action="write")
      do i = 1, gen%fit%iNoParams
        write (2,*) best (i)
      end do
      close (2)
!     call print_gsp(trim(cFdfFile))  
! end simplex fit          
      deallocate (y, p, best)
    case (k_SimplexSA)
      allocate (best(1:gen%fit%iNoParams))
      call SimplexSA (gen, atomic, tb, sol, io)
      do i = 1, gen%fit%neps
        tol (i) = 1.0e25_k_pr
      end do
      copt = UpdateCost (best, gen, atomic, tb, sol, io)
      write (io%uout,*) "simplex step"
      do
        p (1, :) = best (:)
        do i = 2, gen%fit%iNoParams + 1
          do j = 1, gen%fit%iNoParams
            if (i-1 == j) then
              p (i, j) = bounds (j, 1) + Abs (bounds(j, 2)-bounds(j, 1)) * ranmar (sol%seed)
            else
              p (i, j) = p (1, j)
            end if
          end do
        end do
        do i = 1, gen%fit%iNoParams + 1
          y (i) = UpdateCost (p(i, :), gen, atomic, tb, sol, io)
        end do
        call amoeba (p, y, gen%fit%fitTol, UpdateCost, kk, bounds, gen%fit%iter, gen, atomic, tb, sol, io)
        do i = gen%fit%neps, 2, - 1
          tol (i) = tol (i-1)
        end do
        tol (1) = y (1)
        call PrintVectorA (p(1, :), 'current parameters:', .true., .false., io)
        quit = .false.
        if (Abs(copt-tol(1)) < gen%fit%fitTol) then
          quit = .true.
        end if
        if (copt > y(1)) then
          best (:) = p (1, :)
          copt = y (1)
        end if
        write (io%uout,*) "current  ", y (1), "and best cost functions", copt
        do i = 2, gen%fit%neps
          if (Abs(tol(1)-tol(i)) > gen%fit%fitTol) then
            quit = .false.
          end if
        end do
        call PrintVectorA (best, 'best parameters so far:', .true., .false., io)
        if (quit) then
          call PrintVectorA (best, 'optimal parameters:', .true., .false., io)
          exit
        end if
      end do
      call PrintFit (gen, atomic, tb, sol, io)
      open (2, file=trim(cRestart), status='unknown', action="write")
      do i = 1, gen%fit%iNoParams
        write (2,*) best (i)
      end do
      close (2)
      deallocate (y, p, best)
    case (k_SA)
      allocate (best(1:gen%fit%iNoParams))
      call InitSA (gen, atomic, tb, sol, io)
      call PrintFit (gen, atomic, tb, sol, io)
      un = GetUnit ()
      open (un, file=trim(cRestart), status='unknown', action="write")
      do i = 1, gen%fit%iNoParams
        write (un,*) best (i)
      end do
      close (un)
      deallocate (best)
    case (k_TrustRegion)
      allocate (best(1:gen%fit%iNoParams))
      call TrustRegion (gen, atomic, tb, sol, io)
      call PrintFit (gen, atomic, tb, sol, io)
      un = GetUnit ()
      open (un, file=trim(cRestart), status='unknown', action="write")
      do i = 1, gen%fit%iNoParams
        write (un,*) best (i)
      end do
      close (un)
      deallocate (best)
    end select
!
    deallocate (fitData%x, fitData%exper, fitData%fit)
    deallocate (bounds)
    deallocate (tol)
  end subroutine fitting
!
!> \brief initializes the simplex method
!> \author Alin M Elena
!> \date 14/11/07, 13:35:55
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine InitSimplex (gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = 'InitSimplex'
    integer :: i, j, errno
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    allocate (y(1:gen%fit%iNoParams+1), p(1:gen%fit%iNoParams+1, gen%fit%iNoParams))
    if (gen%fit%RestartFit) then
      open (2, file=trim(cRestart), status="old", action="read", iostat=errno)
      if (errno /= 0) then
        call error ("restart file not found!!!", myname, .true., io)
      end if
      do i = 1, gen%fit%iNoParams
        read (2,*) p (1, i)
      end do
      close (2)
      call UpdateParams (p(1, :), gen, atomic, tb, sol, io)
    else
      call SetInitialParams (p(1, :), atomic, tb)
    end if
    call ReadBounds (gen, io)
    do i = 1, gen%fit%iNoParams
      if ((p(1, i) < bounds(i, 1)) .or. (p(1, i) > bounds(i, 2))) then
        write (io%uout, '(a,i0,f16.8,1x,f16.8,1x,f16.8)') "check parameter ", i, bounds (i, 1), p (1, i), bounds (i, 2)
        call error ("the starting value (x) is outside the bounds execution terminated without any"//" optimization. lb(i) < x(i) <&
       &ub(i), i = 1, n.", myname, .true., io)
      end if
    end do
!
    do i = 2, gen%fit%iNoParams + 1
      do j = 1, gen%fit%iNoParams
        if (i-1 == j) then
          p (i, j) = bounds (j, 1) + Abs (bounds(j, 2)-bounds(j, 1)) * ranmar (sol%seed)
!generate a  new vertex in the box           
        else
          p (i, j) = p (1, j)
        end if
      end do
    end do
!
    do i = 1, gen%fit%iNoParams + 1
      y (i) = UpdateCost (p(i, :), gen, atomic, tb, sol, io)
    end do
  end subroutine InitSimplex
!
!> \brief initializes the simulated annealing method
!> \author Alin M Elena
!> \date 14/11/07, 13:35:55
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine InitSA (gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "InitSA"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
    real (k_pr), allocatable :: x (:), xopt (:), c (:), vm (:), fstar (:), xp (:)
    real (k_pr) :: t, eps, rt, fopt
    integer, allocatable :: nacp (:)
    integer :: ns, nt, nfcnev, ier, maxevl, nacc, nobds, i, neps, errno
    logical :: Max
!
    allocate (x(1:gen%fit%iNoParams), xopt(1:gen%fit%iNoParams), c(1:gen%fit%iNoParams), vm(1:gen%fit%iNoParams), &
   & xp(1:gen%fit%iNoParams), nacp(1:gen%fit%iNoParams))
!
!  set input parameters.      
    Max = .false.
    neps = gen%fit%neps
    eps = gen%fit%fitTol
    rt = gen%fit%rt
    ns = gen%fit%ns
    nt = gen%fit%nt
    maxevl = gen%fit%feval
    allocate (fstar(1:neps))
    do i = 1, gen%fit%iNoParams
      c (i) = gen%fit%stepAd
    end do
    call ReadBounds (gen, io)
    if (gen%fit%RestartFit) then
      open (2, file=trim(cRestart), status="old", action="read", iostat=errno)
      if (errno /= 0) then
        call error ("restart file not found!!", myname, .true., io)
      end if
      do i = 1, gen%fit%iNoParams
        read (2,*) x (i)
      end do
      close (2)
      call UpdateParams (x, gen, atomic, tb, sol, io)
    else
      call SetInitialParams (x, atomic, tb)
    end if
    t = gen%fit%temp
    do i = 1, gen%fit%iNoParams
      vm (i) = gen%fit%step
    end do
    write (io%uout, 1000) gen%fit%iNoParams, Max, t, rt, eps, ns, nt, neps, maxevl
1000 format (/, ' simulated annealing ', /, /, ' number of parameters: ', i3, '   maximazation: ', l5, /, ' initial temp: ', g8.2, &
   & '   rt: ', g8.2, '   eps: ', g8.2, /, ' ns: ', i3, '   nt: ', i2, '   neps: ', i2, /, ' maxevl: ', i10)
    call PrintVectorA (x, 'starting values', .true., .false., io)
    call PrintVectorA (vm, 'initial step length', .true., .false., io)
    call PrintVectorA (bounds(:, 1), 'lower bound', .true., .false., io)
    call PrintVectorA (bounds(:, 2), 'upper bound', .true., .false., io)
    call PrintVectorA (c, 'c vector', .true., .false., io)
    write (io%uout, '(a)') "/  ****   end of driver routine output   **** /"
    write (io%uout, '(a)') "****  before call to SimulAnnealing.  ****"
    call SimulAnnealing (gen%fit%iNoParams, x, Max, rt, eps, ns, nt, neps, maxevl, bounds(:, 1), bounds(:, 2), c, t, vm, xopt, &
   & fopt, nacc, nfcnev, nobds, ier, fstar, xp, nacp, UpdateCost, gen, atomic, sol, tb, io)
    write (io%uout, '(/,''  ****   results after sa   ****   '')')
    call PrintVectorA (xopt, 'solution', .true., .false., io)
    call PrintVectorA (vm, 'final step length', .true., .false., io)
    write (io%uout, 1001) fopt, nfcnev, nacc, nobds, t, ier
1001 format (/, ' optimal function value: ', g20.13 /, ' number of function evaluations:     ', i10, /, ' number of accepted evalua&
   &tions:     ', i10, /, ' number of out of bound evaluations: ', i10, /, ' final temp: ', g20.13, '  ier: ', i3)
    best = xopt
    call UpdateParams (xopt, gen, atomic, tb, sol, io)
    open (2, file=trim(cRestart), status='unknown', action="write")
    do i = 1, gen%fit%iNoParams
      write (2,*) x (i)
    end do
    close (2)
!     call print_gsp(trim(cFdfFile))
    deallocate (fstar, x, xopt, c, vm, xp, nacp)
  end subroutine InitSA
!
!> \brief initializes the simplex-simulated annealing method
!> \author Alin M Elena
!> \date 14/11/07, 13:48:55
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine SimplexSA (gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "SimplexSA"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
    real (k_pr) :: opt (1:gen%fit%iNoParams), yb, temperature, copt
    integer :: i, iter, k
    logical :: quit
    real (k_pr), allocatable :: tol (:)
!init simplex        
    call InitSimplex (gen, atomic, tb, sol, io)
    yb = huge (yb)
    allocate (tol(1:gen%fit%neps))
    do i = 1, gen%fit%neps
      tol (i) = 1.0e20_k_pr
    end do
!end init simplex        
    temperature = gen%fit%temp
    write (io%uout,*) "simulated annealing stage"
    write (io%uout,*) "initial values:"
    call PrintVectorA (p(1, :), "initial parameters", .true., .false., io)
    write (io%uout,*) "starting temperature: ", gen%fit%temp
    write (io%uout,*) "maximum simplex iterations: ", gen%fit%iter
    write (io%uout,*) "fit tolerance: ", gen%fit%fitTol
    write (io%uout,*) " temperature reduction factor: ", gen%fit%rt
    copt = UpdateCost (p(1, :), gen, atomic, tb, sol, io)
    write (io%uout,*) "initial cost function: ", copt
    best (:) = p (1, :)
    k = 0
    do
      iter = gen%fit%iter
      k = k + 1
      write (io%uout, '(a,i0)') "Enter simplex step: ", k
      write (io%uout, '(a,f16.8)') "current temperature: ", temperature
      call amebsa (p, y, opt, yb, gen%fit%fitTol, UpdateCost, iter, temperature, bounds, gen, atomic, tb, sol, io)
      write (io%uout, '(a,i0)') "iterations in simplex step: ", iter
      do i = gen%fit%neps, 2, - 1
        tol (i) = tol (i-1)
      end do
      tol (1) = yb
      call PrintVectorA (opt, 'current parameters:', .true., .false., io)
      quit = .false.
      if (Abs(copt-tol(1)) < gen%fit%fitTol) then
        quit = .true.
      end if
      if (copt > yb) then
        best (:) = opt (:)
        copt = yb
      end if
      write (io%uout,*) "current  ", yb, "and best cost functions", copt
      do i = 2, gen%fit%neps
        if (Abs(tol(1)-tol(i)) > gen%fit%fitTol) then
          quit = .false.
        end if
      end do
      call PrintVectorA (best, 'best parameters so far:', .true., .false., io)
      if (quit) then
        call PrintVectorA (opt, 'optimal parameters:', .true., .false., io)
        exit
      end if
      temperature = temperature * gen%fit%rt
      if (temperature < epsilon(temperature)) then
        write (io%uout,*) "temperature under machine precision!!!"
        exit
      end if
    end do
    deallocate (tol)
  end subroutine SimplexSA
!
!> \brief initializes the trust-region method
!> \author Alin M Elena
!> \date 14/11/07, 13:48:55
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine TrustRegion (gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "TrustRegion"
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    integer :: n !! n - number of function variables
    integer :: m !! m - dimension of function value
    real (k_pr) :: eps (1:6)!! precisions for stop-criteria (see manual for more details)
    real (k_pr) :: jac_eps !! jacobi calculation precision
    integer :: rci_request !! reverse communication interface parameter
    real (k_pr), allocatable :: fvec (:)!! function (f(x)) value vector
    real (k_pr), allocatable :: fjac (:, :)!! jacobi matrix
    integer :: iter !! number of iterations
    integer :: st_cr !! number of stop-criterion
    integer :: successful !! controls of rci cycle
    integer :: iter1 !! maximum number of iterations
    integer :: iter2 !! maximum number of iterations of calculation of trial-step
    real (k_pr) :: rs !! initial step bound
    real (k_pr) :: r1, r2 !! initial and final residuals
    integer (kind=8) :: handle !! tr solver handle
    integer :: i, j, errno, un
    real (k_pr), allocatable :: lw (:), up (:)
    real (k_pr), allocatable :: f1 (:), f2 (:)
!
!
    n = gen%fit%iNoParams
    m = fitData%n - 1
    allocate (lw(1:n), up(1:n), fvec(1:m), fjac(1:m, 1:n))
    allocate (f1(1:m), f2(1:m))
!! set precisions for stop-criteria    
    eps (1:6) = gen%fit%fitTol
!! set maximum number of iterations    
    iter1 = gen%fit%iter
!! set maximum number of iterations of calculation of trial-step    
    iter2 = iter1 / 100
!! set initial step bound    
    rs = gen%fit%step
!! precisions for jacobi calculation    
    jac_eps = gen%fit%fitTol
!! set the initial guess    
    if (gen%fit%RestartFit) then
      un = GetUnit ()
      open (un, file=trim(cRestart), status="old", action="read", iostat=errno)
      if (errno /= 0) then
        call error ("restart "//trim(cRestart)//" file not found!!", myname, .true., io)
      end if
      do i = 1, gen%fit%iNoParams
        read (un,*) best (i)
      end do
      close (un)
      call UpdateParams (best, gen, atomic, tb, sol, io)
    else
      call SetInitialParams (best, atomic, tb)
    end if
!! set lower and upper bounds    
    call ReadBounds (gen, io)
    lw (:) = bounds (:, 1)
    up (:) = bounds (:, 2)
!! set initial values    
    do i = 1, m
      fvec (i) = 0.0_k_pr
      do j = 1, n
        fjac (i, j) = 0.0_k_pr
      end do
    end do
!! initialize solver (allocate memory, set initial values)    
!! handle in/out: tr solver handle    
!! n in: number of function variables    
!! m in: dimension of function value    
!! x in: solution vector. contains values x for f(x)    
!! bounds(:,1) in: lower bound    
!! bounds(:,2) in: upper bound    
!! eps in: precisions for stop-criteria    
!! iter1 in: maximum number of iterations    
!! iter2 in: maximum number of iterations of calculation of trial-step    
!! rs in: initial step bound    
    if (dtrnlspbc_init(handle, n, m, best, lw, up, eps, iter1, iter2, rs) /= tr_success) then
!! if function does not complete successfully then print error message    
      call error ("Initialization of optimizer failed!!", myname, .true., io)
    end if
!! set initial rci cycle variables    
    rci_request = 0
    successful = 0
!! rci cycle    
    do while (successful ==  0)
!! call tr solver    
!! handle in/out: tr solver handle    
!! fvec in: vector    
!! fjac in: jacobi matrix    
!! rci_request in/out: return number that denotes next step for performing    
      if (dtrnlspbc_solve(handle, fvec, fjac, rci_request) /= tr_success) then
!! if function does not complete successfully then print error message      
        call error ("Optimizer failed!!", myname, .true., io)
      end if
!! rci_request in/out: return number that denotes next step for performing      
!! according to rci_request value we do next step      
      select case (rci_request)
      case (-1,-2,-3,-4,-5,-6)
!! stop rci cycle          
        successful = 1
      case (1)
!! recalculate function value        
!! m in: dimension of function value        
!! n in: number of function variables        
!! best in: solution vector        
!! fvec out: function value f(x)        
        call ObjectiveFunctions (m, n, best, fvec, gen, atomic, tb, sol, io)
      case (2)
!! compute jacobi matrix        
!! ObjectiveFunctions in: external objective function        
!! n in: number of function variables        
!! m in: dimension of function value        
!! fjac out: jacobi matrix        
!! x in: solution vector        
!! jac_eps in: jacobi calculation precision        
!
        if (ComputeJacobi(n, m, fjac, best, jac_eps, f1, f2, gen, atomic, tb, sol, io) /= tr_success) then
!! if function does not complete successfully then print error message        
          call error ("Jacobian calculation failed!!", myname, .true., io)
        end if
      end select
    end do
!! get solution statuses    
!! handle in: tr solver handle    
!! iter out: number of iterations    
!! st_cr out: number of stop criterion    
!! r1 out: initial residuals    
!! r2 out: final residuals    
    if (dtrnlspbc_get(handle, iter, st_cr, r1, r2) /= tr_success) then
!! if function does not complete successfully then print error message    
      call error ("getting status of optimizer failed!!", myname, .true., io)
    end if
!
!! free handle memory    
    if (dtrnlspbc_delete(handle) /= tr_success) then
!! if function does not complete successfully then print error message      
      call error ("Cleaning the memory of optimizer failed!!", myname, .true., io)
    end if
!
!! if final residual is less than required precision then print pass    
    if (r2 < gen%fit%fitTol) then
      write (io%uout, '(a,i0,a,f16.8)') "Optimizer was successfull in ", iter, " iterations with a final residual of", r2
      write (io%uout, '(a,i0)') "Termination criteria ", st_cr
!! else print failed      
    else
      write (io%uout, '(a,f16.8)') "Optimizer failed with a final residual of ", r2
      write (io%uout, '(a,i0)') "Termination criteria ", st_cr
    end if
!
    deallocate (lw, up, fvec, fjac)
    deallocate (f1, f2)
  end subroutine TrustRegion
!
!
  subroutine ObjectiveFunctions (m, n, x, f, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "ObjectiveFunctions"
    integer, intent (inout) :: n, m
    real (k_pr), intent (inout) :: x (:), f (:)
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
    integer :: i, ilevels
!
    call UpdateParams (x, gen, atomic, tb, sol, io)
    call SinglePoint (io, gen, atomic, tb, sol)
    if ( .not. gen%lIsSCFConverged) then
      f (1:m) = k_infinity
    else
      ilevels = sol%h%dim
      fitData%fit (1:ilevels) = sol%eigenvals(1:ilevels)
!
      do i = 1, ilevels / 2 - 1
        f (i) = fitData%exper(i) - fitData%exper(i+1) - (sol%eigenvals(i)-sol%eigenvals(i+1))
      end do
!
      do i = 1 + ilevels / 2, ilevels - 1
        f (i-1) = fitData%exper(i) - fitData%exper(i+1) - (sol%eigenvals(i)-sol%eigenvals(i+1))
      end do
      f (ilevels-1) = fitData%exper(1+ilevels/2) - fitData%exper(1) - (sol%eigenvals(1+ilevels/2)-sol%eigenvals(1))
! add charges in the objective functions  
      fitData%fit (ilevels+1:ilevels+atomic%atoms%natoms) = atomic%atoms%chrg(1:atomic%atoms%natoms)
      do i = 1, atomic%atoms%natoms
        f (i+ilevels-1) = fitData%exper(i+ilevels) - atomic%atoms%chrg(i)
      end do
    end if
  end subroutine ObjectiveFunctions
!
  integer function ComputeJacobi (n, m, fjac, x, jac_eps, f1, f2, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "ComputeJacobi"
    integer, intent (inout) :: n, m
    real (k_pr), intent (inout) :: x (:), fjac (:, :), f1 (:), f2 (:), jac_eps
    integer (kind=8) :: handle
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
    integer :: successful !! controls of rci cycle
    integer :: rci_request !! reverse communication interface parameter
!
    rci_request = 0
    successful = 0
    ComputeJacobi = tr_success - 1
!
    if (djacobi_init(handle, n, m, x, fjac, jac_eps) /= tr_success) then
! if function does not complete successfully then print error message    
      call error ("allocating memory for the jacobian failed!!", myname, .true., io)
    end if
!
!! rci cycle    
    do while (successful ==  0)
!! call solver    
      if (djacobi_solve(handle, f1, f2, rci_request) /= tr_success) then
!! if function does not complete successfully then print error message      
        call error ("computation of jacobian failed!!", myname, .true., io)
      end if
      if (rci_request == 1) then
!! calculate function value f1 = f(x+eps)      
        call ObjectiveFunctions (m, n, x, f1, gen, atomic, tb, sol, io)
      else if (rci_request == 2) then
!! calculate function value f1 = f(x-eps)        
        call ObjectiveFunctions (m, n, x, f2, gen, atomic, tb, sol, io)
      else if (rci_request == 0) then
!! exit rci cycle        
        successful = 1
      end if
    end do
    if (djacobi_delete(handle) /= tr_success) then
!! if function does not complete successfully then print error message    
      call error ("Cleaning the memory of jacobian failed!!", myname, .true., io)
    end if
    ComputeJacobi = tr_success
  end function ComputeJacobi
!
end module m_Fit
