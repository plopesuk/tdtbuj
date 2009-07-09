!> \brief implements bfgs method for large scale optimization
!> \details (Numerical Recipes 3rd edition chapter 10)
!> \author Alin M. Elena
!> \date 18th of January 2008
module m_BFGS
  use m_Constants
  use m_Useful
  use m_Types
  use m_Gutenberg
#ifdef MKL95
  use mkl95_BLAS, only: dot
#endif
  implicit none
!
  private
!
  public :: driverBFGS
!
contains
!
!
!> \brief driver routine for BFGS optimization
!> \author Alin M Elena
!> \date 30/05/08, 15:22:53
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param func function the function that evaluates the cost
  subroutine driverBFGS (io, gen, atomic, tb, sol, func)
    character (len=*), parameter :: myname = 'driverBFGS'
    type (ioType), intent (inout) :: io
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (modelType), intent (inout) :: tb
!
    interface
      function func (gen, atomic, tb, sol, io, x, f, gradient)
        use m_Constants
        use m_Types
        implicit none
        real (k_pr), dimension (:), intent (in) :: x
        real (k_pr), dimension (:), intent (inout), optional :: gradient
        real (k_pr), intent (inout) :: f
        type (generalType), intent (inout) :: gen
        type (atomicxType), intent (inout) :: atomic
        type (modelType), intent (inout) :: tb
        type (solutionType), intent (inout) :: sol
        type (ioType), intent (inout) :: io
        integer :: func
      end function func
    end interface
!
    integer :: i
    integer :: m, n, res, itermax
    real (k_pr) :: eps, epsg, gtol, xtol, ftol, stpmax
    real (k_pr), allocatable :: x (:)
    integer :: info, atom
!
!
!
!
    write (io%uout, '(/a/)') "--BFGS Optimization Run-----------------------------------------"
    n = atomic%atoms%nmoving * 3
    eps = gen%ftol
    xtol = gen%xtol
    gtol = gen%gtol
    itermax = gen%nsteps
    stpmax = 100.0_k_pr
    info = 0
!
    allocate (x(1:n))
    do i = 1, atomic%atoms%nmoving
      atom = atomic%atoms%id(atomic%atoms%moving(i))
      x (3*(i-1)+1) = atomic%atoms%x(atom)
      x (3*(i-1)+2) = atomic%atoms%y(atom)
      x (3*(i-1)+3) = atomic%atoms%z(atom)
    end do
!
    call bfgs (n, x, itermax, stpmax, eps, xtol, gtol, info, func, gen, atomic, tb, sol, io)
    res = func (gen, atomic, tb, sol, io, x, ftol, x)
    deallocate (x)
!
    open (unit=1, file="new_coords.xyz", status="unknown", action="write")
    call PrintXYZ (1, atomic, .false., "Optimised coordinates")
    close (unit=1)
!
    write (io%uout, '(/a)') "Final forces:"
    call PrintForces (atomic%atoms, io)
    write (io%uout, '(/a)') "Final coordinates:"
    call PrintCoordinates (io, atomic, gen)
    if (info <= 0) then
      write (io%uout, '(/a,i4)') "WARNING: Optimization did not converge.", info
    else
      write (io%uout, '(/a,i4)') "BFGS report: Optimization converged with code: ", info
    end if
! !      call write_fdf_coords()
    write (io%uout, '(/a)') "----------------------------------------------------------------"
  end subroutine driverBFGS
!
!
!
!> \remarks code for info at exit
!> info = -2 error in lineSearch
!>        -1 unexpected error
!>         0 maximum number of iterations has been reached
!>         1 xtol criteria has been satisfied
!>         2 gtol criteria has been satisfied
!
  subroutine bfgs (n, x, itermax, stpmx, eps, xtol, gtol, info, func, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "BFGS"
    integer, intent (inout) :: n, itermax, info
    real (k_pr), intent (inout) :: stpmx, eps, gtol, xtol
    real (k_pr), dimension (:), intent (inout) :: x
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    interface
      function func (gen, atomic, tb, sol, io, x, f, gradient)
        use m_Constants
        use m_Types
        implicit none
        real (k_pr), dimension (:), intent (in) :: x
        real (k_pr), dimension (:), intent (inout), optional :: gradient
        real (k_pr), intent (inout) :: f
        type (generalType), intent (inout) :: gen
        type (atomicxType), intent (inout) :: atomic
        type (modelType), intent (inout) :: tb
        type (solutionType), intent (inout) :: sol
        type (ioType), intent (inout) :: io
        integer :: func
      end function func
    end interface
!
!! local variables now
    real (k_pr) :: den, fac, fad, fae, fp, stpmax, sum, sumdg, sumxi, temp, test, fret, f
    real (k_pr), allocatable :: dg (:), g (:), hdg (:), xnew (:), xi (:), hessin (:, :)
    integer :: i, j, iters, res, infoLine
    allocate (dg(1:n), g(1:n), hdg(1:n), xnew(1:n), xi(1:n), hessin(1:n, 1:n))
    sum = 0.0_k_pr
    info = - 1
    res = func (gen, atomic, tb, sol, io, x, fp, g)
    hessin = 0.0_k_pr
    do i = 1, n
      hessin (i, i) = 1.0_k_pr
      xi (i) = - g (i)
      sum = sum + x (i) * x (i)
    end do
    stpmax = stpmx * Max (Sqrt(sum), real(n, k_pr))
    do iters = 1, itermax
      call lineSearch (x, fp, g, xi, xnew, fret, stpmax, n, xtol, infoLine, func, gen, atomic, tb, sol, io)
      if (infoLine < 0) then
        deallocate (dg, g, hdg, xnew, xi, hessin)
        info = - 2
        return
      end if
      fp = fret
      do i = 1, n
        xi (i) = xnew (i) - x (i)
        x (i) = xnew (i)
      end do
      test = 0.0_k_pr
      do i = 1, n
        temp = Abs (xi(i)) / Max (Abs(x(i)), 1.0_k_pr)
        if (temp > test) then
          test = temp
        end if
      end do
      if (test < xtol) then
        deallocate (dg, g, hdg, xnew, xi, hessin)
        info = 1
        return
      end if
      dg = g
      res = func (gen, atomic, tb, sol, io, x, f, g)
      test = 0.0_k_pr
      den = Max (fret, 1.0_k_pr)
      do i = 1, n
        temp = Abs (g(i)) * Max (Abs(x(i)), 1.0_k_pr) / den
        if (temp > test) then
          test = temp
        end if
      end do
      if (test < gtol) then
        info = 2
        deallocate (dg, g, hdg, xnew, xi, hessin)
        return
      end if
      dg = g - dg
      do i = 1, n
        hdg (i) = 0.0_k_pr
        do j = 1, n
          hdg (i) = hdg (i) + hessin (i, j) * dg (j)
        end do
      end do
      fac = 0.0_k_pr
      fae = 0.0_k_pr
      sumdg = 0.0_k_pr
      sumxi = 0.0_k_pr
      do i = 1, n
        fac = fac + dg (i) * xi (i)
        fae = fae + dg (i) * hdg (i)
        sumdg = sumdg + dg (i) * dg (i)
        sumxi = sumxi + xi (i) * xi (i)
      end do
      if (fac > Sqrt(eps*sumdg*sumxi)) then
        fac = 1.0_k_pr / fac
        fad = 1.0_k_pr / fae
        do i = 1, n
          dg (i) = fac * xi (i) - fad * hdg (i)
        end do
        do i = 1, n
          do j = i, n
            hessin (i, j) = hessin (i, j) + fac * xi (i) * xi (j) - fad * hdg (i) * hdg (j) + fae * dg (i) * dg (j)
            hessin (j, i) = hessin (i, j)
          end do
        end do
      end if
      do i = 1, n
        xi (i) = 0.0_k_pr
        do j = 1, n
          xi (i) = xi (i) - hessin (i, j) * g (j)
        end do
      end do
    end do
    info = 0
    call error ("Maximum number of iterations has been reached. It may be that the solution found is not the best!!!", myname, &
   & .false., io)
    deallocate (dg, g, hdg, xnew, xi, hessin)
  end subroutine bfgs
!
  subroutine lineSearch (xold, fold, g, p, x, f, stpmax, n, xtol, info, func, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "lineSearch"
    integer, intent (inout) :: n, info
    real (k_pr), intent (inout) :: stpmax, xtol, fold, f
    real (k_pr), dimension (:), intent (inout) :: xold, g, p, x
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    interface
      function func (gen, atomic, tb, sol, io, x, f, gradient)
        use m_Constants
        use m_Types
        implicit none
        real (k_pr), dimension (:), intent (in) :: x
        real (k_pr), dimension (:), intent (inout), optional :: gradient
        real (k_pr), intent (inout) :: f
        type (generalType), intent (inout) :: gen
        type (atomicxType), intent (inout) :: atomic
        type (modelType), intent (inout) :: tb
        type (solutionType), intent (inout) :: sol
        type (ioType), intent (inout) :: io
        integer :: func
      end function func
    end interface
!
!!! internal variables
    real (k_pr) :: alf, a, alam, alam2, alamin, b, disc, f2
    real (k_pr) :: rhs1, rhs2, slope, sum, temp, test, tmplam
    integer :: i, res
    logical :: check
#ifndef MKL95
    real(k_pr),external :: ddot
#endif
    alf = 1.0e-4_k_pr
    info = 0
    alam2 = 0.0_k_pr
    f2 = 0.0_k_pr
    slope = 0.0_k_pr
    sum = 0.0_k_pr
    check = .false.
#ifdef MKL95
    sum = Sqrt (dot(p, p))
#else
    sum = Sqrt (ddot(n,p,k_ione, p,k_ione))
#endif
    if (sum > stpmax) then
      p = p * stpmax / sum
    end if
#ifdef MKL95
    slope = dot (g, p)
#else
    slope = ddot (n,g, k_ione,p,k_ione)
#endif
    if (slope >= 0.0_k_pr) then
      call error ("Roundoff problem in "//trim(myname), myname, .false., io)
      info = - 1
      return
    end if
    test = 0.0_k_pr
    do i = 1, n
      temp = Abs (p(i)) / Max (Abs(xold(i)), 1.0_k_pr)
      if (temp > test) then
        test = temp
      end if
    end do
    alamin = xtol / test
    alam = 1.0_k_pr
    do
      x = xold + alam * p
      res = func (gen, atomic, tb, sol, io, x, f)
      if (alam < alamin) then
        x = xold
        check = .true.
        info = 1
        return
      else if (f <= fold+alf*alam*slope) then
        return
      else
        if (Abs(alam-1.0_k_pr) <= tiny(1.0_k_pr)) then
          tmplam = - slope / (2.0_k_pr*(f-fold-slope))
        else
          rhs1 = f - fold - slope * alam
          rhs2 = f2 - fold - slope * alam2
          a = (rhs1/(alam*alam)-rhs2/(alam2*alam2)) / (alam-alam2)
          b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2)) / (alam-alam2)
          if (Abs(a) <= tiny(1.0_k_pr)) then
            tmplam = - slope / (2.0_k_pr*b)
          else
            disc = b * b - 3.0_k_pr * a * slope
            if (disc < 0.0_k_pr) then
              tmplam = 0.5_k_pr * alam
            else if (b <= 0.0_k_pr) then
              tmplam = (-b+Sqrt(disc)) / (3.0_k_pr*a)
            else
              tmplam = - slope / (b+Sqrt(disc))
            end if
          end if
          if (tmplam > 0.5_k_pr*alam) then
            tmplam = 0.5_k_pr * alam
          end if
        end if
      end if
      alam2 = alam
      f2 = f
      alam = Max (tmplam, 0.1_k_pr*alam)
    end do
  end subroutine lineSearch
end module m_BFGS
