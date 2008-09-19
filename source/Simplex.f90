!> \brief implements simplex algorithm for minimization in the spiriti of NR
!> \author Alin M Elena
!> \date 12th of November 2007
module m_Simplex
  use m_Constants
  use m_Types
  use m_Useful
  implicit none
  private
!
  public :: amoeba
!
contains
!> \brief minimization of the function func in n dimensions by the downhill simplex method
!> \details of nelder and mead. the (n + 1)Xn matrix p is input. its n + 1 rows are n-dimensional
!> vectors that are the vertices of the starting simplex. also input is the vector y of length
!> n + 1, whose components must be preinitialized to the values of func evaluated at the
!> n + 1 vertices (rows) of p;and ftol the fractional convergence tolerance to be achieved
!> in the function value (n.b.!). on output, p and y will have been reset to n +1 new points
!> all within ftol of a minimum function value, and iter gives the number of function
!> evaluations taken.
!> \date 12th of November 2007
!> \param p matrix n+1 by n contains the initial simplex
!> \param y array, n+1, contains the cost function for each point in simplex
!> @param func function the cost function
!> \param bounds array array, n, the parameters to be optimized
!> \param iter integer number of iterations done
!> \param itmax integer maximum number of iterations
!> \param ftol real the desired tolerance
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!> \remarks all the job is done by AmoebaPrivate
!> \internal add citation to the original paper
  subroutine amoeba (p, y, ftol, func, iter, bounds, itmax, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "Amoeba"
    integer, intent (out) :: iter
    real (k_pr), intent (in) :: ftol
    real (k_pr), dimension (:), intent (inout) :: y
    real (k_pr), dimension (:, :), intent (inout) :: p
    real (k_pr), dimension (:, :), intent (in) :: bounds
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    integer, intent (in) :: itmax
    interface
      function func (x, gen, atomic, tb, sol, io)
        use m_Constants
        use m_Types
        implicit none
        real (k_pr), dimension (:), intent (in) :: x
        type (generalType), intent (inout) :: gen
        type (atomicxType), intent (inout) :: atomic
        type (modelType), intent (inout) :: tb
        type (solutionType), intent (inout) :: sol
        type (ioType), intent (inout) :: io
        real (k_pr) :: func
      end function func
    end interface
!
    real (k_pr), parameter :: tiny = 1.0e-10
    integer :: ihi, ndim !global variables.
    real (k_pr), dimension (size(p, 2)) :: psum
    call AmoebaPrivate
!
  contains
!
    subroutine AmoebaPrivate
      implicit none
      character (len=*), parameter :: myname = "AmoebaPrivate"
      integer :: i, ilo, inhi, iii
      real (k_pr) :: rtol, ysave, ytry, ytmp
      ndim = AssertEq (size(p, 2), size(p, 1)-1, size(y)-1, io)
      iter = 0
      psum (:) = sum (p(:, :), dim=1)
      do
        ilo = iminloc (y(:))!determine which point is the highest (worst),
        ihi = imaxloc (y(:))!next-highest, and lowest (best).
        ytmp = y (ihi)
        y (ihi) = y (ilo)
        inhi = imaxloc (y(:))
        y (ihi) = ytmp
        rtol = 2.0_k_pr * Abs (y(ihi)-y(ilo)) / (Abs(y(ihi))+Abs(y(ilo))+tiny)
!  rtol=abs(y(ihi)-y(ilo))        
!compute the fractional range from highest to lowest and return if satisfactory.  
        if (rtol < ftol) then !if returning, put best point and value in slot 1
          call swap (y(1), y(ilo))
          call swap (p(1, :), p(ilo, :))
          return
        end if
        if (iter >= itmax) then
          call error ("itmax achieved in amoeba", myname, .false., io)
          return
        end if
! begin a new iteration. first extrapolate by a factor ?1 through the face of the simplex  
! across from the high point, i.e., reflect the simplex from the high point.  
        ytry = amotry (-1.0_k_pr)
        iter = iter + 1
        if (ytry <= y(ilo)) then
!gives a result better than the best point, so  
!try an additional extrapolation by a factor of 2.  
          ytry = amotry (2.0_k_pr)
          iter = iter + 1
        else if (ytry >= y(inhi)) then
!the reflected point is worse than the second  
!highest, so look for an intermediate lower point, i.e., do a one-dimensional contraction.  
          ysave = y (ihi)
          ytry = amotry (0.5_k_pr)
          iter = iter + 1
          if (ytry >= ysave) then
!can't seem to get rid of that high point. better contract around the lowest (best) point.  
            p (:, :) = 0.5_k_pr * (p(:, :)+spread(p(ilo, :), 1, size(p, 1)))
            do i = 1, ndim + 1
              if (i /= ilo) then
                y (i) = func (p(i, :), gen, atomic, tb, sol, io)
!!! check the bounds                
                do iii = 1, ndim
                  if ((p(i, iii) < bounds(iii, 1)) .or. (p(i, iii) > bounds(iii, 2))) then
                    y (i) = k_infinity
                  end if
                end do
              end if
            end do
            iter = iter + ndim !keep track of function evaluations.
            psum (:) = sum (p(:, :), dim=1)
          end if
        end if
      end do !go back for the test of doneness and the next
    end subroutine AmoebaPrivate
!> \brief generates a new point in the simplex
!> \details extrapolates by a factor fac through the face of the simplex across from the high point,
!> tries it, and replaces the high point if the new point is better.
!> \author Alin M Elena
!> \date 12th of November 2007
!> \param fac real the factpor used in extrapolation
    function amotry (fac)
      character (len=*), parameter :: myname = "AmoTry"
      real (k_pr), intent (in) :: fac
      real (k_pr) :: amotry
      integer :: ii
!
      real (k_pr) :: fac1, fac2, ytry
      real (k_pr), dimension (size(p, 2)) :: ptry
      fac1 = (1.0_k_pr-fac) / ndim
      fac2 = fac1 - fac
      ptry (:) = psum (:) * fac1 - p (ihi, :) * fac2
      ytry = func (ptry, gen, atomic, tb, sol, io)!evaluate the function at the trial point.
!!! check the bounds  
      do ii = 1, ndim
        if ((ptry(ii) < bounds(ii, 1)) .or. (ptry(ii) > bounds(ii, 2))) then
          ytry = k_infinity
        end if
      end do
      if (ytry < y(ihi)) then !if it's better than the highest, then replace
        y (ihi) = ytry !the highest.
        psum (:) = psum (:) - p (ihi, :) + ptry (:)
        p (ihi, :) = ptry (:)
      end if
      amotry = ytry
    end function amotry
  end subroutine amoeba
end module m_Simplex
