!> \brief Simplex-Simulated Annealing in the spirit of NR
!> \author Alin M Elena
!> \date 12/11/07, 13:38:02
!
module m_SimplexSA
  use m_Constants
  use m_Types
  use m_Useful
  implicit none
!
  private
!
  public :: amebsa
contains
!
!> \brief Simplex-Simulated Annealing method
!> \author Alin M Elena
!> \date 12/11/07, 13:39:28
!> \param p matrix n+1 by n contains the initial simplex
!> \param y array, n+1, contains the cost function for each point in simplex
!> @param func function the cost function
!> \param bounds array array, n, the parameters to be optimized
!> \param iter integer number of iterations done
!> \param ftol real the desired tolerance
!> \param temptr the temperature
!> \param pb the best point
!> \param yb the cost of the best point
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!> \remarks all the job is done by AmebsaPrivate
!> \internal add citation to the original paper
  subroutine amebsa (p, y, pb, yb, ftol, func, iter, temptr, bounds, gen, atomic, tb, sol, io)
    integer, intent (inout) :: iter
    real (k_pr), intent (inout) :: yb
    real (k_pr), intent (inout) :: ftol, temptr
    real (k_pr), dimension (:), intent (inout) :: y, pb
    real (k_pr), dimension (:, :), intent (inout) :: p
    real (k_pr), dimension (:, :), intent (in) :: bounds
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
!
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
    integer :: ihi, ndim
    real (k_pr) :: yhi
    real (k_pr), dimension (size(p, 2)) :: psum
!
    call AmebsaPrivate
  contains
!
    subroutine AmebsaPrivate
      integer :: i, iii, j, ilo, inhi
      real (k_pr) :: rtol, ylo, ynhi, ysave, ytry
      real (k_pr), dimension (size(y)) :: yt, harvest
      ndim = AssertEq (size(p, 2), size(p, 1)-1, size(y)-1, size(pb), io)
      psum (:) = sum (p(:, :), dim=1)
      do
        do j = 1, ndim + 1
          harvest (j) = ranmar (sol%seed)
        end do
        yt (:) = y (:) - temptr * Log (harvest)
        ilo = iminloc (yt(:))
        ylo = yt (ilo)
        ihi = imaxloc (yt(:))
        yhi = yt (ihi)
        yt (ihi) = ylo
        inhi = imaxloc (yt(:))
        ynhi = yt (inhi)
        rtol = 2.0_k_pr * Abs (yhi-ylo) / (Abs(yhi)+Abs(ylo))
        if (rtol < ftol .or. iter < 0) then
          call swap (y(1), y(ilo))
          call swap (p(1, :), p(ilo, :))
          return
        end if
        ytry = amotsa (-1.0_k_pr)
        iter = iter - 1
        if (ytry <= ylo) then
          ytry = amotsa (2.0_k_pr)
          iter = iter - 1
        else if (ytry >= ynhi) then
          ysave = yhi
          ytry = amotsa (0.5_k_pr)
          iter = iter - 1
          if (ytry >= ysave) then
            p (:, :) = 0.5_k_pr * (p(:, :)+spread(p(ilo, :), 1, size(p, 1)))
            do i = 1, ndim + 1
              if (i /= ilo) then
                y (i) = func (p(i, :), gen, atomic, tb, sol, io)
!check the bounds
                do iii = 1, ndim
                  if ((p(i, iii) < bounds(iii, 1)) .or. (p(i, iii) > bounds(iii, 2))) then
                    y (i) = k_infinity
                  end if
                end do
              end if
            end do
            iter = iter - ndim
            psum (:) = sum (p(:, :), dim=1)
          end if
        end if
      end do
    end subroutine AmebsaPrivate
!> \brief extrapolates a point an checks it
!> \author Alin M Elena
!> \date 12/11/07, 13:35:59
!> \param fac real the factor used for extrapolation
    function amotsa (fac)
      real (k_pr), intent (in) :: fac
      real (k_pr) :: amotsa
      real (k_pr) :: fac1, fac2, yflu, ytry, harv
      real (k_pr), dimension (size(p, 2)) :: ptry
      integer :: i
      fac1 = (1.0_k_pr-fac) / ndim
      fac2 = fac1 - fac
      ptry (:) = psum (:) * fac1 - p (ihi, :) * fac2
      ytry = func (ptry, gen, atomic, tb, sol, io)
      do i = 1, ndim
        if ((ptry(i) < bounds(i, 1)) .or. (ptry(i) > bounds(i, 2))) then
          ytry = k_infinity
        end if
      end do
      if (ytry <= yb) then
        pb (:) = ptry (:)
        yb = ytry
      end if
      harv = ranmar (sol%seed)
      yflu = ytry + temptr * Log (harv)
      if (yflu < yhi) then
        y (ihi) = ytry
        yhi = yflu
        psum (:) = psum (:) - p (ihi, :) + ptry (:)
        p (ihi, :) = ptry (:)
      end if
      amotsa = yflu
    end function amotsa
  end subroutine amebsa
end module m_SimplexSA
