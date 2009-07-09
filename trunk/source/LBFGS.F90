!> \brief implements limited memory bfgs method for large scale optimization
!> \details (citation needed)
!> \author Alin M. Elena
!> \date 18th of January 2008
module m_LBFGS
  use m_Constants
  use m_Useful
  use m_Types
  use m_Gutenberg
#ifdef MKL95
  use mkl95_BLAS, only: dot, axpy
#endif
  implicit none
!
  private
!
  public :: linearBFGS
!
contains
!
!> \brief driver routine for linearBFGS optimization
!> \author Alin M Elena
!> \date 10/11/07, 15:22:53
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tb type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine linearBFGS (io, gen, atomic, tb, sol, func)
    character (len=*), parameter :: myname = 'linearBFGS'
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
    integer :: m, n, res
    real (k_pr) :: epsx, epsf, epsg, gtol, xtol, ftol
    real (k_pr), allocatable :: x (:)
    integer :: info, atom, iprint (1:2), maxfev
!
    write (io%uout, '(/a/)') "--LBFGS Optimization Run-----------------------------------------"
    n = atomic%atoms%nmoving * 3
    m = Min (gen%HessianM, n)
    if (io%verbosity > k_highVerbos) then
      iprint (1) = 1
      iprint (2) = 3
    else
      iprint (1) = - 1
      iprint (2) = 0
    end if
!
    epsf = gen%epsf
    epsg = gen%epsg
    epsx = gen%epsx
    xtol = gen%xtol
    gtol = gen%gtol
    info = 0
    maxfev = gen%maxFeval
    ftol = gen%ftol
    allocate (x(1:n))
!
    do i = 1, atomic%atoms%nmoving
      atom = atomic%atoms%id(atomic%atoms%moving(i))
      x (3*(i-1)+1) = atomic%atoms%x(atom)
      x (3*(i-1)+2) = atomic%atoms%y(atom)
      x (3*(i-1)+3) = atomic%atoms%z(atom)
    end do
    call lbfgs (n, m, x, epsg, epsf, epsx, xtol, gtol, ftol, maxfev, gen%nsteps, iprint, info, func, gen, atomic, tb, sol, io)
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
    if (info < 0) then
      write (io%uout, '(/a,i4)') "WARNING: Optimization did not converge.", info
    else
      write (io%uout, '(/a,i4)') "LBFGS report: Optimization converged with code: ", info
    end if
    write (io%uout, '(/a)') "----------------------------------------------------------------"
  end subroutine linearBFGS
!
!
!> \brief        limited memory bfgs method for large scale optimization
!> \details     this subroutine solves the unconstrained minimization problem
!> \details   min f(x),    x= (x1,x2,...,xn),
!> \details  using the limited memory bfgs method
!> \author jorge nocedal (see module description), Alin M Elena
!> \date   july 1990,  18th of January 2008
!>\param   n       is an integer variable that must be set by the user to the
!>             number of variables. it is not altered by the routine.
!>             restriction: n>0.
!>
!> \param     m       is an integer variable that must be set by the user to
!>             the number of corrections used in the bfgs update. it
!>             is not altered by the routine. values of m less than 3 are
!>             not recommended; large values of m will result in excessive
!>             computing time. 3<= m <=7 is recommended. restriction: m>0.
!>
!> \param    x       is a double precision array of length n. on initial entry
!>             it must be set by the user to the values of the initial
!>             estimate of the solution vector. on exit with iflag=0, it
!>             contains the values of the variables at the best point
!>             found (usually a solution).
!>
!> \param     f       is a double precision variable. before initial entry and on
!>             a re-entry with iflag=1, it must be set by the user to
!>             contain the value of the function f at the point x.
!>
!> \param    g       is a double precision array of length n. before initial
!>             entry and on a re-entry with iflag=1, it must be set by
!>             the user to contain the components of the gradient g at
!>             the point x.
!>
!> \param     diagco  is a logical variable that must be set to .true. if the
!>             user  wishes to provide the diagonal matrix hk0 at each
!>            iteration. otherwise it should be set to .false., in which
!>            case  lbfgs will use a default value described below. if
!>            diagco is set to .true. the routine will return at each
!>            iteration of the algorithm with iflag=2, and the diagonal
!>             matrix hk0  must be provided in the array diag.
!>
!>
!> \param   diag    is a double precision array of length n. if diagco=.true.,
!>            then on initial entry or on re-entry with iflag=2, diag
!>            it must be set by the user to contain the values of the
!>            diagonal matrix hk0.  restriction: all elements of diag
!>            must be positive.
!>
!> \param    iprint  is an integer array of length two which must be set by the
!>            user.
!>
!>            iprint(1) specifies the frequency of the output:
!>               iprint(1) < 0 : no output is generated,
!>               iprint(1) = 0 : output only at first and last iteration,
!>               iprint(1) > 0 : output every iprint(1) iterations.
!>
!>             iprint(2) specifies the type of output generated:
!>                iprint(2) = 0 : iteration count, number of function
!>                                evaluations, function value, norm of the
!>                                gradient, and steplength,
!>                iprint(2) = 1 : same as iprint(2)=0, plus vector of
!>                                variables and  gradient vector at the
!>                               initial point,
!>               iprint(2) = 2 : same as iprint(2)=1, plus vector of
!>                               variables,
!>               iprint(2) = 3 : same as iprint(2)=2, plus gradient vector.
!>
!>
!> \param   eps     is a positive double precision variable that must be set by
!>            the user, and determines the accuracy with which the solution
!>            is to be found. the subroutine terminates when
!>
!>                        ||g|| < eps max(1,||x||),
!>
!>            where ||.|| denotes the euclidean norm.
!>
!> \param   xtol    is a  positive double precision variable that must be set by
!>            the user to an estimate of the machine precision (e.g.
!>            10**(-16) on a sun station 3/60). the line search routine will
!>            terminate if the relative width of the interval of uncertainty
!>            is less than xtol.
!>
!> \param   w       is a double precision array of length n(2m+1)+2m used as
!>            workspace for lbfgs. this array must not be altered by the
!>            user.
!>
!> \param   iflag   is an integer variable that must be set to 0 on initial entry
!>            to the subroutine. a return with iflag<0 indicates an error,
!>            and iflag=0 indicates that the routine has terminated without
!>            detecting errors. on a return with iflag=1, the user must
!>            evaluate the function f and gradient g. on a return with
!>            iflag=2, the user must provide the diagonal matrix hk0.
!>
!>            the following negative values of iflag, detecting an error,
!>            are possible:
!>
!>             iflag=-1  the line search routine mcsrch failed. the
!>                       parameter info provides more detailed information
!>                       (see also the documentation of mcsrch):
!>
!>                      info = 0  improper input parameters.
!>
!>                      info = 2  relative width of the interval of
!>                                uncertainty is at most xtol.
!>
!>                      info = 3  more than n function evaluations were
!>                                required at the present iteration.
!>
!>                      info = 4  the step is too small.
!>
!>                      info = 5  the step is too large.
!>
!>                      info = 6  rounding errors prevent further progress.
!>                                there may not be a step which satisfies
!>                                the sufficient decrease and curvature
!>                                conditions. tolerances may be too small.
!>
!>
!>             iflag=-2  the i-th diagonal element of the diagonal inverse
!>                       hessian approximation, given in diag, is not
!>                       positive.
!>
!>              iflag=-3  improper input parameters for lbfgs (n or m are
!>                        not positive).
!> \remarks the routine is especially
!>     effective on problems involving a large number of variables. in
!>     a typical iteration of this method an approximation hk to the
!>     inverse of the hessian is obtained by applying m bfgs updates to
!>     a diagonal matrix hk0, using information from the previous m steps.
!>     the user specifies the number m, which determines the amount of
!>     storage required by the routine. the user may also provide the
!>     diagonal matrices hk0 if not satisfied with the default choice.
!>    the algorithm is described in "on the limited memory bfgs method
!>    for large scale optimization", by d. liu and j. nocedal,
!>      mathematical programming b 45 (1989) 503-528.
!>
!>      the user is required to calculate the function value f and its
!>      gradient g. in order to allow the user complete control over
!>      these computations, reverse  communication is used. the routine
!>    must be called repeatedly under the control of the parameter
!>     iflag.
!>
!>     the steplength is determined at each iteration by means of the
!>     line search routine mcvsrch, which is a slight modification of
!>     the routine csrch written by more' and thuente.
  subroutine lbfgs (n, m, x, epsg, epsf, epsx, xtol, gtol, ftol, maxfev, maxits, iprint, info, func, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "lbfgs"
    real (k_pr), intent (in) :: xtol, epsg, epsf, epsx, gtol, ftol
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    integer, intent (in) :: maxits, n, m, maxfev, iprint (1:2)
    integer, intent (inout) :: info
    real (k_pr), intent (inout) :: x (:)
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
    real (k_pr) :: gnorm, stp1, stp, ys, yy, sq, yr, beta, xnorm
    integer :: iter, nfun, point, ispt, iypt, bound, npt, cp, i, nfev, inmc, iycn, iscn
    logical :: finish
!   character(len=255) :: saux
    real (k_pr), allocatable :: g (:), diag (:), w (:), tx (:), ta (:), xold (:)
    real (k_pr) :: f, fold, tf, txnorm, v
    integer :: res
    real (k_pr) :: stpmin, stpmax
#ifndef MKL95
    real(k_pr),external :: ddot
#endif

!
    allocate (g(1:n))
    allocate (diag(1:n))
    allocate (w(1:n*(2*m+1)+2*m))
    allocate (xold(1:n))
    allocate (tx(1:n))
    allocate (ta(1:n))
    g = 0.0_k_pr
    diag = 0.0_k_pr
    w = 0.0_k_pr
    tx = 0.0_k_pr
    ta = 0.0_k_pr
    f = 0.0_k_pr
    write (io%uout, '(a,a)') "Energy update from: ", myname
    res = func (gen, atomic, tb, sol, io, x, f, g)
!
    fold = f
    iter = 0
    info = 0
    if (n <= 0 .or. m <= 0 .or. m > n .or. epsg < 0.0_k_pr .or. epsf < 0.0_k_pr .or. epsx < 0.0_k_pr .or. maxits < 0) then
      info = - 1
      return
    end if
    nfun = 1
    point = 0
    finish = .false.
    diag = 1.0_k_pr
    stpmin = 10.0e-20_k_pr
    stpmax = 10.0e20_k_pr
    ispt = n + 2 * m
    iypt = ispt + n * m
    do i = 1, n
      w (ispt+i) = - g (i) * diag (i)
    end do
#ifdef MKL95
    gnorm = Sqrt (dot(g, g))
#else
   gnorm = Sqrt (ddot(n,g,k_ione, g,k_ione))
#endif
    stp1 = 1 / gnorm
!
    do while (.true.)
      xold = x
      iter = iter + 1
      info = 0
      bound = iter - 1
      if (iter /= 1) then
        if (iter > m) then
          bound = m
        end if
#ifdef MKL95
        ys = dot (w(iypt+npt+1:iypt+npt+n), w(ispt+npt+1:ispt+npt+n))
        yy = dot (w(iypt+npt+1:iypt+npt+n), w(iypt+npt+1:iypt+npt+n))
#else
        ys = ddot (n,w(iypt+npt+1:iypt+npt+n),k_ione, w(ispt+npt+1:ispt+npt+n),k_ione)
        yy = ddot (n,w(iypt+npt+1:iypt+npt+n),k_ione, w(iypt+npt+1:iypt+npt+n),k_ione)
#endif
        diag = ys / yy
        cp = point
        if (point == 0) then
          cp = m
        end if
        w (n+cp) = 1 / ys
        w (1:n) = - g (1:n)
        cp = point
        do i = 1, bound
          cp = cp - 1
          if (cp ==-1) then
            cp = m - 1
          end if
#ifdef MKL95
          sq = dot (w(ispt+cp*n+1:ispt+cp*n+n), w(1:n))
#else
          sq = ddot (n,w(ispt+cp*n+1:ispt+cp*n+n),k_ione, w(1:n),k_ione)
#endif
          inmc = n + m + cp + 1
          iycn = iypt + cp * n
          w (inmc) = w (n+cp+1) * sq
#ifdef MKL95
          call axpy (w(iycn+1:iycn+n), w(1:n),-w(inmc))
#else
          call daxpy (n,w(iycn+1:iycn+n), w(1:n),k_ione,-w(inmc),k_ione)
#endif
        end do
        do i = 1, n
          w (i) = diag (i) * w (i)
        end do
        do i = 1, bound
#ifdef MKL95
          yr = dot (w(iypt+cp*n+1:iypt+cp*n+n), w(1:n))
#else
          yr = ddot (n,w(iypt+cp*n+1:iypt+cp*n+n),k_ione, w(1:n),k_ione)
#endif
          beta = w (n+cp+1) * yr
          inmc = n + m + cp + 1
          beta = w (inmc) - beta
          iscn = ispt + cp * n
#ifdef MKL95
          call axpy (w(iscn+1:iscn+n), w(1:n), beta)
#else
          call daxpy (n,w(iscn+1:iscn+n), w(1:n),k_ione, beta,k_ione)
#endif
          cp = cp + 1
          if (cp == m) then
            cp = 0
          end if
        end do
        do i = 1, n
          w (ispt+point*n+i) = w (i)
        end do
      end if
      nfev = 0
      stp = 1
      if (iter == 1) then
        stp = stp1
      end if
      w (1:n) = g (1:n)
      call mcsrch (n, x, f, g, w(ispt+point*n+1:ispt+point*n+n), stp, ftol, xtol, maxfev, info, nfev, diag, gtol, stpmin, stpmax, &
     & func, gen, atomic, tb, sol, io)
      if (info /= 1) then
        if (info == 0) then
          info = - 1
          return
        end if
      end if
      nfun = nfun + nfev
      npt = point * n
!!!! potentially thread unsafe check that ispt/=iypt
      do i = 1, n
        w (ispt+npt+i) = stp * w (ispt+npt+i)
        w (iypt+npt+i) = g (i) - w (i)
      end do
      point = point + 1
      if (point == m) then
        point = 0
      end if
      if ((iter > maxits) .and. (maxits > 0)) then
        info = 5
        return
      end if
      call IterationReport (iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish, io)
#ifdef MKL95
      gnorm = Sqrt (dot(g, g))
#else
      gnorm = Sqrt (ddot(n,g,k_ione, g,k_ione))
#endif
      if (gnorm <= epsg) then
        info = 4
        return
      end if
      tf = Max (Abs(fold), Max(Abs(f), 1.0_k_pr))
      if (Abs(fold-f) <= epsf*tf) then
        info = 1
        return
      end if
!
      tx = xold
      tx = tx - x
#ifdef MKL95
      xnorm = Sqrt (dot(x, x))
      txnorm = Max (xnorm, Sqrt(dot(xold, xold)))
      v = Sqrt (dot(tx, tx))
#else
      xnorm = Sqrt (ddot(n,x,k_ione, x,k_ione))
      txnorm = Max (xnorm, Sqrt(ddot(n,xold,k_ione, xold,k_ione)))
      v = Sqrt (ddot(n,tx,k_ione, tx,k_ione))
#endif
      txnorm = Max (txnorm, 1.0_k_pr)
      if (v <= epsx) then
        info = 2
        return
      end if
      fold = f
      xold = x
    end do
!
    deallocate (g, diag, w, tx, ta, xold)
!
  end subroutine lbfgs
!
!
  subroutine IterationReport (iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish, io)
!
!     -------------------------------------------------------------
!     this routine prints monitoring information. the frequency and
!     amount of output are controlled by iprint.
!     -------------------------------------------------------------
!
    integer, intent (in) :: iprint (2), iter, nfun, n, m
    real (k_pr), intent (in) :: x (:), g (:), f, gnorm, stp
    type (ioType), intent (in) :: io
    logical, intent (in) :: finish
    integer :: i
!
    if (iter == 0) then
      write (io%uout, '(a)') "=============================================================="
      write (io%uout, '(a,i0,a,i0,a)') "n= ", n, " number of corrections= ", m, "initial values"
      write (io%uout, '(a,ES12.4,a,ES12.4)') "f=", f, " gnorm= ", gnorm
      if (iprint(2) >= 1) then
        write (io%uout, '(a)') "vector X:"
        write (io%uout, '(1500f16.8)') (x(i), i=1, n)
        write (io%uout, '(a)') "gradient vector (Forces): "
        write (io%uout, '(1500f16.8)') (g(i), i=1, n)
      end if
      write (io%uout, '(a)') "=============================================================="
      write (io%uout, '(a,4x,a,8x,a,7x,a)') "   i   nfn", "func", "gnorm", "steplength"
    else
      if ((iprint(1) == 0) .and. (iter /= 1 .and. .not. finish)) return
      if (iprint(1) /= 0) then
        if (Mod(iter-1, iprint(1)) == 0 .or. finish) then
          if (iprint(2) > 1 .and. iter > 1) write (io%uout, '(a,4x,a,8x,a,7x,a)') "   i   nfn", "func", "gnorm", "steplength"
          write (io%uout, '(2(i4,1x),3x,3(1f10.3,2x))') iter, nfun, f, gnorm, stp
        else
          return
        end if
      else
        if (iprint(2) > 1 .and. finish) write (io%uout, '(a,4x,a,8x,a,7x,a)') "   i   nfn", "func", "gnorm", "steplength"
        write (io%uout, '(2(i4,1x),3x,3(1f10.3,2x))') iter, nfun, f, gnorm, stp
      end if
      if (iprint(2) == 2 .or. iprint(2) == 3) then
        if (finish) then
          write (io%uout, '(a)') "final point X: "
        else
          write (io%uout, '(a)') "vector X:"
        end if
        write (io%uout, '(1500f16.8)') (x(i), i=1, n)
        if (iprint(2) == 3) then
          write (io%uout, '(a)') "gradient vector (Forces): "
          write (io%uout, '(1500f16.8)') (g(i), i=1, n)
        end if
      end if
      if (finish) write (io%uout, '(a)') " the minimization terminated without detecting errors. iflag = 0"
    end if
  end subroutine IterationReport
!
! !c
! !c     **************************
! !c     line search routine mcsrch
! !c     **************************
! !c
! !c
! !c                     subroutine mcsrch
! !c
! !c     a slight modification of the subroutine csrch of more' and thuente.
! !c     the changes are to allow reverse communication, and do not affect
! !c     the performance of the routine.
! !c
! !c     the purpose of mcsrch is to find a step which satisfies
! !c     a sufficient decrease condition and a curvature condition.
! !c
! !c     at each stage the subroutine updates an interval of
! !c     uncertainty with endpoints stx and sty. the interval of
! !c     uncertainty is initially chosen so that it contains a
! !c     minimizer of the modified function
! !c
! !c          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
! !c
! !c     if a step is obtained for which the modified function
! !c     has a nonpositive function value and nonnegative derivative,
! !c     then the interval of uncertainty is chosen so that it
! !c     contains a minimizer of f(x+stp*s).
! !c
! !c     the algorithm is designed to find a step which satisfies
! !c     the sufficient decrease condition
! !c
! !c           f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),
! !c
! !c     and the curvature condition
! !c
! !c           abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).
! !c
! !c     if ftol is less than gtol and if, for example, the function
! !c     is bounded below, then there is always a step which satisfies
! !c     both conditions. if no step can be found which satisfies both
! !c     conditions, then the algorithm usually stops when rounding
! !c     errors prevent further progress. in this case stp only
! !c     satisfies the sufficient decrease condition.
! !c
! !c     the subroutine statement is
! !c
! !c        subroutine mcsrch(n,x,f,g,s,stp,ftol,xtol, maxfev,info,nfev,wa)
! !c     where
! !c
! !c       n is a positive integer input variable set to the number
! !c         of variables.
! !c
! !c       x is an array of length n. on input it must contain the
! !c         base point for the line search. on output it contains
! !c         x + stp*s.
! !c
! !c       f is a variable. on input it must contain the value of f
! !c         at x. on output it contains the value of f at x + stp*s.
! !c
! !c       g is an array of length n. on input it must contain the
! !c         gradient of f at x. on output it contains the gradient
! !c         of f at x + stp*s.
! !c
! !c       s is an input array of length n which specifies the
! !c         search direction.
! !c
! !c       stp is a nonnegative variable. on input stp contains an
! !c         initial estimate of a satisfactory step. on output
! !c         stp contains the final estimate.
! !c
! !c       ftol and gtol are nonnegative input variables.
! !          termination occurs when the sufficient decrease
! !c         condition and the directional derivative condition are
! !c         satisfied.
! !c
! !c       xtol is a nonnegative input variable. termination occurs
! !c         when the relative width of the interval of uncertainty
! !c         is at most xtol.
! !c
! !c       stpmin and stpmax are nonnegative input variables which
! !c         specify lower and upper bounds for the step. (in this reverse
! !c         communication implementatin they are defined in a common
! !c         statement).
! !c
! !c       maxfev is a positive integer input variable. termination
! !c         occurs when the number of calls to fcn is at least
! !c         maxfev by the end of an iteration.
! !c
! !c       info is an integer output variable set as follows:
! !c
! !c         info = 0  improper input parameters.
! !c
! !c         info =-1  a return is made to compute the function and gradient.
! !c
! !c         info = 1  the sufficient decrease condition and the
! !c                   directional derivative condition hold.
! !c
! !c         info = 2  relative width of the interval of uncertainty
! !c                   is at most xtol.
! !c
! !c         info = 3  number of calls to fcn has reached maxfev.
! !c
! !c         info = 4  the step is at the lower bound stpmin.
! !c
! !c         info = 5  the step is at the upper bound stpmax.
! !c
! !c         info = 6  rounding errors prevent further progress.
! !c                   there may not be a step which satisfies the
! !c                   sufficient decrease and curvature conditions.
! !c                   tolerances may be too small.
! !c
! !c       nfev is an integer output variable set to the number of
! !c         calls to fcn.
! !c
! !c       wa is a work array of length n.
! !c     **********
!
  subroutine mcsrch (n, x, f, g, s, stp, ftol, xtol, maxfev, info, nfev, wa, gtol, stpmin, stpmax, func, gen, atomic, tb, sol, io)
    character (len=*), parameter :: myname = "mcsrch"
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tb
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    integer, intent (in) :: n, maxfev
    integer, intent (inout) :: info, nfev
    real (k_pr), intent (inout) :: f, stp
    real (k_pr), intent (in) :: ftol, gtol, xtol, stpmin, stpmax
    real (k_pr), intent (inout) :: x (:), g (:), wa (:), s (:)
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
!
    integer :: res
    integer :: infoc, j
    logical :: brackt, stage1
    real (k_pr) :: dg, dgm, dginit, dgtest, dgx, dgxm, dgy, dgym
    real (k_pr) :: finit, ftest1, fm, fx, fxm, fy, fym, p5, p66, stx, sty
    real (k_pr) :: stmin, stmax, width, width1, xtrapf, mytemp
!
    p5 = 0.5_k_pr
    p66 = 0.66_k_pr
    xtrapf = 4.0_k_pr
    f = 0.0_k_pr
    write (io%uout, '(a,a)') "Energy update from: ", myname
    res = func (gen, atomic, tb, sol, io, x, f, g)
    infoc = 1
    info = 0
    if (n <= 0 .or. stp <= k_zero .or. ftol < 0.0_k_pr .or. gtol < k_zero .or. xtol < k_zero .or. stpmin < k_zero .or. stpmax < &
   & stpmin .or. maxfev <= 0) then
      return
    end if
    dginit = 0.0_k_pr
    do j = 1, n
      dginit = dginit + g (j) * s (j)
    end do
    if (dginit >= 0) then
      return
    end if
    brackt = .false.
    stage1 = .true.
    nfev = 0
    finit = f
    dgtest = ftol * dginit
    width = stpmax - stpmin
    width1 = width / p5
    do j = 1, n
      wa (j) = x (j)
    end do
    stx = 0.0_k_pr
    fx = finit
    dgx = dginit
    sty = 0.0_k_pr
    fy = finit
    dgy = dginit
    do while (.true.)
      if (brackt) then
        if (stx < sty) then
          stmin = stx
          stmax = sty
        else
          stmin = sty
          stmax = stx
        end if
      else
        stmin = stx
        stmax = stp + xtrapf * (stp-stx)
      end if
      if (stp > stpmax) then
        stp = stpmax
      end if
      if (stp < stpmin) then
        stp = stpmin
      end if
      if (brackt .and. (stp <= stmin .or. stp >= stmax) .or. nfev >= maxfev-1 .or. infoc == 0 .or. brackt .and. stmax-stmin <= &
     & xtol*stmax) then
        stp = stx
      end if
      do j = 1, n
        x (j) = wa (j) + stp * s (j)
      end do
      write (io%uout, '(a,a)') "Energy update from: ", myname
      res = func (gen, atomic, tb, sol, io, x, f, g)
      info = 0
      nfev = nfev + 1
      dg = 0.0_k_pr
      do j = 1, n
        dg = dg + g (j) * s (j)
      end do
      ftest1 = finit + stp * dgtest
      if (brackt .and. (stp <= stmin .or. stp >= stmax) .or. (infoc == 0)) then
        info = 6
      end if
      if ((stp == stpmax) .and. (f <= ftest1) .and. (dg <= dgtest)) then
        info = 5
      end if
      if (stp == stpmin .and. (f > ftest1 .or. dg >= dgtest)) then
        info = 4
      end if
      if (nfev >= maxfev) then
        info = 3
      end if
      if (brackt .and. stmax-stmin <= xtol*stmax) then
        info = 2
      end if
      if (f <= ftest1 .and. Abs(dg) <=-gtol*dginit) then
        info = 1
      end if
      if (info /= 0) then
        return
      end if
      mytemp = ftol
      if (gtol < ftol) then
        mytemp = gtol
      end if
      if (stage1 .and. f <= ftest1 .and. dg >= mytemp*dginit) then
        stage1 = .false.
      end if
      if (stage1 .and. f <= fx .and. f > ftest1) then
        fm = f - stp * dgtest
        fxm = fx - stx * dgtest
        fym = fy - sty * dgtest
        dgm = dg - dgtest
        dgxm = dgx - dgtest
        dgym = dgy - dgtest
        call mcstep (stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc)
        fx = fxm + stx * dgtest
        fy = fym + sty * dgtest
        dgx = dgxm + dgtest
        dgy = dgym + dgtest
      else
        call mcstep (stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc)
      end if
      if (brackt) then
        if (Abs(sty-stx) >= p66*width1) then
          stp = stx + p5 * (sty-stx)
        end if
        width1 = width
        width = Abs (sty-stx)
      end if
    end do
  end subroutine mcsrch
!
! !c
! !c     subroutine mcstep
! !c
! !c     the purpose of mcstep is to compute a safeguarded step for
! !c     a linesearch and to update an interval of uncertainty for
! !c     a minimizer of the function.
! !c
! !c     the parameter stx contains the step with the least function
! !c     value. the parameter stp contains the current step. it is
! !c     assumed that the derivative at stx is negative in the
! !c     direction of the step. if brackt is set true then a
! !c     minimizer has been bracketed in an interval of uncertainty
! !c     with endpoints stx and sty.
! !c
! !c     the subroutine statement is
! !c
! !c       subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
! !c                        stpmin,stpmax,info)
! !c
! !c     where
! !c
! !c       stx, fx, and dx are variables which specify the step,
! !c         the function, and the derivative at the best step obtained
! !c         so far. the derivative must be negative in the direction
! !c         of the step, that is, dx and stp-stx must have opposite
! !c         signs. on output these parameters are updated appropriately.
! !c
! !c       sty, fy, and dy are variables which specify the step,
! !c         the function, and the derivative at the other endpoint of
! !c         the interval of uncertainty. on output these parameters are
! !c         updated appropriately.
! !c
! !c       stp, fp, and dp are variables which specify the step,
! !c         the function, and the derivative at the current step.
! !c         if brackt is set true then on input stp must be
! !c         between stx and sty. on output stp is set to the new step.
! !c
! !c       brackt is a logical variable which specifies if a minimizer
! !c         has been bracketed. if the minimizer has not been bracketed
! !c         then on input brackt must be set false. if the minimizer
! !c         is bracketed then on output brackt is set true.
! !c
! !c       stpmin and stpmax are input variables which specify lower
! !c         and upper bounds for the step.
! !c
! !c       info is an integer output variable set as follows:
! !c         if info = 1,2,3,4,5, then the step has been computed
! !c         according to k_one of the five cases below. otherwise
! !c         info = 0, and this indicates improper input parameters.
! !c
! !c     subprograms called
! !c
! !c       fortran-supplied ... abs,max,min,sqrt
! !c
! !c     argonne national laboratory. minpack project. june 1983
! !c     jorge j. more', david j. thuente
! !c
!       real(k_pr) :: gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta
!       info = 0
! !c
! !c     check the input parameters for errors.
! !c
!       if ((brackt .and. (stp <= min(stx,sty) .or. &
!          stp >= max(stx,sty))) .or. &
!          dx*(stp-stx) >= 0.0_k_pr .or. stpmax < stpmin) return
! !c
! !c     determine if the derivatives have opposite sign.
! !c
!       sgnd = dp*(dx/abs(dx))
! !c
! !c     first case. a higher function value.
! !c     the minimum is bracketed. if the cubic step is closer
! !c     to stx than the quadratic step, the cubic step is taken,
! !c     else the average of the cubic and quadratic steps is taken.
! !c
!
  subroutine mcstep (stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, stmin, stmax, info)
    integer, intent (inout) :: info
    real (k_pr), intent (inout) :: stx, fx, dx, sty, fy, dy, stp, fp, dp
    real (k_pr), intent (in) :: stmin, stmax
    logical, intent (inout) :: brackt
    real (k_pr) :: p, q, r, s, sgnd, stpc, stpf, stpq, theta, p66, gamma
    logical :: bound
!
    p66 = 0.66_k_pr
    info = 0
    if (brackt .and. (stp <= Min(stx, sty) .or. stp >= Max(stx, sty)) .or. dx*(stp-stx) >= k_zero .or. stmax < stmin) then
      return
    end if
    sgnd = dp * (dx/Abs(dx))
    if (fp > fx) then
      info = 1
      bound = .true.
      theta = 3.0_k_pr * (fx-fp) / (stp-stx) + dx + dp
      s = Max (Abs(theta), Max(Abs(dx), Abs(dp)))
      gamma = s * Sqrt ((theta*theta/s/s)-dx/s*(dp/s))
      if (stp < stx) then
        gamma = - gamma
      end if
      p = gamma - dx + theta
      q = gamma - dx + gamma + dp
      r = p / q
      stpc = stx + r * (stp-stx)
      stpq = stx + dx / ((fx-fp)/(stp-stx)+dx) / 2.0_k_pr * (stp-stx)
      if (Abs(stpc-stx) < Abs(stpq-stx)) then
        stpf = stpc
      else
        stpf = stpc + (stpq-stpc) / 2
      end if
      brackt = .true.
    else
      if (sgnd < 0) then
        info = 2
        bound = .false.
        theta = 3.0_k_pr * (fx-fp) / (stp-stx) + dx + dp
        s = Max (Abs(theta), Max(Abs(dx), Abs(dp)))
        gamma = s * Sqrt (theta*theta/s/s-dx/s*(dp/s))
        if (stp > stx) then
          gamma = - gamma
        end if
        p = gamma - dp + theta
        q = gamma - dp + gamma + dx
        r = p / q
        stpc = stp + r * (stx-stp)
        stpq = stp + dp / (dp-dx) * (stx-stp)
        if (Abs(stpc-stp) > Abs(stpq-stp)) then
          stpf = stpc
        else
          stpf = stpq
        end if
        brackt = .true.
!
      else
        if (Abs(dp) < Abs(dx)) then
          info = 3
          bound = .true.
          theta = 3.0_k_pr * (fx-fp) / (stp-stx) + dx + dp
          s = Max (Abs(theta), Max(Abs(dx), Abs(dp)))
          gamma = s * Sqrt (Max(k_zero, theta*theta/s/s-dx/s*(dp/s)))
          if (stp > stx) then
            gamma = - gamma
          end if
          p = gamma - dp + theta
          q = gamma + (dx-dp) + gamma
          r = p / q
          if (r < 0 .and. gamma /= 0) then
            stpc = stp + r * (stx-stp)
          else
            if (stp > stx) then
              stpc = stmax
            else
              stpc = stmin
            end if
          end if
          stpq = stp + dp / (dp-dx) * (stx-stp)
          if (brackt) then
            if (Abs(stp-stpc) < Abs(stp-stpq)) then
              stpf = stpc
            else
              stpf = stpq
            end if
          else
            if (Abs(stp-stpc) > Abs(stp-stpq)) then
              stpf = stpc
            else
              stpf = stpq
            end if
          end if
        else
          info = 4
          bound = .false.
          if (brackt) then
            theta = 3.0_k_pr * (fp-fy) / (sty-stp) + dy + dp
            s = Max (Abs(theta), Max(Abs(dy), Abs(dp)))
            gamma = s * Sqrt (theta*theta/s/s-dy/s*(dp/s))
            if (stp > sty) then
              gamma = - gamma
            end if
            p = gamma - dp + theta
            q = gamma - dp + gamma + dy
            r = p / q
            stpc = stp + r * (sty-stp)
            stpf = stpc
          else
            if (stp > stx) then
              stpf = stmax
            else
              stpf = stmin
            end if
          end if
        end if
      end if
    end if
    if (fp > fx) then
      sty = stp
      fy = fp
      dy = dp
    else
      if (sgnd < k_zero) then
        sty = stx
        fy = fx
        dy = dx
      end if
      stx = stp
      fx = fp
      dx = dp
    end if
    stpf = Min (stmax, stpf)
    stpf = Max (stmin, stpf)
    stp = stpf
    if (brackt .and. bound) then
      if (sty > stx) then
        stp = Min (stx+p66*(sty-stx), stp)
      else
        stp = Max (stx+p66*(sty-stx), stp)
      end if
    end if
  end subroutine mcstep
!
end module m_LBFGS
