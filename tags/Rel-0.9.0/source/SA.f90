!> \brief Simulated Annealing
!> \author Alin M Elena
!> \date 12/11/07, 13:43:16

module m_SA
  use m_Constants
  use m_Useful
  use m_Gutenberg
  use m_Types
  implicit none
  private

  public :: SimulAnnealing

contains
!> \brief Simulated annealing
!> \author Alin M Elena
!> \date 12/11/07, 13:44:21
!> \remarks     Goffe, Ferrier and Rogers, J. Econometrics, 60:1/2,65-100 (Jan/Feb 1994)
!> "Global Optimization of Statistical Functions with Simulated Annealing"
!> \internal add the science citation + document parameters
  subroutine SimulAnnealing(n,x,max,rt,eps,ns,nt,neps,maxevl,lb,ub,c,&
            t,vm,xopt,fopt,nacc,nfcnev,nobds,ier,&
            fstar,xp,nacp,func,gen,atomic,sol,tb,io)
    character(len=*),parameter :: myname="SimulAnnealing"
    real(k_pr)  :: x(:), c(:), vm(:), lb(:),ub(:),fstar(:),&
            xopt(:), xp(:), t, eps, rt, fopt
    integer :: nacp(:), n, ns, nt, neps, nacc, maxevl,&
            nobds, ier, nfcnev
    logical :: max
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(modelType), intent(inout) :: tb
    type(solutionType), intent(inout) :: sol
    type(ioType), intent(inout) :: io

    interface
      function func(x,gen,atomic,tb,sol,io)
        use m_Constants
        use m_Types
        implicit none
          real(k_pr), dimension(:), intent(in) :: x
          type(generalType), intent(inout) :: gen
          type(atomicxType), intent(inout) :: atomic
          type(modelType), intent(inout) :: tb
          type(solutionType), intent(inout) :: sol
          type(ioType), intent(inout) :: io
          real(k_pr) :: func
        end function func
    end interface

      !c  type all internal variables.
    real(k_pr)  :: f, fp, p, pp, ratio
    integer  :: nup, ndown, nrej, nnew, lnobds, h, i, j, m
    logical  :: quit
    character(len=10) :: saux
!  set initial values.
      nacc = 0
      nobds = 0
      nfcnev = 0
      ier = 99
      nacp=0
      xopt(:) = x(:)
      fstar=k_infinity
        if (t <= 0.0_k_pr)then
           write(io%uout,*)"the initial temperature is not positive.  reset the variable t"
           ier = 3
           return
        end if
        do i = 1, n
            if ((x(i) < lb(i)) .or. (x(i) > ub(i))) then
              write(saux,'(i0)')i
              call error("the starting value of parameter "//trim(saux)//" is outside the bounds execution terminated "//&
              "without any optimization. lb(i) < x(i) < ub(i), i = 1, n.",myname,.false.,io)
              ier = 2
              return
            end if
        enddo

      !  evaluate the function with input x and return value as f.
        f=func(x,gen,atomic,tb,sol,io)
        if(.not. max) f = -f
        nfcnev = nfcnev + 1
        fopt = f
        fstar(1) = f
          write(io%uout,*)
          call PrintVectorA(x,"initial x",.true.,.false.,io)
          if (max) then
            write(io%uout,'(a,g25.18)')" )initial f: ", f
          else
            write(io%uout,'(a,g25.18)')" )initial f: ", -f
          end if
    do
        nup = 0
        nrej = 0
        nnew = 0
        ndown = 0
        lnobds = 0

        do m = 1, nt
            do  j = 1, ns
                do  h = 1, n
               !  generate xp, the trial value of x. note use of vm to choose xp.
                    do i = 1, n
                        if (i == h) then
                            xp(i) = x(i) + (ranmar(sol%seed)*2._k_pr- 1._k_pr) * vm(i)
                        else
                            xp(i) = x(i)
                        end if

                  !  if xp is out of bounds, select a point in bounds for the trial.
                        if((xp(i) < lb(i)) .or. (xp(i) > ub(i))) then
                            xp(i) = lb(i) + (ub(i) - lb(i))*ranmar(sol%seed)
                            lnobds = lnobds + 1
                            nobds = nobds + 1
                            call prt3(max,xp,x,f,io)
                        end if
                    enddo

              !  evaluate the function with the trial point xp and return as fp.
                    fp=func(xp,gen,atomic,tb,sol,io)
                    if(.not. max) fp = -fp
                    nfcnev = nfcnev + 1
                    call prt4(max,xp,x,fp,f,io)

               !  if too many function evaluations occur, terminate the algorithm.
                    if(nfcnev >= maxevl) then
                      write(io%uout,'(a)')"/  too many function evaluations; consider  /,  increasing maxevl or "&
                         " eps, or decreasing  /,  nt or rt. these results are likely to be  /,  poor./"
                        if (.not. max) fopt = -fopt
                        ier = 1
                        return
                    end if

               !  accept the new point if the function value increases.
                    if(fp >= f) then
                        write(io%uout,'(''  point accepted'')')
                        x(:) = xp(:)
                        f = fp
                        nacc = nacc + 1
                        nacp(h) = nacp(h) + 1
                        nup = nup + 1

                  !  if greater than any other point, record as new optimum.
                        if (fp > fopt) then
                            write(io%uout,'(''  new optimum'')')
                            xopt(:) = xp(:)
                            fopt = fp
                            nnew = nnew + 1
                        end if

                  !  if the point is lower, use the metropolis criteria to decide on
                  !  acceptance or rejection.
                    else
                        p = exprep((fp - f)/t)
                        pp = ranmar(sol%seed)
                        if (pp < p) then
                            if (max) then
                              write(io%uout,*) "though lower, point accepted"
                            else
                              write(io%uout,*) "though higher, point accepted"
                            end if
                                x(:) = xp(:)
                            f = fp
                            nacc = nacc + 1
                            nacp(h) = nacp(h) + 1
                            ndown = ndown + 1
                        else
                            nrej = nrej + 1
                            if (max) then
                              write(io%uout,'(''  lower point rejected'')')
                            else
                              write(io%uout,'(''  higher point rejected'')')
                            end if
                        end if
                    end if

                enddo
            enddo

         !  adjust vm so that approximately half of all evaluations are accepted.
            do i = 1, n
                ratio = real(nacp(i),k_pr) /real(ns,k_pr)
                if (ratio > .6_k_pr) then
                    vm(i) = vm(i)*(1._k_pr + c(i)*(ratio - .6_k_pr)/.4_k_pr)
                else if (ratio < .4_k_pr) then
                    vm(i) = vm(i)/(1._k_pr + c(i)*((.4_k_pr - ratio)/.4_k_pr))
                end if
                if (vm(i) > (ub(i)-lb(i))) then
                    vm(i) = ub(i) - lb(i)
                end if
            enddo
            call prt8(vm,xopt,x,io)
            nacp = 0
        enddo

        call prt9(max,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew,io)
!  check termination criteria.
        quit = .false.
        fstar(1) = f
        if ((fopt - fstar(1)) <= eps) quit = .true.
        do  i = 1, neps
            if (abs(f - fstar(i)) > eps) quit = .false.
        enddo

      !  terminate sa if appropriate.
        if (quit) then
            x(:) = xopt(:)
            ier = 0
            if (.not. max) fopt = -fopt
            write(io%uout,*)"  sa achieved termination criteria. ier = 0. "
            return
        end if

      ! if termination criteria is not met, prepare for another loop.
        t = rt*t
        do  i = neps, 2, -1
            fstar(i) = fstar(i-1)
        enddo
        f = fopt
        x(:) = xopt(:)

        if (t<epsilon(t)) then
            ier=111
            write(io%uout,*)"temperature has reached machine precision !!!"
            exit
        endif
      !  loop again.
    enddo

    end subroutine SimulAnnealing

  subroutine prt3(max,xp,x,f,io)
    real(k_pr), intent(in) :: xp(:), x(:), f
    logical, intent(in) ::  max
    type(ioType), intent(inout) :: io
        write(io%uout,'(''  '')')
        call PrintVectorA(x,'current x',.true.,.false.,io)
        if (max) then
            write(io%uout,'(a,g25.18)')"  current f: ", f
        else
            write(io%uout,'(a,g25.18)') -f
        end if
          call PrintVectorA(xp,'trial x',.true.,.false.,io)
        write(io%uout,*)"  point rejected since out of bounds "
    end subroutine prt3

    subroutine prt4(max,xp,x,fp,f,io)
      real(k_pr),intent(in) ::  xp(:), x(:), fp, f
      logical,intent(in) ::  max
      type(ioType), intent(inout) :: io

      write(io%uout,*)
      call PrintVectorA(x,'current x',.true.,.false.,io)
      if (max) then
        write(io%uout,'(a,g25.18)')" current f: ", f
        call PrintVectorA(xp,'trial x',.true.,.false.,io)
        write(io%uout,'(a,g25.18)')" resulting f: ", fp
      else
        write(io%uout,'(a,g25.18)')" current f: ", -f
        call PrintVectorA(xp,'trial x',.true.,.false.,io)
        write(io%uout,'(a,g25.18)')" resulting f: ", -fp
      end if
    end subroutine prt4

    subroutine prt8(vm,xopt,x,io)
        type(ioType), intent(inout) :: io
        real(k_pr) :: vm(:), xopt(:), x(:)


        write(io%uout,*)" intermediate results after step length adjustment"
          call PrintVectorA(vm,'new step length (vm)',.true.,.false.,io)
          call PrintVectorA(xopt,'current optimal x',.true.,.false.,io)
          call PrintVectorA(x,'current x',.true.,.false.,io)
        write(io%uout,*)
    end subroutine prt8

    subroutine prt9(max,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew,io)
        type(ioType), intent(inout) :: io
        real(k_pr) ::          xopt(:), vm(:), t, fopt
        integer::  nup, ndown, nrej, lnobds, nnew, totmov
        logical ::  max

        totmov = nup + ndown + nrej

        write(io%uout,*)" intermediate results before next temperature reduction"
        write(io%uout,'(a,g12.5)')" current temperature:            ", t
        if (max) then
            write(io%uout,'(a,g12.5)')" max function value so far:  ", fopt
            write(io%uout,'(a,g12.5)')"  total moves:                ", totmov
            write(io%uout,'(a,i8)')"     uphill:                  ", nup
            write(io%uout,'(a,i8)') "    accepted downhill:       ",ndown
            write(io%uout, '(a,i8)')  "  rejected downhill:       ", nrej
            write(io%uout,'(a,i8)')  "out of bounds trials:       ",lnobds
            write(io%uout,'(a,i8)')"  new maxima this temperature:", nnew
        else
            write(io%uout,'(a,g12.5)')" min function value so far:  ", -fopt
            write(io%uout,'(a,g12.5)')"  total moves:                ", totmov
            write(io%uout,'(a,i8)')"     uphill:                  ", nup
            write(io%uout,'(a,i8)') "    accepted downhill:       ",ndown
            write(io%uout, '(a,i8)')  "  rejected downhill:       ", nrej
            write(io%uout,'(a,i8)')  "out of bounds trials:       ",lnobds
            write(io%uout,'(a,i8)')"  new minima this temperature:", nnew
        end if
          call PrintVectorA(xopt,'current optimal x',.true.,.false.,io)
          call PrintVectorA(vm,'step length (vm)',.true.,.false.,io)
        write(io%uout,'('' '')')
  end subroutine prt9

end module m_SA