!> \brief deals with the general tight binding functions
!> \author Alin M Elena
!> \date 01/11/07, 11:48:23
!
module m_TightBinding
  use m_Constants
  use m_Types
  use m_Useful
  use m_TailFunctions
  use m_Gutenberg
  use m_LinearAlgebra
  implicit none
!
  private
!
  public :: SetSolutionSpace
  public :: rad
  public :: RadP
  public :: RadPp
  public :: rep
  public :: RepP
  public :: RepPp
  public :: onsite
  public :: embedding
  public :: EmbeddingP
  public :: EmbeddingPp
  public :: InitMagneticMoment
  public :: ComputeMagneticMoment
  public :: SetTails
  public :: ComputeBondCurrents
  public :: ComputeBondCurrentsOnOrbitals
  public :: UpdateNeighboursList
!
contains
!> \brief Initializes the solution space
!> \details allocates memory for all the variables necessary for the calculation. It will
!> also call the routines that precompute some values (eg factorial, Gaunt coefficients)
!> \author Alin M Elena
!> \date 01/11/07, 11:49:38
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!
  subroutine SetSolutionSpace (ioLoc, genLoc, atomic, tbMod, sol)
    character (len=*), parameter :: sMyName = "SetSolutionSpace"
    type (ioType), intent (inout) :: ioLoc
    type (generalType), intent (inout) :: genLoc
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tbMod
    type (solutionType), intent (inout) :: sol
    integer :: l, n, m, z, info, i
    real (k_pr) :: cr
!
    call cpu_time (cr)
    if ( .not. sol%h%created) then
      call CreateMatrix (sol%h, atomic%basis%norbitals, .true.)
      call CreateMatrix (sol%forceOpX, sol%h%dim, .true.)
      call CreateMatrix (sol%forceOpY, sol%h%dim, .true.)
      call CreateMatrix (sol%forceOpZ, sol%h%dim, .true.)
    end if
    call ZeroMatrix (sol%h,ioLoc)
    call ZeroMatrix (sol%forceOpX,ioLoc)
    call ZeroMatrix (sol%forceOpY,ioLoc)
    call ZeroMatrix (sol%forceOpZ,ioLoc)
    n = atomic%basis%norbitals
    if ( .not. sol%eigenvecs%created) then
      allocate (sol%eigenvals(1:atomic%basis%norbitals))
      call CreateMatrix (sol%eigenvecs, atomic%basis%norbitals, .true.)
    end if
    sol%eigenvals (1:n) = 0.0_k_pr
    call ZeroMatrix (sol%eigenvecs, ioLoc)
!
    if (genLoc%scf) then
      if ( .not. sol%hin%created) then
        call CreateMatrix (sol%hin, atomic%basis%norbitals, .true.)
        call CreateMatrix (sol%h2, atomic%basis%norbitals, .true.)
        allocate (sol%potential(1:atomic%atoms%natoms))
        allocate (sol%field(1:atomic%atoms%natoms, 1:3))
      end if
      call ZeroMatrix (sol%hin,ioLoc)
      call ZeroMatrix (sol%h2,ioLoc)
      sol%potential (1:atomic%atoms%natoms) = 0.0_k_pr
      sol%field (1:atomic%atoms%natoms, 1:3) = 0.0_k_pr
    end if
!
    if (genLoc%spin) then
!spin down
      if ( .not. sol%hdown%created) then
        call CreateMatrix (sol%hdown, atomic%basis%norbitals/2, .true.)
      end if
      call ZeroMatrix (sol%hdown, ioLoc)
      if ( .not. sol%hup%created) then
        call CreateMatrix (sol%hup, atomic%basis%norbitals/2, .true.)
      end if
      call ZeroMatrix (sol%hup, ioLoc)
    end if
    if ( .not. sol%rho%created) then
      call CreateMatrix (sol%rho, atomic%basis%norbitals, .true.)
    end if
    call ZeroMatrix (sol%rho, ioLoc)
    l = maxval (tbMod%hopping(:, :)%l2)
    z = 6 * (l+1) + 1
    if (genLoc%smearMethod == k_smMP) then
      z = Max (z, genLoc%MPN)
    end if
    if ( .not. associated(sol%fact)) then
      allocate (sol%fact(0:z))
    end if
    sol%fact (0:z) = 0.0_k_pr
    call initFact (z, sol%fact)
    call PrintVectorP (sol%fact, "Factorial Values", .true., .true., ioLoc)
    if ((genLoc%scf) .and. (genLoc%electrostatics == k_electrostaticsMultipoles)) then
      call InitGaunt (2*l+1, 2*l+1, 4*l+2, sol, ioLoc)
    end if
    if ( .not. allocated(sol%n0)) then
      allocate (sol%n0(1:atomic%basis%norbitals))
    end if
    sol%n0 (1:n) = 0.0_k_pr
    call Buildn0 (genLoc, atomic, sol)
    if (ioLoc%Verbosity >= k_highVerbos) then
      write (ioLoc%uout, '(a)') "Diagonal part of the reference density:"
      write (ioLoc%uout,*)
      do l = 1, n
        write (ioLoc%uout, '(f12.8,a1,i0,a1,i0,a1,1x)', advance="no") sol%n0(l), "(", l, ",", l, ")"
      end do
      write (ioLoc%uout,*)
      write (ioLoc%uout, '(a)') "End reference density"
    end if
    l = sol%h%dim * (sol%h%dim-1) / 2
    allocate (sol%density(1:l+sol%h%dim))
    sol%density (1:l+sol%h%dim) = 0.0_k_pr
    call SetTails (ioLoc, genLoc, atomic, tbMod, sol)
    call rmarin (Int((cr-Int(cr))*3132), Int(genLoc%ranseed*30081), sol%seed, ioLoc)
    if ( .not. allocated(sol%buff%dins)) then
      n = atomic%basis%norbitals
      m = sol%hup%dim
      allocate (sol%buff%dins(1:l+n, 1:genLoc%scfMixn))
      allocate (sol%buff%douts(1:l+n, 1:genLoc%scfMixn))
      allocate (sol%buff%res(1:l+n, 1:genLoc%scfMixn))
      allocate (sol%buff%densityin(1:l+n))
      allocate (sol%buff%densityout(1:l+n))
      allocate (sol%buff%densitynext(1:l+n))
      allocate (sol%buff%tmpA(1:m))
      call CreateMatrix (sol%buff%tmpB, m, .true.)
      allocate (sol%buff%f(1:sol%eigenvecs%dim))
      allocate (sol%buff%g(1:sol%eigenvecs%dim))
      allocate (sol%buff%a(1:sol%eigenvecs%dim, 1:sol%eigenvecs%dim))
      allocate (sol%buff%pos1(1:sol%rho%dim), sol%buff%pos2(1:sol%rho%dim))
      allocate (sol%buff%nstart(1:atomic%basis%norbitals))
      allocate (sol%buff%itmp(1:atomic%atoms%natoms))
      allocate (sol%Distances(1:atomic%atoms%natoms, 1:atomic%atoms%natoms))
    end if
    call ZeroMatrix (sol%buff%tmpB, ioLoc)
    sol%buff%dins (1:l+n, 1:genLoc%scfMixn) = 0.0_k_pr
    sol%buff%douts (1:l+n, 1:genLoc%scfMixn) = 0.0_k_pr
    sol%buff%res (1:l+n, 1:genLoc%scfMixn) = 0.0_k_pr
    sol%buff%densityin (1:l+n) = 0.0_k_pr
    sol%buff%densityout (1:l+n) = 0.0_k_pr
    sol%buff%densitynext (1:l+n) = 0.0_k_pr
    sol%buff%tmpA (1:m) = 0.0_k_pr
    sol%buff%f (1:n) = 0.0_k_pr
    sol%buff%g (1:n) = 0.0_k_pr
    sol%buff%a (1:n, 1:n) = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
    sol%buff%nstart (1:n) = 0.0_k_pr
    sol%buff%pos1 (1:n) = 0
    sol%buff%pos2 (1:n) = 0
    sol%buff%itmp = 0
! allocate space for the smearing method
    if (genLoc%smearMethod == k_smMP) then
      if ( .not. associated(sol%hermite)) then
        allocate (sol%hermite(0:2*genLoc%MPN+1))
      else
        sol%hermite = 0.0_k_pr
      end if
    end if
    call BuildNeighboursList (atomic, sol, tbMod, ioLoc)
    if ( .not. allocated(sol%CurrentMatrix)) then
      allocate (sol%CurrentMatrix(1:atomic%basis%norbitals, 1:atomic%basis%norbitals))
      allocate (sol%CurrentMatrix2(1:atomic%atoms%natoms, 1:atomic%atoms%natoms))
    end if
    sol%CurrentMatrix = 0.0_k_pr
    sol%CurrentMatrix2 = 0.0_k_pr
    if ((genLoc%runType == k_runEhrenfestDamped) .or. (genLoc%runType == k_runEhrenfest)) then
      if ( .not. sol%rhodot%created) then
        call CreateMatrix (sol%rhodot, sol%h%dim, .true.)
        call CreateMatrix (sol%deltaRho, sol%h%dim, .true.)
        call CreateMatrix (sol%rhoold, sol%h%dim, .true.)
        call CreateMatrix (sol%rhonew, sol%h%dim, .true.)
        call CreateMatrix (sol%rho0, sol%h%dim, .true.)
      else
        call ZeroMatrix (sol%rhodot, ioLoc)
        call ZeroMatrix (sol%deltaRho, ioLoc)
        call ZeroMatrix (sol%rhonew, ioLoc)
        call ZeroMatrix (sol%rhoold, ioLoc)
        call ZeroMatrix (sol%rho0, ioLoc)
      end if
    end if
!
    if ( .not. sol%buff%h%created) then
      if (genLoc%spin) then
        call CreateMatrix (sol%buff%h, sol%h%dim/2, .true.)
      else
        call CreateMatrix (sol%buff%h, sol%h%dim, .true.)
      end if
    else
      call ZeroMatrix (sol%buff%h, ioLoc)
    end if
!
!
    if (( .not. genLoc%compElec) .and. (genLoc%electrostatics == k_electrostaticsMultipoles)) then
      if ( .not. allocated(sol%delq)) then
! this is the first time
        allocate (sol%delq(1:atomic%atoms%natoms), sol%vs(1:atomic%atoms%natoms), stat=info)
        if (info /= 0) then
          call error ("not enough memory to precompute electrostatic, please change to ComputeMultipolesOnFly T", sMyName, .true., &
         & ioLoc)
        end if
        do i = 1, atomic%atoms%natoms
          z = GetLmax (atomic%atoms%sp(i), atomic%speciesBasis, atomic%species)
          l = 2 * z
          sol%delq(i)%dim = l
          allocate (sol%delq(i)%a(1:l*l+2*l+1), stat=info)
          if (info /= 0) then
            call error ("not enough memory to precompute electrostatic, please change to ComputeMultipolesOnFly T", sMyName, &
           & .true., ioLoc)
          end if
          sol%delq(i)%a = 0.0_k_pr
          l = 2 * z + 1
          sol%delq(i)%dim = l
          allocate (sol%vs(i)%a(1:l*l+2*l+1), stat=info)
          if (info /= 0) then
            call error ("not enough memory to precompute electrostatic, please change to ComputeMultipolesOnFly T", sMyName, &
           & .true., ioLoc)
          end if
          sol%vs(i)%a = 0.0_k_pr
        end do
      end if
    end if

    if (.not. allocated(sol%sk%wignerD)) then
      l= maxval (atomic%speciesBasis(:, :)%l)
      allocate(sol%sk%wignerD(0:l,-l:l,0:l))
    endif
      sol%sk%wignerD=0.0_k_pr
!
  end subroutine SetSolutionSpace
!
!
!> \brief sets the tails parameters
!> \author Alin M Elena
!> \date 02/11/07, 16:18:21
!> \param ioLoc type(ioType) contains all the info about I/O files
!> \param genLoc type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains information about the tight binding model parameters
!> \param sol type(solutionType) contains information about the solution space
!
  subroutine SetTails (ioLoc, genLoc, atomic, tbMod, sol)
    character (len=*), parameter :: sMyName = "setTails"
    type (ioType), intent (inout) :: ioLoc
    type (generalType), intent (inout) :: genLoc
    type (atomicxType), intent (inout) :: atomic
    type (modelType), intent (inout) :: tbMod
    type (solutionType), intent (inout) :: sol
!
    integer :: i, j
    real (k_pr) :: f, fp, fpp
    integer :: k, k1, k2
!
!    calculate the tail function parameters
    do i = 1, atomic%species%nspecies
      do j = 1, atomic%species%nspecies
! repulsions
        f = RepNoTail (tbMod%hopping(i, j)%d1, atomic%species%id(i), atomic%species%id(j), tbMod, genLoc)
        fp = RepPNoTail (tbMod%hopping(i, j)%d1, atomic%species%id(i), atomic%species%id(j), tbMod, genLoc)
        fpp = RepPpNoTail (tbMod%hopping(i, j)%d1, atomic%species%id(i), atomic%species%id(j), tbMod, genLoc)
        tbMod%hopping(i, j)%repTail = makeTailType (f, fp, fpp, tbMod%hopping(i, j)%d1, tbMod%hopping(i, j)%dcut, ioLoc)
!
! hoppings
        if ( .not. allocated(tbMod%hopping(i, j)%hmnTail)) then
          allocate (tbMod%hopping(i, j)%hmnTail(0:tbMod%hopping(i, j)%l1, 0:tbMod%hopping(i, j)%l2, 0:tbMod%hopping(i, j)%ll))
        end if
        do k = 0, tbMod%hopping(i, j)%l1
          do k1 = 0, tbMod%hopping(i, j)%l2
            do k2 = 0, Min (k, k1)
!               write(io%uout,*) &
!               trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2)
              f = tbMod%hopping(i, j)%a(k, k1, k2) * RadNoTail (tbMod%hopping(i, j)%r1, atomic%species%id(i), atomic%species%id(j), &
             & genLoc, tbMod)
              fp = tbMod%hopping(i, j)%a(k, k1, k2) * RadPNoTail (tbMod%hopping(i, j)%r1, atomic%species%id(i), &
             & atomic%species%id(j), genLoc, tbMod)
              fpp = tbMod%hopping(i, j)%a(k, k1, k2) * RadPpNoTail (3, 3, tbMod%hopping(i, j)%r1, atomic%species%id(i), &
             & atomic%species%id(j), genLoc, tbMod)
              tbMod%hopping(i, j)%hmnTail(k, k1, k2) = makeTailType (f, fp, fpp, tbMod%hopping(i, j)%r1, tbMod%hopping(i, j)%rcut, &
             & ioLoc)
            end do
          end do
        end do
!
      end do
    end do
!
    call PrintTail (atomic, tbMod, ioLoc)
  end subroutine SetTails
!
!
!==Gaunt coefficients====================================================
!> \brief returns the Gaunt coefficient for l1,m1,l2,m2,l3,m3
!> \author Alin M Elena
!> \date 01/11/07, 17:15:07
!> \param l1,m1,l2,m2,l3,m3 angular momentum and quantum magnetic numbers
!> \param io type(ioType) contains all the info about I/O files
!> \param sol type(solutionType) contains information about the solution space
  real (k_pr) function cgaunt (l1, m1, l2, m2, l3, m3, sol, io)
    character (len=*), parameter :: myname = "cgaunt"
    integer, intent (in) :: l1, m1, l2, m2, l3, m3
    type (ioType), intent (in) :: io
    type (solutionType), intent (in) :: sol
    real (k_pr) :: sq
    integer :: g
    character (len=k_ml) :: saux
!     check for validity
    if (l1 < 0 .or. l2 < 0 .or. l3 < 0) then
      write (saux, '(a,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') "error in the input l-values in gaunt", l1, m1, l2, m2, l3, m3
      call error (trim(saux), myname, .true., io)
    end if
    if (Abs(m1) > l1 .or. Abs(m2) > l2 .or. Abs(m3) > l3) then
      write (saux, '(a,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0)') "error in the input m-values in gaunt", l1, m1, l2, m2, l3, m3
      call error (trim(saux), myname, .true., io)
    end if
!
    if ((Abs(l1-l2) <= l3) .and. (l3 <= l1+l2) .and. (Mod(l1+l2+l3, 2) == 0) .and. (m1 == m3-m2)) then
!
      g = (l1+l2+l3) / 2
      sq = real ((2*l1+1)*(2*l2+1), k_pr)
      sq = Sqrt (sq/(4.0_k_pr*k_pi)) * sol%fact(g) / (sol%fact(g-l1)*sol%fact(g-l2)*sol%fact(g-l3))
      cgaunt = (-1) ** (g-l3) * sq * cgc (l1, m1, l2, m2, l3, m3, sol) * dl (l1, l2, l3, sol)
    else
      cgaunt = 0.0_k_pr
    end if
!
  end function cgaunt
!
!> \brief the "heart" of cgaunt
!> \author Alin M Elena
!> \date 01/11/07, 17:20:53
!> \param a,x,b,y,c,z
!> \param sol type(solutionType) contains information about the solution space
  real (k_pr) function cgc (a, x, b, y, c, z, sol)
    integer, intent (in) :: a, x, b, y, c, z
    type (solutionType), intent (in) :: sol
    integer :: kmin, kmax, k
    real (k_pr) :: sum, gamma
!
    gamma = Sqrt (sol%fact(a+x)*sol%fact(b+y)*sol%fact(c+z)*sol%fact(a-x)*sol%fact(b-y)*sol%fact(c-z)*real(2*c+1, k_pr))
    sum = 0.0_k_pr
    kmin = Max (b-c-x, a+y-c, 0)
    kmax = Min (b+y, a-x, a+b-c)
    do k = kmin, kmax
      sum = sum + (-1) ** k / (sol%fact(k)*sol%fact(a+b-c-k)*sol%fact(a-x-k)*sol%fact(b+y-k)*sol%fact(c-b+x+k)*sol%fact(c-a-y+k))
    end do
    cgc = dl (a, b, c, sol) * gamma * sum
  end function cgc
!
!
  complex (k_pr) function ulmmiu (m, miu)
    integer :: m, miu, z
    real (k_pr) :: sq2, r, i
!
    sq2 = dsqrt (2.0_k_pr)
    z = 0
    r = d (m, z) * d (miu, z) + h (miu) * (d(m, miu)+(-1)**m*d(m,-miu)) / sq2
    i = h (-miu) * (-d(m,-miu)+d(m, miu)*(-1)**m) / sq2
!
    ulmmiu = cmplx (r, i, k_pr)
  end function ulmmiu
!
!> \brief Builds the Gaunt coefficients tables
!> \author Alin M Elena
!> \param l1m,l2m,l3m integers maxima ls for angular momenta
!> \param io type(ioType) contains all the info about I/O files
!> \param sol type(solutionType) contains information about the solution space
  subroutine InitGaunt (l1m, l2m, l3m, sol, io)
    character (len=*), parameter :: myname = "InitGaunt"
    integer, intent (in) :: l1m, l2m, l3m
    type (ioType), intent (in) :: io
    type (solutionType), intent (inout) :: sol
    integer :: n1, n2, n3
    integer :: i, j, k, l1, l2, l3, m1, m2, m3, ii, kk, jj, p
    complex (k_pr) :: aux, tt
!
    p = 0
    n1 = l1m * l1m + 2 * l1m + 1
    n2 = l2m * l2m + 2 * l2m + 1
    n3 = l3m * l3m + 2 * l3m + 1
    allocate (sol%gcoeff(1:n1, 1:n2, 1:n3), sol%rgc(1:n1, 1:n2, 1:n3))
!
    sol%gcoeff = 0.0_k_pr
    sol%rgc = 0.0_k_pr
!
    do l1 = 0, l1m
      do m1 = - l1, l1
        do l2 = 0, l2m
          do m2 = - l2, l2
            do l3 = 0, l3m
              do m3 = - l3, l3
                sol%gcoeff (idx(l1, m1), idx(l2, m2), idx(l3, m3)) = cgaunt (l1, m1, l2, m2, l3, m3, sol, io)
              end do
            end do
          end do
        end do
      end do
    end do
!
    if (io%Verbosity >= k_highVerbos) then
      write (io%uout, '(/a)') "===Gaunt Coefficients and Real Gaunt coefficients=================="
      write (io%uout,*) "  Gaunt R-Gaunt  l1 m1 l2 m2 l3 m3 i  j  k"
    end if
    do i = 1, n1
      call ridx (i, l1, m1)
      do j = 1, n2
        call ridx (j, l2, m2)
        do k = 1, n3
          call ridx (k, l3, m3)
          tt = cmplx (0.0_k_pr, 0.0_k_pr, k_pr)
!
          do ii = - l1, l1
            do jj = - l2, l2
              do kk = - l3, l3
                if (sol%gcoeff(idx(l1, ii), idx(l2, jj), idx(l3, kk)) /= k_zero) then
                  aux = ulmmiu (ii, m1)
                  tt = tt + (-1) ** kk * ulmmiu (ii, m1) * ulmmiu (jj, m2) * ulmmiu (-kk, m3) * sol%gcoeff(idx(l1, ii), idx(l2, &
                 & jj), idx(l3, kk))
                end if
              end do
            end do
          end do
          sol%rgc (i, j, k) = real (tt, k_pr)
          if (io%Verbosity >= k_highVerbos) then
            if ((sol%gcoeff(i, j, k) /= k_zero) .or. (sol%rgc(i, j, k) /= k_zero)) then
              write (io%uout, '(2F8.4,6I3,3I3)') sol%gcoeff(i, j, k), sol%rgc(i, j, k), l1, m1, l2, m2, l3, m3, i, j, k
              p = p + 1
              if (Mod(p, 40) == 0) write (io%uout,*) "  Gaunt R-Gaunt  l1 m1 l2 m2 l3 m3 i  j  k"
            end if
          end if
        end do
      end do
    end do
    if (io%Verbosity >= k_highVerbos) then
      write (io%uout, '(/a/)') "==End Gaunt Coefficients============================================="
    end if
  end subroutine InitGaunt
! =====End Gaunt Coefficients ================================================================================
! ====Radial dependence=========================================================================================
  function rad (r, sp1, sp2, gen, tb, l1, l2, m)
! Gives you the radial dependence for elemental overlaps
! for specie sp1, sp2
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'rad'
!--subroutine parameters -------------------------!
    real (k_pr), intent (inout) :: r
    integer, intent (inout) :: sp1, sp2
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
    integer, intent (in) :: l1, l2, m
    real (k_pr) :: rad
!-------------------------------------------------!
    if (r <= (tb%hopping(sp1, sp2)%r1)) then
      rad = tb%hopping(sp1, sp2)%a(l1, l2, m) * RadNoTail (r, sp1, sp2, gen, tb)
    else if (r <= (tb%hopping(sp1, sp2)%rcut)) then
      rad = tailFunction (r, tb%hopping(sp1, sp2)%hmnTail(l1, l2, m))
    else
      rad = 0.0_k_pr
    end if
  end function rad
!
  function RadNoTail (r, sp1, sp2, gen, tb)
!Gives you the radial dependence for elemental overlaps
! for specie sp1, sp2
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'rad'
!--subroutine parameters -------------------------!
    real (k_pr), intent (in) :: r
    integer, intent (in) :: sp1, sp2
    type (modelType), intent (in) :: tb
    type (generalType), intent (in) :: gen
!
    real (k_pr) :: RadNoTail
!-------------------------------------------------!
    real (k_pr) :: r0, n, nc, rc
    select case (gen%bond)
    case (k_bondGSP)
      r0 = tb%hopping(sp1, sp2)%r0
      n = tb%hopping(sp1, sp2)%n
      nc = tb%hopping(sp1, sp2)%nc
      rc = tb%hopping(sp1, sp2)%rc
      RadNoTail = (r0/r) ** n * Exp (n*(-(r/rc)**nc+(r0/rc)**nc))
    case (k_bondHarrison)
      n = tb%hopping(sp1, sp2)%n
      RadNoTail = k_hbar ** 2 / k_me / r ** (n)
    end select
!
  end function RadNoTail
!
!
  function RadP (r, sp1, sp2, gen, tb, l1, l2, m)
!Gives you the radial dependence for elemental overlaps
! for specie sp1, sp2
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RadP'
!--subroutine parameters -------------------------!
    real (k_pr), intent (inout) :: r
    integer :: sp1, sp2!, alpha
    integer, intent (in) :: l1, l2, m
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
    real (k_pr) :: RadP
!-------------------------------------------------!
!
!
    if (r <= tb%hopping(sp1, sp2)%r1) then
      RadP = tb%hopping(sp1, sp2)%a(l1, l2, m) * RadPNoTail (r, sp1, sp2, gen, tb)
    else if (r <= tb%hopping(sp1, sp2)%rcut) then
      RadP = TailFunctionP (r, tb%hopping(sp1, sp2)%hmnTail(l1, l2, m))
    else
      RadP = 0.0_k_pr
    end if
!
!
  end function RadP
!
!
  function RadPNoTail (r, sp1, sp2, gen, tb)
!Gives you the radial dependence for elemental overlaps
! for specie sp1, sp2
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RadP_notail'
!--subroutine parameters -------------------------!
    real (k_pr), intent (inout) :: r
    integer, intent (inout) :: sp1, sp2
    !integer, intent (in) :: alpha
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
    real (k_pr) :: RadPNoTail
!-------------------------------------------------!
    real (k_pr) :: r0, n, nc, rc
!
    select case (gen%bond)
    case (k_bondGSP)
      r0 = tb%hopping(sp1, sp2)%r0
      n = tb%hopping(sp1, sp2)%n
      nc = tb%hopping(sp1, sp2)%nc
      rc = tb%hopping(sp1, sp2)%rc
      RadPNoTail = - RadNoTail (r, sp1, sp2, gen, tb) * (1.0_k_pr+nc*(r/rc)**nc) * (n/r)
    case (k_bondHarrison)
      n = tb%hopping(sp1, sp2)%n
      RadPNoTail = - n / r ** (n+1) * k_hbar ** 2 / k_me
    end select
  end function RadPNoTail
!
  function RadPp (alpha, beta, r, sp1, sp2, gen, tb, l1, l2, m)
!Gives you the radial dependence for elemental overlaps
! for specie sp1, sp2
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RadPpNoTail'
!--subroutine parameters -------------------------!
    real (k_pr), intent (inout) :: r
    integer, intent (inout) :: sp1, sp2, alpha, beta
    type (modelType), intent (inout) :: tb
    integer, intent (in) :: l1, l2, m
    type (generalType), intent (inout) :: gen
    real (k_pr) :: RadPp
!-------------------------------------------------!
!
    if (r <= tb%hopping(sp1, sp2)%r1) then
      RadPp = tb%hopping(sp1, sp2)%a(l1, l2, m) * RadPpNoTail (alpha, beta, r, sp1, sp2, gen, tb)
    else if (r <= tb%hopping(sp1, sp2)%rcut) then
      if ((alpha == 3) .and. (beta == 3)) then
        RadPp = TailFunctionPp (r, tb%hopping(sp1, sp2)%hmnTail(l1, l2, m))
      else
        RadPp = 0.0_k_pr
      end if
    else
      RadPp = 0.0_k_pr
    end if
!
!
!
  end function RadPp
!
  function RadPpNoTail (alpha, beta, r, sp1, sp2, gen, tb)
!Gives you the radial dependence for elemental overlaps
! for specie sp1, sp2
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RadPpNoTail'
!--subroutine parameters -------------------------!
    real (k_pr), intent (inout) :: r
    integer, intent (inout) :: sp1, sp2
    integer, intent (in) :: alpha, beta
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
    real (k_pr) :: RadPpNoTail
!-------------------------------------------------!
    real (k_pr) :: r0, n, nc, rc
!
    select case (gen%bond)
    case (k_bondGSP)
      r0 = tb%hopping(sp1, sp2)%r0
      n = tb%hopping(sp1, sp2)%n
      nc = tb%hopping(sp1, sp2)%nc
      rc = tb%hopping(sp1, sp2)%rc
      if ((alpha == 3) .and. (beta == 3)) then
        RadPpNoTail = RadNoTail (r, sp1, sp2, gen, tb) * &
       & ((-n/r-n*nc*(r/rc)**nc/r)**2-(-n/(r*r)+nc*(nc-1.0_k_pr)*n/(r*r)*(r/rc)**nc))
      else
        RadPpNoTail = 0.0_k_pr
      end if
    case (k_bondHarrison)
      if ((alpha == 3) .and. (beta == 3)) then
        n = tb%hopping(sp1, sp2)%n
        RadPpNoTail = n * (n+1.0_k_pr) / r ** (n+2) * k_hbar ** 2 / k_me
      else
        RadPpNoTail = 0.0_k_pr
      end if
    end select
  end function RadPpNoTail
! ============end radial dependence==============================================================================
! ====repulsive term======================================================================================
  function rep (r, sp1, sp2, tb, gen)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'rep'
!--subroutine parameters -------------------------!
    real (k_pr) :: rep
    real (k_pr), intent (in) :: r
    integer, intent (in) :: sp1, sp2
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
!--internal variables ----------------------------!
!    real(k_pr) :: phi0,d0,m,dc,mc
!-------------------------------------------------!
!
!
    if (r <= tb%hopping(sp1, sp2)%d1) then
      rep = RepNoTail (r, sp1, sp2, tb, gen)
    else if (r <= tb%hopping(sp1, sp2)%dcut) then
      rep = tailFunction (r, tb%hopping(sp1, sp2)%repTail)
    else
      rep = 0.0_k_pr
    end if
!
  end function rep
!
!
  function RepNoTail (r, sp1, sp2, tb, gen)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RepNoTail'
!--subroutine parameters -------------------------!
    real (k_pr) :: RepNoTail
    real (k_pr), intent (in) :: r
    integer, intent (in) :: sp1, sp2
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
!--internal variables ----------------------------!
    real (k_pr) :: phi0, d0, m, dc, mc
!-------------------------------------------------!
    select case (gen%bond)
    case (k_bondGSP)
      phi0 = tb%hopping(sp1, sp2)%phi0
      d0 = tb%hopping(sp1, sp2)%d0
      m = tb%hopping(sp1, sp2)%m
      mc = tb%hopping(sp1, sp2)%mc
      dc = tb%hopping(sp1, sp2)%dc
      RepNoTail = phi0 * (d0/r) ** m * Exp (m*(-(r/dc)**mc+(d0/dc)**mc))
!
    case (k_bondHarrison)
      phi0 = tb%hopping(sp1, sp2)%phi0
      m = tb%hopping(sp1, sp2)%m
      RepNoTail = phi0 / r ** m
    end select
  end function RepNoTail
!
  function RepP (r, sp1, sp2, tb, gen)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RepP'
!--subroutine parameters -------------------------!
    real (k_pr) :: RepP
    real (k_pr), intent (in) :: r
    integer, intent (in) :: sp1, sp2
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
!--internal variables ----------------------------!
!    real(k_pr) :: phi0,d0,m,dc,mc
!-------------------------------------------------!
!
    if (r <= tb%hopping(sp1, sp2)%d1) then
      RepP = RepPNoTail (r, sp1, sp2, tb, gen)
    else if (r <= tb%hopping(sp1, sp2)%dcut) then
      RepP = TailFunctionP (r, tb%hopping(sp1, sp2)%repTail)
    else
      RepP = 0.0_k_pr
    end if
!
  end function RepP
!
  function RepPNoTail (r, sp1, sp2, tb, gen)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RepPNoTail'
!--subroutine parameters -------------------------!
    real (k_pr) :: RepPNoTail
    real (k_pr), intent (in) :: r
    integer, intent (in) :: sp1, sp2
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
!--internal variables ----------------------------!
    real (k_pr) :: phi0, d0, m, dc, mc
!-------------------------------------------------!
    select case (gen%bond)
    case (k_bondGSP)
      phi0 = tb%hopping(sp1, sp2)%phi0
      d0 = tb%hopping(sp1, sp2)%d0
      m = tb%hopping(sp1, sp2)%m
      mc = tb%hopping(sp1, sp2)%mc
      dc = tb%hopping(sp1, sp2)%dc
!
      RepPNoTail = RepNoTail (r, sp1, sp2, tb, gen) * (-m/r-m*mc*(r/dc)**mc/r)
    case (k_bondHarrison)
!
      phi0 = tb%hopping(sp1, sp2)%phi0
      m = tb%hopping(sp1, sp2)%m
!
      RepPNoTail = - m * phi0 / r ** (m+1)
    end select
  end function RepPNoTail
!
!
  function RepPp (r, sp1, sp2, tb, gen)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RepPp'
!--subroutine parameters -------------------------!
    real (k_pr) :: RepPp
    real (k_pr), intent (in) :: r
    integer, intent (in) :: sp1, sp2
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
!--internal variables ----------------------------!
!    real(k_pr) :: phi0,d0,m,dc,mc
!-------------------------------------------------!
!
    if (r <= tb%hopping(sp1, sp2)%d1) then
      RepPp = RepPpNoTail (r, sp1, sp2, tb, gen)
    else if (r <= tb%hopping(sp1, sp2)%dcut) then
      RepPp = TailFunctionPp (r, tb%hopping(sp1, sp2)%repTail)
    else
      RepPp = 0.0_k_pr
    end if
!
  end function RepPp
!
!
  function RepPpNoTail (r, sp1, sp2, tb, gen)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'RepPpNoTail'
!--subroutine parameters -------------------------!
    real (k_pr) :: RepPpNoTail
    real (k_pr), intent (in) :: r
    integer, intent (in) :: sp1, sp2
    type (modelType), intent (inout) :: tb
    type (generalType), intent (inout) :: gen
!--internal variables ----------------------------!
    real (k_pr) :: phi0, d0, m, dc, mc
!-------------------------------------------------!
    select case (gen%bond)
    case (k_bondGSP)
      phi0 = tb%hopping(sp1, sp2)%phi0
      d0 = tb%hopping(sp1, sp2)%d0
      m = tb%hopping(sp1, sp2)%m
      mc = tb%hopping(sp1, sp2)%mc
      dc = tb%hopping(sp1, sp2)%dc
!
      RepPpNoTail = RepNoTail (r, sp1, sp2, tb, gen) * ((-m/r-m*mc*(r/dc)**mc/r)**2-(-m/(r*r)+mc*(mc-1.0_k_pr)*m/(r*r)*(r/dc)**mc))
    case (k_bondHarrison)
!
      phi0 = tb%hopping(sp1, sp2)%phi0
      m = tb%hopping(sp1, sp2)%m
!
      RepPpNoTail = m * (m+1.0_k_pr) * phi0 / r ** (m+2)
    end select
  end function RepPpNoTail
! ========end  repulsive term ============================================================================
!======== onsite term======================
  function onsite (orba, orbb, tb)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = "onsite"
!--subroutine parameters -------------------------!
    real (k_pr) :: onsite
    type (orbitalType), intent (in) :: orba, orbb
    type (modelType), intent (in) :: tb
!--internal variables ----------------------------!
!
    if ((orba%l == orbb%l) .and. (orba%m == orbb%m)) then
      onsite = tb%hopping(orba%sp, orba%sp)%eps(orba%l)
    else
      onsite = 0.0_k_pr
    end if
  end function onsite
!=====end onsite term ======================
!========embedding term=========================================
  function embedding (x, sp, tb)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = "embedding"
!--subroutine parameters -------------------------!
    real (k_pr) :: embedding
    real (k_pr), intent (in) :: x
    integer, intent (in) :: sp
    type (modelType), intent (in) :: tb
!--internal variables ----------------------------!
    real (k_pr) :: a1, a2, a3, a4
!
    a1 = tb%hopping(sp, sp)%a1
    a2 = tb%hopping(sp, sp)%a2
    a3 = tb%hopping(sp, sp)%a3
    a4 = tb%hopping(sp, sp)%a4
!
    embedding = a1 * x + a2 * x * x + a3 * x * x * x + a4 * x * x * x * x
!
  end function embedding
!
  function EmbeddingP (x, sp, tb)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = "embedding"
!--subroutine parameters -------------------------!
    real (k_pr) :: EmbeddingP
    real (k_pr), intent (in) :: x
    integer, intent (in) :: sp
    type (modelType), intent (in) :: tb
!--internal variables ----------------------------!
    real (k_pr) :: a1, a2, a3, a4
!
    a1 = tb%hopping(sp, sp)%a1
    a2 = tb%hopping(sp, sp)%a2
    a3 = tb%hopping(sp, sp)%a3
    a4 = tb%hopping(sp, sp)%a4
    EmbeddingP = a1 + 2.0_k_pr * a2 * x + 3.0_k_pr * a3 * x * x + 4.0_k_pr * a4 * x * x * x
!
  end function EmbeddingP
!
  function EmbeddingPp (x, sp, tb)
!--subroutine name--------------------------------!
    character (len=*), parameter :: myname = 'EmbeddingPp'
!--subroutine parameters -------------------------!
    real (k_pr) :: EmbeddingPp
    real (k_pr), intent (in) :: x
    integer, intent (in) :: sp
    type (modelType), intent (in) :: tb
!--internal variables ----------------------------!
    real (k_pr) :: a2, a3, a4
!
!
    a2 = tb%hopping(sp, sp)%a2
    a3 = tb%hopping(sp, sp)%a3
    a4 = tb%hopping(sp, sp)%a4
!
    EmbeddingPp = 2.0_k_pr * a2 + 6.0_k_pr * a3 * x + 12.0_k_pr * a4 * x ** 2
!
  end function EmbeddingPp
!
!========end embedding term ====================================
!> \brief initializes the magnetic moment
!> \details it just computes the difference between the spin up and spin down charges
!> \author Alin M Elena
!> \date 08/11/07, 11:39:49
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
  subroutine InitMagneticMoment (atomic)
    character (len=*), parameter :: myname = "InitMagneticMoment"
    type (atomicxType), intent (inout) :: atomic
!
    real (k_pr) :: tmom, lmom
    integer :: i, k
    tmom = 0.0_k_pr
!
    do i = 1, atomic%atoms%natoms
      lmom = 0.0_k_pr
      do k = 1, atomic%atoms%norbs(i) / 2
        lmom = lmom - (atomic%basis%orbitals(atomic%atoms%orbs(i, k))%occup-atomic%basis%orbitals(atomic%atoms%orbs(i, &
       & k+atomic%atoms%norbs(i)/2))%occup)
      end do
      atomic%atoms%MagMom (i) = lmom
      tmom = tmom + atomic%atoms%MagMom(i)
    end do
    atomic%atoms%MagneticMoment = tmom
  end subroutine InitMagneticMoment
!
!> \brief computes the n0
!> \details we suppose that we start with a diagonal "density matrix"
!> \author Alin M Elena
!> \date 08/11/07, 10:41:28
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
  subroutine Buildn0 (gen, atomic, sol)
    character (len=*), parameter :: myname = 'Buildn0'
    type (generalType), intent (inout) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    integer :: i, k, j, m
!
    if (gen%spin) then
      m = atomic%basis%norbitals / 2
      do i = 1, atomic%atoms%natoms
        do k = 1, atomic%species%norbs(atomic%atoms%sp(i)) / 2
          j = atomic%atoms%orbs(i, k)
          if (gen%SymRefRho) then
            sol%n0 (j) = (atomic%basis%orbitals(j)%occup+atomic%basis%orbitals(j+m)%occup) / 2.0_k_pr - gen%netCharge / &
           & atomic%basis%norbitals
            sol%n0 (j+m) = sol%n0(j)
          else
            sol%n0 (j) = atomic%basis%orbitals(j)%occup - gen%netCharge / atomic%basis%norbitals
            sol%n0 (j+m) = atomic%basis%orbitals(j+m)%occup - gen%netCharge / atomic%basis%norbitals
          end if
        end do
      end do
    else
      do i = 1, atomic%atoms%natoms
        do k = 1, atomic%species%norbs(atomic%atoms%sp(i))
          j = atomic%atoms%orbs(i, k)
          sol%n0 (j) = atomic%basis%orbitals(j)%occup - gen%netCharge / atomic%basis%norbitals
        end do
      end do
    end if
  end subroutine Buildn0
!
!> \brief computes the magnetic pmoment for each atom
!> \details internally all the "hard" work is done by LocalMoment
!> \author Alin M Elena
!> \date 08/11/07, 11:41:42
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param io type(ioType) contains all the info about I/O files
  subroutine ComputeMagneticMoment (gen, atomic, sol, io)
    character (len=*), parameter :: myname = 'ComputeMagneticMoment'
    type (generalType), intent (in) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    integer :: i
    real (k_pr) :: tmom
    tmom = 0.0_k_pr
    if (gen%spin) then
      do i = 1, atomic%atoms%natoms
        atomic%atoms%MagMom (i) = LocalMoment (atomic%atoms%id(i), atomic, .false., io, sol)
        tmom = tmom + atomic%atoms%MagMom(i)
      end do
      atomic%atoms%MagneticMoment = tmom
    else
      atomic%atoms%MagMom = 0.0_k_pr
      atomic%atoms%MagneticMoment = 0.0_k_pr
    end if
  end subroutine ComputeMagneticMoment
!
  subroutine BuildNeighboursList (atomic, sol, tb, io)
    character (len=*), parameter :: myname = "BuildNeighboursList"
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    type (modelType), intent (inout) :: tb
    integer :: i, k
!
    call CountNeighbours (atomic, sol, tb, io)
    do i = 1, atomic%atoms%ncurrent
      k = atomic%atoms%current(i)
      atomic%atoms%neighbours(k)%n = sol%buff%itmp(k)
      allocate (atomic%atoms%neighbours(k)%a(1:sol%buff%itmp(k)))
      call GetNeighbours (k, atomic, sol, tb, io)
    end do
    call PrintNeighbours (atomic, io)
  end subroutine BuildNeighboursList
!
  subroutine CountNeighbours (atomic, sol, tb, io)
    character (len=*), parameter :: myname = "CountNeighbours"
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    type (modelType), intent (inout) :: tb
    integer :: i, j, sp1, sp2
!
    sol%buff%itmp = 0
    call ComputeEuclideanMatrix (atomic%atoms, io, sol%Distances)
    if (atomic%atoms%ncurrent == atomic%atoms%natoms) then
      do i = 1, atomic%atoms%ncurrent - 1
        sp1 = atomic%atoms%sp(atomic%atoms%current(i))
        do j = i + 1, atomic%atoms%natoms
          sp2 = atomic%atoms%sp(j)
          if (sol%Distances(atomic%atoms%current(i), j) <= tb%hopping(sp1, sp2)%rcut) then
            sol%buff%itmp (atomic%atoms%current(i)) = sol%buff%itmp(atomic%atoms%current(i)) + 1
            sol%buff%itmp (j) = sol%buff%itmp(j) + 1
          end if
        end do
      end do
    else
      do i = 1, atomic%atoms%ncurrent
        sp1 = atomic%atoms%sp(atomic%atoms%current(i))
        do j = 1, atomic%atoms%natoms
          sp2 = atomic%atoms%sp(j)
          if ((atomic%atoms%current(i) /= j) .and. (sol%Distances(atomic%atoms%current(i), j) <= tb%hopping(sp1, sp2)%rcut)) then
            sol%buff%itmp (atomic%atoms%current(i)) = sol%buff%itmp(atomic%atoms%current(i)) + 1
          end if
        end do
      end do
    end if
  end subroutine CountNeighbours
!
  subroutine GetNeighbours (at, atomic, sol, tb, io)
    character (len=*), parameter :: myname = "GetNeighbours"
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    type (modelType), intent (inout) :: tb
    integer, intent (inout) :: at
    integer :: i, j, sp1, sp2, k
!
!
    sp1 = atomic%atoms%sp(at)
    k = 1
    do j = 1, atomic%atoms%natoms
      sp2 = atomic%atoms%sp(j)
      if ((at /= j) .and. (sol%Distances(at, j) <= tb%hopping(sp1, sp2)%rcut) .and. k <= atomic%atoms%neighbours(at)%n) then
        atomic%atoms%neighbours(at)%a(k) = j
        k = k + 1
      end if
    end do
!
  end subroutine GetNeighbours
!
  subroutine UpdateNeighboursList (atomic, sol, tb, io)
    character (len=*), parameter :: myname = "UpdateNeighboursList"
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    type (modelType), intent (inout) :: tb
    integer :: i, k
!
    call CountNeighbours (atomic, sol, tb, io)
    do i = 1, atomic%atoms%ncurrent
      k = atomic%atoms%current(i)
      if (sol%buff%itmp(k) > atomic%atoms%neighbours(k)%n) then
        atomic%atoms%neighbours(k)%n = sol%buff%itmp(k)
        deallocate (atomic%atoms%neighbours(k)%a)
        allocate (atomic%atoms%neighbours(k)%a(1:sol%buff%itmp(k)))
      else
        atomic%atoms%neighbours(k)%n = sol%buff%itmp(k)
      end if
      call GetNeighbours (k, atomic, sol, tb, io)
    end do
    if (io%Verbosity >= k_highVerbos) then
      call PrintNeighbours (atomic, io)
    end if
  end subroutine UpdateNeighboursList
!
  real (k_pr) function ComputeBondCurrentsOnAtom (atp, at, gen, atomic, sol, io, lOnOrbitals)
    character (len=*), parameter :: myname = 'ComputeBondCurrentsOnAtom'
    type (generalType), intent (in) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    integer, intent (in) :: at, atp
    logical, intent (in) :: lOnOrbitals
!
    integer :: m, n, i, j
    real (k_pr) :: aux, aux2
!
    aux = 0.0_k_pr
    do n = 1, atomic%species%norbs(atomic%atoms%sp(atp))
      i = atomic%atoms%orbs(atp, n)
      do m = 1, atomic%species%norbs(atomic%atoms%sp(at))
        j = atomic%atoms%orbs(at, m)
        aux2 = real (sol%h%a(j, i)) * aimag (sol%rho%a(i, j))
        aux = aux + aux2
        if (lOnOrbitals) then
          sol%CurrentMatrix (i, j) = aux2
        end if
      end do
    end do
    ComputeBondCurrentsOnAtom = 2.0_k_pr * aux * k_e / k_hbar
  end function ComputeBondCurrentsOnAtom
!
!
  subroutine ComputeBondCurrents (gen, atomic, sol, io, lOnOrbitals)
    character (len=*), parameter :: myname = 'ComputeBondCurrents'
    type (generalType), intent (in) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    logical, intent (in) :: lOnOrbitals
    integer :: i, j, at, atp
    real (k_pr) :: inpn, aux, aux2
!

    do i = 1, atomic%atoms%ncurrent
      atp = atomic%atoms%current(i)
      aux = 0.0_k_pr
      do j = 1, atomic%atoms%neighbours(atp)%n
        at = atomic%atoms%neighbours(atp)%a(j)
        inpn = ComputeBondCurrentsOnAtom (atp, at, gen, atomic, sol, io, lOnOrbitals)
        sol%CurrentMatrix2 (atp, at) = inpn
        aux = aux + inpn
      end do
      sol%CurrentMatrix2 (atp, atp) = aux
    end do
  end subroutine ComputeBondCurrents
!
!
  real (k_pr) function ComputeBondCurrentsOnOrbital (atp, at, gen, atomic, sol, io, l, m)
    character (len=*), parameter :: myname = 'ComputeBondCurrentsOnOrbital'
    type (generalType), intent (in) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    integer, intent (inout) :: at, atp
    integer, intent (inout) :: l, m
!
    integer :: i, j
    real (k_pr) :: aux, aux2
!
    if (GetLmax(atomic%atoms%sp(atp), atomic%speciesBasis, atomic%species) < l) then
      ComputeBondCurrentsOnOrbital = 0.0_k_pr
      return
    end if
    if (GetLmax(atomic%atoms%sp(at), atomic%speciesBasis, atomic%species) < l) then
      ComputeBondCurrentsOnOrbital = 0.0_k_pr
      return
    end if
    i = getOrbitalIndex (atp, atomic, l, m)
    j = getOrbitalIndex (at, atomic, l, m)
    if (gen%spin) then
      aux = sol%CurrentMatrix (i, j) + sol%CurrentMatrix(i+atomic%basis%norbitals/2, j+atomic%basis%norbitals/2)
    else
      aux = sol%CurrentMatrix (i, j)
    end if
    ComputeBondCurrentsOnOrbital = 2.0_k_pr * aux * k_e / k_hbar
  end function ComputeBondCurrentsOnOrbital
!
  subroutine ComputeBondCurrentsOnOrbitals (gen, atomic, sol, io, l, m)
    character (len=*), parameter :: myname = 'ComputeBondCurrentsOnOrbitals'
    type (generalType), intent (in) :: gen
    type (atomicxType), intent (inout) :: atomic
    type (solutionType), intent (inout) :: sol
    type (ioType), intent (inout) :: io
    integer, intent (inout) :: l, m
    integer :: i, j, at, atp
    real (k_pr) :: inpn, aux
!
    do i = 1, atomic%atoms%ncurrent
      atp = atomic%atoms%current(i)
      aux = 0.0_k_pr
      do j = 1, atomic%atoms%neighbours(atp)%n
        at = atomic%atoms%neighbours(atp)%a(j)
        inpn = ComputeBondCurrentsOnOrbital (atp, at, gen, atomic, sol, io, l, m)
        sol%CurrentMatrix2 (atp, at) = inpn
        aux = aux + inpn
      end do
      sol%CurrentMatrix2 (atp, atp) = aux
    end do
  end subroutine ComputeBondCurrentsOnOrbitals
end module m_TightBinding
