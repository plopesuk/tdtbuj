!> \brief subprograms of general use
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
module m_Useful
  use m_Constants
  use m_Types

  implicit none
  private
  public :: cstr
  public :: error
  public :: GetUnit
  public :: isUnique
  public :: isInList
  public :: DateAndTime
  public :: ccvar
  public :: ccnlm
  public :: rmarin
  public :: ranmar
  public :: InitFact
  public :: h
  public :: d
  public :: GetZ
  public :: GetSpecie
  public :: LocalMoment
  public :: GetLmax
  public :: AtomDistance
  public :: AtomDistVec
  public :: Distance
  public :: Fermi
  public :: ka
  public :: aidx
  public :: norm
  public :: VectorModulus
  public :: iMaxLoc
  public :: iMinLoc
  public :: Swap
  public :: AssertEq
  public :: ExpRep
  public :: ComputeEuclideanMatrix
  public :: InitializeHermite
  public :: occupMP
  public :: sn
  public :: MarzariF
  public :: MarzariS
  public :: LMax
  public :: getUnits
  interface Swap
    module procedure SwapScalar,SwapVector
  end interface

  interface AssertEq
    module procedure AssertEq3, AssertEq4
  end interface

contains

!> \brief logical function compares two strings
!> \details not case sensitive returns true if the strings are the same false otherwise
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \param str1, str2 character(len=*) that get compared

  logical function cstr(str1,str2)
    character(len=*), parameter :: myname = 'cstr'
    character(len=*) :: str1, str2
    integer :: len1,len2,i,s

    len1=len(str1)
    len2=len(str2)
    cstr=.false.
    if (len1==len2) then
      s=0
      do i=1,len1
        s=s+abs(up(str1(i:i))-up(str2(i:i)))
      enddo
      if (s==0) cstr=.true.
    endif
  end function cstr

!> \brief integer function returns the ascii code of a letter
!> \details always will be the ascii code of the capital letter
!> as long as the characters are grouped in the set consecutive it should work
!> (A-Z a-z)
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th of January, 2006
!> \param a character
  integer function up(a)
    character(len=1) :: a
    if (iachar(a)>iachar("Z")) then
      up=iachar(a)-(iachar("a")-iachar("A"))
    else
      up=iachar(a)
    endif
  end function up

!> \brief Prints an error/warning message and aborts the programs if necessary
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th January 2006
!> \param message the message that will be displayed
!> \param routine the name of the caller
!> \param ioLoc I/O details
!> \param critical logical if true aborts the program
!> \remarks  20th of January, 2007, by Alin M Elena (Queen's University Belfast), added ioLoc parameter
  subroutine error(message, routine, critical,ioLoc)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: routine
    logical, intent(in)   :: critical
    type(ioType), intent(in) :: ioLoc

    if (critical) then
      write(ioLoc%uout,*) &
        "Critical error in subroutine: ", routine
    else
      write(ioLoc%uout,*) &
        "Error message from subroutine: ", routine
    end if

    write(ioLoc%uout,*) routine,": ", message

    if (critical) then
      write(ioLoc%uout,*)routine,": User stop."
      write(*,*)routine," User stop."
      stop
    endif

  end subroutine error

!> \brief gives an integer that can be used as a unit number in a open command
!> \details this function should be used to avoid using the same unit in different open commands
!> \author Alin M Elena
!> \date 14-15th of January, 2006
!> \warning using a save attribute makes the function thread unsafe, however
!> I do not expect any problems from it (you do not open so many files after all)
  integer function GetUnit()
    character(len=*), parameter :: myname = 'GetUnit()'
    integer,save ::ustart=9

    ustart=ustart+1
    GetUnit=ustart

  end function GetUnit

!> \brief checks if the elements of an integer array are unique
!> \author Alin M Elena
!> \date 29/10/07, 19:32:47
!> \param x integer array
  logical function isUnique(x)
    character(len=*), parameter :: sMyName="isUnique"
    integer,intent(in) :: x(:)
    integer :: i,j
    integer :: n
    logical :: aux

    n=size(x)
    aux=.false.
    if (n/=1) then
      main: do i=1,n-1
        do j=i+1,n
          if (x(i) == x(j)) then
            aux = .true.
            exit main
          endif
        enddo
      enddo main
    endif
    isUnique=aux
  end function isUnique

!> \brief checks if the elements of an integer is in a list
!> \author Alin M Elena
!> \date 29/10/07, 23:50:47
!> \param x integer array
!> \param y integer number to be searched
!> \param pos integer optional returns the position where y is found or -1 if is not found
  logical function isInList(y,x,pos)
    character(len=*), parameter :: sMyName="isInList"
    integer, intent(in) :: x(:)
    integer,intent(in) :: y
    integer, intent(inout), optional :: pos

    integer :: i
    integer :: n
    logical :: aux

    n=size(x)
    if (present(pos)) pos = -1
    aux=.false.
    if (n>=1) then
      do i=1,n
        if (y == x(i)) then
          aux = .true.
          if (present(pos)) pos = i
        endif
      enddo
    endif
    isInList=aux
  end function isInList

!> \brief return the date and time in a human format
!> \author Alin M Elena
!> \date 31/10/07, 09:48:57
!> \param date character contains the date in the format dd-mm-yyyy
!> \param time character contains the time in the format hh:mm:ss.mmmgoo
  subroutine DateAndTime(date,time)
    character(len=*), parameter :: sMyName="DateAndTime"
    character(len=*), intent(out) :: date, time
    character(len=8) :: dt
    character(len=10) :: tm

    call date_and_time(dt,tm)
    date(7:10)=dt(1:4)
    date(3:3)="-"
    date(4:5)=dt(5:6)
    date(6:6)="-"
    date(1:2)=dt(7:8)
    time(1:2)=tm(1:2)
    time(3:3)=":"
    time(4:5)=tm(3:4)
    time(6:6)=":"
    time(7:12)=tm(5:10)
  end subroutine DateAndTime

!> \brief creates a string having as input two integers and a string by concatenation and removing blanks
!> \author Alin M. Elena, Queen's University Belfast, UK
!> \date 20th of January, 2005
!> \param i,j integers to be concatenated
!> \param strr string
!> \remarks rewrote on 31st of October 2007

  character(len=30) function ccvar(i,j,strr)
    character(len=*),parameter :: myname="ccvar"
    integer, intent(in) :: i,j
    character(len=10), intent(in) :: strr
    character(len=k_ml)   :: dump

    write(dump,'(a,i0,i0)')trim(strr),i,j
    ccvar=trim(dump)
  end function ccvar

!> \brief creates a string having as input five integers by concatenation and removing blanks
!> \author Alin M. Elena, Queen's University Belfast, UK
!> \date 20th of January, 2005
!> \param i,j,k,k1,k2 integers to be concatenated
!> \remarks rewrote on 31st of October 2007
  character(len=30) function ccnlm(i,j,k,k1,k2)
    character(len=*),parameter :: myname="ccnlm"
    integer, intent(in) :: i,j,k,k1,k2
    character(len=30) ::dump

    write(dump,'(a,i0,a,i0,a,i0,a,i0,i0)')"l",k,"l",k1,"m",k2,"_",i,j
    ccnlm=trim(dump)
  end function ccnlm

!> \brief Initializes the random number generator
!> \details returns an array of "seeds"
!> \author netlib
!> \param ij,kl starting seeds
!> \param u real array of "seeds
!> \param io type(ioType) i/o units
!> \remarks
!> here is a sequence of code that will Initialize the generator
!> \code
!> real(k_pr) :: cr,seed(97)
!>  call cpu_time(cr)
!>  call rmarin(int((cr-int(cr))*3132),general%ranseed*30081),seed,io)
!> \endcode
  subroutine rmarin(ij,kl,u,io)
    character(len=*),parameter :: myname="rmarin"
    real(k_pr),intent(inout) ::  u(:)
    integer, intent(in) :: ij,kl
    type(ioType),intent(in) :: io
    character(len=k_ml) :: saux
    integer :: i,j,k,l,ii,jj,m
    real(k_pr) :: s,t

    if( ij < 0  .or.  ij > 31328  .or. &
        kl < 0  .or.  kl > 30081 ) then
      write(saux,'(a,i0,a,i0)')"first seed must be between 0 and 31328 the second seed must be between 0 and  30081! you supplied "&
         ,ij," and ",kl
      call error(trim(saux),myname,.true.,io)
    endif
    i = mod(ij/177, 177) + 2
    j = mod(ij    , 177) + 2
    k = mod(kl/169, 178) + 1
    l = mod(kl,     169)
    do ii = 1, 97
        s = 0.0_k_pr
        t = 0.5_k_pr
        do jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) >= 32) then
              s = s + t
            endif
            t = 0.5_k_pr * t
        enddo
        u(ii) = s
    enddo

  end subroutine rmarin

!> \brief generates a random number
!> \author netlib
!> \param u(:) real array of "seeds"
  real(k_pr) function ranmar(u)
    real(k_pr), intent(inout) :: u(:)
    real(k_pr) :: c, cd, cm,uni
    integer :: i97, j97


    c = 362436.0_k_pr / 16777216.0_k_pr
    cd = 7654321.0_k_pr / 16777216.0_k_pr
    cm = 16777213.0_k_pr /16777216.0_k_pr
    i97 = 97
    j97 = 33
    uni = u(i97) - u(j97)
    if( uni < 0.0_k_pr) uni = uni + 1.0_k_pr
    u(i97) = uni
    i97 = i97 - 1
    if(i97 == 0) i97 = 97
    j97 = j97 - 1
    if(j97 == 0) j97 = 97
    c = c - cd
    if( c < 0.0_k_pr ) c = c + cm
    uni = uni - c
    if( uni < 0.0_k_pr ) uni = uni + 1.0_k_pr
    ranmar = uni
  end function ranmar

!> \brief computes all the factorials up to n
!> \author Alin M Elena
!> \date 01/11/07, 14:19:01
!> \param n integer
!> \param fact real(k_pr) will contain the factorials on position 0 -> 0! ... n->n!

  subroutine InitFact(n,fact)
    character(len=*), parameter :: myname="fact"
    integer, intent(in) :: n
    real(k_pr), intent(inout) :: fact(0:n)
    integer :: i

    fact(0)=1.0_k_pr
    do i=1,n
      fact(i)=fact(i-1)*real(i,k_pr)
    enddo

  end subroutine InitFact
!> \brief implements delta function
!> \details \f[ \delta_{ab}= \begin{cases} 1 & a=b \\
!> 0 & a\neq b
!> \end{cases}
!> \f]
!> \author Alin M Elena
!> \date 01/11/07, 15:07:30
!> \param a,b  integers

  real(k_pr) function d(a,b)
    integer :: a,b
    d=0.0_k_pr
    if (a==b) then
      d=1.0_k_pr
    else
      d=0.0_k_pr
    endif
  end function d

!> \brief implements an discrete step function (Heaviside)
!> \details \f[ \theta(a)= \begin{cases} 1 & a > 0 \\
!> 0 & a\leq 0
!> \end{cases}
!> \f]
!> \author Alin M Elena
!> \date 01/11/07, 15:07:30
!> \param a  integers
!> \param eq optional logical its presence defines the behaviour in 0
!> if present and true the function is 1 in 0, if absent is 0
  real(k_pr) function h(a,eq)
    character(len=*),parameter :: mynale="h"
    integer,intent(in) :: a
    logical, intent(in),optional ::eq

    if (present(eq)) then
      if (a>=0) then
        h=1.0_k_pr
      else
        h=0.0_k_pr
      endif
    else
      if (a>0) then
        h=1.0_k_pr
      else
        h=0.0_k_pr
      endif
    endif
  end function h

!> \brief returns the atomic no Z
!> \author Alin M Elena
!> \date 30th of October 2007
!> \param atomic type(atomicxType) data about atoms
!> \param atom integer the id of an atom
  integer function GetZ(atomic,atom)
    character(len=*), parameter :: sMyName="GetZ"
    type(atomicxType), intent(in) :: atomic
    integer, intent(in) :: atom

    GetZ=atomic%species%z(GetSpecie(atomic,atom))
  end function GetZ

!> \brief returns the id of the specie of an atom
!> \author Alin M Elena
!> \date 30th of October 2007
!> \param atomic type(atomicxType) data about atoms
!> \param atom integer the id of an atom
  integer function GetSpecie(atomic,atom)
    character(len=*), parameter :: sMyName="GetSpecie"
    type(atomicxType), intent(in) :: atomic
    integer, intent(in) :: atom

    GetSpecie=atomic%atoms%sp(atom)
  end function GetSpecie

!> \brief computes the local magnetic moment from the excess density for an atom
!> \author Alin M Elena
!> \date 3rd November 2007
!> \param io type(ioType) contains all the info about I/O files
!> \param atomic type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param at integer the atom
!> \param show logical if true and verbosity high enough Prints the moemnt info
  real(k_pr) function LocalMoment(at,atomic,show,io,sol)
    character(len=*), parameter :: myname="LocalMoment"

    integer, intent(inout),optional :: at
    logical,intent(in)  :: show
    type(atomicxType),intent(inout) :: atomic
    type(ioType),intent(in) :: io
    type(solutionType), intent(in) :: sol
    integer :: from,to,k,j,m
    real(k_pr) :: m_down,m_up,aux,mu,md

    j=0
    m_down=0.0_k_pr
    m_up=0.0_k_pr
    mu=0.0_k_pr
    md=0.0_k_pr
    LocalMoment=0.0_k_pr
    aux=0.0_k_pr
    m=atomic%basis%norbitals*(atomic%basis%norbitals-1)/2
    do k=0,GetLmax(atomic%atoms%sp(at),atomic%speciesBasis,atomic%species)
    ! spin down
      from=atomic%atoms%orbs(at,1)+j
      to=from+2*k
      m_down=sum(sol%density(from+m:to+m))+sum(sol%n0(from:to))
      md=md+m_down
    !spin up
      from=atomic%atoms%orbs(at,1)+j+atomic%basis%norbitals/2
      to=from+2*k
      m_up=sum(sol%density(from+m:to+m))+sum(sol%n0(from:to))
      mu=mu+m_up
      aux=aux+(m_up-m_down)
      if (show .and. io%verbosity > k_highVerbos) then
        write(io%uout,"(a,i5,a,i5,a,f16.8,a,f16.8,a,f16.8,a)")&
            "local moment on atom ",at," for l= ",k,": ",m_up-m_down,"(",m_up,"-",m_down,")"
      endif
      j=j+2*k+1
    enddo
    LocalMoment=aux
    if (show .and. io%verbosity > k_highVerbos) then
      write(io%uout,"(a,i5,a,f16.8,a,f16.8,a,f16.8,a)")&
          "local moment on atom ",at,": ",LocalMoment,"(",mu,"-",md,")"
    endif
  end function LocalMoment

!> \brief returns the maximum l for an atom
!> \author Alin M Elena
!> \date 03/11/07, 11:03:47
!> \param at integer the atom
!> \param atomic type(atomicx) all the info about atoms
!> \param specBas type(orbitalType) all the info about species basis
!> \param spec type(speciesType) all the info about species
  integer function GetLmax(sp,specBasis,spec)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'GetLmax'
    !--subroutine parameters -------------------------!
    integer, intent(in) ::  sp
    type(orbitalType),intent(in) :: specBasis(:,:)
    type(speciesType), intent(in) :: spec
! for an atom at gives you the maximum l (quantum orbital number)
! in the basis
    !--internal variables ----------------------------!
    integer :: i,lmax
    lmax=0
!    do i=1,atomic%atoms%norbs(at)
! we suppose that the basis is the same for spin up and spin down
!      if (lmax<atomic%basis%orbitals(atomic%atoms%orbs(at,i))%l) &
!        lmax=atomic%basis%orbitals(atomic%atoms%orbs(at,i))%l
!    enddo

    do  i=1,spec%norbs(sp)
      if (lmax<specBasis(sp,i)%l) &
        lmax=specBasis(sp,i)%l
    enddo
    GetLmax=lmax
  end function GetLmax

!> \brief returns the maximum l for an atom
!> \author Alin M Elena
!> \date 25/11/08, 23:03:47
!> \param at integer the atom
!> \param specBas type(orbitalType) all the info about species basis
!> \param spec type(speciesType) all the info about species
  integer function Lmax(specBas,spec)
    !--subroutine name--------------------------------!
    character(len=*), parameter :: myname = 'Lmax'
    !--subroutine parameters -------------------------!
    type(orbitalType),intent(in) :: specBas(:,:)
    type(speciesType), intent(in) :: spec
! for an atom at gives you the maximum l (quantum orbital number)
! in the basis
    !--internal variables ----------------------------!
    integer :: i,lm,tmp
    lm=0
    do i=1,spec%nspecies
! we suppose that the basis is the same for spin up and spin down
      tmp=GetLMax(i,specBas,spec)
      if (lm<tmp) lm=tmp
    enddo
    Lmax=lm
  end function Lmax



!> \brief returns the distance between two atoms
!> \details also returns the cosine directions l,m,n
!> \author Alin M Elena
!> \date 07/11/07, 11:02:05
!> \param l,m,n reals output cosine directions
!> \param r real the distance between atoms
!> \param atoms type(atomicType) atoms data
!> \param j,i integers the ids of the atoms

  subroutine AtomDistance(atoms,j,i,r,l,m,n)
    character(len=*), parameter :: myname = 'AtomDistance'
    !--subroutine parameters--------------------------!
    integer, intent(in)  :: i,j
    type(atomicType),intent(in) :: atoms
    real(k_pr), intent(inout) :: r,l,m,n

    l = atoms%x(j) - atoms%x(i)
    m = atoms%y(j) - atoms%y(i)
    n = atoms%z(j) - atoms%z(i)
    r= sqrt(l*l+m*m+n*n)
    l=l/r
    m=m/r
    n=n/r
  end subroutine AtomDistance


!> \brief returns the distance between two atoms
!> \details also returns the componets of the vector bewtween them
!> \author Alin M Elena
!> \date 07/11/07, 16:50:05
!> \param rx,ry,rz reals output components of the r vector
!> \param r real the distance between atoms
!> \param atoms type(atomicType) atoms data
!> \param j,i integers the ids of the atoms

  subroutine AtomDistVec(atoms,j,i,r,rx,ry,rz)
    character(len=*), parameter :: myname = 'AtomDistVec'
    !--subroutine parameters--------------------------!
    integer, intent(in)  :: i,j
    type(atomicType),intent(in) :: atoms
    real(k_pr), intent(inout) :: r,rx,ry,rz

    rx = atoms%x(j) - atoms%x(i)
    ry = atoms%y(j) - atoms%y(i)
    rz = atoms%z(j) - atoms%z(i)
    r= sqrt(rx*rx+ry*ry+rz*rz)
  end subroutine AtomDistVec

!> \brief returns the distance between two atoms
!> \author Alin M Elena
!> \date 07/11/07, 17:02:05
!> \param atoms type(atomicType) atoms data
!> \param j,i integers the ids of the atoms
  real(k_pr) function Distance(atoms,j,i)
    character(len=*), parameter :: myname = 'Distance'
    !--subroutine parameters--------------------------!
    integer, intent(in)  :: i,j
    type(atomicType),intent(in) :: atoms

    Distance = sqrt((atoms%x(j) - atoms%x(i))*(atoms%x(j) - atoms%x(i))+&
                    (atoms%y(j) - atoms%y(i))*(atoms%y(j) - atoms%y(i))+&
                     (atoms%z(j) - atoms%z(i))*(atoms%z(j) - atoms%z(i)))
  end function Distance

!> \brief implements FermiDirac distribution
!> \author Alin M Elena
!> \date 07/11/07, 14:49:02
!> \param temp real the temperature
!> \param energy real the energy
!> \param mu real chemical potential
!>\remarks
!> \f[ f_{FD}=\cfrac{1}{e^{\cfrac{\epsilon - \mu}{k_BT}}+1} \f]
   function fermi(temp,energy,mu)
    !--subroutine name--------------------------------!
      character(len=*), parameter :: myname = 'fermi'
      real(k_pr) :: temp
      real(k_pr) :: energy
      real(k_pr) :: mu
      real(k_pr) :: fermi
    !--internal variables ----------------------------!
      real(k_pr) :: expo
    !-------------------------------------------------!
         expo = (energy-mu)/(temp*k_kb)
         if (expo<-100.0_k_pr) then
            fermi = 1.0_k_pr
         elseif (expo>100.0_k_pr) then
            fermi = 0.0_k_pr
         else
            fermi = 1.0_k_pr / (exp(expo)+1.0_k_pr)
         endif

   end function fermi

  real(k_pr) function occupMP(gen,sol,energy)
    character(len=*),parameter :: myname="occupMP"
    type(generalType), intent(inout) :: gen
    type(solutionType), intent(inout) :: sol
    real(k_pr), intent(in) :: energy
    real(k_pr) :: aux,expo,sum
    integer :: i

    aux=(energy-gen%electronicMU)/gen%mpW
    call InitializeHermite(aux,2*gen%mpN+1,sol%hermite)
    if (aux*aux<-100.0_k_pr) then
      expo= 1.0_k_pr
    elseif (aux*aux>100.0_k_pr) then
      expo = 0.0_k_pr
    else
      expo = exp(-aux*aux)
    endif

    sum=0.5_k_pr*derfc(aux)
    do i=1,gen%mpN
      sum=sum+(-1)**i/(sol%fact(i)*4.0_k_pr**i*sqrt(k_pi))*expo*sol%hermite(2*i-1)
    enddo
    occupMP=sum
  end function occupMP

  real(k_pr) function derfc(x)
    character(len=*), parameter :: myname = 'derfc'
    real(k_pr), intent(in) :: x
    real(k_pr) :: z, t
    z = abs(x)
    t = 1.0_k_pr / (1.0_k_pr + 0.5_k_pr * z)

    derfc = t*exp(-(z*z)-1.26551223_k_pr+t*(1.00002368_k_pr+t*(0.37409196_k_pr+&
      t*(0.09678418_k_pr+t*(-0.18628806_k_pr+ &
      t*(0.27886807_k_pr+t*(-1.13520398_k_pr+ &
      t*(1.48851587_k_pr+t*(-0.82215223_k_pr+t*.17087277_k_pr)))))))))

    if (x<0.0_k_pr) derfc=2.0_k_pr-derfc

   end function derfc

!> \brief tells you where the line i stops \see aidx function
!> \author Alin M Elena
!> \date 08/11/07, 12:21:21
!> \param i integer the line
!> \param n integer the dimension of the matrix
  integer function ka(i,n)
    character(len=*), parameter :: myname = 'ka'
    integer :: i,n

    ka=n*(i-1)-(i-1)*i/2
  end function ka

!> \brief gives you the position in an array of the i,j matrix element
!> \details having a symmetric square matrix of dimension n, we want to pack
!> in an array the upper triangular part by line followed by the diagonal part
!>  this function will give you the position in this array where the i,j element
!>  of the original matrix is stored
!> \author Alin M Elena
!> \date 08/11/07, 12:16:19
!> \param i,j integers coordinates of the matrix element
!> \param n integer the dimension of the square matrix
  integer function aidx(i,j,n)
    character(len=*), parameter :: myname = 'aidx'
    integer, intent(in) ::  i,j,n
    if (j>i) then
      aidx=ka(i,n)+j-i
    elseif (j<i) then
      aidx=ka(j,n)+i-j
    else
      aidx=n*(n-1)/2+i
    endif
  end function aidx



!> computes a norm to estimate errors
!> \param approx is a real(k_pr) array containing the energy
!> \param exact eal(k_pr) array containing the forces
!> \param dh real represents the step
!> \details for each point from lbound(approx,1)+1 to ubound(approx,1)-1
!> computes the numeric derivative of approx and return the maximum error according to
!> formula
!> \f[
!> error=\max(|approx_i-exact_i|_{i=1,n})
!> \f]
!> divided by the modulus of the exact force in that point
!> The result is in percents.
!> author Alin M. Elena (Belfast)
!> \date 25th of April 2006
  real(k_pr) function norm(approx,exact,dh)
    character(len=*), parameter ::   myname="norm"
    real(k_pr), intent(in) :: approx(:), exact(:),dh
    real(k_pr) :: sum,work
    integer :: i

    sum=0.0_k_pr
!compute linf error or minmax error:sum
!returns what percent sum represents from the exact value
    do i=2,size(approx)-1
      work=(approx(i+1)-approx(i-1))/(2.0_k_pr*dh)
      if (abs(exact(i)-work)>sum) then
        sum=abs(exact(i)-work)
        norm=sum/abs(exact(i))*100.0_k_pr
      endif
    enddo
  end function norm

!> \brief computes the modulus of vecotr
!> \author Alin M Elena
!> \date 10/11/07, 13:25:51
!> \param x,y,z reals cartesian coordinates of the vector
  real(k_pr) function VectorModulus(x,y,z)
    character(len=*), parameter :: myname = 'VectorModulus'
    real(k_pr) ::x,y,z
    VectorModulus=sqrt(x*x+y*y+z*z)
  end function VectorModulus
!> \brief returns the position in an array of the maximum
!> \author Alin M Elena
!> \date 12/11/07, 13:25:12
!> \param arr theh array
  function imaxloc(arr)
    real(k_pr), dimension(:), intent(in) :: arr
    integer :: imaxloc
    integer, dimension(1) :: imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  end function imaxloc

!> \brief returns the position in an array of the minimum
!> \author Alin M Elena
!> \date 12/11/07, 13:26:12
!> \param arr theh array
  function iminloc(arr)
    real(k_pr), dimension(:), intent(in) :: arr
    integer, dimension(1) :: imin
    integer :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc

!> \brief asserts if 3 integers are equal
!> \author Alin M Elena
!> \date 12/11/07, 13:23:47
!> \param n1,n2,n3 the integers
!> \param io type(ioType) i/o units

  integer function AssertEq3(n1,n2,n3,io)
    character(len=*), parameter :: myname="AssertEq3"
    integer, intent(in) :: n1,n2,n3
    type(ioType), intent(inout) :: io

    if (n1 == n2 .and. n2 == n3) then
      AssertEq3=n1
    else
     call  error("AssertEq failed:",myname,.true.,io)
   end if
  end function AssertEq3

!> \brief asserts if 4 integers are equal
!> \author Alin M Elena
!> \date 12/11/07, 13:23:47
!> \param n1,n2,n3,n4 the integers
!> \param io type(ioType) i/o units
  integer function AssertEq4(n1,n2,n3,n4,io)
    character(len=*), parameter :: myname="AssertEq4"
    integer, intent(in) :: n1,n2,n3,n4
    type(ioType), intent(inout) :: io
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
      AssertEq4=n1
    else
      call  error("AssertEq failed:",myname,.true.,io)
    end if
  end function AssertEq4

!> \brief swaps two scalars
!> \author Alin M Elena
!> \date 12/11/07, 13:23:20
!> \param a,b reals the scalars
  subroutine SwapScalar(a,b)
    real(k_pr), intent(inout) :: a,b
    real(k_pr) :: dum
    dum=a
    a=b
    b=dum
  end subroutine SwapScalar

!> \brief swaps two vectors
!> \author Alin M Elena
!> \date 12/11/07, 13:22:37
!> \param a,b the vectors
  subroutine SwapVector(a,b)
    real(k_pr), dimension(:), intent(inout) :: a,b
    real(k_pr), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine SwapVector
!> \brief retunr the exponential of a number
!> \details it takes care that the exponential has a resonable big or small value
!> \author Alin M Elena
!> \date 12/11/07, 15:14:16
!> \param rdum real the exponent of the exponential

  real(k_pr) function  ExpRep(rdum)
    real(k_pr),intent(in) ::  rdum
    if (rdum > 100._k_pr) then
      ExpRep = 3.69+50_k_pr
    else if (rdum < -100._k_pr) then
      ExpRep = 0.0_k_pr
    else
        ExpRep = exp(rdum)
    end if
  end function ExpRep

  subroutine ComputeEuclideanMatrix(atoms,io,euclidDistances)
    character(len=*), parameter :: myname="ComputeEuclideanMatrix"
    type(ioType), intent(in) :: io
    type(atomicType), intent(in) :: atoms
    real(k_pr), intent(inout),optional :: euclidDistances(:,:)
    real(k_pr),allocatable :: tmp(:,:)
    integer :: i,j

    if (.not.present(euclidDistances)) then
      allocate(tmp(1:atoms%natoms,1:atoms%natoms))
      tmp=0.0_k_pr
    endif
    if (.not.present(euclidDistances)) then
      do i=1,atoms%natoms-1
        do j=i+1,atoms%natoms
          tmp(i,j)=Distance(atoms,i,j)
          tmp(j,i)=tmp(i,j)
        enddo
      enddo
    else
      do i=1,atoms%natoms-1
        do j=i+1,atoms%natoms
          euclidDistances(i,j)=Distance(atoms,i,j)
          euclidDistances(j,i)=euclidDistances(i,j)
        enddo
      enddo
    endif
    if (.not.present(euclidDistances)) then
      write(io%uout,*) "==Euclidean Distances Matrix======"
      write(io%uout,'(7x)',advance="no")
      do i=1,atoms%natoms
         write(io%uout,'(i12,1x)',advance="no") i
      enddo
      write(io%uout,*)
      do i=1,atoms%natoms
         write(io%uout,'(i6,1x)',advance="no") i
         do j=1,atoms%natoms
            write(io%uout,'(f12.6,1x)',advance="no") tmp(i,j)
         enddo
         write(io%uout,'(1x,i0)')i
      enddo
      write(io%uout,'(7x)',advance="no")
      do i=1,atoms%natoms
         write(io%uout,'(i12,1x)',advance="no") i
      enddo
      write(io%uout,*)
      write(io%uout,*)"==================================="
      deallocate(tmp)
    endif
  end subroutine ComputeEuclideanMatrix

  subroutine InitializeHermite(x, n, h)
      character(len=*), parameter :: myname = 'InitializeHermite'
      real(k_pr), intent(inout)  :: h(0:n)
      real(k_pr), intent(in) :: x
      integer, intent(in) :: n

      integer :: i
      h(0)=1.0_k_pr
      h(1)=2.0_k_pr*x
      do i=2,n
         h(i)=2.0_k_pr*x*h(i-1)-2.0_k_pr*real(i-1,k_pr)*h(i-2)
      enddo
   end subroutine InitializeHermite

  real(k_pr) function sn(x,n,sol)
    character(len=*), parameter :: myname = 'sn'
    real(k_pr), intent(in) :: x
    integer, intent(in) :: n
    type(solutionType), intent(inout) :: sol
    real(k_pr) :: expo

    call InitializeHermite(x,2*n+1,sol%hermite)
    if (x*x<-100.0_k_pr) then
        expo= 1.0_k_pr
    elseif (x*x>100.0_k_pr) then
        expo = 0.0_k_pr
    else
        expo = exp(-x*x)
    endif
    sn=0.5_k_pr*expo*sol%hermite(2*n)*(-1)**n/(sol%fact(n)*4**n*sqrt(k_pi))

   end function sn

  real(k_pr) function MarzariF(x)
    character(len=*), parameter :: myname="MarzariF"
    real(k_pr), intent(in) :: x
    real(k_pr) :: aux,expo

    aux= (1.0_k_pr/sqrt(2.0_k_pr)-x)*(1.0_k_pr/sqrt(2.0_k_pr)-x)
    if (aux<-100.0_k_pr) then
        expo= 1.0_k_pr
    elseif (aux>100.0_k_pr) then
        expo = 0.0_k_pr
    else
        expo = exp(-aux)
    endif
    MarzariF=expo/sqrt(2.0_k_pr*k_pi)+0.5_k_pr*derfc(1.0_k_pr/sqrt(2.0_k_pr)-x)
  end function MarzariF

  real(k_pr) function MarzariS(x)
    character(len=*), parameter :: myname="MarzariS"
    real(k_pr), intent(in) :: x
    real(k_pr) :: aux,expo

    aux= (1.0_k_pr/sqrt(2.0_k_pr)-x)*(1.0_k_pr/sqrt(2.0_k_pr)-x)
    if (aux<-100.0_k_pr) then
        expo= 1.0_k_pr
    elseif (aux>100.0_k_pr) then
        expo = 0.0_k_pr
    else
        expo = exp(-aux)
    endif
    MarzariS=expo/sqrt(k_pi)+0.0_k_pr*(1.0_k_pr-sqrt(2.0_k_pr)*x)
  end function MarzariS

  character(len=50) function getUnits(gen)
    character(len=*),parameter :: myName="getUnits"
    type(generalType), intent(in) :: gen

    getUnits="Unknown system of units, please update "//myName
    select case(gen%units)
      case(k_unitsEV)
        getUnits="eV units (eV-Angstrom)"
      case(k_unitsAU)
        getUnits="atomic units (Hartree-Angstrom)"
      case(k_unitsSI)
        getUnits="SI units (J-m)"
      case(k_unitsARU)
        getUnits="atomic Rydberg units (Rydberg-Bohr)"
    end select
  end function getUnits
end module m_Useful
