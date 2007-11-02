!> \brief subprograms of general use
!> \author Alin M. Elena (Queen's University Belfast) 
!> \date 14-15th of January, 2006
module useful
  use constants
  use types

  implicit none
  private
  public :: cstr
  public :: error
  public :: get_unit
  public :: isUnique
  public :: isInList
  public :: dateandtime
  public :: ccvar
  public :: ccnlm
  public :: rmarin
  public :: ranmar
  public :: init_fact
  public :: h
  public :: d

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

!> \brief prints an error/warning message and aborts the programs if necessary
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th January 2006  
!> \param message the message that will be displayed
!> \param routine the name of the caller
!> \param io_loc I/O details
!> \param critical logical if true aborts the program
!> \remarks  20th of January, 2007, by Alin M Elena (Queen's University Belfast), added io_loc parameter

! provides an elegant way of handling different errors
  subroutine error(message, routine, critical,io_loc)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: routine
    logical, intent(in)   :: critical
    type(io_type), intent(in) :: io_loc
    
    if (critical) then
      write(io_loc%uout,*) &
        "Critical error in subroutine: ", routine
    else
      write(io_loc%uout,*) &
        "Error message from subroutine: ", routine
    end if
    
    write(io_loc%uout,*) routine,": ", message
    
    if (critical) then
      write(io_loc%uout,*)routine,": User stop."
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
  integer function get_unit()
    character(len=*), parameter :: myname = 'get_unit()'
    integer,save ::ustart=9

    ustart=ustart+1
    get_unit=ustart

  end function get_unit

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
  subroutine dateandtime(date,time)
    character(len=*), parameter :: sMyName="dateandtime"
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
  end subroutine dateandtime

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
    character(len=ml)   :: dump

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

!> \brief initializes the random number generator
!> \details returns an array of "seeds"
!> \author netlib
!> \param ij,kl starting seeds
!> \param u real array of "seeds
!> \param io type(io_type) i/o units
!> \remarks
!> here is a sequence of code that will initialize the generator
!> \code
!> real(pr) :: cr,seed(97)
!>  call cpu_time(cr)
!>  call rmarin(int((cr-int(cr))*3132),general%ranseed*30081),seed)
!> \endcode
  subroutine rmarin(ij,kl,u,io)
    character(len=*),parameter :: myname="rmarin"
    real(pr),intent(out) ::  u(:)
    integer, intent(in) :: ij,kl
    type(io_type),intent(in) :: io

    integer :: i,j,k,l,ii,jj,m
    real(pr) :: s,t

    if( ij < 0  .or.  ij > 31328  .or. &
        kl < 0  .or.  kl > 30081 ) then
      call error("first seed must be between 0 and 31328 the second seed must be between 0 and  30081",&
                 myname,.true.,io)
    endif
    i = mod(ij/177, 177) + 2
    j = mod(ij    , 177) + 2
    k = mod(kl/169, 178) + 1
    l = mod(kl,     169)
    do ii = 1, 97
        s = 0.0_pr
        t = 0.5_pr
        do jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) >= 32) then
              s = s + t
            endif
            t = 0.5_pr * t
        enddo
        u(ii) = s
    enddo

  end subroutine rmarin

!> \brief generates a random number
!> \author netlib
!> \param u(:) real array of "seeds"
  real(pr) function ranmar(u)
    real(pr), intent(inout) :: u(:)
    real(pr) :: c, cd, cm,uni
    integer :: i97, j97


    c = 362436.0_pr / 16777216.0_pr
    cd = 7654321.0_pr / 16777216.0_pr
    cm = 16777213.0_pr /16777216.0_pr
    i97 = 97
    j97 = 33
    uni = u(i97) - u(j97)
    if( uni < 0.0_pr) uni = uni + 1.0_pr
    u(i97) = uni
    i97 = i97 - 1
    if(i97 == 0) i97 = 97
    j97 = j97 - 1
    if(j97 == 0) j97 = 97
    c = c - cd
    if( c < 0.0_pr ) c = c + cm
    uni = uni - c
    if( uni < 0.0_pr ) uni = uni + 1.0_pr
    ranmar = uni
  end function ranmar

!> \brief computes all the factorials up to n
!> \author Alin M Elena
!> \date 01/11/07, 14:19:01
!> \param n integer 
!> \param fact real(pr) will contain the factorials on position 0 -> 0! ... n->n!

  subroutine init_fact(n,fact)
    character(len=*), parameter :: myname="fact"
    integer, intent(in) :: n
    real(pr), intent(inout) :: fact(0:n)
    integer :: i

    fact(0)=1.0_pr
    do i=1,n
      fact(i)=fact(i-1)*real(i,pr)
    enddo

  end subroutine init_fact
!> \brief implements delta function
!> \details \f[ \delta_{ab}= \begin{cases} 1 & a=b \\
!> 0 & a\neq b
!> \end{cases}
!> \f]
!> \author Alin M Elena
!> \date 01/11/07, 15:07:30
!> \param a,b  integers

  real(pr) function d(a,b)
    integer :: a,b
    d=0.0_pr
    if (a==b) then
      d=1.0_pr
    else
      d=0.0_pr
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
  real(pr) function h(a,eq)
    character(len=*),parameter :: mynale="h"
    integer,intent(in) :: a
    logical, intent(in),optional ::eq

    if (present(eq)) then
      if (a>=0) then
        h=1.0_pr
      else
        h=0.0_pr
      endif
    else
      if (a>0) then
        h=1.0_pr
      else
        h=0.0_pr
      endif
    endif
  end function h


end module useful

