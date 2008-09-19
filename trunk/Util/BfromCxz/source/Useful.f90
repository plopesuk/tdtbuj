module m_Useful
  implicit none
  private
!
  integer, parameter, public :: k_pr = kind (1.0d0)!< define precision for reals
  integer, parameter, public :: k_mw = 40 !< number of character per word
  integer, parameter, public :: k_ml = 255 !< number of character per line
  real (k_pr), parameter, public :: k_pi = 3.14159265358979323846264338327950_k_pr !< \f$ \pi \f$
  real (k_pr), parameter, public :: k_mju0 = 1.6729405d-04 ! mju0 in Atomic Rydberg Units
  real (k_pr), parameter, public :: k_b2tesla = 3.3241346D+05 ! atomic magnetic field to tesla
!
!
!> data structure that contains the info about a frame
  type, private :: frameType
    character (len=k_mw) :: frameLabel
    real (k_pr) :: timestamp
    real (k_pr), allocatable, dimension (:) :: x, y, z, dx, dy, dz, q, curr1, curr2, bondCurrent
    character (len=2), allocatable, dimension (:) :: element
    real (k_pr) :: tdx, tdy, tdz, dipole, charge, projDipole, b (1:3)
  end type frameType
!
!> data structure that contains the info about frames and atoms
  type, public :: cxzType
    integer :: natoms = 0 !< number of atoms in each frame
    integer :: nframes = 0 !< number of frames in the file
    integer :: nring = 0 !< number of atoms in the ring
    integer, allocatable :: ring (:)!< the atoms that make the ring
    type (frameType), allocatable :: frames (:)!< contains the info about frames
  end type cxzType
!
  public :: error
  public :: GetUnit
  public :: isInList
  public :: norm
contains
!
!> \brief Prints an error/warning message and aborts the programs if necessary
!> \author Alin M. Elena (Queen's University Belfast)
!> \date 14-15th January 2006
!> \param message the message that will be displayed
!> \param routine the name of the caller
!> \param uout integer  unit number for output
!> \param critical logical if true aborts the program
  subroutine error (message, routine, critical, uout)
    character (len=*), intent (in) :: message
    character (len=*), intent (in) :: routine
    logical, intent (in) :: critical
    integer, intent (in) :: uout
!
    if (critical) then
      write (uout,*) "Critical error in subroutine: ", routine
    else
      write (uout,*) "Error message from subroutine: ", routine
    end if
!
    write (uout,*) routine, ": ", message
!
    if (critical) then
      write (uout,*) routine, ": User stop."
      write (*,*) routine, " User stop."
      stop
    end if
  end subroutine error
!
!> \brief gives an integer that can be used as a unit number in a open command
!> \details this function should be used to avoid using the same unit in different open commands
!> \author Alin M Elena
!> \date 14-15th of January, 2006
!> \warning using a save attribute makes the function thread unsafe, however
!> I do not expect any problems from it (you do not open so many files after all)
  integer function GetUnit ()
    character (len=*), parameter :: myname = 'GetUnit()'
    integer, save :: ustart = 9
!
    ustart = ustart + 1
    GetUnit = ustart
!
  end function GetUnit
!
!
!> \brief checks if the elements of an integer is in a list
!> \author Alin M Elena
!> \date 29/10/07, 23:50:47
!> \param x integer array
!> \param y integer number to be searched
!> \param pos integer optional returns the position where y is found or -1 if is not found
  logical function isInList (y, x, pos)
    character (len=*), parameter :: sMyName = "isInList"
    integer, intent (in) :: x (:)
    integer, intent (in) :: y
    integer, intent (inout), optional :: pos
!
    integer :: i
    integer :: n
    logical :: aux
!
    n = size (x)
    if (present(pos)) pos = - 1
    aux = .false.
    if (n >= 1) then
      do i = 1, n
        if (y == x(i)) then
          aux = .true.
          if (present(pos)) pos = i
        end if
      end do
    end if
    isInList = aux
  end function isInList
!
!
!> \brief computes the norm of a vector
!> \author Alin M Elena
!> \date 18/08/2008, 10:41
!> \param x real array
!
  real (k_pr) function norm (x)
    character (len=*), parameter :: myname = "norm"
    real (k_pr), intent (in) :: x (:)
    norm = Sqrt (dot_product(x, x))
  end function norm
!
end module m_Useful
