module m_Gutenberg
  use m_Useful
  implicit none
  private
!
  public :: PrintCurrents
  public :: PrintMagneticField
!
contains
!
  subroutine PrintCurrents (cxz, unitout)
    character (len=*), parameter :: myname = "PrintCurrents"
    type (cxzType), intent (inout) :: cxz
    integer, intent (inout) :: unitout
!
    character (len=k_ml) :: fmtContainer, saux
    integer :: i
!
    write (unitout, '(a1,a16)', advance="no") "#", "Time "
    do i = 1, cxz%nring
      write (unitout, '(a4,i8,a1,a2,a1)', advance="no") "Curr", i + 1, "-", trim (cxz%frames(1)%element(i)), ";"
    end do
    do i = 1, cxz%nring
      write (unitout, '(a4,i8,a1,a2,a1)', advance="no") "GCur", i + 1 + cxz%nring, "-", trim (cxz%frames(1)%element(i)), ";"
    end do
    write (unitout,*)
    write (unitout,*) "#to get the atom number for Current component substract 1 from last number after Curr"
    write (unitout, '(a,i0,a)') "#to get the atom number for Gamma Current component substract ", 1 + cxz%nring, " from last number&
   & after GCur"
    write (fmtContainer, '("(f16.8,1x,",i0,"g,",i0,"g)")') cxz%natoms, cxz%natoms
    do i = 1, cxz%nframes
      write (unitout, trim(fmtContainer)) cxz%frames(i)%timestamp, cxz%frames(i)%curr1, cxz%frames(i)%curr2
    end do
  end subroutine PrintCurrents
!
!
  subroutine PrintMagneticField (cxz, unitout)
    character (len=*), parameter :: myname = "PrintMagneticField"
    type (cxzType), intent (inout) :: cxz
    integer, intent (inout) :: unitout
    integer :: i
!
    do i = 1, cxz%nframes
      write (unitout, "(5g)") cxz%frames(i)%timestamp, cxz%frames(i)%b, norm (cxz%frames(i)%b)
    end do
!
  end subroutine PrintMagneticField
!
!
end module m_Gutenberg
