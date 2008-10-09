module m_Gutenberg
  use m_Useful
  implicit none
  private
!
  public :: PrintCurrents
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
    write (unitout, '(a,i0,1x,a1,2a2,a,3f12.8)') "#Bond current between atom ", cxz%atom1, "(", cxz%element1, ") ", "position: ", &
   & cxz%x1, cxz%y1, cxz%z1
    write (unitout, '(a,i0,1x,a1,2a2,a,3f12.8)') "#and atom ", cxz%atom2, "(", cxz%element2, ") ", "position: ", cxz%x2, cxz%y2, &
   & cxz%z2
    write (unitout, '(a1,2a16)') "#", "Time ", "BondCurrent"
    do i = 1, cxz%nframes
      write (unitout, '(2g)') cxz%frames(i)%timestamp, cxz%frames(i)%bondCurrent
    end do
  end subroutine PrintCurrents
!
!
end module m_Gutenberg
