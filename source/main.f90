!> \mainpage TDTB+UJ Manual
!> Time-Dependent Tight-Binding+UJ
!> \brief main program Time-Dependent Tight-Binding+UJ
!> \author Alin M Elena
!> \date 14th-15th January 2006

program tbuj
  use constants
  use types
  use read_data
  use parser, only : end_parse
  implicit none


  integer :: narguments
  character(len=mw) :: arg
  type(io_type) :: io_info
  type(general_type) :: general
  type(atomicx_type) :: atomicx

  narguments=iargc()
  call cpu_time(general%time%start)

! it reads the name of the input as inline argument
! if there is none the default name is inp

  if (narguments==1) then
    call getarg(1,arg)
    io_info%inp_file=arg
  else
    io_info%inp_file="inp"
  endif

  call initialize(io_info,general,atomicx)

!!!!!! closes all the units associated with blocks and deallocated the trees for tokens and blocks
  call end_parse

  call cpu_time(general%time%end)
  write(io_info%udeb,'(a,f16.6,a)')"Program has run for "&
    ,general%time%end-general%time%start," seconds"
  write(io_info%uout,'(a,f16.6,a)')"Program has run for "&
    ,general%time%end-general%time%start," seconds"
  call close_io_general(io_info)

end program tbuj
