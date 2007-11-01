!> \mainpage TDTB+UJ Manual
!> Time-Dependent Tight-Binding+UJ
!> \brief main program Time-Dependent Tight-Binding+UJ
!> \author Alin M Elena


program tbuj
  use constants
  use types
  use read_data
  use useful, only : dateandtime
  implicit none


  integer :: narguments
  character(len=mw) :: arg
  type(io_type) :: io_info
  type(general_type) :: general
  type(atomicx_type) :: atomicx
  type(model_type) :: tb_model
  character(len=10) :: dt
  character(len=12) :: tm

  call dateandtime(dt,tm)
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

  call initialize(io_info,general,atomicx,tb_model)

!!!!!! closes all the units associated with blocks and deallocated the trees for tokens and blocks
  

  call clean_memory(atomicx,general,tb_model)
  call cpu_time(general%time%end)
  write(io_info%uout,'(a,a,a,a)',advance="no")"Program has started at ",dt," ",tm
  call dateandtime(dt,tm)
  write(io_info%uout,'(a,a,a,a)')" ended at ",dt," ",tm
  write(io_info%udeb,'(a,f16.6,a)')"Program has run for "&
    ,general%time%end-general%time%start," seconds"
  write(io_info%uout,'(a,f0.6,a)')"Program has run for "&
    ,general%time%end-general%time%start," seconds"

  call close_io_general(io_info)

end program tbuj
