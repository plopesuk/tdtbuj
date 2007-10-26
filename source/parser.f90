!> deals with all subroutines related with the parsing
!> \author Alin M. Elena
!> \date 14th of January 2006
!> \remark 17th of July 2007 removed the save attribute from variables
module parser


  use constants
  use useful
  use types

  implicit none
  private
  public :: parse_file,end_parse
  public :: get_logical,get_string,get_real,get_integer,get_block


!> data structure of a node in the parsing tree
type, private :: names
  integer :: lines  !< the number of lines useful in the case of blocks
  character(len=mw) :: nam !< name of the entity read in file
  character(len=mw) :: value !< the value associated with the node
  integer  :: ut !< unit number where write the info to be parsed (write(ut,*))
  type(names),pointer :: next=>null() !< next element in tree
end type names

!> data structure of info about the parsed files
	type, private :: info_parse
	  integer :: comments=0 !< number of comments parsed
	  integer :: empty=0 !< number of empty lines
	  integer :: includ=0 !<  number of include statements
	  integer :: tokens=0 !< number of tokens
	  integer :: blocks=0 !< number of blocks
	  integer :: lines=0 !< number of lines
	  integer :: blocklines=0 !< number of lines in blocks
    end type info_parse

    
    
  type(names),pointer:: bnames,tnames,currentb,currentt !< blocks and tokens lists for the files plus the current ones

contains

!>  \brief opens the input and error files
!>  \details allocates the lists for blocks and tokens
!>  starts the parsing
!> \author Alin M. Elena (Belfast)              
!> \date 14th-15th of January 2006
!> \param io_loc is type(io_type) (see types::io_type)
!> \remark 20th of January 2007 added the parameter \em io_loc


  subroutine parse_file(io_loc)
    character(len=*), parameter :: myname = 'parse_file'
    integer :: errno
    type(io_type), intent(inout) :: io_loc	
	type(info_parse) :: ireport

!allocate the lists for name of the blocks and tokens
    allocate(bnames)
    allocate(tnames)
    call parse(io_loc%uinp,io_loc, ireport)

    write(io_loc%uerr,'(a)') "This is not an error"
    write(io_loc%uerr,'(a)')"1. Blocks labels"
    call print_name(bnames,io_loc)
    write(io_loc%uerr,'(a)')"2. Tokens labels"
    call print_name(tnames,io_loc)
    write(io_loc%uerr,'(a)')"3. General info"
    write(io_loc%uerr,'(a,i5)') "No of lines read",ireport%lines
    write(io_loc%uerr,'(a,i5)') "No of comment lines",ireport%comments
    write(io_loc%uerr,'(a,i5)') "No of empty lines",ireport%empty
    write(io_loc%uerr,'(a,i5)') "No of included files",ireport%includ
    write(io_loc%uerr,'(a,i5)',advance='no') "No of blocks",ireport%blocks
    write(io_loc%uerr,'(a,i5,a)') " containing ",ireport%blocklines," lines"
    write(io_loc%uerr,'(a,i5)') "No of tokens",ireport%tokens

    call finalize(io_loc)
  end subroutine parse_file

 
!> \brief    recursive subroutine which effectively parses the file
!> \details    associated with unit=nounit
!>          if an include statement is found the routine is called again to parse the new file
!> \author Alin M. Elena (Belfast)
!> \date 14th-15th of January 2006
!> \param nounit integer, represents the unit number of a file which was previous opened and now
!> is parsed 
!> \param  io_loc type(io_type) (see io_type)
!> \param  ireport type(info_parse),  keeps the info about the parsed files (see types::info_parse type)
!> \remark 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added parameter \em io_loc


  recursive subroutine parse(nounit,io_loc, ireport)
    character(len=*), parameter :: myname = 'parse'
    type(io_type), intent(inout) :: io_loc	
    integer,intent(in) :: nounit	
    type(info_parse), intent(inout) :: ireport

    character(len=ml) :: line,lineaux,storeline,lineaux2
    character(len=mw) :: filename,nam,blockname,endblockname,tokename
    integer :: errno,noinc,lineno,blines,nt
    logical :: ex,op,bool

    lineno=0
    inquire(unit=nounit,name=nam)

    do while (get_line(nounit,line))
      lineaux=adjustl(line)
      lineno=lineno+1
      if (lineaux(1:1)=="!" .or. lineaux(1:1)=="#".or.lineaux(1:2)=="//") then
        ireport%comments=ireport%comments+1
    ! just count the comments nothing else
      elseif (len_trim(line)==0) then
        ireport%empty=ireport%empty+1
    ! just count empty lines
      elseif (cstr(lineaux(1:7),"include")) then
        ! defines the behaviour for an include statement
        ireport%includ=ireport%includ+1
        read(lineaux(8:255),*,iostat=errno)filename

        if (errno/=0) call parse_err("invalid name after include",myname,nam,lineno,io_loc)
        inquire(file=trim(filename),exist=ex,opened=op)        
        if (.not.ex) call parse_err("can not open file after include: "//trim(filename),myname,nam,lineno,io_loc)
        if (op) call parse_err("file has been already opened (two indentical lines in input)",&
          myname,nam,lineno,io_loc)

        noinc=get_unit()
        open(unit=noinc,file=trim(filename),status="old",iostat=errno)
        call parse(noinc,io_loc,ireport)

      elseif (cstr(lineaux(1:5),"block")) then
! defines the behaviour of a block - endblock statement

        ireport%blocks=ireport%blocks+1
        read(lineaux(6:ml),*,iostat=errno)blockname
        if (errno/=0) call parse_err("invalid name after block statement",myname,nam,lineno,io_loc)
        bool= .true.
        blines=0
!  create a scratch file
        nt=get_unit()
        filename=trim(blockname)//".blk"
        open(unit=nt,file=trim(filename),status="replace",action="write")
        do while (bool)
          bool= get_line(nounit,line)
          lineno=lineno+1
          lineaux=adjustl(line)
! take some action if we have a endblock       
          if (cstr(lineaux(1:8),"endblock")) then
            ! check the name
            read(lineaux(9:ml),*,iostat=errno)endblockname
            if (errno/=0) call parse_err("invalid name after endblock statement",myname,nam,lineno,io_loc)
! closing the wrong block
            if (.not.cstr(trim(endblockname),trim(blockname))) &
              call parse_err("closing wrong block, expected endblock "//trim(blockname)//&
              " found endblock "//trim(endblockname),myname,nam,lineno,io_loc)
            exit
          endif
! parse the content of the block line by line and get rid of commented and empty lines           
          if (lineaux(1:1)=="!" .or. lineaux(1:1)=="#".or.lineaux(1:2)=="//") then
            ireport%comments=ireport%comments+1
          elseif (len_trim(line)==0) then
            ireport%empty=ireport%empty+1
          else
! count only not obviuos without value lines
! get rid of any inline comment
            storeline=parse_line(lineaux)
            write(nt,'(a)')trim(storeline)
            blines=blines+1
          endif     
        end do

! check if we have reached end-of-file
        if (.not.bool) call parse_err("missing endblock "//trim(blockname),myname,nam,-1,io_loc)
! check empty blocks and give an warning, empty blocks are ignored
        if (blines==0)  then
          call parse_war("detected empty(probably only comments and empty lines) block "&
            //trim(blockname),myname,nam,lineno,io_loc)
! empty files we delete them immediately    
          close(unit=nt,status="delete")
        endif

! check the uniquness of the block  
        ireport%blocklines=blines+ireport%blocklines   
        if (blines/=0) then
          close(nt)   
! check for the block in the existent list
          if (ireport%blocks==1) then
            currentb=>bnames 
            bnames%nam=trim(blockname)
            bnames%lines=blines
            bnames%ut=nt
            bnames%value=trim(filename)
          else
            if (find_name(trim(blockname),bnames)) &
              call parse_err("found block "//trim(blockname)//" duplicated",&
              myname,nam,lineno,io_loc)
            call add_name(trim(blockname),currentb,blines)
            currentb%ut=nt
            currentb%value=trim(filename)
          endif
        endif 
      elseif (cstr(lineaux(1:8),"endblock")) then
! endblock without block    
        call parse_err("endblock without block",myname,nam,lineno,io_loc)
      else
        ireport%tokens=ireport%tokens+1
        read(lineaux,*,iostat=errno)tokename
        lineaux2=parse_line(lineaux)
        read(lineaux2(len_trim(tokename)+1:ml),*,iostat=errno)storeline
        if (errno/=0) storeline=''

! check for the token in the existent list and add it if is new
        if (ireport%tokens==1) then
          currentt=>tnames 
          tnames%nam=trim(tokename)
        else
          if (find_name(trim(tokename),tnames)) &
            call parse_err("found token "//trim(tokename)//" duplicated",&
            myname,nam,lineno,io_loc)
          call add_name(trim(tokename),currentt)
        endif 
        currentt%value=trim(storeline)
      endif

    enddo



    ireport%lines=ireport%lines+lineno
  end subroutine parse


!> \brief logical function reads a line from a unit
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param uno integer, unit number from where to read the line
!> \param line character(len=*), the line that was read
!> \return .true. if successfull in reading the line, .false. otherwise
!> \remarks
              
  logical function get_line(uno,line)
    character(len=*), parameter :: myname = 'get_line'
    character(len=ml), intent(out) :: line 
    integer, intent(in) ::uno
    integer :: errno
    inquire(unit=uno,iostat=errno)
    get_line=.false.
    if (errno/=0) then
      write(*,*)"Unexpected error opening the input file(s)"
      stop
    endif
    read(uno,fmt='(a)',iostat=errno)line
    if (errno==0) get_line=.true.

  end function get_line

!> \brief adds a new node in the list
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param word character(len=*), value for field nam of the node
!> \param current pointer to the last node in the list
!> \param lines integer, optional value for field lines 
!> \remarks


  subroutine add_name(word,current,lines)
    character(len=*), parameter :: myname = 'add_name'
    character(len=*),intent(in) :: word 
    integer,intent(in),optional :: lines
    type(names), pointer,intent(inout) :: current
    type(names),pointer :: node

    allocate(node)
    node%nam=trim(word)
    if(present(lines))  node%lines=lines
    current%next=>node
    current=>node

  end subroutine add_name

!> \brief logical function finds a field in a list starting at root
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param word character(len=*), value of the field nam to be searched for
!> \param root type(names), pointer, starting point in search
!> \param loc tppe(names), pointer, optional returns the location where the info was found.
!> \remarks

  logical function find_name(word,root,loc)
    character(len=*), parameter :: myname = 'find_name'
    character(len=*),intent(in) :: word
    type(names), pointer,intent(in) :: root
    type(names),pointer,intent(out),optional :: loc
    type(names), pointer :: current

    current=>root
    find_name=.false.
    do while(associated(current))
      if (cstr(trim(word),trim(current%nam))) then
        find_name=.true.
        if (present(loc)) loc=>current
        exit
      endif
      current=>current%next

    end do

  end function find_name


!> \brief prints at a specified unit the nam field  from a list
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param root type(names), pointer  starting node in the printing list 
!> \param io_loc type(io_type) conatins the unit to print (see types::io_type)
!> \remarks


  subroutine print_name(root,io_loc)
    character(len=*), parameter :: myname = 'print_name'
    type(names), pointer,intent(in) :: root
    type(io_type), intent(inout) :: io_loc		
    type(names),pointer :: current

    current=>root
    do while(associated(current))
      write(io_loc%uerr,'(a)') trim(current%nam)
      current=>current%next
    end do

  end subroutine print_name

!> \brief recursive function that deallocates all the nodes starting with root
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param root type(names), pointer the starting node
!> \remarks

  recursive subroutine delete_list(root)
    character(len=*), parameter :: myname = 'delete_list'
    type(names), pointer :: root
    type(names),pointer :: current

    current=>root%next

    if (associated(current)) then
      call delete_list(current)
    else
      deallocate(root)
    endif 

  end subroutine delete_list

!> \brief parses a line removing comments
!> \details everything after ! \# or // is considered a comment
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param line character(len=*) the line to be parsed
!> \remarks

  character(ml) function parse_line(line)
    character(len=*), parameter :: myname = 'parse_line' 
    character(len=*), intent(in) :: line
    integer ::p(1:3),pos


    p(1)=scan(line,"!")
    p(2)=scan(line,"#")
    p(3)=index(line,"//")
    where(p==0) p=ml+1
    pos=minval(p)-1
    if (pos/=ml) then
      parse_line=line(1:pos)
    else
      parse_line=line
    endif

  end function parse_line

!> \brief closes all the units that where opened for reading or writing during the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param io_loc type(io_type) contain the staring unit for the list of temporary units (see types::io_type)
!> \remarks

        
  subroutine finalize(io_loc)
    character(len=*), parameter :: myname = 'finalize'
    type(io_type), intent(inout) :: io_loc	
    integer :: i
    logical :: op
    character(len=mw) :: filename


    do i=io_loc%uerr+1,get_unit()-1
      inquire(unit=i,opened=op,name=filename)
      if (op) close(unit=i)
    enddo

  end subroutine finalize

!> \brief prints an error message and aborts the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param message the error message
!> \param routine the subprogram that generated the error
!> \param filename raeding this \em filename the error occured
!> \param lineno the line number that generated the error
!> \param io_loc type(io_type) (see types::io_type)
!> \remarks

          
  subroutine parse_err(message,routine,filename,lineno,io_loc)
    character(len=*), parameter :: myname = 'parse_err'
    character(len=*),intent(in) :: message,routine,filename
    integer, intent(in) :: lineno
    type(io_type), intent(inout) :: io_loc

    write(io_loc%uerr,'(a,a)')"Error: ",trim(message)
    write(io_loc%uerr,'(a,a)')"Routine: ",trim(routine)
    write(io_loc%uerr,'(a,a)')"File: ", trim(filename)
    write(io_loc%uerr,'(a,i7)')"Line number: ",lineno
    write(io_loc%uerr,'(a)')"User stop"
    write(*,'(a,a)')"User stop, but probably this is not what you want, check file: ",io_loc%inp_err
    stop

  end subroutine parse_err

!> \brief prints a warning message occured in the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param message the error message
!> \param routine the subprogram that generated the warning
!> \param filename raeding this \em filename the warning occured
!> \param lineno the line number that generated the warning
!> \param io_loc type(io_type) (see types::io_type)
!> \remarks
  subroutine parse_war(message,routine,filename,lineno,io_loc)
    character(len=*), parameter :: myname = 'parse_war'
    character(len=*),intent(in) :: message,routine,filename
    integer, intent(in) :: lineno
    type(io_type), intent(inout) :: io_loc	

    write(io_loc%uerr,'(a,a)')"Warning: ",trim(message)
    write(io_loc%uerr,'(a,a)')"Routine: ",trim(routine)
    write(io_loc%uerr,'(a,a)')"File: ", trim(filename)
    write(io_loc%uerr,'(a,i7)')"Line number: ",lineno
    write(io_loc%uerr,'(a)')"User warning"
    write(*,'(a,a)')"User warning, check file: ",trim(io_loc%inp_err)
  end subroutine parse_war


!> \brief  returns a logical value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to .true. if no dflt parameter is present
!> output is put in the units indicated by \em io_loc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param io_loc type(io_type) (see types::io_type)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em io_loc parameter


  logical function get_logical(io_loc,label,dflt)
    character(len=*), parameter :: myname = 'get_logical'
    character(len=*),intent(in) :: label
    type(io_type),intent(inout) :: io_loc
    logical, intent(in),optional :: dflt
    type(names),pointer :: found
!default value is set to .true. if no dflt parameter is present


    if (find_name(trim(label),tnames,found)) then
      if ((cstr(trim(found%value),"yes")).or.(cstr(trim(found%value),"true")) &
        .or.(cstr(trim(found%value),".true.")).or.(cstr(trim(found%value),"t"))&
        .or.(cstr(trim(found%value),"y")).or.(cstr(trim(found%value),"1"))) then 
        get_logical=.true.
      elseif ((cstr(trim(found%value),"no")).or.(cstr(trim(found%value),"false")) &
        .or.(cstr(trim(found%value),".false.")).or.(cstr(trim(found%value),"f")) &
        .or.(cstr(trim(found%value),"n")).or.(cstr(trim(found%value),"0"))) then 
        get_logical=.false.
      else
        call error("Unsuported value: "//&
          trim(label)//" "//trim(found%value),myname,.true.,io_loc)
      endif
      write(io_loc%udeb,'(a,2x,l)') trim(label),get_logical
    else
      if (present(dflt)) then
        get_logical=dflt
      else
        get_logical=.true.
      endif 
      write(io_loc%udeb,'(a,2x,l,a)') trim(label),get_logical, " default" 
    endif

  end function get_logical

!> \brief  returns a string value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to empty string if no dflt parameter is present
!> output is put in the units indicated by \em io_loc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param io_loc type(io_type) (see types::io_type)
!> \param print optional logical determines if the token is printed or not. if the parameter is missing is assumed .false.
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em io_loc parameter
  
character(len=mw) function get_string(io_loc,label,dflt,print)
    character(len=*), parameter :: myname = 'get_string'
    type(io_type),intent(inout) :: io_loc
    character(len=*),intent(in) :: label
    character(len=*), intent(in),optional :: dflt
    logical, intent(in), optional :: print

    type(names),pointer :: found
! default value if no dflt parameter is present is empty string
! so be sure that you provide a default parameter
    if (find_name(trim(label),tnames,found)) then
      get_string=trim(found%value)
      if (present(print).and.(.not.print)) then 
      else 
        write(io_loc%udeb,'(a,2x,a)') trim(label),trim(get_string)
      endif   
    else
      if (present(dflt)) then
        get_string=trim(dflt)
      else
        get_string=''
      endif
      if (present(print).and.(.not.print)) then 
      else   
        write(io_loc%udeb,'(a,2x,a,a)') trim(label),trim(get_string), " default" 
      endif   
    endif  
  end function get_string


!> \brief  returns a real(kind=pr) value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to 0.0_pr if no dflt parameter is present
!> output is put in the units indicated by \em io_loc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param io_loc type(io_type) (see types::io_type)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em io_loc parameter


  real(pr) function get_real(io_loc,label,dflt)
    character(len=*), parameter :: myname = 'get_real'
    type(io_type),intent(inout) :: io_loc
    character(len=*),intent(in) :: label
    real(pr), intent(in),optional :: dflt
    type(names),pointer :: found
    character(len=mw) :: aux
    integer :: errno
! default value if no dflt parameter is present is 0.0_pr

    if (find_name(trim(label),tnames,found)) then
      aux=trim(found%value)
      read(aux,*,iostat=errno)get_real
      if (errno/=0) call error("wrong value "//trim(found%value)&
        //" suplied for label "//trim(label),myname,.true.,io_loc)
      write(io_loc%udeb,'(a,2x,f32.16)') trim(label),get_real
    else
      if (present(dflt)) then
        get_real=dflt
      else
        get_real=0.0_pr
      endif
      write(io_loc%udeb,'(a,2x,f32.16,a)') trim(label),get_real, " default" 
    endif

  end function get_real


!> \brief  returns an integer value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to 0 if no dflt parameter is present
!> output is put in the units indicated by \em io_loc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param io_loc type(io_type) (see types::io_type)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em io_loc parameter


  integer function get_integer(io_loc,label,dflt)
    character(len=*), parameter :: myname = 'get_integer'
    type(io_type),intent(inout) :: io_loc
    character(len=*),intent(in) :: label
    integer, intent(in),optional :: dflt
    type(names),pointer :: found
    character(len=mw) :: aux
    integer :: errno
! default value if no dflt parameter is present is 0

    if (find_name(trim(label),tnames,found)) then
      aux=trim(found%value)
      read(aux,*,iostat=errno)get_integer
      if (errno/=0) call error("wrong value "//trim(found%value)&
        //" suplied for label "//trim(label),myname,.true.,io_loc)
      write(io_loc%udeb,'(a,2x,i0)') trim(label),get_integer
    else
      if (present(dflt)) then
        get_integer=dflt
      else
        get_integer=0
      endif
      write(io_loc%udeb,'(a,2x,i0,a)') trim(label),get_integer, " default" 
    endif

  end function get_integer


!> \brief  returns a logical value to indicate if it found in the input file(s) a valid block associated with the token \em label
!> \details  if no token is found the value returned is .false.
!> a unit from where the content of the block cand be read is returned via nt output is put in the units indicated by io_loc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param nt unit number from where to read the block data
!> \param io_loc type(io_type) (see types::io_type)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em io_loc parameter


  logical function get_block(io_loc,label,nt)
    character(len=*), parameter :: myname = 'get_block'
    type(io_type),intent(inout) :: io_loc	
    character(len=*),intent(in) :: label
    integer, intent(out) :: nt
    type(names),pointer :: found
    character(len=mw) :: line


    if (find_name(trim(label),bnames,found)) then
      get_block=.true.
      nt=found%ut
      line=found%value  
      open(unit=nt,action="read",status="old",file=trim(line))
      write(io_loc%udeb,'(a,2x,l)') trim(label),get_block
    else
      call error("Block "//trim(label)//"was not found",myname,.true.,io_loc)
    endif

  end function get_block

!> \brief ends the parsing process and cuts the \em tree (deallocates the memory used during the parsing process)
!> \author Alin M Elena
!> \date 14th of January 2006
!> \warning after the call to this function none of the get_* function will work
!> all the input has to be done before the call to it.

  subroutine end_parse
    character(len=*), parameter :: myname = 'end_parse'

! to be called only when there is nothing to be read
    type(names),pointer :: current
    integer :: nt
    character(len=mw) :: filename

    current=>bnames
    do while(associated(current))
      nt = current%ut
      filename=current%value
      close(nt)
      open(nt,file=trim(filename),status="old",action="read")
      close(nt,status="delete")
      current=>current%next
    end do

    call delete_list(tnames)
    call delete_list(bnames)
  end subroutine end_parse

end module parser
