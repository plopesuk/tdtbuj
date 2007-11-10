!> \brief deals with all subroutines related with the parsing
!> \author Alin M. Elena
!> \date 14th of January 2006
!> \remark 17th of July 2007 removed the save attribute from variables
!> \todo suspected memory leak
module m_Parser
  use m_Constants
  use m_Useful
  use m_Types
  implicit none
  private
  public :: ParseFile,EndParse
  public :: GetLogical,GetString,GetReal,GetInteger,GetBlock

!> data structure of a node in the parsing tree
  type, private :: names
    integer :: lines=0  !< the number of lines m_Useful in the case of blocks
    character(len=k_mw) :: nam="" !< name of the entity read in file
    character(len=k_mw) :: value="" !< the value associated with the node
    integer  :: ut !< unit number where write the info to be parsed (write(ut,*))
    type(names),pointer :: next=>null() !< next element in tree
  end type names

!> data structure of info about the parsed files
  type, private :: infoParseType
    integer :: comments=0 !< number of comments parsed
    integer :: empty=0 !< number of empty lines
    integer :: includ=0 !<  number of include statements
    integer :: tokens=0 !< number of tokens
    integer :: blocks=0 !< number of blocks
    integer :: lines=0 !< number of lines
    integer :: blocklines=0 !< number of lines in blocks
  end type infoParseType

  type(names),pointer:: bnames,tnames,currentb,currentt !< blocks and tokens lists for the files plus the current ones

contains
!>  \brief opens the input and error files
!>  \details allocates the lists for blocks and tokens
!>  starts the parsing
!> \author Alin M. Elena (Belfast)              
!> \date 14th-15th of January 2006
!> \param ioLoc is type(ioType) (see m_Types::ioType)
!> \remark 20th of January 2007 added the parameter \em ioLoc
  subroutine ParseFile(ioLoc)
    character(len=*), parameter :: myname = 'ParseFile'
    integer :: errno
    type(ioType), intent(inout) :: ioLoc
    type(infoParseType) :: ireport

!allocate the lists for name of the blocks and tokens
    allocate(bnames)
    allocate(tnames)
    call parse(ioLoc%uinp,ioLoc, ireport)

    write(ioLoc%uerr,'(a)') "This is not an error"
    write(ioLoc%uerr,'(a)')"1. Blocks labels"
    call PrintName(bnames,ioLoc)
    write(ioLoc%uerr,'(a)')"2. Tokens labels"
    call PrintName(tnames,ioLoc)
    write(ioLoc%uerr,'(a)')"3. General info"
    write(ioLoc%uerr,'(a,i5)') "No of lines read",ireport%lines
    write(ioLoc%uerr,'(a,i5)') "No of comment lines",ireport%comments
    write(ioLoc%uerr,'(a,i5)') "No of empty lines",ireport%empty
    write(ioLoc%uerr,'(a,i5)') "No of included files",ireport%includ
    write(ioLoc%uerr,'(a,i5)',advance='no') "No of blocks",ireport%blocks
    write(ioLoc%uerr,'(a,i5,a)') " containing ",ireport%blocklines," lines"
    write(ioLoc%uerr,'(a,i5)') "No of tokens",ireport%tokens

    call finalize(ioLoc)
  end subroutine ParseFile

!> \brief    recursive subroutine which effectively parses the file
!> \details    associated with unit=nounit
!>          if an include statement is found the routine is called again to parse the new file
!> \author Alin M. Elena (Belfast)
!> \date 14th-15th of January 2006
!> \param nounit integer, represents the unit number of a file which was previous opened and now
!> is parsed 
!> \param  ioLoc type(ioType) (see ioType)
!> \param  ireport type(infoParseType),  keeps the info about the parsed files (see m_Types::infoParseType type)
!> \remark 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added parameter \em ioLoc

  recursive subroutine parse(nounit,ioLoc, ireport)
    character(len=*), parameter :: myname = 'parse'
    type(ioType), intent(inout) :: ioLoc
    integer,intent(in) :: nounit
    type(infoParseType), intent(inout) :: ireport

    character(len=k_ml) :: line,lineaux,storeline,lineaux2,filename,nam
    character(len=k_mw) :: blockname,endblockname,tokename
    integer :: errno,noinc,lineno,blines,nt
    logical :: ex,op,bool

    lineno=0
    inquire(unit=nounit,name=nam)
    do while (GetLine(nounit,line))
      lineaux=adjustl(line)
      lineno=lineno+1
      if (lineaux(1:1)=="!" .or. lineaux(1:1)=="#".or.lineaux(1:2)=="//") then
        ireport%comments=ireport%comments+1
    ! just count the comments nothing else
      elseif (len(trim(line))==0) then
        ireport%empty=ireport%empty+1
    ! just count empty lines
      elseif (cstr(lineaux(1:7),"include")) then
        ! defines the behaviour for an include statement
        ireport%includ=ireport%includ+1
        read(lineaux(8:255),*,iostat=errno)filename

        if (errno/=0) call ParseErr("invalid name after include",myname,nam,lineno,ioLoc)
        inquire(file=trim(filename),exist=ex,opened=op)        
        if (.not.ex) call ParseErr("can not open file after include: "//trim(filename),myname,nam,lineno,ioLoc)
        if (op) call ParseErr("file has been already opened (two indentical lines in input)",&
          myname,nam,lineno,ioLoc)

        noinc=GetUnit()
        open(unit=noinc,file=trim(filename),status="old",iostat=errno)
        call parse(noinc,ioLoc,ireport)

      elseif (cstr(lineaux(1:5),"block")) then
! defines the behaviour of a block - endblock statement
        ireport%blocks=ireport%blocks+1
        read(lineaux(6:k_ml),*,iostat=errno)blockname
        if (errno/=0) call ParseErr("invalid name after block statement",myname,nam,lineno,ioLoc)
        bool= .true.
        blines=0
!  create a scratch file
        nt=GetUnit()
        filename=trim(blockname)//".blk"
        open(unit=nt,file=trim(filename),status="replace",action="write")
        do while (bool)
          bool= GetLine(nounit,line)
          lineno=lineno+1
          lineaux=adjustl(line)
! take some action if we have a endblock       
          if (cstr(lineaux(1:8),"endblock")) then
            ! check the name
            read(lineaux(9:k_ml),*,iostat=errno)endblockname
            if (errno/=0) call ParseErr("invalid name after endblock statement",myname,nam,lineno,ioLoc)
! closing the wrong block
            if (.not.cstr(trim(endblockname),trim(blockname))) &
              call ParseErr("closing wrong block, expected endblock "//trim(blockname)//&
              " found endblock "//trim(endblockname),myname,nam,lineno,ioLoc)
            exit
          endif
! parse the content of the block line by line and get rid of commented and empty lines           
          if (lineaux(1:1)=="!" .or. lineaux(1:1)=="#".or.lineaux(1:2)=="//") then
            ireport%comments=ireport%comments+1
          elseif (len(trim(line))==0) then
            ireport%empty=ireport%empty+1
          else
! count only not obviuos without value lines
! get rid of any inline comment
            storeline=ParseLine(lineaux)
            write(nt,'(a)')trim(storeline)
            blines=blines+1
          endif
        end do
! check if we have reached end-of-file
        if (.not.bool) call ParseErr("missing endblock "//trim(blockname),myname,nam,-1,ioLoc)
! check empty blocks and give an warning, empty blocks are ignored
        if (blines==0)  then
          call ParseWar("detected empty(probably only comments and empty lines) block "&
            //trim(blockname),myname,nam,lineno,ioLoc)
! empty files we delete them immediately
        endif
! check the uniquness of the block
        ireport%blocklines=blines+ireport%blocklines

        close(nt)
! check for the block in the existent list
        if (ireport%blocks==1) then
          currentb=>bnames
          bnames%nam=trim(blockname)
          bnames%lines=blines
          bnames%ut=nt
          bnames%value=trim(filename)
        else
          if (FindName(trim(blockname),bnames)) &
            call ParseErr("found block "//trim(blockname)//" duplicated",&
              myname,nam,lineno,ioLoc)
          call AddName(trim(blockname),currentb,blines)
          currentb%ut=nt
          currentb%value=trim(filename)
        endif
      elseif (cstr(lineaux(1:8),"endblock")) then
! endblock without block    
        call ParseErr("endblock without block",myname,nam,lineno,ioLoc)
      else
        ireport%tokens=ireport%tokens+1
        read(lineaux,*,iostat=errno)tokename
        lineaux2=ParseLine(lineaux)
        read(lineaux2(len(trim(tokename))+1:k_ml),*,iostat=errno)storeline
        if (errno/=0) storeline=''
! check for the token in the existent list and add it if is new
        if (ireport%tokens==1) then
          currentt=>tnames 
          tnames%nam=trim(tokename)
        else
          if (FindName(trim(tokename),tnames)) &
            call ParseErr("found token "//trim(tokename)//" duplicated",&
            myname,nam,lineno,ioLoc)
          call AddName(trim(tokename),currentt)
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

  logical function GetLine(uno,line)
    character(len=*), parameter :: myname = 'GetLine'
    character(len=k_ml), intent(out) :: line 
    integer, intent(in) ::uno
    integer :: errno
    inquire(unit=uno,iostat=errno)
    GetLine=.false.
    if (errno/=0) then
      write(*,*)"Unexpected error opening the input file(s)"
      stop
    endif
    read(uno,fmt='(a)',iostat=errno)line
    if (errno==0) GetLine=.true.
  end function GetLine

!> \brief adds a new node in the list
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param word character(len=*), value for field nam of the node
!> \param current pointer to the last node in the list
!> \param lines integer, optional value for field lines 

  subroutine AddName(word,current,lines)
    character(len=*), parameter :: myname = 'AddName'
    character(len=*),intent(in) :: word 
    integer,intent(in),optional :: lines
    type(names),pointer :: current
    type(names),pointer :: node

    allocate(node)
    node%nam=trim(word)
    if(present(lines))  node%lines=lines
    current%next=>node
    current=>node

  end subroutine AddName

!> \brief logical function finds a field in a list starting at root
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param word character(len=*), value of the field nam to be searched for
!> \param root type(names), pointer, starting point in search
!> \param loc tppe(names), pointer, optional returns the location where the info was found.
  logical function FindName(word,root,loc)
    character(len=*), parameter :: myname = 'FindName'
    character(len=*),intent(in) :: word
    type(names), pointer :: root
    type(names), pointer, optional :: loc
    type(names), pointer :: current

    current=>root
    FindName=.false.
    do while(associated(current))
      if (cstr(trim(word),trim(current%nam))) then
        FindName=.true.
        if (present(loc)) loc=>current
        exit
      endif
      current=>current%next
    end do
  end function FindName

!> \brief Prints at a specified unit the nam field  from a list
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param root type(names), pointer  starting node in the Printing list 
!> \param ioLoc type(ioType) conatins the unit to Print (see m_Types::ioType)
!> \remarks
  subroutine PrintName(root,ioLoc)
    character(len=*), parameter :: myname = 'PrintName'
    type(names), pointer :: root
    type(ioType), intent(inout) :: ioLoc
    type(names),pointer :: current

    current=>root
    do while(associated(current))
      write(ioLoc%uerr,'(a)') trim(current%nam)
      current=>current%next
    end do
  end subroutine PrintName

!> \brief recursive function that deallocates all the nodes starting with root
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param root type(names), pointer the starting node
!> \remarks
  recursive subroutine DeleteList(root)
    character(len=*), parameter :: myname = 'DeleteList'
    type(names), pointer :: root
    type(names), pointer :: current

    current=>root%next
    if (associated(current)) then
      call DeleteList(current)
    else
      deallocate(root)
    endif 
  end subroutine DeleteList

!> \brief parses a line removing comments
!> \details everything after ! \# or // is considered a comment
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param line character(len=*) the line to be parsed
  character(k_ml) function ParseLine(line)
    character(len=*), parameter :: myname = 'ParseLine' 
    character(len=*), intent(in) :: line
    integer ::p(1:3),pos

    p(1)=scan(line,"!")
    p(2)=scan(line,"#")
    p(3)=index(line,"//")
    where(p==0) p=k_ml+1
    pos=minval(p)-1
    if (pos/=k_ml) then
      ParseLine=line(1:pos)
    else
      ParseLine=line
    endif
  end function ParseLine

!> \brief closes all the units that where opened for reading or writing during the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param ioLoc type(ioType) contain the staring unit for the list of temporary units (see m_Types::ioType)
!> \remarks
  subroutine finalize(ioLoc)
    character(len=*), parameter :: myname = 'finalize'
    type(ioType), intent(inout) :: ioLoc
    integer :: i
    logical :: op
    character(len=k_mw) :: filename

    do i=ioLoc%uerr+1,GetUnit()-1
      inquire(unit=i,opened=op,name=filename)
      if (op) close(unit=i)
    enddo
  end subroutine finalize

!> \brief Prints an error message and aborts the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param message the error message
!> \param routine the subprogram that generated the error
!> \param filename raeding this \em filename the error occured
!> \param lineno the line number that generated the error
!> \param ioLoc type(ioType) (see m_Types::ioType)
  subroutine ParseErr(message,routine,filename,lineno,ioLoc)
    character(len=*), parameter :: myname = 'ParseErr'
    character(len=*),intent(in) :: message,routine,filename
    integer, intent(in) :: lineno
    type(ioType), intent(inout) :: ioLoc

    write(ioLoc%uerr,'(a,a)')"Error: ",trim(message)
    write(ioLoc%uerr,'(a,a)')"Routine: ",trim(routine)
    write(ioLoc%uerr,'(a,a)')"File: ", trim(filename)
    write(ioLoc%uerr,'(a,i7)')"Line number: ",lineno
    write(ioLoc%uerr,'(a)')"User stop"
    write(*,'(a,a)')"User stop, but probably this is not what you want, check file: ",ioLoc%inpErr
    stop
  end subroutine ParseErr

!> \brief Prints a warning message occured in the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param message the error message
!> \param routine the subprogram that generated the warning
!> \param filename raeding this \em filename the warning occured
!> \param lineno the line number that generated the warning
!> \param ioLoc type(ioType) (see m_Types::ioType)

  subroutine ParseWar(message,routine,filename,lineno,ioLoc)
    character(len=*), parameter :: myname = 'ParseWar'
    character(len=*),intent(in) :: message,routine,filename
    integer, intent(in) :: lineno
    type(ioType), intent(inout) :: ioLoc

    write(ioLoc%uerr,'(a,a)')"Warning: ",trim(message)
    write(ioLoc%uerr,'(a,a)')"Routine: ",trim(routine)
    write(ioLoc%uerr,'(a,a)')"File: ", trim(filename)
    write(ioLoc%uerr,'(a,i7)')"Line number: ",lineno
    write(ioLoc%uerr,'(a)')"User warning"
    write(*,'(a,a)')"User warning, check file: ",trim(ioLoc%inpErr)
  end subroutine ParseWar

!> \brief  returns a logical value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to .true. if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  logical function GetLogical(ioLoc,label,dflt)
    character(len=*), parameter :: myname = 'GetLogical'
    character(len=*),intent(in) :: label
    type(ioType),intent(inout) :: ioLoc
    logical, intent(in),optional :: dflt
    type(names),pointer :: found
!default value is set to .true. if no dflt parameter is present

    if (FindName(trim(label),tnames,found)) then
      if ((cstr(trim(found%value),"yes")).or.(cstr(trim(found%value),"true")) &
        .or.(cstr(trim(found%value),".true.")).or.(cstr(trim(found%value),"t"))&
        .or.(cstr(trim(found%value),"y")).or.(cstr(trim(found%value),"1"))) then 
        GetLogical=.true.
      elseif ((cstr(trim(found%value),"no")).or.(cstr(trim(found%value),"false")) &
        .or.(cstr(trim(found%value),".false.")).or.(cstr(trim(found%value),"f")) &
        .or.(cstr(trim(found%value),"n")).or.(cstr(trim(found%value),"0"))) then 
        GetLogical=.false.
      else
        call error("Unsuported value: "//&
          trim(label)//" "//trim(found%value),myname,.true.,ioLoc)
      endif
      write(ioLoc%udeb,'(a,2x,l1)') trim(label),GetLogical
    else
      if (present(dflt)) then
        GetLogical=dflt
      else
        GetLogical=.true.
      endif 
      write(ioLoc%udeb,'(a,2x,l1,a)') trim(label),GetLogical, " default"
    endif
  end function GetLogical

!> \brief  returns a string value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to empty string if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \param Print optional logical determines if the token is Printed or not. if the parameter is missing is assuk_med .false.
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
character(len=k_mw) function GetString(ioLoc,label,dflt,Print)
    character(len=*), parameter :: myname = 'GetString'
    type(ioType),intent(inout) :: ioLoc
    character(len=*),intent(in) :: label
    character(len=*), intent(in),optional :: dflt
    logical, intent(in), optional :: Print
    type(names),pointer :: found
! default value if no dflt parameter is present is empty string
! so be sure that you provide a default parameter
    if (FindName(trim(label),tnames,found)) then
      GetString=trim(found%value)
      if (present(Print).and.(.not.Print)) then 
      else 
        write(ioLoc%udeb,'(a,2x,a)') trim(label),trim(GetString)
      endif   
    else
      if (present(dflt)) then
        GetString=trim(dflt)
      else
        GetString=''
      endif
      if (present(Print).and.(.not.Print)) then 
      else   
        write(ioLoc%udeb,'(a,2x,a,a)') trim(label),trim(GetString), " default"
      endif   
    endif  
  end function GetString

!> \brief  returns a real(kind=k_pr) value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to 0.0_k_pr if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  real(k_pr) function GetReal(ioLoc,label,dflt)
    character(len=*), parameter :: myname = 'GetReal'
    type(ioType),intent(inout) :: ioLoc
    character(len=*),intent(in) :: label
    real(k_pr), intent(in),optional :: dflt
    type(names),pointer :: found
    character(len=k_mw) :: aux
    integer :: errno
! default value if no dflt parameter is present is 0.0_k_pr
    if (FindName(trim(label),tnames,found)) then
      aux=trim(found%value)
      read(aux,*,iostat=errno)GetReal
      if (errno/=0) call error("wrong value "//trim(found%value)&
        //" suplied for label "//trim(label),myname,.true.,ioLoc)
      write(ioLoc%udeb,'(a,2x,f32.16)') trim(label),GetReal
    else
      if (present(dflt)) then
        GetReal=dflt
      else
        GetReal=0.0_k_pr
      endif
      write(ioLoc%udeb,'(a,2x,f32.16,a)') trim(label),GetReal, " default"
    endif
  end function GetReal

!> \brief  returns an integer value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to 0 if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  integer function GetInteger(ioLoc,label,dflt)
    character(len=*), parameter :: myname = 'GetInteger'
    type(ioType),intent(inout) :: ioLoc
    character(len=*),intent(in) :: label
    integer, intent(in),optional :: dflt
    type(names),pointer :: found
    character(len=k_mw) :: aux
    integer :: errno
! default value if no dflt parameter is present is 0

    if (FindName(trim(label),tnames,found)) then
      aux=trim(found%value)
      read(aux,*,iostat=errno)GetInteger
      if (errno/=0) call error("wrong value "//trim(found%value)&
        //" suplied for label "//trim(label),myname,.true.,ioLoc)
      write(ioLoc%udeb,'(a,2x,i0)') trim(label),GetInteger
    else
      if (present(dflt)) then
        GetInteger=dflt
      else
        GetInteger=0
      endif
      write(ioLoc%udeb,'(a,2x,i0,a)') trim(label),GetInteger, " default"
    endif
  end function GetInteger

!> \brief  returns a logical value to indicate if it found in the input file(s) a valid block associated with the token \em label
!> \details  if no token is found the value returned is .false.
!> a unit from where the content of the block cand be read is returned via nt output is put in the units indicated by ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param nt unit number from where to read the block data
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  logical function GetBlock(ioLoc,label,nt)
    character(len=*), parameter :: myname = 'GetBlock'
    type(ioType),intent(inout) :: ioLoc
    character(len=*),intent(in) :: label
    integer, intent(out) :: nt
    type(names),pointer :: found
    character(len=k_mw) :: line

    if (FindName(trim(label),bnames,found)) then
      GetBlock=.true.
      nt=found%ut
      line=found%value  
      open(unit=nt,action="read",status="old",file=trim(line))
      write(ioLoc%udeb,'(a,2x,l1)') trim(label),GetBlock
    else
      write(ioLoc%udeb,'(a)')"Block "//trim(label)//" was not found"
      GetBlock=.false.
    endif
  end function GetBlock

!> \brief ends the parsing process and cuts the \em tree (deallocates the memory used during the parsing process)
!> \author Alin M Elena
!> \date 14th of January 2006
!> \warning after the call to this function none of the get_* function will work
!> all the input has to be done before the call to it.
  subroutine EndParse
    character(len=*), parameter :: myname = 'EndParse'

! to be called only when there is nothing to be read
    type(names),pointer :: current
    integer :: nt
    character(len=k_mw) :: filename

    current=>bnames
    do while(associated(current))
      nt = current%ut
      filename=current%value
      close(nt)
      open(nt,file=trim(filename),status="old",action="read")
      close(nt,status="delete")
      current=>current%next
    end do
    call DeleteList(tnames)
    call DeleteList(bnames)
  end subroutine EndParse

end module m_Parser
