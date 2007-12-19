!> \mainpage TDTB+UJ Manual
!> Time-Dependent Tight-Binding+UJ \n
!>  \section intro Introduction
!>  The tight binding method is implemented in the spirit of Mike Finnis' book\n
!>  \em Interatomic \em forces \em in \em condensed \em matter, Oxford University Press, 2003\n
!> \section run Run methods
!> \section inp Input Variables
!> - \subpage ioVars "I/O Variables"
!> - \subpage control "Program Flow Control Variables"
!> \section format Format specifiers
!> - \subpage atoms "Atoms"
!> - \subpage species "Atomic Species Data"
!> - \subpage basis "Basis Set"
!> - \subpage delta "Multipoles parameters(Delta Block)"
!> \section tb Tight Binding Parameters
!>
!>  A pdf version of this documentation can be found at http://titus.phy.qub.ac.uk/Programs/TDTB_UJ/refman.pdf
!> \page code Coding Style Document
!> \section general General Style
!>- You \em SHOULD use the Object Oriented features of Fotran as much as possible.
!>  If you are constantly casting things or doing IsKindOf calls on objects, then your design is probably wrong.
!>- You should protect your data. All member variables \em SHOULD be protected and only accessible
!> (where apropriate) via methods. Furthermore, you \em SHOULD avoid exposing details of implementation that
!> might change. Supply a feature set of methods that could, if necessary,
!> have an entirely different implementation or platform
!>- When allocating memory on the heap you \em SHOULD  make sure it is explicit in the object design where objects need
!> to be deleted and whose responsibility it is. Use try/catch blocks where you are allocating memory
!> on the heap and clean up after exceptions.
!>- You \em SHOULD program defensively. Anticipate breakages.
!> Always test your code, especially with boundary conditions.
!> Consider writing tests before you write your code (the XP way).
!>- You \em SHOULD program defensively. Anticipate breakages.
!> Always test your code, especially with boundary conditions.
!> Consider writing tests before you write your code (the XP way).
!> \section file Files Organisation
!>- You \em SHOULD use one file per module where practical.
!> Small helper modules that are 'local' to a major module can be incorporated
!> in the same file but you should use separate files if the code gets too large.
!>-  You \em SHOULD name source files after the module name they implement and omit the leading 'm_'.
!> For example, the module m_Widget will reside in the file Widget.f90
!> \section notation Hungarian notation
!>- Some Hungarian notation to indicate the type of a variable in its name \em may be used.
!> Names start with one or more lower-case letters to indicate the type of variable. For example:
!> -# iCount integer
!> -# lCompleted logical
!> -# cNameOfTag character
!> -# pThing pointer (can be combined with the others to indicate a pointer to a specific type).
!> -# ppThing pointer to a pointer.
!>- Otherwise, you \em SHOULD use studly caps to separate words, not underscores. The first letter \em may be in caps.
!> If type identifier is present the letter after \em SHOULD be in caps.
!> \section globals Global variables and constants
!>- Constants \em SHOULD be all implemented in a separate module (eg. m_Constants) and \em SHOULD be prefixed by k_ (eg. k_pi),
!>- You \em SHOULD NOT use global variables (with obvious exceptions, like the sole instance of your primary object).
!> They are a source of bugs, and can usually be avoided.
!> If you do need to use global variables, you \em SHOULD prefix their names with 'g_' (for instance, g_iTotalHits).
!>\section module Module comments
!>- Each module \em SHOULD have a similar block comment preceding with Author, Date, and
!> description. The description should give any specific instantiation or initialisation notes,
!> and where possible, usage examples. You \em SHOULD use doxygen style comments.
!> \section moremodule Module naming convention
!>- Module names \em SHOULD start with 'm_' and use studly caps to break up words. Use descriptive module
!> names like m_ReadData or m_SlaterKoster
!>- You \em SHOULD NOT use other underscores than "m_" in module names.
!> \section Variables
!>- Variables of a module or subprograme \em SHOULD use the conventions above. For example:
!> -# integer NoOfOrbitals
!> -# real ElectronicEntropy
!>- There is no specific convention for local variables but the name \em SHOULD be descriptive.
!>- You \em SHOULD always initialise variables to a known initial value.
!>\section procedures Subprograms naming
!>- Subprograms \em SHOULD always start with a capital and \em SHOULD use capitals to brake words.
!> For example, if you have subprograms to get or set properties of your modules, use \b GetSomeProperty
!> and \b SetSomeProperty. For logical properties, \b IsSomeState is a better choice than \b Get.
!>\section data Data types naming
!>- User defined data types \em SHOULD use capitals to brake words and \em SHOULD use 'Type' as a suffix.(eg OrbitalsType)
!> \section statements Statements
!>- Each line \em SHOULD contain at most one statement do not use multi-statement lines.
!> An if statement \em SHOULD always have a then.Eg. \n
!> if (SomeLogicalExpression) then \n
!>   call PrintMe() \n
!> endif \n
!> rather than: \n
!> if (SomeLogicalExpression) call PrintMe() \n
!>- When using select case statements, you \em SHOULD explicitly comment any deliberate fall-through behaviour (where you deliberately omit the default case from the end of a case clause).
!> \section comment Comments
!>- You \em SHOULD use exclamation sign (!) for comments. You \em SHOULD put comments on lines starting in the 0 column and
!> avoid long lines. Use comments to explain why an action/procedure is required. What you are doing
!> \em SHOULD be obvious from the code itself.
!>- For commenting modules/subprograms/types definitions, you \em SHOULD use doxygen block comments.
!>- A comment should preferably look like this
!> \verbatim
!> \brief short description of the action
!> \details (optional) more details and notes
!> \author The Name
!> \date creation date (optional) time
!> \remarks (optional) extra comments
!> \endverbatim
!> \section whites White space and indentation
!>-  You \em SHOULD use spaces rather than tabs (2 spaces for a tab stops are recommended). More editors will
!> have options to assist with this.
!>- Blank lines \em SHOULD never appear more than one.
!>- You \em SHOULD avoid lines longer than 132 characters.
!>- You \em SHOULD break up long lists of parameters, indenting the broken lines to the right of the start of the first line.
!>  Call Procedure(BigName, AnotherParameter, GoshThisIsLong) \n
!> becomes: \n
!>  Call Procedure(BigName,\n
!>   AnotherParameter,\n
!>   GoshThisIsLong);
!>- You \em SHOULD use 2 spaces to indent the code. The end of a block should be on the same column with the start block.
!> The content of a block should be indented 2 spase with respect to the belonging code. Eg.
!> \verbatim
!> do i=1,n
!>   if (lCondition) then
!>     a=3
!>   endif
!> enddo
!> a=5
!> \endverbatim

!> \brief main program Time-Dependent Tight-Binding+UJ
!> \author Alin M Elena
program tbuj
  use m_Constants
  use m_Types
  use m_ReadData
  use m_Useful, only : DateAndTime, error
  use m_TightBinding
  use m_DriverRoutines
  use m_Fit
  implicit none


  integer :: narguments
  character(len=k_mw) :: arg
  character(len=10) :: dt
  character(len=12) :: tm
!
  type(ioType) :: ioInfo
  type(generalType) :: general
  type(atomicxType) :: atomicx
  type(modelType) :: tbModel
! solution spece variable
  type(solutionType) :: SolSpace


  call DateAndTime(dt,tm)
  narguments=iargc()
  call cpu_time(general%time%start)

! it reads the name of the input as inline argument
! if there is none the default name is inp

  if (narguments==1) then
    call getarg(1,arg)
    ioInfo%inpFile=arg
  else
    ioInfo%inpFile="inp"
  endif

  call Initialize(ioInfo,general,atomicx,tbModel)
  call SetSolutionSpace(ioInfo,general,atomicx,tbModel,SolSpace)

  if (general%runType == k_runSp) then
    call SinglePoint(ioInfo,general,atomicx,tbModel,SolSpace)
  elseif (general%runType == k_runGeomBFGS) then
!     call bfgs
  elseif (general%runType == k_runBO) then
    call BornOppenheimerDynamics(ioInfo,general,atomicx,tbModel,SolSpace)
  elseif (general%runType == k_runEhrenfest) then
    call EhrenfestDynamics(ioInfo,general,atomicx,tbModel,SolSpace)
  elseif (general%runType==k_runFit) then
    call Fitting(ioInfo,general,atomicx,tbModel,SolSpace)
  elseif (general%runType==k_runForceTest) then
!     call forceTest
  elseif (general%runType==k_runForceTestx) then
!     call forceTestx
  elseif (general%runType==k_runForceTesty) then
!     call forceTesty
  elseif (general%runType==k_runForceTestz) then
!     call forceTestz
  elseif (general%runType==k_runEhrenfestDamped) then
    call EhrenfestDynamicsDamped(ioInfo,general,atomicx,tbModel,SolSpace)
  elseif (general%runType==k_runFragments) then
!     call fragments_k_run
  elseif (general%runType==k_runSpecial) then
!     call special_k_run
  else
    call error("RunType not implemented",'TDTB+UJ',.true.,ioInfo)
  endif

!!!!!! closes all the units associated with blocks and deallocated the trees for tokens and blocks
  call CleanMemory(ioInfo,atomicx,general,tbModel,SolSpace)
  call cpu_time(general%time%end)
  write(ioInfo%uout,'(a,a,a,a)',advance="no")"Program TDTB+UJ has started at ",dt," ",tm
  call DateAndTime(dt,tm)
  write(ioInfo%uout,'(a,a,a,a)')" ended at ",dt," ",tm
  write(ioInfo%udeb,'(a,f16.6,a)')"Program TDTB+UJ has run for "&
    ,general%time%end-general%time%start," seconds"
  write(ioInfo%uout,'(a,f0.6,a)')"Program TDTB+UJ has run for "&
    ,general%time%end-general%time%start," seconds"

  call CloseIoGeneral(ioInfo)

end program tbuj
