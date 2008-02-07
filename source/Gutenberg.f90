!> \brief Prints complicated structures
!> \details all the Prints that get used more than once or are larger than a line
!> should be here
!> \author Alin M Elena
!> \date 02/11/07, 18:13:15
module m_Gutenberg
  use m_Constants
  use m_Types
  use m_Useful
  implicit none

  private

  public :: PrintTail
  public :: PrintAtoms
  public :: PrintSpecies
  public :: PrintBasis
  public :: PrintTbGSP
  public :: PrintTbHarrison
  public :: PrintDelta
  public :: PrintMagneticMoment
  public :: PrintMatrix
  public :: PrintVector
  public :: PrintForces
  public :: PrintCharges
  public :: PrintDipoles
  public :: PrintAtomMatrix
  public :: PrintCoordinates
  public :: PrintVelocities
  public :: PrintXYZ
contains

!> \brief Prints the tail parameters
!> \author Alin M Elena
!> \date 02/11/07, 09:00:08
!> \param io type(ioType) i/o units
!> \param tbmod type(modelType) model parameters
!> \param atomix type(atomicxType) atomic details
  subroutine PrintTail(atomix,tbmod,io)
    character(len=*), parameter :: sMyName="PrintTail"
    type(ioType), intent(in) :: io
    type(atomicxType), intent(in) :: atomix
    type(modelType),intent(in) :: tbmod

    integer :: i,j
    
    integer :: k,k1,k2
    write(io%uout,'(a)')"==Tail functions========================================================================================"
    write(io%uout,'(a)')"==Repulsive terms tail parameters======================================================"
        write(io%uout,'(2a4,3a20,2a10)')"Sp1","Sp2","a","b","c","r_1","r_cut"
    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        write(io%uout,'(2i4,3f20.8,2f10.4)')i,j,tbMod%hopping(i,j)%repTail%a,&
          tbMod%hopping(i,j)%repTail%b,tbMod%hopping(i,j)%repTail%c,&
          tbMod%hopping(i,j)%repTail%rIn,tbMod%hopping(i,j)%repTail%rOut
      enddo
    enddo
        write(io%uout,'(a)')"==End Repulsive terms tail parameters================================================"
    write(io%uout,'(a)')"==Radial terms tail parameters====================================================================="
        write(io%uout,'(a10,9x,3a20,2a10)')"hopping","a","b","c","r_1","r_cut"
    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
              write(io%uout,'(a9,f10.4,3f20.8,2f10.4)')trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2),&
                  tbMod%hopping(i,j)%hmnTail(k,k1,k2)%a,&
                  tbMod%hopping(i,j)%hmnTail(k,k1,k2)%b,tbMod%hopping(i,j)%hmnTail(k,k1,k2)%c,&
                  tbMod%hopping(i,j)%hmnTail(k,k1,k2)%rIn,tbMod%hopping(i,j)%hmnTail(k,k1,k2)%rOut
            enddo
          enddo
        enddo
      enddo
    enddo
    write(io%uout,'(a)')"==End Radial terms tail parameters================================================================="
    write(io%uout,'(a)')"==End Tail functions===================================================================================="
  end subroutine PrintTail

!> \brief Prints coordinates, velocities and type of atoms
!> \author Alin M Elena
!> \date 30/10/07, 13:22:04
!> \param io type(ioType) i/o units
!> \param general type(generalType) general data
!> \param atomix type(atomicxType) contains info about atoms

  subroutine PrintAtoms(io,general,atomix)
    character(len=*), parameter :: sMyName="PrintAtoms"
    type(ioType), intent(inout) :: io
    type(generalType), intent(in)  :: general
    type(atomicxType), intent(in) :: atomix
    integer :: i

    call PrintCoordinates(io,atomix)
    if (general%ReadVelocity) then
      call PrintVelocities(io,atomix)
    endif
    write(io%uout,'(a)')"=AtomicLists=========================================="
    if(atomix%atoms%nscf/=0) then
      write(io%uout,'(a)')"SCFAtoms:"
      do i=1,atomix%atoms%nscf
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%scf(i),"(",symbol(GetZ(atomix,atomix%atoms%scf(i))),") "
      enddo
      write(io%uout,*)
    endif
    if(atomix%atoms%nmoving/=0) then
      write(io%uout,'(a)')"Moving Atoms:"
      do i=1,atomix%atoms%nmoving
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%moving(i),"(",symbol(GetZ(atomix,atomix%atoms%moving(i))),") "
      enddo
      write(io%uout,*)
    endif
    write(io%uout,'(a,i0)')"Atoms in Donor Group(NDonor): "&
        ,atomix%atoms%ndonor
    if(atomix%atoms%ndonor/=0) then
      write(io%uout,'(a)')"Donor Atoms:"
      do i=1,atomix%atoms%ndonor
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%donor(i),"(",symbol(GetZ(atomix,atomix%atoms%donor(i))),") "
      enddo
      write(io%uout,*)
    endif
    write( io%uout,'(a,i0)')"Atoms in Acceptor Group(NAcceptor): "&
        ,atomix%atoms%nacceptor
    if(atomix%atoms%nacceptor/=0) then
      write(io%uout,'(a)')"Acceptor Atoms:"
      do i=1,atomix%atoms%nacceptor
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%acceptor(i),"(",symbol(GetZ(atomix,atomix%atoms%acceptor(i))),") "
      enddo
      write(io%uout,*)
    endif
    write(io%uout,'(a,i0)')"Atoms in Spacer Group(NSpacer): "&
        ,atomix%atoms%nspacer
    if(atomix%atoms%nspacer/=0) then
      write(io%uout,'(a)')"Spacer Atoms:"
      do i=1,atomix%atoms%nspacer
        write(io%uout,'(i0,a,a,a)',advance="no")atomix%atoms%spacer(i),"(",symbol(GetZ(atomix,atomix%atoms%spacer(i))),") "
      enddo
      write(io%uout,*)
    endif
    write(io%uout,'(a)')"=EndLists============================================="
  end subroutine PrintAtoms

!> \brief Prints velocities of atoms
!> \author Alin M Elena
!> \date 30/10/07, 13:22:04
!> \param io type(ioType) i/o units
!> \param atomix type(atomicxType) contains info about atoms
  subroutine PrintVelocities(io,atomix)
    character(len=*), parameter :: sMyName="PrintVelocities"
    type(ioType), intent(inout) :: io
    type(atomicxType), intent(in) :: atomix
    integer :: i
    write(io%uout,'(a)')  "=StartVelocitiesData==========================================="
    write(io%uout,'(a)')  "      Id  Sp El              VX              VY              VZ"
    write(io%uout,'(a)')  "---------------------------------------------------------------"
    do i=1,atomix%atoms%natoms
      write(io%uout,'(i8,i4,1x,a2,3f16.8,5x)')atomix%atoms%id(i),atomix%atoms%sp(i),symbol(GetZ(atomix,i)),&
            atomix%atoms%vx(i),atomix%atoms%vy(i),atomix%atoms%vz(i)
    enddo
    write(io%uout,'(a)')  "=EndVelocitiesData============================================="

  end subroutine PrintVelocities

!> \brief Prints coordinates of atoms
!> \author Alin M Elena
!> \date 30/10/07, 13:22:04
!> \param io type(ioType) i/o units
!> \param atomix type(atomicxType) contains info about atoms
  subroutine PrintCoordinates(io,atomix)
    character(len=*), parameter :: sMyName="PrintCoordinates"
    type(ioType), intent(inout) :: io
    type(atomicxType), intent(in) :: atomix
    integer :: i

    write(io%uout,'(a)')  "=StartAtomsData==============================================================================="
    write(io%uout,'(a)')  "      Id  Sp El               X               Y               Z            Bias IsSCF IsMoving"
    write(io%uout,'(a)')"----------------------------------------------------------------------------------------------"
    do i=1,atomix%atoms%natoms
      write(io%uout,'(i8,i4,1x,a2,4f16.8,5x,l1,8x,l1)')atomix%atoms%id(i),atomix%atoms%sp(i),symbol(GetZ(atomix,i)),&
            atomix%atoms%x(i),atomix%atoms%y(i),atomix%atoms%z(i),atomix%atoms%bias(i),&
            atomix%atoms%isscf(i),atomix%atoms%ismoving(i)
    enddo
    write(io%uout,'(a)')   "=EndAtomsData================================================================================="
  end subroutine PrintCoordinates

!> \brief Prints the info about species
!> \details at high verbosity some intristing info is Printed about each specie
!> \author Alin M Elena
!> \date 29/10/07, 19:59:04
!> \param io type(ioType) i/o units
!> \param specs type(speciesType) species data

  subroutine PrintSpecies(io,specs)
    character(len=*), parameter :: sMyName="PrintSpecies"
    type(ioType), intent(in) :: io
    type(speciesType), intent(inout) :: specs
    integer :: i
    write(io%uout,'(a)')"=SpeciesData============================================================================"
    write(io%uout,'(a)')"  Id  Z el    Zval            Mass          Ulocal          Jlocal          Uinter Norbs"
    write(io%uout,'(a)')"----------------------------------------------------------------------------------------"
    do i=1,specs%nspecies
      write(io%uout,'(i4,i3,1x,a,f8.4,4f16.8,i5)')specs%id(i),specs%z(i),symbol(specs%z(i)),specs%zval(i),specs%mass(i),&
              specs%ulocal(i),specs%jlocal(i),specs%uinter(i),specs%norbs(i)
      specs%mass(i)=specs%mass(i)*k_amuToInternal
    enddo
    write(io%uout,'(a)')"=EndSpeciesData========================================================================="

    if (io%verbosity>=k_highVerbos) then
      do i=1, specs%nspecies
        write(io%uout,'(a)')"========================================================"
        write(io%uout,'(a,i0,a)')"Specie: ",specs%id(i)," detailed info"
        write(io%uout,'(a,a,a,a)')"name: ",trim(ElName(specs%z(i)))," symbol: ",trim(symbol(specs%z(i)))
        write(io%uout,'(a,i0,a,f0.3,a,f0.3,a,f0.3,a)')"Z=",specs%z(i)," Zvalence=",specs%zval(i),&
          " mass=",specs%mass(i)," (",weight(specs%z(i)),") a.u."
        write(io%uout,'(a,a,a,i0,a,a)')"period: ",trim(period(specs%z(i)))," group: ",group(specs%z(i)),&
          " electronic config:",trim(ElConfig(specs%z(i)))
        write(io%uout,'(a,f0.3,a,f0.3,a)')"van der Waals radius=",VdWRadius(specs%z(i)),&
          "A covalent radius=",CovalentRadius(specs%z(i)),"A"
      enddo
      write(io%uout,'(a)')"========================================================"
    endif
  end subroutine PrintSpecies


!> \brief Prints the basis set information
!> \author Alin M Elena
!> \date 30/10/07, 17:56:02
!> \param io type(ioType) i/o units
!> \param gen type(generalType) general data
!> \param atomix type(atomicxType) contains info about atoms
  subroutine PrintBasis(io,gen,atomix)
    character(len=*), parameter :: sMyName="PrintBasis"
    type(ioType), intent(inout) :: io
    type(generalType), intent(in)  :: gen
    type(atomicxType), intent(inout) :: atomix
    integer :: i,j,sp
    character(len=1) :: spin=""
    character(len=2) :: el=""
    write(io%uout,'(a)')"==BasisSetInfo==================================="
    write(io%uout,'(a)')"==Listed by species==================="
    write(io%uout,'(a)')"  #Orb  Sp El  N  L  M Spin Occupation"
    do i=1,atomix%species%nspecies
      sp=atomix%species%id(i)
      do j=1,atomix%species%norbs(i)
        if (gen%spin) then
          if (atomix%speciesBasis(i,j)%spin) then
            spin="U"
          else
            spin="D"
          endif
        endif
        write(io%uout,'(i6,i4,1x,a,3i3,4x,a,1x,f0.8)')j,atomix%speciesBasis(i,j)%sp,&
              symbol(atomix%species%z(atomix%speciesBasis(i,j)%sp)),&
            atomix%speciesBasis(i,j)%n,atomix%speciesBasis(i,j)%l,atomix%speciesBasis(i,j)%m,&
                spin,atomix%speciesBasis(i,j)%occup
      enddo
    enddo
    write(io%uout,'(a)')"==End Listed by species==============="
    write(io%uout,'(a)')"==Listed by atoms============================"
    write(io%uout,'(a)')"  #Orb  Atom  Sp  El  N  L  M Spin Occupation"
    spin=""
    el="el"
    do i=1,atomix%basis%norbitals
      if (atomix%basis%orbitals(i)%spin) then
        spin="U"
      else
        spin="D"
      endif
      el=symbol(GetZ(atomix,atomix%basis%orbitals(i)%atom))
      write(io%uout,'(2i6,i4,a4,3i3,1x,a4,1x,f0.8)')i,atomix%basis%orbitals(i)%atom,atomix%basis%orbitals(i)%sp,el,&
            atomix%basis%orbitals(i)%n,atomix%basis%orbitals(i)%l,atomix%basis%orbitals(i)%m,spin,&
            atomix%basis%orbitals(i)%occup
    enddo
    write(io%uout,'(a)')"==End Listed by atoms========================"
    write(io%uout,'(a)')"==EndBasisSetInfo================================"
  end subroutine PrintBasis

!> \brief Prints the parameters of the gsp model that will be used
!> \author Alin M Elena
!> \date 31/10/07, 16:55:21
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine PrintTbGSP(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="PrintTbGSP"
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomix
    type(modelType), intent(inout) :: tbMod
    integer :: i,j,i1,k,k1,k2
    character(len=k_mw) :: read_var

    write(io%uout,'(a)')"==GSP Model======================================="
    write(io%uout,*) &
        "Using tight binding model as described in:"
    write(io%uout,*) &
        "A.P. Horsfield,P.D. Godwin,D.G. Pettifor,A.P. Sutton"
    write(io%uout,*) &
        "Computational materials synthesis. I. A tight-binding schek_me"
    write(io%uout,*) &
        "for hydrocarbons - Phys. Rev. B 54, 22, 15 773"
    write(io%uout,'(a)')"===================================================="

    do i=1,atomix%species%nspecies
      do j=1,atomix%species%nspecies
        write(io%uout,'(a,2i2)')"# Species",i,j
        if (i.eq.j) then
          do i1=0,atomix%species%norbs(i)-1
            read_var="eps"
            read_var=trim(ccvar(i,i,read_var))
            write(io%uout,*) trim(ccvar(i1,i1,read_var)),tbMod%hopping(i,j)%eps(i1)
          enddo
          if (gen%embedding) then
            read_var="a1"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a1
            read_var="a2"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a2
            read_var="a3"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a3
            read_var="a4"
            write(io%uout,*) &
              trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a4
          endif
        endif
        read_var="phi0"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%phi0
        read_var="r0"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%r0
        read_var="rc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%rc
        read_var="r1"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%r1
        read_var="rcut"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%rcut
        read_var="n"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%n
        read_var="nc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),  tbMod%hopping(i,j)%nc
        read_var="d0"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%d0
        read_var="dc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%dc
        read_var="d1"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%d1
        read_var="dcut"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%dcut
        read_var="m"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),  tbMod%hopping(i,j)%m
        read_var="mc"
        write(io%uout,*) &
          trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%mc
        do k=0,tbMod%hopping(i,j)%l1
          do k1=0,tbMod%hopping(i,j)%l2
            do k2=0,min(k,k1)
              write(io%uout,*) &
              trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2)
            enddo
          enddo
        enddo
      enddo
    enddo
    write(io%uout,*) &
      "========================Finish========================="
  end subroutine PrintTbGSP


  !> \brief Prints the initialization routine for the Harrison tight binding model
!> \author Alin M Elena
!> \date 31/10/07, 23:20:15
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters
  subroutine PrintTbHarrison(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="PrintTbHarrison"
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomix
    type(modelType), intent(inout) :: tbMod

    integer :: i,j,i1,k,k1,k2
    character(len=k_mw) :: read_var

    write(io%uout,*)"==Harrison Parameters================================================"
    write(io%uout,*) &
        "Using Harrison type tight binding model "
      do i=1,atomix%species%nspecies
        do j=1,atomix%species%nspecies
          write(io%uout,'(a,2i2)')"#",i,j
          if (i==j) then
            do i1=0,atomix%species%norbs(i)-1
              read_var="eps"
              read_var=trim(ccvar(i,i,read_var))
              write(io%uout,*) trim(ccvar(i1,i1,read_var)),tbMod%hopping(i,j)%eps(i1)
            enddo
            if (gen%embedding) then
              read_var="a1"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a1
              read_var="a2"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a2
              read_var="a3"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a3
              read_var="a4"
              write(io%uout,*) &
                trim(ccvar(i,j,read_var)), tbMod%hopping(i,j)%a4
            endif
          endif
          read_var="phi0"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%phi0

          read_var="r1"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%r1
          read_var="rcut"
          write(io%uout,*) &
              trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%rcut
          read_var="n"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%n
          read_var="d1"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%d1
          read_var="dcut"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%dcut
          read_var="m"
          write(io%uout,*) &
            trim(ccvar(i,j,read_var)),tbMod%hopping(i,j)%m
          do k=0,tbMod%hopping(i,j)%l1
            do k1=0,tbMod%hopping(i,j)%l2
              do k2=0,min(k,k1)
                write(io%uout,*) &
                trim(ccnlm(i,j,k,k1,k2)),tbMod%hopping(i,j)%a(k,k1,k2)
              enddo
            enddo
          enddo
        enddo
      enddo

      write(io%uout,*)"==End Harrison Parameters============================================"
  end subroutine PrintTbHarrison

!> \brief Prints the delta block
!> \details the block is read only if multipoles k_electrostatics is asked for
!> \author Alin M Elena
!> \date 31/10/07, 17:52:20
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomix type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param tbMod type(modelType) contains iformation about the tight binding model parameters

  subroutine PrintDelta(io,gen,atomix,tbMod)
    character(len=*), parameter :: sMyName="PrintDelta"
    type(ioType), intent(inout) :: io
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomix
    type(modelType), intent(inout) :: tbMod
    integer :: i,j,k,l

    write(io%uout,*) &
      '==DeltaPol Block====================================================='
    write(io%uout,*) "Note that: delta (l,l1,l2) with |l1-l2|<=l<=|l1+l2|"
    do i = 1,atomix%species%nspecies
      write(io%uout,'(a,i2)') "Species ",tbMod%delta(i)%sp
      do j=0,tbMod%delta(i)%l
        do k=j,tbMod%delta(i)%l
          if (k==j) then
            write(io%uout,'(100(a,3i3,a,f0.8))')(' delta(',l,j,k,') = ', tbMod%delta(i)%d(l,j,k),l=abs(j-k),abs(j+k))
          else
            write(io%uout,'(100(a,3i3,a,3i3,a,f14.8))')(' delta(',l,j,k,') = delta(',l,k,j,') = ', &
                  tbMod%delta(i)%d(l,j,k),l=abs(j-k),abs(j+k))
          endif
        enddo
      enddo
    enddo
    write(io%uout,*) &
        '==End DeltaPol Block================================================='
  end subroutine PrintDelta

!> \brief Prints the magnetic moment
!> \author Alin M Elena
!> \date 03/11/07, 12:48:24
!> \param io type(ioType) contains all the info about I/O files
!> \param atomic type(atomicType) contains all info about the atoms and basis set and some parameters
!> \param sol type(solutionType) contains information about the solution space
!> \param blocal logical \internal complete the description of parameters
  subroutine PrintMagneticMoment(atomic,sol,blocal,io)
    character(len=*), parameter :: myname="PrintMagneticMoment"

    type(atomicxType), intent(inout) :: atomic
    type(solutionType), intent(inout) :: sol
    logical, intent(in) :: blocal
    type(ioType), intent(in) :: io

    integer :: i
    real(k_pr) :: lm
    write(io%uout,"(a)")"Atom | Specie | Magnetic Moment |"
    if (blocal) then
      do i=1,atomic%atoms%natoms
        write(io%uout,"(i4,1x,i6,1x,f16.8)")i,atomic%atoms%sp(i),atomic%atoms%MagMom(i)
      enddo
      do i=1,atomic%atoms%natoms
        lm=LocalMoment(atomic%atoms%id(i),atomic,.true.,io,sol)
      enddo
    else
      do i=1,atomic%atoms%natoms
        write(io%uout,"(i4,1x,i6,1x,f16.8)")i,atomic%atoms%sp(i),atomic%atoms%MagMom(i)
      enddo
    endif
    write(io%uout,"(a,f16.8)")"Total magnetic moment: ",atomic%atoms%MagneticMoment
  end subroutine PrintMagneticMoment


!> \brief prints a matrix from matrixType
!> \author Alin M Elena
!> \date 07/11/07, 11:24:28
!> \param mat type(matrixType) the matrixType
!> \param label character a label that should be added as a caption
!> \param io type(ioType) i/o units
!> \param isComplex logical optional on/off complex part if present and true complex part is printed
  subroutine PrintMatrix(mat,label,io,isComplex)
    character(len=*), parameter :: myname = 'PrintMatrix'
    type(matrixType), intent(in) :: mat
    character(len=*), intent(in) :: label
    type(ioType), intent(in) :: io
    logical, intent(in),optional :: isComplex
    integer :: i,j,m


      write(io%uout,*) "==",trim(label),"======"
      write(io%uout,*) '==Real part==='
      m=mat%dim
      write(io%uout,'(7x)',advance="no")
      do i=1,m
         write(io%uout,'(i12,1x)',advance="no") i
      enddo
      write(io%uout,*)
      do i=1,m
         write(io%uout,'(i6,1x)',advance="no") i
         do j=1,m
            write(io%uout,'(f12.6,1x)',advance="no") real(mat%a(i,j))
         enddo
         write(io%uout,'(1x,i0)')i
      enddo
      write(io%uout,'(7x)',advance="no")
      do i=1,m
         write(io%uout,'(i12,1x)',advance="no") i
      enddo
      write(io%uout,*)
      if (present(isComplex).and.isComplex) then
        write(io%uout,*) '==Imaginary part==='
        write(io%uout,'(7x)',advance="no")
        do i=1,m
          write(io%uout,'(i12,1x)',advance="no") i
        enddo
        write(io%uout,*)
        do i=1,m
          write(io%uout,'(i6,1x)',advance="no") i
          do j=1,m
              write(io%uout,'(f12.6,1x)',advance="no") aimag(mat%a(i,j))
          enddo
          write(io%uout,'(1x,i0)')i
        enddo
        write(io%uout,'(7x)',advance="no")
        do i=1,m
          write(io%uout,'(i12,1x)',advance="no") i
        enddo
        write(io%uout,*)
      endif
      write(io%uout,*) "==End ",trim(label),"======"
   end subroutine PrintMatrix

!> \brief prints a vector
!> \details it can be printed on row or column, with index or without
!> \author Alin M Elena
!> \date 07/11/07, 13:13:52
!> \param vec real array to be printed
!> \param label character label to be printed
!> \param lIsRow logical if true is printed on row else as a column
!> \param lIsIndex logical if true the index is printed else no index
!> \param io type(ioType) i/o units
  subroutine PrintVector(vec,label,lIsRow,lIsIndex,io)
    character(len=*), parameter :: sMyName="PrintVector"
    real(k_pr), intent(in) :: vec(:)
    character(len=*), intent(in) :: label
    logical, intent(in) :: lIsRow, lIsIndex
    type(ioType),intent(in) :: io
    integer :: i

    write(io%uout,'(a,a,a)')"==",trim(label),"==="
    if (lisRow .and. lIsIndex) then
      do i=1,size(vec)
        write(io%uout,'(a,1x,i0,1x,f0.8,a)',advance="no")"(",i,vec(i),")"
      enddo
        write(io%uout,*)
    elseif (lisRow .and. (.not.lIsIndex)) then
      do i=1,size(vec)
        write(io%uout,'(f0.8,1x)',advance="no")vec(i)
      enddo
        write(io%uout,*)
    elseif ((.not.lisRow) .and. (.not.lIsIndex)) then
      do i=1,size(vec)
        write(io%uout,'(f16.8)')vec(i)
      enddo
    else
      do i=1,size(vec)
        write(io%uout,'(i5,f16.8,1x)')i,vec(i)
      enddo
    endif
    write(io%uout,'(a,a,a)')"==End ",trim(label),"==="
  end subroutine PrintVector


!> \brief prints the forces
!> \author Alin M Elena
!> \date 07/11/07, 17:18:42
!> \param io type(ioType) model parameters
!> \param atomic type(atomicType) atoms details
  subroutine PrintForces(atomic,io)
    character(len=*), parameter :: myname = 'PrintForces'
    type(atomicType), intent(inout) :: atomic
    type(ioType), intent(in) :: io
    integer :: i
    !-------------------------------------------------!

    write(io%uout,*) &
      "==atomic forces================================================="
    write(io%uout,*) &
      "    i      FX         FY         FZ      species   IsMoving"
    write(io%uout,*) &
      "----------------------------------------------------------------"
    do i=1,atomic%natoms
      write(io%uout,'(i6,3f11.5,i7,11x,l1)') &
      i,atomic%fx(i),atomic%fy(i),atomic%fz(i),atomic%sp(i),atomic%isMoving(i)
    end do
    write(io%uout,*) &
      "----------------------------------------------------------------"
    write(io%uout,'(a,3f11.5)')'Sum   ',sum(atomic%fx(1:atomic%natoms)),&
      sum(atomic%fy(1:atomic%natoms)),sum(atomic%fz(1:atomic%natoms))
    write(io%uout,*) &
      "==End atomic forces ============================================"

  end subroutine PrintForces

!> \brief prints the charges and excess charges
!> \author Alin M Elena
!> \date 08/11/07, 09:57:05
!> \param io type(ioType) contains all the info about I/O files
!> \param gen type(generalType) contains the info needed by the program to k_run
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param i integer optional, if present only the data about i atom will be printed
  recursive subroutine PrintCharges(gen,atomic,io,i)
    character(len=*), parameter :: myname = 'PrintCharges'
    integer, optional,intent(in) :: i
    type(generalType), intent(inout) :: gen
    type(atomicxType), intent(inout) :: atomic
    type(ioType), intent(inout) :: io
    integer :: m,n,j
    real(k_pr) :: ch=0.0_k_pr

    n=atomic%basis%norbitals
    m=n*(n-1)/2
    if (present(i)) then
      write(io%uout,'(i5,1x,i8,1x,f16.8,1x,f16.8)')i,atomic%atoms%sp(i),atomic%atoms%chrg(i),&
          atomic%atoms%chrg(i)+atomic%species%zval(atomic%atoms%sp(i))
    else
      write(io%uout,'(a)')"Atom | Specie | Charge execess |    Total Charge |"
      ch=0.0_k_pr
        do j=1,atomic%atoms%natoms
          call PrintCharges(gen,atomic,io,j)
          ch=ch+atomic%atoms%chrg(j)
        enddo
      write(io%uout,'(a,f16.8)') &
          "Excess charge on all atoms ",ch
    endif
   end subroutine PrintCharges

!> \brief prints the dipole moment
!> \author Alin M Elena
!> \date 08/11/07, 10:47:22
!> \param io type(ioType) contains all the info about I/O files
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
  subroutine PrintDipoles(atomic,io)
    character(len=*), parameter :: myname = 'PrintDipoles'
    type(atomicxType), intent(inout) :: atomic
    type(ioType), intent(inout) :: io
    real(k_pr) :: u2SI,u2D, ax,ay,az,ld
    integer :: at
    ax=0.0_k_pr
    ay=0.0_k_pr
    az=0.0_k_pr
    u2SI=k_length2SI*k_charge2SI
    u2D=u2SI/k_debye2SI
!! point charges dipoles
    write(io%uout,'(a)')" Dipoles from point charges =========================================="
    write(io%uout,'(a)')"  Atom |Specie|  Dipole x    |  Dipole y   |    Dipole z     | Units |"
    do at=1, atomic%atoms%natoms
      write(io%uout,'(2i7,3f16.8,a)')at,atomic%atoms%sp(at),atomic%atoms%x(at)*atomic%atoms%chrg(at)*u2D,&
            atomic%atoms%y(at)*atomic%atoms%chrg(at)*u2D,atomic%atoms%z(at)*atomic%atoms%chrg(at)*u2D," Debye"
      ax=ax+atomic%atoms%x(at)*atomic%atoms%chrg(at)*u2D
      ay=ay+atomic%atoms%y(at)*atomic%atoms%chrg(at)*u2D
      az=az+atomic%atoms%z(at)*atomic%atoms%chrg(at)*u2D
    enddo
    ld=sqrt(ax*ax+ay*ay+az*az)
    write(io%uout,'(a,4f16.8,a)')"dipole moment ",ax,ay,az,ld," Debye"
  !  avoid division by zero
    if (abs(ld) > epsilon(ax)) then
      write(io%uout,'(a,f9.6,1x,f9.6,1x,f9.6,a)')"dipole moment orientation (",ax/ld,ay/ld,az/ld, ") unit vector"
    endif
    write(io%uout,'(a)')" ====================================================================="
  !! dipoles in local units
    write(io%uout,'(a)')"  Atom |Specie|  Dipole x    |  Dipole y   |    Dipole z     | Units |"
    do at=1, atomic%atoms%natoms
      write(io%uout,'(2i7,3g16.8,a)')at,atomic%atoms%sp(at),atomic%atoms%dx(at),atomic%atoms%dy(at),atomic%atoms%dz(at)
    enddo
    write(io%uout,'(a,3g16.8)')"dipole moment ",atomic%atoms%tdipx,atomic%atoms%tdipy,atomic%atoms%tdipz
  !! dipoles in SI units
        write(io%uout,'(a)')"  Atom |Specie|  Dipole x    |  Dipole y   |    Dipole z     | Units |"
    do at=1, atomic%atoms%natoms
      write(io%uout,'(2i7,3g16.8,a)')at,atomic%atoms%sp(at),atomic%atoms%dx(at)*u2SI,atomic%atoms%dy(at)*u2SI,&
          atomic%atoms%dz(at)*u2SI," Cm"
    enddo
    write(io%uout,'(a,3g16.8,a)')"dipole moment ",atomic%atoms%tdipx*u2SI,atomic%atoms%tdipy*u2SI,&
      atomic%atoms%tdipz*u2SI, " Cm"
  !! dipoles in Debye
    write(io%uout,'(a)')"  Atom |Specie|  Dipole x    |  Dipole y   |    Dipole z     | Units |"
    do at=1, atomic%atoms%natoms
      write(io%uout,'(2i7,3f16.8,a)')at,atomic%atoms%sp(at),atomic%atoms%dx(at)*u2D,atomic%atoms%dy(at)*u2D,&
            atomic%atoms%dz(at)*u2D," Debye"
    enddo
    ld=sqrt(atomic%atoms%tdipx*u2D*atomic%atoms%tdipx*u2D+atomic%atoms%tdipy*u2D*atomic%atoms%tdipy*u2D+&
      atomic%atoms%tdipz*u2D*atomic%atoms%tdipz*u2D)
    write(io%uout,'(a,4f16.8,a)')"dipole moment ",atomic%atoms%tdipx*u2D,atomic%atoms%tdipy*u2D,&
      atomic%atoms%tdipz*u2D,ld, " Debye"
    if (abs(ld) > epsilon(ax)) then
      write(io%uout,'(a,f9.6,1x,f9.6,1x,f9.6,a)')"dipole moment orientation(",&
          atomic%atoms%tdipx*u2D/ld,atomic%atoms%tdipy*u2D/ld,atomic%atoms%tdipz*u2D/ld, ") unit vector"
    endif
  end subroutine PrintDipoles


!> \brief selects from a matrix type(matrixType) only data belonging to a certain atoms
!> \author Alin M Elena
!> \date 09/11/07, 10:25:27
!> \param mat type(matrixType) the matrixType
!> \param label character a label that should be added as a caption
!> \param io type(ioType) i/o units
!> \param isComplex logical on/off complex part if present and true complex part is printed
!> \param i integer the atom
!> \param atomic type(atomicxType) atomic details
  subroutine PrintAtomMatrix(i,atomic,mat,label,io,isComplex)
    character(len=*), parameter :: myname = 'PrintAtomMatrix'
    character(len=*), intent(in) :: label
    integer,intent(in) :: i
    type(matrixType), intent(in) :: mat
    logical, intent(in) :: isComplex
    type(ioType), intent(in) :: io
    type(atomicxType), intent(in) :: atomic
    integer :: u,v,m,n,k1,k2

    write(io%uout,'(a,a,a,i0)') "==Start ",trim(label),"====== for atom ",i,"========="
    write(io%uout,*) 'Real part'
    m= atomic%basis%norbitals
    n= atomic%species%norbs(atomic%atoms%sp(i))

    write(io%uout,'(4x)',advance="no")
    do u=1,n
        write(io%uout,'(i12,1x)',advance="no") atomic%atoms%orbs(i,u)
    enddo
    write(io%uout,*)
    do u=1,n
        k1=atomic%atoms%orbs(i,u)
        write(io%uout,'(i4,1x)',advance="no") k1
        do v=1,n
          k2=atomic%atoms%orbs(i,v)
          write(io%uout,'(f12.6,1x)',advance="no") real(mat%a(k1,k2))
        enddo
        write(io%uout,'(1x,i0)')u
    enddo
    write(io%uout,'(4x)',advance="no")
    do u=1,n
        write(io%uout,'(i12,1x)',advance="no") u
    enddo
    write(io%uout,*)
    if (isComplex) then
      write(io%uout,*) 'Imaginary part'
      write(io%uout,'(4x)',advance="no")
      do u=1,n
        write(io%uout,'(i12,1x)',advance="no") atomic%atoms%orbs(i,u)
      enddo
      write(io%uout,*)
      do u=1,n
        k1=atomic%atoms%orbs(i,u)
        write(io%uout,'(i4,1x)',advance="no") k1
        do v=1,n
          k2=atomic%atoms%orbs(i,v)
          write(io%uout,'(f12.6,1x)',advance="no") aimag(mat%a(k1,k2))
        enddo
        write(io%uout,'(1x,i0)')u
      enddo
      write(io%uout,'(4x)',advance="no")
      do u=1,n
        write(io%uout,'(i12,1x)',advance="no") u
      enddo
      write(io%uout,*)
    endif
    write(io%uout,*) "==End ",trim(label),"======"
    end subroutine PrintAtomMatrix

!> \brief prints to the unit the atoms in xyz format
!> \author Alin M Elena
!> \date 15/11/07, 10:15:42
!> \param unit output unit
!> \param atomic type(atomicxType) contains all info about the atoms and basis set and some parameters
!> \param lIsVelocity logical if set to true prints the velocities also
!> \param label characters if present will be printed in the comment line
  subroutine PrintXYZ(unit,atomic,lIsVelocity,label)
    character(len=*), parameter :: sMyName="PrintXYZ"
    type(atomicxType), intent(in) :: atomic
    logical, intent(in) :: lIsVelocity
    character(len=*), intent(in),optional :: label
    integer, intent(in) :: unit
    integer :: i
    write(unit,'(i0)')atomic%atoms%natoms
    if (present(label)) then
      write(unit,'(a)')trim(label)
    else
      write(unit,'(a)')""
    endif
    do i=1,atomic%atoms%natoms
    if (lIsVelocity) then
      write(unit,'(a2,1x,6f16.8)') symbol(atomic%atoms%sp(i)),atomic%atoms%x(i),atomic%atoms%y(i),atomic%atoms%z(i),&
        atomic%atoms%vx(i),atomic%atoms%vy(i),atomic%atoms%vz(i)
    else
      write(unit,'(a2,1x,3f16.8)')  &
        symbol(atomic%species%z(atomic%atoms%sp(i))),atomic%atoms%x(i),atomic%atoms%y(i),atomic%atoms%z(i)
    endif
    enddo
  end subroutine PrintXYZ
end module m_Gutenberg