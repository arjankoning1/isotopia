   subroutine prodres(iz, ia, is)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Residual production cross sections
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-12-19   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              sgl, &     ! single precision kind
!              Atarget, & ! mass number of target nucleus
!              Eback, &   ! lower end of energy range in MeV for isotope
!              Ebeam, &   ! incident energy in MeV for isotope production
!              iso, &     ! counter for isotope
!              isotope, & ! isotope of natural element
!              k0, &      ! index of incident particle
!              libname, & ! library name
!              Nenrp, &   ! number of incident energies for residual prod
!              Niso, &    ! number of isotopes produced after irradiation
!              nuc, &     ! symbol of nucleus
!              numen, &   ! number of energies
!              parsym, &  ! symbol of particle
!              path, &    ! directory containing files to be read
!              xsrp, &    ! residual production cross section in mb
!              Ztarget    ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  logical            :: flagiso        !
  logical            :: flagpositive   ! flag for existence of non-zero cross sections
  logical            :: lexist         ! logical to determine existence
  character(len=7)   :: ZAstring       !
  character(len=132) :: nucfile         ! file with cross sections
  character(len=132) :: csfile         ! file with cross sections
  character(len=132) :: pfile          ! parameter file
  character(len=132) :: string         ! line with parameter value
  character(len=16) :: reaction   ! reaction
  character(len=132) :: topline    ! topline
  character(len=16) :: col(2)     ! header
  character(len=16) :: un(2)     ! header
  character(len=80) :: quantity   ! quantity
  integer            :: i              ! counter
  integer            :: Ncol           ! 
  integer            :: istat          ! 
  integer            :: indent
  integer            :: ia             ! mass number from abundance table
  integer            :: iE             ! energy counter
  integer            :: is             ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: iz             ! charge number of residual nucleus
  integer            :: nen            ! energy counter
  real(sgl)          :: E              ! incident energy
  real(sgl)          :: xs             ! help variable
!
! ******************* Read non-elastic cross sections ******************
!
! For convenience in later loops, we store the non-elastic cross section in the 0th element of xsrp. In the next loop, we subtract
! the inelastic cross section from this, i.e. xsrp will contain all non-elastic cross sections other than inelastic.
!
  indent = 0
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  un(2)='mb'
  Ncol=2
  quantity='Cross section'
  Ein = 0.
  xsrp = 0.
  if (iz == 0 .and. ia == 0 .and. is ==  -1) then
    Enon = 0.
    xsnon = 0.
    Nennon = 0
    pfile = trim(xspath)//'xs/'//parsym(k0)//'-'//trim(nuclide)// '-MT005.'//trim(libname)
    open (unit = 1, status = 'unknown', file = trim(pfile))
    iE = 0
    do
      read(1, '(a80)', iostat = istat) string
      if (istat == -1) exit
      if (string(1:1) == '#') cycle
      read(string, * ) E, xs
      iE = iE + 1
      if (iE > numen) then
        write(*, '(" ISOTOPIA-error: too many incident energies: increase numen in A0_isotopia_mod.f90")')
        stop
      endif
      Ein(iE) = E
      xsrp(iE) = xs
      Enon(iE) = E
      xsnon(iE) = xs
    enddo
    close (unit = 1)
    Nenrp = iE
    Nennon = iE
    rpexist(0, 0, -1) = .true.
    if (flagcross) then
      reaction='('//ptype0//',non)'
      csfile = parsym(k0)//'-'//trim(nuclide)//'.non'
      open (unit = 1, status = 'unknown', file = csfile)
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(indent,topline,source,user,date,oformat)
      call write_target(indent)
      call write_reaction(indent,reaction,0.d0,0.d0,3,3)
      call write_quantity(indent,quantity)
      call write_datablock(indent,Ncol,Nennon,col,un)
      do i = 1, Nennon
        write(1, '(2es15.6)') Enon(i), xsnon(i)
      enddo
      close(1)
    endif
    return
  endif
!
! **************** Read residual production cross sections *************
!
  rpexist(iz, ia, is) = .false.
  Nenrp = 0
  flagiso = .false.
  Nisomer(iz, ia) = -1
  ZAstring = '000000 '
  write(ZAstring(1:3), '(i3.3)') iz
  write(ZAstring(4:6), '(i3.3)') ia
  if (is == 0) ZAstring = trim(ZAstring)//'g'
  if (is == 1) ZAstring = trim(ZAstring)//'m'
  Nisomer(iz, ia) = max(Nisomer(iz, ia), is)
!
! Read cross sections and check for positive values in energy range of interest.
!
  reaction='('//ptype0//',x)'
  if (xsfile(iz, ia, is) /= ' ') then
    nucfile = 'file.loc'
    csfile = xsfile(iz, ia, is)
  else
    nucfile = parsym(k0)//'-'// trim(nuclide)//'-rp'//trim(ZAstring)//'.'//trim(libname)
    csfile = trim(xspath)//'residual/'//trim(nucfile)
  endif
  inquire (file = csfile, exist = lexist)
  if (lexist) then
    rpexist(iz, ia, is) = .true.
    if (is >= 0) flagiso = .true.
    flagpositive = .false.
    iE = 0
    open (unit = 1, status = 'unknown', file = csfile)
    if (flagcross) open (unit = 2, status = 'unknown', file = nucfile)
    do
      read(1, '(a80)', iostat = istat) string
      if (istat == -1) exit
      if (flagcross) write(2, '(a)') trim(string)
      if (string(1:1) == '#') cycle
      read(string, * , iostat = istat) E, xs
      if (istat > 0) then
        write(*, '(" ISOTOPIA-error: Problem in cross section file ", a, " around E = ", es12.5, " xs = ", es12.5)') &
 &        trim(csfile), E, xs
        stop
      endif
      if (E >= Eback .and. E <= Ebeam .and. xs > 0.) flagpositive = .true.
      iE = iE + 1
      if (iE > numen) then
        write(*, '(" ISOTOPIA-error: too many incident energies: increase numen in A0_isotopia_mod")')
        stop
      endif
      if (iE > 1) then
        if (E < Ein(iE - 1)) then
          write(*, '(" ISOTOPIA-error: incident energies in cross section file not in increasing order at E = ", es12.5)') E
          stop
        endif
      endif
      Ein(iE) = E
      xsrp(iE) = xs
    enddo
    close (unit = 1)
    if (flagcross) close (unit = 1)
!
! Subtract inelastic cross section from non-elastic cross section
!
    Nenrp = iE
    if ( .not. flagpositive) rpexist(iz, ia, is) = .false.
    if (iz == Ztarget .and. ia == Atarget .and. is <= 0) then
      do i = 1, Nennon
        E = Enon(i)
        if (E < Ein(1)) cycle
        if (E > Ein(Nenrp)) cycle
        call locate(Ein, 1, Nenrp, E, nen)
        if (nen == 0) cycle
        xsnon(i) = xsnon(i) - xsrp(nen)
      enddo
      if (flagcross) then
        csfile = parsym(k0)//'-'//trim(nuclide)//'-non.xs'
        open (unit = 1, status = 'unknown', file = csfile)
        reaction='('//ptype0//',non)'
        topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,0.d0,0.d0,3,3)
        call write_quantity(indent,quantity)
        call write_datablock(indent,Ncol,Nennon,col,un)
        do i = 1, Nennon
          write(1, '(2es15.6)') Enon(i), xsnon(i)
        enddo
        close(1)
      endif
    endif
  endif
  return
end subroutine prodres
! Copyright A.J. Koning 2023
