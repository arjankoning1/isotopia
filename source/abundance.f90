subroutine abundance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Natural abundances
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              sgl, &       ! single precision kind
!              abun, &      ! natural abundance
!              inline, &    ! input line
!              iso, &       ! counter for isotope
!              isonum, &    ! number of isotopes in element
!              isotope, &   ! isotope of natural element
!              natstring, & ! string extension for file names
!              nlines, &    ! number of input lines
!              nuc, &       ! symbol of nucleus
!              numiso, &    ! maximum number of isotopes per element
!              path, &      ! directory containing files to be read
!              Ztarget      ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist    ! logical to determine existence
  character(len=1)   :: ch        ! character
  character(len=4)   :: abchar    ! help variable
  character(len=132) :: abfile    ! isotopic abundance file
  integer            :: i         ! counter
  integer            :: istat     ! level
  integer            :: i2        ! value
  integer            :: ia        ! mass number from abundance table
  real(sgl)          :: ab        ! isotopic abundance
  real(sgl)          :: abtot     ! summed abundances for normalization
!      INTEGER*4 getcwd, status
!      character*64 dirname
!
! ****************************** Abundances ****************************
!
! Note that for non-natural elements we take the longest-lived isotope
! as default.
!
! 1. Isotopic abundances from user file
!
!      status = getcwd(dirname)
!      write(*,*) 'start'
!      write(*,*) dirname
  do i = 1, nlines
    if (inline(i)(1:9) == 'abundance') then
      do i2 = 11, 132
        ch = inline(i)(i2:i2)
        if (ch /= ' ') then
          read(inline(i)(i2:132), * , iostat = istat) abfile
          if (istat > 0) then
            write(*,'(" ISOTOPIA-error: No natural isotopes for this element, the mass keyword must be different from 0")')
            stop
          endif
          goto 100
        endif
      enddo
    endif
  enddo
!
! 2. Isotopic abundances from abundance directory
!
  abchar = 'z   '
  write(abchar(2:4), '(i3.3)') Ztarget
  abfile = trim(path)//'files/abundance/'//abchar
!
! Read abundances from file
!
100 inquire (file = abfile, exist = lexist)
  if ( .not. lexist) then
    write(*, '(" ISOTOPIA-error: No natural isotopes for this element, the mass keyword must be different from 0")')
    stop
  endif
  open (unit = 2, status = 'old', file = abfile)
  i = 1
  do
    read(2, '(4x, i4, f11.6)', iostat = istat) ia, ab
    if (istat > 0) then
      write(*, '(" ISOTOPIA-error: Format error in abundance file ", a)') trim(abfile)
      stop
    endif
    if (istat == -1) exit
    if (ia == 0) cycle
    isotope(i) = ia
    abun(i) = 0.01 * ab
    i = i + 1
    if (i > numiso) then
      write(*, '(" ISOTOPIA-error: Too many isotopes in abundance file", a, " Increase numiso in A0_isotopia_mod.f90")') &
 &      trim(abfile)
      stop
    endif
  enddo
  close (unit = 2)
  isonum = i - 1
!
! Normalize abundances to 1.
!
  abtot = 0.
  do i = 1, isonum
    abtot = abtot + abun(i)
  enddo
  do i = 1, isonum
    abun(i) = abun(i) / abtot
  enddo
  write(*, '(/" Calculation for multi-isotope case"/)')
  write(*, '("  Isotope Abundance"/)')
  do i = 1, isonum
    write(*, '(2x, i3, a2, f11.6)') isotope(i), nuc(Ztarget), abun(i)
  enddo
!
! Create file name extensions
!
  do i = 1, isonum
    write(natstring(i)(1:4), '(".", i3.3)') isotope(i)
  enddo
  return
end subroutine abundance
! Copyright A.J. Koning 2023
