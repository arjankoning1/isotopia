subroutine input
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read user input
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
  use A1_error_handling_mod
!
!              activity, &    ! activity of produced isotope in MBq
!              Area, &        ! target area in cm^2
!              Atarget, &     ! mass number of target nucleus
!              crosspath, &   ! directory containing cross sections
!              Eback, &       ! lower end of energy range in MeV for isotope
!              Ebeam, &       ! incident energy in MeV for isotope production
!              flagnatural, & ! flag for calculation of natural element
!              Ibeam, &       ! beam current in mA for isotope production
!              inline, &      ! input line
!              iso, &         ! counter for isotope
!              nlines, &      ! number of input lines
!              nuc, &         ! symbol of nucleus
!              numA, &        ! number of masses
!              numZ, &        ! number of elements
!              path, &        ! directory containing files to be read
!              ptype0, &      ! type of incident particle
!              radiounit, &   ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              rhotarget, &   ! target material density
!              Starget, &     ! symbol of target nucleus
!              Tco, &         ! cooling time per unit
!              Tcool, &       ! cooling time per unit cooling time unit (y, d, h, m, s)
!              Tir, &         ! irradiation time per unit
!              Tirrad, &      ! irradiation time per unit irradiation time unit
!              unitTcool, &   ! cooling time unit (y, d, h, m, s)
!              unitTirrad, &  ! irradiation time unit (y, d, h, m, s)
!              yield, &       ! yield of produced isotope in MBq / (mA.h)
!              yieldunit, &   ! unit for isotope yield: num (number), mug, mg, g, or kg
!              Ztarget        ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch          ! character
  character(len=132) :: line        ! line
  character(len=132) :: key         ! keyword
  character(len=132) :: val         ! value or string
  character(len=132) :: word(40)    ! words on input line
  integer            :: i           ! counter
  integer            :: ix
  integer            :: i2          ! value
  integer            :: ia          ! mass number from abundance table
  integer            :: is          ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: istat       ! logical for file access
  integer            :: iz          ! charge number of residual nucleus
  integer            :: k           ! designator for particle
  integer            :: lental      ! length of string
!
! ************ Read first set of variables from input lines ************
!
! 1. Initializations
!
  flagnatural = .false.
  ptype0 = ' '
  Ztarget = 0
  Starget = '  '
  Atarget = -1
  Zdepth = 10
  Adepth = 20
  Ebeam = -1.
  Eback = -1.
  Ibeam = 1.
  radiounit = 'gbq'
  yieldunit = 'num'
  Area = 1.
  Tirrad = 0
  unitTirrad = ' '
  Tcool = 0
  unitTcool = ' '
  Tirrad(1) = 1
  unitTirrad(1) = 'd'
  Tcool(1) = 1
  unitTcool(1) = 'd'
  rhotarget = -1.
  flagdecay = .true.
  flagZAoutput = .false.
  flagcross = .false.
  xsfile = ' '
  source ='ISOTOPIA'
  oformat ='YANDF-0.4'
!
! 2. Read input
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified and the corresponding values are read.
! Erroneous input is immediately checked. The keywords and number of values on each line are retrieved from the input.
!
Loop1:  do i = 1, nlines
    line = inline(i)
    call getkeywords(inline(i), word)
    key = word(1)
    val = word(2)
    ch = word(2)(1:1)
    if (key == 'projectile') then
      ptype0 = ch
      cycle
    endif
    if (key == 'element') then
      if (ch >= '0' .and. ch <= '9') then
        read(val, * , iostat = istat) Ztarget
        if (istat /= 0) call read_error(line, istat)
        if (Ztarget < 1 .or. Ztarget > numZ)  call read_error(line, istat)
        Starget = nuc(Ztarget)
        cycle
      else
        read(val, '(a2)', iostat = istat) Starget
        if (istat /= 0) call read_error(line, istat)
        Starget(1:1) = char(ichar(Starget(1:1)) - 32)
        cycle
      endif
    endif
    if (key == 'mass') then
      read(val, * , iostat = istat) Atarget
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'zdepth') then
      read(val, * , iostat = istat) Zdepth
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'adepth') then
      read(val, * , iostat = istat) Adepth
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ebeam') then
      read(val, * , iostat = istat) Ebeam
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ibeam') then
      read(val, * , iostat = istat) Ibeam
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'eback') then
      read(val, * , iostat = istat) Eback
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'area') then
      read(val, * , iostat = istat) area
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'radiounit') then
      read(val, * , iostat = istat) radiounit
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'yieldunit') then
      read(val, * , iostat = istat) yieldunit
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'zaoutput') then
      if (ch == 'n') flagZAoutput = .false.
      if (ch == 'y') flagZAoutput = .true.
      if (ch /= 'y' .and. ch /= 'n') cycle
      cycle
    endif
    if (key == 'outcross') then
      if (ch == 'n') flagcross = .false.
      if (ch == 'y') flagcross = .true.
      if (ch /= 'y' .and. ch /= 'n') cycle
      cycle
    endif
    if (key == 'tirrad') then
      Tirrad(1) = 0
      unitTirrad(1) = ' '
      do k = 1, 5
        read(word(2*k), '(i9)', iostat = istat) Tirrad(k)
        if (istat < 0) call read_error(line, istat)
        if (istat > 0) cycle Loop1
        read(word(2*k+1), '(a1)', iostat = istat) unitTirrad(k)
        if (istat /= 0) call read_error(line, istat)
      enddo
      cycle
    endif
    if (key == 'rho') then
      read(val, * , iostat = istat) rhotarget
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'tcool') then
      Tcool(1) = 0
      unitTcool(1) = ' '
      do k = 1, 5
        read(word(2 * k), '(i9)' , iostat = istat) Tcool(k)
        if (istat < 0) call read_error(line, istat)
        if (istat > 0) cycle Loop1
        read(word(2 * k + 1), '(a1)' , iostat = istat) unitTcool(k)
        if (istat /= 0) call read_error(line, istat)
      enddo
      cycle
    endif
    if (key == 'decay') then
      if (ch == 'n') flagdecay = .false.
      if (ch == 'y') flagdecay = .true.
      if (ch /= 'y' .and. ch /= 'n') cycle
      cycle
    endif
    if (key == 'crosspath') then
      lental = 0
      do i2 = 11, 132
        ch = inline(i)(i2:i2)
        if (ch /= ' ') then
          lental = lental + 1
          crosspath(lental:lental) = ch
        endif
      enddo
      cycle
    endif
    if (key == 'xsfile') then
      read(word(2), * , iostat = istat) iz
      if (istat /= 0) call read_error(line, istat)
      read(word(3), * , iostat = istat) ia
      if (istat /= 0) call read_error(line, istat)
      if (iz < 0 .or. iz > numZ .or. ia < 0 .or. ia > numA) then
        write(*, '(" ISOTOPIA-error: Z, A index out of range ", a)') trim(inline(i))
        stop
      endif
      is = -1
      read(word(5), * , iostat = istat) is
      if (istat < 0) call read_error(line, istat)
      if (istat > 0) cycle Loop1
      if (is <  -1 .or. is > 1)  call read_error(line, istat)
      xsfile(iz, ia, is) = word(4)
      cycle
    endif
    if (key == 'source') then
      ix=index(line,'source')+7
      source=trim(adjustl(line(ix:132)))
      cycle
    endif
    if (key == 'user') then
      ix=index(line,'user')+5
      user=trim(adjustl(line(ix:132)))
      cycle
    endif
    if (key == 'format') then
      ix=index(line,'format')+7
      oformat=trim(adjustl(line(ix:132)))
      cycle
    endif
  enddo Loop1
  return
  stop
end subroutine input
! Copyright A.J. Koning 2023
