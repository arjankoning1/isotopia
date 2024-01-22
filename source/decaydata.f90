subroutine decaydata
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Decay data (half lives)
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
!              daysec, &    ! number of seconds in a day
!              hoursec, &   ! number of seconds in an hour
!              iso, &       ! counter for isotope
!              isonum, &    ! number of isotopes in element
!              isotope, &   ! isotope of natural element
!              k0, &        ! index of incident particle
!              lambda, &    ! decay rate per isotope
!              minutesec, & ! number of seconds in a minute
!              nuc, &       ! symbol of nucleus
!              parA, &      ! mass number of particle
!              path, &      ! directory containing files to be read
!              rtyp, &      ! type of beta decay, beta - : 1 , beta + : 2 (from ENDF format)
!              Td, &        ! half life per time unit
!              Thalf, &     ! half life of nuclide in sec.
!              yearsec      ! number of seconds in a year
!
! *** Declaration of local data
!
  implicit none
  character(len=4)   :: decaystring !
  character(len=132) :: decayfile   ! decay data file
  character(len=80)  :: string      ! line with parameter value
  integer            :: Ahigh       !
  integer            :: i           ! level
  integer            :: istat       ! 
  integer            :: ia          ! mass number from abundance table
  integer            :: iline       ! number of lines
  integer            :: is          ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: iz          ! charge number of residual nucleus
  integer            :: NC          ! number of lines for MF/MT number
  real(sgl)          :: rt          ! radii of projectile and target
  real(sgl)          :: TT          ! help variable
!
! ************** Find half lives for all involved nuclides *************
!
! Read decay constants from JEFF-3.1.1 radioactive decay data library
!
  decaystring = 'z000'
  if (isonum > 1) then
    Ahigh = isotope(isonum) + parA(k0)
  else
    Ahigh = Acomp
  endif
  lambda = 0.
  rtyp = 0
  Thalf = 1.e30
  Td = 0.
  do iz = Zcomp + 1, Zcomp - Zdepth - 1, -1
    write(decaystring(2:4), '(i3.3)') iz
    decayfile = trim(decaypath) //decaystring
    open (unit = 1, status = 'unknown', file = decayfile)
Loop1: do
      read(1, '(a80)', iostat = istat) string
      if (istat == -1) exit Loop1
      if (string(72:80) /= '1451    5') cycle
      read(string(8:10), '(i3)') ia
      if (ia < Acomp - Adepth) cycle
      if (ia > Ahigh) exit Loop1
      is = -1
      if (string(11:11) == 'M') is = 1
      if (string(11:11) == 'N') is = 2
      do
        read(1, '(a80)', iostat = istat) string
        if (istat == -1) exit Loop1
        if (string(72:80) == '8457    2') exit
      enddo
      read(string(1:11), '(e11.5)') Thalf(iz, ia, is)
      if (flagdecay) then
        read(string(45:55), '(i11)') NC
        iline = 1 + (NC - 1) / 6
        do i = 1, iline
          read(1, '(a80)', iostat = istat) string
          if (istat == -1) exit Loop1
        enddo
        read(1, '(a80)', iostat = istat) string
        if (istat == -1) exit Loop1
        read(1, '(a80)', iostat = istat) string
        if (istat == -1) exit Loop1
        read(string(1:11), '(e11.5)') rt
        rtyp(iz, ia, is) = int(rt)
      endif
      if (Thalf(iz, ia, is) == 0.) Thalf(iz, ia, is) = 1.e30
!
! Write half life in years, days, etc.
!
      TT = Thalf(iz, ia, is)
      Td(iz, ia, is, 1) = int(TT / yearsec)
      TT = TT - Td(iz, ia, is, 1) * yearsec
      Td(iz, ia, is, 2) = int(TT / daysec)
      TT = TT - Td(iz, ia, is, 2) * daysec
      Td(iz, ia, is, 3) = int(TT / hoursec)
      TT = TT - Td(iz, ia, is, 3) * hoursec
      Td(iz, ia, is, 4) = int(TT / minutesec)
      TT = TT - Td(iz, ia, is, 4) * minutesec
      Td(iz, ia, is, 5) = int(TT)
!
! Calculate decay rates
! Nuclides with a lifetime longer 1.e17 sec are considered stable
!
      if (Thalf(iz, ia, is) > 1.e17) then
        lambda(iz, ia, is) = 0.
      else
        lambda(iz, ia, is) = log(2.) / Thalf(iz, ia, is)
      endif
      if (is == -1) then
        rtyp(iz, ia, 0) = rtyp(iz, ia, -1)
        lambda(iz, ia, 0) = lambda(iz, ia, -1)
        Thalf(iz, ia, 0) = Thalf(iz, ia, -1)
        do i = 1, 5
          Td(iz, ia, 0, i) = Td(iz, ia, -1, i)
        enddo
      endif
    enddo Loop1
    close (unit = 1)
    do ia = Ahigh, Acomp - Adepth, -1
      if (Thalf(iz, ia, 1) > 1.e17) then
        rtyp(iz, ia, 1) = rtyp(iz, ia, -1)
        lambda(iz, ia, 1) = lambda(iz, ia, -1)
        Thalf(iz, ia, 1) = Thalf(iz, ia, -1)
        do i = 1, 5
          Td(iz, ia, 1, i) = Td(iz, ia, -1, i)
        enddo
      endif
    enddo
  enddo
  return
end subroutine decaydata
! Copyright A.J. Koning 2023
