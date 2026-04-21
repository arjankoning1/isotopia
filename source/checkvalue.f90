subroutine checkvalue
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Check for errors in values
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-04-09   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              abun, &        ! natural abundance
!              Area, &        ! target area in cm^2
!              Atarget, &     ! mass number of target nucleus
!              Eback, &       ! lower end of energy range in MeV for isotope
!              Ebeam, &       ! incident energy in MeV for isotope production
!              flagnatural, & ! flag for calculation of natural element
!              Ibeam, &       ! beam current in mA for isotope production
!              iso, &         ! counter for isotope
!              isotope, &     ! isotope of natural element
!              k0, &          ! index of incident particle
!              nuc, &         ! symbol of nucleus
!              numA, &        ! number of masses
!              numZ, &        ! number of elements
!              parsym, &      ! symbol of particle
!              ptype0, &      ! type of incident particle
!              rho_target, &   ! target material density
!              Starget, &     ! symbol of target nucleus
!              Tco, &         ! cooling time per unit
!              Tcool, &       ! cooling time per unit cooling time unit (y, d, h, m, s)
!              Tir, &         ! irradiation time per unit
!              Tirrad, &      ! irradiation time per unit irradiation time unit
!              unitTcool, &   ! cooling time unit (y, d, h, m, s)
!              unitTirrad, &  ! irradiation time unit (y, d, h, m, s)
!              Ztarget        ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  integer :: iz                ! charge number of residual nucleus
  integer :: k                 ! designator for particle
  integer :: type              ! particle type
!
! All parameters need to fall within certain ranges. These ranges are
! specified in this subroutine and in the manual.
!
! ******************* Check for wrong input variables ******************
!
! While checking, we also set a few parameters on the basis of the input
!
  if (ptype0 == ' ') then
    write(*, '(" ISOTOPIA-error: projectile must be given")')
    stop
  endif
  do type = 0, 6
    if (ptype0 == parsym(type)) then
      k0 = type
      goto 100
    endif
  enddo
  write(*, '(" ISOTOPIA-error: Wrong symbol for projectile: ", a1)') ptype0
  stop
100 if (Ztarget == 0 .and. Starget == '  ') then
    write(*, '(" ISOTOPIA-error: element must be given")')
    stop
  endif
  if (Ztarget == 0) then
    do iz = 1, numZ
      if (nuc(iz) == Starget) then
        Ztarget = iz
        goto 200
      endif
    enddo
    write(*, '(" ISOTOPIA-error: Wrong symbol for element: ", a2)') Starget
    stop
  endif
!
! A calculation for a natural element is specified by target mass 0
!
! abundance  : subroutine for natural abundances
!
200 if (Atarget ==  -1) then
    write(*, '(" ISOTOPIA-error: target mass must be given")')
    stop
  endif
  if (Atarget == 0) then
    flagnatural = .true.
    if (iso == 1) call abundance
    Atarget = isotope(iso)
  endif
  if (Atarget <= 5 .or. Atarget > numA) then
    write(*, '(" ISOTOPIA-error: 5 < Target mass < = ", i3)') numA
    stop
  endif
  if (k0 /= 1) then
    if (Ebeam ==  -1.) then
      write(*, '(" ISOTOPIA-error: incident energy must be given")')
      stop
    endif
    if (Ebeam <= 0..or.Ebeam > 250.) then
      write(*, '(" ISOTOPIA-error: 0 < Ebeam < 250 MeV")')
      stop
    endif
    if (k0 > 1) then
      if (Eback ==  -1.) then
        Eback = max(Ebeam - 5., 0.1)
      else
        if (Eback <= 0..or.Eback > 250.) then
          write(*, '(" ISOTOPIA-error: 0 < Eback < 250 MeV")')
          stop
        endif
      endif
      if (Eback >= Ebeam) then
        write(*, '(" ISOTOPIA-error: Ebeam must be larger than Eback")')
        stop
      endif
    endif
  endif
  Zdepth = min(Zdepth, Ztarget)
  Adepth = min(Adepth, Atarget)
  if (Zdepth < 0 .or. Zdepth > Ztarget) then
    write(*, '(" ISOTOPIA-error: 0 < = Zdepth < = Ztarget")')
    stop
  endif
  if (Adepth < 0 .or. Adepth > Atarget) then
    write(*, '(" ISOTOPIA-error: 0 < = Adepth < = Atarget")')
    stop
  endif
  if (Ibeam <= 0. .or. Ibeam > 10000.) then
    write(*, '(" ISOTOPIA-error: 0 < = Ibeam < 10000 mA")')
    stop
  endif
  if (Area <= 0. .or. Area > 10000.) then
    write(*, '(" ISOTOPIA-error: 0 < = Area < 10000 cm^2")')
    stop
  endif
  if (thickness <= 0. .or. thickness > 10000.) then
    write(*, '(" ISOTOPIA-error: 0 < = thickness < 10000 cm")')
    stop
  endif
  do k = 1, 5
    if (Tirrad(k) < 0 .or. Tirrad(k) >= 1000000) then
      write(*, '(" ISOTOPIA-error: 0 < = Tirrad < 1.e6")')
      stop
    endif
    if (Tcool(k) < 0 .or. Tcool(k) >= 1000000) then
      write(*, '(" ISOTOPIA-error: 0 < = Tcool < 1.e6")')
      stop
    endif
  enddo
  do k = 1, 5
    if (unitTirrad(k) /= ' ' .and. unitTirrad(k) /= 'y' .and. unitTirrad(k) /= 'd' .and. unitTirrad(k) /= 'h' .and. &
      unitTirrad(k) /= 'm' .and. unitTirrad(k) /= 's') then
      write(*, '(" ISOTOPIA-error: wrong unit for Tirrad = ", i9)') Tirrad(k)
      stop
    endif
    if (unitTcool(k) /= ' ' .and. unitTcool(k) /= 'y' .and. unitTcool(k) /= 'd' .and. unitTcool(k) /= 'h' .and. &
      unitTcool(k) /= 'm' .and. unitTcool(k) /= 's') then
      write(*, '(" ISOTOPIA-error: wrong unit for Tcool = ", i9)') Tcool(k)
      stop
    endif
  enddo
  if (rho_target /=  - 1. .and. (rho_target <= 0. .or. rho_target > 100.)) then
    write(*, '(" ISOTOPIA-error: 0 < rhotarget <= 100.")')
    stop
  endif
  if (targetmass /= -1. .and. (targetmass <=  0. .or. targetmass > 1.e9)) then
    write(*, '(" ISOTOPIA-error: 0 < targetmass <= 1.e9")')
    stop
  endif
  if (fluxtotal <=  0. .or. fluxtotal > 1.e38) then
    write(*, '(" ISOTOPIA-error: 0 < fluxtotal <= 1.e38")')
    stop
  endif
  if (fgamma <  0. .or. fgamma > 1.) then
    write(*, '(" ISOTOPIA-error: 0 < fgamma <= 1.")')
    stop
  endif
  if (yieldunit /= 'mug' .and. yieldunit /= 'mg' .and. yieldunit /= 'g' .and. yieldunit /= 'kg') then
    write(*, '(" ISOTOPIA-error: yieldunit must be equal to mug, mg, g or kg")')
    stop
  endif
  if (radiounit /= 'bq' .and. radiounit /= 'kbq' .and. radiounit /= 'mbq' .and. radiounit /= 'gbq' .and. &
 &  radiounit /= 'ci' .and. radiounit /= 'kci' .and. radiounit /= 'mci') then
    write(*, '(" ISOTOPIA-error: radiounit must be equal to bq, kbq, mbq, gbq, ci, kci or mci")')
    stop
  endif
  if (timeunit /= 's' .and. timeunit /= 'h' .and. timeunit /= 'd') then
    write(*, '(" ISOTOPIA-error: timeunit must be equal to s, h or d")')
    stop
  endif
  if (currentunit /= 'mua' .and. currentunit /= 'ma' .and. currentunit /= 'a') then
    write(*, '(" ISOTOPIA-error: currentunit must be equal to mua, ma or a")')
    stop
  endif
  return
end subroutine checkvalue
! Copyright A.J. Koning 2026
