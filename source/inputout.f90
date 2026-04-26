subroutine inputout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write input parameters
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-04-26   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              activity, &   ! activity of produced isotope in MBq
!              Area, &       ! target area in cm^2
!              Atarget, &    ! mass number of target nucleus
!              Eback, &      ! lower end of energy range in MeV for isotope
!              Ebeam, &      ! incident energy in MeV for isotope production
!              Ibeam, &      ! beam current for isotope production
!              inline, &     ! input line
!              iso, &        ! counter for isotope
!              isotope, &    ! isotope of natural element
!              nlines, &     ! number of input lines
!              nuc, &        ! symbol of nucleus
!              ptype0, &     ! type of incident particle
!              radiounit, &  ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              rho_target, &  ! target material density
!              Starget, &    ! symbol of target nucleus
!              Tco, &        ! cooling time per unit
!              Tcool, &      ! cooling time per unit cooling time unit (y, d, h, m, s)
!              Tir, &        ! irradiation time per unit
!              Tirrad, &     ! irradiation time per unit irradiation time unit
!              unitTcool, &  ! cooling time unit (y, d, h, m, s)
!              unitTirrad, & ! irradiation time unit (y, d, h, m, s)
!              yield, &      ! yield of produced isotope 
!              yieldunit     ! unit for isotope yield: num (number), mug, mg, g, or kg
!
! *** Declaration of local data
!
  implicit none
  integer :: i              ! level
  character(len=1) :: yesno ! y or n function
!
! *************************** Code and version *************************
!
  write(*, '(/"    ISOTOPIA-2.23   (Version: April 26, 2026)"/)')
  write(*, '(10x, " Simulation of medical isotope production")')
  write(*, '(/" Copyright (C) 2026  A.J. Koning")')
  write(*, '(/" User: ",a)') trim(user)
  write(*, '(" Date: ",a)') trim(date)
!
! ************************** User input file ***************************
!
  write(*, '(/" ########## USER INPUT ##########")')
  write(*, '(/" USER INPUT FILE"/)')
  do i = 1, nlines
    write(*, '(1x, a)') trim(inline(i))
  enddo
!
! ********* All possible input parameters including defaults ***********
!
  write(*, '(/" USER INPUT FILE + DEFAULTS"/)')
  write(*, '(" Keyword           Value   Variable     Explanation"/)')
  write(*, '(" user      ", a, t28, "user         user for this calculation")') trim(user)
  write(*, '(" source    ", a, t28, "source       source for this calculation")') trim(source)
  write(*, '(" format    ", a, t28, "oformat      format for output")') trim(oformat)
  write(*, '(" crosspath ", a, "   crosspath    directory containing cross sections ")') trim(crosspath)
!
! Nuclear reaction
!
  write(*, '(" #"/" # Nuclear reaction"/" #")')
  write(*, '(" projectile          ", a1, "     ptype0       type of incident particle")') ptype0
  write(*, '(" element            ", a2, "     Starget      symbol of target nucleus")') Starget
  write(*, '(" mass              ", i3, "     Atarget      mass number of target nucleus")') Atarget
!
! Accelerator
!
  if (k0 /= 1) then
    write(*, '(" #"/" # Accelerator "/" #")')
    write(*, '(" Ebeam             ", f7.3, " Ebeam        incident energy in MeV")') Ebeam
    if (k0 > 1) write(*, '(" Eback             ", f7.3, " Eback        lower end of energy range in MeV")') Eback
    write(*, '(" Ibeam             ", f7.3, " Ibeam_input  beam current in ",a)') Ibeam_input, trim(cstr)
  endif
!
! Output
!
  write(*, '(" #"/" # Output "/" #")')
  write(*, '(" radiounit             ", a3, " radiounit    unit for radioactivity")') radiounit
  write(*, '(" yieldunit (massunit)  ", a3, " yieldunit    unit for isotope yield")') yieldunit
  write(*, '(" ZAoutput              ", a1, "   flagZAoutput output files per Z, A instead of El, A")') yesno(flagZAoutput)
!
! Target
!
  write(*, '(" #"/" # Target"/" #")')
  do i = 1, 5
    if (Tirrad(i) > 0) write(*, '(" Tirrad      ", i9, "     Tirrad       ", a1, " irradiation time")') Tirrad(i), unitTirrad(i)
  enddo
  write(*, '(" Area              ", f7.3, " Area         target area in cm^2")') Area
  do i = 1, 5
    if (Tcool(i) > 0) write(*, '(" Tcool       ", i9, "     Tcool        ", a1, " cooling time")') Tcool(i), unitTcool(i)
  enddo
  if (rho_target /= -1.) write(*, '(" rho               ", f7.3, " rho_target   target density [g/cm^3] ")') rho_target
  if (targetmass_input /= -1.) write(*, '(" targetmass    ", f7.3, " targetmass_input target mass ")') targetmass_input
  if (k0 <= 1) write(*, '(" fluxtotal ", es15.6, " fluxtotal    total flux")') fluxtotal
  write(*, '(" selfshield            ", a1, " flagselfshield flag for self-shielding ")') yesno(flagselfshield)
!
! Cross section and decay data
!
  write(*, '(" #"/" # Cross section and decay data"/" #")')
  write(*, '(" Zdepth            ", i3, "     Zdepth       depth to which Z numbers are scanned for cross sections")') Zdepth
  write(*, '(" Adepth            ", i3, "     Adepth       depth to which A numbers are scanned for cross sections")') Adepth
  write(*, '(" decay               ", a1, "     flagdecay    flag to include decay of nuclides (n only for diagnosis) ")') &
 &  yesno(flagdecay)
  write(* , *)
  return
end subroutine inputout
! Copyright A.J. Koning 2026
