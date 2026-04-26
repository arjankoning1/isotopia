subroutine conversion
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Conversion factors 
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-04-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
! *** Declaration of local data
!
  implicit none
! 
! Normalization and strings for output
!
  rstr = 'Bq'
  ystr = 'g'
  cstr = 'mA'
  tstr = 'h'
  rfac = 1.
  mfac = 1.
  tfac = 1.
  cfac = 1.
  if (radiounit == 'bq') then
    rstr = 'Bq'
    rfac = 1.
  endif
  if (radiounit == 'kbq') then
    rstr = 'KBq'
    rfac = 1.e-3
  endif
  if (radiounit == 'mbq') then
    rstr = 'MBq'
    rfac = 1.e-6
  endif
  if (radiounit == 'gbq') then
    rstr = 'GBq'
    rfac = 1.e-9
  endif
  if (radiounit == 'ci') then
    rstr = 'Ci'
    rfac = 1./3.7e10
  endif
  if (radiounit == 'kci') then
    rstr = 'KCi'
    rfac = 1./3.7e13
  endif
  if (radiounit == 'mci') then
    rstr = 'mCi'
    rfac = 1./3.7e7
  endif
  if (yieldunit == 'g') then
    ystr = 'g'
    mfac = 1.
  endif
  if (yieldunit == 'mg') then
    ystr = 'mg'
    mfac = 1.e3
  endif
  if (yieldunit == 'mug') then
    ystr = 'mug'
    mfac = 1.e6
  endif
  if (yieldunit == 'kg') then
    ystr = 'kg'
    mfac = 1.e-3
  endif
  if (timeunit == 's') then
    tstr = 's'
    tfac = 1.
  endif
  if (timeunit == 'h') then
    tstr = 'h'
    tfac = 3600.
  endif
  if (timeunit == 'd') then
    tstr = 'd'
    tfac = 86400.
  endif
  if (currentunit == 'ma') then
    cstr = 'mA'
    cfac = 1.e-3
  endif
  if (currentunit == 'mua') then
    cstr = 'muA'
    cfac = 1.e-6
  endif
  if (currentunit == 'a') then
    cstr = 'A'
    cfac = 1.
  endif
  Ibeam = Ibeam_input * cfac
  if (targetmass_input /= -1.) targetmass = targetmass_input / mfac
  return
end subroutine conversion
! Copyright A.J. Koning 2026
