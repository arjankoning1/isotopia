subroutine production
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculate isotope production
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-01-04   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
! reactionrates: subroutine to calculate reaction rates
! equation: subroutine for production and depletion equations
! activities: subroutine to calculate activities
! prodout  : subroutine for output
!
  call reactionrates
  call equations
  call activities
  call prodout
  return
end subroutine production
! Copyright A.J. Koning 2026
