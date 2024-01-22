subroutine production
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculate isotope production
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
! prodrates: subroutine to calculate reaction rates
! prodyield: subroutine to calculate production yields
! prodout  : subroutine for output
!
  call prodrates
  call prodyield
  call prodout
  return
end subroutine production
! Copyright A.J. Koning 2023
