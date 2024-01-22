subroutine isotopiainitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Initialization of nuclear reaction info
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              iso        ! counter for isotope
!
! *** Declaration of local data
!
  implicit none
!
! ******************* Initialization ***********************************
!
! reacinitial: subroutine for initialization of nuclear reaction info
! decaydata  : subroutine for decay data (half lives)
!
  call reacinitial
  if (iso == 1) call decaydata
  return
end subroutine isotopiainitial
! Copyright A.J. Koning 2023
