subroutine isotopiainput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Input
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
! readinput   : subroutine to read input lines
! input       : subroutine to read input for keywords
! checkkeyword: subroutine to check for errors in keywords
! checkvalue  : subroutine to check for errors in values
! checkfiles  : subroutine to check for errors in external data files
! inputout    : subroutine to write input parameters
!
  call readinput
  call input
  call checkkeyword
  call checkvalue
  call checkfiles
  call inputout
  return
end subroutine isotopiainput
! Copyright A.J. Koning 2023
