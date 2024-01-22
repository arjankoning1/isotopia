subroutine checkkeyword
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Check for errors in keywords
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              abun, &      ! natural abundance
!              crosspath, & ! directory containing cross sections
!              inline, &    ! input line
!              nlines, &    ! number of input lines
!              path, &      ! directory containing files to be read
!              radiounit, & ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              yield, &     ! yield of produced isotope in MBq / (mA.h)
!              yieldunit    ! unit for isotope yield: num (number), mug, mg, g, or kg
!
! *** Declaration of local data
!
  implicit none

  integer, parameter :: numkey=24         ! number of keywords
  integer            :: i
  integer            :: j
  character(len=132) :: key               ! keyword
  character(len=132) :: keyword(numkey)   ! keyword
  character(len=132) :: word(40)          ! words on input line
!
! Although it is difficult to prevent the user from all possible input
! errors, we can check for the use of wrong keywords and for unphysical
! values for most of the input variables.
!
! *********************** Check for wrong keywords *********************
!
! ISOTOPIA will stop if a keyword is incorrect
!
  data (keyword(i), i = 1, numkey) / ' ', 'abundance', 'adepth', 'area', 'crosspath', 'decay', &
 &  'eback', 'ebeam', 'element', 'format', 'ibeam', 'mass', 'outcross', 'projectile', 'radiounit', 'rho', 'source', 'tirrad', &
 &  'tcool', 'user', 'xsfile', 'yieldunit', 'zaoutput', 'zdepth'/
!
! A keyword can be de-activated by putting a # in front of it.
! All first words of the input lines are checked against the list of keywords.
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified.
!
Loop1:  do i = 1, nlines
    call getkeywords(inline(i), word)
    key = word(1)
    if (key(1:1) == '#') cycle
    do j = 1, numkey
      if (keyword(j) == key) cycle Loop1
    enddo
    write(*, '(/" ISOTOPIA-error: Wrong keyword: ", a20)') key
    stop
  enddo Loop1
  return
end subroutine checkkeyword
! Copyright A.J. Koning 2023
