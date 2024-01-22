subroutine convert(i)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Convert input line from upper case to lowercase
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
! use A0_isotopia_mod, only: & ! General module with all global variables
!              abun, &      ! natural abundance
!              crosspath, & ! directory containing cross sections
!              inline, &    ! input line
!              path         ! directory containing files to be read
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: str    ! input line
  integer            :: i      ! level
  integer            :: k      ! designator for particle
!
! ************** Convert uppercase to lowercase characters *************
!
! For easy handling of all the input parameters, the whole input is
! converted to lowercase characters, with the exception of filenames or
! other character strings.
!
  str = inline(i)
  do k = 1, 132
    if (inline(i)(k:k) >= 'A' .and. inline(i)(k:k) <= 'Z') inline(i)(k:k) = char(ichar(inline(i)(k:k)) + 32)
  enddo
  do k = 0, 60
    if (inline(i)(k+1:k+6) == 'xsfile') then
      inline(i)(k + 7:132) = str(k + 7:132)
      return
    endif
    if (inline(i)(k+1:k+9) == 'abundance') then
      inline(i)(k + 10:132) = str(k + 10:132)
      return
    endif
    if (inline(i)(k+1:k+9) == 'crosspath') then
      inline(i)(k + 10:132) = str(k + 10:132)
      return
    endif
  enddo
  return
end subroutine convert
! Copyright A.J. Koning 2023
