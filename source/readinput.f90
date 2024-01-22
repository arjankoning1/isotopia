subroutine readinput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read user input
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              inline, &   ! input line
!              iso, &      ! counter for isotope
!              nlines, &   ! number of input lines
!              numlines    ! number of input lines
!
! *** Declaration of local data
!
  implicit none
  integer :: i              ! counter
  integer :: istat          ! help variable
!
! ************************** User Input ********************************
!
! We read the complete input file first as a set of character strings.
! The actual keywords will be read from these later on.
!
  if (iso /= 1) return
!
! If command line parameters are present
!
  nlines = command_argument_count()
  if (nlines > 0) then
!
! loop over options
!
   do i=1,nlines
     call get_command_argument(i,inline(i))
   enddo
  else
!
! If command line parameters are NOT present. Standard input file, e.g. isotopia < isotopia.inp
!
    i = 1
    do
      read(*, '(a132)', iostat = istat) inline(i)
      if (istat == -1) exit
      i = i + 1
      if (i > numlines) then
        write(*, '(" ISOTOPIA-error: Number of input lines exceeds", i5)') numlines
        write(*, '(" numlines in isotopia.cmb should be increased")')
        stop
      endif
    enddo
    nlines = i - 1
  endif
!
! ************** Convert uppercase to lowercase characters *************
!
! For easy handling of all the input parameters, the whole input is converted to lowercase characters, 
! with the exception of filenames or other character strings.
!
! convert: subroutine to convert input line from upper case to lowercase
!
  do i = 1, nlines
    call convert(i)
  enddo
  return
end subroutine readinput
! Copyright A.J. Koning 2023
