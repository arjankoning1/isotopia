subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Machine dependent statements
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2024-12-19   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              crosspath, & ! directory containing cross sections
!              iso, &       ! counter for isotope
!              isonum, &    ! number of isotopes in element
!              libname, &   ! library name
!              natstring, & ! string extension for file names
!              path         ! directory containing files to be read
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: codedir   ! code directory
  character(len=132) :: basedir   ! base directory
  integer            :: i         ! level
  integer            :: ix        ! level
  integer            :: year      ! year
  integer            :: month     ! month
  integer            :: day       ! day
  integer            :: values(8) ! date and time values
!
! ************************ Set directories *****************************
!
  codedir = '/Users/koning/isotopia/'
  path = trim(codedir)
  ix = index(codedir,'/isotopia/')
  basedir = codedir(1:ix)
!
! ************************ Set nuclear data directory ******************
!
  crosspath = trim(basedir)//'isotopia.libs/'
  libname = 'iaea.2024'
!
! ************************ Set counter for isotope *********************
!
  iso = 1
  do i = 1, numiso
    natstring(i) = '    '
  enddo
  isonum = 1
!
! ************************ Set date and user ***************************
!
  call date_and_time(VALUES=values)
  year=values(1)
  month=values(2)
  day=values(3)
  date='xxxx-xx-xx'
  write(date(1:4),'(i4.4)') year
  write(date(6:7),'(i2.2)') month
  write(date(9:10),'(i2.2)') day
  user = 'Arjan Koning'
  return
end subroutine machine
! Copyright A.J. Koning 2024
