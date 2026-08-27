subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Machine dependent statements
!
! Revision    Date      Author      Quality  Description
! =====================================================
!    1     2024-12-19   A.J. Koning    A     Original code
!    2     2026-08-27   A.J. Koning    A     Runtime definition of ISOTOPIA directory and user
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_isotopia_mod
!
  implicit none
  logical             :: lexist
  character(len=1024) :: code_dir
  character(len=1024) :: base_dir
  character(len=1024) :: isotopia_dir
  character(len=1024) :: isotopia_user
  integer             :: envstat
  integer             :: i
  integer             :: n
  integer             :: year
  integer             :: month
  integer             :: day
  integer             :: values(8)
!
! ************************ Set directories *****************************
!
  call get_environment_variable('ISOTOPIA_DIR', isotopia_dir, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    code_dir = trim(isotopia_dir)
  else
    code_dir = '/path/to/isotopia/'
  endif
!
! Remove a trailing slash, if present, and determine the parent directory.
! By default isotopia.libs is a sibling directory of ISOTOPIA.
!
  i = len_trim(code_dir)
  if (i > 1) then
    if (code_dir(i:i) == '/') code_dir = code_dir(:i - 1)
  endif
  i = scan(trim(code_dir), '/', back=.true.)
  if (i > 0) then
    base_dir = code_dir(:i)
  else
    base_dir = './'
  endif
!
  path = trim(code_dir)//'/'
  crosspath = trim(base_dir)//'isotopia.libs/'
  libname = 'iaea.2026'
!
! Test accessibility of the internal ISOTOPIA data files.
!
  inquire (file=trim(path)//'files/abundance/z001', exist=lexist)
  if (.not. lexist) then
    write(*, '(a)') 'ISOTOPIA error: data files not found.'
    write(*, '(2a)') 'Expected file: ', trim(path)//'files/abundance/z001'
    write(*, '(a)') 'Set the ISOTOPIA_DIR environment variable:'
    write(*, '(a)') '  export ISOTOPIA_DIR=/path/to/isotopia'
    write(*, '(a)') 'Alternatively, edit code_dir in source/machine.f90'
    write(*, '(a)') 'and rebuild ISOTOPIA.'
    error stop 77
  endif
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
!
  call get_environment_variable('ISOTOPIA_USER', isotopia_user, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    user = trim(isotopia_user)
  else
    user = 'Unknown User'
  endif
  return
end subroutine machine
! Copyright A.J. Koning 2026
