subroutine timer
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Output of execution time
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
  use A0_isotopia_mod
!
!              sgl                  ! single precision kind
!
! *** Declaration of local data
!
  integer   :: hour                 ! number of hours
  integer   :: hundred              ! number of 1/100th of seconds
  integer   :: minute               ! number of minutes
  integer   :: second               ! number of seconds
  real(sgl) :: etime                ! time function
  real(sgl) :: tarray(2)            ! help variable
  real(sgl) :: time                 ! time
!
! ****** Get elapsed time in seconds from beginning of execution *******
!
! etime  : time function
!
! The returned time should be "charge time" (e.g., cp+pp+sys). This
! could be machine dependent.
!
  time = etime(tarray)
  hour = int(time / 3600.)
  minute = int((time - hour * 3600) / 60.)
  second = int(time - hour * 3600 - minute * 60)
  hundred = int(100 * (time - int(time)))
!      write(*,'(/" Execution time:",i3," hours ",i2," minutes ",i2,".",
!     +  i2.2," seconds ")') hour,minute,second,hundred
!      write(*,*)
!      write(*,*) " End of ISOTOPIA isotopic run"
  return
end subroutine timer
! Copyright A.J. Koning 2023
