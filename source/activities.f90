subroutine activities
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculate activity
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-04-01   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              sgl, &                             ! single precision kind
!              dbl                                ! double precision kind
!              activity, &                        ! activity of produced isotope in MBq
!              Atarget, &                         ! mass number of target nucleus
!              daysec, &                          ! number of seconds in a day
!              hoursec, &                         ! number of seconds in an hour
!              Niso, &                            ! number of isotopes produced after irradiation
!              Nisorel, &                         ! fraction of number of produced isotopes per ele
!              Nisotot, &                         ! number of elemental isotopes produced after irr
!              N_0, &                             ! number of original target atoms
!              Ntime, &                           ! number of time points
!              nuc, &                             ! symbol of nucleus
!              numtime, &                         ! number of time points
!              radiounit, &                       ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              rtyp, &                            ! type of beta decay, beta - : 1 , beta + : 2 (from ENDF format)
!              Tco, &                             ! cooling time per unit
!              Tmax, &                            ! irradiation time with maximal yield
!              Tmaxactivity, &                    ! time of maximum activity of produced isoto
!              Tp, &                              ! irradiation time with maximal yield per time unit
!              yearsec, &                         ! number of seconds in a year
!              yield, &                           ! yield of produced isotope in MBq / (mA.h)
!              Ztarget                            ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: ia      ! mass number from abundance table
  integer   :: is      ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: it      ! counter for tritons
  integer   :: itm     ! maxmum time
  integer   :: iz      ! charge number of residual nucleus
  real(sgl) :: acmax   ! maximum activity
  real(sgl) :: activity_factor
!
! Calculate activities from decay constants and number of nuclides
!
  if (k0 <= 1) then
    activity_factor = M_target
  else
    activity_factor = Ibeam
  endif
  activity = 0.
  yield = 0.
  specactivity = 0.
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        do it = 1, numtime
          activity(iz, ia, is, it) = dble(lambda(iz, ia, is)) * Niso(iz, ia, is, it)
        enddo
        do it = 1, Ntime
          yield(iz, ia, is, it) = max( (activity(iz, ia, is, it) - activity(iz, ia, is, it - 1)) / &
 &          (activity_factor * dble(Tgrid(it) - Tgrid(it - 1))), 0.)
        enddo
      enddo
    enddo
  enddo
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        acmax = 0.
        do it = 1, numtime
!
! MV correct the case is=-1  when an isomer is present.
! The condition on tmax is to guarantee that if the isomer activity goes to zero, the code is still executed.
!
          itm = int(Tmaxactivity(iz,ia,-1))
          if (is == 1 .and. Niso(iz,ia,1,itm) > 0.) then
            Niso(iz,ia,-1,it) = Niso(iz,ia,0,it) + Niso(iz,ia,1,it)
            activity(iz,ia,-1,it) = activity(iz,ia,0,it) + activity(iz,ia,1,it)
          endif
          if (is == 1 .and. activity(iz,ia,-1,it) > activity(iz,ia,-1,itm)) Tmaxactivity(iz,ia,-1) = it
! MV end
          if (activity(iz, ia, is, it) > acmax) then
            Tmaxactivity(iz, ia, is) = it
            acmax = activity(iz, ia, is, it)
          endif
          specactivity(iz, ia, is, it) = activity(iz, ia, is, it) / M_target
        enddo
      enddo
    enddo
  enddo
  return
end subroutine activities
! Copyright A.J. Koning 2026
