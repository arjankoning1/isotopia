subroutine prodyield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculate production yields
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              sgl, &                             ! single precision kind
!              dbl                                ! double precision kind
!              activity, &                        ! activity of produced isotope in MBq
!              Atarget, &                         ! mass number of target nucleus
!              avogadro, &                        ! Avogadro's number
!              daysec, &                          ! number of seconds in a day
!              dbl, selected_real_kind(15, 307) & ! double precision kind
!              hoursec, &                         ! number of seconds in an hour
!              Ibeam, &                           ! beam current in mA for isotope production
!              iso, &                             ! counter for isotope
!              isotope, &                         ! isotope of natural element
!              lambda, &                          ! decay rate per isotope
!              minutesec, &                       ! number of seconds in a minute
!              Niso, &                            ! number of isotopes produced after irradiation
!              Nisorel, &                         ! fraction of number of produced isotopes per ele
!              Nisotot, &                         ! number of elemental isotopes produced after irr
!              Ntar0, &                           ! number of original target atoms
!              Ntime, &                           ! number of time points
!              nuc, &                             ! symbol of nucleus
!              numtime, &                         ! number of time points
!              prate, &                           ! production rate per isotope
!              radiounit, &                       ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              rhotarget, &                       ! target material density
!              rtyp, &                            ! type of beta decay, beta - : 1 , beta + : 2 (from ENDF format)
!              Tco, &                             ! cooling time per unit
!              Tgrid, &                           ! time
!              Tir, &                             ! irradiation time per unit
!              Tmax, &                            ! irradiation time with maximal yield
!              Tmaxactivity, &                    ! time of maximum activity of produced isoto
!              Tp, &                              ! irradiation time with maximal yield per time unit
!              Vtar, &                            ! active target volume
!              yearsec, &                         ! number of seconds in a year
!              yield, &                           ! yield of produced isotope in MBq / (mA.h)
!              yieldunit, &                       ! unit for isotope yield: num (number), mug, mg, g, or kg
!              Ztarget                            ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  integer   :: ia      ! mass number from abundance table
  integer   :: is      ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: isob    ! counter
  integer   :: it      ! counter for tritons
  integer   :: itm     ! maxmum time
  integer   :: iz      ! charge number of residual nucleus
  integer   :: Ncool   ! number of cooling time steps
  integer   :: Zparent ! Z of parent isotope
  real(sgl) :: dT      ! time step
  real(sgl) :: acmax   ! maximum activity
  real(sgl) :: N0      ! number of isotopes
  real(sgl) :: rfac    ! conversion factor for radioactivity
  real(sgl) :: yfac    ! conversion factor for isotope yield
  real(dbl) :: C1      ! constant
  real(dbl) :: CP1     ! constant
  real(dbl) :: denom   ! help variable
  real(dbl) :: denomD  ! help variable
  real(dbl) :: exp1    ! exponent
  real(dbl) :: exp2    ! exponent
  real(dbl) :: expo1   ! exponent
  real(dbl) :: expo2   ! exponent
  real(dbl) :: lamD    ! decay rate per isotope
  real(dbl) :: lamPD   ! decay rate for parent isotope
  real(dbl) :: prate0  ! production rate for all isotopes
  real(dbl) :: pratei  ! production rate per isotope
  real(dbl) :: prateP  ! production rate
  real(dbl) :: t1      ! help variable
  real(dbl) :: t2      ! help variable
  real(dbl) :: term    ! help variable
  real(dbl) :: TT      ! help variable
!
! ************ Initial condition for irradiation ***********************
!
  Ntar0 = avogadro / Atarget * rhotarget * Vtar
!
! ******************** Set time grid ***********************************
!
  Ntime = numtime / 2
  dT = Tir / Ntime
  do it = 0, Ntime
    Tgrid(it) = it * dT
  enddo
  if (Tco > 0.) then
    Ncool = numtime - Ntime
    dT = Tco / Ncool
    do it = Ntime + 1, numtime
      Tgrid(it) = Tir + (it - Ntime) * dT
    enddo
  endif
!
! ******************** Activity yield in MBq ***************************
!
  prate0 = dble(prate(0, 0, -1))
  do iz = Zcomp, Zcomp - Zdepth, -1
    do it=0,numtime
      Nisotot(iz,it)=0.
      Nisototnat(iz,it)=0.
    enddo
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, 1
        Yexist(iz,ia,is)=.false.
        if ( .not. rpexist(iz, ia, is)) cycle
        Tmax(iz,ia,is)=0.
        Tmaxactivity(iz,ia,is)=0.
        do it=1,5
          Tp(iz,ia,is,it)=0.
        enddo
        do it=0,numtime
          Niso(iz,ia,is,it)=0.
          activity(iz,ia,is,it)=0.
          yield(iz,ia,is,it)=0.
          Nisorel(iz,ia,is,it)=0.
          Nisonat(iz,ia,is,it)=0.
          activitynat(iz,ia,is,it)=0.
          yieldnat(iz,ia,is,it)=0.
          Nisorelnat(iz,ia,is,it)=0.
        enddo
        if (iz == Ztarget .and. ia == Atarget .and. is ==  -1) then
          Niso(iz, ia, is, 0) = dble(Ntar0)
          Nisotot(iz, 0) = Nisotot(iz, 0) + Niso(iz, ia, is, 0)
        endif
        pratei = dble(prate(iz, ia, is))
        lamD = dble(lambda(iz, ia, is))
        C1 = dble(Ntar0) * pratei
        denomD = lamD - prate0
        do it = 1, numtime
          TT = dble(Tgrid(it))
!
! Depletion of target
!
          if (iz == Ztarget .and. ia == Atarget .and. is ==  -1) then
            if (it <= Ntime) then
              Niso(iz, ia, is, it) = dble(Ntar0) * exp( -prate0 * TT)
            else
              Niso(iz, ia, is, it) = Niso(iz, ia, is, Ntime)
            endif
          else
!
! Production and decay of other isotopes
!
! 1. Production directly from target
!
            if (it <= Ntime) then
              t1 = prate0 * TT
              expo1 = exp( - t1)
              t2 = lamD * TT
              expo2 = exp( - t2)
              if (denomD /= 0) then
                exp1 = expo1 / denomD - expo2 / denomD
                Niso(iz, ia, is, it) = C1 * exp1
              endif
            else
              t2 = lamD * (TT - Tir)
              expo2 = exp( - t2)
              term = Niso(iz, ia, is, Ntime) * expo2
              Niso(iz, ia, is, it) = term
            endif
!
! 2. Production from decay of other isotope
!
            if (flagdecay) then
              do isob = -1, 1
                Zparent = iz + isob
                if (Zparent < 0) cycle
                if ((isob ==  -1 .and. rtyp(Zparent, ia, -1) == 1) .or. (isob == 1 .and. rtyp(Zparent, ia, -1) == 2)) then
                  lamPD = dble(lambda(Zparent, ia, is))
                  if (it <= Ntime) then
                    prateP = dble(prate(Zparent, ia, is))
                    if (prateP > 0..and.lamPD > 0.) then
                      CP1 = dble(Ntar0) * prateP
                      t1 = lamPD * TT
                      expo1 = exp( -t1)
                      denom = lamD - lamPD
                      exp2 = expo1 / denom - expo2 / denom
                      term = lamPD * CP1 / (lamPD - prate0) * (exp1 - exp2)
                      Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + term
                    endif
                  else
!
! Cooling only
!
                    N0 = Niso(Zparent, ia, is, Ntime)
                    t1 = lamPD * (TT - Tir)
                    exp1 = exp( -t1)
                    t2 = lamD * (TT - Tir)
                    exp2 = exp( -t2)
                    denom = lamD - lamPD
                    if (denom /= 0.) then
                      term = N0 * lamPD / denom * (exp1 - exp2)
                      Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + term
                    endif
                  endif
                endif
              enddo
            endif
            activity(iz, ia, is, it) = lamD * Niso(iz, ia, is, it) * 1.e-6
            if (it <= Ntime) yield(iz, ia, is, it) = max( (activity(iz, ia, is, it) - activity(iz, ia, is, it - 1)) / &
 &            (Ibeam * dble(Tgrid(it) - Tgrid(it - 1))), 0.)
          endif
          if (Niso(iz, ia, is, it) > 0.) Yexist(iz, ia, is) = .true.
          Nisotot(iz, it) = Nisotot(iz, it) + Niso(iz, ia, is, it)
        enddo
        if ( .not. Yexist(iz, ia, is)) cycle
        if (lamD > 0..and.prate0 > 0.) then
          Tmax(iz, ia, is) = log(lamD / prate0) / (lamD - prate0)
        else
          Tmax(iz, ia, is) = 1.e30
        endif
!
! Write irradiation time with maximum yield in years, days, etc.
!
        TT = Tmax(iz, ia, is)
        Tp(iz, ia, is, 1) = int(TT / yearsec)
        TT = TT - Tp(iz, ia, is, 1) * yearsec
        Tp(iz, ia, is, 2) = int(TT / daysec)
        TT = TT - Tp(iz, ia, is, 2) * daysec
        Tp(iz, ia, is, 3) = int(TT / hoursec)
        TT = TT - Tp(iz, ia, is, 3) * hoursec
        Tp(iz, ia, is, 4) = int(TT / minutesec)
        TT = TT - Tp(iz, ia, is, 4) * minutesec
        Tp(iz, ia, is, 5) = int(TT)
      enddo
    enddo
Loop1: do ia = Acomp, Acomp - Adepth, -1
      do is = -1, 1
        if ( .not. Yexist(iz, ia, is)) cycle Loop1
        do it = 0, numtime
          if (Nisotot(iz, it) /= 0.) Nisorel(iz, ia, is, it) = Niso(iz, ia, is, it) / Nisotot(iz, it)
        enddo
      enddo
    enddo Loop1
  enddo
!
! Transform quantities to user-dependent units
!
  rfac = 1.
  yfac = 1.
  if (radiounit == 'bq') rfac = 1.e6
  if (radiounit == 'kbq') rfac = 1.e3
  if (radiounit == 'gbq') rfac = 1.e-3
  if (radiounit == 'ci') rfac = 1./3.7e4
  if (radiounit == 'kci') rfac = 1./3.7e7
  if (radiounit == 'mci') rfac = 1./3.7e1
  do iz = Zcomp, Zcomp - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      if (yieldunit == 'g') yfac = real(ia)/avogadro
      if (yieldunit == 'mug') yfac = real(ia)/avogadro*1.e6
      if (yieldunit == 'mg') yfac = real(ia)/avogadro*1.e3
      if (yieldunit == 'kg') yfac = real(ia)/avogadro*1.e-3
      do is = -1, Nisomer(iz, ia)
        acmax = 0.
        do it = 1, numtime
          activity(iz, ia, is, it) = rfac * activity(iz, ia, is, it)
          yield(iz, ia, is, it) = rfac * yield(iz, ia, is, it)
          Niso(iz, ia, is, it) = yfac * Niso(iz, ia, is, it)
!
! MV correct the case is=-1  when an isomer is present.
! The condition on tmax is to guarantee that if the isomer activity goes to zero, the code is still executed.
!
          itm = int(Tmaxactivity(iz,ia,-1))
          if (is.eq.1 .and. Niso(iz,ia,1,itm) > 0.) then
            Niso(iz,ia,-1,it) = Niso(iz,ia,0,it) + Niso(iz,ia,1,it)
            activity(iz,ia,-1,it) = activity(iz,ia,0,it) + activity(iz,ia,1,it)
          endif
          if (is == 1 .and. activity(iz,ia,-1,it) > activity(iz,ia,-1,itm)) Tmaxactivity(iz,ia,-1) = it
! MV end
          if (activity(iz, ia, is, it) > acmax) then
            Tmaxactivity(iz, ia, is) = it
            acmax = activity(iz, ia, is, it)
          endif
        enddo
      enddo
    enddo
  enddo
  return
end subroutine prodyield
! Copyright A.J. Koning 2023
