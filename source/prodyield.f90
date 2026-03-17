subroutine prodyield
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculate production yields
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-02-26   A.J. Koning    A     Original code
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
!              N_0, &                             ! number of original target atoms
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
  integer   :: is_p    ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: isob    ! counter
  integer   :: it      ! counter for tritons
  integer   :: itm     ! maxmum time
  integer   :: iz      ! charge number of residual nucleus
  integer   :: Ncool   ! number of cooling time steps
  integer   :: Zparent ! Z of parent isotope
  real(sgl) :: dT      ! time step
  real(sgl) :: acmax   ! maximum activity
  real(sgl) :: N_p     ! number of isotopes for parent
  real(sgl) :: rfac    ! conversion factor for radioactivity
  real(sgl) :: yfac    ! conversion factor for isotope yield
  real(dbl) :: enum_i  ! constant
  real(dbl) :: enum_p  ! constant
  real(dbl) :: denom_i ! help variable
  real(dbl) :: denom_p ! help variable
  real(dbl) :: denom_ip! help variable
  real(dbl) :: exp_cp  ! exponent
  real(dbl) :: exp_ci  ! exponent
  real(dbl) :: exp_i   ! exponent
  real(dbl) :: exp_c   ! exponent
  real(dbl) :: exp_p   ! exponent
  real(dbl) :: exp_T   ! exponent
  real(dbl) :: lambda_i! decay rate per isotope
  real(dbl) :: lambda_p! decay rate for parent isotope
  real(dbl) :: R_T     ! production rate for target
  real(dbl) :: R_Ti    ! production rate per isotope
  real(dbl) :: R_Tp    ! production rate for parent
  real(dbl) :: term    ! help variable
  real(dbl) :: term_pi ! help variable
  real(dbl) :: term_Ti ! help variable
  real(dbl) :: TT      ! help variable
!
! ************ Initial condition for irradiation ***********************
!
  N_0 = avogadro / Atarget * Mtar
!
! ******************** Set time grid ***********************************
!
! Ntime is number of time points until end of irration, i.e. before cooling
! time points bewteen Ntime and numtime are for cooling only
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
! ******************** Activity ****************************************
!
  R_T = dble(prate(0, 0, -1))
  do iz = Zcomp, Zcomp - Zdepth, -1
    do it=0,numtime
      Nisotot(iz,it)=0.
      Nisototnat(iz,it)=0.
    enddo
    do ia = Acomp, Acomp - Adepth, -1
!print*,iz,ia
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
          Niso(iz, ia, is, 0) = dble(N_0)
          Nisotot(iz, 0) = Nisotot(iz, 0) + Niso(iz, ia, is, 0)
        endif
        R_Ti = dble(prate(iz, ia, is))
        lambda_i = dble(lambda(iz, ia, is))
        enum_i = dble(N_0) * R_Ti
        denom_i = lambda_i - R_T
        do it = 1, numtime
          TT = dble(Tgrid(it))
!
! Depletion of target
!
          if (iz == Ztarget .and. ia == Atarget .and. is ==  -1) then
            if (it <= Ntime) then
              Niso(iz, ia, is, it) = dble(N_0) * exp( -R_T * TT)
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
              exp_T = exp( -R_T * TT)
              exp_i = exp( -lambda_i * TT)
              if (denom_i /= 0) then
                term_Ti = exp_T / denom_i - exp_i / denom_i
                Niso(iz, ia, is, it) = enum_i * term_Ti
!if (iz == 53 .and. ia == 122 .and. is == -1) print*," A ", TT, Niso(iz, ia, is, it)
              endif
            else
!
! Cooling
!
              exp_c = exp( -lambda_i * (TT - Tir))
              Niso(iz, ia, is, it) = Niso(iz, ia, is, Ntime) * exp_c
!if (iz == 53 .and. ia == 122 .and. is == -1) print*," B ", TT, Niso(iz, ia, is, it)
            endif
!
! 2. Production from decay of other isotope
!
            if (flagdecay) then
              do isob = -1, 1
                Zparent = iz + isob
                if (Zparent < 0) cycle
!
! Account for isomeric decay of parent
!
                do is_p = -1, -1
                  if ((isob ==  -1 .and. rtyp(Zparent, ia, is_p) == 1) .or. (isob == 1 .and. rtyp(Zparent, ia, is_p) == 2)) then
                    lambda_p = dble(lambda(Zparent, ia, is_p))
                    denom_ip = lambda_i - lambda_p
                    if (it <= Ntime) then
                      R_Tp = dble(prate(Zparent, ia, is_p))
                      if (R_Tp > 0. .and. lambda_p > 0.) then
                        enum_p = dble(N_0) * R_Tp
                        denom_p = lambda_p - R_T
                        exp_p = exp( -lambda_p * TT)
                        term_pi = exp_p / denom_ip - exp_i / denom_ip
                        term = lambda_p * enum_p / denom_p * (term_Ti - term_pi)
                        Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + term
!if (iz == 53 .and. ia == 122 .and. is == -1) print*," C ", TT, Niso(iz, ia, is, it),term
                      endif
                    else
!
! Cooling only
!
                      N_p = Niso(Zparent, ia, is_p, Ntime)
                      exp_cp = exp( -lambda_p * (TT - Tir))
                      exp_ci = exp( -lambda_i * (TT - Tir))
                      if (denom_ip /= 0.) then
                        term = N_p * lambda_p / denom_ip * (exp_cp - exp_ci)
                        Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + term
!if (iz == 53 .and. ia == 122 .and. is == -1) print*," D ", TT, Niso(iz, ia, is, it),term,Tir
                      endif
                    endif
                  endif
                enddo
              enddo
            endif
            activity(iz, ia, is, it) = lambda_i * Niso(iz, ia, is, it) * 1.e-6
            if (it <= Ntime) yield(iz, ia, is, it) = max( (activity(iz, ia, is, it) - activity(iz, ia, is, it - 1)) / &
 &            (Ibeam * dble(Tgrid(it) - Tgrid(it - 1))), 0.)
          endif
          if (Niso(iz, ia, is, it) > 0.) Yexist(iz, ia, is) = .true.
          Nisotot(iz, it) = Nisotot(iz, it) + Niso(iz, ia, is, it)
        enddo
        if ( .not. Yexist(iz, ia, is)) cycle
        if (lambda_i > 0..and. R_T > 0.) then
          Tmax(iz, ia, is) = log(lambda_i / R_T) / (lambda_i - R_T)
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
