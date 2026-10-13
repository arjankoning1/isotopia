subroutine equations
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Production and depletion equations
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-07-09   A.J. Koning    A     Original code
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
!              Ibeam, &                           ! beam current for isotope production
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
!              reaction_rate, &                   ! reaction rate per isotope
!              radiounit, &                       ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              ndk, &                             ! number of decay modes (from ENDF format)
!              dkZ, &                             ! Z of daughter per decay mode (-1: no daughter, e.g. SF)
!              dkA, &                             ! A of daughter per decay mode
!              dkis, &                            ! isomeric state of daughter per decay mode (RFS)
!              dkbr, &                            ! branching ratio per decay mode
!              rpexist, &                         ! logical for residual production c.s.
!              Tco, &                             ! cooling time per unit
!              Tgrid, &                           ! time
!              Tir, &                             ! irradiation time per unit
!              Tmax, &                            ! irradiation time with maximal yield
!              Tmaxactivity, &                    ! time of maximum activity of produced isoto
!              Tp, &                              ! irradiation time with maximal yield per time unit
!              V_target, &                            ! active target volume
!              yearsec, &                         ! number of seconds in a year
!              yield, &                           ! yield of produced isotope in MBq / (mA.h)
!              Ztarget                            ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: maxfeed=5*5*3*numdk ! maximum number of feeding decay modes per nuclide
  integer   :: Ap      ! A of parent isotope
  integer   :: ia      ! mass number from abundance table
  integer   :: ip      ! counter for feeding decay mode
  integer   :: is      ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: is_p    ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: is_p1   ! first isomeric state of parent
  integer   :: is_p2   ! last isomeric state of parent
  integer   :: it      ! counter for tritons
  integer   :: iz      ! charge number of residual nucleus
  integer   :: k       ! counter for decay mode
  integer   :: Ncool   ! number of cooling time steps
  integer   :: Nfeed   ! number of feeding decay modes
  integer   :: Zp      ! Z of parent isotope
  integer, dimension(maxfeed) :: Afeed  ! A of parent per feeding decay mode
  integer, dimension(maxfeed) :: isfeed ! isomeric state of parent per feeding decay mode
  integer, dimension(maxfeed) :: Zfeed  ! Z of parent per feeding decay mode
  real(dbl), dimension(maxfeed) :: brfeed ! branching ratio per feeding decay mode
  real(sgl) :: dT      ! time step
  real(sgl) :: N_p     ! number of isotopes for parent
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
  real(dbl) :: br      ! branching ratio
  real(dbl) :: TT      ! help variable
  real(dbl) :: tau     ! help variable
  real(dbl) :: eps_decay ! help variable
!
! ************ Initial condition for irradiation ***********************
!
  N_0 = avogadro / Atarget * M_target
  if (iso == 1) then
    N_0_nat = abun(1) * N_0
  else
    N_0_nat = N_0_nat + abun(iso) * N_0
  endif
!
! ******************** Set time grid ***********************************
!
! Ntime is number of time points until end of irration, i.e. before cooling
! time points between Ntime and numtime are for cooling only
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
  R_T = dble(reaction_rate(0, 0, -1))
!
! The equations are solved in three passes over all nuclides:
! 1. irradiation: feeding only depends on the direct production rate of the parent, so the order of nuclides is irrelevant
! 2. cooling: feeding depends on the parent inventory at the end of irradiation, which is complete after pass 1
! 3. totals, existence flags and time of maximum yield
!
! ******************** Pass 1: irradiation *****************************
!
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
    do it=0,numtime
      Nisotot(iz,it)=0.
    enddo
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        Yexist(iz,ia,is)=.false.
        do it=0,numtime
          Niso(iz,ia,is,it)=0.
          Nisorel(iz,ia,is,it)=0.
        enddo
        call feeding
        if (iz /= Zcomp + 1 .and. .not. rpexist(iz, ia, is) .and. Nfeed == 0) cycle
        Tmax(iz,ia,is)=0.
        Tmaxactivity(iz,ia,is)=0.
        do it=1,5
          Tp(iz,ia,is,it)=0.
        enddo
        if (iso == 1) Ynatexist(iz,ia,is)=.false.
        if (iz == Ztarget .and. ia == Atarget .and. is ==  -1) then
          Niso(iz, ia, is, 0) = dble(N_0)
          Nisotot(iz, 0) = Nisotot(iz, 0) + Niso(iz, ia, is, 0)
        endif
        R_Ti = dble(reaction_rate(iz, ia, is))
        lambda_i = dble(lambda(iz, ia, is))
        enum_i = dble(N_0) * R_Ti
        denom_i = lambda_i - R_T
        do it = 1, Ntime
          TT = dble(Tgrid(it))
!
! Depletion of target
!
          if (iz == Ztarget .and. ia == Atarget .and. is ==  -1) then
            Niso(iz, ia, is, it) = dble(N_0) * exp( -R_T * TT)
          else
!
! Production and decay of other isotopes
!
! 1. Production directly from target
!
            exp_T = exp( -R_T * TT)
            exp_i = exp( -lambda_i * TT)
            eps_decay = 1.d-12 * max(abs(lambda_i), abs(R_T), 1.d-30)
            if (abs(denom_i) > eps_decay) then
              term_Ti = (exp_T - exp_i) / denom_i
            else
              term_Ti = TT * exp_i
            endif
            Niso(iz, ia, is, it) = enum_i * term_Ti
!
! 2. Production from decay of other isotopes, weighted with the branching ratio of each decay mode
!
            do ip = 1, Nfeed
              Zp = Zfeed(ip)
              Ap = Afeed(ip)
              is_p = isfeed(ip)
              br = brfeed(ip)
              lambda_p = dble(lambda(Zp, Ap, is_p))
              denom_ip = lambda_i - lambda_p
              R_Tp = dble(reaction_rate(Zp, Ap, is_p))
              if (R_Tp > 0. .and. lambda_p > 0.) then
                exp_p = exp( -lambda_p * TT)
                eps_decay = 1.d-12 * max(abs(lambda_i), abs(lambda_p), 1.d-30)
                if (abs(denom_ip) > eps_decay) then
                  term_pi = (exp_p - exp_i) / denom_ip
                else
                  term_pi = TT * exp_p
                endif
                enum_p = dble(N_0) * R_Tp
                denom_p = lambda_p - R_T
                if (abs(denom_p) > eps_decay) then
                  term = lambda_p * enum_p / denom_p * (term_Ti - term_pi)
                else
                  if (abs(denom_i) > eps_decay) then
                    term = R_T * enum_p * (TT * exp_T * denom_i - exp_T + exp_i) / (denom_i * denom_i)
                  else
                    term = 0.5d0 * R_T * enum_p * TT * TT * exp_T
                  endif
                endif
                Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + br * term
              endif
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
!
! ******************** Pass 2: cooling *********************************
!
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        call feeding
        if (iz /= Zcomp + 1 .and. .not. rpexist(iz, ia, is) .and. Nfeed == 0) cycle
        lambda_i = dble(lambda(iz, ia, is))
        do it = Ntime + 1, numtime
          TT = dble(Tgrid(it))
          if (iz == Ztarget .and. ia == Atarget .and. is ==  -1) then
            Niso(iz, ia, is, it) = Niso(iz, ia, is, Ntime)
          else
            exp_c = exp( -lambda_i * (TT - Tir))
            Niso(iz, ia, is, it) = Niso(iz, ia, is, Ntime) * exp_c
            do ip = 1, Nfeed
              Zp = Zfeed(ip)
              Ap = Afeed(ip)
              is_p = isfeed(ip)
              br = brfeed(ip)
              lambda_p = dble(lambda(Zp, Ap, is_p))
              denom_ip = lambda_i - lambda_p
              N_p = Niso(Zp, Ap, is_p, Ntime)
              exp_cp = exp( -lambda_p * (TT - Tir))
              exp_ci = exp( -lambda_i * (TT - Tir))
              tau = TT - Tir
              eps_decay = 1.d-12 * max(abs(lambda_i), abs(lambda_p), 1.d-30)
              if (abs(denom_ip) > eps_decay) then
                term = N_p * lambda_p / denom_ip * (exp_cp - exp_ci)
              else
                term = N_p * lambda_p * tau * exp_cp
              endif
              Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + br * term
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
!
! ******************** Pass 3: totals and time of maximum yield ********
!
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        call feeding
        if (iz /= Zcomp + 1 .and. .not. rpexist(iz, ia, is) .and. Nfeed == 0) cycle
        do it = 1, numtime
          if (Niso(iz, ia, is, it) > 0.) then
            Yexist(iz, ia, is) = .true.
            Ynatexist(iz, ia, is) = .true.
          endif
          Nisotot(iz, it) = Nisotot(iz, it) + Niso(iz, ia, is, it)
        enddo
        if ( .not. Yexist(iz, ia, is)) cycle
        lambda_i = dble(lambda(iz, ia, is))
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
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        if (Yexist(iz, ia, is)) then
          do it = 0, numtime
            if (Nisotot(iz, it) /= 0.) Nisorel(iz, ia, is, it) = Niso(iz, ia, is, it) / Nisotot(iz, it)
          enddo
        endif
      enddo
    enddo
  enddo
  return
contains
!
! ******************** Feeding decay modes *****************************
!
  subroutine feeding
!
! Collect all parent decay modes that end in nuclide (iz, ia, is). Parents lie within Z-2..Z+2 and A..A+4 of the daughter.
! For a parent with isomeric states, the states 0..Nisomer are used, not the total, to avoid double counting.
! The total of the daughter (is=-1) is fed by all decay modes, an isomeric state only by modes with that final state (RFS).
! Internal transition within the same nuclide only feeds a different isomeric state, not the total.
! The target is only depleted, it is not fed by decay.
!
    Nfeed = 0
    if (.not. flagdecay) return
    if (iz == Ztarget .and. ia == Atarget .and. is == -1) return
    do Zp = max(iz - 2, Zcomp + 1 - Zdepth, 0), min(iz + 2, Zcomp + 1)
      do Ap = ia, min(ia + 4, Acomp)
        if (Nisomer(Zp, Ap) >= 0) then
          is_p1 = 0
          is_p2 = Nisomer(Zp, Ap)
        else
          is_p1 = -1
          is_p2 = -1
        endif
        do is_p = is_p1, is_p2
          if (Zp == iz .and. Ap == ia .and. (is == -1 .or. is_p == is)) cycle
          do k = 1, ndk(Zp, Ap, is_p)
            if (dkZ(Zp, Ap, is_p, k) /= iz .or. dkA(Zp, Ap, is_p, k) /= ia) cycle
            if (is /= -1 .and. dkis(Zp, Ap, is_p, k) /= is) cycle
            if (dkbr(Zp, Ap, is_p, k) <= 0.) cycle
            Nfeed = Nfeed + 1
            Zfeed(Nfeed) = Zp
            Afeed(Nfeed) = Ap
            isfeed(Nfeed) = is_p
            brfeed(Nfeed) = dble(dkbr(Zp, Ap, is_p, k))
          enddo
        enddo
      enddo
    enddo
  end subroutine feeding
end subroutine equations
! Copyright A.J. Koning 2026
