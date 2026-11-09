subroutine equations
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Production and depletion equations
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-07-09   A.J. Koning    A     Original code
!    2     2026-11-09   A.J. Koning    A     Two-pass solution: all end-of-beam inventories before cooling
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
!              rtyp, &                            ! type of beta decay, beta - : 1 , beta + : 2 (from ENDF format)
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
  integer   :: ia      ! mass number from abundance table
  integer   :: is      ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: is_p    ! isotope counter: -1=total, 0=ground state 1=isomer
  integer   :: isob    ! counter
  integer   :: it      ! counter for tritons
  integer   :: iz      ! charge number of residual nucleus
  integer   :: Ncool   ! number of cooling time steps
  integer   :: Zparent ! Z of parent isotope
  real(sgl) :: dT      ! time step
  real(dbl) :: N_p     ! number of isotopes for parent at end of irradiation
  integer   :: ipass   ! pass: 1 = irradiation (all nuclides), 2 = cooling (all nuclides)
  integer   :: it1     ! first time point of current pass
  integer   :: it2     ! last time point of current pass
  logical   :: parent_ok ! parent inventory at end of irradiation is available
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
! The nuclides are processed from high to low Z. For beta- decay the parent (Z-1) is therefore only reached after its daughter.
! The cooling inventory of a daughter needs the parent inventory at the end of irradiation, so the calculation is done in two
! passes over all nuclides:
!   pass 1: irradiation (time points 1..Ntime) for every nuclide, which fixes Niso(.,.,.,Ntime) for all of them
!   pass 2: cooling (time points Ntime+1..numtime) for every nuclide, using the end-of-beam parent inventories from pass 1
! Passes added after error report of Sam Ogunnubi
!
  R_T = dble(reaction_rate(0, 0, -1))
  do ipass = 1, 2
    if (ipass == 1) then
      it1 = 1
      it2 = Ntime
    else
      it1 = Ntime + 1
      it2 = numtime
    endif
    do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
      if (ipass == 1) then
        do it=0,numtime
          Nisotot(iz,it)=0.
        enddo
      endif
      do ia = Acomp, Acomp - Adepth, -1
        do is = -1, Nisomer(iz, ia)
          if (ipass == 1) Yexist(iz,ia,is)=.false.
          if (iz /= Zcomp + 1 .and. .not. rpexist(iz, ia, is)) cycle
          if (ipass == 1) then
            Tmax(iz,ia,is)=0.
            Tmaxactivity(iz,ia,is)=0.
            do it=1,5
              Tp(iz,ia,is,it)=0.
            enddo
            if (iso == 1) Ynatexist(iz,ia,is)=.false.
            do it=0,numtime
              Niso(iz,ia,is,it)=0.
              Nisorel(iz,ia,is,it)=0.
            enddo
            if (iz == Ztarget .and. ia == Atarget .and. is ==  -1) then
              Niso(iz, ia, is, 0) = dble(N_0)
              Nisotot(iz, 0) = Nisotot(iz, 0) + Niso(iz, ia, is, 0)
            endif
          endif
          R_Ti = dble(reaction_rate(iz, ia, is))
          lambda_i = dble(lambda(iz, ia, is))
          enum_i = dble(N_0) * R_Ti
          denom_i = lambda_i - R_T
          do it = it1, it2
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
                eps_decay = 1.d-12 * max(abs(lambda_i), abs(R_T), 1.d-30)
                if (abs(denom_i) > eps_decay) then
                  term_Ti = (exp_T - exp_i) / denom_i
                else
                  term_Ti = TT * exp_i
                endif
                Niso(iz, ia, is, it) = enum_i * term_Ti
              else
!
! Cooling: decay of the inventory present at the end of irradiation
!
                exp_c = exp( -lambda_i * (TT - Tir))
                Niso(iz, ia, is, it) = Niso(iz, ia, is, Ntime) * exp_c
              endif
!
! 2. Production from decay of other isotope
!
              if (flagdecay) then
                do isob = -1, 1, 2
                  Zparent = iz + isob
                  if (Zparent < 0) cycle
!
! Account for isomeric decay of parent
!
                  do is_p = -1, Nisomer(Zparent, ia)
                    if ((isob ==  -1 .and. rtyp(Zparent, ia, is_p) == 1) .or. &
 &                    (isob == 1 .and. rtyp(Zparent, ia, is_p) == 2)) then
                      lambda_p = dble(lambda(Zparent, ia, is_p))
                      denom_ip = lambda_i - lambda_p
                      if (it <= Ntime) then
                        R_Tp = dble(reaction_rate(Zparent, ia, is_p))
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
                          Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + term
                        endif
                      else
!
! Cooling only: ingrowth from the parent inventory at the end of irradiation (computed in pass 1).
! The parent is only used if it was part of this calculation, otherwise Niso could hold data from a previous call.
!
                        parent_ok = (Zparent <= Zcomp + 1 .and. Zparent >= Zcomp + 1 - Zdepth)
                        if (parent_ok) parent_ok = Yexist(Zparent, ia, is_p)
                        if (parent_ok) then
                          N_p = dble(Niso(Zparent, ia, is_p, Ntime))
                          exp_cp = exp( -lambda_p * (TT - Tir))
                          exp_ci = exp( -lambda_i * (TT - Tir))
                          tau = TT - Tir
                          eps_decay = 1.d-12 * max(abs(lambda_i), abs(lambda_p), 1.d-30)
                          if (abs(denom_ip) > eps_decay) then
                            term = N_p * lambda_p / denom_ip * (exp_cp - exp_ci)
                          else
                            term = N_p * lambda_p * tau * exp_cp
                          endif
                          Niso(iz, ia, is, it) = Niso(iz, ia, is, it) + term
                        endif
                      endif
                    endif
                  enddo
                enddo
              endif
            endif
            if (Niso(iz, ia, is, it) > 0.) then
              Yexist(iz, ia, is) = .true.
              Ynatexist(iz, ia, is) = .true.
            endif
            Nisotot(iz, it) = Nisotot(iz, it) + Niso(iz, ia, is, it)
          enddo
          if (ipass == 2) cycle
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
    enddo
  enddo
!
! Relative isotope fractions per element, after all time points are known
!
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
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
end subroutine equations
! Copyright A.J. Koning 2026
