subroutine prodrates
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculate reaction rates
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-03-09   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              sgl, &       ! single precision kind
!              Area, &      ! target area in cm^2
!              Eback, &     ! lower end of energy range in MeV for isotope
!              Ebeam, &     ! incident energy in MeV for isotope production
!              heat, &      ! produced heat
!              Ibeam, &     ! beam current in mA for isotope production
!              iso, &       ! counter for isotope
!              isotope, &   ! isotope of natural element
!              k0, &        ! index of incident particle
!              Mtar, &      ! active target mass
!              Nenrp, &     ! number of incident energies for residual prod
!              nuc, &       ! symbol of nucleus
!              numen, &     ! number of energies
!              parZ, &      ! charge number of particle
!              prate, &     ! production rate per isotope
!              projnum, &   ! number of incident particles [s^ -1]
!              qelem, &     ! elementary charge in C
!              rhotarget, & ! target material density
!              targetdx, &  ! effective thickness of target
!              Vtar, &      ! active target volume
!              xsrp         ! residual production cross section in mb
!
! *** Declaration of local data
!
  implicit none
  character(len=3)   :: Estring          ! string
  character(len=132) :: string           ! line
  character(len=132) :: specfile        ! 
  integer, parameter :: numint=10000     ! number of integration points
  integer            :: ia               ! mass number from abundance table
  integer            :: is               ! isotope counter: -1=total, 0=ground state 1=isomer
  integer            :: iz               ! charge number of residual nucleus
  integer            :: N                ! neutron number of residual nucleus
  integer            :: nE               ! number of energies
  integer            :: nen              ! energy counter
  integer            :: iE
  integer            :: istat
  integer            :: NintE            ! number of integration points
  real(sgl)          :: dE               ! help variable
  real(sgl)          :: dEdx             ! stopping power
  real(sgl)          :: dxdE(numint)     ! 1/stopping power
  real(sgl)          :: phi(numint)      ! spectrum
  real(sgl)          :: dxdEsum          ! integration sum
  real(sgl)          :: E                ! incident energy
  real(sgl)          :: Ea               ! start energy of local adjustment
  real(sgl)          :: Eb               ! end energy of local adjustment
  real(sgl)          :: Egrid(0:numen)   !
  real(sgl)          :: dum
  real(sgl)          :: Eint(0:numint)   ! energy on integration grid
  real(sgl)          :: dEint(numint)    ! delta energy on integration grid
  real(sgl)          :: ratesum          ! integration sum
  real(sgl)          :: xs               ! help variable
  real(sgl)          :: xsa              ! help variable
  real(sgl)          :: xsb              ! help variable
!
! **** Photons: read Bremsstrahlung spectrum for incident electron energy
!
  Mtar = targetmass
  dEint = 0.
  Eint = 0.
  phi = 0.
  if (k0 == 0) then
    iE = min(200, int(Ebeam))
    iE = max(iE, 1)
    Estring = '000'
    write(Estring(1:3), '(i3.3)') iE
    specfile = trim(path)//'files/photon_spectra/brems_'//Estring//'MeV.txt'
    open (unit = 1, status = 'unknown', file = specfile)
    nen = 0
    do
      read(1, '(a)', iostat = istat) string
      if (istat == -1) exit
      if (string(1:1) == '#') cycle
      nen = nen + 1
      read(string, *) Eint(nen), phi(nen)
    enddo
    close(1)
    NintE = nen
    do nen = 2, nintE
      dEint(nen-1) = Eint(nen) - Eint(nen-1)
    enddo
    dEint(NintE) = dEint(NintE-1)
  endif
!
! **** Neutrons: read neutron spectrum 
!
  if (k0 == 1) then
    iE = min(200, int(Ebeam))
    iE = max(iE, 1)
    Estring = '000'
    write(Estring(1:3), '(i3.3)') iE
    specfile = trim(path)//'files/neutron_spectra/n-hfr.dat'
    open (unit = 1, status = 'unknown', file = specfile)
    nen = 0
    do
      read(1, '(a)', iostat = istat) string
      if (istat == -1) exit
      if (string(1:1) == '#') cycle
      nen = nen + 1
      read(string, *) iE, Ea, Eb, dum,  phi(nen)
      Eint(nen) = 0.5 * (Ea + Eb) * 1.e-6
      dEint(nen) = (Ea - Eb) * 1.e-6
    enddo
    close(1)
    NintE = nen
  endif
!
! **** Charged particles: Determine integration grid and stopping power *
!
! stoppingpower: subroutine to calculate stopping power
!
  if (k0 > 1) then
    NintE = 100
    dE = (Ebeam - Eback) / NintE
    dxdEsum = 0.
    do nE = 1, NintE
      Eint(nE) = Eback + (nE - 0.5) * dE
      call stoppingpower(Eint(nE), dEdx)
      dxdE(nE) = 1. / dEdx
      dxdEsum = dxdEsum + dxdE(nE)
    enddo
    targetdx = dxdEsum * dE
    Vtar = Area * targetdx
    Mtar = rhotarget * Vtar
    heat = Ibeam * (Ebeam - Eback)
    projnum = Ibeam / (1000. * parZ(k0) * qelem)
  endif
!
! ********************* Calculate reaction rates ***********************
!
! locate   : subroutine to find value in ordered table
! pol1     : subroutine for interpolation of first order
!
  Egrid = 0.
  call prodres(0, 0, -1)
  do iz = Zcomp + 1, 0, -1
    do ia = Acomp, 0, -1
      do is = -1, 1
        prate(iz, ia, is) = 0.
        if (iz == 0 .and. ia == 0 .and. is ==  -1) then
          Nenrp = Nennon
          xsrp = xsnon
          goto 140
        endif
        if (iz < Zcomp - Zdepth .or. ia < Acomp - Adepth) cycle
140     ratesum = 0.
        call prodres(iz, ia, is)
        if ( .not. rpexist(iz, ia, is)) cycle
        N = Nenrp
        do nen = 1, N
          Egrid(nen) = Ein(nen)
        enddo
        do nE = 1, NintE
          E = Eint(nE)
          if (E < Egrid(1)) cycle
          if (E > Egrid(N)) exit
          call locate(Egrid, 1, N, E, nen)
          if (nen == 0) cycle
          Ea = Egrid(nen)
          Eb = Egrid(nen + 1)
          xsa = xsrp(nen)
          xsb = xsrp(nen + 1)
          call pol1(Ea, Eb, xsa, xsb, E, xs)
          if (k0 == 0) ratesum = ratesum + phi(nE) * xs * dEint(nE)
          if (k0 == 1) ratesum = ratesum + phi(nE) * xs
          if (k0 > 1) ratesum = ratesum + dxdE(nE) * xs
        enddo
        if (k0 <= 1) then
          prate(iz, ia, is) = fluxtotal  * ratesum * 1.e-27
        else
          prate(iz, ia, is) = projnum / Vtar * ratesum * dE * 1.e-27
        endif
      enddo
    enddo
  enddo
  return
end subroutine prodrates
! Copyright A.J. Koning 2026
