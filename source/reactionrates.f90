subroutine reactionrates
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculate reaction rates
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-04-06   A.J. Koning    A     Original code
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
!              M_target, &      ! active target mass
!              Nenrp, &     ! number of incident energies for residual prod
!              nuc, &       ! symbol of nucleus
!              numen, &     ! number of energies
!              parZ, &      ! charge number of particle
!              reaction_rate, & ! reaction rate per isotope
!              projnum, &   ! number of incident particles [s^ -1]
!              qelem, &     ! elementary charge in C
!              rho_target, & ! target material density
!              Leff, &  ! effective path length or thickness of target
!              V_target, &      ! active target volume
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
  real(sgl)          :: S                ! stopping power
  real(sgl)          :: phi(numint)      ! spectrum
  real(sgl)          :: phisum           ! integration sum
  real(sgl)          :: E                ! incident energy
  real(sgl)          :: Ea               ! start energy of local adjustment
  real(sgl)          :: Eb               ! end energy of local adjustment
  real(sgl)          :: Egrid(0:numen)   !
  real(sgl)          :: dum
  real(sgl)          :: Eint(0:numint)   ! energy on integration grid
  real(sgl)          :: dEint(numint)    ! delta energy on integration grid
  real(sgl)          :: ratesum          ! integration sum
  real(sgl)          :: number_density   ! number densitym
  real(sgl)          :: xs               ! help variable
  real(sgl)          :: xsa              ! help variable
  real(sgl)          :: xsb              ! help variable
!
! **** Photons: read Bremsstrahlung spectrum for incident electron energy
!
  number_density = 0.
  M_target = targetmass
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
      dEint(nen) = abs(Ea - Eb) * 1.e-6
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
    phisum = 0.
    do nE = 1, NintE
      Eint(nE) = Eback + (nE - 0.5) * dE
      dEint(nE) = dE
      call stoppingpower(Eint(nE), S)
      if (S /= 0.) phi(nE) = 1. / S
      phisum = phisum + phi(nE) * dE
    enddo
    Leff = phisum
    V_target = Area * Leff
    heat = Ibeam * (Ebeam - Eback)
    projnum = Ibeam / (1000. * parZ(k0) * qelem)
    number_density = projnum / V_target
  endif
  if (k0 <= 1) then
    if (targetmass /= -1.) then
      V_target = targetmass / rho_target 
      thickness = V_target / Area
    else
      V_target = Area * thickness
    endif
  endif
  M_target = rho_target * V_target
!
! ********************* Calculate reaction rates ***********************
!
! locate   : subroutine to find value in ordered table
! pol1     : subroutine for interpolation of first order
!
  Egrid = 0.
  call crosssections(0, 0, -1)
  do iz = Zcomp + 1, 0, -1
    do ia = Acomp, 0, -1
      do is = -1, 1
        reaction_rate(iz, ia, is) = 0.
        if ((iz < Zcomp - Zdepth .or. ia < Acomp - Adepth) .and. .not. (iz == 0 .and. ia == 0 .and. is ==  -1)) cycle
        if (iz == 0 .and. ia == 0 .and. is ==  -1) then
          Nenrp = Nennon
          xsrp = xsnon
        endif
        ratesum = 0.
        call crosssections(iz, ia, is)
        if ( .not. rpexist(iz, ia, is)) cycle
        N = Nenrp
        do nen = 1, N
          Egrid(nen) = Ein(nen)
        enddo
        do nE = 1, NintE
          E = Eint(nE)
          if (E < Egrid(1)) cycle
          if (E >= Egrid(N)) cycle
          call locate(Egrid, 1, N, E, nen)
          if (nen == 0) cycle
          Ea = Egrid(nen)
          Eb = Egrid(nen + 1)
          xsa = xsrp(nen)
          xsb = xsrp(nen + 1)
          call pol1(Ea, Eb, xsa, xsb, E, xs)
          if (k0 == 1) then
            ratesum = ratesum + phi(nE) * xs
          else
            ratesum = ratesum + phi(nE) * xs * dEint(nE)
          endif
        enddo
!
! Conversion of mb to cm^2: 1.e-27
!
        if (k0 <= 1) then
          reaction_rate(iz, ia, is) = fluxtotal  * ratesum * 1.e-27
        else
          reaction_rate(iz, ia, is) = number_density * ratesum * 1.e-27
        endif
      enddo
    enddo
  enddo
  return
end subroutine reactionrates
! Copyright A.J. Koning 2026
