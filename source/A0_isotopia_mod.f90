module A0_isotopia_mod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: General module with all global variables
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-04-26  A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Definition of single and double precision variables
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
!
!-----------------------------------------------------------------------------------------------------------------------------------
! All global dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: numpar=6         ! number of particles
  integer, parameter :: numZ=124         ! number of elements
  integer, parameter :: numisom=2        ! number of isomers
  integer, parameter :: numiso=20        ! maximum number of isotopes per element
  integer, parameter :: numlines=100     ! number of input lines
  integer, parameter :: numA=339         ! number of masses
  integer, parameter :: numen=1000000    ! number of energies
  integer, parameter :: numtime=100      ! number of time points
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading input lines
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132), dimension(numlines) :: inline  ! input line
  integer                                 :: nlines  ! number of input lines
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=1), dimension(-1:numisom) :: isochar     ! symbol of isomer
  character(len=1), dimension(0:numpar)   :: parsym      ! symbol of particle
  character(len=2), dimension(numZ)       :: nuc         ! symbol of nucleus
  character(len=8), dimension(0:numpar)   :: parname     ! name of particle
  character(len=4), dimension(numiso)     :: natstring   ! string extension for file names
  integer                                 :: iso         ! counter for isotope
  integer, dimension(0:numpar)            :: parA        ! mass number of particle
  integer, dimension(0:numpar)            :: parZ        ! charge number of particle
  real(dbl)                               :: amu         ! atomic mass unit in MeV
  real(sgl)                               :: emass       ! electron mass
  real(sgl)                               :: avogadro    ! Avogadro's number
  real(dbl), dimension(0:numpar)          :: parmass     ! mass of particle in a.m.u.
  real(sgl)                               :: qelem       ! elementary charge in C
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for abundance
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                      :: isonum  ! number of isotopes in element
  integer, dimension(numiso)   :: isotope ! isotope of natural element
  real(sgl), dimension(numiso) :: abun    ! natural abundance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for main input
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical            :: flagnatural                    ! flag for calculation of natural element
  logical            :: flagdecay                      ! flag to include decay from parent nuclide
  logical            :: flagcross                      ! flag for output of cross sections
  logical            :: flagZAoutput                   ! flag for output per Z,A
  logical            :: flagselfshield                 ! flag for self-shielding  
  character(len=1)   :: ptype0                         ! type of incident particle
  character(len=2)   :: Starget                        ! symbol of target nucleus
  character(len=132) :: xsfile(numZ, numA, -1:numisom) ! cross section file
  character(len=132) :: source                         ! source of data
  character(len=132) :: oformat                        ! format of data
  character(len=132) :: user                           ! user of data
  integer            :: Atarget                        ! mass number of target nucleus
  integer            :: k0                             ! index of incident particle
  integer            :: Ztarget                        ! charge number of target nucleus
  integer            :: Adepth                         ! depth of A
  integer            :: Zdepth                         ! depth of Z
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for files (machine)
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=10)  :: date     ! date
  character(len=20)  :: libname   ! library name
  character(len=132) :: crosspath ! directory containing cross sections
  character(len=132) :: path      ! directory containing files to be read
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for checkfiles
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=5)   :: nuclide   ! nuclide
  character(len=6)   :: targetnuclide ! target nuclide
  character(len=132) :: decaypath ! directory with decay data
  character(len=132) :: xspath    ! directory with cross section data
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for prodres
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer                                        :: Nennon     ! number of energies for nonelastic n c.s.
  integer                                        :: Nentot     ! number of energies for total c.s.
  integer, dimension(0:numZ, 0:numA)             :: Nisomer    ! number of energies for residual produci n c.s.
  logical, dimension(0:numZ, 0:numA, -1:numisom) :: rpexist    ! logical for residual production c.s.
  real(sgl), dimension(0:numen)                  :: Ein        ! incident energy
  real(sgl), dimension(0:numen)                  :: Enon       ! incident energy
  real(sgl), dimension(0:numen)                  :: xsnon      ! nonelastic cross section
  real(sgl), dimension(0:numen)                  :: Etot       ! incident energy
  real(sgl), dimension(0:numen)                  :: xstot      ! total cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for medical isotope production
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=1), dimension(5) :: unitTcool  ! cooling time unit (y,d,h,m,s)
  character(len=1), dimension(5) :: unitTirrad ! irradiation time unit (y,d,h,m,s)
  character(len=3)               :: radiounit  ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
  character(len=3)               :: yieldunit  ! unit for isotope yield: num (number), mug, mg, g, or kg
  character(len=3)               :: timeunit   ! unit for time for production output: s, h or d
  character(len=3)               :: currentunit ! unit for accelerator current: muA, mA, or A
  integer, dimension(5)          :: Tcool      ! cooling time per unit cooling time unit (y,d,h,m,s)
  integer, dimension(5)          :: Tirrad     ! irradiation time per unit irradiation time unit
  real(sgl)                      :: Area       ! target area in cm^2
  real(sgl)                      :: thickness  ! target thickness in cm
  real(sgl)                      :: Eback      ! lower end of energy range in MeV for isotope
  real(sgl)                      :: Ebeam      ! incident energy in MeV for isotope production
  real(sgl)                      :: Ibeam_input! beam current in mA for isotope production
  real(sgl)                      :: Ibeam      ! beam current in mA for isotope production
  real(sgl)                      :: rho_target  ! target material density
  real(sgl)                      :: targetmass ! target mass in grams
  real(sgl)                      :: targetmass_input ! target mass in grams
  real(sgl)                      :: fluxtotal  ! total flux
  real(sgl)                      :: fgamma     ! electron-to-photon conversion effiency
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for decay data
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, dimension(-1:numZ, 0:numA, -1:numisom)    :: rtyp      ! type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
  integer, dimension(-1:numZ, 0:numA, -1:numisom, 5) :: Td        ! half life per time unit
  integer                                            :: Zcomp     ! compound nucleus Z
  integer                                            :: Acomp     ! compound nucleus A
  real(sgl)                                          :: daysec    ! number of seconds in a day
  real(sgl)                                          :: hoursec   ! number of seconds in an hour
  real(sgl), dimension(-1:numZ, 0:numA, -1:numisom)  :: lambda    ! decay rate per isotope
  real(sgl)                                          :: minutesec ! number of seconds in a minute
  real(sgl), dimension(-1:numZ, 0:numA, -1:numisom)  :: Thalf     ! half life of nuclide in sec.
  real(sgl)                                          :: yearsec   ! number of seconds in a year
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for isotope production
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(0:numZ, 0:numA, -1:numisom)           :: Yexist       ! logical for yield
  integer                                                  :: Nenrp        ! number of incident energies for residual prod
  real(sgl)                                                :: Leff         ! effective thickness of target
  real(sgl)                                                :: V_target         ! active target volume
  real(sgl)                                                :: M_target         ! active target mass
  real(sgl)                                                :: projnum      ! number of incident particles [s^-1]
  real(sgl)                                                :: heat         ! produced heat
  real(sgl)                                                :: Eaverage    ! average energy of spectrum
  real(sgl), dimension(-1:numZ, -1:numA, -1:numisom)       :: reaction_rate ! reaction rate per isotope
  real(sgl)                                                :: selfshield_av ! average self-shielding factor
  real(sgl), dimension(-1:numZ, -1:numA, -1:numisom)       :: sacs         ! spectrum averaged cross section
  real(sgl), dimension(0:numen)                            :: Erp          ! incident energy
  real(sgl), dimension(0:numen)                            :: xsrp         ! residual production cross section in mb
  integer                                                  :: Ntime        ! number of time points
  integer, dimension(0:numZ,0:numA,-1:numisom)             :: Tmaxactivity ! time of maximum activity of produced isoto
  integer, dimension(0:numZ,0:numA,-1:numisom,5)           :: Tp           ! irradiation time with maximal yield per time unit
  real(sgl)                                                :: N_0          ! number of original target atoms
  real(sgl)                                                :: N_0_nat      ! number of original target atoms
  real(sgl), dimension(0:numtime)                          :: Tgrid        ! time
  real(sgl)                                                :: Tir          ! irradiation time per unit
  real(sgl)                                                :: Tco          ! cooling time per unit
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: Niso         ! number of isotopes produced after irradiation
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: activity     ! activity of produced isotope in MBq
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: Nisonat      ! number of isotopes produced after irradiation
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: activitynat  ! activity of produced isotope in MBq
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: specactivitynat  ! activity of produced isotope in MBq
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: yield        ! yield of produced isotope in MBq/(mA.h)
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: yieldnat     ! yield of produced isotope in MBq/(mA.h)
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: Nisorel      ! fraction of number of produced isotopes per ele
  real(sgl), dimension(0:numZ, 0:numtime)                  :: Nisotot      ! number of elemental isotopes produced after irr
  real(sgl), dimension(0:numZ, 0:numtime)                  :: Nelrel       ! relative amount of produced element
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: Nisorelnat   ! fraction of number of produced isotopes per ele
  real(sgl), dimension(0:numZ, 0:numtime)                  :: Nisototnat   ! number of elemental isotopes produced after irr
  real(sgl), dimension(0:numZ,0:numA,-1:numisom)           :: Tmax         ! irradiation time with maximal yield
  real(sgl), dimension(0:numZ,0:numA,-1:numisom,0:numtime) :: specactivity ! specific activity of produced isotope
  real(sgl), dimension(-1:numZ, -1:numA, -1:numisom)       :: reaction_ratenat ! reaction rate
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for conversion
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=3)               :: rstr       ! string
  character(len=3)               :: ystr       ! string
  character(len=3)               :: tstr       ! string
  character(len=3)               :: cstr       ! string
  real(sgl) :: rfac    ! conversion factor for radioactivity
  real(sgl) :: mfac    ! conversion factor for mass
  real(sgl) :: tfac    ! conversion factor for time
  real(sgl) :: cfac    ! conversion factor for current
end module A0_isotopia_mod
! Copyright A.J. Koning 2026
