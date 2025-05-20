  subroutine constants
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Constants and basic properties of particles
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-05-20   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_isotopia_mod
!
!              amu, &      ! atomic mass unit in MeV
!              avogadro, & ! Avogadro's number
!              iso, &      ! counter for isotope
!              isochar, &  ! symbol of isomer
!              nuc, &      ! symbol of nucleus
!              numiso, &   ! maximum number of isotopes per element
!              numisom, &  ! number of isomers
!              numZ, &     ! number of elements
!              parA, &     ! mass number of particle
!              parmass, &  ! mass of particle in a.m.u.
!              parname, &  ! name of particle
!              parsym, &   ! symbol of particle
!              parZ, &     ! charge number of particle
!              qelem       ! elementary charge in C
!
! *** Declaration of local data
!
  implicit none
!
! ****************** General properties of particles *******************
!
!          photon  = 0
!          neutron = 1
!          proton  = 2
!          deuteron= 3
!          triton  = 4
!          helium-3= 5
!          alpha   = 6
!
  parname = (/'gamma   ', 'neutron ', 'proton  ', 'deuteron', 'triton  ', 'helium-3', 'alpha   '/)
  parsym =  (/'g', 'n', 'p', 'd', 't', 'h', 'a'/)
  parZ =    (/ 0, 0, 1, 1, 1, 2, 2 /)
  parA =    (/ 0, 1, 1, 2, 3, 3, 4 /)
  parmass = (/ 0., 1.00866491574, 1.00727646662, 2.01355321275, 3.01550071621, 3.01493224717, 4.0015061791 /)
!
! ************************ Nuclear symbols *****************************
!
  nuc = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
    'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
    'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
    'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
    'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', &
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og', 'B9', 'C0', 'C1', 'C2', 'C3', 'C4'/)
  isochar =                               (/ ' ', 'g', 'm', 'n'/)
!
! *********************** Fundamental constants ************************
!
  amu = 931.49410242
  emass = 0.510998950
  avogadro = 6.02214076e23
  qelem = 1.602176634e-19
end subroutine constants
! Copyright A.J. Koning 2025
