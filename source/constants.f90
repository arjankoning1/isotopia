  subroutine constants
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Constants and basic properties of particles
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2024-12-29   A.J. Koning    A     Original code
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
  parmass = (/ 0., 1.00866491590, 1.007825031898, 2.014101777844, 3.01604928132, 3.01602932197, 4.00260325413 /)
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
  amu = 931.49386
  emass = 0.510999
  avogadro = 6.0221367e23
  qelem = 1.6021773e-19
end   subroutine constants
! Copyright A.J. Koning 2023
