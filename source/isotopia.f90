program isotopia
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Main program - Prediction of radio-isotope production
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-05-20   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
!   |-------------------------------------------------------|
!   |                 Arjan Koning                          |
!   |                                                       |
!   | Email: A.Koning@@iaea.org                             |
!   |-------------------------------------------------------|
!
! MIT License
!
! Copyright (c) 2025 Arjan Koning
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! *** Use data from other modules
!
  use A0_isotopia_mod
!
! use A1_set_main_mod, only: & ! Variables for main input
!              flagnatural    ! flag for calculation of natural element
!
! ********** Input, initialization and calculation of yields ***********
!
! machine        : subroutine for machine dependent statements
! constants      : subroutine for constants
! isotopiainput  : subroutine for input
! isotopiainitial: subroutine for initialization
! production     : subroutine to calculate the isotope production
! flagnatural    : flag for calculation of natural element
! natural        : subroutine for calculation of natural element
!
  call machine
  call constants
  call isotopiainput
  call isotopiainitial
  call production
  if (flagnatural) call natural
  call timer
end program isotopia
! Copyright A.J. Koning 2025
