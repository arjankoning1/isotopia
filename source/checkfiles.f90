subroutine checkfiles
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Check existence of files
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2022-09-18   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              Atarget, &   ! mass number of target nucleus
!              crosspath, & ! directory containing cross sections
!              k0, &        ! index of incident particle
!              libname, &   ! library name
!              parsym, &    ! symbol of particle
!              path, &      ! directory containing files to be read
!              Starget, &   ! symbol of target nucleus
!              Ztarget      ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist        ! logical to determine existence
  character(len=3)   :: Astr          !
  character(len=4)   :: decaystring   !
  character(len=16)  :: nucpath       !
  character(len=132) :: pfile         ! parameter file
!
! ************************ Check cross section files *******************
!
  Astr = '   '
  write(Astr(1:3), '(i3.3)') Atarget
  nuclide = trim(Starget) //Astr
  Astr = '   '
  write(Astr, '(i3)') Atarget
  targetnuclide = trim(Starget) // adjustl(Astr)
  nucpath = parsym(k0)//'/'//trim(nuclide)//'/'
  xspath = trim(crosspath)//trim(nucpath)//trim(libname)//'/tables/'
  pfile = trim(xspath)//'xs/'//parsym(k0)//'-'//trim(nuclide)// '-MT005.'//trim(libname)
  inquire (file = trim(pfile), exist = lexist)
  if ( .not. lexist) then
    write(*, '(" ISOTOPIA-error: Cross section file does not exist: " , a)') trim(pfile)
    stop
  endif
!
! ************************ Check decay files ***************************
!
  decaypath = trim(path)//'files/decay/'
  decaystring = 'z000'
  write(decaystring(2:4), '(i3.3)') Ztarget
  inquire (file = trim(decaypath) //decaystring, exist = lexist)
  if ( .not. lexist) then
    write(*, '(" ISOTOPIA-error: Decay data file does not exist: ", a)') trim(decaypath)//decaystring
    stop
  endif
  return
end subroutine checkfiles
! Copyright A.J. Koning 2023
