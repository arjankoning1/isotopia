subroutine prodout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write output
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
!
!              activity, &     ! activity of produced isotope in MBq
!              Area, &         ! target area in cm^2
!              Atarget, &      ! mass number of target nucleus
!              Eback, &        ! lower end of energy range in MeV for isotope
!              Ebeam, &        ! incident energy in MeV for isotope production
!              heat, &         ! produced heat
!              Ibeam, &        ! beam current in mA for isotope production
!              iso, &          ! counter for isotope
!              isochar, &      ! symbol of isomer
!              isotope, &      ! isotope of natural element
!              k0, &           ! index of incident particle
!              lambda, &       ! decay rate per isotope
!              Mtar, &         ! active target mass
!              natstring, &    ! string extension for file names
!              Niso, &         ! number of isotopes produced after irradiation
!              Nisorel, &      ! fraction of number of produced isotopes per ele
!              Ntar0, &        ! number of original target atoms
!              Ntime, &        ! number of time points
!              nuc, &          ! symbol of nucleus
!              numtime, &      ! number of time points
!              parsym, &       ! symbol of particle
!              prate, &        ! production rate per isotope
!              projnum, &      ! number of incident particles [s^ -1]
!              ptype0, &       ! type of incident particle
!              radiounit, &    ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              rhotarget, &    ! target material density
!              targetdx, &     ! effective thickness of target
!              Tco, &          ! cooling time per unit
!              Tcool, &        ! cooling time per unit cooling time unit (y, d, h, m, s)
!              Td, &           ! half life per time unit
!              Tgrid, &        ! time
!              Thalf, &        ! half life of nuclide in sec.
!              Tir, &          ! irradiation time per unit
!              Tirrad, &       ! irradiation time per unit irradiation time unit
!              Tmax, &         ! irradiation time with maximal yield
!              Tmaxactivity, & ! time of maximum activity of produced isoto
!              Tp, &           ! irradiation time with maximal yield per time unit
!              Vtar, &         ! active target volume
!              yield, &        ! yield of produced isotope in MBq / (mA.h)
!              yieldunit, &    ! unit for isotope yield: num (number), mug, mg, g, or kg
!              Ztarget         ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  character(len=3)  :: Astr        !
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=13) :: state       ! state of final nuclide
  character(len=15) :: Yfile       ! file with production yields
  character(len=38) :: halflife    ! half life
  character(len=38) :: maxprod     ! maximum production
  character(len=16) :: reaction   ! reaction
  character(len=15) :: col(6)    ! header
  character(len=15) :: un(6)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: string    !
  character(len=132) :: topline    ! topline
  integer           :: Ncol       ! counter
  integer           :: ia          ! mass number from abundance table
  integer           :: is          ! isotope counter: -1=total, 0=ground state 1=isomer
  integer           :: it          ! counter for tritons
  integer           :: iz          ! charge number of residual nucleus
  integer           :: k           ! designator for particle
  real(sgl)         :: Th          ! time in hours
!
! ************************* Main output ********************************
!
  rstr = 'MBq'
  ystr = '   '
  if (radiounit == 'bq') rstr = ' Bq'
  if (radiounit == 'kbq') rstr = 'KBq'
  if (radiounit == 'gbq') rstr = 'GBq'
  if (radiounit == 'ci') rstr = ' Ci'
  if (radiounit == 'kci') rstr = 'KCi'
  if (radiounit == 'mci') rstr = 'mCi'
  if (yieldunit == 'g') ystr = '  g'
  if (yieldunit == 'mug') ystr = 'mug'
  if (yieldunit == 'mg') ystr = ' mg'
  if (yieldunit == 'kg') ystr = ' kg'
  write(*, '(/" Summary of isotope production for ", a1, " + ", a/)') ptype0, trim(targetnuclide)
  string=''
  write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tirrad(k), k = 1, 5) 
  write(*, '(" Maximal irradiation time: ", a)') trim(string)
  write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tcool(k), k = 1, 5) 
  write(*, '(" Cooling time: ", a)') trim(string)
  write(*, '(" E-beam [MeV]:", es15.6)') Ebeam
  write(*, '(" E-back [MeV]:", es15.6)') Eback
  write(*, '(" Beam current [mA]: ", f12.3, " mA")') Ibeam
  write(*, '(" Target material density [g/cm^3]:", es15.6)') rhotarget
  write(*, '(" Target area [cm^2]:", es15.6)') Area
  write(*, '(" Effective target thickness [cm]:", es15.6)') targetdx
  write(*, '(" Effective target volume [cm^3]:", es15.6)') Vtar
  write(*, '(" Effective target mass [g]:", es15.6)') Mtar
  write(*, '(" Number of target atoms: ", es15.6)') Ntar0
  write(*, '(" Number of incident particles [s^-1]:", es15.6)') projnum
  write(*, '(" Produced heat in target [kW]:", es15.6)') heat
  write(*, '(/" (Maximum) production and decay rates per isotope"/)')
  write(*, '(" Total production rate [s^-1]:", es15.6/)') prate(0, 0, -1)
  write(*, '("#  Nuc     Production rate Decay rate     Activity       #isotopes      Yield          Isotopic", &
 &  "                   Half-life               Time of maximum production")')
  write(*, '("#             [s^-1]         [s^-1]         [", a3, "]          [", a3, "]        [", a3, "/mAh]      fraction")') &
 & rstr, ystr, rstr
  do iz = Zcomp, Zcomp - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        if ( .not. Yexist(iz, ia, is)) cycle
        it = int(Tmaxactivity(iz, ia, is))
        halflife = '                                      '
        if (Thalf(iz, ia, is) > 1.e17) then
          write(halflife, '(a15)') '       stable  '
        else
          write(halflife, '(i11, " y ", i3, " d ", i3, " h ", i3, " m ", i3, " s ")') (Td(iz, ia, is, k), k = 1, 5)
        endif
        maxprod = '                                      '
        if (Tmax(iz, ia, is) > 1.e17) then
          write(maxprod, '(a15)') '       infinite'
        else
          write(maxprod, '(i11, " y ", i3, " d ", i3, " h ", i3, " m ", i3, " s ")') (Tp(iz, ia, is, k), k = 1, 5)
        endif
        if (is.eq.-1 .and. Niso(iz,ia,0,it) > 0. .and. Niso(iz,ia,1,it) > 0.) then
          halflife = '                                      '
          maxprod = '                                      '
        endif
        write(*, '(1x, a2, i4, 1x, a1, 5es15.6, f10.5, 2a38)') nuc(iz), ia, isochar(is), prate(iz, ia, is), &
 &        lambda(iz, ia, is), activity(iz, ia, is, it), Niso(iz, ia, is, it), &
 &        yield(iz, ia, is, 1), Nisorel(iz, ia, is, it), halflife, maxprod
      enddo
    enddo
  enddo
!
! Output to files per residual product
!
  reaction='('//ptype0//',x)'
  quantity='Isotope production'
  write(*,'()')
  do iz = Zcomp, Zcomp - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        if ( .not. Yexist(iz, ia, is)) cycle
        if (flagZAoutput) then
          Yfile = 'Y000000.tot'//natstring(iso)
          write(Yfile(2:7), '(2i3.3)') iz, ia
          if (is >= 0) Yfile(9:11) = 'L00'
          if (is >= 1) write(Yfile(10:11), '(i2.2)') is
        else
          Astr = '000'
          write(Astr(1:3), '(i3.3)') ia
          Yfile = trim(nuc(iz)) //Astr
          if (is == 0) Yfile = trim(Yfile)//'g'
          if (is == 1) Yfile = trim(Yfile)//'m'
          Yfile = trim(Yfile)//'.act'//natstring(iso)
        endif
        if (is >= 1) then
          state = ' Isomer=     '
          write(state(10:10), '(a1)') isochar(is)
        else
          state = ' Ground state'
        endif
        massstring='   '
        write(massstring,'(i3)') ia
        finalnuclide=trim(nuc(iz))//adjustl(massstring)//isochar(is)
        write(*,'("file: ",a)') trim(Yfile)
        open (unit = 1, file = Yfile, status = 'replace')
        topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
        call write_header(topline,source,user,date,oformat)
        call write_target
        call write_reaction(reaction,0.d0,0.d0,6,5)
        call write_residual(iz,ia,finalnuclide)
        write(1,'("# parameters:")')
        call write_real(2,'Beam current [mA]',Ibeam)
        call write_real(2,'E-Beam [MeV]',Ebeam)
        call write_real(2,'E-Back [MeV]',Eback)
        string='Initial production rate [s^-1]'
        call write_real(2,string,prate(iz, ia, is))
        string='Decay rate [s^-1]'
        call write_real(2,string,lambda(iz, ia, is))
        string='Initial production yield ['//rstr//'/mAh]'
        call write_real(2,string,yield(iz, ia, is, 1))
        string='Total activity at EOI ['//rstr//']'
        call write_real(2,string,activity(iz, ia, is, Ntime))
        string=''
        write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tirrad(k), k = 1, 5)
        call write_char(2,'Irradiation time',string)
        write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tcool(k), k = 1, 5)
        call write_char(2,'Cooling time',string)
        if (Thalf(iz, ia, is) > 1.e17) then
          string='stable'
        else
          write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') &
 & (Td(iz, ia, is, k), k = 1, 5)
        endif
        call write_char(2,'Half-life',string)
        if (Tmax(iz, ia, is) > 1.e17) then
          string='infinity'
        else
          write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")')  &
 & (Tp(iz, ia, is, k), k = 1, 5)
        endif
        call write_char(2,'Maximum production at',string)
        un = ''
        col(1)='Time'
        un(1) = 'sec'
        col(2)='Activity'
        un(2) = trim(rstr)
        col(3)='Isotopes'
        un(3) = trim(ystr)
        col(4)='Yield'
        un(4) = trim(rstr)
        col(5)='Isotopic frac.'
        col(6)='Time'
        un(1) = 'h'
        Ncol=6
        call write_datablock(quantity,Ncol,numtime,col,un)
        do it = 1, numtime
          Th = Tgrid(it)/hoursec
          write(1, '(6es15.6)') Tgrid(it), activity(iz, ia, is, it), Niso(iz, ia, is, it), &
 &          yield(iz, ia, is, it), Nisorel(iz, ia, is, it), Th
        enddo
        close (unit = 1)
      enddo
    enddo
  enddo
  write(*, '(/"  End of ISOTOPIA calculation for ", a1, " + ", a)') ptype0, trim(targetnuclide)
  return
end subroutine prodout
! Copyright A.J. Koning 2023
