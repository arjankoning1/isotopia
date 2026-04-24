subroutine prodout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write output
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2026-04-23   A.J. Koning    A     Original code
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
!              M_target, &         ! active target mass
!              natstring, &    ! string extension for file names
!              Niso, &         ! number of isotopes produced after irradiation
!              Nisorel, &      ! fraction of number of produced isotopes per ele
!              N_0, &          ! number of original target atoms
!              Ntime, &        ! number of time points
!              nuc, &          ! symbol of nucleus
!              numtime, &      ! number of time points
!              parsym, &       ! symbol of particle
!              reaction_rate, &  ! reaction rate per isotope
!              projnum, &      ! number of incident particles [s^ -1]
!              ptype0, &       ! type of incident particle
!              radiounit, &    ! unit for radioactivity: Bq, kBq, MBq, Gbq, mCi,
!              rho_target, &    ! target material density
!              Leff, &     ! effective thickness of target
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
!              V_target, &         ! active target volume
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
  character(len=12) :: yieldstring     
  character(len=13) :: state       ! state of final nuclide
  character(len=15) :: Yfile       ! file with production yields
  character(len=38) :: halflife    ! half life
  character(len=38) :: maxprod     ! maximum production
  character(len=16) :: reaction   ! reaction
  character(len=15) :: col(7)    ! header
  character(len=15) :: un(7)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: string    !
  character(len=132) :: topline    ! topline
  integer           :: Ncol       ! counter
  integer           :: ia          ! mass number from abundance table
  integer           :: is          ! isotope counter: -1=total, 0=ground state 1=isomer
  integer           :: it          ! counter for tritons
  integer           :: iz          ! charge number of residual nucleus
  integer           :: k           ! designator for particle
  integer           :: indent
  integer           :: id2
  integer           :: id4
  real(sgl)         :: Th          ! time in hours
  real(sgl)         :: act_out          ! help variable
  real(sgl)         :: specact_out          ! help variable
  real(sgl)         :: yield_out          ! help variable
  real(sgl)         :: Niso_out          ! help variable
  real(sgl)         :: Nisorel_out          ! help variable
  real(sgl)         :: yfac
!
! ************************* Main output ********************************
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
! 
! Normalization and strings for output
!
  call conversion
  if (k0 > 1) then
    yfac =  rfac * cfac * tfac
  else
    yfac =  rfac * mfac * tfac
  endif
  write(*, '(/" Summary of isotope production for ", a1, " + ", a/)') ptype0, trim(targetnuclide)
  string=''
  write(string, '(es15.6," (",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds)")') &
 &  Tir, (Tirrad(k), k = 1, 5) 
  write(*, '(" Maximal irradiation time [s]: ", a)') trim(string)
  write(string, '(es15.6," (",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds)")') &
 &  Tco, (Tcool(k), k = 1, 5) 
  write(*, '(" Cooling time [s]: ", a)') trim(string)
  if (k0 > 1) then
    write(*, '(" E-beam [MeV]:", es15.6)') Ebeam
    write(*, '(" E-back [MeV]:", es15.6)') Eback
  endif
  if (k0 /= 1) write(*, '(" Beam current [mA]: ", f12.3)') Ibeam
  write(*, '(" Target material density [g/cm^3]:", es15.6)') rho_target
  write(*, '(" Target area [cm^2]:", es15.6)') Area
  if (k0 > 1) then
    write(*, '(" Effective target thickness [cm]:", es15.6)') Leff
  else
    write(*, '(" Target thickness [cm]:", es15.6)') thickness
  endif
  write(*, '(" Effective target volume [cm^3]:", es15.6)') V_target
  write(*, '(" Effective target mass [",a,"]:", es15.6)') trim(ystr), M_target * mfac
  write(*, '(" Average energy of spectrum [MeV]: ", es15.6)') Eaverage
  write(*, '(" Number of target atoms: ", es15.6)') N_0
  if (k0 > 1) then
    write(*, '(" Number of incident particles [s^-1]:", es15.6)') projnum
    write(*, '(" Produced heat in target [kW]:", es15.6)') heat
  else
    write(*, '(" Total flux [cm^-2s^-1]:", es15.6)') fluxtotal
  endif
  write(*, '(/" (Maximum) production and decay rates per isotope"/)')
  write(*, '(" Total production rate [s^-1]:", es15.6/)') reaction_rate(0, 0, -1)
  write(*, '("#  Nuc       Activity     Spec_activity    Prod_rate      #Isotopes Isot_fraction ", &
 &  " Reac_constant   Decay_constant            Half-life                      Time of maximum production")')
  if (k0 > 1) then
    yieldstring = trim(rstr)//'/('//trim(cstr)//'.'//trim(tstr)//')'
  else
    yieldstring = trim(rstr)//'/('//trim(ystr)//'.'//trim(tstr)//')'
  endif
  write(*, '("#              [", a, "]          [", a, "/", a, "]        [", a, "]        []            []", &
 &  "        [s^-1]         [s^-1]       ")') trim(rstr), trim(rstr), trim(ystr), trim(yieldstring)
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, Nisomer(iz, ia)
        if ( .not. Yexist(iz, ia, is)) cycle
        it = int(Tmaxactivity(iz, ia, is))
        if (it == 0) it = Ntime
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
        if (is == -1 .and. Niso(iz,ia,0,it) > 0. .and. Niso(iz,ia,1,it) > 0.) then
          halflife = '                                      '
          maxprod = '                                      '
        endif
        write(*, '(1x, a2, i4, 1x, a1, 4es15.6, f10.5, 2es15.6, 2a38)') nuc(iz), ia, isochar(is), &
 &         activity(iz, ia, is, it) * rfac, specactivity(iz, ia, is, it) * rfac * mfac, &
 &        yield(iz, ia, is, 1) * yfac, Niso(iz, ia, is, it), Nisorel(iz, ia, is, it), &
 &        reaction_rate(iz, ia, is), lambda(iz, ia, is), halflife, maxprod
      enddo
    enddo
  enddo
!
! Output of files per residual product
!
  reaction='('//ptype0//',x)'
  quantity='Isotope production'
  write(*,'()')
  do iz = Zcomp + 1, Zcomp + 1 - Zdepth, -1
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
        finalnuclide=trim(nuc(iz))//trim(adjustl(massstring))//isochar(is)
        write(*,'("file: ",a)') trim(Yfile)
        open (unit = 1, file = Yfile, status = 'replace')
        topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
        call write_header(indent,topline,source,user,date,oformat)
        call write_target(indent)
        call write_reaction(indent,reaction,0.d0,0.d0,0,0)
        call write_residual(id2,iz,ia,finalnuclide)
        call write_char(id2,'parameters','')
        if (k0 > 1) then
          call write_real(id4,'E-Beam [MeV]',Ebeam)
          call write_real(id4,'E-Back [MeV]',Eback)
        endif
        if (k0 /= 1) call write_real(id4,'Beam current [mA]',Ibeam)
        string='Reaction constant [s^-1]'
        call write_real(id4,string,reaction_rate(iz, ia, is))
        string='Decay constant [s^-1]'
        call write_real(id4,string,lambda(iz, ia, is))
        string='Spectrum averaged cross section [mb]'
        call write_real(id4,string,sacs(iz, ia, is))
        if (k0 == 1) then
          string='Average self-shielding factor'
          call write_real(id4,string,selfshield(iz, ia, is))
        endif
        string='Initial production rate ['//trim(yieldstring)//']'
        call write_real(id4,string,yield(iz, ia, is, 1) * yfac)
        string='Total activity at EOI ['//trim(rstr)//']'
        call write_real(id4,string,activity(iz, ia, is, Ntime) * rfac)
        string='Specific activity at EOI ['//trim(rstr)//'/'//trim(ystr)//']'
        call write_real(id4,string,specactivity(iz, ia, is, Ntime) * rfac * mfac)
        string=''
        write(string, '(es15.6," (",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds)")') &
 &        Tir, (Tirrad(k), k = 1, 5)
        call write_char(id4,'Irradiation time [s]',string)
        write(string, '(es15.6," (",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds)")') &
 &        Tco, (Tcool(k), k = 1, 5)
        call write_char(id4,'Cooling time [s]',string)
        if (Thalf(iz, ia, is) > 1.e17) then
          string='stable'
        else
          write(string, '(es15.6," (",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds)")') &
 &          Thalf(iz, ia, is), (Td(iz, ia, is, k), k = 1, 5)
        endif
        call write_char(id4,'Half-life [s]',string)
        if (Tmax(iz, ia, is) > 1.e17) then
          string='infinity'
        else
          write(string, '(es15.6," (",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds)")')  &
 &         Tmax(iz, ia, is), (Tp(iz, ia, is, k), k = 1, 5)
        endif
        call write_char(id4,'Maximum production [s]',string)
        un = ''
        col(1)='Time'
        un(1) = 's'
        col(2)='Activity'
        un(2) = trim(rstr)
        col(3)='Spec_activity'
        un(3) = trim(rstr)//'/'//trim(ystr)
        col(4)='Prod_rate'
        un(4) = trim(yieldstring)
        col(5)='#Isotopes'
        col(6)='Isot_fraction'
        col(7)='Time'
        un(7) = 'h'
        Ncol=7
        call write_quantity(id2,quantity)
        call write_datablock(id2,Ncol,numtime,col,un)
        do it = 1, numtime
          act_out = activity(iz, ia, is, it) * rfac
          specact_out = specactivity(iz, ia, is, it) * rfac * mfac
          yield_out = yield(iz, ia, is, it) * yfac
          Niso_out = Niso(iz, ia, is, it)
          Nisorel_out = Nisorel(iz, ia, is, it)
          Th = Tgrid(it)/hoursec
          write(1, '(7es15.6)') Tgrid(it), act_out, specact_out, yield_out, Niso_out, Nisorel_out, Th
        enddo
        close (unit = 1)
      enddo
    enddo
  enddo
  write(*, '(/,"  End of ISOTOPIA calculation for ", a1, " + ", a)') ptype0, trim(targetnuclide)
  return
end subroutine prodout
! Copyright A.J. Koning 2026
