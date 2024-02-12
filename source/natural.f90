subroutine natural
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Calculation for natural elements
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2023-02-25   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
  use A0_isotopia_mod
  use A1_error_handling_mod
!
!              sgl, &       ! single precision kind
!              abun, &      ! natural abundance
!              activity, &  ! activity of produced isotope in MBq
!              Eback, &     ! lower end of energy range in MeV for isotope
!              Ebeam, &     ! incident energy in MeV for isotope production
!              Ibeam, &     ! beam current in mA for isotope production
!              iso, &       ! counter for isotope
!              isochar, &   ! symbol of isomer
!              isonum, &    ! number of isotopes in element
!              isotope, &   ! isotope of natural element
!              k0, &        ! index of incident particle
!              lambda, &    ! decay rate per isotope
!              natstring, & ! string extension for file names
!              Niso, &      ! number of isotopes produced after irradiation
!              Nisorel, &   ! fraction of number of produced isotopes per ele
!              Nisotot, &   ! number of elemental isotopes produced after irr
!              Ntar0, &     ! number of original target atoms
!              Ntime, &     ! number of time points
!              nuc, &       ! symbol of nucleus
!              numA, &      ! number of masses
!              numtime, &   ! number of time points
!              numZ, &      ! number of elements
!              parsym, &    ! symbol of particle
!              prate, &     ! production rate per isotope
!              Tco, &       ! cooling time per unit
!              Tcool, &     ! cooling time per unit cooling time unit (y, d, h, m, s)
!              Td, &        ! half life per time unit
!              Tgrid, &     ! time
!              Thalf, &     ! half life of nuclide in sec.
!              Tir, &       ! irradiation time per unit
!              Tirrad, &    ! irradiation time per unit irradiation time unit
!              Tp, &        ! irradiation time with maximal yield per time unit
!              yield, &     ! yield of produced isotope in MBq / (mA.h)
!              Ztarget      ! charge number of target nucleus
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist                            ! logical to determine existence
  logical           :: resexist(0:numZ, 0:numA, -1:1)    ! logical to determine existence of residual production file
  character(len=3)  :: Astr                              !
  character(len=3)  :: massstring !
  character(len=6)  :: finalnuclide !
  character(len=13) :: state                             ! state of final nuclide
  character(len=15) :: Yfile0                            !
  character(len=20) :: Yf
  character(len=15) :: Yfile(0:numZ, 0:numA, -1:1)       ! file with production yields
  character(len=40) :: form1                             ! format
  character(len=40) :: form2                             ! format
  character(len=16) :: reaction   ! reaction
  character(len=16) :: col(6)    ! header
  character(len=16) :: un(6)    ! header
  character(len=80) :: quantity   ! quantity
  character(len=132) :: string    !
  character(len=132) :: topline    ! topline
  character(len=132) :: line        !
  integer           :: Ncol       ! counter
  integer           :: ia                                ! mass number from abundance table
  integer           :: is                                ! isotope counter: -1=total, 0=ground state 1=isomer
  integer           :: it                                ! counter for tritons
  integer           :: iz                                ! charge number of residual nucleus
  integer           :: k                                 ! designator for particle
  integer           :: istat                             ! istat
  real(sgl)         :: Nis                               ! number of produced isotopes
  real(sgl)         :: Tdum                              !
  real(sgl)         :: Y                                 ! product yield (in ENDF-6 format)
  real(sgl)         :: Ym                                !
  real(sgl)         :: Th                                !
!
! **************** Create runs and directories per isotope *************
!
! isotopiainput  : subroutine for input
! isotopiainitial: subroutine for initialization
! production     : subroutine to calculate the isotope production
!
! For mono-isotopic nuclides we are done.
!
  if (isonum == 1) return
!
! do a full ISOTOPIA calculation for each isotope
!
  do iso = 2, isonum
    call isotopiainput
    call isotopiainitial
    call production
  enddo
!
! ******** Merge output files in results for natural elements **********
!
  quantity='Isotope production'
  reaction='('//ptype0//',x)'
  do iz = Zcomp, Zcomp - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, 1
        resexist(iz, ia, is) = .false.
        do iso = 1, isonum
          if (flagZAoutput) then
            Yfile0 = 'Y000000.tot'
            write(Yfile0(2:7), '(2i3.3)') iz, ia
            if (is == 0) Yfile0(9:11) = 'L00'
            if (is == 1) Yfile0(9:11) = 'L01'
          else
            write(Astr(1:3), '(i3.3)') ia
            Yfile0 = trim(nuc(iz)) //Astr
            if (is == 0) Yfile0 = trim(Yfile0)//'g'
            if (is == 1) Yfile0 = trim(Yfile0)//'m'
            Yfile0 = trim(Yfile0)//'.act'
          endif
          if (is >= 1) then
            state = ' Isomer=     '
            write(state(10:10), '(a1)') isochar(is)
          else
            state = ' Ground state'
          endif
          inquire (file = trim(Yfile0) //natstring(iso), exist = lexist)
          if (lexist) then
            resexist(iz, ia, is) = .true.
            Yfile(iz, ia, is) = trim(Yfile0)
            Yf=trim(Yfile0)//natstring(iso)
            open (2, status = 'old', file = Yf)
            do 
              read(2,'(a)',iostat = istat) line
              if (istat == -1) exit
              if (line(1:9) == '# Entries') then
                read(line(42:47),'(i6)', iostat = istat) Ntime
                if (istat /= 0) call read_error(Yf, istat)
                read(2,'()')
                do it = 1, Ntime
                  read(2, * , iostat = istat) Tdum, Y, Nis, Ym
                  if (istat < 0) call read_error(Yf, istat)
                  if (istat > 0) cycle
                  activitynat(iz, ia, is, it) = activitynat(iz, ia, is, it) + abun(iso) * Y
                  yieldnat(iz, ia, is, it) = yieldnat(iz, ia, is, it) + abun(iso) * Ym
                  Nisonat(iz, ia, is, it) = Nisonat(iz, ia, is, it) + abun(iso) * Nis
                  Nisototnat(iz, it) = Nisototnat(iz, it) + abun(iso) * Nis
                enddo
                exit
              endif
            enddo
            close (unit = 2)
          endif
        enddo
      enddo
    enddo
Loop1: do ia = Acomp, Acomp - Adepth, -1
      do is = -1, 1
        if ( .not. resexist(iz, ia, is)) cycle Loop1
        do it = 1, Ntime
          if (Nisototnat(iz, it) /= 0.) Nisorelnat(iz, ia, is, it) = Nisonat(iz, ia, is, it) / Nisototnat(iz, it)
        enddo
      enddo
    enddo Loop1
    do it = 1, Ntime
      Nelrel(iz, it) = Nisototnat(iz, it) / Ntar0
    enddo
  enddo
  write(*,'()')
  do iz = Zcomp, Zcomp - Zdepth, -1
    do ia = Acomp, Acomp - Adepth, -1
      do is = -1, 1
        if (resexist(iz, ia, is)) then
          massstring='   '
          write(massstring,'(i3)') ia
          finalnuclide=trim(nuc(iz))//adjustl(massstring)//isochar(is)
          write(*,'("file: ",a)') trim(Yfile(iz, ia, is))
          open (1, status = 'unknown', file = Yfile(iz, ia, is))
          topline=trim(targetnuclide)//trim(reaction)//trim(finalnuclide)//' '//trim(quantity)
          call write_header(topline,source,user,date,oformat)
          call write_target
          call write_reaction(reaction,0.d0,0.d0,6,5)
          call write_residual(iz,ia,finalnuclide)
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
          string='      years     days    hours  minutes seconds'
          string=''
          write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tirrad(k), k = 1, 5)
          call write_char(2,'Irradiation time',string)
          write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') (Tcool(k), k = 1, 5)
          call write_char(2,'Cooling time',string)
          if (Thalf(iz, ia, is) > 1.e17) then
            string='stable'
          else
            write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")') &
 &   (Td(iz, ia, is, k), k = 1, 5)
          endif
          call write_char(2,'Half-life',string)
          if (Tmax(iz, ia, is) > 1.e17) then
            string='infinity'
          else
            write(string, '(" ",i6, " years ", i3, " days", i3, " hours", i3, " minutes", i3, " seconds ")')  &
 &   (Tp(iz, ia, is, k), k = 1, 5)
          endif
          call write_char(2,'Maximum production at',string)
          un = ''
          col(1)='Time'
          un(1)='sec'
          col(2)='Activity'
          un(2)= trim(rstr)
          col(3)='Isotopes'
          un(3)= trim(ystr)
          col(4)='Yield'
          un(4)= trim(rstr)
          col(5)='Isotopic_frac.'
          col(6)='Time'
          un(6)='h'
          Ncol=6
          call write_datablock(quantity,Ncol,numtime,col,un)
          do it = 1, numtime
            Th = Tgrid(it)/hoursec
            write(1, '(6es15.6)') Tgrid(it),  activitynat(iz, ia, is, it), Nisonat(iz, ia, is, it), &
 &            yieldnat(iz, ia, is, it), Nisorelnat(iz, ia, is, it), Th
          enddo
          close (unit = 1)
        endif
      enddo
    enddo
  enddo
  form1='("# Time [h]   ",xx(a2,"      "))'
  write(form1(18:19), '(i2)') Zdepth+1
  form2='(f8.1,xxf8.5)'
  write(form2(7:8), '(i2)') Zdepth+1
  open (2, status = 'unknown', file = 'elements.out')
  write(2, fmt = form1) (nuc(iz), iz = Zcomp, Zcomp - Zdepth, -1)
  do it = 1, Ntime
    write(2, fmt = form2) Tgrid(it), (Nelrel(iz, it), iz = Zcomp, Zcomp - Zdepth, -1)
  enddo
  write( * , * ) " End of ISOTOPIA for natural target"
  return
end subroutine natural
! Copyright A.J. Koning 2023
