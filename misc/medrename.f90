program medrename
!
! Read IAEA medical isotope files provided by Marco Verpelli and put them into exact libraries format
!
  character*132 string(5)
  do i = 1,5
    read(*,'(a)') string(i)
  enddo
  read(string(4),'(14x,i6)') N
  write(*,'(a)') trim(string(1))
  write(*,'("# IAEA-MED         EVAL-2024")')
  write(*,'("# uncertainties: y")')
  write(*,'("# # energies =",i6)') N
  write(*,'("#  E(MeV)          xs(mb)         dxs(mb)")')
  do
    read(*, *, iostat = istat) E, xs, dxs
    if (istat == -1) exit
    write(*,'(es12.5,2es15.5)') E, xs, dxs
  enddo
end program medrename
