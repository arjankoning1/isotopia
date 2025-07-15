program medmerge
!
! Merge IAEA medical isotope file with TENDL for smooth file up to 200 MeV
!
  implicit none
  integer, parameter :: numen=40000
  integer i, Niaea, Ntendl, k, N, Nhead
  integer Nlow, Nup, nen
  integer istat,keyix
  character*132 iaeafile,tendlfile,string(100),str,key
  real Eiaea(numen), xsiaea(numen), dxsiaea(numen)
  real Etendl(numen), xstendl(numen), dxstendl(numen)
  real E(numen), xs(numen), dxs(numen)
  real Elow, Eup, Rlow, Rup
  string = ''
  dxstendl = 0.
  read(*,'(a)') iaeafile 
  read(*,'(a)') tendlfile 
  open (unit=1,status='unknown',file=iaeafile)
  do i = 1, 5
    read(1,'(a)') string(i)
  enddo
  i=1
  do
    read(1,*, iostat = istat) Eiaea(i), xsiaea(i), dxsiaea(i)
    if (istat == -1) then
      Niaea = i - 1
      exit
    else
      i=i+1
    endif
  enddo
  close(1)
  open (unit=2,status='unknown',file=tendlfile)
  Nhead = 0
Loop1:  do
    read(2,'(a)', iostat = istat) str
    if (istat == -1) exit
    Nhead = Nhead + 1
    string(Nhead) = str
    if (str(1:1) /= '#') then
      Nhead = Nhead - 1
      backspace 2
      i=1
      do
        read(2,*, iostat = istat) Etendl(i), xstendl(i)
        if (istat == -1) then
          Ntendl = i - 1
          exit Loop1
        else
          i=i+1
        endif
      enddo
    endif
  enddo Loop1
  close(2)
  Elow = Eiaea(1)
  Eup = Eiaea(Niaea)
  nen = 0
  Rlow = 1.
  do i = 1, Ntendl
    if (Etendl(i).ge.Elow) then
      Nlow = i - 1
      if (xstendl(i) > 0.) Rlow = xsiaea(1) / xstendl(i)
      exit
    endif
  enddo
  Rup = 1.
  do i = Ntendl, 1, -1
    if (Etendl(i).le.Eup) then
      Nup = i + 1
      if (xstendl(i) > 0.) Rup = xsiaea(Niaea) / xstendl(i)
      exit
    endif
  enddo
  k = 0
  do i = 1, Nlow
    k = k + 1
    E(k) = Etendl(i)
    xs(k) = Rlow * xstendl(i)
    dxs(k) = dxstendl(i)
  enddo
  do i = 1, Niaea
    k = k + 1
    E(k) = Eiaea(i)
    xs(k) = xsiaea(i)
    dxs(k) = dxsiaea(i)
  enddo
  do i = Nup , Ntendl
    k = k + 1
    E(k) = Etendl(i)
    xs(k) = Rup * xstendl(i)
    dxs(k) = dxstendl(i)
  enddo
  N = k
  do i =1, Nhead
    key='source'
    keyix=index(string(i),trim(key))
    if (keyix > 0) write(string(i)(keyix+len_trim(key)+2:80),'(a)', iostat = istat) 'IAEA medical isotope consortium'
    key='library'
    keyix=index(string(i),trim(key))
    if (keyix > 0) write(string(i)(keyix+len_trim(key)+2:80),'(a)', iostat = istat) 'iaea.2024 + TENDL-2023'
    key='author'
    keyix=index(string(i),trim(key))
    if (keyix > 0) write(string(i)(keyix+len_trim(key)+2:80),'(a)', iostat = istat) 'IAEA medical isotope consortium'
    key='year'
    keyix=index(string(i),trim(key))
    if (keyix > 0) write(string(i)(keyix+len_trim(key)+2:80),'(a)', iostat = istat) '2024'
    key='columns'
    keyix=index(string(i),trim(key))
    if (keyix > 0) write(string(i)(keyix+len_trim(key)+2:80),'(a)', iostat = istat) '3'
    key='entries'
    keyix=index(string(i),trim(key))
    if (keyix > 0) write(string(i)(keyix+len_trim(key)+2:keyix+len_trim(key)+8),'(i7)', iostat = istat) N
  enddo
  string(Nhead-1) = '##       E             xs            dxs'
  string(Nhead)   = '##     [MeV]          [mb]           [mb]'
  do i =1, Nhead
    write(*,'(a)', iostat = istat) trim(string(i))
  enddo
  do i = 1, N
    write(*,'(3es15.6)') E(i), xs(i), dxs(i)
  enddo
end program medmerge
