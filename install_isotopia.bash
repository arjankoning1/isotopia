#!/bin/bash
unalias -a
echo
echo "       ISOTOPIA installation (Version December 29 2024) (C) Copyright 2024 Arjan Koning All Rights Reserved"
echo
echo " Two ways to use this script:"
echo
echo " install_isotopia.bash 'Arjan Koning' 'iaea.2026'"
echo "     Replace  Arjan Koning by your own name"
echo "     Replace  iaea.2026 by the nuclear data library being used"
echo
echo " or"
echo
echo " install_isotopia.bash"
echo "     after which you will be prompted to input your name"
echo "     and then the nuclear data library name e.g. iaea.2026"
echo
if [ $# -eq 2 ] ; then
  yourname=$1
  libname=$2
else
  echo 'Enter your name (which will appear in the output files): '
  read yourname
  echo 'Enter the nuclear data library being used: '
  read libname
fi
echo ${yourname}
echo ${libname}
pfile=path_change.bash
if [ -e $pfile ] ; then
  sed "s/user=.*/user='${yourname}'/" $pfile > tmp
  sed "s/user=.*/user='${yourname}'/" $pfile > tmp
  mv -f tmp $pfile
  sed "s/libname=.*/libname='${libname}'/" $pfile > tmp
  mv -f tmp $pfile
fi
chmod a+x $pfile
./code_build.bash isotopia
