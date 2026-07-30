#!/bin/bash
unalias -a
echo
echo "       path_change.bash (Version December 29 2024) (C) Copyright 2024 Arjan Koning All Rights Reserved"
echo
#
# Script to change the code directory in machine.f or machine.f90
#
# Set user
# set nuclear data library
#
user='Arjan Koning'
libname='iaea.2026'
#
#  Ensure that the code directory is changed by replacing
#    the path name in subroutine machine.f or machine.f90
#
codedir=`pwd`'/'
#
if [ -e machine.f ] ; then
  ffile=machine.f
fi
if [ -e machine.f90 ] ; then
  ffile=machine.f90
fi
if [ -e $ffile ] ; then
  sed "s/  user =.*/  user = '${user}'/" $ffile > tmp
  mv -f tmp $ffile
  echo "If not done, change the 'user' variable in the path_change.bash script by your own name "
  echo 'Right now it is: ' $user
  echo
  newdir=`echo $codedir | sed 's/\//\\\\\//g' | sed 's/source.*$//g'`
  sed "s/  codedir.*/  codedir = '${newdir}'/" $ffile > tmp
  mv -f tmp $ffile
  echo 'Replaced code directory in ' $ffile' by:'
  grep '  codedir =' $ffile
  echo 
  sed "s/  libname.*/  libname = '${libname}'/" $ffile > tmp
  mv -f tmp $ffile
  echo 'Replaced iaea library name in ' $ffile' by:' $libname
  echo
fi
