#!/bin/bash
Thome=/Users/koning/
libdir=$Thome'libraries/'
for file in `ls -1 *iaea.med`; do
  proj=${file:0:1}
  nuc0=${file:2:6}
  nuc=`echo $nuc0 | sed 's/\-.*//'`
  tfile=`echo $file | sed 's/iaea\.med/tendl\.2023/'`
  tendlfile=$libdir$proj'/'$nuc'/tendl.2023/tables/residual/'$tfile
  if [ -e $tendlfile ] ; then
    cat > inp << EOI
$file
$tendlfile
EOI
    echo $file
    echo $tendlfile
    ffile=`echo $file | sed 's/iaea\.med/iaea\.2024/'`
    ~/isotopia/misc/medmerge < inp > $ffile
# mv inp inp-$proj$nuc
    rm inp
  fi
done
