#!/bin/bash
for file in `ls -1`; do
  echo $file
  ~/isotopia/misc/medrename < $file > tmp
  mv tmp $file
done
