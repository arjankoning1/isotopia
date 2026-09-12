#!/usr/bin/env bash
# Compare the EOI Mo-99 activity with the 0 mm Ta total in TECDOC-2051 Table 2.
set -euo pipefail

case_root=$(cd "$(dirname "$0")" && pwd)
isotopia_bin=${ISOTOPIA_BIN:-isotopia}
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

printf '%-24s %-10s %18s %18s %16s\n' case product 'ISOTOPIA EOI [GBq]' 'TECDOC EOI [GBq]' 'ISOTOPIA/TECDOC'
while IFS=$'\t' read -r case_name product reference scope; do
  [[ $case_name == \#* || -z $case_name ]] && continue
  work="$scratch/$case_name"
  mkdir "$work"
  cp "$case_root/$case_name/isotopia.inp" "$work/"
  (cd "$work" && "$isotopia_bin" < isotopia.inp > isotopia.out)
  activity=$(awk '/Total activity at EOI \[GBq\]/{print $NF; exit}' "$work/$product.act")
  ratio=$(awk -v isotope_activity="$activity" -v tecdoc_activity="$reference" 'BEGIN {printf "%.6f", isotope_activity / tecdoc_activity}')
  printf '%-24s %-10s %18.8g %18s %16s  (%s)\n' "$case_name" "$product" "$activity" "$reference" "$ratio" "$scope"
done < "$case_root/cases.tsv"
