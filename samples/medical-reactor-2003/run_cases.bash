#!/usr/bin/env bash
# Compare final requested activities with the operational values in TECDOC-1340.
set -euo pipefail

case_root=$(cd "$(dirname "$0")" && pwd)
isotopia_bin=${ISOTOPIA_BIN:-isotopia}
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

printf '%-18s %-10s %18s %18s %16s\n' case product 'ISOTOPIA [GBq]' 'TECDOC [GBq]' 'ISOTOPIA/TECDOC'
while IFS=$'\t' read -r case_name product reference scope; do
  [[ $case_name == \#* || -z $case_name ]] && continue
  work="$scratch/$case_name"
  mkdir "$work"
  cp "$case_root/$case_name/isotopia.inp" "$work/"
  (cd "$work" && "$isotopia_bin" < isotopia.inp > isotopia.out)
  activity=$(awk '$1 ~ /^[0-9.]+E[+-][0-9]+$/ {value = $2} END {print value}' "$work/$product.act")
  ratio=$(awk -v isotope_activity="$activity" -v tecdoc_activity="$reference" 'BEGIN {printf "%.6f", isotope_activity / tecdoc_activity}')
  printf '%-18s %-10s %18.8g %18s %16s  (%s)\n' "$case_name" "$product" "$activity" "$reference" "$ratio" "$scope"
done < "$case_root/cases.tsv"
