#!/usr/bin/env bash
# Print the ISOTOPIA initial yield in the same GBq/C convention as TECDOC-1211.
set -euo pipefail

case_root=$(cd "$(dirname "$0")" && pwd)
isotopia_bin=${ISOTOPIA_BIN:-isotopia}
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

printf '%-20s %-10s %16s %18s %16s\n' case product 'ISOTOPIA GBq/C' 'TECDOC GBq/C' 'ISOTOPIA/TECDOC'
while IFS=$'\t' read -r case_name product table reference; do
  [[ $case_name == \#* || -z $case_name ]] && continue
  work="$scratch/$case_name"
  mkdir "$work"
  cp "$case_root/$case_name/isotopia.inp" "$work/"
  (cd "$work" && "$isotopia_bin" < isotopia.inp > isotopia.out)
  initial=$(awk '/Initial production rate \[GBq\/\(mA.h\)\]/{print $NF; exit}' "$work/$product.act")
  yield_per_c=$(awk -v value="$initial" 'BEGIN {printf "%.8g", value / 3.6}')
  ratio=''
  if [[ $reference =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    ratio=$(awk -v isotope_yield="$yield_per_c" -v tecdoc_yield="$reference" 'BEGIN {printf "%.6f", isotope_yield / tecdoc_yield}')
  fi
  printf '%-20s %-10s %16s %18s %16s  (%s)\n' "$case_name" "$product" "$yield_per_c" "$reference" "$ratio" "$table"
done < "$case_root/cases.tsv"
