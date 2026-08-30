#!/usr/bin/env bash

set -euo pipefail

isotopia_dir=$(cd "$(dirname "$0")" && pwd)
source_dir="$isotopia_dir/source"

if [[ ! -d "$source_dir" ]]; then
  echo "ISOTOPIA installation error: source directory not found:" >&2
  echo "  $source_dir" >&2
  exit 1
fi

if [[ ! -f "$source_dir/Makefile" ]]; then
  echo "ISOTOPIA installation error: Makefile not found:" >&2
  echo "  $source_dir/Makefile" >&2
  exit 1
fi

data_file="$isotopia_dir/files/abundance/z001"

if [[ ! -f "$data_file" ]]; then
  echo "ISOTOPIA installation error: data files missing or incomplete:" >&2
  echo "  $data_file" >&2
  exit 1
fi

echo
echo "Installing ISOTOPIA"
echo "Installation directory: $isotopia_dir"
echo

make -C "$source_dir" clean
make -C "$source_dir" all "$@"

isotopia_exe="$isotopia_dir/bin/isotopia"

if [[ ! -x "$isotopia_exe" ]]; then
  echo "ISOTOPIA installation error: executable not created:" >&2
  echo "  $isotopia_exe" >&2
  exit 1
fi

echo
echo "ISOTOPIA executable:"
echo "  $isotopia_exe"
echo
echo "If not already done, add the following lines to your shell configuration:"
echo
echo "  export ISOTOPIA_DIR=\"$isotopia_dir\""
echo "  export PATH=\"\$ISOTOPIA_DIR/bin:\$PATH\""
echo "  export ISOTOPIA_USER=\"Your Name\""
echo
echo "By default, the cross-section database is expected at:"
echo "  $(dirname "$isotopia_dir")/isotopia.libs"
echo
echo "Alternatively, edit code_dir in source/machine.f90 and rebuild ISOTOPIA."
echo
