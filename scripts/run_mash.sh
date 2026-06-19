#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: run_mash.sh <genome_list.txt> <output_prefix>" >&2
  exit 2
fi
LIST="$1"
PREFIX="$2"
if ! command -v mash >/dev/null 2>&1; then
  echo "ERROR: mash is not available on PATH" >&2
  exit 127
fi
mash sketch -k 21 -s 10000 -p 8 -l "$LIST" -o "$PREFIX"
mash dist "${PREFIX}.msh" "${PREFIX}.msh" > "${PREFIX}.dist.tsv"
