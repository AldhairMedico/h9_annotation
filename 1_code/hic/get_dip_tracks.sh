#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Build diploid bedgraph tracks by merging hap1 + hap2 while keeping ONE header.
# Script location: h9_annotation/1_code/hic/get_dip_tracks.sh
# Looks for inputs in: h9_annotation/2_data/2.2_processed
#
# File pattern:
#   ASM_HAP.ALN.RES.TRACKNAME.bedgraph
# Example:
#   H9_T2T_v0.1_hap1.bwa.5000.coding_cov.bedgraph
# Output:
#   H9_T2T_v0.1_dip.bwa.5000.coding_cov.bedgraph
# -----------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/../../2_data/2.2_processed"

ASM="${ASM:-H9_T2T_v0.1}"
ALN="${ALN:-bwa}"
TRACK="${TRACK:-coding_cov}"

HAPS=(hap1 hap2)
RESOLUTIONS=(5000 10000 20000 50000 100000 200000)

die() { echo "[ERROR] $*" >&2; exit 1; }

[[ -d "${DATA_DIR}" ]] || die "Data directory not found: ${DATA_DIR}"

for res in "${RESOLUTIONS[@]}"; do
  f1="${DATA_DIR}/${ASM}_${HAPS[0]}.${ALN}.${res}.${TRACK}.bedgraph"
  f2="${DATA_DIR}/${ASM}_${HAPS[1]}.${ALN}.${res}.${TRACK}.bedgraph"
  out="${DATA_DIR}/${ASM}_dip.${ALN}.${res}.${TRACK}.bedgraph"

  [[ -s "${f1}" ]] || die "Missing/empty input: ${f1}"
  [[ -s "${f2}" ]] || die "Missing/empty input: ${f2}"

  h1="$(head -n 1 "${f1}")"
  h2="$(head -n 1 "${f2}")"
  if [[ "${h1}" != "${h2}" ]]; then
    echo "[WARN] Header mismatch for res=${res} track=${TRACK}; keeping hap1 header" >&2
  fi

  tmp="$(mktemp "${out}.tmp.XXXXXX")"
  {
    printf '%s\n' "${h1}"
    tail -n +2 "${f1}"
    tail -n +2 "${f2}"
  } > "${tmp}"

  mv -f "${tmp}" "${out}"
  echo "[OK] Wrote: ${out}"
done
