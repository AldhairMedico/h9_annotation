#!/bin/bash
set -euo pipefail

########################################
# User config
########################################

# Get script directory and derive repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Working directory containing teloscope outputs (flat structure)
WD="${REPO_ROOT}/data/2_processed"

# Tag for combined output files
RUN_TAG="asms_x1_TTAGGG_v1.3"

OUTPUT_BED="${WD}/${RUN_TAG}.terminal_telomeres.bed"
PATH_SUMMARY_OUT="${WD}/${RUN_TAG}.path_summary.tsv"
METRICS_OUT="${WD}/${RUN_TAG}.assembly_metrics.tsv"

########################################
# Collect per-assembly BED + report files,
# deduplicating with preference for ".chr"
########################################

declare -A BED_FOR_KEY      # canonical_id → bed path
declare -A LABEL_FOR_KEY    # canonical_id → label used in "asm" column
declare -A CHR_FLAG_FOR_KEY # canonical_id → 1 if ".chr" chosen, else 0

echo "[INFO] Searching for terminal_telomeres.bed files in ${WD}..."

if [[ ! -d "$WD" ]]; then
  echo "[ERROR] Working directory not found: $WD"
  exit 1
fi

while IFS= read -r bed_file; do
  bn=$(basename "$bed_file")
  # Strip trailing teloscope suffix
  bn_trim="${bn%_terminal_telomeres.bed}"   # e.g. GCA_...genomic.chr.fna
  # Canonical ID: same assembly with/without ".chr" collapse to the same key
  canonical="${bn_trim//.chr./.}"           # e.g. ...genomic.fna
  chr_flag=0
  [[ "$bn_trim" == *".chr."* ]] && chr_flag=1

  prev=${BED_FOR_KEY[$canonical]-}

  if [[ -z "$prev" ]]; then
    # First time we see this assembly
    BED_FOR_KEY["$canonical"]="$bed_file"
    LABEL_FOR_KEY["$canonical"]="$bn_trim"
    CHR_FLAG_FOR_KEY["$canonical"]=$chr_flag
  else
    # Already have a representation; prefer ".chr" over non-chr
    prev_chr=${CHR_FLAG_FOR_KEY[$canonical]:-0}
    if (( chr_flag > prev_chr )); then
      BED_FOR_KEY["$canonical"]="$bed_file"
      LABEL_FOR_KEY["$canonical"]="$bn_trim"
      CHR_FLAG_FOR_KEY["$canonical"]=$chr_flag
    fi
  fi
done < <(find "$WD" -maxdepth 1 -type f -name '*_terminal_telomeres.bed' | sort || true)

echo "[INFO] Found ${#BED_FOR_KEY[@]} unique assemblies after .chr dedup."

########################################
# Init combined outputs
########################################

: > "$OUTPUT_BED"
: > "$PATH_SUMMARY_OUT"
: > "$METRICS_OUT"

########################################
# Helper: sorted iteration over assemblies
########################################

keys=( "${!BED_FOR_KEY[@]}" )
IFS=$'\n' sorted_keys=( $(printf '%s\n' "${keys[@]}" | sort) )
unset IFS

########################################
# Build combined BED + reports
########################################

for key in "${sorted_keys[@]}"; do
  bed_file="${BED_FOR_KEY[$key]}"
  asm_label="${LABEL_FOR_KEY[$key]}"     # label used as "asm" / column name
  bed_dir=$(dirname "$bed_file")
  bed_bn=$(basename "$bed_file")
  asm_base="${bed_bn%_terminal_telomeres.bed}"   # base name before teloscope suffix
  report_file="${bed_dir}/${asm_base}.telo.report"

  echo "[INFO] Assembly: ${asm_label}"
  echo "       BED   : ${bed_file}"
  echo "       REPORT: ${report_file}"

  ######################################
  # 1) BED: vertical concat with asm label
  ######################################
  if [[ -f "$bed_file" ]]; then
    awk -v OFS='\t' -v asm="$asm_label" '{print $0, asm}' "$bed_file" >> "$OUTPUT_BED"
  else
    echo "[WARN] Missing BED for ${asm_label}: ${bed_file}"
  fi

  ######################################
  # 2) .telo.report: Path Summary (vertical) + metrics (horizontal)
  ######################################
  if [[ ! -f "$report_file" ]]; then
    echo "[WARN] Missing report for ${asm_label}: ${report_file}"
    continue
  fi

  # 2.1 Path Summary → vertical
  if [[ ! -s "$PATH_SUMMARY_OUT" ]]; then
    # First assembly: keep header, prepend "asm"
    awk -v asm="$asm_label" '
      BEGIN { in_sec = 0 }
      /^\+\+\+ Path Summary Report \+\+\+/ { in_sec = 1; next }
      /^\+\+\+/ && in_sec { in_sec = 0 }
      in_sec && NF > 0 {
        if ($1 == "pos") {
          print "asm\t" $0
        } else {
          print asm "\t" $0
        }
      }
    ' "$report_file" >> "$PATH_SUMMARY_OUT"
  else
    # Subsequent assemblies: skip header line ("pos ...")
    awk -v asm="$asm_label" '
      BEGIN { in_sec = 0 }
      /^\+\+\+ Path Summary Report \+\+\+/ { in_sec = 1; next }
      /^\+\+\+/ && in_sec { in_sec = 0 }
      in_sec && NF > 0 && $1 != "pos" {
        print asm "\t" $0
      }
    ' "$report_file" >> "$PATH_SUMMARY_OUT"
  fi

  # 2.2 Remaining sections → horizontal metrics
  metrics_tmp=$(mktemp)

  awk '
    BEGIN { seen_path = 0; metrics = 0 }
    # Mark that we have seen the Path Summary section
    /^\+\+\+ Path Summary Report \+\+\+/ { seen_path = 1; next }
    # First +++ after Path Summary marks start of the "metrics" block
    seen_path && /^\+\+\+/ && !metrics { metrics = 1; next }
    metrics {
      if (/^\+\+\+/) next      # skip further section headers
      if (NF == 0) next        # skip blank lines
      pos = index($0, ":")
      if (pos <= 0) next       # require "label: value"
      label = substr($0, 1, pos - 1)
      value = substr($0, pos + 1)
      gsub(/^[ \t]+/, "", value)
      print label "\t" value
    }
  ' "$report_file" > "$metrics_tmp"

  if [[ -s "$metrics_tmp" ]]; then
    if [[ ! -s "$METRICS_OUT" ]]; then
      # First assembly: create header + metric/value table
      {
        echo -e "metric\t${asm_label}"
        cat "$metrics_tmp"
      } > "$METRICS_OUT"
    else
      # Subsequent assemblies: paste new value column
      col_file=$(mktemp)
      {
        echo "${asm_label}"
        cut -f2 "$metrics_tmp"
      } > "$col_file"

      tmp_out=$(mktemp)
      paste "$METRICS_OUT" "$col_file" > "$tmp_out"
      mv "$tmp_out" "$METRICS_OUT"
      rm -f "$col_file"
    fi
  else
    echo "[WARN] No metrics parsed from ${asm_label}"
  fi

  rm -f "$metrics_tmp"
done

echo "[DONE] Combined BED          : ${OUTPUT_BED}"
echo "[DONE] Combined Path Summary : ${PATH_SUMMARY_OUT}"
echo "[DONE] Combined Metrics      : ${METRICS_OUT}"
