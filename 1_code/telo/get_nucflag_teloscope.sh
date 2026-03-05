#!/usr/bin/env bash
set -euo pipefail

# get_nucflag_teloscope.sh
# ------------------------
# Intersect H9 terminal telomere BED with NucFlag annotations
# (25 kbp chromosome ends evaluated by NucFlag) to produce
# per-base overlap of each telomere region with NucFlag categories.
#
# Output per dataset (HiFi / ONT):
#   {NUCFLAG_DIR}/nucflag_teloscope_{dataset}.bed
#   Columns: telo_chr  telo_start  telo_end  telo_length  arm
#            nf_chr  nf_start  nf_end  nf_name  nf_score  nf_strand
#            nf_thickStart  nf_thickEnd  nf_itemRgb  nf_zscore  nf_af
#            overlap_bp
#
# Usage: bash 1_code/telo/get_nucflag_teloscope.sh

# --- Configuration ---
NUCFLAG_DIR="2_data/2.2_processed/nucflag"
TELOSCOPE_DIR="2_data/2.2_processed/25.12.10_teloscope_compiled"
TERMINAL_BED="${TELOSCOPE_DIR}/H9_T2T_v0.1_dip.fasta_terminal_telomeres.bed"

DATASETS=("hifi" "ont")
NUCFLAG_FILES=("nucflag_telo_hifi.tsv" "nucflag_telo_ont.tsv")

# --- Preflight ---
[[ -d "$NUCFLAG_DIR" ]]    || { echo "Error: ${NUCFLAG_DIR} not found"  >&2; exit 1; }
[[ -f "$TERMINAL_BED" ]]   || { echo "Error: ${TERMINAL_BED} not found" >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found" >&2; exit 1; }

# --- Intersection ---
for i in "${!DATASETS[@]}"; do
    DATASET="${DATASETS[$i]}"
    NUCFLAG="${NUCFLAG_DIR}/${NUCFLAG_FILES[$i]}"
    OUT="${NUCFLAG_DIR}/nucflag_teloscope_${DATASET}.bed"

    [[ -f "$NUCFLAG" ]] || { echo "Error: ${NUCFLAG} not found" >&2; exit 1; }

    echo "=== ${DATASET}: intersecting telomeres with NucFlag ==="

    # 1. Extract telomere regions as BED
    #    Input has 10 cols: chr start end length arm fwdCounts revCounts canCounts nonCanCounts chrSize
    #    Output 5 cols:     chr start end length arm
    TMP_TELO=$(mktemp)
    awk -F'\t' 'BEGIN{OFS="\t"} { print $1, $2, $3, $4, $5 }' "$TERMINAL_BED" \
      | sort -k1,1 -k2,2n > "$TMP_TELO"

    N_TELO=$(wc -l < "$TMP_TELO")
    echo "  Telomere regions: ${N_TELO}"

    # 2. bedtools intersect -wo: report overlap (bp) for each pair
    #    Output: telo (5 cols) + nucflag (11 cols) + overlap_bp (1 col) = 17 cols
    #    NucFlag header (#chrom…) is stripped via grep before passing to bedtools
    bedtools intersect -a "$TMP_TELO" -b <(grep -v '^#' "$NUCFLAG") -wo > "$OUT"

    N_ROWS=$(wc -l < "$OUT")
    echo "  ${OUT} (${N_ROWS} intersection rows)"

    rm -f "$TMP_TELO"
done

echo ""
echo "Done."
