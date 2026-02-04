#!/usr/bin/env bash
set -euo pipefail

# get_gene_coverage.sh
# --------------------
# 1. Extract protein-coding gene BED from H9-T2T GTF annotations
#    (GTF has transcript/exon/CDS features only — no "gene" rows)
# 2. Compute coding coverage per Hi-C bin using bedtools
#
# Output per haplotype + resolution:
#   {BINS_DIR}/{PREFIX}.{RES}.coding_cov.bedgraph
#   Columns: chrom  start  end  weight  coding
#
# Usage: bash 1_code/hic/get_gene_coverage.sh

# --- Configuration ---
ANNOT_DIR="2_data/2.2_processed/25.12.16_gene_annotations"
BINS_DIR="2_data/2.2_processed"
HAPLOTYPES=("HAP1" "HAP2")

# --- Preflight ---
[[ -d "$ANNOT_DIR" ]] || { echo "Error: ${ANNOT_DIR} not found" >&2; exit 1; }
command -v awk     >/dev/null 2>&1 || { echo "Error: awk not found"      >&2; exit 1; }
command -v gunzip  >/dev/null 2>&1 || { echo "Error: gunzip not found"   >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found" >&2; exit 1; }

# --- Step 1: Extract coding gene BEDs ---
declare -A CODING_BEDS

for HAP in "${HAPLOTYPES[@]}"; do
    GTF_GZ="${ANNOT_DIR}/H9_${HAP}.gtf.gz"
    CODING_BED="${ANNOT_DIR}/H9_${HAP}.genes.coding.bed"

    [[ -f "$GTF_GZ" ]] || { echo "Error: ${GTF_GZ} not found" >&2; exit 1; }

    echo "=== ${HAP}: extracting coding genes ==="

    # Filter protein_coding transcripts, collapse by gene_id → gene-level BED
    # (GTF is 1-based; subtract 1 from start for 0-based BED)
    gunzip -c "$GTF_GZ" \
      | awk -F'\t' -v OFS='\t' '
          $3 == "transcript" {
              bt = $9; sub(/.*gene_biotype "/, "", bt); sub(/".*/, "", bt)
              if (bt != "protein_coding") next
              gid = $9; sub(/.*gene_id "/, "", gid); sub(/".*/, "", gid)
              key = $1 SUBSEP gid
              if (!(key in starts) || $4 < starts[key]) { starts[key] = $4; chroms[key] = $1 }
              if (!(key in ends)   || $5 > ends[key])   ends[key] = $5
          }
          END {
              for (k in starts) print chroms[k], starts[k]-1, ends[k]
          }
        ' | sort -k1,1 -k2,2n > "$CODING_BED"

    NGENES=$(wc -l < "$CODING_BED")
    echo "  ${CODING_BED} (${NGENES} coding genes)"

    CODING_BEDS["${HAP}"]="$CODING_BED"
done

# --- Step 2: Compute coding coverage per bin ---
echo ""
echo "=== Computing coding coverage per bin ==="

for TSV in "${BINS_DIR}"/H9_T2T_v0.1_hap*.bwa.*.bins.tsv; do
    [[ -f "$TSV" ]] || continue

    BASENAME=$(basename "$TSV")
    BASE="${BASENAME%.bins.tsv}"
    OUT="${BINS_DIR}/${BASE}.coding_cov.bedgraph"

    # Determine haplotype from filename
    if   [[ "$BASENAME" == *hap1* ]]; then HAP="HAP1"
    elif [[ "$BASENAME" == *hap2* ]]; then HAP="HAP2"
    else echo "  Skipping ${BASENAME}: unknown haplotype" >&2; continue
    fi

    CODING_BED="${CODING_BEDS[$HAP]}"
    echo "  ${BASENAME} → $(basename "$OUT")  [${HAP}]"

    # Headerless BED from bins
    TMP_BED=$(mktemp)
    tail -n +2 "$TSV" | cut -f1-3 > "$TMP_BED"

    # bedtools coverage: col 7 = fraction of bin covered by coding genes
    TMP_COV=$(mktemp)
    bedtools coverage -a "$TMP_BED" -b "$CODING_BED" > "$TMP_COV"

    # Extract coding fraction
    TMP_FRAC=$(mktemp)
    awk '{ print $7 }' "$TMP_COV" > "$TMP_FRAC"

    # Original data columns (skip header)
    TMP_DATA=$(mktemp)
    tail -n +2 "$TSV" | cut -f1-4 > "$TMP_DATA"

    # Write output
    printf 'chrom\tstart\tend\tweight\tcoding\n' > "$OUT"
    paste "$TMP_DATA" "$TMP_FRAC" >> "$OUT"

    rm -f "$TMP_BED" "$TMP_COV" "$TMP_FRAC" "$TMP_DATA"
done

echo ""
echo "Done."