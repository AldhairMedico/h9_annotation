#!/usr/bin/env bash
set -euo pipefail

# flip_E1_collate.sh
# -------------------
# Flips cooltools E1 eigenvector tracks per-chromosome based on the Spearman
# correlation with the coding_cov phasing track.  Chromosomes whose
# spearman_rho is negative get their E1 sign inverted (×-1), so that
# positive E1 → A (active) compartment across all chromosomes.
#
# Outputs per haplotype (dip, hap1, hap2):
#   {base}_{hap}.{res}.AB.flipped.bedgraph   bedGraph
#   {base}_{hap}.{res}.AB.flipped.bw   bigWig
#   {base}_{hap}.{res}.AB.flipped.bed  merged A/B compartment BED
#
# Collated outputs (hap1 + hap2 combined):
#   {base}_collated.{res}.AB.flipped.bedgraph / .bw / .bed
#
# Required inputs (per haplotype × resolution):
#   - cooltools vecs.tsv  (2_data/2.2_processed/)
#   - coding_cov summary  (3_figures/3.1_draft/26.02.05_hic/)

# ── Configuration ─────────────────────────────────────────────────────────
HAPS=("dip" "hap1" "hap2")
RESOLUTIONS=(5000 10000 20000 50000 100000 200000)
BASE="H9_T2T_v0.1"
METHOD="bwa"

# ── Directories ───────────────────────────────────────────────────────────
wd="/mnt/d/research/h9_annotation"
DATA_DIR="${wd}/2_data/2.2_processed"
FIG_DIR="${wd}/3_figures/3.1_draft/26.02.05_hic"
OUT_DIR="${DATA_DIR}/26.02.07_cooltools_flipped_E1"

# ── Helper: derive chrom.sizes from a bedGraph (preserves chrom order) ────
derive_sizes() {
  local bg="$1" sizes="$2"
  awk -F'\t' '{
    if (!($1 in s)) { order[++n] = $1 }
    if ($3+0 > s[$1]) s[$1] = $3+0
  } END {
    for (i = 1; i <= n; i++) print order[i], s[order[i]]
  }' "$bg" > "$sizes"
}

# ── Helper: generate merged A/B compartment BED from bedGraph ─────────────
#   A (E1 > 0) → #E38AAA (pink/salmon)
#   B (E1 ≤ 0) → #6F9DD0 (steel blue)
generate_ab_bed() {
  local bg="$1" bed="$2"
  awk 'BEGIN {
    OFS = "\t"
    last_chr = ""; last_start = 0; last_end = 0; last_comp = ""
  }
  {
    comp = ($4 > 0) ? "A" : "B"
    if ($1 == last_chr && comp == last_comp && $2+0 == last_end) {
      last_end = $3+0
    } else {
      if (last_comp != "") {
        color = (last_comp == "A") ? "#E38AAA" : "#6F9DD0"
        print last_chr, last_start, last_end, last_comp, 0, ".", last_start, last_end, color
      }
      last_chr = $1; last_start = $2+0; last_end = $3+0; last_comp = comp
    }
  }
  END {
    if (last_comp != "") {
      color = (last_comp == "A") ? "#E38AAA" : "#6F9DD0"
      print last_chr, last_start, last_end, last_comp, 0, ".", last_start, last_end, color
    }
  }' "$bg" > "$bed"
}

# ── Per-haplotype processing ──────────────────────────────────────────────
for HAP in "${HAPS[@]}"; do
  for RES in "${RESOLUTIONS[@]}"; do
    VECS="${DATA_DIR}/${BASE}_${HAP}.${METHOD}.${RES}.mcool.compartments.cis.vecs.tsv"
    STATS_TSV="${FIG_DIR}/${BASE}_${HAP}_${RES}/${BASE}_${HAP}_${RES}_coding_cov_summary_stats.tsv"

    OUT_BG="${OUT_DIR}/${BASE}_${HAP}.${RES}.AB.flipped.bedgraph"
    OUT_BW="${OUT_DIR}/${BASE}_${HAP}.${RES}.AB.flipped.bw"
    OUT_BED="${OUT_DIR}/${BASE}_${HAP}.${RES}.AB.flipped.bed"
    SIZES_TMP="${OUT_DIR}/.${BASE}_${HAP}.${RES}.chrom.sizes.tmp"

    echo "━━ ${HAP} @ ${RES} ━━"
    echo "  Stats : ${STATS_TSV}"
    echo "  Vecs  : ${VECS}"

    # 1) Build bedGraph: read summary stats for flip decisions, then
    #    process vecs.tsv with per-chromosome conditional sign inversion.
    awk -F'\t' '
    BEGIN { OFS = "\t" }
    # Pass 1 – summary stats: if spearman_rho < 0, flag chromosome for flip
    NR == FNR {
      if (FNR > 1 && $3 != "") { flip[$2] = ($3 + 0 < 0) ? 1 : 0 }
      next
    }
    # Pass 2 – vecs.tsv: emit bedGraph lines (skip header & empty/NaN E1)
    FNR > 1 && NF >= 5 && $5 != "" && $5 != "nan" && $5 != "NaN" {
      val = $5 + 0
      if (flip[$1]) val = -val
      print $1, $2, $3, val
    }
    ' "$STATS_TSV" "$VECS" > "$OUT_BG"

    # 2) Derive chrom.sizes from the bedGraph & convert → bigWig
    derive_sizes "$OUT_BG" "$SIZES_TMP"
    echo "  → ${OUT_BW}"
    bedGraphToBigWig "$OUT_BG" "$SIZES_TMP" "$OUT_BW"
    rm -f "$SIZES_TMP"

    # 3) Merged A/B compartment BED
    echo "  → ${OUT_BED}"
    generate_ab_bed "$OUT_BG" "$OUT_BED"
  done
done

# ── Collate hap1 + hap2 ──────────────────────────────────────────────────
for RES in "${RESOLUTIONS[@]}"; do
  HAP1_BG="${OUT_DIR}/${BASE}_hap1.${RES}.AB.flipped.bedgraph"
  HAP2_BG="${OUT_DIR}/${BASE}_hap2.${RES}.AB.flipped.bedgraph"

  COLL_BG="${OUT_DIR}/${BASE}_collated.${RES}.AB.flipped.bedgraph"
  COLL_BW="${OUT_DIR}/${BASE}_collated.${RES}.AB.flipped.bw"
  COLL_BED="${OUT_DIR}/${BASE}_collated.${RES}.AB.flipped.bed"
  COLL_SIZES="${OUT_DIR}/.${BASE}_collated.${RES}.chrom.sizes.tmp"

  echo "━━ collate @ ${RES} ━━"

  # 1) Concatenate hap1 + hap2 bedGraphs
  cat "$HAP1_BG" "$HAP2_BG" > "$COLL_BG"

  # 2) Derive chrom.sizes & convert → bigWig
  derive_sizes "$COLL_BG" "$COLL_SIZES"
  echo "  → ${COLL_BW}"
  bedGraphToBigWig "$COLL_BG" "$COLL_SIZES" "$COLL_BW"
  rm -f "$COLL_SIZES"

  # 3) Merged A/B compartment BED
  echo "  → ${COLL_BED}"
  generate_ab_bed "$COLL_BG" "$COLL_BED"
done

echo "Done. All tracks flipped, collated, and A/B BED generated."
