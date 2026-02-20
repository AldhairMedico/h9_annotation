#!/bin/bash
#SBATCH --job-name=topsicle_h9
#SBATCH --output=%x_%j_%a_%A.out
#SBATCH --error=%x_%j_%a_%A.err
#SBATCH --partition=vgl_a
#SBATCH --account=jarv_condo_bank
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32

set -euo pipefail

# ---- compact base ----
WD="/lustre/fs5/vgl/scratch/amedico/h9_annotation"
RAW_DIR="${WD}/2_data/2.1_raw"
INPUT_DIR="${RAW_DIR}/hifi"
OUTPUT_DIR="${WD}/2_data/2.2_processed/topsicle_hifi_r4k_t25k"

mkdir -p "$OUTPUT_DIR"
cd "$WD"

# ---- params (override with env vars if needed) ----
THREADS="${THREADS:-32}"

# Topsicle wants telomere repeat 5'->3'. For human it suggests CCCTAA.
# Override if you want TTAGGG:  PATTERN=TTAGGG sbatch run_topsicle_ont.sh
PATTERN="${PATTERN:-CCCTAA}"

MIN_SEQ_LEN="${MIN_SEQ_LEN:-4000}"
WINDOW_SIZE="${WINDOW_SIZE:-100}"
TRIMFIRST="${TRIMFIRST:-100}"
MAXLEN_TELO="${MAXLEN_TELO:-25000}"
CUTOFF="${CUTOFF:-0.7}"

# Optional flags
PLOT="${PLOT:-1}"         # 1 => --plot
OVERRIDE="${OVERRIDE:-1}" # 1 => --override

# ---- sanity ----
command -v topsicle >/dev/null 2>&1 || { echo "ERROR: topsicle not found in PATH"; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: input dir not found: $INPUT_DIR"; exit 1; }
ls -1 "$INPUT_DIR"/*.fastq.gz >/dev/null 2>&1 || { echo "ERROR: no *.fastq.gz in $INPUT_DIR"; exit 1; }

# ---- run (single call on folder) ----
cmd=(topsicle
  -i "$INPUT_DIR"
  -o "$OUTPUT_DIR"
  --pattern "$PATTERN"
  --minSeqLength "$MIN_SEQ_LEN"
  --cutoff "$CUTOFF"
  --windowSize "$WINDOW_SIZE"
  --trimfirst "$TRIMFIRST"
  --maxlengthtelo "$MAXLEN_TELO"
  -t "$THREADS"
)

(( PLOT == 1 )) && cmd+=(--plot)
(( OVERRIDE == 1 )) && cmd+=(--override)

echo "[topsicle] $(date) :: ${cmd[*]}"
"${cmd[@]}"
