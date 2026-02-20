#!/bin/bash
#SBATCH --job-name=telogator2_h9
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
INPUT_DIR="${WD}/2_data/2.1_raw/hifi"
OUTPUT_DIR="${WD}/2_data/2.2_processed/telogator2_hifi_n3_fixed_feb19"
TELOGATOR2_DIR="/lustre/fs5/vgl/scratch/amedico/tools/telogator2"
TELOGATOR2_REF="${WD}/2_data/2.2_processed/H9_T2T_v0.1_dip.telogator2.ref.fixed.fa"

mkdir -p "$OUTPUT_DIR"
cd "$WD"

# ---- params ----
READ_TYPE="${READ_TYPE:-hifi}"
MIN_CLUSTER_READS="${MIN_CLUSTER_READS:-3}"

# ---- inputs (array-safe) ----
INPUTS=(
  "${INPUT_DIR}/m84091_241205_180120_s2.hifi_reads.bc1003.fastq.gz"
  "${INPUT_DIR}/m84091_241209_213135_s3.hifi_reads.bc1003.fastq.gz"
  "${INPUT_DIR}/m84091_241209_233055_s4.hifi_reads.bc1003.fastq.gz"
)

# ---- run ----
python "${TELOGATOR2_DIR}/telogator2.py" \
  -i "${INPUTS[@]}" \
  -o "$OUTPUT_DIR" \
  -t "$TELOGATOR2_REF"
  -r "$READ_TYPE" \
  -n "$MIN_CLUSTER_READS" \
  -p 32
