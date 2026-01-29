#!/bin/bash
#SBATCH --job-name=telogator2_h9
#SBATCH --output=%x_%j_%a_%A.out
#SBATCH --error=%x_%j_%a_%A.err
#SBATCH --partition=vgl_c
#SBATCH --account=jarv_condo_bank
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32

set -euo pipefail

# ---- compact base ----
WD="/lustre/fs5/vgl/scratch/amedico/h9_annotation"
INPUT_DIR="${WD}/2_data/2.1_raw/ont"
OUTPUT_DIR="${WD}/2_data/2.2_processed/telogator2_ont_n3"
TELOGATOR2_DIR="/lustre/fs5/vgl/scratch/amedico/tools/telogator2"

mkdir -p "$OUTPUT_DIR"
cd "$WD"

# ---- params ----
READ_TYPE="${READ_TYPE:-ont}"
MIN_CLUSTER_READS="${MIN_CLUSTER_READS:-3}"

# ---- inputs (array-safe) ----
INPUTS=(
  "${INPUT_DIR}/06_21_23_06_28_23_R1041_H9_UL.dorado_0.8.1_5mC_5hmC.fastq.gz"
  "${INPUT_DIR}/11_20_2024_R1041_UL_H9.dorado_0.8.1_sup.5mC_5hmC.fastq.gz"
)

# ---- run ----
python "${TELOGATOR2_DIR}/telogator2.py" \
  -i "${INPUTS[@]}" \
  -o "$OUTPUT_DIR" \
  -r "$READ_TYPE" \
  -n "$MIN_CLUSTER_READS" \
  -p 32
