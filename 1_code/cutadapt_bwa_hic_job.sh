#!/usr/bin/env bash
#SBATCH --job-name=h9_hic_trim+bwa
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=vgl_a
#SBATCH --account=jarv_condo_bank
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32

set -euo pipefail

# -----------------------------
# Paths (H9 project)
# -----------------------------
WD="/lustre/fs5/vgl/scratch/amedico/h9_annotation"

RAW_DIR="${WD}/2_data/2.1_raw"
OUT_PROC="${WD}/2_data/2.2_processed"
OUT_TRIM="${OUT_PROC}/hic_cutadapt"
OUT_BWA="${OUT_PROC}/bwa_hic"

R1_IN="${RAW_DIR}/H9_p36_S2_R1_001.fastq.gz"
R2_IN="${RAW_DIR}/H9_p36_S2_R2_001.fastq.gz"

HAP1_REF="${RAW_DIR}/H9_T2T_v0.1_hap1.fasta"
HAP2_REF="${RAW_DIR}/H9_T2T_v0.1_hap2.fasta"

R1_TRIM="${OUT_TRIM}/H9_p36_S2.R1.cutadapt.trimmed.fastq.gz"
R2_TRIM="${OUT_TRIM}/H9_p36_S2.R2.cutadapt.trimmed.fastq.gz"

mkdir -p "${OUT_TRIM}" "${OUT_BWA}"

# -----------------------------
# Basic checks
# -----------------------------
for f in "${R1_IN}" "${R2_IN}" "${HAP1_REF}" "${HAP2_REF}"; do
  [[ -s "${f}" ]] || { echo "ERROR: missing/empty file: ${f}" >&2; exit 1; }
done
for exe in cutadapt bwa samtools; do
  command -v "${exe}" >/dev/null 2>&1 || { echo "ERROR: ${exe} not found in PATH" >&2; exit 1; }
done

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-32}"
THREADS_MAP=$(( THREADS_TOTAL / 2 ))
[[ "${THREADS_MAP}" -ge 1 ]] || THREADS_MAP=1
THREADS_SORT=$(( THREADS_MAP / 2 ))
[[ "${THREADS_SORT}" -ge 1 ]] || THREADS_SORT=1

echo "[INFO] Host: $(hostname)"
echo "[INFO] Job: ${SLURM_JOB_ID:-NA}"
echo "[INFO] Start: $(date)"
echo "[INFO] Threads total=${THREADS_TOTAL} map(each)=${THREADS_MAP} sort(each)=${THREADS_SORT}"

# -----------------------------
# 1) Cutadapt trim (parallel via -j)
# -----------------------------
echo "[INFO] Trimming Hi-C reads with cutadapt..."
cutadapt -j "${THREADS_TOTAL}" \
  --error-rate=0.1 \
  --times=1 \
  --overlap=3 \
  --action=trim \
  --cut=5 -U 5 \
  --minimum-length=1 \
  -o "${R1_TRIM}" \
  -p "${R2_TRIM}" \
  "${R1_IN}" "${R2_IN}"

# -----------------------------
# Helpers
# -----------------------------
need_bwa_index() {
  local ref="$1"
  local exts=(amb ann bwt pac sa)
  for e in "${exts[@]}"; do
    [[ -f "${ref}.${e}" ]] || return 0
  done
  return 1
}

index_ref() {
  local ref="$1"
  if need_bwa_index "${ref}"; then
    echo "[INFO] bwa index: ${ref}"
    bwa index "${ref}"
  else
    echo "[INFO] bwa index exists: ${ref}"
  fi
}

align_ref() {
  local ref="$1"
  local outprefix="$2"

  echo "[INFO] Mapping -> ${ref}"
  bwa mem -SP5M -M -t "${THREADS_MAP}" "${ref}" "${R1_TRIM}" "${R2_TRIM}" | \
    samtools view -bS - | \
    samtools sort -@ "${THREADS_SORT}" -o "${outprefix}.sorted.bam"

  samtools index -@ "${THREADS_SORT}" "${outprefix}.sorted.bam"
  samtools flagstat -@ "${THREADS_SORT}" "${outprefix}.sorted.bam" > "${outprefix}.flagstat.txt"
  echo "[INFO] Done: ${outprefix}.sorted.bam"
}

# -----------------------------
# 2) Index both haplotypes in parallel
# -----------------------------
index_ref "${HAP1_REF}" &
index_ref "${HAP2_REF}" &
wait

# -----------------------------
# 3) Align to hap1 and hap2 in parallel
# -----------------------------
align_ref "${HAP1_REF}" "${OUT_BWA}/H9_T2T_v0.1_hap1.hic" &
align_ref "${HAP2_REF}" "${OUT_BWA}/H9_T2T_v0.1_hap2.hic" &
wait

echo "[INFO] All done: $(date)"
