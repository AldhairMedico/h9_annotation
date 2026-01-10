#!/bin/bash
#SBATCH --job-name=aws_h9
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

# Working directory
WD=/lustre/fs5/vgl/scratch/amedico/h9_annotation
cd "$WD"

OUTDIR="/lustre/fs5/vgl/scratch/amedico/h9_annotation/2_data/2.1_raw"
LOGDIR="$WD/2_data/2.1_raw/_logs"
mkdir -p "$OUTDIR" "$LOGDIR"

# If aws is in a conda env, activate it here (edit if needed)
# module load mambaforge
# source activate awscli

# -------- robustness knobs --------
export AWS_MAX_ATTEMPTS=20
export AWS_RETRY_MODE=adaptive

# Multipart tuning: higher = more parallel parts (good on fat pipes)
# 128MB chunks is a sane HPC default; tune if needed.
aws configure set default.s3.multipart_threshold 128MB
aws configure set default.s3.multipart_chunksize 128MB
aws configure set default.s3.max_concurrent_requests 32
aws configure set default.s3.max_queue_size 10000

# List of URLs
URLS=(
  "https://genomeark.s3.amazonaws.com/species/Homo_sapiens/H9/genomic_data/pacbio_hifi/m84091_241205_180120_s2.hifi_reads.bc1003.fastq.gz"
  "https://genomeark.s3.amazonaws.com/species/Homo_sapiens/H9/genomic_data/pacbio_hifi/m84091_241209_213135_s3.hifi_reads.bc1003.fastq.gz"
  "https://genomeark.s3.amazonaws.com/species/Homo_sapiens/H9/genomic_data/pacbio_hifi/m84091_241209_233055_s4.hifi_reads.bc1003.fastq.gz"
  "https://genomeark.s3.amazonaws.com/species/Homo_sapiens/H9/genomic_data/ont/11_20_2024_R1041_UL_H9.dorado_0.8.1_sup.5mC_5hmC.fastq.gz"
  "https://genomeark.s3.amazonaws.com/species/Homo_sapiens/H9/genomic_data/ont/H9_ONT_SUL.fastq.gz"
  "https://genomeark.s3.amazonaws.com/species/Homo_sapiens/H9/genomic_data/ont/06_21_23_06_28_23_R1041_H9_UL.dorado_0.8.1_5mC_5hmC.fastq.gz"
)

# -------- helpers --------
url_to_s3() {
  # https://genomeark.s3.amazonaws.com/<key>  -> s3://genomeark/<key>
  local url="$1"
  local key="${url#https://genomeark.s3.amazonaws.com/}"
  echo "s3://genomeark/${key}"
}

ready_marker() {
  local f="$1"
  echo "${f}.done"
}

validate_gz() {
  local f="$1"
  # gzip integrity check (fast, no full decompression to disk)
  gzip -t "$f"
}

download_one() {
  local url="$1"
  local s3uri
  s3uri="$(url_to_s3 "$url")"

  local fname
  fname="$(basename "$url")"

  local out="$OUTDIR/$fname"
  local part="$out.part"
  local done
  done="$(ready_marker "$out")"

  if [[ -f "$done" && -s "$out" ]]; then
    echo "READY (already): $out"
    return 0
  fi

  echo "START: $s3uri -> $out"

  # Use a .part target so incomplete files are obvious.
  # aws s3 cp will resume only if file exists AND size matches expectation poorly; so we keep retries + re-run safety.
  # We also remove .part ONLY after validation.
  local attempt=0
  local max_attempts=10

  while (( attempt < max_attempts )); do
    attempt=$((attempt + 1))
    echo "Attempt $attempt/$max_attempts: $fname"

    # Download (multipart+parallel handled by aws settings above)
    # --no-sign-request works for public buckets (GenomeArk is public)
    aws s3 cp --no-sign-request "$s3uri" "$part" --only-show-errors || true

    # Basic checks: exists + non-empty
    if [[ -s "$part" ]]; then
      # Validate gzip integrity (only for .gz)
      if [[ "$part" == *.gz ]]; then
        if validate_gz "$part"; then
          mv -f "$part" "$out"
          date -Is > "$done"
          echo "READY: $out"
          return 0
        else
          echo "WARN: gzip test failed for $fname (will retry)"
        fi
      else
        mv -f "$part" "$out"
        date -Is > "$done"
        echo "READY: $out"
        return 0
      fi
    fi

    # If we got here, attempt failed or validation failed
    # Keep partial file for potential resume behavior and for debugging
    sleep $(( 30 * attempt ))
  done

  echo "ERROR: Failed to download $fname after $max_attempts attempts" >&2
  return 1
}

export -f url_to_s3 ready_marker validate_gz download_one
export OUTDIR LOGDIR

# Run all 6 concurrently, using up to 6 jobs (or fewer if you want to throttle)
# Each aws cp will itself use multipart concurrency; don’t crank this too high.
printf "%s\n" "${URLS[@]}" | \
  xargs -n 1 -P 6 -I {} bash -lc 'download_one "$@"' _ {}

echo "ALL DONE. Completed markers:"
ls -lh "$OUTDIR"/*.done 2>/dev/null || true
