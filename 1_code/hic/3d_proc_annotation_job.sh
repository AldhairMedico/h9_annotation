#!/bin/bash
#SBATCH --job-name=3d_proc_annotation
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

# Activate environment with pairtools and cooltools installed
# conda activate cooltools_env

# Working directory: SLURM_SUBMIT_DIR on HPC, script-relative locally
WD=$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}/../.." && pwd)
cd "$WD"

# Which steps to run (1-3: preprocessing, 4-15: downstream). Separated by spaces.
STEPS=(2 3 4 5 6 7 8 9 10 11 12a 12b 12c 13 14 15)
RESOLUTIONS=(5000 10000 20000 50000 100000 200000)

# Force overwrite existing files (set to true to rerun steps even if output exists)
FORCE_OVERWRITE=true

# Helper: only run if $1 in STEPS
function run_step() {
  local s="$1"
  for step in "${STEPS[@]}"; do
    [[ "$step" == "$s" ]] && return 0
  done
  return 1
}

# Helper: check if file exists (skip step if exists and not forcing overwrite)
function should_skip() {
  local file="$1"
  if [[ -s "${file}" ]] && [[ "${FORCE_OVERWRITE}" == "false" ]]; then
    return 0
  fi
  return 1
}

# Prefix and params (HAP: hap1, hap2, or dip)
HAP="dip"
ASM="H9_T2T_v0.1_${HAP}"
ALIGNER="bwa"
PREFIX="${ASM}.${ALIGNER}"
WALKS_POLICY="mask"
BIN_SIZE=5000  # For cload (must be <= finest RESOLUTIONS entry)

# Derive comma-separated resolution list for cooler zoomify (single source of truth)
ZOOMIFY_RES_CSV=$(IFS=,; echo "${RESOLUTIONS[*]}")

# Validate: every RESOLUTIONS entry must be >= BIN_SIZE (zoomify can't go finer than input)
for _r in "${RESOLUTIONS[@]}"; do
  (( _r >= BIN_SIZE )) || {
    echo "ERROR: resolution ${_r} < BIN_SIZE ${BIN_SIZE}; cload base is too coarse" >&2
    exit 1
  }
done

# Paths
RAW_DIR="${WD}/2_data/2.1_raw"
PROC_DIR="${WD}/2_data/2.2_processed"
TRIM_DIR="${PROC_DIR}/hic_cutadapt"
BAM_DIR="${PROC_DIR}/bwa_hic"
ASM_FASTA="${RAW_DIR}/${ASM}.fasta"

# Input reads
R1_IN="${RAW_DIR}/H9_p36_S2_R1_001.fastq.gz"
R2_IN="${RAW_DIR}/H9_p36_S2_R2_001.fastq.gz"

# Trimmed reads
R1_TRIM="${TRIM_DIR}/H9_p36_S2.R1.cutadapt.trimmed.fastq.gz"
R2_TRIM="${TRIM_DIR}/H9_p36_S2.R2.cutadapt.trimmed.fastq.gz"

mkdir -p "${TRIM_DIR}" "${BAM_DIR}"

# Basic checks
for f in "${R1_IN}" "${R2_IN}" "${ASM_FASTA}"; do
  [[ -s "${f}" ]] || { echo "ERROR: missing/empty file: ${f}" >&2; exit 1; }
done
#for exe in cutadapt bwa samtools gfastats pairtools cooler cooltools; do
#  command -v "${exe}" >/dev/null 2>&1 || { echo "ERROR: ${exe} not found in PATH" >&2; exit 1; }
#done

# Thread allocation
THREADS_TOTAL=32
THREADS_MAP=32
THREADS_SORT=32

echo "[INFO] Host: $(hostname)"
echo "[INFO] Job: ${SLURM_JOB_ID:-NA}"
echo "[INFO] Start: $(date)"
echo "[INFO] Processing: ${PREFIX}"
echo "[INFO] Threads total=${THREADS_TOTAL} map=${THREADS_MAP} sort=${THREADS_SORT}"
echo "[INFO] FORCE_OVERWRITE=${FORCE_OVERWRITE}"

# =============================================================================
# PREPROCESSING STEPS (1-3)
# =============================================================================

# 1) Cutadapt trim
if run_step 1; then
  if should_skip "${R1_TRIM}" && should_skip "${R2_TRIM}"; then
    echo "Step 1: trimmed reads exist, skipping!"
  else
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
    echo "Step 1: cutadapt done!"
  fi
fi

# 2) BWA index
if run_step 2; then
  INDEX_EXISTS=true
  for ext in amb ann bwt pac sa; do
    [[ -f "${ASM_FASTA}.${ext}" ]] || INDEX_EXISTS=false
  done
  if [[ "${INDEX_EXISTS}" == "true" ]] && [[ "${FORCE_OVERWRITE}" == "false" ]]; then
    echo "Step 2: bwa index exists, skipping!"
  else
    echo "[INFO] bwa index: ${ASM_FASTA}"
    bwa index "${ASM_FASTA}"
    echo "Step 2: bwa index done!"
  fi
fi

# 3) BWA mem + sort + index + flagstat
if run_step 3; then
  SORTED_BAM="${BAM_DIR}/${ASM}.hic.sorted.bam"
  if should_skip "${SORTED_BAM}"; then
    echo "Step 3: sorted BAM exists, skipping!"
  else
    echo "[INFO] Mapping -> ${ASM_FASTA}"
    bwa mem -SP5M -M -t "${THREADS_MAP}" "${ASM_FASTA}" "${R1_TRIM}" "${R2_TRIM}" | \
      samtools view -bS - | \
      samtools sort -@ "${THREADS_SORT}" -o "${SORTED_BAM}"
    samtools index -@ "${THREADS_SORT}" "${SORTED_BAM}"
    samtools flagstat -@ "${THREADS_SORT}" "${SORTED_BAM}" > "${BAM_DIR}/${ASM}.hic.flagstat.txt"
    echo "Step 3: bwa mem + sort done!"
  fi
fi

# =============================================================================
# 3D ANNOTATION STEPS (4-15)
# =============================================================================

# 4) Generate chrom.sizes from assembly
if run_step 4; then
  CHROM_SIZES="${PROC_DIR}/${ASM}.chrom.sizes"
  if should_skip "${CHROM_SIZES}"; then
    echo "Step 4: chrom.sizes exists, skipping!"
  else
    samtools faidx "${ASM_FASTA}"
    cut -f1,2 "${ASM_FASTA}.fai" > "${CHROM_SIZES}"
    echo "Step 4: chrom.sizes done!"
  fi
fi

# 5) Name-sort BAM for pairtools
if run_step 5; then
  NAMESORTED_BAM="${BAM_DIR}/${ASM}.hic.namesorted.bam"
  SORTED_BAM="${BAM_DIR}/${ASM}.hic.sorted.bam"
  if should_skip "${NAMESORTED_BAM}"; then
    echo "Step 5: namesorted BAM exists, skipping!"
  else
    samtools sort -n -@ "${THREADS_SORT}" ${SORTED_BAM} -o ${NAMESORTED_BAM}
    echo "Step 5: namesort done!"
  fi
fi

# 6) parse → ${PREFIX}.pairs.gz
if run_step 6; then
  PAIRS_GZ="${PROC_DIR}/${PREFIX}.pairs.gz"
  if should_skip "${PAIRS_GZ}"; then
    echo "Step 6: pairs.gz exists, skipping!"
  else
    pairtools parse \
      --chroms-path   ${PROC_DIR}/${ASM}.chrom.sizes \
      --assembly      ${ASM} \
      --walks-policy  ${WALKS_POLICY} \
      --add-columns   mapq \
      --drop-sam --drop-seq \
      --output        ${PAIRS_GZ} \
      --output-stats  ${PROC_DIR}/${PREFIX}.pairs.parse.stats.txt \
      --nproc-in      $SLURM_CPUS_PER_TASK \
      --nproc-out     $SLURM_CPUS_PER_TASK \
      ${BAM_DIR}/${ASM}.hic.namesorted.bam
    echo "Step 6: parse done!"
  fi
fi

# 7) sort → ${PREFIX}.pairs.sorted.gz
if run_step 7; then
  PAIRS_SORTED="${PROC_DIR}/${PREFIX}.pairs.sorted.gz"
  if should_skip "${PAIRS_SORTED}"; then
    echo "Step 7: pairs.sorted.gz exists, skipping!"
  else
    pairtools sort \
      --output        ${PAIRS_SORTED} \
      --nproc-in      $SLURM_CPUS_PER_TASK \
      --nproc-out     $SLURM_CPUS_PER_TASK \
      ${PROC_DIR}/${PREFIX}.pairs.gz
    echo "Step 7: sort done!"
  fi
fi

# 8) dedup → ${PREFIX}.pairs.sorted.dedup.gz
if run_step 8; then
  PAIRS_DEDUP="${PROC_DIR}/${PREFIX}.pairs.sorted.dedup.gz"
  if should_skip "${PAIRS_DEDUP}"; then
    echo "Step 8: pairs.sorted.dedup.gz exists, skipping!"
  else
    pairtools dedup \
      --max-mismatch 0 \
      --mark-dups \
      --output        ${PAIRS_DEDUP} \
      --output-stats  ${PROC_DIR}/${PREFIX}.pairs.sorted.dedup.stats.txt \
      --nproc-in      $SLURM_CPUS_PER_TASK \
      --nproc-out     $SLURM_CPUS_PER_TASK \
      ${PROC_DIR}/${PREFIX}.pairs.sorted.gz
    echo "Step 8: dedup done!"
  fi
fi

# 9) cload pairs → ${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
if run_step 9; then
  COOL_FILE="${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool"
  if should_skip "${COOL_FILE}"; then
    echo "Step 9: cool file exists, skipping!"
  else
    gunzip -c ${PROC_DIR}/${PREFIX}.pairs.sorted.dedup.gz | \
    cooler cload pairs \
      -c1 2 -p1 3 -c2 4 -p2 5 \
      --assembly ${ASM} \
      ${PROC_DIR}/${ASM}.chrom.sizes:${BIN_SIZE} \
      - \
      ${COOL_FILE}
    echo "Step 9: cload pairs done!"
  fi
fi

# 10) balance
if run_step 10; then
  COOL_FILE="${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool"
  cooler balance \
    -p $SLURM_CPUS_PER_TASK \
    ${COOL_FILE}
  echo "Step 10: balance done!"
fi

# 11) zoomify → multi-resolution .mcool
if run_step 11; then
  MCOOL_FILE="${PROC_DIR}/${PREFIX}.${BIN_SIZE}.mcool"
  if should_skip "${MCOOL_FILE}"; then
    echo "Step 11: mcool file exists, skipping!"
  else
    cooler zoomify \
      --nproc $SLURM_CPUS_PER_TASK \
      --out ${MCOOL_FILE} \
      --resolutions ${ZOOMIFY_RES_CSV} \
      --balance \
      ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
    echo "Step 11: zoomify done!"
  fi
fi

# 12a) Dump bin coords
if run_step 12a; then
  for RES in "${RESOLUTIONS[@]}"; do
    BINS_TSV="${PROC_DIR}/${PREFIX}.${RES}.bins.tsv"
    if should_skip "${BINS_TSV}"; then
      echo "Step 12a: ${RES} bp bins.tsv exists, skipping!"
    else
      cooler dump \
        --join \
        --table bins \
        --header \
        ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.mcool::resolutions/${RES} \
        > ${BINS_TSV}
      echo "Step 12a: dump ${RES} bp bins done!"
    fi
  done
fi

# 12b) Compute GC per bin
if run_step 12b; then
  for RES in "${RESOLUTIONS[@]}"; do
    GC_BEDGRAPH="${PROC_DIR}/${PREFIX}.${RES}.gc.bedgraph"
    if should_skip "${GC_BEDGRAPH}"; then
      echo "Step 12b: ${RES} bp gc.bedgraph exists, skipping!"
    else
      cooltools genome gc \
        ${PROC_DIR}/${PREFIX}.${RES}.bins.tsv \
        ${ASM_FASTA} \
        > ${GC_BEDGRAPH}
      echo "Step 12b: GC ${RES} bp bedgraph done!"
    fi
  done
fi

# 12c) Eigendecompose + phase to GC content
if run_step 12c; then
  for RES in "${RESOLUTIONS[@]}"; do
    COMPARTMENTS_PREFIX="${PROC_DIR}/${PREFIX}.${RES}.pairs.sorted.dedup.cool.compartments"
    if should_skip "${COMPARTMENTS_PREFIX}.cis.vecs.tsv"; then
      echo "Step 12c: ${RES} bp compartments exists, skipping!"
    else
      cooltools eigs-cis \
        --out-prefix ${COMPARTMENTS_PREFIX} \
        --bigwig \
        --phasing-track ${PROC_DIR}/${PREFIX}.${RES}.gc.bedgraph \
        ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.mcool::resolutions/${RES}
      echo "Step 12c: eigs-cis ${RES} bp (gc) done!"
    fi
  done
fi

# 13) cis-expected (QC)
if run_step 13; then
  EXPECTED_CIS="${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.expected_cis.tsv"
  if should_skip "${EXPECTED_CIS}"; then
    echo "Step 13: expected_cis.tsv exists, skipping!"
  else
    cooltools expected-cis \
      -p $SLURM_CPUS_PER_TASK \
      -o ${EXPECTED_CIS} \
      ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
    echo "Step 13: expected-cis done!"
  fi
fi

# 14) insulation (TAD boundaries at 200 kb)
if run_step 14; then
  INSULATION="${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.insulation.tsv"
  if should_skip "${INSULATION}"; then
    echo "Step 14: insulation.tsv exists, skipping!"
  else
    cooltools insulation \
      -p $SLURM_CPUS_PER_TASK \
      --bigwig \
      -o ${INSULATION} \
      ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool \
      200000
    echo "Step 14: insulation done!"
  fi
fi

# 15) saddle plot (QC)
if run_step 15; then
  SADDLE_PREFIX="${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.saddle"
  if should_skip "${SADDLE_PREFIX}.pdf"; then
    echo "Step 15: saddle plot exists, skipping!"
  else
    cooltools saddle \
      --out-prefix ${SADDLE_PREFIX} \
      --fig pdf \
      --fig svg \
      --fig png \
      --qrange 0.02 0.98 \
      ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool \
      ${PROC_DIR}/${PREFIX}.${RESOLUTIONS[0]}.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv::E1 \
      ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.expected_cis.tsv::balanced.avg
    echo "Step 15: saddle done!"
  fi
fi

echo "[INFO] All done: $(date)"
