#!/bin/bash
#SBATCH --job-name=3d_annotation
#SBATCH --output=%x_%j_%a_%A.out
#SBATCH --error=%x_%j_%a_%A.err
#SBATCH --partition=vgl_a
#SBATCH --account=jarv_condo_bank
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32

# Activate environment with pairtools and cooltools installed
# conda activate cooltools_env

# Working directory
WD=/lustre/fs5/vgl/scratch/amedico/h9_annotation
cd $WD

# Which steps to run (0: chrom.sizes, 1–2: parse & sort, 4–11: downstream). Separated by spaces.
STEPS=(0 1 2 4 5 6 7 8a 8b 8c 9 10 11)
RESOLUTIONS=(20000 50000 100000)

# Helper: only run if $1 in STEPS
function run_step() {
  local s="$1"
  for step in "${STEPS[@]}"; do
    [[ "$step" == "$s" ]] && return 0
  done
  return 1
}

# Prefix and params (HAP: hap1, hap2, or dip)
HAP="hap2"
ASM="H9_T2T_v0.1_${HAP}"
ALIGNER="bwa"
PREFIX="${ASM}.${ALIGNER}"
WALKS_POLICY="mask"
BIN_SIZE=10000  # For cload

# Paths
RAW_DIR="${WD}/2_data/2.1_raw"
PROC_DIR="${WD}/2_data/2.2_processed"
BAM_DIR="${PROC_DIR}/bwa_hic"
ASM_FASTA="${RAW_DIR}/${ASM}.fasta"
INPUT_BAM="${BAM_DIR}/${ASM}.hic.sorted.bam"

echo "Processing: ${PREFIX}"

# 0) Generate chrom.sizes from assembly
if run_step 0; then
  gfastats ${ASM_FASTA} --path-report | tail -n +2 | cut -f1,4 > ${PROC_DIR}/${ASM}.chrom.sizes
  echo "Step 0: chrom.sizes done!"
fi

# 1) parse → ${PREFIX}.pairs.gz
# NOTE: Input BAM must be name-sorted. If coordinate-sorted, run:
  samtools sort -n -@ $SLURM_CPUS_PER_TASK ${INPUT_BAM} -o ${BAM_DIR}/${ASM}.hic.namesorted.bam
if run_step 1; then
  pairtools parse \
    --chroms-path   ${PROC_DIR}/${ASM}.chrom.sizes \
    --assembly      ${ASM} \
    --walks-policy  ${WALKS_POLICY} \
    --add-columns   mapq \
    --drop-sam --drop-seq \
    --output        ${PROC_DIR}/${PREFIX}.pairs.gz \
    --output-stats  ${PROC_DIR}/${PREFIX}.pairs.parse.stats.txt \
    --nproc-in      $SLURM_CPUS_PER_TASK \
    --nproc-out     $SLURM_CPUS_PER_TASK \
    ${BAM_DIR}/${ASM}.hic.namesorted.bam
  echo "Step 1: parse done!"
fi

# 2) sort → ${PREFIX}.pairs.sorted.gz
if run_step 2; then
  pairtools sort \
    --output        ${PROC_DIR}/${PREFIX}.pairs.sorted.gz \
    --nproc-in      $SLURM_CPUS_PER_TASK \
    --nproc-out     $SLURM_CPUS_PER_TASK \
    ${PROC_DIR}/${PREFIX}.pairs.gz
  echo "Step 2: sort done!"
fi


# 3) restrict → ${PREFIX}.pairs.sorted.restricted.gz (commented out)
#if run_step 3; then
#  pairtools restrict \
#    --frags       ${ASM}.arima.fragments.bed \
#    --output      ${PREFIX}.pairs.sorted.restricted.gz \
#    --nproc-in    $SLURM_CPUS_PER_TASK \
#    --nproc-out   $SLURM_CPUS_PER_TASK \
#    ${PREFIX}.pairs.sorted.gz
#  echo "Step 3: restrict done!"
#fi

# 4) dedup → ${PREFIX}.pairs.sorted.dedup.gz
if run_step 4; then
  pairtools dedup \
    --max-mismatch 0 \
    --mark-dups \
    --output        ${PROC_DIR}/${PREFIX}.pairs.sorted.dedup.gz \
    --output-stats  ${PROC_DIR}/${PREFIX}.pairs.sorted.dedup.stats.txt \
    --nproc-in      $SLURM_CPUS_PER_TASK \
    --nproc-out     $SLURM_CPUS_PER_TASK \
    ${PROC_DIR}/${PREFIX}.pairs.sorted.gz
  echo "Step 4: dedup done!"
fi

# 5) cload pairs → ${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
if run_step 5; then
  gunzip -c ${PROC_DIR}/${PREFIX}.pairs.sorted.dedup.gz | \
  cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    --assembly ${ASM} \
    ${PROC_DIR}/${ASM}.chrom.sizes:${BIN_SIZE} \
    - \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
  echo "Step 5: cload pairs done!"
fi

# 6) balance
if run_step 6; then
  cooler balance \
    -p $SLURM_CPUS_PER_TASK \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
  echo "Step 6: balance done!"
fi

# 7) zoomify → multi-resolution .mcool
if run_step 7; then
  cooler zoomify \
    --nproc $SLURM_CPUS_PER_TASK \
    --out ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.mcool \
    --resolutions 10000,20000,50000,100000,200000 \
    --balance \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
  echo "Step 7: zoomify done!"
fi

# 8a) Dump bin coords
if run_step 8a; then
  for RES in "${RESOLUTIONS[@]}"; do
    cooler dump \
      --join \
      --table bins \
      --header \
      ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.mcool::resolutions/${RES} \
      > ${PROC_DIR}/${PREFIX}.${RES}.bins.tsv
    echo "Step 8a: dump ${RES} bp bins done!"
  done
fi

# 8b) Compute GC per bin
if run_step 8b; then
  for RES in "${RESOLUTIONS[@]}"; do
    cooltools genome gc \
      ${PROC_DIR}/${PREFIX}.${RES}.bins.tsv \
      ${ASM_FASTA} \
      > ${PROC_DIR}/${PREFIX}.${RES}.gc.bedgraph
    echo "Step 8b: GC ${RES} bp bedgraph done!"
  done
fi

# 8c) Eigendecompose + phase to GC content
if run_step 8c; then
  for RES in "${RESOLUTIONS[@]}"; do
    cooltools eigs-cis \
      --out-prefix ${PROC_DIR}/${PREFIX}.${RES}.pairs.sorted.dedup.cool.compartments \
      --bigwig \
      --phasing-track ${PROC_DIR}/${PREFIX}.${RES}.gc.bedgraph \
      ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.mcool::resolutions/${RES}
    echo "Step 8c: eigs-cis ${RES} bp (gc) done!"
  done
fi


# 9) cis-expected (QC)
if run_step 9; then
  cooltools expected-cis \
    -p $SLURM_CPUS_PER_TASK \
    -o ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.expected_cis.tsv \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool
  echo "Step 9: expected-cis done!"
fi

# 10) insulation (TAD boundaries at 200 kb)
if run_step 10; then
  cooltools insulation \
    -p $SLURM_CPUS_PER_TASK \
    --bigwig \
    -o ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.insulation.tsv \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool \
    200000
  echo "Step 10: insulation done!"
fi

# 11) saddle plot (QC)
if run_step 11; then
  cooltools saddle \
    --out-prefix ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.saddle \
    --fig pdf \
    --fig svg \
    --fig png \
    --qrange 0.02 0.98 \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv::E1 \
    ${PROC_DIR}/${PREFIX}.${BIN_SIZE}.pairs.sorted.dedup.cool.expected_cis.tsv::balanced.avg
  echo "Step 11: saddle done!"
fi
