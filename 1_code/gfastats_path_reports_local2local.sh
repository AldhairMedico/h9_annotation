#!/bin/bash
set -e

# Paths
wd="/mnt/d/research/h9_annotation"
in_dir="$wd/2_data/2.1_raw"
out_dir="$wd/2_data/2.2_processed/gfastats"

mkdir -p "$out_dir"

# Function to process a single assembly
# Runs 3 gfastats jobs in parallel (1 thread each)
process_assembly() {
  local input_file="$1"
  local basename=$(basename "$input_file")

  # Extract assembly name (remove .gz if present, then remove .fasta/.fna/.fa extension)
  local asm_name="${basename%.gz}"
  asm_name="${asm_name%.fasta}"
  asm_name="${asm_name%.fna}"
  asm_name="${asm_name%.fa}"

  local asm_out_dir="$out_dir/$asm_name"

  # Check if all expected outputs already exist - skip if so
  if [[ -f "$asm_out_dir/gfastats_${asm_name}.tsv" ]] && \
     [[ -f "$asm_out_dir/gfastatsNxContig_${asm_name}.tsv" ]] && \
     [[ -f "$asm_out_dir/gfastatsNxScaffold_${asm_name}.tsv" ]]; then
    echo "Skipping (outputs exist): $asm_name"
    return 0
  fi

  mkdir -p "$asm_out_dir"

  echo "Processing: $asm_name"

  # Run 3 gfastats jobs in parallel using background processes
  # Job 1: Basic stats
  (
    printf "$asm_name\t" > "$asm_out_dir/gfastats_${asm_name}.tsv"
    gfastats -t "$input_file" -j 1 | cut -f2 | sed -z 's/\n/\t/g; s/.$//' >> "$asm_out_dir/gfastats_${asm_name}.tsv"
    echo "" >> "$asm_out_dir/gfastats_${asm_name}.tsv"
  ) &
  pid1=$!

  # Job 2: Nx contig stats
  (
    printf "$asm_name\t" > "$asm_out_dir/gfastatsNxContig_${asm_name}.tsv"
    gfastats "$input_file" -j 1 -s c | sort -nrk2 | awk 'BEGIN{pos=0}{total+=$2; size[pos] = $2; cum_size[pos++] = total}END{if(total>0){for (p = 0; p < pos; p++) {printf size[p]","cum_size[p]/total"\t"}}; printf "\n"}' >> "$asm_out_dir/gfastatsNxContig_${asm_name}.tsv"
  ) &
  pid2=$!

  # Job 3: Nx scaffold stats
  (
    printf "$asm_name\t" > "$asm_out_dir/gfastatsNxScaffold_${asm_name}.tsv"
    gfastats "$input_file" -j 1 -s s | sort -nrk2 | awk 'BEGIN{pos=0}{total+=$2; size[pos] = $2; cum_size[pos++] = total}END{if(total>0){for (p = 0; p < pos; p++) {printf size[p]","cum_size[p]/total"\t"}}; printf "\n"}' >> "$asm_out_dir/gfastatsNxScaffold_${asm_name}.tsv"
  ) &
  pid3=$!

  # Wait for all 3 jobs to complete before moving to next assembly
  wait $pid1 $pid2 $pid3

  echo "Completed: $asm_name"
}

# Initialize output files with headers
rm -f "$out_dir/gfastats.tsv" "$out_dir/gfastatsNxContig.tsv" "$out_dir/gfastatsNxScaffold.tsv"
printf 'Assembly\t# scaffolds\tTotal scaffold length\tAverage scaffold length\tScaffold N50\tScaffold auN\tScaffold L50\tLargest scaffold\tSmallest scaffold\t# contigs\tTotal contig length\tAverage contig length\tContig N50\tContig auN\tContig L50\tLargest contig\tSmallest contig\t# gaps in scaffolds\tTotal gap length in scaffolds\tAverage gap length in scaffolds\tGap N50 in scaffolds\tGap auN in scaffolds\tGap L50 in scaffolds\tLargest gap in scaffolds\tSmallest gap in scaffolds\tBase composition (A\tGC content\t# soft-masked bases\t# segments\tTotal segment length\tAverage segment length\t# gaps\t# paths\n' > "$out_dir/gfastats.tsv"

# Find and process all genome files sequentially
# Supports .fa, .fasta, .fna with or without .gz compression
shopt -s nullglob
for genome_file in "$in_dir"/*.fa "$in_dir"/*.fasta "$in_dir"/*.fna "$in_dir"/*.fa.gz "$in_dir"/*.fasta.gz "$in_dir"/*.fna.gz; do
  # Skip .chr versions if non-.chr version exists (we want unplaced contigs/scaffolds)
  basename_file=$(basename "$genome_file")
  if [[ "$basename_file" == *".chr."* ]]; then
    non_chr_file="${genome_file/.chr./.}"
    if [[ -f "$non_chr_file" ]]; then
      echo "Skipping .chr version (non-.chr exists): $basename_file"
      continue
    fi
  fi

  process_assembly "$genome_file"
done
shopt -u nullglob

# Combine all individual results into final output files
echo "Combining results..."
cat "$out_dir"/*/gfastats_*.tsv >> "$out_dir/gfastats.tsv" 2>/dev/null || true
cat "$out_dir"/*/gfastatsNxContig_*.tsv > "$out_dir/gfastatsNxContig.tsv" 2>/dev/null || true
cat "$out_dir"/*/gfastatsNxScaffold_*.tsv > "$out_dir/gfastatsNxScaffold.tsv" 2>/dev/null || true

echo "Done! Results in: $out_dir"