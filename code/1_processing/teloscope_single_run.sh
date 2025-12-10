raw_data_dir="/mnt/d/research/h9_annotation/data/1_raw/"
processed_data_dir="/mnt/d/research/h9_annotation/data/2_processed/"

# Run teloscope on haploid assembly
assembly_name="GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna"
assembly_path="$raw_data_dir/$assembly_name"
assembly_report="$processed_data_dir/$assembly_name.telo.report"

touch "$assembly_report"
teloscope -f "$assembly_path" -o "$processed_data_dir" -j 2 -x 1 -gremi > "$assembly_report"

# # Run teloscope on diploid assembly
# assembly_name_1="I002Cv0.7.hap1.fasta.gz"
# assembly_path_1="$raw_data_dir/$assembly_name_1"
# assembly_report_1="$processed_data_dir/$assembly_name_1.telo.report"

# touch "$assembly_report_1"
# teloscope -f "$assembly_path_1" -o "$processed_data_dir" -j 2 -x 1 -gremi > "$assembly_report_1"

# assembly_name_2="I002Cv0.7.hap2.fasta.gz"
# assembly_path_2="$raw_data_dir/$assembly_name_2"
# assembly_report_2="$processed_data_dir/$assembly_name_2.telo.report"

# touch "$assembly_report_2"
# teloscope -f "$assembly_path_2" -o "$processed_data_dir" -j 2 -x 1 -gremi > "$assembly_report_2"