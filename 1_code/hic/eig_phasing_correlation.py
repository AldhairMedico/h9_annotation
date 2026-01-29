#!/usr/bin/env python3
"""
Generate per-chromosome gc cov vs. E1 scatter plots and genome-wide correlation summaries.

Run from the directory containing your *.cis.vecs.tsv and *.gc.bedgraph files.  
Outputs will be written on figures_wd

processed_wd: /mnt/d/research/h9_annotation/2_data/2.2_processed
figures_wd:  /mnt/d/research/h9_annotation/3_figures/26.01.27_hic

(base) aldhair@jackalienware:/mnt/d/research/h9_annotation/2_data/2.2_processed$ pwd
/mnt/d/research/h9_annotation/2_data/2.2_processed
(base) aldhair@jackalienware:/mnt/d/research/h9_annotation/2_data/2.2_processed$ ls -lah
drwxrwxrwx 1 aldhair aldhair  512 Jan 23 11:52 .
drwxrwxrwx 1 aldhair aldhair  512 Dec 29 10:05 ..
-rwxrwxrwx 1 aldhair aldhair  51K Jan 23 11:43 H9_T2T_v0.1_dip.bwa.pairs.parse.stats.txt
-rwxrwxrwx 1 aldhair aldhair  51K Jan 23 11:43 H9_T2T_v0.1_dip.bwa.pairs.sorted.dedup.stats.txt
-rwxrwxrwx 1 aldhair aldhair  573 Jan 23 11:43 H9_T2T_v0.1_dip.chrom.sizes
-rwxrwxrwx 1 aldhair aldhair  40M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.10000.pairs.sorted.dedup.cool.expected_cis.tsv
-rwxrwxrwx 1 aldhair aldhair  24M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.10000.pairs.sorted.dedup.cool.insulation.tsv
-rwxrwxrwx 1 aldhair aldhair 2.9M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.10000.pairs.sorted.dedup.cool.insulation.tsv.200000.bw
-rwxrwxrwx 1 aldhair aldhair 1.2M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.100000.bins.tsv
-rwxrwxrwx 1 aldhair aldhair 1.4M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.100000.gc.bedgraph
-rwxrwxrwx 1 aldhair aldhair 463K Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.bw
-rwxrwxrwx 1 aldhair aldhair 2.1K Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.lam.txt
-rwxrwxrwx 1 aldhair aldhair 3.0M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
-rwxrwxrwx 1 aldhair aldhair 5.5M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.20000.bins.tsv
-rwxrwxrwx 1 aldhair aldhair 6.6M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.20000.gc.bedgraph
-rwxrwxrwx 1 aldhair aldhair 1.4M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.20000.pairs.sorted.dedup.cool.compartments.cis.bw
-rwxrwxrwx 1 aldhair aldhair 2.1K Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.20000.pairs.sorted.dedup.cool.compartments.cis.lam.txt
-rwxrwxrwx 1 aldhair aldhair  15M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.20000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
-rwxrwxrwx 1 aldhair aldhair 2.3M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.50000.bins.tsv
-rwxrwxrwx 1 aldhair aldhair 2.7M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.50000.gc.bedgraph
-rwxrwxrwx 1 aldhair aldhair 902K Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.50000.pairs.sorted.dedup.cool.compartments.cis.bw
-rwxrwxrwx 1 aldhair aldhair 2.1K Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.50000.pairs.sorted.dedup.cool.compartments.cis.lam.txt
-rwxrwxrwx 1 aldhair aldhair 6.0M Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.50000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
-rwxrwxrwx 1 aldhair aldhair  23K Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.pairs.parse.stats.txt
-rwxrwxrwx 1 aldhair aldhair  23K Jan 23 11:43 H9_T2T_v0.1_hap1.bwa.pairs.sorted.dedup.stats.txt
-rwxrwxrwx 1 aldhair aldhair  465 Jan 23 11:43 H9_T2T_v0.1_hap1.chrom.sizes
-rwxrwxrwx 1 aldhair aldhair  39M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.10000.pairs.sorted.dedup.cool.expected_cis.tsv
-rwxrwxrwx 1 aldhair aldhair  24M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.10000.pairs.sorted.dedup.cool.insulation.tsv
-rwxrwxrwx 1 aldhair aldhair 2.8M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.10000.pairs.sorted.dedup.cool.insulation.tsv.200000.bw
-rwxrwxrwx 1 aldhair aldhair 1.2M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.100000.bins.tsv
-rwxrwxrwx 1 aldhair aldhair 1.4M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.100000.gc.bedgraph
-rwxrwxrwx 1 aldhair aldhair 461K Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.bw
-rwxrwxrwx 1 aldhair aldhair 2.1K Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.lam.txt
-rwxrwxrwx 1 aldhair aldhair 3.0M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
-rwxrwxrwx 1 aldhair aldhair 5.5M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.20000.bins.tsv
-rwxrwxrwx 1 aldhair aldhair 6.6M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.20000.gc.bedgraph
-rwxrwxrwx 1 aldhair aldhair 1.4M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.20000.pairs.sorted.dedup.cool.compartments.cis.bw
-rwxrwxrwx 1 aldhair aldhair 2.1K Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.20000.pairs.sorted.dedup.cool.compartments.cis.lam.txt
-rwxrwxrwx 1 aldhair aldhair  15M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.20000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
-rwxrwxrwx 1 aldhair aldhair 2.3M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.50000.bins.tsv
-rwxrwxrwx 1 aldhair aldhair 2.7M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.50000.gc.bedgraph
-rwxrwxrwx 1 aldhair aldhair 898K Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.50000.pairs.sorted.dedup.cool.compartments.cis.bw
-rwxrwxrwx 1 aldhair aldhair 2.1K Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.50000.pairs.sorted.dedup.cool.compartments.cis.lam.txt
-rwxrwxrwx 1 aldhair aldhair 6.0M Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.50000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
-rwxrwxrwx 1 aldhair aldhair  23K Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.pairs.parse.stats.txt
-rwxrwxrwx 1 aldhair aldhair  23K Jan 23 11:43 H9_T2T_v0.1_hap2.bwa.pairs.sorted.dedup.stats.txt
-rwxrwxrwx 1 aldhair aldhair  465 Jan 23 11:43 H9_T2T_v0.1_hap2.chrom.sizes

"""

# import os
# import glob
# import pandas as pd
# import matplotlib.pyplot as plt

# def make_plots():
#     processed_wd = "/mnt/d/research/h9_annotation/2_data/2.2_processed"
#     figures_wd = "/mnt/d/research/h9_annotation/3_figures/26.01.27_hic"
#     os.makedirs(figures_wd, exist_ok=True)

#     # find all cis.vecs.tsv files
#     vecs_files = glob.glob(os.path.join(processed_wd, "*.cis.vecs.tsv"))

#     for vecs_path in vecs_files:
#         basename = os.path.basename(vecs_path)
#         # e.g. H9_T2T_v0.1_hap1.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
#         pre, post = basename.split(".bwa.", 1)
#         hap = pre              # H9_T2T_v0.1_hap1
#         res = post.split(".")[0]  # 100000

#         # corresponding gc bedgraph
#         gc_path = os.path.join(processed_wd, f"{hap}.bwa.{res}.gc.bedgraph")
#         if not os.path.exists(gc_path):
#             print(f"[!] Missing GC track: {gc_path}, skipping {hap} @ {res}")
#             continue

#         # prepare output directory
#         out_dir = os.path.join(figures_wd, f"{hap}_{res}")
#         os.makedirs(out_dir, exist_ok=True)

#         # load data
#         df_eigs = pd.read_csv(vecs_path, sep="\t")
#         df_gc = pd.read_csv(
#             gc_path, sep="\t", header=None,
#             names=["chrom", "start", "end", "gc"],
#             dtype={"chrom": str, "start": int, "end": int, "gc": float},
#         )

#         # merge on chrom/start/end
#         df = pd.merge(
#             df_eigs[["chrom", "start", "end", "E1"]],
#             df_gc[["chrom", "start", "end", "gc"]],
#             on=["chrom", "start", "end"],
#             how="inner",
#         ).dropna(subset=["E1", "gc"])

#         # compute per-chrom correlations and scatter plots
#         corrs = {}
#         for chrom, grp in df.groupby("chrom"):
#             r = grp["E1"].corr(grp["gc"])
#             corrs[chrom] = r

#             fig, ax = plt.subplots(figsize=(4, 3))
#             ax.scatter(grp["gc"] * 100, grp["E1"], s=5, alpha=0.5)
#             ax.axhline(0, color="gray", linewidth=1)
#             ax.set_title(f"{hap} · {chrom} (r={r:.2f})")
#             ax.set_xlabel("GC content (%)")
#             ax.set_ylabel("E1")
#             fig.tight_layout()
#             fig.savefig(os.path.join(out_dir, f"{chrom}.png"), dpi=600)
#             plt.close(fig)

#         # summary barplot of all correlations
#         chs = sorted(corrs)
#         vals = [corrs[c] for c in chs]
#         fig, ax = plt.subplots(figsize=(8, 4))
#         ax.bar(chs, vals)
#         ax.axhline(0, color="gray", linewidth=1)
#         ax.set_xticks(range(len(chs)))
#         ax.set_xticklabels(chs, rotation=90, fontsize=6)
#         ax.set_ylabel("Pearson r (E1 vs GC)")
#         ax.set_title(f"{hap} ({res} bp) genome-wide E1/GC correlation")
#         fig.tight_layout()
#         fig.savefig(os.path.join(out_dir, f"{hap}_{res}_correlation_summary.png"), dpi=600)
#         plt.close(fig)

# if __name__ == "__main__":
#     make_plots()
