#!/usr/bin/env python3
"""
Generate per-chromosome gc cov vs. E1 scatter plots and genome-wide correlation summaries.

Run from the directory containing your *.cis.vecs.tsv and *.gc.bedgraph files.  
Outputs will be written on figures_wd

"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def make_plots():
    # Paths - relative to repo root
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

    processed_wd = os.path.join(REPO_ROOT, "2_data", "2.2_processed")
    figures_wd = os.path.join(REPO_ROOT, "3_figures", "3.1_draft", "26.01.30_hic")
    os.makedirs(figures_wd, exist_ok=True)

    # find all cis.vecs.tsv files
    vecs_files = glob.glob(os.path.join(processed_wd, "*.cis.vecs.tsv"))

    for vecs_path in vecs_files:
        basename = os.path.basename(vecs_path)
        # e.g. H9_T2T_v0.1_hap1.bwa.100000.pairs.sorted.dedup.cool.compartments.cis.vecs.tsv
        pre, post = basename.split(".bwa.", 1)
        hap = pre              # H9_T2T_v0.1_hap1
        res = post.split(".")[0]  # 100000

        # corresponding gc bedgraph
        gc_path = os.path.join(processed_wd, f"{hap}.bwa.{res}.gc.bedgraph")
        if not os.path.exists(gc_path):
            print(f"[!] Missing GC track: {gc_path}, skipping {hap} @ {res}")
            continue

        # prepare output directory
        out_dir = os.path.join(figures_wd, f"{hap}_{res}")
        os.makedirs(out_dir, exist_ok=True)

        # load data
        df_eigs = pd.read_csv(vecs_path, sep="\t")
        df_gc = pd.read_csv(gc_path, sep="\t")
        df_gc.rename(columns={"GC": "gc"}, inplace=True)

        # merge on chrom/start/end
        df = pd.merge(
            df_eigs[["chrom", "start", "end", "E1"]],
            df_gc[["chrom", "start", "end", "gc"]],
            on=["chrom", "start", "end"],
            how="inner",
        ).dropna(subset=["E1", "gc"])

        # compute per-chrom correlations and scatter plots
        corrs = {}
        for chrom, grp in df.groupby("chrom"):
            r = grp["E1"].corr(grp["gc"])
            corrs[chrom] = r

            fig, ax = plt.subplots(figsize=(4, 3))
            ax.scatter(grp["gc"] * 100, grp["E1"], s=5, alpha=0.5)
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_title(f"{hap} · {chrom} (r={r:.2f})")
            ax.set_xlabel("GC content (%)")
            ax.set_ylabel("E1")
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"{chrom}.png"), dpi=600)
            plt.close(fig)

        # summary barplot of all correlations
        chs = sorted(corrs)
        vals = [corrs[c] for c in chs]
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(chs, vals)
        ax.axhline(0, color="gray", linewidth=1)
        ax.set_xticks(range(len(chs)))
        ax.set_xticklabels(chs, rotation=90, fontsize=6)
        ax.set_ylabel("Pearson r (E1 vs GC)")
        ax.set_title(f"{hap} ({res} bp) genome-wide E1/GC correlation")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"{hap}_{res}_correlation_summary.png"), dpi=600)
        plt.close(fig)

if __name__ == "__main__":
    make_plots()
