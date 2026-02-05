#!/usr/bin/env python3
"""
Generate per-chromosome track vs. E1 scatter plots and genome-wide correlation summaries.

"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt


# track name -> (bedgraph suffix, x-axis label, display scale factor)
TRACKS = {
    "gc": ("gc.bedgraph", "GC content (%)", 100),
    "coding_cov": ("coding_cov.bedgraph", "Coding gene coverage", 1),
}


def make_plots():
    # Paths - relative to repo root
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

    processed_wd = os.path.join(REPO_ROOT, "2_data", "2.2_processed")
    figures_wd = os.path.join(REPO_ROOT, "3_figures", "3.1_draft", "26.02.04_hic")
    os.makedirs(figures_wd, exist_ok=True)

    # find all cis.vecs.tsv files
    vecs_files = glob.glob(os.path.join(processed_wd, "*.mcool.compartments.cis.vecs.tsv"))

    for vecs_path in vecs_files:
        basename = os.path.basename(vecs_path)
        # e.g. H9_T2T_v0.1_hap1.bwa.100000.mcool.compartments.cis.vecs.tsv
        pre, post = basename.split(".bwa.", 1)
        hap = pre
        res = post.split(".")[0]

        # load eigenvector data once per resolution
        df_eigs = pd.read_csv(vecs_path, sep="\t")

        # prepare output directory
        out_dir = os.path.join(figures_wd, f"{hap}_{res}")
        os.makedirs(out_dir, exist_ok=True)

        for track_name, (suffix, xlabel, scale) in TRACKS.items():
            track_path = os.path.join(processed_wd, f"{hap}.bwa.{res}.{suffix}")
            if not os.path.exists(track_path):
                print(f"[!] Missing {track_name} track: {track_path}, skipping {hap} @ {res}")
                continue

            df_track = pd.read_csv(track_path, sep="\t")
            # normalise the value column (4th) to the track name
            value_col = df_track.columns[3]
            if value_col != track_name:
                df_track.rename(columns={value_col: track_name}, inplace=True)

            # merge on chrom/start/end
            df = pd.merge(
                df_eigs[["chrom", "start", "end", "E1"]],
                df_track[["chrom", "start", "end", track_name]],
                on=["chrom", "start", "end"],
                how="inner",
            ).dropna(subset=["E1", track_name])

            # per-chrom correlations and scatter plots
            corrs = {}
            for chrom, grp in df.groupby("chrom"):
                r = grp["E1"].corr(grp[track_name])
                corrs[chrom] = r

                fig, ax = plt.subplots(figsize=(3, 3))
                ax.scatter(grp[track_name] * scale, grp["E1"], s=5, alpha=0.5)
                ax.axhline(0, color="gray", linewidth=1)
                ax.set_xlabel(xlabel)
                ax.set_ylabel("E1")
                fig.tight_layout()
                fig.savefig(os.path.join(out_dir, f"{chrom}_{track_name}.png"), dpi=600)
                plt.close(fig)

            # summary barplot of all correlations
            chs = sorted(corrs)
            vals = [corrs[c] for c in chs]
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.bar(chs, vals)
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_xticks(range(len(chs)))
            ax.set_xticklabels(chs, rotation=90, fontsize=6)
            ax.set_ylabel(f"Pearson r (E1 vs {track_name})")
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"{hap}_{res}_{track_name}_correlation_summary.png"), dpi=600)
            plt.close(fig)


if __name__ == "__main__":
    make_plots()
