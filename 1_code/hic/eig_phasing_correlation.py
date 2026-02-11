#!/usr/bin/env python3
"""
Generate per-chromosome track vs. E1 scatter plots and genome-wide correlation summaries.

"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr


# track name -> (bedgraph suffix, x-axis label, display scale factor)
TRACKS = {
    "gc": ("gc.bedgraph", "GC content (%)", 100),
    "coding_cov": ("coding_cov.bedgraph", "Coding gene coverage", 1),
}


def compute_delta(e1_values, coding_values):
    """Delta = median(coding | top 20% E1) - median(coding | bottom 20% E1)."""
    e1 = np.asarray(e1_values, dtype=float)
    cod = np.asarray(coding_values, dtype=float)
    if len(e1) < 10:
        return np.nan
    q80, q20 = np.percentile(e1, 80), np.percentile(e1, 20)
    top, bot = e1 >= q80, e1 <= q20
    if top.sum() < 2 or bot.sum() < 2:
        return np.nan
    return float(np.median(cod[top]) - np.median(cod[bot]))


def make_plots():
    # Paths - relative to repo root
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

    processed_wd = os.path.join(REPO_ROOT, "2_data", "2.2_processed")
    figures_wd = os.path.join(REPO_ROOT, "3_figures", "3.1_draft", "26.02.08_hic")
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
            # normalise the value column (5th, index 4) to the track name
            value_col = df_track.columns[4]
            if value_col != track_name:
                df_track.rename(columns={value_col: track_name}, inplace=True)

            # merge on chrom/start/end
            df = pd.merge(
                df_eigs[["chrom", "start", "end", "E1"]],
                df_track[["chrom", "start", "end", track_name]],
                on=["chrom", "start", "end"],
                how="inner",
            ).dropna(subset=["E1", track_name])

            # First pass: compute correlations and p-values for all chromosomes
            chrom_stats = {}
            chrom_data = {}
            for chrom, grp in df.groupby("chrom"):
                chrom_data[chrom] = grp
                pearson_r, pearson_p = pearsonr(grp[track_name], grp["E1"])
                spearman_rho, spearman_p = spearmanr(grp[track_name], grp["E1"])
                chrom_stats[chrom] = {
                    "pearson_r": pearson_r,
                    "pearson_p": pearson_p,
                    "spearman_rho": spearman_rho,
                    "spearman_p": spearman_p,
                }

                if track_name == "coding_cov":
                    delta = compute_delta(grp["E1"].values, grp[track_name].values)
                    flip_recommended = (
                        True if (not np.isnan(delta) and delta < 0)
                        else (False if not np.isnan(delta) else np.nan)
                    )
                    chrom_stats[chrom]["delta"] = delta
                    chrom_stats[chrom]["flip_recommended"] = flip_recommended

            chroms_sorted = sorted(chrom_stats.keys())

            # Generate scatter plots with legends
            for chrom in chroms_sorted:
                grp = chrom_data[chrom]
                stats = chrom_stats[chrom]

                fig, ax = plt.subplots(figsize=(4, 4))
                ax.scatter(grp[track_name] * scale, grp["E1"], s=5, alpha=0.5)
                ax.axhline(0, color="gray", linewidth=1)
                ax.set_xlabel(xlabel)
                ax.set_ylabel("E1")

                # Add legend with correlation statistics
                stats_text = (
                    f"ρ = {stats['spearman_rho']:.3f}, p = {stats['spearman_p']:.2e}\n"
                    f"r = {stats['pearson_r']:.3f}, p = {stats['pearson_p']:.2e}"
                )
                ax.text(
                    0.95, 0.95, stats_text,
                    transform=ax.transAxes,
                    fontsize=7,
                    verticalalignment="top",
                    horizontalalignment="right",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgrey", alpha=0.9),
                )

                fig.tight_layout()
                fig.savefig(os.path.join(out_dir, f"{chrom}_{track_name}.png"), dpi=600)
                plt.close(fig)

            # Summary barplot: Spearman correlations
            spearman_vals = [chrom_stats[c]["spearman_rho"] for c in chroms_sorted]
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.bar(chroms_sorted, spearman_vals)
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_xticks(range(len(chroms_sorted)))
            ax.set_xticklabels(chroms_sorted, rotation=90, fontsize=6)
            ax.set_ylabel(f"Spearman ρ (E1 vs {track_name})")
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"{hap}_{res}_{track_name}_spearman_summary.png"), dpi=600)
            plt.close(fig)

            # Summary barplot: Pearson correlations
            pearson_vals = [chrom_stats[c]["pearson_r"] for c in chroms_sorted]
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.bar(chroms_sorted, pearson_vals)
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_xticks(range(len(chroms_sorted)))
            ax.set_xticklabels(chroms_sorted, rotation=90, fontsize=6)
            ax.set_ylabel(f"Pearson r (E1 vs {track_name})")
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"{hap}_{res}_{track_name}_pearson_summary.png"), dpi=600)
            plt.close(fig)

            # Delta barplot (coding_cov only)
            if track_name == "coding_cov":
                delta_vals = [chrom_stats[c].get("delta", 0) for c in chroms_sorted]
                flip_flags = [chrom_stats[c].get("flip_recommended", False) for c in chroms_sorted]
                colors = [
                    "#CCCCCC" if (isinstance(f, float) and np.isnan(f))
                    else "#E38AAA" if f
                    else "#6F9DD0"
                    for f in flip_flags
                ]
                fig, ax = plt.subplots(figsize=(8, 4))
                ax.bar(chroms_sorted, delta_vals, color=colors)
                ax.axhline(0, color="gray", linewidth=1)
                ax.set_xticks(range(len(chroms_sorted)))
                ax.set_xticklabels(chroms_sorted, rotation=90, fontsize=6)
                ax.set_ylabel("Δ (median coding top 20% E1 − bottom 20% E1)")
                fig.tight_layout()
                fig.savefig(os.path.join(out_dir, f"{hap}_{res}_delta_summary.png"), dpi=600)
                plt.close(fig)

            # Export summary statistics to TSV
            summary_rows = []
            for chrom in chroms_sorted:
                stats = chrom_stats[chrom]
                row = {
                    "track": track_name,
                    "chrom": chrom,
                    "spearman_rho": stats["spearman_rho"],
                    "spearman_p": stats["spearman_p"],
                    "pearson_r": stats["pearson_r"],
                    "pearson_p": stats["pearson_p"],
                }
                if track_name == "coding_cov":
                    row["delta"] = stats.get("delta", "")
                    row["flip_recommended"] = stats.get("flip_recommended", "")
                summary_rows.append(row)
            df_summary = pd.DataFrame(summary_rows)
            summary_path = os.path.join(out_dir, f"{hap}_{res}_{track_name}_summary_stats.tsv")
            df_summary.to_csv(summary_path, sep="\t", index=False)


if __name__ == "__main__":
    make_plots()
