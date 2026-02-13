#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
telomere_lengths_comparison.py

Using Teloscope terminal telomere annotations.

This script:
  1) Rainclouds (half-violin + thin box + jitter) of telomere length per assembly
     in both vertical and horizontal layouts; outliers shown hollow and excluded
     from stats.
  2) Outputs pairwise tests (U test) with BH & Bonferroni p-adj to TSV.

All PNG saved at 600 dpi; PNG + PDF go to plots_dir.
"""

import os
import math
from itertools import combinations
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Global font: Arial 10
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

# stats
try:
    from scipy.stats import mannwhitneyu, pearsonr, gaussian_kde
except Exception:
    mannwhitneyu = None
    pearsonr = None
    gaussian_kde = None

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------

script_dir  = os.path.dirname(os.path.abspath(__file__))
base_dir    = os.path.dirname(os.path.dirname(script_dir))
data_dir    = os.path.join(base_dir, "2_data", "2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures", "3.1_draft")
run_label   = "v1"

plots_dir = os.path.join(figures_dir, "26.01.29_telomeres")
os.makedirs(plots_dir, exist_ok=True)

INPUT_BED = os.path.join(
    data_dir, "25.12.10_teloscope_compiled",
    "25.12.10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed",
)
OUT_DIR = plots_dir

# Figure basenames
FIG_VERTICAL_BASENAME  = "raincloud_by_assembly_vertical"
FIG_HORIZONTAL_BASENAME = "raincloud_by_assembly_horizontal"

PAIRWISE_PVALS_TSV = "telomere_length_pairwise_pvalues.tsv"  # written into OUT_DIR

# Explicit mapping of raw sequence file path -> concise assembly label
ASSEMBLY_LABELS: Dict[str, str] = {
    "GCA_009914755.4_T2T-CHM13v2.0_genomic.chr.fna": "CHM13",
    "GWHGEYB00000000.1.genome.fasta.gz": "YAO pat",
    "GWHGEYC00000000.1.genome.fasta.gz": "YAO mat",
    "GCA_018852605.3_hg002v1.1.pat_genomic.fna": "HG002 pat",
    "GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna": "HG002 mat",
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.chr.fna": "RPE1 hap1",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.chr.fna": "RPE1 hap2",
    "I002Cv0.7.hap1.fasta.gz": "I002C hap1",
    "I002Cv0.7.hap2.fasta.gz": "I002C hap2",
    "H9_T2T_v0.1_hap1.fasta": "H9 hap1",
    "H9_T2T_v0.1_hap2.fasta": "H9 hap2",
}

# X order
ASSEMBLY_ORDER: List[str] = [
    "CHM13",
    "HG002 mat",
    "HG002 pat",
    "I002C hap1",
    "I002C hap2",
    "YAO mat",
    "YAO pat",
    "RPE1 hap1",
    "RPE1 hap2",
    "H9 hap1",
    "H9 hap2",
]

# Palette (user-updated; do not change)
PALETTE: Dict[str, str] = {
    "CHM13"     : "#F0E442",  # Yellow
    "YAO pat" : "#FF5353",
    "YAO mat" : "#FFA5A5",
    "HG002 mat" : "#0072B2",  # Okabe–Ito Blue
    "HG002 pat" : "#56B4E9",  # Okabe–Ito Sky Blue
    "RPE1 hap1" : "#984EA3",  # Purple
    "RPE1 hap2" : "#CC79A7",  # Magenta
    "I002C hap1": "#D55E00",  # Vermillion
    "I002C hap2": "#E69F00",  # Orange
    "H9 hap1"   : "#00796B",  # Teal 700
    "H9 hap2"   : "#80CBC4",  # Teal 200
}

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

# --------------------------------------------------------------------------
# PARSE INPUT
# --------------------------------------------------------------------------

def parse_teloscope_terminal_file(path: str) -> pd.DataFrame:
    """
    Parse Teloscope terminal telomere BED-like file with 11/12 columns.

    Returns columns:
      header,start,end,length,label,fwdCounts,revCounts,canCounts,nonCanCounts,chrSize,
      assembly,assembly_label,tel_length_kbp,canonical_pct
    """
    records = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 11:
                header, s, e, length, label, fwd, rev, can_, noncan, chrsize, asm = parts
            elif len(parts) == 12 and parts[1] == "Chromosome":
                header, _chrom, s, e, length, label, fwd, rev, can_, noncan, chrsize, asm = parts
            else:
                continue

            try:
                start = int(s); end = int(e); length_bp = int(length)
                fwd_counts = int(fwd); rev_counts = int(rev)
                can_counts = int(can_); noncan_counts = int(noncan)
                chr_size = int(chrsize)
            except ValueError:
                continue

            if length_bp <= 0:
                continue

            tel_kbp = length_bp / 1000.0
            canonical_pct = (can_counts * 600.0) / float(length_bp) if length_bp > 0 else np.nan
            if not math.isnan(canonical_pct):
                canonical_pct = max(0.0, min(100.0, canonical_pct))

            asm_label = ASSEMBLY_LABELS.get(asm, asm)

            records.append({
                "header": header,
                "start": start,
                "end": end,
                "length": length_bp,
                "label": str(label),
                "fwdCounts": fwd_counts,
                "revCounts": rev_counts,
                "canCounts": can_counts,
                "nonCanCounts": noncan_counts,
                "chrSize": chr_size,
                "assembly": asm,
                "assembly_label": asm_label,
                "tel_length_kbp": tel_kbp,
                "canonical_pct": canonical_pct,
            })
    if not records:
        raise RuntimeError(f"No valid records parsed from {path}")
    return pd.DataFrame.from_records(records)

# --------------------------------------------------------------------------
# OUTLIER HANDLING (exclude from stats; plot hollow)
# --------------------------------------------------------------------------

def split_outliers_iqr(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (inliers, outliers) using 1.5*IQR rule."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return x, np.array([])
    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    lo = q1 - 1.5 * iqr
    hi = q3 + 1.5 * iqr
    inliers = x[(x >= lo) & (x <= hi)]
    outliers = x[(x < lo) | (x > hi)]
    return inliers, outliers

# --------------------------------------------------------------------------
# RAINCLOUD HELPERS ─ half-violin + thin box + jitter (no seaborn)
# --------------------------------------------------------------------------

def draw_half_violin(ax, center_x, vals, color, width_violin=0.32, horizontal=False):
    if gaussian_kde is None:
        return
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 2 or np.max(vals) == np.min(vals):
        return
    kde = gaussian_kde(vals, bw_method=0.3)
    t = np.linspace(vals.min(), vals.max(), 300)
    dens = kde(t)
    if dens.max() > 0:
        scale = width_violin / dens.max()
        if horizontal:
            y_top = center_x + dens * scale
            ax.fill_between(t, center_x, y_top, alpha=0.6, linewidth=0, color=color, zorder=1)
            ax.plot(t, y_top, linewidth=1.0, color=color, zorder=2)
        else:
            x_right = center_x + dens * scale
            ax.fill_betweenx(t, center_x, x_right, alpha=0.6, linewidth=0, color=color, zorder=1)
            ax.plot(x_right, t, linewidth=1.0, color=color, zorder=2)

def draw_thin_box(ax, center_x, vals, color, box_width=0.12, horizontal=False):
    if len(vals) == 0:
        return
    bp = ax.boxplot(
        [vals],
        positions=[center_x],
        widths=box_width,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        vert=not horizontal,
        zorder=3,
    )
    for patch in bp['boxes']:
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_linewidth(0.6)
    for key in ('whiskers', 'caps', 'medians'):
        for obj in bp[key]:
            obj.set_color('black')
            obj.set_linewidth(0.6)

# --------------------------------------------------------------------------
# Pairwise tests (exclude outliers) + p-adj
# --------------------------------------------------------------------------

def adjust_pvals_bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR."""
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    order = np.argsort(pvals)
    ranked = np.empty(n, dtype=float)
    for i, idx in enumerate(order, start=1):
        ranked[idx] = pvals[idx] * n / i
    ranked[order[::-1]] = np.minimum.accumulate(ranked[order[::-1]])
    return np.clip(ranked, 0, 1)

def _inlier_stats(v: np.ndarray) -> Dict[str, float]:
    """Return mean, median, sd, min, max for an inlier array."""
    if v.size == 0:
        return {"mean": np.nan, "median": np.nan, "sd": np.nan, "min": np.nan, "max": np.nan}
    return {
        "mean": float(np.mean(v)),
        "median": float(np.median(v)),
        "sd": float(np.std(v, ddof=1)) if v.size > 1 else 0.0,
        "min": float(np.min(v)),
        "max": float(np.max(v)),
    }

def compute_pairwise_pvalues(df: pd.DataFrame, assemblies: List[str], out_dir: str, outfile: str) -> pd.DataFrame:
    rows = []
    inliers_by_asm = {}
    for a in assemblies:
        vals = df.loc[df["assembly_label"] == a, "tel_length_kbp"].dropna().values
        inliers, _ = split_outliers_iqr(vals)
        inliers_by_asm[a] = inliers

    method_str = "Mann-Whitney U (two-sided), outliers excluded (IQR)"
    print(f"[INFO] Statistical method: {method_str}")

    for a1, a2 in combinations(assemblies, 2):
        v1 = inliers_by_asm[a1]
        v2 = inliers_by_asm[a2]
        s1 = _inlier_stats(v1)
        s2 = _inlier_stats(v2)
        if len(v1) == 0 or len(v2) == 0:
            pval = np.nan
        else:
            if mannwhitneyu is None:
                pval = np.nan
            else:
                try:
                    _stat, pval = mannwhitneyu(v1, v2, alternative="two-sided")
                except Exception:
                    pval = np.nan
        rows.append({
            "assembly1": a1,
            "assembly2": a2,
            "mean1": s1["mean"], "median1": s1["median"], "sd1": s1["sd"],
            "min1": s1["min"], "max1": s1["max"],
            "mean2": s2["mean"], "median2": s2["median"], "sd2": s2["sd"],
            "min2": s2["min"], "max2": s2["max"],
            "p_value": pval,
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        pvals = out["p_value"].to_numpy()
        with np.errstate(invalid="ignore"):
            out["p_adj_BH"] = adjust_pvals_bh(np.where(np.isfinite(pvals), pvals, 1.0))
            m = out.shape[0] if out.shape[0] > 0 else 1
            out["p_adj_bonf"] = np.minimum(1.0, pvals * m)
        ensure_dir(out_dir)
        out_path = os.path.join(out_dir, outfile)
        out.to_csv(out_path, sep="\t", index=False)
        print(f"[SAVE] Pairwise p-values -> {out_path}")
    return out

# --------------------------------------------------------------------------
# Significance bars
# --------------------------------------------------------------------------

def add_significance_bars(ax, xlabels: List[str], sig_df: pd.DataFrame, alpha=0.05, use_col="p_adj_bonf", font_size=6):
    """
    Draw bracket-style significance bars for pairs with adjusted p < alpha.
    """
    if sig_df is None or sig_df.empty:
        return
    pos = {lab: i+1 for i, lab in enumerate(xlabels)}
    sdf = sig_df.copy()
    if use_col not in sdf.columns:
        return
    sdf = sdf[np.isfinite(sdf[use_col]) & (sdf[use_col] < alpha)]
    if sdf.empty:
        return
    sdf = sdf.sort_values(use_col)
    y_min, y_max = ax.get_ylim()
    height = y_max
    h_step = (y_max - y_min) * 0.06
    used_spans = []
    for _, row in sdf.iterrows():
        a1, a2, p_adj = row["assembly1"], row["assembly2"], row[use_col]
        if a1 not in pos or a2 not in pos:
            continue
        x1, x2 = pos[a1], pos[a2]
        if x1 == x2:
            continue
        x_lo, x_hi = sorted([x1, x2])
        cur_h = height
        for (lo, hi, h) in used_spans:
            if not (x_hi < lo or x_lo > hi):
                cur_h = max(cur_h, h + h_step)
        ax.plot([x_lo, x_lo, x_hi, x_hi], [cur_h, cur_h + h_step*0.35, cur_h + h_step*0.35, cur_h],
                color="black", linewidth=0.8)
        ax.text((x_lo + x_hi)/2.0, cur_h + h_step*0.45, f"{p_adj:.2e}", ha="center", va="bottom", fontsize=font_size)
        used_spans.append((x_lo, x_hi, cur_h + h_step*0.35))
        height = max(height, cur_h + h_step)
    ax.set_ylim(y_min, height + h_step*1.2)

# --------------------------------------------------------------------------
# RAINCLOUD: Vertical layout
# --------------------------------------------------------------------------

def make_plot_raincloud_vertical(df: pd.DataFrame, out_dir: str, basename: str) -> List[str]:
    candidate = [a for a in ASSEMBLY_ORDER if a in set(df["assembly_label"].unique())]
    if not candidate:
        print("[WARN] No assemblies; skipping vertical raincloud.")
        return []

    rng = np.random.default_rng(6)
    jitter_sd = 0.03
    x_offset = 0.20
    violin_w = 0.28
    point_size = 4

    fig, ax = plt.subplots(figsize=(2, 2), dpi=600)

    xs, xticklabels, plotted = [], [], []
    for i, asm in enumerate(candidate, start=1):
        sub = df[df["assembly_label"] == asm]
        vals_all = sub["tel_length_kbp"].dropna().values
        inliers, outliers = split_outliers_iqr(vals_all)

        if inliers.size == 0 and outliers.size == 0:
            continue

        color = PALETTE.get(asm, "#BBBBBB")
        draw_half_violin(ax, i, inliers, color, width_violin=violin_w, horizontal=False)

        if inliers.size > 0:
            draw_thin_box(ax, i, inliers, color, horizontal=False)
            x_jit = rng.normal(loc=i - x_offset, scale=jitter_sd, size=inliers.size)
            ax.scatter(x_jit, inliers, s=point_size, linewidths=0.3, facecolors=color, edgecolors="black", zorder=4, alpha=0.95)

        if outliers.size > 0:
            x_jit_o = rng.normal(loc=i - x_offset, scale=jitter_sd, size=outliers.size)
            ax.scatter(x_jit_o, outliers, s=point_size, linewidths=0.4, facecolors='none', edgecolors=color, zorder=4)

        xs.append(i); xticklabels.append(asm); plotted.append(asm)

    ax.set_xticks(xs, xticklabels, rotation=45, ha="right")
    ax.set_ylabel("Telomere length (kbp)")
    ax.grid(axis="y", linestyle=":", linewidth=0.4, alpha=0.5)

    for spine in ax.spines.values():
        spine.set_linewidth(0.4)
    ax.tick_params(axis="both", which="both", length=2, width=0.4, direction="out")

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    plt.close(fig)

    return plotted

# --------------------------------------------------------------------------
# RAINCLOUD: Horizontal layout
# --------------------------------------------------------------------------

def make_plot_raincloud_horizontal(df: pd.DataFrame, out_dir: str, basename: str) -> List[str]:
    candidate = [a for a in ASSEMBLY_ORDER if a in set(df["assembly_label"].unique())]
    if not candidate:
        print("[WARN] No assemblies; skipping horizontal raincloud.")
        return []

    rng = np.random.default_rng(6)
    jitter_sd = 0.03
    y_offset = 0.20
    violin_w = 0.28
    point_size = 4

    fig, ax = plt.subplots(figsize=(2, 2), dpi=600)

    ys, yticklabels, plotted = [], [], []
    for i, asm in enumerate(candidate, start=1):
        sub = df[df["assembly_label"] == asm]
        vals_all = sub["tel_length_kbp"].dropna().values
        inliers, outliers = split_outliers_iqr(vals_all)

        if inliers.size == 0 and outliers.size == 0:
            continue

        color = PALETTE.get(asm, "#BBBBBB")
        draw_half_violin(ax, i, inliers, color, width_violin=violin_w, horizontal=True)

        if inliers.size > 0:
            draw_thin_box(ax, i, inliers, color, horizontal=True)
            y_jit = rng.normal(loc=i - y_offset, scale=jitter_sd, size=inliers.size)
            ax.scatter(inliers, y_jit, s=point_size, linewidths=0.3, facecolors=color, edgecolors="black", zorder=4, alpha=0.95)

        if outliers.size > 0:
            y_jit_o = rng.normal(loc=i - y_offset, scale=jitter_sd, size=outliers.size)
            ax.scatter(outliers, y_jit_o, s=point_size, linewidths=0.4, facecolors='none', edgecolors=color, zorder=4)

        ys.append(i); yticklabels.append(asm); plotted.append(asm)

    ax.set_yticks(ys, yticklabels)
    ax.set_xlabel("Telomere length (kbp)")
    ax.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.5)

    for spine in ax.spines.values():
        spine.set_linewidth(0.4)
    ax.tick_params(axis="both", which="both", length=2, width=0.4, direction="out")

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    plt.close(fig)

    return plotted

# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------

def main() -> None:
    ensure_dir(OUT_DIR)
    df = parse_teloscope_terminal_file(INPUT_BED)

    # Vertical raincloud
    plotted = make_plot_raincloud_vertical(df, OUT_DIR, FIG_VERTICAL_BASENAME)

    # Horizontal raincloud
    make_plot_raincloud_horizontal(df, OUT_DIR, FIG_HORIZONTAL_BASENAME)

    # Pairwise p-values TSV only
    compute_pairwise_pvalues(df, plotted, OUT_DIR, PAIRWISE_PVALS_TSV)

if __name__ == "__main__":
    main()
