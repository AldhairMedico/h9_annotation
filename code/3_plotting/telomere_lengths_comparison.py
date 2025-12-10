#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
telomere_lengths_comparison.py

Using Teloscope terminal telomere annotations from:
  data/2_processed/asms_x1_TTAGGG_v1.3.terminal_telomeres.bed

This script:
  1) Fig 1: Paired scatter H9 hap1 vs hap2 telomere lengths (per chr arm),
     with chromosome labels and Pearson R² + p-value (bottom-right).
  2) Fig 2A: Rainclouds (half-violin + thin box + jitter) of telomere length per assembly,
     using the given palette; outliers shown hollow and excluded from stats; outputs
     pairwise tests (U test) with BH & Bonferroni p-adj.
     Fig 2B: Pairwise matrix — upper triangle shows Δ(mean Kbp) as text; lower triangle
     shows -log10(p-adj Bonf.) heatmap.
  3) Fig 3: Telomere length (Kbp) vs canonical proportion (%) scatter
     (square, y-axis starts at min(canonical)-1, legend inside).
  4) Fig 4: H9-only (hap1 & hap2) scatter of telomere length vs canonical proportion
     with chromosome name labels highlighted.

All PNG saved at 600 dpi; PNG + PDF go to ./figures.
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
from matplotlib.lines import Line2D

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

WD = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root
INPUT_BED = os.path.join(WD, "data", "2_processed", "asms_x1_TTAGGG_v1.3.terminal_telomeres.bed")
OUT_DIR = os.path.join(WD, "figures")

# Figure basenames (v5)
FIG1_BASENAME  = "scatter_H9_hap1_vs_hap2_v5"
FIG2A_BASENAME = "raincloud_by_assembly_v5"
FIG2B_BASENAME = "pairwise_heatmap_v5"
FIG3A_BASENAME = "scatter_by_assembly_v5"
FIG4_BASENAME  = "scatter_H9_length_vs_canonical_v5"

PAIRWISE_PVALS_TSV = "telomere_length_pairwise_pvalues.tsv"  # written into OUT_DIR

# Explicit mapping of raw sequence file path -> concise assembly label
ASSEMBLY_LABELS: Dict[str, str] = {
    "GCA_000001405.29_GRCh38.p14_genomic.chr.fna": "GRCh38",
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
    "GRCh38",
    "CHM13",
    "YAO mat",
    "YAO pat",
    "HG002 mat",
    "HG002 pat",
    "RPE1 hap1",
    "RPE1 hap2",
    "I002C hap1",
    "I002C hap2",
    "H9 hap1",
    "H9 hap2",
]

# Palette (user-updated; do not change)
PALETTE: Dict[str, str] = {
    "GRCh38"    : "#999999",  # Neutral Gray
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

def draw_half_violin(ax, center_x, vals, color, width_violin=0.32):
    if gaussian_kde is None:
        return
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 2 or np.max(vals) == np.min(vals):
        return
    kde = gaussian_kde(vals, bw_method=0.3)
    y = np.linspace(vals.min(), vals.max(), 300)
    dens = kde(y)
    if dens.max() > 0:
        scale = width_violin / dens.max()
        x_right = center_x + dens * scale
        ax.fill_betweenx(y, center_x, x_right, alpha=0.6, linewidth=0, color=color, zorder=1)
        ax.plot(x_right, y, linewidth=1.0, color=color, zorder=2)

def draw_thin_box(ax, center_x, vals, color, box_width=0.20):
    if len(vals) == 0:
        return
    bp = ax.boxplot(
        [vals],
        positions=[center_x],
        widths=box_width,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        zorder=3,
    )
    for patch in bp['boxes']:
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_linewidth(1.0)
    for key in ('whiskers', 'caps', 'medians'):
        for obj in bp[key]:
            obj.set_color('black')
            obj.set_linewidth(1.0)

# --------------------------------------------------------------------------
# FIG 1: H9 hap1 vs hap2 paired scatter with chr labels + R², p
# --------------------------------------------------------------------------

def make_plot_h9_hap1_vs_hap2(df: pd.DataFrame, out_dir: str, basename: str) -> None:
    sub = df[df["assembly_label"].isin(["H9 hap1", "H9 hap2"])].copy()
    if sub.empty:
        print("[WARN] No H9; skipping Fig1.")
        return

    sub["hap"] = np.where(sub["assembly_label"].eq("H9 hap1"), "hap1", "hap2")

    def get_chrom(h: str) -> str:
        for tag in ["_hap1", "_hap2"]:
            if h.endswith(tag):
                return h[:-len(tag)]
        return h
    sub["chrom"] = sub["header"].astype(str).map(get_chrom)
    sub["arm"] = sub["label"].astype(str).str.replace(r"[^pq]", "", regex=True).str[:1]
    sub["tel_id"] = sub["chrom"] + "_" + sub["arm"]

    pivot = (
        sub.pivot_table(index=["chrom", "arm", "tel_id"], columns="hap", values="tel_length_kbp", aggfunc="mean")
        .reset_index()
        .dropna(subset=["hap1", "hap2"])
    )
    if pivot.empty:
        print("[WARN] No matched H9 telomeres; skipping Fig1.")
        return

    # Pearson R^2 and p
    r2, pval = np.nan, np.nan
    if pearsonr is not None:
        r, p = pearsonr(pivot["hap1"].values, pivot["hap2"].values)
        r2, pval = r * r, p

    fig, ax = plt.subplots(figsize=(4.0, 4.0), dpi=600)

    # color by arm; label each point with chromosome number only (no 'chr', no arm)
    arm_colors = {"p": "#1f77b4", "q": "#ff7f0e"}
    for arm_val, grp in pivot.groupby("arm"):
        c = arm_colors.get(arm_val, "#7f7f7f")
        ax.scatter(grp["hap1"], grp["hap2"], s=30, alpha=0.66, edgecolor="black", linewidth=0.3, label=f"{arm_val}", c=c)
        for _, row in grp.iterrows():
            num = str(row["chrom"])
            if num.startswith("chr"):
                num = num[3:]
            ax.text(row["hap1"] + 0.1, row["hap2"] + 0.1, num, fontsize=5,  # moved farther, smaller font
                    va="bottom", ha="left", color=c)

    max_len = float(max(pivot["hap1"].max(), pivot["hap2"].max()))
    ax.plot([0, max_len], [0, max_len], linestyle="--", color="grey", linewidth=1.0)

    ax.set_xlabel("Hap1 telomere length (Kbp)")
    ax.set_ylabel("Hap2 telomere length (Kbp)")
    ax.set_xlim(0, max_len * 1.05)
    ax.set_ylim(0, max_len * 1.05)

    # === SCALE TOGGLE === set both to 'linear' or 'log' as needed
    # ax.set_xscale('linear'); ax.set_yscale('linear')
    # ax.set_xscale('log');    ax.set_yscale('log')

    # Legend: same font for title and entries
    ax.legend(title="Arm", frameon=True, loc="upper left", fontsize=8, title_fontsize=8)

    # R^2 and p (bottom-right inside)
    if np.isfinite(r2) and np.isfinite(pval):
        txt = f"R² = {r2:.3f}\np = {pval:.3g}"
        ax.text(0.98, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, bbox=dict(facecolor="white", edgecolor="none", alpha=0.6))

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", which="both", length=3, width=1.0, direction="out")

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    plt.close(fig)

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

def compute_pairwise_pvalues(df: pd.DataFrame, assemblies: List[str], out_dir: str, outfile: str) -> pd.DataFrame:
    rows = []
    inliers_by_asm = {}
    for a in assemblies:
        vals = df.loc[df["assembly_label"] == a, "tel_length_kbp"].dropna().values
        inliers, _ = split_outliers_iqr(vals)
        inliers_by_asm[a] = inliers

    for a1, a2 in combinations(assemblies, 2):
        v1 = inliers_by_asm[a1]
        v2 = inliers_by_asm[a2]
        if len(v1) == 0 or len(v2) == 0:
            pval = np.nan; method = "NA (insufficient data)"
        else:
            if mannwhitneyu is None:
                pval = np.nan; method = "Mann-Whitney U not available"
            else:
                try:
                    _stat, pval = mannwhitneyu(v1, v2, alternative="two-sided")
                    method = "Mann-Whitney U (two-sided), outliers excluded (IQR)"
                except Exception as e:
                    pval = np.nan; method = f"Error: {e}"
        rows.append({"assembly1": a1, "assembly2": a2, "p_value": pval, "method": method})

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
# FIG 2A: Rainclouds per assembly + Bonferroni bars limited to H9 pairs
# --------------------------------------------------------------------------

def make_plot_assemblies_raincloud(df: pd.DataFrame, out_dir: str, basename: str) -> List[str]:
    candidate = [a for a in ASSEMBLY_ORDER if a in set(df["assembly_label"].unique())]
    if not candidate:
        print("[WARN] No assemblies; skipping Fig2A.")
        return []

    rng = np.random.default_rng(6)
    jitter_sd = 0.04
    x_offset = 0.25
    violin_w = 0.32
    point_size = 14

    # Taller figure to accommodate brackets
    fig, ax = plt.subplots(figsize=(max(6.5, 0.5 * len(candidate) + 2), 6.8), dpi=600)

    xs, xticklabels, plotted = [], [], []
    for i, asm in enumerate(candidate, start=1):
        sub = df[df["assembly_label"] == asm]
        vals_all = sub["tel_length_kbp"].dropna().values
        inliers, outliers = split_outliers_iqr(vals_all)

        if inliers.size == 0 and outliers.size == 0:
            continue

        color = PALETTE.get(asm, "#BBBBBB")
        draw_half_violin(ax, i, inliers, color, width_violin=violin_w)

        if inliers.size > 0:
            draw_thin_box(ax, i, inliers, color)
            x_jit = rng.normal(loc=i - x_offset, scale=jitter_sd, size=inliers.size)
            ax.scatter(x_jit, inliers, s=point_size, linewidths=0.6, facecolors=color, edgecolors="black", zorder=4, alpha=0.95)

        if outliers.size > 0:
            x_jit_o = rng.normal(loc=i - x_offset, scale=jitter_sd, size=outliers.size)
            ax.scatter(x_jit_o, outliers, s=point_size, linewidths=0.8, facecolors='none', edgecolors=color, zorder=4)

        xs.append(i); xticklabels.append(asm); plotted.append(asm)

    ax.set_xticks(xs, xticklabels, rotation=45, ha="right")
    ax.set_ylabel("Telomere length (Kbp)")
    # ax.set_xlabel("Assembly")
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.5)

    # Pairwise p-vals for plotted groups only
    sig_df_all = compute_pairwise_pvalues(df, plotted, out_dir, PAIRWISE_PVALS_TSV)

    # # Only Bonferroni-significant pairs that involve H9 hap1/hap2
    # if sig_df_all is not None and not sig_df_all.empty:
    #     h9_set = {"H9 hap1", "H9 hap2"}
    #     mask = sig_df_all.apply(lambda r: (r["assembly1"] in h9_set) or (r["assembly2"] in h9_set), axis=1)
    #     sig_df_h9 = sig_df_all[mask]
    #     add_significance_bars(ax, xticklabels, sig_df_h9, alpha=0.05, use_col="p_adj_bonf", font_size=5)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", which="both", length=3, width=1.0, direction="out")

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    plt.close(fig)

    return plotted

# --------------------------------------------------------------------------
# FIG 2B: Pairwise upper Δmean (text) / lower -log10(padj) heatmap
# --------------------------------------------------------------------------

def make_plot_pairwise_heatmap(df: pd.DataFrame, assemblies: List[str], out_dir: str, basename: str) -> None:
    if not assemblies:
        return
    # Compute inlier means (consistent with stats)
    means = {}
    for a in assemblies:
        vals = df.loc[df["assembly_label"] == a, "tel_length_kbp"].dropna().values
        inl, _ = split_outliers_iqr(vals)
        means[a] = float(np.mean(inl)) if inl.size else np.nan

    # Get p-adj (Bonf) matrix
    sig_df = compute_pairwise_pvalues(df, assemblies, out_dir, PAIRWISE_PVALS_TSV)
    n = len(assemblies)
    padj = np.full((n, n), np.nan, dtype=float)
    for _, r in sig_df.iterrows():
        a1, a2 = r["assembly1"], r["assembly2"]
        if a1 in assemblies and a2 in assemblies and np.isfinite(r["p_adj_bonf"]):
            i, j = assemblies.index(a1), assemblies.index(a2)
            padj[i, j] = r["p_adj_bonf"]
            padj[j, i] = r["p_adj_bonf"]

    # Lower triangle: -log10(padj)
    with np.errstate(divide="ignore"):
        neglog = -np.log10(padj)
    # mask upper triangle (including diagonal) for heatmap
    mask_lower = np.triu(np.ones_like(neglog, dtype=bool))
    neglog_masked = np.ma.masked_where(mask_lower, neglog)

    # Build figure
    fig, ax = plt.subplots(figsize=(max(6.0, 0.5*n + 2), max(6.0, 0.5*n + 2)), dpi=600)
    im = ax.imshow(neglog_masked, cmap="viridis", interpolation="nearest")

    # Upper triangle text: Δmean(row - col) in Kbp
    for i in range(n):
        for j in range(n):
            if j > i:
                m = means.get(assemblies[i], np.nan) - means.get(assemblies[j], np.nan)
                if np.isfinite(m):
                    ax.text(j, i, f"{m:.1f}", ha="center", va="center", fontsize=7, color="black")

    # Ticks/labels
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(assemblies, rotation=45, ha="right")
    ax.set_yticklabels(assemblies)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("-log10(p-adj Bonf.)")

    # Grid lines
    ax.set_xlim(-0.5, n-0.5); ax.set_ylim(n-0.5, -0.5)
    ax.grid(False)

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    plt.close(fig)

# --------------------------------------------------------------------------
# FIG 3A: scatter (square, y>=min-1)   |   FIG 3B/4: single robust contours
# --------------------------------------------------------------------------

def legend_inside_br(ax, labels_colors):
    handles = [Line2D([0],[0], marker='o', linestyle='None', color=c, markerfacecolor=c,
                      markeredgecolor='black', markeredgewidth=0.4, label=lab)
               for lab, c in labels_colors]
    ax.legend(handles=handles, frameon=True, fontsize=7, loc="lower right")

def make_plot_length_vs_canonical_scatter(df: pd.DataFrame, out_dir: str, basename: str) -> None:
    present = []
    for a in ASSEMBLY_ORDER:
        if a in set(df["assembly_label"].unique()):
            present.append(a)
    for a in df["assembly_label"].unique():
        if a not in present:
            present.append(a)
    if not present:
        return

    fig, ax = plt.subplots(figsize=(5.0, 5.0), dpi=600)  # square
    labels_colors = []
    for asm in present:
        sub = df[df["assembly_label"] == asm]
        if sub.empty:
            continue
        col = PALETTE.get(asm, "#BBBBBB")
        labels_colors.append((asm, col))
        ax.scatter(sub["tel_length_kbp"], sub["canonical_pct"], s=16, alpha=0.85,
                   edgecolor="black", linewidth=0.2, c=col)

    ax.set_xlabel("Telomere length (Kbp)")
    ax.set_ylabel("Canonical proportion (%)")
    y_min_scatter = float(np.nanmin(df["canonical_pct"])) - 1.0
    ax.set_ylim(y_min_scatter, 100)

    # === SCALE TOGGLE === set both to 'linear' or 'log' as needed
    # ax.set_xscale('linear'); ax.set_yscale('linear')
    # ax.set_xscale('log');    ax.set_yscale('log')

    ax.grid(True, which="major", axis="both", linewidth=0.6, alpha=0.35, color="lightgrey")
    legend_inside_br(ax, labels_colors)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", which="both", length=3, width=1.0, direction="out")

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    plt.close(fig)

# --------------------------------------------------------------------------
# FIG 4: H9-only scatter with chromosome labels (length vs canonical %)
# --------------------------------------------------------------------------

def make_plot_h9_length_vs_canonical_scatter(df: pd.DataFrame, out_dir: str, basename: str) -> None:
    """
    Scatter plot of telomere length (Kbp) vs canonical proportion (%)
    for H9 hap1 and hap2 only, with chromosome name labels on each point.
    """
    sub = df[df["assembly_label"].isin(["H9 hap1", "H9 hap2"])].copy()
    if sub.empty:
        print("[WARN] No H9 data; skipping H9 length vs canonical scatter.")
        return

    # Extract chromosome name from header (e.g., "chr1_hap1" -> "chr1")
    def get_chrom(h: str) -> str:
        for tag in ["_hap1", "_hap2"]:
            if h.endswith(tag):
                return h[:-len(tag)]
        return h
    sub["chrom"] = sub["header"].astype(str).map(get_chrom)
    # Extract arm (p or q)
    sub["arm"] = sub["label"].astype(str).str.replace(r"[^pq]", "", regex=True).str[:1]
    # Create label like "1p" or "Xq"
    sub["chr_label"] = sub["chrom"].str.replace("chr", "") + sub["arm"]

    fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=600)

    labels_colors = []
    for asm in ["H9 hap1", "H9 hap2"]:
        asm_sub = sub[sub["assembly_label"] == asm]
        if asm_sub.empty:
            continue
        col = PALETTE.get(asm, "#BBBBBB")
        labels_colors.append((asm, col))
        ax.scatter(asm_sub["tel_length_kbp"], asm_sub["canonical_pct"], s=25, alpha=0.8,
                   edgecolor="black", linewidth=0.4, c=col, zorder=3)
        # Add chromosome labels
        for _, row in asm_sub.iterrows():
            ax.text(row["tel_length_kbp"] + 0.15, row["canonical_pct"] + 0.15,
                    row["chr_label"], fontsize=5, va="bottom", ha="left", color=col, zorder=4)

    ax.set_xlabel("Telomere length (Kbp)")
    ax.set_ylabel("Canonical proportion (%)")
    ax.set_xlim(0, 25)
    y_min_scatter = min(102, float(np.nanmin(sub["canonical_pct"])) - 1.0)
    ax.set_ylim(y_min_scatter, 102)
    ax.grid(True, which="major", axis="both", linewidth=0.6, alpha=0.35, color="lightgrey")

    legend_inside_br(ax, labels_colors)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", which="both", length=3, width=1.0, direction="out")

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    plt.close(fig)

# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------

def main() -> None:
    ensure_dir(OUT_DIR)
    df = parse_teloscope_terminal_file(INPUT_BED)

    # Fig 1: H9 hap1 vs hap2, paired by chr arm, labels + R²/p
    make_plot_h9_hap1_vs_hap2(df, OUT_DIR, FIG1_BASENAME)

    # Fig 2A: rainclouds per assembly + Bonferroni bars limited to H9 pairs
    plotted = make_plot_assemblies_raincloud(df, OUT_DIR, FIG2A_BASENAME)
    # Fig 2B: pairwise matrix (using the same plotted assemblies)
    make_plot_pairwise_heatmap(df, plotted, OUT_DIR, FIG2B_BASENAME)

    # Fig 3: scatter (square), y>=min-1
    make_plot_length_vs_canonical_scatter(df, OUT_DIR, FIG3A_BASENAME)

    # Fig 4: H9-only scatter with chr labels (length vs canonical %)
    make_plot_h9_length_vs_canonical_scatter(df, OUT_DIR, FIG4_BASENAME)

if __name__ == "__main__":
    main()
