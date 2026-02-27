#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h9_tl_haps.py

Using Teloscope terminal telomere annotations from:
  data/2_processed/asms_x1_TTAGGG_v1.3.terminal_telomeres.bed

This script:
  1) Fig 1: Paired scatter H9 hap1 vs hap2 telomere lengths (per chr arm),
     with chromosome labels and Pearson R² + p-value (bottom-right).

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

# stats
try:
    from scipy.stats import mannwhitneyu, pearsonr, spearmanr, wilcoxon, shapiro, gaussian_kde, ttest_ind, ttest_rel
except Exception:
    mannwhitneyu = None
    pearsonr = None
    spearmanr = None
    wilcoxon = None
    shapiro = None
    gaussian_kde = None
    ttest_ind = None
    ttest_rel = None

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(os.path.dirname(script_dir))  # repo root
data_dir = os.path.join(base_dir, "2_data", "2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures", "3.1_draft")

INPUT_BED = os.path.join(
    data_dir, "25.12.10_teloscope_compiled",
    "25.12.10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed"
)
OUT_DIR = os.path.join(figures_dir, "26.02.25_telomere_stats")

# Figure basenames
FIG1_BASENAME  = "scatter_H9_hap1_vs_hap2_pearson"
FIG2_BASENAME  = "scatter_H9_hap1_vs_hap2_spearman"

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
# STATISTICAL ANALYSIS: normality, pairing, and haplotype differences
# --------------------------------------------------------------------------

def perform_statistical_analysis(df: pd.DataFrame, out_dir: str) -> None:
    """
    Comprehensive statistical analysis of H9 telomeres:
    1. Report n explicitly for each haplotype
    2. Identify matched vs unmatched ends
    3. Normality tests (Shapiro-Wilk)
    4. Paired analysis (Wilcoxon signed-rank test)
    """
    sub = df[df["assembly_label"].isin(["H9 hap1", "H9 hap2"])].copy()
    if sub.empty:
        print("[WARN] No H9 data for statistical analysis.")
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
    
    # Count ends per haplotype
    n_hap1 = len(sub[sub["hap"] == "hap1"])
    n_hap2 = len(sub[sub["hap"] == "hap2"])
    
    # Create pivot to identify matched pairs
    pivot = (
        sub.pivot_table(index=["chrom", "arm", "tel_id"], columns="hap", 
                       values="tel_length_kbp", aggfunc="mean")
        .reset_index()
    )
    
    # Identify matched and unmatched ends
    matched = pivot.dropna(subset=["hap1", "hap2"])
    hap1_only = pivot[pivot["hap2"].isna() & pivot["hap1"].notna()]
    hap2_only = pivot[pivot["hap1"].isna() & pivot["hap2"].notna()]
    
    n_matched_pairs = len(matched)
    
    # Prepare output file
    stats_file = os.path.join(out_dir, "H9_telomere_statistics.txt")
    ensure_dir(out_dir)
    
    with open(stats_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("H9 TELOMERE LENGTH STATISTICAL ANALYSIS\n")
        f.write("=" * 80 + "\n\n")
        
        # Sample sizes
        f.write("1. SAMPLE SIZES\n")
        f.write("-" * 80 + "\n")
        f.write(f"Hap1 total chromosome ends: {n_hap1}\n")
        f.write(f"Hap2 total chromosome ends: {n_hap2}\n")
        f.write(f"Matched pairs (chr-end present in both haplotypes): {n_matched_pairs}\n\n")
        
        if len(hap1_only) > 0:
            f.write("Hap1-only ends (missing in hap2):\n")
            for _, row in hap1_only.iterrows():
                f.write(f"  {row['tel_id']}: {row['hap1']:.2f} kbp\n")
            f.write("\n")
        
        if len(hap2_only) > 0:
            f.write("Hap2-only ends (missing in hap1):\n")
            for _, row in hap2_only.iterrows():
                f.write(f"  {row['tel_id']}: {row['hap2']:.2f} kbp\n")
            f.write("\n")
        
        # Descriptive statistics
        f.write("\n2. DESCRIPTIVE STATISTICS\n")
        f.write("-" * 80 + "\n")
        hap1_vals = sub[sub["hap"] == "hap1"]["tel_length_kbp"].values
        hap2_vals = sub[sub["hap"] == "hap2"]["tel_length_kbp"].values
        
        f.write(f"Hap1 (n={len(hap1_vals)}): "
                f"mean={np.mean(hap1_vals):.2f}, "
                f"median={np.median(hap1_vals):.2f}, "
                f"SD={np.std(hap1_vals, ddof=1):.2f}, "
                f"range=[{np.min(hap1_vals):.2f}, {np.max(hap1_vals):.2f}] kbp\n")
        
        f.write(f"Hap2 (n={len(hap2_vals)}): "
                f"mean={np.mean(hap2_vals):.2f}, "
                f"median={np.median(hap2_vals):.2f}, "
                f"SD={np.std(hap2_vals, ddof=1):.2f}, "
                f"range=[{np.min(hap2_vals):.2f}, {np.max(hap2_vals):.2f}] kbp\n")
        
        # Normality tests
        f.write("\n3. NORMALITY TESTS (Shapiro-Wilk)\n")
        f.write("-" * 80 + "\n")
        
        if shapiro is not None:
            if len(hap1_vals) >= 3:
                w1, p1 = shapiro(hap1_vals)
                f.write(f"Hap1: W={w1:.4f}, p={p1:.4g} (n={len(hap1_vals)})\n")
            else:
                f.write(f"Hap1: insufficient data (n={len(hap1_vals)})\n")
            
            if len(hap2_vals) >= 3:
                w2, p2 = shapiro(hap2_vals)
                f.write(f"Hap2: W={w2:.4f}, p={p2:.4g} (n={len(hap2_vals)})\n")
            else:
                f.write(f"Hap2: insufficient data (n={len(hap2_vals)})\n")
        else:
            f.write("scipy.stats.shapiro not available\n")
        
        # Paired analysis
        f.write("\n4. PAIRED ANALYSIS (matched chromosome ends)\n")
        f.write("-" * 80 + "\n")
        
        if n_matched_pairs > 0:
            matched_hap1 = matched["hap1"].values
            matched_hap2 = matched["hap2"].values
            differences = matched_hap1 - matched_hap2
            
            f.write(f"Matched pairs: n={n_matched_pairs}\n")
            f.write(f"Mean difference (hap1 - hap2): {np.mean(differences):.2f} kbp\n")
            f.write(f"Median difference: {np.median(differences):.2f} kbp\n")
            f.write(f"SD of differences: {np.std(differences, ddof=1):.2f} kbp\n\n")
            
            # Pearson correlation
            if pearsonr is not None:
                r, p = pearsonr(matched_hap1, matched_hap2)
                f.write(f"Pearson correlation: r={r:.4f}, R²={r*r:.4f}, p={p:.4g}\n")
            
            # Spearman correlation
            if spearmanr is not None:
                rho, p_spear = spearmanr(matched_hap1, matched_hap2)
                f.write(f"Spearman correlation: ρ={rho:.4f}, p={p_spear:.4g}\n\n")
            
            # Wilcoxon signed-rank test (paired, non-parametric)
            if wilcoxon is not None and n_matched_pairs >= 1:
                try:
                    stat, p_wilc = wilcoxon(matched_hap1, matched_hap2)
                    f.write(f"Wilcoxon signed-rank test: statistic={stat:.2f}, p={p_wilc:.4g} (n pairs={n_matched_pairs})\n")
                    f.write("Interpretation: Tests if median difference between paired observations is zero\n\n")
                except Exception as e:
                    f.write(f"Wilcoxon test failed: {e}\n\n")
            else:
                f.write("Insufficient data for Wilcoxon signed-rank test\n\n")
            
            # Paired t-test (parametric alternative to Wilcoxon)
            if ttest_rel is not None and n_matched_pairs >= 2:
                try:
                    stat_t, p_t = ttest_rel(matched_hap1, matched_hap2)
                    f.write(f"Paired t-test: t={stat_t:.4f}, p={p_t:.4g} (n pairs={n_matched_pairs})\n")
                    f.write("Interpretation: Tests if mean difference between paired observations is zero\n\n")
                except Exception as e:
                    f.write(f"Paired t-test failed: {e}\n\n")
        else:
            f.write("No matched pairs available for paired analysis.\n\n")
        
        # Independent samples analysis (treating haplotypes as separate groups)
        f.write("\n5. INDEPENDENT SAMPLES ANALYSIS (ignoring pairing)\n")
        f.write("-" * 80 + "\n")
        f.write("NOTE: This approach ignores the natural pairing of chromosome ends.\n")
        f.write("The paired analysis above is more appropriate for this dataset.\n\n")
        
        if len(hap1_vals) > 1 and len(hap2_vals) > 1:
            # Welch's t-test (does not assume equal variances or pairing)
            if ttest_ind is not None:
                try:
                    stat_welch, p_welch = ttest_ind(hap1_vals, hap2_vals, equal_var=False)
                    df_approx = len(hap1_vals) + len(hap2_vals) - 2  # approximate df
                    f.write(f"Welch's t-test: t={stat_welch:.4f}, p={p_welch:.4g}, df≈{df_approx}\n")
                    f.write(f"  (n_hap1={len(hap1_vals)}, n_hap2={len(hap2_vals)})\n")
                    f.write("Interpretation: Tests if means differ between independent groups\n\n")
                except Exception as e:
                    f.write(f"Welch's t-test failed: {e}\n\n")
            
            # Mann-Whitney U test (non-parametric independent samples)
            if mannwhitneyu is not None:
                try:
                    stat_mw, p_mw = mannwhitneyu(hap1_vals, hap2_vals, alternative='two-sided')
                    f.write(f"Mann-Whitney U test: U={stat_mw:.2f}, p={p_mw:.4g}\n")
                    f.write(f"  (n_hap1={len(hap1_vals)}, n_hap2={len(hap2_vals)})\n")
                    f.write("Interpretation: Tests if distributions differ between independent groups\n")
                except Exception as e:
                    f.write(f"Mann-Whitney U test failed: {e}\n")
        else:
            f.write("Insufficient data for independent samples tests\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    print(f"[OK] Statistical analysis saved to {stats_file}")

# --------------------------------------------------------------------------
# FIG 1: H9 hap1 vs hap2 paired scatter with chr labels + Pearson R², p
# --------------------------------------------------------------------------

def make_plot_h9_hap1_vs_hap2_pearson(df: pd.DataFrame, out_dir: str, basename: str) -> None:
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
    
    n_pairs = len(pivot)

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

    ax.set_xlabel("Hap1 telomere length (kbp)")
    ax.set_ylabel("Hap2 telomere length (kbp)")
    ax.set_xlim(0, max_len * 1.05)
    ax.set_ylim(0, max_len * 1.05)

    # === SCALE TOGGLE === set both to 'linear' or 'log' as needed
    # ax.set_xscale('linear'); ax.set_yscale('linear')
    # ax.set_xscale('log');    ax.set_yscale('log')

    # Legend: same font for title and entries
    ax.legend(title="Arm", frameon=True, loc="upper left", fontsize=8, title_fontsize=8)

    # R^2 and p (bottom-right inside)
    if np.isfinite(r2) and np.isfinite(pval):
        txt = f"R² = {r2:.3f}\np = {pval:.3g}\nn = {n_pairs}"
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
# FIG 2: H9 hap1 vs hap2 paired scatter with chr labels + Spearman ρ, p
# --------------------------------------------------------------------------

def make_plot_h9_hap1_vs_hap2_spearman(df: pd.DataFrame, out_dir: str, basename: str) -> None:
    sub = df[df["assembly_label"].isin(["H9 hap1", "H9 hap2"])].copy()
    if sub.empty:
        print("[WARN] No H9; skipping Fig2 (Spearman).")
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
        print("[WARN] No matched H9 telomeres; skipping Fig2 (Spearman).")
        return

    # Spearman ρ and p
    rho, pval = np.nan, np.nan
    if spearmanr is not None:
        rho, pval = spearmanr(pivot["hap1"].values, pivot["hap2"].values)

    fig, ax = plt.subplots(figsize=(4.0, 4.0), dpi=600)

    # color by arm; label each point with chromosome number
    arm_colors = {"p": "#1f77b4", "q": "#ff7f0e"}
    for arm_val, grp in pivot.groupby("arm"):
        c = arm_colors.get(arm_val, "#7f7f7f")
        ax.scatter(grp["hap1"], grp["hap2"], s=30, alpha=0.66, edgecolor="black", linewidth=0.3, label=f"{arm_val}", c=c)
        for _, row in grp.iterrows():
            num = str(row["chrom"])
            if num.startswith("chr"):
                num = num[3:]
            ax.text(row["hap1"] + 0.1, row["hap2"] + 0.1, num, fontsize=5,
                    va="bottom", ha="left", color=c)

    max_len = float(max(pivot["hap1"].max(), pivot["hap2"].max()))
    ax.plot([0, max_len], [0, max_len], linestyle="--", color="grey", linewidth=1.0)

    ax.set_xlabel("Hap1 telomere length (kbp)")
    ax.set_ylabel("Hap2 telomere length (kbp)")
    ax.set_xlim(0, max_len * 1.05)
    ax.set_ylim(0, max_len * 1.05)

    # Legend
    ax.legend(title="Arm", frameon=True, loc="upper left", fontsize=8, title_fontsize=8)

    # Spearman ρ and p (bottom-right inside)
    if np.isfinite(rho) and np.isfinite(pval):
        txt = f"ρ = {rho:.3f}\np = {pval:.3g}\nn = {len(pivot)}"
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
# MAIN
# --------------------------------------------------------------------------

def main() -> None:
    ensure_dir(OUT_DIR)
    df = parse_teloscope_terminal_file(INPUT_BED)

    # Statistical analysis
    perform_statistical_analysis(df, OUT_DIR)

    # Fig 1: H9 hap1 vs hap2, paired by chr arm, labels + Pearson R²/p
    make_plot_h9_hap1_vs_hap2_pearson(df, OUT_DIR, FIG1_BASENAME)
    
    # Fig 2: H9 hap1 vs hap2, paired by chr arm, labels + Spearman ρ/p
    make_plot_h9_hap1_vs_hap2_spearman(df, OUT_DIR, FIG2_BASENAME)

if __name__ == "__main__":
    main()
