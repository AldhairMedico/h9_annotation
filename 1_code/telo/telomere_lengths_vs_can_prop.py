#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
telomere_lengths_comparison.py

Using Teloscope terminal telomere annotations from:
  2_data/2.2_processed/25.12.10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed

This script:
  1) Fig 1: Telomere length (Kbp) vs canonical proportion (%) scatter
     (square, y-axis starts at min(canonical)-1, legend inside).
  2) Fig 2: H9-only (hap1 & hap2) scatter of telomere length vs canonical proportion
     with chromosome name labels highlighted.

All PNG saved at 600 dpi; PNG + PDF go to ./figures.
"""

import os
import math
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from adjustText import adjust_text

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------

# Get the repo root directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

INPUT_BED = os.path.join(REPO_ROOT, "2_data", "2.2_processed", "25.12.10_teloscope_compiled", "25.12.10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed")

# Output dir
OUT_DIR = os.path.join(REPO_ROOT, "3_figures", "3.1_draft", "26.01.29_telomeres")
os.makedirs(OUT_DIR, exist_ok=True)

# Figure basenames (v5)
FIG1_BASENAME = "scatter_by_assembly"
FIG2_BASENAME = "scatter_H9_length_vs_canonical"

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
# FIG 1: scatter (square, y>=min-1)
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
    fig.savefig(os.path.join(out_dir, f"{basename}.svg"))
    plt.close(fig)

# --------------------------------------------------------------------------
# FIG 2: H9-only scatter with chromosome labels (length vs canonical %)
# --------------------------------------------------------------------------

# Exclusion lists for trimmed/removed telomeres (not included in regression)
EXCLUDE_HAP1 = {"21p", "15p", "12p", "12q", "16p"}
EXCLUDE_HAP2 = {"13p", "15q", "2q", "21p"}

def make_plot_h9_length_vs_canonical_scatter(df: pd.DataFrame, out_dir: str, basename: str) -> None:
    """
    Scatter plot of telomere length (Kbp) vs canonical proportion (%)
    for H9 hap1 and hap2 only, with chromosome name labels on each point.
    Excluded telomeres are shown as hollow points and not included in regression.
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

    # Determine which points are excluded based on haplotype
    def is_excluded(row) -> bool:
        if row["assembly_label"] == "H9 hap1":
            return row["chr_label"] in EXCLUDE_HAP1
        elif row["assembly_label"] == "H9 hap2":
            return row["chr_label"] in EXCLUDE_HAP2
        return False
    sub["excluded"] = sub.apply(is_excluded, axis=1)

    fig, ax = plt.subplots(figsize=(5, 5), dpi=600)

    texts = []  # For adjustText
    labels_colors = []
    for asm in ["H9 hap1", "H9 hap2"]:
        asm_sub = sub[sub["assembly_label"] == asm]
        if asm_sub.empty:
            continue
        col = PALETTE.get(asm, "#BBBBBB")
        labels_colors.append((asm, col))

        # Included points (filled)
        included = asm_sub[~asm_sub["excluded"]]
        if not included.empty:
            ax.scatter(included["tel_length_kbp"], included["canonical_pct"], s=30, alpha=0.85,
                       edgecolor="black", linewidth=0.5, c=col, zorder=3)

        # Excluded points (hollow with grey border)
        excluded = asm_sub[asm_sub["excluded"]]
        if not excluded.empty:
            ax.scatter(excluded["tel_length_kbp"], excluded["canonical_pct"], s=30, alpha=0.85,
                       edgecolor="grey", linewidth=0.5, facecolors="none", zorder=2)

        # Add chromosome labels for all points
        for _, row in asm_sub.iterrows():
            txt = ax.text(row["tel_length_kbp"], row["canonical_pct"],
                          row["chr_label"], fontsize=7, va="center", ha="center",
                          color=col if not row["excluded"] else "grey", zorder=4)
            texts.append(txt)

    # Adjust text positions to avoid overlapping
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.3),
                expand_points=(1.8, 1.8), force_text=(0.6, 0.6))

    # Regression on included points only
    included_all = sub[~sub["excluded"]].dropna(subset=["tel_length_kbp", "canonical_pct"])
    if len(included_all) >= 3:
        x_reg = included_all["tel_length_kbp"].values
        y_reg = included_all["canonical_pct"].values
        # Spearman correlation
        rho, p_val = stats.spearmanr(x_reg, y_reg)
        # Linear regression for the line
        slope, intercept = np.polyfit(x_reg, y_reg, 1)
        x_line = np.linspace(0, 25, 100)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, color="black", linestyle="--", linewidth=1.2, alpha=0.7, zorder=1)
        # Display Spearman rho and p-value at the top
        ax.text(0.05, 0.95, f"Spearman ρ = {rho:.3f}\np-adj = {p_val:.2e}",
                transform=ax.transAxes, fontsize=8, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="lightgrey", alpha=0.9))

    ax.set_xlabel("Telomere length (Kbp)", fontsize=11)
    ax.set_ylabel("Canonical proportion (%)", fontsize=11)
    # ax.set_xscale("log")
    ax.set_xlim(0, 23)
    y_min_scatter = min(102, float(np.nanmin(sub["canonical_pct"])) - 1.0)
    ax.set_ylim(y_min_scatter, 102)
    ax.tick_params(axis="both", which="major", labelsize=10, length=4, width=1.0, direction="out")
    ax.tick_params(axis="both", which="minor", labelsize=8, length=2, width=0.8, direction="out")
    ax.grid(True, which="major", axis="both", linewidth=0.6, alpha=0.35, color="lightgrey")

    legend_inside_br(ax, labels_colors)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    plt.tight_layout()
    ensure_dir(out_dir)
    fig.savefig(os.path.join(out_dir, f"{basename}.png"), dpi=600)
    fig.savefig(os.path.join(out_dir, f"{basename}.pdf"))
    fig.savefig(os.path.join(out_dir, f"{basename}.svg"))
    plt.close(fig)

# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------

def main() -> None:
    ensure_dir(OUT_DIR)
    df = parse_teloscope_terminal_file(INPUT_BED)

    # Fig 1: scatter (square), y>=min-1
    make_plot_length_vs_canonical_scatter(df, OUT_DIR, FIG1_BASENAME)

    # Fig 2: H9-only scatter with chr labels (length vs canonical %)
    make_plot_h9_length_vs_canonical_scatter(df, OUT_DIR, FIG2_BASENAME)

if __name__ == "__main__":
    main()
