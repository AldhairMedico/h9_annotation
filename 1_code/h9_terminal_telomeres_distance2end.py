#!/usr/bin/env python3
"""
Terminal telomere analysis for H9 assemblies.

Goal
--------------
- Explore telomere length vs distance to chromosome end.
- Generate scatterplots colored by assembly using the palette.
"""

import os
import re
from typing import Dict
import numpy as np
import pandas as pd

# --- headless-friendly backend (must be set before importing pyplot) ---
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.transforms import blended_transform_factory

# ───────────────────────────────────────────────────────────────────────────────
# Directory structure
# ───────────────────────────────────────────────────────────────────────────────
# Determine base_dir dynamically from script location
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(os.path.dirname(script_dir))
data_dir = os.path.join(base_dir, "2_data", "2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures", "3.1_draft")
run_label = "v3"  # Label for output filenames

# Create subdirectories for organized output
distance2end_dir = os.path.join(figures_dir, "25.12.14_distance2end")
chr_positioning_dir = os.path.join(figures_dir, "25.12.14_chr_positioning")
os.makedirs(distance2end_dir, exist_ok=True)
os.makedirs(chr_positioning_dir, exist_ok=True)

terminal_telomeres = os.path.join(data_dir, "25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed")
terminal_telomeres_extended = os.path.join(data_dir, "25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres_extended.bed")
terminal_telomeres_filtered = os.path.join(data_dir, "25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres_filtered.bed")

# ───────────────────────────────────────────────────────────────────────────────
# Chromosome name mapping for GRCh38 and RPE1
# ───────────────────────────────────────────────────────────────────────────────
GRCH38_CHR_MAP = {
    "CM000663.2": "chr1", "CM000664.2": "chr2", "CM000665.2": "chr3",
    "CM000666.2": "chr4", "CM000667.2": "chr5", "CM000668.2": "chr6",
    "CM000669.2": "chr7", "CM000670.2": "chr8", "CM000671.2": "chr9",
    "CM000672.2": "chr10", "CM000673.2": "chr11", "CM000674.2": "chr12",
    "CM000675.2": "chr13", "CM000676.2": "chr14", "CM000677.2": "chr15",
    "CM000678.2": "chr16", "CM000679.2": "chr17", "CM000680.2": "chr18",
    "CM000681.2": "chr19", "CM000682.2": "chr20", "CM000683.2": "chr21",
    "CM000684.2": "chr22", "CM000685.2": "chrX", "CM000686.2": "chrY",
}

RPE1_HAP2_CHR_MAP = {
    "CM116119.1": "chr1", "CM116120.1": "chr2", "CM116121.1": "chr3",
    "CM116122.1": "chr4", "CM116123.1": "chr5", "CM116124.1": "chr6",
    "CM116125.1": "chr7", "CM116126.1": "chr8", "CM116127.1": "chr9",
    "CM116128.1": "chr10", "CM116129.1": "chr11", "CM116130.1": "chr12",
    "CM116131.1": "chr13", "CM116132.1": "chr14", "CM116133.1": "chr15",
    "CM116134.1": "chr16", "CM116135.1": "chr17", "CM116136.1": "chr18",
    "CM116137.1": "chr19", "CM116138.1": "chr20", "CM116139.1": "chr21",
    "CM116140.1": "chr22", "CM116141.1": "chrX",
}

# RPE1 hap1 uses CM116095.1 to CM116117.1
RPE1_HAP1_CHR_MAP = {
    "CM116095.1": "chr1", "CM116096.1": "chr2", "CM116097.1": "chr3",
    "CM116098.1": "chr4", "CM116099.1": "chr5", "CM116100.1": "chr6",
    "CM116101.1": "chr7", "CM116102.1": "chr8", "CM116103.1": "chr9",
    "CM116104.1": "chr10", "CM116105.1": "chr11", "CM116106.1": "chr12",
    "CM116107.1": "chr13", "CM116108.1": "chr14", "CM116109.1": "chr15",
    "CM116110.1": "chr16", "CM116111.1": "chr17", "CM116112.1": "chr18",
    "CM116113.1": "chr19", "CM116114.1": "chr20", "CM116115.1": "chr21",
    "CM116116.1": "chr22", "CM116117.1": "chrX", "CM116118.1": "chrY",
}


def get_display_chr_name(chr_name: str, assembly_label: str) -> str:
    """Convert chromosome name to display format based on assembly."""
    # For H9: remove _hap1/_hap2 suffix
    if assembly_label.startswith("H9"):
        return re.sub(r'_hap[12]$', '', chr_name)

    # For GRCh38: map accession to chr name
    if assembly_label == "GRCh38":
        return GRCH38_CHR_MAP.get(chr_name, chr_name)

    # For RPE1 hap1: map accession to chr name
    if assembly_label == "RPE1 hap1":
        return RPE1_HAP1_CHR_MAP.get(chr_name, chr_name)

    # For RPE1 hap2: map accession to chr name
    if assembly_label == "RPE1 hap2":
        return RPE1_HAP2_CHR_MAP.get(chr_name, chr_name)

    return chr_name


# Read with flexible columns (some assemblies have extra "Chromosome" field)
# Handle mixed column counts (11 or 12 columns) by reading line by line
rows = []
with open(terminal_telomeres, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) == 11:
            rows.append(fields)
        elif len(fields) == 12:
            # Remove the extra "Chromosome" field (index 1)
            rows.append([fields[0]] + fields[2:])

df = pd.DataFrame(rows, columns=[
    "chr", "start", "end", "length", "label",
    "fwdCounts", "revCounts", "canCounts", "nonCanCounts",
    "chr_length", "assembly"
])

# Convert numeric columns
numeric_cols = ["start", "end", "length", "fwdCounts", "revCounts",
                "canCounts", "nonCanCounts", "chr_length"]
for col in numeric_cols:
    df[col] = pd.to_numeric(df[col])

# Palette (user-updated; do not change)
PALETTE: Dict[str, str] = {
    "GRCh38"    : "#999999",  # Neutral Gray
    "CHM13"     : "#F0E442",  # Yellow
    "YAO mat"   : "#FF5353",  # Red
    "YAO pat"   : "#FFA5A5",  # Light Red
    "HG002 mat" : "#0072B2",  # Okabe-Ito Blue
    "HG002 pat" : "#56B4E9",  # Okabe-Ito Sky Blue
    "RPE1 hap1" : "#984EA3",  # Purple
    "RPE1 hap2" : "#CC79A7",  # Magenta
    "I002C hap1": "#D55E00",  # Vermillion
    "I002C hap2": "#E69F00",  # Orange
    "H9 hap1"   : "#00796B",  # Teal 700
    "H9 hap2"   : "#80CBC4",  # Teal 200
}

# Map raw assembly names to display names for palette lookup
ASSEMBLY_MAP: Dict[str, str] = {
    "GCA_000001405.29_GRCh38.p14_genomic.chr.fna": "GRCh38",
    "GCA_009914755.4_T2T-CHM13v2.0_genomic.chr.fna": "CHM13",
    "GCA_018852605.3_hg002v1.1.pat_genomic.fna": "HG002 pat",
    "GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna": "HG002 mat",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.chr.fna": "RPE1 hap2",
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.chr.fna": "RPE1 hap1",
    "GWHDOOG00000000.genome.chr.fasta.gz": "YAO pat",
    "GWHDQZJ00000000.genome.chr.fasta.gz": "YAO mat",
    "GWHGEYB00000000.1.genome.fasta.gz": "YAO pat",
    "GWHGEYC00000000.1.genome.fasta.gz": "YAO mat",
    "H9_T2T_v0.1_hap1.fasta": "H9 hap1",
    "H9_T2T_v0.1_hap2.fasta": "H9 hap2",
    "I002Cv0.7.hap1.fasta.gz": "I002C hap1",
    "I002Cv0.7.hap2.fasta.gz": "I002C hap2",
}

# Map assembly to display name
df["assembly_label"] = df["assembly"].map(ASSEMBLY_MAP)

# Create display chromosome name
df["chr_display"] = df.apply(lambda row: get_display_chr_name(row["chr"], row["assembly_label"]), axis=1)

# Define assembly order for consistent legend ordering
ASSEMBLY_ORDER = [
    "GRCh38", "CHM13", "HG002 mat", "HG002 pat",
    "I002C hap1", "I002C hap2", "YAO mat", "YAO pat",
    "RPE1 hap1", "RPE1 hap2", "H9 hap1", "H9 hap2"
]

# ───────────────────────────────────────────────────────────────────────────────
# Compute distance to end and chromosome positioning
# ───────────────────────────────────────────────────────────────────────────────
# distance2end: minimum of "start" vs ("chr_length" - "end")
df["distance2end"] = np.minimum(df["start"], df["chr_length"] - df["end"])
df["distance2end_pct"] = (df["distance2end"] / df["chr_length"]) * 100

# chr_position: position along chromosome (0% = start, 100% = end)
df["chr_position"] = (df["start"] + df["end"]) / 2
df["chr_position_pct"] = ((df["start"] + df["end"]) / 2 / df["chr_length"]) * 100

# canonical_pct: canonical proportion as percentage (canCounts * 6 / length)
df["canonical_pct"] = (df["canCounts"] * 6 / df["length"]) * 100

# ───────────────────────────────────────────────────────────────────────────────
# Save extended telomere data with new columns
# ───────────────────────────────────────────────────────────────────────────────
df.to_csv(terminal_telomeres_extended, sep="\t", index=False, header=False)
print(f"[OK] Extended telomere data saved to {terminal_telomeres_extended}")

# ───────────────────────────────────────────────────────────────────────────────
# Save filtered telomere data (distance2end <= 100)
# ───────────────────────────────────────────────────────────────────────────────
df_filtered = df[df["distance2end"] <= 100]
df_filtered.to_csv(terminal_telomeres_filtered, sep="\t", index=False, header=False)
print(f"[OK] Filtered telomere data ({len(df_filtered)} records with distance2end <= 100) saved to {terminal_telomeres_filtered}")

# ───────────────────────────────────────────────────────────────────────────────
# Common settings for distance2end plots
# ───────────────────────────────────────────────────────────────────────────────
CUTOFF_BP = 100
CUTOFF_PCT = 0.1

# Setup for chromosome positioning plots
bar_height = 0.6  # compact
cmap = plt.cm.cividis

# Normalizers for heatmap coloring
norm_length = Normalize(vmin=df["length"].min() / 1000, vmax=df["length"].max() / 1000)
norm_canonical = Normalize(vmin=df["canonical_pct"].min(), vmax=df["canonical_pct"].max())

# Chromosome sort order (natural sort: chr1, chr2, ..., chr22, chrX, chrY)
CHR_ORDER = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

def chr_sort_key(chr_name: str) -> int:
    """Sort key for chromosome names."""
    if chr_name in CHR_ORDER:
        return CHR_ORDER.index(chr_name)
    return 100  # Unknown chromosomes go last


# ═══════════════════════════════════════════════════════════════════════════════
# DISTANCE TO END PLOTS (Figures 1-4) -> distance2end_dir
# Consistent aesthetics: black text, fontsize=8, linewidth=1
# ═══════════════════════════════════════════════════════════════════════════════

# Calculate outliers (same for both bp and pct based plots)
n_below_bp = len(df[df["distance2end"] <= CUTOFF_BP])
n_above_bp = len(df[df["distance2end"] > CUTOFF_BP])
outliers_bp = df[df["distance2end"] > CUTOFF_BP].copy()

n_below_pct = len(df[df["distance2end_pct"] <= CUTOFF_PCT])
n_above_pct = len(df[df["distance2end_pct"] > CUTOFF_PCT])
outliers_pct = df[df["distance2end_pct"] > CUTOFF_PCT].copy()

# ───────────────────────────────────────────────────────────────────────────────
# Plot 1a: End distance (log bp) vs Telomere length - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    subset = subset[subset["distance2end"] > 0]
    ax.scatter(
        subset["distance2end"],
        subset["length"] / 1000.0,
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.axvline(x=CUTOFF_BP, color="black", linestyle="--", linewidth=1)
ax.text(0.02, 0.98, f"n={n_below_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', color='black')
ax.text(0.98, 0.98, f"n={n_above_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right', color='black')
for _, row in outliers_bp.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (row["distance2end"], row["length"] / 1000.0),
        fontsize=5, alpha=0.7, color='black',
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xscale("log")
ax.set_xlabel("Distance to end (bp)")
ax.set_ylabel("Telomere length (Kbp)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"1a_Telomere_vs_distance2end_log_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 1b: End distance (log bp) vs Telomere length - Heatmap
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
df_pos = df[df["distance2end"] > 0].copy()
h, xedges, yedges, im = ax.hist2d(
    np.log10(df_pos["distance2end"]),
    df_pos["length"] / 1000.0,
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.axvline(x=np.log10(CUTOFF_BP), linestyle="--", linewidth=1, color='black')
ax.text(0.02, 0.98, f"n={n_below_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top')
ax.text(0.98, 0.98, f"n={n_above_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right')
for _, row in outliers_bp.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (np.log10(row["distance2end"]), row["length"] / 1000.0),
        fontsize=5, alpha=0.7,
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xlabel("Distance to end (log10 bp)")
ax.set_ylabel("Telomere length (Kbp)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"1b_Telomere_vs_distance2end_log_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 2a: End distance (% chr) vs Telomere length - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    ax.scatter(
        subset["distance2end_pct"],
        subset["length"] / 1000.0,
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.axvline(x=CUTOFF_PCT, color="black", linestyle="--", linewidth=1)
ax.text(0.02, 0.98, f"n={n_below_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', color='black')
ax.text(0.98, 0.98, f"n={n_above_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right', color='black')
for _, row in outliers_pct.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (row["distance2end_pct"], row["length"] / 1000.0),
        fontsize=5, alpha=0.7, color='black',
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xlabel("Distance to end (% chr)")
ax.set_ylabel("Telomere length (Kbp)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"2a_Telomere_vs_distance2end_pct_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 2b: End distance (% chr) vs Telomere length - Heatmap
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
h, xedges, yedges, im = ax.hist2d(
    df["distance2end_pct"],
    df["length"] / 1000.0,
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.axvline(x=CUTOFF_PCT, linestyle="--", linewidth=1, color='black')
ax.text(0.02, 0.98, f"n={n_below_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top')
ax.text(0.98, 0.98, f"n={n_above_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right')
for _, row in outliers_pct.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (row["distance2end_pct"], row["length"] / 1000.0),
        fontsize=5, alpha=0.7,
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xlabel("Distance to end (% chr)")
ax.set_ylabel("Telomere length (Kbp)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"2b_Telomere_vs_distance2end_pct_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 3a: End distance (log bp) vs Canonical proportion - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    subset = subset[subset["distance2end"] > 0]
    ax.scatter(
        subset["distance2end"],
        subset["canonical_pct"],
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.axvline(x=CUTOFF_BP, color="black", linestyle="--", linewidth=1)
ax.text(0.02, 0.98, f"n={n_below_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', color='black')
ax.text(0.98, 0.98, f"n={n_above_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right', color='black')
for _, row in outliers_bp.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (row["distance2end"], row["canonical_pct"]),
        fontsize=5, alpha=0.7, color='black',
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xscale("log")
ax.set_xlabel("Distance to end (bp)")
ax.set_ylabel("Canonical proportion (%)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"3a_Canonical_vs_distance2end_log_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 3b: End distance (log bp) vs Canonical proportion - Heatmap
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
df_pos = df[df["distance2end"] > 0].copy()
h, xedges, yedges, im = ax.hist2d(
    np.log10(df_pos["distance2end"]),
    df_pos["canonical_pct"],
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.axvline(x=np.log10(CUTOFF_BP), linestyle="--", linewidth=1, color='black')
ax.text(0.02, 0.98, f"n={n_below_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top')
ax.text(0.98, 0.98, f"n={n_above_bp}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right')
for _, row in outliers_bp.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (np.log10(row["distance2end"]), row["canonical_pct"]),
        fontsize=5, alpha=0.7,
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xlabel("Distance to end (log10 bp)")
ax.set_ylabel("Canonical proportion (%)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"3b_Canonical_vs_distance2end_log_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 4a: End distance (% chr) vs Canonical proportion - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    ax.scatter(
        subset["distance2end_pct"],
        subset["canonical_pct"],
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.axvline(x=CUTOFF_PCT, color="black", linestyle="--", linewidth=1)
ax.text(0.02, 0.98, f"n={n_below_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', color='black')
ax.text(0.98, 0.98, f"n={n_above_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right', color='black')
for _, row in outliers_pct.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (row["distance2end_pct"], row["canonical_pct"]),
        fontsize=5, alpha=0.7, color='black',
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xlabel("Distance to end (% chr)")
ax.set_ylabel("Canonical proportion (%)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"4a_Canonical_vs_distance2end_pct_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 4b: End distance (% chr) vs Canonical proportion - Heatmap
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
h, xedges, yedges, im = ax.hist2d(
    df["distance2end_pct"],
    df["canonical_pct"],
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.axvline(x=CUTOFF_PCT, linestyle="--", linewidth=1, color='black')
ax.text(0.02, 0.98, f"n={n_below_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top')
ax.text(0.98, 0.98, f"n={n_above_pct}", transform=ax.transAxes,
        fontsize=8, verticalalignment='top', horizontalalignment='right')
for _, row in outliers_pct.iterrows():
    ax.annotate(
        f"{row['chr_display']}",
        (row["distance2end_pct"], row["canonical_pct"]),
        fontsize=5, alpha=0.7,
        xytext=(3, 3), textcoords="offset points"
    )
ax.set_xlabel("Distance to end (% chr)")
ax.set_ylabel("Canonical proportion (%)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(distance2end_dir, f"4b_Canonical_vs_distance2end_pct_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)