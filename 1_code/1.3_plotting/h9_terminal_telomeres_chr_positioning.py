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
df["distance2end"] = np.minimum(df["start"], df["chr_length"] - df["end"])
df["distance2end_pct"] = (df["distance2end"] / df["chr_length"]) * 100

df["chr_position"] = (df["start"] + df["end"]) / 2
df["chr_position_pct"] = ((df["start"] + df["end"]) / 2 / df["chr_length"]) * 100

df["canonical_pct"] = (df["canCounts"] * 6 / df["length"]) * 100

# ───────────────────────────────────────────────────────────────────────────────
# Save extended telomere data with new columns
# ───────────────────────────────────────────────────────────────────────────────
df.to_csv(terminal_telomeres_extended, sep="\t", index=False, header=False)
print(f"[OK] Extended telomere data saved to {terminal_telomeres_extended}")

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
cmap = plt.cm.viridis

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
# CHROMOSOME POSITIONING PLOTS -> chr_positioning_dir
# ═══════════════════════════════════════════════════════════════════════════════

# Constants for arm filtering
ARM_CUTOFF_BP = 25000  # 25 Kbp from ends
ARM_CUTOFF_PCT = 0.05  # 0.05% from ends
ARM_CUTOFF_KB = ARM_CUTOFF_BP / 1000  # 25 Kbp

# Separate p-arm and q-arm telomeres
df_p_arm = df[df["label"] == "p"].copy()
df_q_arm = df[df["label"] == "q"].copy()

# --- sort order: assembly first, then chr_name (natural), then start ---
ASM_SORT = {asm: i for i, asm in enumerate(ASSEMBLY_ORDER)}
df_p_arm["assembly_sort"] = df_p_arm["assembly_label"].map(ASM_SORT).fillna(999).astype(int)
df_q_arm["assembly_sort"] = df_q_arm["assembly_label"].map(ASM_SORT).fillna(999).astype(int)

df_p_arm["chr_sort"] = df_p_arm["chr_display"].apply(chr_sort_key)
df_q_arm["chr_sort"] = df_q_arm["chr_display"].apply(chr_sort_key)

df_p_sorted = df_p_arm.sort_values(["assembly_sort", "chr_sort", "chr_display", "start"]).reset_index(drop=True)
df_q_sorted = df_q_arm.sort_values(["assembly_sort", "chr_sort", "chr_display", "start"]).reset_index(drop=True)

# --- build aligned p/q tables so both panels have identical row count/order ---
KEY_COLS = ["assembly_sort", "assembly_label", "chr_sort", "chr_display", "chr_length"]

base_rows = (
    pd.concat([df_p_sorted[KEY_COLS], df_q_sorted[KEY_COLS]], ignore_index=True)
      .drop_duplicates(subset=KEY_COLS)
      .sort_values(["assembly_sort", "chr_sort", "chr_display"])
      .reset_index(drop=True)
)

df_p_plot = base_rows.merge(df_p_sorted, on=KEY_COLS, how="left")
df_q_plot = base_rows.merge(df_q_sorted, on=KEY_COLS, how="left")

df_p_plot["is_placeholder"] = df_p_plot["start"].isna()
df_q_plot["is_placeholder"] = df_q_plot["start"].isna()

_fill_cols = ["start", "end", "length", "canonical_pct"]
df_p_plot[_fill_cols] = df_p_plot[_fill_cols].fillna(0)
df_q_plot[_fill_cols] = df_q_plot[_fill_cols].fillna(0)

_max_n = len(base_rows)

# --- smaller axis fonts for figs 5–6 ---
AX_TICK_FONTSIZE = 7
AX_LABEL_FONTSIZE = 8
TITLE_FONTSIZE = 12

def style_xaxis(ax):
    ax.tick_params(axis="x", labelsize=AX_TICK_FONTSIZE)
    ax.xaxis.label.set_size(AX_LABEL_FONTSIZE)
    ax.title.set_size(TITLE_FONTSIZE)

def add_assembly_bar(ax, df_sorted: pd.DataFrame, ylim, palette: Dict[str, str]) -> None:
    """
    Left-side annotation strip:
      - left: color (palette)
      - right: assembly name (right-aligned), not drawn on top of color
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(*ylim)

    color_w = 0.18  # fraction reserved for color
    # per-row color + white label background
    for i, asm in enumerate(df_sorted["assembly_label"].tolist()):
        ax.barh(i, color_w, left=0, height=bar_height, color=palette.get(asm, "#CCCCCC"), edgecolor="none")
        ax.barh(i, 1 - color_w, left=color_w, height=bar_height, color="white", edgecolor="none")

    # block labels (one per contiguous assembly block), right-aligned in the white region
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    asm_list = df_sorted["assembly_label"].tolist()
    if asm_list:
        start = 0
        curr = asm_list[0]
        for i in range(1, len(asm_list) + 1):
            nxt = asm_list[i] if i < len(asm_list) else None
            if nxt != curr:
                mid = (start + (i - 1)) / 2
                ax.text(
                    0.98, mid, str(curr),
                    transform=trans, va="center", ha="right",
                    fontsize=AX_LABEL_FONTSIZE, color="black"
                )
                start = i
                curr = nxt

    ax.axis("off")

# ───────────────────────────────────────────────────────────────────────────────
# Plot 5A: Relative positioning near ends (%chr)
# ───────────────────────────────────────────────────────────────────────────────
fig, (ax_asm, ax_p, ax_q) = plt.subplots(
    1, 3, figsize=(6.45, 4), sharey=True, constrained_layout=True,
    gridspec_kw={"width_ratios": [0.40, 1, 1]}
)

# p-arm
for i, (_, row) in enumerate(df_p_plot.iterrows()):
    ax_p.barh(i, ARM_CUTOFF_PCT, height=bar_height, color="#E0E0E0", edgecolor="none")
    if not row["is_placeholder"]:
        start_pct = row["start"] / row["chr_length"] * 100
        end_pct = row["end"] / row["chr_length"] * 100
        if end_pct > start_pct:
            ax_p.barh(i, end_pct - start_pct, left=start_pct, height=bar_height,
                      color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_p.set_xlim(0, ARM_CUTOFF_PCT)
ax_p.set_ylim(-0.5, _max_n - 0.5)
ax_p.set_xlabel("Distance to end (%chr)")
ax_p.set_title("p-arm")
ax_p.set_yticks([])
style_xaxis(ax_p)
sns.despine(ax=ax_p, top=True, right=True, left=True)

# q-arm
for i, (_, row) in enumerate(df_q_plot.iterrows()):
    ax_q.barh(i, ARM_CUTOFF_PCT, left=100 - ARM_CUTOFF_PCT, height=bar_height, color="#E0E0E0", edgecolor="none")
    if not row["is_placeholder"]:
        start_pct = row["start"] / row["chr_length"] * 100
        end_pct = row["end"] / row["chr_length"] * 100
        if end_pct > start_pct:
            ax_q.barh(i, end_pct - start_pct, left=start_pct, height=bar_height,
                      color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_q.set_xlim(100 - ARM_CUTOFF_PCT, 100)
ax_q.set_ylim(-0.5, _max_n - 0.5)
ax_q.set_xlabel("Distance to end (%chr)")
ax_q.set_title("q-arm")
ax_q.set_yticks([])
style_xaxis(ax_q)
sns.despine(ax=ax_q, top=True, right=True, left=True)

add_assembly_bar(ax_asm, base_rows, ylim=(-0.5, _max_n - 0.5), palette=PALETTE)

for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(chr_positioning_dir, f"5A_chr_positioning_pct_{run_label}.{ext}"),
                dpi=600 if ext == "png" else None)
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 5B: Absolute near ends (Kbp) (q-arm shown as 25->0 Kbp distance-to-end)
# ───────────────────────────────────────────────────────────────────────────────
fig, (ax_asm, ax_p, ax_q) = plt.subplots(
    1, 3, figsize=(6.45, 4), sharey=True, constrained_layout=True,
    gridspec_kw={"width_ratios": [0.40, 1, 1]}
)

for i, (_, row) in enumerate(df_p_plot.iterrows()):
    ax_p.barh(i, ARM_CUTOFF_KB, height=bar_height, color="#E0E0E0", edgecolor="none")
    if not row["is_placeholder"]:
        start_kb = row["start"] / 1000.0
        end_kb = row["end"] / 1000.0
        if end_kb > start_kb:
            ax_p.barh(i, end_kb - start_kb, left=start_kb, height=bar_height,
                      color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_p.set_xlim(0, ARM_CUTOFF_KB)
ax_p.set_ylim(-0.5, _max_n - 0.5)
ax_p.set_xlabel("Distance to end (Kbp)")
ax_p.set_title("p-arm")
ax_p.set_yticks([])
style_xaxis(ax_p)
sns.despine(ax=ax_p, top=True, right=True, left=True)

for i, (_, row) in enumerate(df_q_plot.iterrows()):
    ax_q.barh(i, ARM_CUTOFF_KB, left=0, height=bar_height, color="#E0E0E0", edgecolor="none")
    if not row["is_placeholder"]:
        d_end_kb = (row["chr_length"] - row["end"]) / 1000.0
        d_start_kb = (row["chr_length"] - row["start"]) / 1000.0
        if d_start_kb > d_end_kb:
            ax_q.barh(i, d_start_kb - d_end_kb, left=d_end_kb, height=bar_height,
                      color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_q.set_xlim(ARM_CUTOFF_KB, 0)
ax_q.set_ylim(-0.5, _max_n - 0.5)
ax_q.set_xlabel("Distance to end (Kbp)")
ax_q.set_title("q-arm")
ax_q.set_yticks([])
style_xaxis(ax_q)
sns.despine(ax=ax_q, top=True, right=True, left=True)

add_assembly_bar(ax_asm, base_rows, ylim=(-0.5, _max_n - 0.5), palette=PALETTE)

for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(chr_positioning_dir, f"5B_chr_positioning_Kbp_{run_label}.{ext}"),
                dpi=600 if ext == "png" else None)
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 6A: Heatmap - Telomere length × %chr near ends
# p: 0..0.05 ; q: 99.95..100 (NOT inverted)
# ───────────────────────────────────────────────────────────────────────────────
fig, (ax_asm, ax_p, ax_q) = plt.subplots(
    1, 3, figsize=(6.45, 4), sharey=True, constrained_layout=True,
    gridspec_kw={"width_ratios": [0.40, 1, 1]}
)

for i, (_, row) in enumerate(df_p_sorted.iterrows()):
    ax_p.barh(i, ARM_CUTOFF_PCT, height=bar_height, color="#E0E0E0", edgecolor="none")
    start_pct = row["start"] / row["chr_length"] * 100
    end_pct = row["end"] / row["chr_length"] * 100
    ax_p.barh(i, end_pct - start_pct, left=start_pct, height=bar_height,
              color=cmap(norm_length(row["length"] / 1000.0)), edgecolor="none")

ax_p.set_xlim(0, ARM_CUTOFF_PCT)
ax_p.set_ylim(-0.5, _max_n - 0.5)
ax_p.set_xlabel("Distance to end (%chr)")
ax_p.set_title("p-arm")
ax_p.set_yticks([])
style_xaxis(ax_p)
sns.despine(ax=ax_p, top=True, right=True, left=True)

for i, (_, row) in enumerate(df_q_sorted.iterrows()):
    ax_q.barh(i, ARM_CUTOFF_PCT, left=100 - ARM_CUTOFF_PCT, height=bar_height, color="#E0E0E0", edgecolor="none")
    start_pct = row["start"] / row["chr_length"] * 100
    end_pct = row["end"] / row["chr_length"] * 100
    ax_q.barh(i, end_pct - start_pct, left=start_pct, height=bar_height,
              color=cmap(norm_length(row["length"] / 1000.0)), edgecolor="none")

ax_q.set_xlim(100 - ARM_CUTOFF_PCT, 100)
ax_q.set_ylim(-0.5, _max_n - 0.5)
ax_q.set_xlabel("Distance to end (%chr)")
ax_q.set_title("q-arm")
ax_q.set_yticks([])
style_xaxis(ax_q)
sns.despine(ax=ax_q, top=True, right=True, left=True)

_bar_df = df_p_sorted if len(df_p_sorted) >= len(df_q_sorted) else df_q_sorted
add_assembly_bar(ax_asm, _bar_df, ylim=(-0.5, _max_n - 0.5), palette=PALETTE)

sm = ScalarMappable(cmap=cmap, norm=norm_length)
sm.set_array([])
cb = fig.colorbar(sm, ax=[ax_p, ax_q], fraction=0.04, pad=0.02)
cb.set_label("Telomere length (Kbp)", fontsize=AX_LABEL_FONTSIZE)
cb.ax.tick_params(labelsize=AX_TICK_FONTSIZE)

for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(chr_positioning_dir, f"6A_Length_heatmap_pct_{run_label}.{ext}"),
                dpi=600 if ext == "png" else None)
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 6B: Heatmap - Telomere length × Kbp near ends
# p: 0..25 ; q: distance-to-q-end, axis 25..0 (inverted)
# ───────────────────────────────────────────────────────────────────────────────
fig, (ax_asm, ax_p, ax_q) = plt.subplots(
    1, 3, figsize=(6.45, 4), sharey=True, constrained_layout=True,
    gridspec_kw={"width_ratios": [0.40, 1, 1]}
)

for i, (_, row) in enumerate(df_p_sorted.iterrows()):
    ax_p.barh(i, ARM_CUTOFF_KB, height=bar_height, color="#E0E0E0", edgecolor="none")
    start_kb = row["start"] / 1000.0
    end_kb = row["end"] / 1000.0
    ax_p.barh(i, end_kb - start_kb, left=start_kb, height=bar_height,
              color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_p.set_xlim(0, ARM_CUTOFF_KB)
ax_p.set_ylim(-0.5, _max_n - 0.5)
ax_p.set_xlabel("Distance to end (Kbp)")
ax_p.set_title("p-arm")
ax_p.set_yticks([])
style_xaxis(ax_p)
sns.despine(ax=ax_p, top=True, right=True, left=True)

# q-arm: distance from chromosome end (Kbp)
for i, (_, row) in enumerate(df_q_sorted.iterrows()):
    ax_q.barh(i, ARM_CUTOFF_KB, left=0, height=bar_height, color="#E0E0E0", edgecolor="none")
    d_end_kb = (row["chr_length"] - row["end"]) / 1000.0
    d_start_kb = (row["chr_length"] - row["start"]) / 1000.0
    ax_q.barh(i, d_start_kb - d_end_kb, left=d_end_kb, height=bar_height,
              color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_q.set_xlim(ARM_CUTOFF_KB, 0)  # 25 -> 0 (expected)
ax_q.set_ylim(-0.5, _max_n - 0.5)
ax_q.set_xlabel("Distance to end (Kbp)")
ax_q.set_title("q-arm")
ax_q.set_yticks([])
style_xaxis(ax_q)
sns.despine(ax=ax_q, top=True, right=True, left=True)

_bar_df = df_p_sorted if len(df_p_sorted) >= len(df_q_sorted) else df_q_sorted
add_assembly_bar(ax_asm, _bar_df, ylim=(-0.5, _max_n - 0.5), palette=PALETTE)

sm = ScalarMappable(cmap=cmap, norm=norm_length)
sm.set_array([])
cb = fig.colorbar(sm, ax=[ax_p, ax_q], fraction=0.04, pad=0.02)
cb.set_label("Telomere length (Kbp)", fontsize=AX_LABEL_FONTSIZE)
cb.ax.tick_params(labelsize=AX_TICK_FONTSIZE)

for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(chr_positioning_dir, f"6B_Length_heatmap_Kbp_{run_label}.{ext}"),
                dpi=600 if ext == "png" else None)
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 6C: Heatmap - Canonical proportion × %chr near ends
# p: 0..0.05 ; q: 99.95..100 (NOT inverted)
# ───────────────────────────────────────────────────────────────────────────────
fig, (ax_asm, ax_p, ax_q) = plt.subplots(
    1, 3, figsize=(6.45, 4), sharey=True, constrained_layout=True,
    gridspec_kw={"width_ratios": [0.40, 1, 1]}
)

for i, (_, row) in enumerate(df_p_sorted.iterrows()):
    ax_p.barh(i, ARM_CUTOFF_PCT, height=bar_height, color="#E0E0E0", edgecolor="none")
    start_pct = row["start"] / row["chr_length"] * 100
    end_pct = row["end"] / row["chr_length"] * 100
    ax_p.barh(i, end_pct - start_pct, left=start_pct, height=bar_height,
              color=cmap(norm_canonical(row["canonical_pct"])), edgecolor="none")

ax_p.set_xlim(0, ARM_CUTOFF_PCT)
ax_p.set_ylim(-0.5, _max_n - 0.5)
ax_p.set_xlabel("Distance to end (%chr)")
ax_p.set_title("p-arm")
ax_p.set_yticks([])
style_xaxis(ax_p)
sns.despine(ax=ax_p, top=True, right=True, left=True)

for i, (_, row) in enumerate(df_q_sorted.iterrows()):
    ax_q.barh(i, ARM_CUTOFF_PCT, left=100 - ARM_CUTOFF_PCT, height=bar_height, color="#E0E0E0", edgecolor="none")
    start_pct = row["start"] / row["chr_length"] * 100
    end_pct = row["end"] / row["chr_length"] * 100
    ax_q.barh(i, end_pct - start_pct, left=start_pct, height=bar_height,
              color=cmap(norm_canonical(row["canonical_pct"])), edgecolor="none")

ax_q.set_xlim(100 - ARM_CUTOFF_PCT, 100)
ax_q.set_ylim(-0.5, _max_n - 0.5)
ax_q.set_xlabel("Distance to end (%chr)")
ax_q.set_title("q-arm")
ax_q.set_yticks([])
style_xaxis(ax_q)
sns.despine(ax=ax_q, top=True, right=True, left=True)

_bar_df = df_p_sorted if len(df_p_sorted) >= len(df_q_sorted) else df_q_sorted
add_assembly_bar(ax_asm, _bar_df, ylim=(-0.5, _max_n - 0.5), palette=PALETTE)

sm_can = ScalarMappable(cmap=cmap, norm=norm_canonical)
sm_can.set_array([])
cb = fig.colorbar(sm_can, ax=[ax_p, ax_q], fraction=0.04, pad=0.02)
cb.set_label("Canonical proportion (%)", fontsize=AX_LABEL_FONTSIZE)
cb.ax.tick_params(labelsize=AX_TICK_FONTSIZE)

for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(chr_positioning_dir, f"6C_Canonical_heatmap_pct_{run_label}.{ext}"),
                dpi=600 if ext == "png" else None)
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 6D: Heatmap - Canonical proportion × Kbp near ends
# p: 0..25 ; q: distance-to-q-end, axis 25..0 (inverted)
# ───────────────────────────────────────────────────────────────────────────────
fig, (ax_asm, ax_p, ax_q) = plt.subplots(
    1, 3, figsize=(6.45, 4), sharey=True, constrained_layout=True,
    gridspec_kw={"width_ratios": [0.40, 1, 1]}
)

for i, (_, row) in enumerate(df_p_sorted.iterrows()):
    ax_p.barh(i, ARM_CUTOFF_KB, height=bar_height, color="#E0E0E0", edgecolor="none")
    start_kb = row["start"] / 1000.0
    end_kb = row["end"] / 1000.0
    ax_p.barh(i, end_kb - start_kb, left=start_kb, height=bar_height,
              color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_p.set_xlim(0, ARM_CUTOFF_KB)
ax_p.set_ylim(-0.5, _max_n - 0.5)
ax_p.set_xlabel("Distance to end (Kbp)")
ax_p.set_title("p-arm")
ax_p.set_yticks([])
style_xaxis(ax_p)
sns.despine(ax=ax_p, top=True, right=True, left=True)

for i, (_, row) in enumerate(df_q_sorted.iterrows()):
    ax_q.barh(i, ARM_CUTOFF_KB, left=0, height=bar_height, color="#E0E0E0", edgecolor="none")
    d_end_kb = (row["chr_length"] - row["end"]) / 1000.0
    d_start_kb = (row["chr_length"] - row["start"]) / 1000.0
    ax_q.barh(i, d_start_kb - d_end_kb, left=d_end_kb, height=bar_height,
              color=PALETTE[row["assembly_label"]], edgecolor="none", alpha=0.9)

ax_q.set_xlim(ARM_CUTOFF_KB, 0)
ax_q.set_ylim(-0.5, _max_n - 0.5)
ax_q.set_xlabel("Distance to end (Kbp)")
ax_q.set_title("q-arm")
ax_q.set_yticks([])
style_xaxis(ax_q)
sns.despine(ax=ax_q, top=True, right=True, left=True)

_bar_df = df_p_sorted if len(df_p_sorted) >= len(df_q_sorted) else df_q_sorted
add_assembly_bar(ax_asm, _bar_df, ylim=(-0.5, _max_n - 0.5), palette=PALETTE)

sm_can = ScalarMappable(cmap=cmap, norm=norm_canonical)
sm_can.set_array([])
cb = fig.colorbar(sm_can, ax=[ax_p, ax_q], fraction=0.04, pad=0.02)
cb.set_label("Canonical proportion (%)", fontsize=AX_LABEL_FONTSIZE)
cb.ax.tick_params(labelsize=AX_TICK_FONTSIZE)

for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(chr_positioning_dir, f"6D_Canonical_heatmap_Kbp_{run_label}.{ext}"),
                dpi=600 if ext == "png" else None)
plt.close(fig)

print(f"[OK] Distance2end plots (1-4) saved to {distance2end_dir}")
print(f"[OK] Chr positioning plots (5-6) saved to {chr_positioning_dir}")
