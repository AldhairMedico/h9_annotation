#!/usr/bin/env python3
"""
Terminal telomere analysis for H9 assemblies.

Goal
--------------
- Explore telomere length vs distance to chromosome end.
- Generate scatterplots colored by assembly using the palette.
"""

import os
from typing import Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ───────────────────────────────────────────────────────────────────────────────
# Directory structure
# ───────────────────────────────────────────────────────────────────────────────
base_dir = "/mnt/d/research/h9_annotation"
data_dir = os.path.join(base_dir, "2_data/2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures/3.1_draft")
run_label = "v1"  # Label for output filenames

terminal_telomeres = os.path.join(data_dir, "25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed")
terminal_telomeres_extended = os.path.join(data_dir, "25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres_extended.bed")
terminal_telomeres_filtered = os.path.join(data_dir, "25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres_filtered.bed")

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
    "YAO pat"   : "#FF5353",
    "YAO mat"   : "#FFA5A5",
    "HG002 mat" : "#0072B2",  # Okabe–Ito Blue
    "HG002 pat" : "#56B4E9",  # Okabe–Ito Sky Blue
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

# Define assembly order for consistent legend ordering
ASSEMBLY_ORDER = [
    "GRCh38", "CHM13", "HG002 mat", "HG002 pat",
    "YAO mat", "YAO pat", "RPE1 hap1", "RPE1 hap2",
    "I002C hap1", "I002C hap2", "H9 hap1", "H9 hap2"
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
# Plot 1: End distance (log) vs Telomere length - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    # Filter out zero/negative distances for log scale
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
ax.set_xscale("log")
ax.set_xlabel("Distance to end (bp)")
ax.set_ylabel("Telomere length (Kbp)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"1a_Telomere_vs_distance2end_log_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 1b: End distance (log) vs Telomere length - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
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
ax.set_xlabel("Distance to end (log₁₀ bp)")
ax.set_ylabel("Telomere length (Kbp)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"1b_Telomere_vs_distance2end_log_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 2a: End distance (% chr) vs Telomere length - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
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
ax.set_xlabel("Distance to end (%chr)")
ax.set_ylabel("Telomere length (Kbp)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"2a_Telomere_vs_distance2end_pct_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 2b: End distance (% chr) vs Telomere length - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
h, xedges, yedges, im = ax.hist2d(
    df["distance2end_pct"],
    df["length"] / 1000.0,
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.set_xlabel("Distance to end (%chr)")
ax.set_ylabel("Telomere length (Kbp)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"2b_Telomere_vs_distance2end_pct_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 3a: Chromosome positioning vs Telomere length - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    ax.scatter(
        subset["chr_position_pct"],
        subset["length"] / 1000.0,
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.set_xlabel("Chromosome positioning (%chr)")
ax.set_ylabel("Telomere length (Kbp)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"3a_Telomere_vs_chr_position_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 3b: Chromosome positioning vs Telomere length - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
h, xedges, yedges, im = ax.hist2d(
    df["chr_position_pct"],
    df["length"] / 1000.0,
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.set_xlabel("Chromosome positioning (%chr)")
ax.set_ylabel("Telomere length (Kbp)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"3b_Telomere_vs_chr_position_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 4a: End distance (log) vs Canonical proportion - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    # Filter out zero/negative distances for log scale
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
ax.set_xscale("log")
ax.set_xlabel("Distance to end (bp)")
ax.set_ylabel("Canonical proportion (%)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"4a_Canonical_vs_distance2end_log_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 4b: End distance (log) vs Canonical proportion - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
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
ax.set_xlabel("Distance to end (log₁₀ bp)")
ax.set_ylabel("Canonical proportion (%)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"4b_Canonical_vs_distance2end_log_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 5a: End distance (% chr) vs Canonical proportion - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
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
ax.set_xlabel("Distance to end (%chr)")
ax.set_ylabel("Canonical proportion (%)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"5a_Canonical_vs_distance2end_pct_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 5b: End distance (% chr) vs Canonical proportion - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
h, xedges, yedges, im = ax.hist2d(
    df["distance2end_pct"],
    df["canonical_pct"],
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.set_xlabel("Distance to end (%chr)")
ax.set_ylabel("Canonical proportion (%)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"5b_Canonical_vs_distance2end_pct_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 6a: Chromosome positioning vs Canonical proportion - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    ax.scatter(
        subset["chr_position_pct"],
        subset["canonical_pct"],
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.set_xlabel("Chromosome positioning (%chr)")
ax.set_ylabel("Canonical proportion (%)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"6a_Canonical_vs_chr_position_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 6b: Chromosome positioning vs Canonical proportion - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
h, xedges, yedges, im = ax.hist2d(
    df["chr_position_pct"],
    df["canonical_pct"],
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.set_xlabel("Chromosome positioning (%chr)")
ax.set_ylabel("Canonical proportion (%)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"6b_Canonical_vs_chr_position_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 7a: Chromosome position (log bp) vs Telomere length - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    # Filter out zero/negative positions for log scale
    subset = subset[subset["chr_position"] > 0]
    ax.scatter(
        subset["chr_position"],
        subset["length"] / 1000.0,
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.set_xscale("log")
ax.set_xlabel("Chromosome position (bp)")
ax.set_ylabel("Telomere length (Kbp)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"7a_Telomere_vs_chr_position_log_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 7b: Chromosome position (log bp) vs Telomere length - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
df_pos = df[df["chr_position"] > 0].copy()
h, xedges, yedges, im = ax.hist2d(
    np.log10(df_pos["chr_position"]),
    df_pos["length"] / 1000.0,
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.set_xlabel("Chromosome position (log₁₀ bp)")
ax.set_ylabel("Telomere length (Kbp)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"7b_Telomere_vs_chr_position_log_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 8a: Chromosome position (log bp) vs Canonical proportion - Scatter by assembly
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
for asm in ASSEMBLY_ORDER:
    subset = df[df["assembly_label"] == asm]
    if subset.empty:
        continue
    # Filter out zero/negative positions for log scale
    subset = subset[subset["chr_position"] > 0]
    ax.scatter(
        subset["chr_position"],
        subset["canonical_pct"],
        c=PALETTE[asm],
        label=asm,
        alpha=0.6,
        s=10,
        edgecolors="none"
    )
ax.set_xscale("log")
ax.set_xlabel("Chromosome position (bp)")
ax.set_ylabel("Canonical proportion (%)")
ax.legend(title="Assembly", fontsize=7, title_fontsize=8, markerscale=2)
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"8a_Canonical_vs_chr_position_log_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

# ───────────────────────────────────────────────────────────────────────────────
# Plot 8b: Chromosome position (log bp) vs Canonical proportion - Heatmap (all data)
# ───────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5.5))
df_pos = df[df["chr_position"] > 0].copy()
h, xedges, yedges, im = ax.hist2d(
    np.log10(df_pos["chr_position"]),
    df_pos["canonical_pct"],
    bins=100,
    cmap="hot",
    cmin=1
)
cb = plt.colorbar(im, ax=ax)
cb.set_label("Count")
ax.set_xlabel("Chromosome position (log₁₀ bp)")
ax.set_ylabel("Canonical proportion (%)")
sns.despine(ax=ax, top=True, right=True)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(
        os.path.join(figures_dir, f"8b_Canonical_vs_chr_position_log_heatmap_{run_label}.{ext}"),
        dpi=600 if ext == "png" else None
    )
plt.close(fig)

print(f"[OK] Plots saved to {figures_dir}")
