#!/usr/bin/env python3
"""
H9 terminal telomere lengths by chromosome and haplotype.

Generates two compact mirrored barplots of telomere length (kbp):
  1. Horizontal (← p-arm | chr labels | q-arm →)
  2. Vertical   (↑ p-arm | chr labels | q-arm ↓)

One tick per chromosome; hap1/hap2 are side-by-side thin bars
distinguished by colour (legend).
"""

import os
from typing import Dict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ───────────────────────────────────────────────────────────────────────────────
# Paths
# ───────────────────────────────────────────────────────────────────────────────
script_dir  = os.path.dirname(os.path.abspath(__file__))
base_dir    = os.path.dirname(os.path.dirname(script_dir))
data_dir    = os.path.join(base_dir, "2_data", "2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures", "3.1_draft")
run_label   = "v1"

plots_dir = os.path.join(figures_dir, "26.01.29_telomeres")
os.makedirs(plots_dir, exist_ok=True)

terminal_telomeres = os.path.join(
    data_dir, "25.12.10_teloscope_compiled",
    "25.12.10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed",
)

# ───────────────────────────────────────────────────────────────────────────────
# Load data & filter to H9
# ───────────────────────────────────────────────────────────────────────────────
ASSEMBLY_MAP: Dict[str, str] = {
    "H9_T2T_v0.1_hap1.fasta": "H9 hap1",
    "H9_T2T_v0.1_hap2.fasta": "H9 hap2",
}

rows = []
with open(terminal_telomeres, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) == 11:
            rows.append(fields)
        elif len(fields) == 12:
            rows.append([fields[0]] + fields[2:])

df = pd.DataFrame(rows, columns=[
    "chr", "start", "end", "length", "label",
    "fwdCounts", "revCounts", "canCounts", "nonCanCounts",
    "chr_length", "assembly",
])

for col in ("start", "end", "length", "chr_length"):
    df[col] = pd.to_numeric(df[col])

df["assembly_label"] = df["assembly"].map(ASSEMBLY_MAP)
df = df.loc[df["assembly_label"].notna()].copy()

df["chr_display"] = df["chr"].str.replace(r"_hap[12]$", "", regex=True)
df["hap"]         = df["assembly_label"].str.replace("H9 ", "")
df["length_kbp"]  = df["length"] / 1000.0

# ───────────────────────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────────────────────
PALETTE = {"hap1": "#00796B", "hap2": "#80CBC4"}

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def _chr_key(c: str) -> int:
    return CHR_ORDER.index(c) if c in CHR_ORDER else 100


df_p = df.loc[df["label"] == "p"]
df_q = df.loc[df["label"] == "q"]

EXTS  = ("png", "pdf", "svg")
TICK  = 5
LABEL = 6
TITLE = 7

LEGEND_HANDLES = [
    Patch(fc=PALETTE["hap1"], label="hap1", linewidth=0),
    Patch(fc=PALETTE["hap2"], label="hap2", linewidth=0),
]
LEGEND_KW = dict(
    handles=LEGEND_HANDLES, fontsize=TICK, loc="upper left",
    frameon=False, handlelength=0.8, handleheight=0.6, handletextpad=0.4,
    borderpad=0.2, labelspacing=0.3,
)

x_max = float(np.ceil(df["length_kbp"].max() / 5) * 5)

# Global thin spines / ticks
plt.rcParams.update({
    "axes.linewidth":    0.4,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
    "xtick.major.size":  2.0,
    "ytick.major.size":  2.0,
})

# ───────────────────────────────────────────────────────────────────────────────
# Helpers
# ───────────────────────────────────────────────────────────────────────────────


def _save(fig, stem: str) -> None:
    for ext in EXTS:
        fig.savefig(
            os.path.join(plots_dir, f"{stem}.{ext}"),
            dpi=600 if ext == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(f"[OK] {stem}")


def _box(ax):
    for sp in ax.spines.values():
        sp.set_edgecolor("black")
        sp.set_linewidth(0.4)
        sp.set_visible(True)


def _len(arm_df: pd.DataFrame, hap: str, chrom: str) -> float:
    """Telomere length (kbp) for a given haplotype + chromosome."""
    m = (arm_df["hap"] == hap) & (arm_df["chr_display"] == chrom)
    v = arm_df.loc[m, "length_kbp"]
    return v.iloc[0] if len(v) else 0.0


# ───────────────────────────────────────────────────────────────────────────────
# Build chromosome list (sorted, deduplicated)
# ───────────────────────────────────────────────────────────────────────────────
chroms = sorted(df["chr_display"].unique(), key=_chr_key)
n_chr  = len(chroms)

# Short chromosome labels for vertical plot (1, 2, … 22, X)
chr_short = [c.replace("chr", "") for c in chroms]

# Grouped-bar geometry  (two thin bars per chromosome)
bw     = 0.30           # bar half-width (thin)
offsets = {"hap1": -bw / 2 - 0.02, "hap2": bw / 2 + 0.02}
pad    = 0.5            # padding so first/last chr aren't clipped

# ═══════════════════════════════════════════════════════════════════════════════
# Horizontal layout  (← p-arm | chr labels | q-arm →)
# ═══════════════════════════════════════════════════════════════════════════════
fig, (ax_p, ax_q) = plt.subplots(
    1, 2, figsize=(2, 4), sharey=True, constrained_layout=True,
)

yy = np.arange(n_chr)

for hap in ("hap1", "hap2"):
    p_vals = [_len(df_p, hap, c) for c in chroms]
    q_vals = [_len(df_q, hap, c) for c in chroms]
    ax_p.barh(yy + offsets[hap], p_vals, height=bw,
              color=PALETTE[hap], edgecolor="none", alpha=0.9)
    ax_q.barh(yy + offsets[hap], q_vals, height=bw,
              color=PALETTE[hap], edgecolor="none", alpha=0.9)

# p-arm: bars grow LEFT
ax_p.set_xlim(x_max, 0)
ax_p.set_ylim(n_chr - 1 + pad, -pad)  # tight, inverted
ax_p.set_yticks(yy)
ax_p.set_yticklabels(chroms, fontsize=TICK)
ax_p.yaxis.tick_right()
ax_p.yaxis.set_label_position("right")
ax_p.set_xlabel("Telomere length (kbp)", fontsize=LABEL)
ax_p.set_title("p-arm", fontsize=TITLE)
ax_p.tick_params(axis="x", labelsize=TICK)
# ax_p.tick_params(axis="y", length=0)  # no tick marks on shared axis
ax_p.legend(**LEGEND_KW)
_box(ax_p)

# q-arm: bars grow RIGHT
ax_q.set_xlim(0, x_max)
ax_q.set_ylim(n_chr - 1 + pad, -pad)
ax_q.set_xlabel("Telomere length (kbp)", fontsize=LABEL)
ax_q.set_title("q-arm", fontsize=TITLE)
ax_q.tick_params(axis="x", labelsize=TICK)
_box(ax_q)

_save(fig, f"h9_telo_lengths_horizontal_{run_label}")

# ═══════════════════════════════════════════════════════════════════════════════
# Vertical layout  (↑ p-arm | chr labels | q-arm ↓)
# ═══════════════════════════════════════════════════════════════════════════════
fig, (ax_p, ax_q) = plt.subplots(
    2, 1, figsize=(4, 2), sharex=True, constrained_layout=True,
)

xx = np.arange(n_chr)

for hap in ("hap1", "hap2"):
    p_vals = [_len(df_p, hap, c) for c in chroms]
    q_vals = [_len(df_q, hap, c) for c in chroms]
    ax_p.bar(xx + offsets[hap], p_vals, width=bw,
             color=PALETTE[hap], edgecolor="none", alpha=0.9)
    ax_q.bar(xx + offsets[hap], q_vals, width=bw,
             color=PALETTE[hap], edgecolor="none", alpha=0.9)

# p-arm: bars extend UP
ax_p.set_xlim(-pad, n_chr - 1 + pad)
ax_p.set_ylim(0, x_max)
ax_p.set_ylabel("Telomere length (kbp)", fontsize=LABEL)
ax_p.set_title("p-arm", fontsize=TITLE)
ax_p.tick_params(axis="y", labelsize=TICK)
ax_p.legend(**LEGEND_KW)
_box(ax_p)

# q-arm: bars extend DOWN
ax_q.set_xlim(-pad, n_chr - 1 + pad)
ax_q.set_ylim(x_max, 0)
ax_q.set_ylabel("Telomere length (kbp)", fontsize=LABEL)
ax_q.set_title("q-arm", fontsize=TITLE)
ax_q.tick_params(axis="y", labelsize=TICK)
ax_q.set_xticks(xx)
ax_q.set_xticklabels(chr_short, fontsize=TICK, rotation=0)
_box(ax_q)

_save(fig, f"h9_telo_lengths_vertical_{run_label}")

print("[DONE] H9 telomere length plots generated.")
