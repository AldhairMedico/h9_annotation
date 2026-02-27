#!/usr/bin/env python3
"""
H9 terminal telomere lengths with NucFlag orthogonal validation.

Generates:
  Plot 1 – Per-chromosome stacked bars (telomere length coloured by NucFlag
           category), mirrored p / q layout.  Variants: {ONT, HiFi} × {H, V}.
  Plot 2 – Summary stacked bars per arm (one tick per haplotype, total telomere
           length coloured by NucFlag category).  Variants: {ONT, HiFi} × {H, V}.
  Plot 3 – Condensed summary (arms combined, one bar per haplotype).
           Variants: {ONT, HiFi} × {H, V}.
  Plot 4 – a) NucFlag composition heatmap (fraction per chr × category),
           b) ONT-vs-HiFi proportion scatter,
           c) Chromosome-type category proportions (metacentric / submetacentric
              / acrocentric; grouped bar).
  Stats  – Comprehensive TSV with per-chr/arm/hap breakdowns, enrichment tests,
           chromosome-type comparisons, hypergeometric tests, ONT-vs-HiFi.

Requires running  1_code/telo/get_nucflag_teloscope.sh  first to produce:
    2_data/2.2_processed/nucflag/nucflag_teloscope_hifi.bed
    2_data/2.2_processed/nucflag/nucflag_teloscope_ont.bed
"""

import os
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcolors

# ═══════════════════════════════════════════════════════════════════════════════
# Paths
# ═══════════════════════════════════════════════════════════════════════════════
script_dir  = os.path.dirname(os.path.abspath(__file__))
base_dir    = os.path.dirname(os.path.dirname(script_dir))
data_dir    = os.path.join(base_dir, "2_data", "2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures", "3.1_draft")
run_label   = "v1"

plots_dir = os.path.join(figures_dir, "26.02.25_telomeres_nucflag")
os.makedirs(plots_dir, exist_ok=True)

nucflag_dir = os.path.join(data_dir, "nucflag")

terminal_telomeres = os.path.join(
    data_dir, "25.12.10_teloscope_compiled",
    "25.12.10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed",
)

NUCFLAG_INTERSECTION = {
    "HiFi": os.path.join(nucflag_dir, "nucflag_teloscope_hifi.bed"),
    "ONT":  os.path.join(nucflag_dir, "nucflag_teloscope_ont.bed"),
}

stats_out = os.path.join(plots_dir, f"h9_nucflag_telomere_stats_{run_label}.tsv")

# ═══════════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════════
ASSEMBLY_MAP: Dict[str, str] = {
    "H9_T2T_v0.1_hap1.fasta": "H9 hap1",
    "H9_T2T_v0.1_hap2.fasta": "H9 hap2",
}

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
ACROCENTRIC    = {"chr13", "chr14", "chr15", "chr21", "chr22"}
METACENTRIC    = {"chr1", "chr3", "chr16", "chr19", "chr20"}
SUBMETACENTRIC = {"chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                  "chr10", "chr11", "chr12", "chr17", "chr18", "chrX"}
# Note: humans have no telocentric chromosomes

def _chr_key(c: str) -> int:
    return CHR_ORDER.index(c) if c in CHR_ORDER else 100

# NucFlag categories (union of HiFi + ONT observed labels)
NUCFLAG_CATS = [
    "correct",
    "misjoin",
    "collapse",
    "false_dup",
    "het_or_mismap",
    "low_quality",
    "deletion",
    "insertion",
    "mismatch",
    "softclip",
    "homopolymer",
    "dinucleotide",
    "simple_repeat",
    "other_repeat",
    "scaffold",
]

# Publication-quality palette (Cell journal style)
NUCFLAG_COLORS: Dict[str, str] = {
    "correct":       "#2D8E47",  # deep green
    "misjoin":       "#CB3335",  # crimson
    "collapse":      "#7B5EA7",  # violet
    "false_dup":     "#8B6D52",  # sienna
    "het_or_mismap": "#C77CBA",  # orchid
    "low_quality":   "#E69F00",  # amber (colorblind-safe)
    "deletion":      "#3C7DC4",  # cobalt
    "insertion":     "#00B6AD",  # teal
    "mismatch":      "#D4A535",  # goldenrod
    "softclip":      "#5A6F80",  # pewter
    "homopolymer":   "#7CB950",  # leaf green
    "dinucleotide":  "#B5CC3E",  # chartreuse
    "simple_repeat": "#E07B42",  # burnt orange
    "other_repeat":  "#47A697",  # sea green
    "scaffold":      "#999999",  # gray
    "uncovered":     "#D4D4D4",  # light gray
}

HAP_ORDER = ["hap1", "hap2"]

EXTS = ("png", "pdf", "svg")

# Cell journal aesthetics (Arial, minimal spines, white background)
plt.rcParams.update({
    "font.family":       "Arial",
    "font.size":         7,
    "axes.linewidth":    0.4,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
    "xtick.major.size":  2.0,
    "ytick.major.size":  2.0,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "savefig.facecolor": "white",
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
})

TICK  = 6   # ticks, arm labels, legend, other text
LABEL = 7   # X / Y axis labels

# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

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
    """Cell-style axes: only left + bottom spines."""
    for side in ("left", "bottom"):
        ax.spines[side].set_edgecolor("black")
        ax.spines[side].set_linewidth(0.4)
        ax.spines[side].set_visible(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Load telomere data  (same logic as original script)
# ═══════════════════════════════════════════════════════════════════════════════
rows: list = []
with open(terminal_telomeres, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) == 11:
            rows.append(fields)
        elif len(fields) == 12:
            rows.append([fields[0]] + fields[2:])

df_telo = pd.DataFrame(rows, columns=[
    "chr", "start", "end", "length", "label",
    "fwdCounts", "revCounts", "canCounts", "nonCanCounts",
    "chr_length", "assembly",
])
for col in ("start", "end", "length", "chr_length"):
    df_telo[col] = pd.to_numeric(df_telo[col])

df_telo["assembly_label"] = df_telo["assembly"].map(ASSEMBLY_MAP)
df_telo = df_telo.loc[df_telo["assembly_label"].notna()].copy()
df_telo["chr_display"] = df_telo["chr"].str.replace(r"_hap[12]$", "", regex=True)
df_telo["hap"]         = df_telo["assembly_label"].str.replace("H9 ", "")
df_telo["length_kbp"]  = df_telo["length"] / 1000.0
df_telo["arm"]         = df_telo["label"]  # p / q

# Chromosome list (deduplicated, sorted)
chroms   = sorted(df_telo["chr_display"].unique(), key=_chr_key)
n_chr    = len(chroms)
chr_short = [c.replace("chr", "") for c in chroms]

print(f"[INFO] Loaded {len(df_telo)} H9 telomere entries across {n_chr} chromosomes.")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Load NucFlag intersection data
# ═══════════════════════════════════════════════════════════════════════════════
def load_nucflag_intersection(path: str) -> pd.DataFrame:
    """Parse bedtools intersect -wo output (11 columns)."""
    col_names = [
        "telo_chr", "telo_start", "telo_end", "telo_length", "arm", "assembly",
        "nf_chr", "nf_start", "nf_end", "nf_category", "overlap_bp",
    ]
    df = pd.read_csv(path, sep="\t", header=None, names=col_names)
    for c in ("telo_start", "telo_end", "telo_length", "nf_start", "nf_end", "overlap_bp"):
        df[c] = pd.to_numeric(df[c])
    df["chr_display"] = df["telo_chr"].str.replace(r"_hap[12]$", "", regex=True)
    df["hap"] = df["assembly"].map({
        "H9_T2T_v0.1_hap1.fasta": "hap1",
        "H9_T2T_v0.1_hap2.fasta": "hap2",
    })
    df["overlap_kbp"] = df["overlap_bp"] / 1000.0
    return df


nf_data: Dict[str, pd.DataFrame] = {}
for ds, path in NUCFLAG_INTERSECTION.items():
    nf_data[ds] = load_nucflag_intersection(path)
    print(f"[INFO] {ds}: {len(nf_data[ds])} intersection rows loaded.")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Build per-telomere NucFlag breakdown
# ═══════════════════════════════════════════════════════════════════════════════
def build_nucflag_matrix(nf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot to one row per (chr_display, arm, hap) with columns for each
    NucFlag category giving overlap_kbp.  An 'uncovered' column captures
    any telomere length not accounted for by NucFlag.
    """
    agg = (
        nf_df.groupby(["chr_display", "arm", "hap", "nf_category"])["overlap_bp"]
        .sum()
        .reset_index()
    )
    piv = agg.pivot_table(
        index=["chr_display", "arm", "hap"],
        columns="nf_category",
        values="overlap_bp",
        fill_value=0,
    ).reset_index()

    # Total covered
    cat_cols = [c for c in piv.columns if c not in ("chr_display", "arm", "hap")]
    piv["covered_bp"] = piv[cat_cols].sum(axis=1)

    # Get telomere length for each entry to compute uncovered
    telo_len = (
        nf_df.groupby(["chr_display", "arm", "hap"])["telo_length"]
        .first()
        .reset_index()
    )
    piv = piv.merge(telo_len, on=["chr_display", "arm", "hap"], how="left")
    piv["uncovered"] = (piv["telo_length"] - piv["covered_bp"]).clip(lower=0)

    # Convert to kbp
    all_cats = cat_cols + ["uncovered"]
    for c in all_cats:
        piv[c] = piv[c] / 1000.0
    piv["telo_length_kbp"] = piv["telo_length"] / 1000.0

    return piv


nf_matrix: Dict[str, pd.DataFrame] = {}
for ds, nf_df in nf_data.items():
    nf_matrix[ds] = build_nucflag_matrix(nf_df)
    print(f"[INFO] {ds}: NucFlag matrix  {nf_matrix[ds].shape}")


# Categories actually present across all datasets
all_present_cats: List[str] = []
for cat in NUCFLAG_CATS + ["uncovered"]:
    for mat in nf_matrix.values():
        if cat in mat.columns and mat[cat].sum() > 0:
            if cat not in all_present_cats:
                all_present_cats.append(cat)
            break

print(f"[INFO] NucFlag categories present: {all_present_cats}")

# ═══════════════════════════════════════════════════════════════════════════════
# Legend handles (shared)
# ═══════════════════════════════════════════════════════════════════════════════
def _legend_handles(cats: List[str]) -> List[Patch]:
    return [
        Patch(fc=NUCFLAG_COLORS.get(c, "#D4D4D4"), label=c, linewidth=0)
        for c in cats
    ]


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Per-chromosome stacked bars  (mirrored p / q)
# ═══════════════════════════════════════════════════════════════════════════════
bw      = 0.30
offsets = {"hap1": -bw / 2 - 0.02, "hap2": bw / 2 + 0.02}
pad     = 0.5


def _get_cat_val(mat: pd.DataFrame, hap: str, chrom: str, arm: str, cat: str) -> float:
    m = (mat["hap"] == hap) & (mat["chr_display"] == chrom) & (mat["arm"] == arm)
    rows = mat.loc[m]
    if rows.empty or cat not in rows.columns:
        return 0.0
    return float(rows[cat].iloc[0])


def _compute_xmax(mat: pd.DataFrame, cats: List[str]) -> float:
    """Compute nice upper limit for the stacked bar axis."""
    if "telo_length_kbp" in mat.columns:
        raw_max = mat["telo_length_kbp"].max()
    else:
        raw_max = mat[cats].sum(axis=1).max()
    return float(np.ceil(raw_max / 5) * 5)


def plot1_horizontal(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    """Horizontal mirrored layout:  ← p-arm | chr labels | q-arm →"""
    fig, (ax_p, ax_q) = plt.subplots(
        1, 2, figsize=(2, 4), sharey=True, constrained_layout=True,
    )
    yy = np.arange(n_chr)
    x_max = _compute_xmax(mat, cats)

    for arm_label, ax, invert in [("p", ax_p, True), ("q", ax_q, False)]:
        for hap in HAP_ORDER:
            bottoms = np.zeros(n_chr)
            for cat in cats:
                vals = np.array([_get_cat_val(mat, hap, c, arm_label, cat) for c in chroms])
                ax.barh(
                    yy + offsets[hap], vals, left=bottoms, height=bw,
                    color=NUCFLAG_COLORS.get(cat, "#D4D4D4"),
                    edgecolor="none", alpha=0.9,
                )
                bottoms += vals
        if invert:
            ax.set_xlim(x_max, 0)
        else:
            ax.set_xlim(0, x_max)
        ax.set_ylim(n_chr - 1 + pad, -pad)
        ax.set_xlabel("Telomere length (kbp)", fontsize=LABEL)
        ax.set_title(f"{arm_label}-arm", fontsize=TICK)
        ax.tick_params(axis="x", labelsize=TICK)
        _box(ax)

    ax_p.set_yticks(yy)
    ax_p.set_yticklabels(chroms, fontsize=TICK)
    ax_p.yaxis.tick_right()
    ax_p.yaxis.set_label_position("right")

    # Legend (outside, below figure)
    handles = _legend_handles(cats) + [
        Patch(fc="white", ec="black", lw=0.4, label="hap1 (left bar)"),
        Patch(fc="white", ec="black", lw=0.4, label="hap2 (right bar)"),
    ]
    fig.legend(handles=handles, fontsize=TICK - 1,
              loc="lower center", bbox_to_anchor=(0.5, -0.08),
              ncol=min(len(handles), 4), frameon=False,
              handlelength=0.8, handleheight=0.6,
              handletextpad=0.4, borderpad=0.2, labelspacing=0.25)

    _save(fig, f"fig1_h9_nucflag_horizontal_{ds}_{run_label}")


def plot1_vertical(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    """Vertical mirrored layout:  ↑ p-arm | chr labels | q-arm ↓"""
    fig, (ax_p, ax_q) = plt.subplots(
        2, 1, figsize=(4, 2), sharex=True, constrained_layout=True,
    )
    xx = np.arange(n_chr)
    x_max = _compute_xmax(mat, cats)

    for arm_label, ax, invert in [("p", ax_p, False), ("q", ax_q, True)]:
        for hap in HAP_ORDER:
            bottoms = np.zeros(n_chr)
            for cat in cats:
                vals = np.array([_get_cat_val(mat, hap, c, arm_label, cat) for c in chroms])
                ax.bar(
                    xx + offsets[hap], vals, bottom=bottoms, width=bw,
                    color=NUCFLAG_COLORS.get(cat, "#D4D4D4"),
                    edgecolor="none", alpha=0.9,
                )
                bottoms += vals
        if invert:
            ax.set_ylim(x_max, 0)
        else:
            ax.set_ylim(0, x_max)
        ax.set_xlim(-pad, n_chr - 1 + pad)
        ax.set_ylabel("Telomere length (kbp)", fontsize=LABEL)
        ax.set_title(f"{arm_label}-arm", fontsize=TICK)
        ax.tick_params(axis="y", labelsize=TICK)
        _box(ax)

    ax_q.set_xticks(xx)
    ax_q.set_xticklabels(chr_short, fontsize=TICK, rotation=0)

    # Legend (outside, below figure)
    handles = _legend_handles(cats) + [
        Patch(fc="white", ec="black", lw=0.4, label="hap1 (left bar)"),
        Patch(fc="white", ec="black", lw=0.4, label="hap2 (right bar)"),
    ]
    fig.legend(handles=handles, fontsize=TICK - 1,
              loc="lower center", bbox_to_anchor=(0.5, -0.12),
              ncol=min(len(handles), 6), frameon=False,
              handlelength=0.8, handleheight=0.6,
              handletextpad=0.4, borderpad=0.2, labelspacing=0.25)

    _save(fig, f"fig1_h9_nucflag_vertical_{ds}_{run_label}")


print("\n── Plot 1: per-chromosome stacked bars ──")
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    plot1_horizontal(mat, cats_ds, ds)
    plot1_vertical(mat, cats_ds, ds)

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Summary stacked bars  (one tick per haplotype)
# ═══════════════════════════════════════════════════════════════════════════════
def _hap_totals(mat: pd.DataFrame, cats: List[str], arm: str) -> Dict[str, np.ndarray]:
    """Total kbp per NucFlag category aggregated across all chromosomes."""
    out = {}
    for hap in HAP_ORDER:
        sub = mat.loc[(mat["hap"] == hap) & (mat["arm"] == arm)]
        vals = np.array([sub[c].sum() if c in sub.columns else 0.0 for c in cats])
        out[hap] = vals
    return out


def plot2_horizontal(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    fig, (ax_p, ax_q) = plt.subplots(
        1, 2, figsize=(3, 1.8), sharey=True, constrained_layout=True,
    )
    yy = np.arange(len(HAP_ORDER))
    bar_h = 0.5

    for arm_label, ax, invert in [("p", ax_p, True), ("q", ax_q, False)]:
        totals = _hap_totals(mat, cats, arm_label)
        for i, hap in enumerate(HAP_ORDER):
            left = 0.0
            for j, cat in enumerate(cats):
                v = totals[hap][j]
                ax.barh(
                    yy[i], v, left=left, height=bar_h,
                    color=NUCFLAG_COLORS.get(cat, "#D4D4D4"),
                    edgecolor="none", alpha=0.9,
                )
                left += v

        x_max_s = max(sum(totals[h]) for h in HAP_ORDER) * 1.05
        if invert:
            ax.set_xlim(x_max_s, 0)
        else:
            ax.set_xlim(0, x_max_s)
        ax.set_ylim(len(HAP_ORDER) - 0.5 - 0.3, -0.5 + 0.3)
        ax.set_xlabel("Total length (kbp)", fontsize=LABEL)
        ax.set_title(f"{arm_label}-arm", fontsize=TICK)
        ax.tick_params(axis="x", labelsize=TICK)
        _box(ax)

    ax_p.set_yticks(yy)
    ax_p.set_yticklabels(HAP_ORDER, fontsize=TICK)
    ax_p.yaxis.tick_right()
    ax_p.yaxis.set_label_position("right")

    handles = _legend_handles(cats)
    fig.legend(handles=handles, fontsize=TICK - 1,
              loc="lower center", bbox_to_anchor=(0.5, -0.18),
              ncol=min(len(handles), 4), frameon=False,
              handlelength=0.6, handleheight=0.5,
              handletextpad=0.3, borderpad=0.2, labelspacing=0.2)

    _save(fig, f"fig2_h9_nucflag_summary_horizontal_{ds}_{run_label}")


def plot2_vertical(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    fig, (ax_p, ax_q) = plt.subplots(
        2, 1, figsize=(1.8, 3), sharex=True, constrained_layout=True,
    )
    xx = np.arange(len(HAP_ORDER))
    bar_w = 0.5

    for arm_label, ax, invert in [("p", ax_p, False), ("q", ax_q, True)]:
        totals = _hap_totals(mat, cats, arm_label)
        for i, hap in enumerate(HAP_ORDER):
            bottom = 0.0
            for j, cat in enumerate(cats):
                v = totals[hap][j]
                ax.bar(
                    xx[i], v, bottom=bottom, width=bar_w,
                    color=NUCFLAG_COLORS.get(cat, "#D4D4D4"),
                    edgecolor="none", alpha=0.9,
                )
                bottom += v

        y_max_s = max(sum(totals[h]) for h in HAP_ORDER) * 1.05
        if invert:
            ax.set_ylim(y_max_s, 0)
        else:
            ax.set_ylim(0, y_max_s)
        ax.set_xlim(-0.5, len(HAP_ORDER) - 0.5)
        ax.set_ylabel("Total length (kbp)", fontsize=LABEL)
        ax.set_title(f"{arm_label}-arm", fontsize=TICK)
        ax.tick_params(axis="y", labelsize=TICK)
        _box(ax)

    ax_q.set_xticks(xx)
    ax_q.set_xticklabels(HAP_ORDER, fontsize=TICK, rotation=0)

    handles = _legend_handles(cats)
    fig.legend(handles=handles, fontsize=TICK - 1,
              loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, handlelength=0.6, handleheight=0.5,
              handletextpad=0.3, borderpad=0.2, labelspacing=0.2)

    _save(fig, f"fig2_h9_nucflag_summary_vertical_{ds}_{run_label}")


print("\n── Plot 2: summary stacked bars ──")
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    plot2_horizontal(mat, cats_ds, ds)
    plot2_vertical(mat, cats_ds, ds)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Condensed summary  (one bar per haplotype, arms combined)
# ═══════════════════════════════════════════════════════════════════════════════
def _hap_totals_combined(mat: pd.DataFrame, cats: List[str]) -> Dict[str, np.ndarray]:
    """Total kbp per NucFlag category aggregated across all chromosomes + arms."""
    out = {}
    for hap in HAP_ORDER:
        sub = mat.loc[mat["hap"] == hap]
        vals = np.array([sub[c].sum() if c in sub.columns else 0.0 for c in cats])
        out[hap] = vals
    return out


def plot3_horizontal(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    """Horizontal condensed: two bars (hap1, hap2), stacked by NucFlag."""
    fig, ax = plt.subplots(figsize=(3, 1.2), constrained_layout=True)
    yy = np.arange(len(HAP_ORDER))
    bar_h = 0.5
    totals = _hap_totals_combined(mat, cats)

    for i, hap in enumerate(HAP_ORDER):
        left = 0.0
        for j, cat in enumerate(cats):
            v = totals[hap][j]
            ax.barh(
                yy[i], v, left=left, height=bar_h,
                color=NUCFLAG_COLORS.get(cat, "#D4D4D4"),
                edgecolor="none", alpha=0.9,
            )
            left += v

    ax.set_yticks(yy)
    ax.set_yticklabels(HAP_ORDER, fontsize=TICK)
    ax.set_xlabel("Total telomere length (kbp)", fontsize=LABEL)
    ax.tick_params(axis="x", labelsize=TICK)
    _box(ax)

    handles = _legend_handles(cats)
    fig.legend(handles=handles, fontsize=TICK - 1,
              loc="lower center", bbox_to_anchor=(0.5, -0.25),
              ncol=min(len(handles), 4), frameon=False,
              handlelength=0.6, handleheight=0.5,
              handletextpad=0.3, borderpad=0.2, labelspacing=0.2)

    _save(fig, f"fig3_h9_nucflag_condensed_horizontal_{ds}_{run_label}")


def plot3_vertical(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    """Vertical condensed: two bars (hap1, hap2), stacked by NucFlag."""
    fig, ax = plt.subplots(figsize=(1.2, 3), constrained_layout=True)
    xx = np.arange(len(HAP_ORDER))
    bar_w = 0.5
    totals = _hap_totals_combined(mat, cats)

    for i, hap in enumerate(HAP_ORDER):
        bottom = 0.0
        for j, cat in enumerate(cats):
            v = totals[hap][j]
            ax.bar(
                xx[i], v, bottom=bottom, width=bar_w,
                color=NUCFLAG_COLORS.get(cat, "#D4D4D4"),
                edgecolor="none", alpha=0.9,
            )
            bottom += v

    ax.set_xticks(xx)
    ax.set_xticklabels(HAP_ORDER, fontsize=TICK, rotation=0)
    ax.set_ylabel("Total telomere length (kbp)", fontsize=LABEL)
    ax.tick_params(axis="y", labelsize=TICK)
    _box(ax)

    handles = _legend_handles(cats)
    fig.legend(handles=handles, fontsize=TICK - 1,
              loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, handlelength=0.6, handleheight=0.5,
              handletextpad=0.3, borderpad=0.2, labelspacing=0.2)

    _save(fig, f"fig3_h9_nucflag_condensed_vertical_{ds}_{run_label}")


print("\n── Plot 3: condensed summary bars ──")
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    plot3_horizontal(mat, cats_ds, ds)
    plot3_vertical(mat, cats_ds, ds)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 4 — Enrichment / comparison panels  (ex-Plot 3)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 4: enrichment & comparisons ──")

# ── 4a  NucFlag composition heatmap (fraction per chr) ──────────────────────
def plot4a_heatmap(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    """Heatmap: rows = chr×arm, columns = NucFlag categories (fraction)."""
    # Sum across haplotypes per chr × arm
    grp = mat.groupby(["chr_display", "arm"])[cats].sum()
    grp["total"] = grp.sum(axis=1)
    frac = grp[cats].div(grp["total"], axis=0).fillna(0)

    # Order rows: p then q for each chr
    row_order = []
    for c in chroms:
        for a in ("p", "q"):
            if (c, a) in frac.index:
                row_order.append((c, a))
    frac = frac.loc[row_order]
    row_labels = [f"{c} {a}" for c, a in row_order]

    # Filter out empty columns
    nonzero_cats = [c for c in cats if frac[c].sum() > 0]
    heat = frac[nonzero_cats].values

    fig, ax = plt.subplots(figsize=(2, 4), constrained_layout=True)
    im = ax.imshow(heat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(len(nonzero_cats)))
    ax.set_xticklabels(nonzero_cats, fontsize=TICK - 1, rotation=60, ha="right")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=TICK - 1)
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=TICK - 1)
    cbar.set_label("Fraction", fontsize=LABEL)
    _box(ax)
    _save(fig, f"fig4_h9_nucflag_heatmap_{ds}_{run_label}")


# ── 4b  ONT vs HiFi proportion scatter ─────────────────────────────────────
def plot4b_scatter() -> None:
    """Scatter: fraction of each NucFlag category — ONT vs HiFi."""
    # Pool across all chr×arm×hap
    fracs = {}
    for ds, mat in nf_matrix.items():
        cats_ds = [c for c in all_present_cats if c in mat.columns]
        sums = mat[cats_ds].sum()
        total = sums.sum()
        fracs[ds] = (sums / total).reindex(all_present_cats, fill_value=0)

    fig, ax = plt.subplots(figsize=(3, 3), constrained_layout=True)
    for cat in all_present_cats:
        x = fracs.get("HiFi", pd.Series(dtype=float)).get(cat, 0)
        y = fracs.get("ONT", pd.Series(dtype=float)).get(cat, 0)
        ax.scatter(x, y, color=NUCFLAG_COLORS.get(cat, "#D4D4D4"), s=30,
                   edgecolors="black", linewidths=0.3, zorder=3)
        ax.annotate(cat, (x, y), fontsize=TICK - 1.5, textcoords="offset points",
                    xytext=(3, 3))

    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], ls="--", lw=0.5, color="grey", zorder=1)
    ax.set_xlabel("HiFi fraction", fontsize=LABEL)
    ax.set_ylabel("ONT fraction", fontsize=LABEL)
    ax.tick_params(labelsize=TICK)
    _box(ax)
    _save(fig, f"fig4_h9_nucflag_ont_vs_hifi_{run_label}")


# ── 4c  Chromosome-type category proportions ─────────────────────────────
CHRTYPE_COLORS = {
    "metacentric":    "#3C7DC4",  # cobalt
    "submetacentric": "#E69F00",  # amber
    "acrocentric":    "#CB3335",  # crimson
}
# Note: humans have no telocentric chromosomes.

def _classify_chr(c: str) -> str:
    if c in ACROCENTRIC:
        return "acrocentric"
    if c in METACENTRIC:
        return "metacentric"
    if c in SUBMETACENTRIC:
        return "submetacentric"
    return "other"

def plot4c_chrtype(mat: pd.DataFrame, cats: List[str], ds: str) -> None:
    """Grouped bar: NucFlag category proportions by chromosome type
    (metacentric, submetacentric, acrocentric)."""
    mat = mat.copy()
    mat["chr_type"] = mat["chr_display"].map(_classify_chr)

    chr_types = ["metacentric", "submetacentric", "acrocentric"]
    cats_use = [c for c in cats if c in mat.columns and mat[c].sum() > 0]

    fracs = {}
    for ct in chr_types:
        sub = mat.loc[mat["chr_type"] == ct]
        sums = sub[cats_use].sum()
        total = sums.sum()
        fracs[ct] = sums / total if total > 0 else sums * 0

    x = np.arange(len(cats_use))
    n_types = len(chr_types)
    w = 0.8 / n_types

    fig, ax = plt.subplots(figsize=(0.5 * len(cats_use) + 1.5, 3), constrained_layout=True)
    for i, ct in enumerate(chr_types):
        offset = (i - (n_types - 1) / 2) * w
        ax.bar(x + offset, [fracs[ct].get(c, 0) for c in cats_use], width=w,
               color=CHRTYPE_COLORS[ct], alpha=0.85, label=ct.capitalize())

    ax.set_xticks(x)
    ax.set_xticklabels(cats_use, fontsize=TICK - 1, rotation=60, ha="right")
    ax.set_ylabel("Fraction of total telomere length", fontsize=LABEL)
    ax.tick_params(axis="y", labelsize=TICK)
    ax.legend(fontsize=TICK, frameon=False)
    _box(ax)
    _save(fig, f"fig4_h9_nucflag_chrtype_{ds}_{run_label}")


for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    plot4a_heatmap(mat, cats_ds, ds)
    plot4c_chrtype(mat, cats_ds, ds)

if len(nf_matrix) == 2:
    plot4b_scatter()


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Comprehensive statistics  →  TSV
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Statistics ──")

stat_rows: list = []


def _add(section: str, dataset: str, key: str, value) -> None:
    stat_rows.append({"section": section, "dataset": dataset, "key": key, "value": value})


# ── 4a  Per chr × arm × hap breakdown ──────────────────────────────────────
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns]
    for _, row in mat.iterrows():
        prefix = f"{row['chr_display']}_{row['arm']}_{row['hap']}"
        tlen = row.get("telo_length_kbp", 0)
        _add("per_telomere", ds, f"{prefix}_telo_length_kbp", f"{tlen:.3f}")
        for cat in cats_ds:
            v = row.get(cat, 0)
            if v > 0:
                _add("per_telomere", ds, f"{prefix}_{cat}_kbp", f"{v:.3f}")
                _add("per_telomere", ds, f"{prefix}_{cat}_frac", f"{v / tlen:.4f}" if tlen else "NA")

# ── 4b  Global category summary ────────────────────────────────────────────
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns]
    total_kbp = mat[cats_ds].sum().sum()
    for cat in cats_ds:
        s = mat[cat].sum()
        _add("global_summary", ds, f"{cat}_total_kbp", f"{s:.3f}")
        _add("global_summary", ds, f"{cat}_fraction", f"{s / total_kbp:.6f}" if total_kbp else "NA")

# ── 4c  p-arm vs q-arm comparison (Mann-Whitney U) ────────────────────────
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    for cat in cats_ds:
        p_vals = mat.loc[mat["arm"] == "p", cat].values
        q_vals = mat.loc[mat["arm"] == "q", cat].values
        if len(p_vals) > 1 and len(q_vals) > 1:
            u_stat, p_val = sp_stats.mannwhitneyu(p_vals, q_vals, alternative="two-sided")
            _add("arm_comparison", ds, f"{cat}_U", f"{u_stat:.1f}")
            _add("arm_comparison", ds, f"{cat}_pvalue", f"{p_val:.4e}")
            _add("arm_comparison", ds, f"{cat}_p_median_kbp", f"{np.median(p_vals):.3f}")
            _add("arm_comparison", ds, f"{cat}_q_median_kbp", f"{np.median(q_vals):.3f}")

# ── 4d  Acrocentric vs non-acrocentric (Mann-Whitney U) ───────────────────
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    mat_c = mat.copy()
    mat_c["is_acro"] = mat_c["chr_display"].isin(ACROCENTRIC)
    for cat in cats_ds:
        acro = mat_c.loc[mat_c["is_acro"], cat].values
        non  = mat_c.loc[~mat_c["is_acro"], cat].values
        if len(acro) > 1 and len(non) > 1:
            u_stat, p_val = sp_stats.mannwhitneyu(acro, non, alternative="two-sided")
            _add("acro_comparison", ds, f"{cat}_U", f"{u_stat:.1f}")
            _add("acro_comparison", ds, f"{cat}_pvalue", f"{p_val:.4e}")
            _add("acro_comparison", ds, f"{cat}_acro_median_kbp", f"{np.median(acro):.3f}")
            _add("acro_comparison", ds, f"{cat}_non_median_kbp",  f"{np.median(non):.3f}")

# ── 4e  Hypergeometric enrichment — is each error category enriched in
#         telomeric ends relative to the full 25 kbp NucFlag window?  ───────
for ds in nf_data:
    nf_raw = nf_data[ds]
    mat    = nf_matrix[ds]
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]

    # Full NucFlag 25kbp data (all categories)
    nf_full_path = NUCFLAG_INTERSECTION[ds].replace("nucflag_teloscope_", "nucflag_telo_").replace(".bed", ".tsv")
    if os.path.isfile(nf_full_path):
        nf_full = pd.read_csv(nf_full_path, sep="\t", header=None,
                              names=["chr", "start", "end", "category"])
        nf_full["length"] = nf_full["end"] - nf_full["start"]
        total_nf_bp = int(nf_full["length"].sum())

        # Total telomere bp
        total_telo_bp = int(nf_raw["overlap_bp"].sum())

        for cat in cats_ds:
            # Category bp in full NucFlag windows (population successes)
            K = int(nf_full.loc[nf_full["category"] == cat, "length"].sum())
            # Category bp in telomere intersections (observed successes)
            k = int(nf_raw.loc[nf_raw["nf_category"] == cat, "overlap_bp"].sum())
            N = total_nf_bp      # population size
            n = total_telo_bp    # number of draws

            if N > 0 and K > 0 and n > 0:
                # P(X >= k) under hypergeometric
                pval = sp_stats.hypergeom.sf(k - 1, N, K, n)
                expected = n * K / N
                fold = k / expected if expected > 0 else float("inf")
                _add("hypergeom_enrichment", ds, f"{cat}_k", str(k))
                _add("hypergeom_enrichment", ds, f"{cat}_K", str(K))
                _add("hypergeom_enrichment", ds, f"{cat}_N", str(N))
                _add("hypergeom_enrichment", ds, f"{cat}_n", str(n))
                _add("hypergeom_enrichment", ds, f"{cat}_expected", f"{expected:.1f}")
                _add("hypergeom_enrichment", ds, f"{cat}_fold_enrichment", f"{fold:.3f}")
                _add("hypergeom_enrichment", ds, f"{cat}_pvalue", f"{pval:.4e}")
    else:
        _add("hypergeom_enrichment", ds, "status", f"full NucFlag file not found: {nf_full_path}")

# ── 4f  ONT vs HiFi category proportions (Fisher exact 2x2) ───────────────
if len(nf_matrix) == 2:
    mats = {ds: m for ds, m in nf_matrix.items()}
    all_cats_both = sorted(
        set(c for m in mats.values() for c in m.columns if c in all_present_cats),
        key=lambda c: all_present_cats.index(c) if c in all_present_cats else 999,
    )
    for cat in all_cats_both:
        vals = {}
        totals = {}
        for ds, m in mats.items():
            vals[ds] = m[cat].sum() * 1000 if cat in m.columns else 0  # back to bp for counts
            cats_ds = [c for c in all_present_cats if c in m.columns]
            totals[ds] = m[cats_ds].sum().sum() * 1000

        a = int(round(vals.get("HiFi", 0)))
        b = int(round(totals.get("HiFi", 0) - vals.get("HiFi", 0)))
        c = int(round(vals.get("ONT", 0)))
        d = int(round(totals.get("ONT", 0) - vals.get("ONT", 0)))

        if min(a + c, a + b) > 0:
            odds, pval = sp_stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
            _add("ont_vs_hifi", "both", f"{cat}_odds_ratio", f"{odds:.4f}")
            _add("ont_vs_hifi", "both", f"{cat}_fisher_pvalue", f"{pval:.4e}")
            _add("ont_vs_hifi", "both", f"{cat}_hifi_frac", f"{a / (a + b):.6f}" if (a + b) else "NA")
            _add("ont_vs_hifi", "both", f"{cat}_ont_frac",  f"{c / (c + d):.6f}" if (c + d) else "NA")

# ── 4g  Per-chromosome type grouping (metacentric / submetacentric / acrocentric)
CHR_TYPE = {}
for c in chroms:
    CHR_TYPE[c] = _classify_chr(c)

for ds, mat in nf_matrix.items():
    mat_c = mat.copy()
    mat_c["chr_type"] = mat_c["chr_display"].map(CHR_TYPE)
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    for ctype in sorted(mat_c["chr_type"].unique()):
        sub = mat_c.loc[mat_c["chr_type"] == ctype]
        total = sub[cats_ds].sum().sum()
        for cat in cats_ds:
            s = sub[cat].sum()
            _add("chr_type_breakdown", ds, f"{ctype}_{cat}_kbp", f"{s:.3f}")
            _add("chr_type_breakdown", ds, f"{ctype}_{cat}_frac", f"{s / total:.6f}" if total else "NA")

# ── 4h  Hap1 vs hap2 per-category comparison (Wilcoxon signed-rank) ───────
for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    for cat in cats_ds:
        h1 = mat.loc[mat["hap"] == "hap1"].set_index(["chr_display", "arm"])[cat]
        h2 = mat.loc[mat["hap"] == "hap2"].set_index(["chr_display", "arm"])[cat]
        common = h1.index.intersection(h2.index)
        if len(common) >= 5:
            v1 = h1.loc[common].values
            v2 = h2.loc[common].values
            diff = v1 - v2
            if np.any(diff != 0):
                stat, pval = sp_stats.wilcoxon(v1, v2, alternative="two-sided")
                _add("hap_comparison", ds, f"{cat}_wilcoxon_stat", f"{stat:.1f}")
                _add("hap_comparison", ds, f"{cat}_wilcoxon_pvalue", f"{pval:.4e}")
                _add("hap_comparison", ds, f"{cat}_hap1_median", f"{np.median(v1):.3f}")
                _add("hap_comparison", ds, f"{cat}_hap2_median", f"{np.median(v2):.3f}")

# ── Write stats TSV ────────────────────────────────────────────────────────
df_stats = pd.DataFrame(stat_rows)
df_stats.to_csv(stats_out, sep="\t", index=False)
print(f"[OK] {stats_out}  ({len(df_stats)} rows)")

print("\n[DONE] H9 telomere NucFlag plots & statistics generated.")
