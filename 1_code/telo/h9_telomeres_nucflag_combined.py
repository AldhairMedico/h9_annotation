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

plots_dir = os.path.join(figures_dir, "26.03.01_telomeres_nucflag")
os.makedirs(plots_dir, exist_ok=True)

nucflag_dir = os.path.join(data_dir, "nucflag")

teloscope_dir = os.path.join(data_dir, "25.12.10_teloscope_asm")
terminal_telomeres = {
    "hap1": os.path.join(teloscope_dir, "H9_T2T_v0.1_hap1.fasta_terminal_telomeres.bed"),
    "hap2": os.path.join(teloscope_dir, "H9_T2T_v0.1_hap2.fasta_terminal_telomeres.bed"),
}

NUCFLAG_INTERSECTION = {
    "HiFi": os.path.join(nucflag_dir, "nucflag_teloscope_hifi.bed"),
    "ONT":  os.path.join(nucflag_dir, "nucflag_teloscope_ont.bed"),
}

stats_out = os.path.join(plots_dir, f"h9_nucflag_telomere_stats_{run_label}.tsv")

# ═══════════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════════
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

# Cell journal aesthetics (Arial / DejaVu Sans fallback, minimal spines, white bg)
_font_family = "Arial"
try:
    from matplotlib.font_manager import findfont, FontProperties
    if findfont(FontProperties(family="Arial")) == findfont(FontProperties()):
        _font_family = "DejaVu Sans"  # Arial not installed; use fallback
except Exception:
    _font_family = "DejaVu Sans"

plt.rcParams.update({
    "font.family":       _font_family,
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
# 1. Load telomere data  (per-hap teloscope BEDs, 10 columns, no assembly)
# ═══════════════════════════════════════════════════════════════════════════════
rows: list = []
for hap_label, telo_path in terminal_telomeres.items():
    with open(telo_path, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) == 10:
                rows.append(fields)

df_telo = pd.DataFrame(rows, columns=[
    "chr", "start", "end", "length", "label",
    "fwdCounts", "revCounts", "canCounts", "nonCanCounts",
    "chr_length",
])
for col in ("start", "end", "length", "chr_length"):
    df_telo[col] = pd.to_numeric(df_telo[col])

df_telo["chr_display"] = df_telo["chr"].str.replace(r"_hap[12]$", "", regex=True)
df_telo["hap"]         = df_telo["chr"].str.extract(r"_(hap[12])$")[0]
df_telo = df_telo.loc[df_telo["hap"].notna()].copy()
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
    """Parse bedtools intersect -wo output (17 columns).

    Teloscope BED (5 cols): chr, start, end, length, arm
    NucFlag BED9+2 (11 cols): chrom, chromStart, chromEnd, name,
        score, strand, thickStart, thickEnd, itemRgb, zscore, af
    Overlap (1 col): overlap_bp
    """
    col_names = [
        "telo_chr", "telo_start", "telo_end", "telo_length", "arm",
        "nf_chr", "nf_start", "nf_end", "nf_category",
        "nf_score", "nf_strand", "nf_thickStart", "nf_thickEnd",
        "nf_itemRgb", "nf_zscore", "nf_af",
        "overlap_bp",
    ]
    df = pd.read_csv(path, sep="\t", header=None, names=col_names)
    for c in ("telo_start", "telo_end", "telo_length", "nf_start", "nf_end", "overlap_bp"):
        df[c] = pd.to_numeric(df[c])
    df["chr_display"] = df["telo_chr"].str.replace(r"_hap[12]$", "", regex=True)
    df["hap"] = df["telo_chr"].str.extract(r"_(hap[12])$")[0]
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
# 4. Comprehensive statistics  →  separate TSV files + console summaries
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STATISTICS — H9 Telomere × NucFlag Analysis")
print("=" * 80)

stats_dir = os.path.join(plots_dir, f"stats_{run_label}")
os.makedirs(stats_dir, exist_ok=True)


def _write_section(filename: str, df: pd.DataFrame, description: str) -> None:
    """Write a stats DataFrame to TSV and print a descriptive console summary."""
    path = os.path.join(stats_dir, filename)
    df.to_csv(path, sep="\t", index=False)
    print(f"\n[SAVED] {path}  ({len(df)} rows)")
    print(f"        {description}")


# ── 4a  Per chr × arm × hap breakdown ──────────────────────────────────────
print("\n" + "─" * 80)
print("4a. PER-TELOMERE BREAKDOWN")
print("    NucFlag category composition for every individual telomere")
print("    (chr × arm × hap), with absolute length (kbp) and fraction.")
print("─" * 80)

for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns]
    telo_rows: list = []
    for _, row in mat.iterrows():
        tlen = row.get("telo_length_kbp", 0)
        entry = {
            "chromosome": row["chr_display"],
            "arm": row["arm"],
            "haplotype": row["hap"],
            "telomere_length_kbp": f"{tlen:.3f}",
        }
        for cat in cats_ds:
            v = row.get(cat, 0)
            entry[f"{cat}_kbp"] = f"{v:.3f}"
            entry[f"{cat}_frac"] = f"{v / tlen:.4f}" if tlen else "NA"
        telo_rows.append(entry)

    df_telo_stats = pd.DataFrame(telo_rows)
    _write_section(
        f"4a_per_telomere_{ds}.tsv", df_telo_stats,
        f"Per-telomere NucFlag breakdown for {ds}: "
        f"{len(df_telo_stats)} telomeres across {mat['chr_display'].nunique()} chromosomes, "
        f"{len(cats_ds)} NucFlag categories.",
    )

    # Console summary: total telomere length per haplotype
    for hap in HAP_ORDER:
        sub = mat.loc[mat["hap"] == hap]
        total = sub["telo_length_kbp"].sum() if "telo_length_kbp" in sub.columns else sub[cats_ds].sum().sum()
        n_telo = len(sub)
        print(f"    {ds} {hap}: {n_telo} telomeres, total length = {total:.1f} kbp")

# ── 4b  Global category summary ────────────────────────────────────────────
print("\n" + "─" * 80)
print("4b. GLOBAL CATEGORY SUMMARY")
print("    Total kbp and fraction of telomere length occupied by each")
print("    NucFlag category, aggregated across all chromosomes and haplotypes.")
print("─" * 80)

for ds, mat in nf_matrix.items():
    cats_ds = [c for c in all_present_cats if c in mat.columns]
    total_kbp = mat[cats_ds].sum().sum()
    summary_rows: list = []
    print(f"\n    {ds}  (total covered telomere length: {total_kbp:.1f} kbp)")
    print(f"    {'Category':<20s} {'Total (kbp)':>12s} {'Fraction':>10s}")
    print(f"    {'─' * 20} {'─' * 12} {'─' * 10}")
    for cat in cats_ds:
        s = mat[cat].sum()
        frac = s / total_kbp if total_kbp else 0
        summary_rows.append({
            "dataset": ds,
            "category": cat,
            "total_kbp": f"{s:.3f}",
            "fraction": f"{frac:.6f}" if total_kbp else "NA",
        })
        print(f"    {cat:<20s} {s:>12.3f} {frac:>10.4f}")

    df_global = pd.DataFrame(summary_rows)
    _write_section(
        f"4b_global_summary_{ds}.tsv", df_global,
        f"Global NucFlag category summary for {ds}: fraction of total "
        f"telomere length ({total_kbp:.1f} kbp) per category.",
    )

# ── 4e  Hypergeometric enrichment — is each error category enriched in
#         telomeric ends relative to the full 25 kbp NucFlag window?  ───────
print("\n" + "─" * 80)
print("4e. HYPERGEOMETRIC ENRICHMENT TEST")
print("    Tests whether each NucFlag category is enriched in telomeric")
print("    regions compared to the full terminal 25 kbp NucFlag window.")
print("    fold > 1 = enriched in telomeres; fold < 1 = depleted.")
print("─" * 80)

hyper_rows: list = []
for ds in nf_data:
    nf_raw = nf_data[ds]
    mat    = nf_matrix[ds]
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]

    # Full NucFlag 25kbp data (all categories)
    nf_full_path = NUCFLAG_INTERSECTION[ds].replace("nucflag_teloscope_", "nucflag_telo_").replace(".bed", ".tsv")
    if os.path.isfile(nf_full_path):
        nf_full = pd.read_csv(nf_full_path, sep="\t", comment="#",
                              header=None,
                              names=["chrom", "chromStart", "chromEnd", "name",
                                     "score", "strand", "thickStart", "thickEnd",
                                     "itemRgb", "zscore", "af"])
        nf_full["chromStart"] = pd.to_numeric(nf_full["chromStart"])
        nf_full["chromEnd"]   = pd.to_numeric(nf_full["chromEnd"])
        nf_full["length"] = nf_full["chromEnd"] - nf_full["chromStart"]
        nf_full["category"] = nf_full["name"]
        total_nf_bp = int(nf_full["length"].sum())

        # Total telomere bp
        total_telo_bp = int(nf_raw["overlap_bp"].sum())

        print(f"\n    {ds}  (NucFlag window: {total_nf_bp:,} bp, telomere overlap: {total_telo_bp:,} bp)")
        print(f"    {'Category':<20s} {'obs (bp)':>12s} {'exp (bp)':>12s} {'fold':>8s} {'p-value':>12s}")
        print(f"    {'─' * 20} {'─' * 12} {'─' * 12} {'─' * 8} {'─' * 12}")

        for cat in cats_ds:
            K = int(nf_full.loc[nf_full["category"] == cat, "length"].sum())
            k = int(nf_raw.loc[nf_raw["nf_category"] == cat, "overlap_bp"].sum())
            N = total_nf_bp
            n = total_telo_bp

            if N > 0 and K > 0 and n > 0:
                pval = sp_stats.hypergeom.sf(k - 1, N, K, n)
                expected = n * K / N
                fold = k / expected if expected > 0 else float("inf")
                sig = " *" if pval < 0.05 else ""
                hyper_rows.append({
                    "dataset": ds, "category": cat,
                    "observed_bp": str(k),
                    "expected_bp": f"{expected:.1f}",
                    "population_category_bp_K": str(K),
                    "population_total_bp_N": str(N),
                    "telomere_total_bp_n": str(n),
                    "fold_enrichment": f"{fold:.3f}",
                    "p_value": f"{pval:.4e}",
                    "significant_0.05": "yes" if pval < 0.05 else "no",
                })
                print(f"    {cat:<20s} {k:>12,d} {expected:>12,.1f} {fold:>8.2f} {pval:>12.4e}{sig}")
    else:
        print(f"\n    {ds}: full NucFlag file not found: {nf_full_path}")

if hyper_rows:
    _write_section(
        "4e_hypergeometric_enrichment.tsv", pd.DataFrame(hyper_rows),
        "Hypergeometric enrichment test: is each NucFlag category over-/under-represented "
        "in telomeric ends vs the full 25 kbp NucFlag terminal window?",
    )

# ── 4f  ONT vs HiFi category proportions (Fisher exact 2x2) ───────────────
print("\n" + "─" * 80)
print("4f. ONT vs HiFi CATEGORY PROPORTIONS  (Fisher exact test)")
print("    Compares the fraction of each NucFlag category between ONT and")
print("    HiFi datasets using a 2×2 contingency table (Fisher exact test).")
print("─" * 80)

fisher_rows: list = []
if len(nf_matrix) == 2:
    mats = {ds: m for ds, m in nf_matrix.items()}
    all_cats_both = sorted(
        set(c for m in mats.values() for c in m.columns if c in all_present_cats),
        key=lambda c: all_present_cats.index(c) if c in all_present_cats else 999,
    )
    print(f"\n    {'Category':<20s} {'HiFi frac':>10s} {'ONT frac':>10s} {'OR':>8s} {'p-value':>12s}")
    print(f"    {'─' * 20} {'─' * 10} {'─' * 10} {'─' * 8} {'─' * 12}")
    for cat in all_cats_both:
        vals = {}
        totals = {}
        for ds, m in mats.items():
            vals[ds] = m[cat].sum() * 1000 if cat in m.columns else 0
            cats_ds = [c for c in all_present_cats if c in m.columns]
            totals[ds] = m[cats_ds].sum().sum() * 1000

        a = int(round(vals.get("HiFi", 0)))
        b = int(round(totals.get("HiFi", 0) - vals.get("HiFi", 0)))
        c = int(round(vals.get("ONT", 0)))
        d = int(round(totals.get("ONT", 0) - vals.get("ONT", 0)))

        if min(a + c, a + b) > 0:
            odds, pval = sp_stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
            hifi_f = a / (a + b) if (a + b) else 0
            ont_f  = c / (c + d) if (c + d) else 0
            sig = " *" if pval < 0.05 else ""
            fisher_rows.append({
                "category": cat,
                "hifi_bp": str(a), "hifi_other_bp": str(b),
                "ont_bp": str(c), "ont_other_bp": str(d),
                "hifi_fraction": f"{hifi_f:.6f}" if (a + b) else "NA",
                "ont_fraction": f"{ont_f:.6f}" if (c + d) else "NA",
                "odds_ratio": f"{odds:.4f}",
                "p_value": f"{pval:.4e}",
                "significant_0.05": "yes" if pval < 0.05 else "no",
            })
            print(f"    {cat:<20s} {hifi_f:>10.4f} {ont_f:>10.4f} {odds:>8.3f} {pval:>12.4e}{sig}")

if fisher_rows:
    _write_section(
        "4f_ont_vs_hifi_fisher.tsv", pd.DataFrame(fisher_rows),
        "Fisher exact test comparing NucFlag category proportions between "
        "ONT and HiFi datasets (2×2 contingency tables, bp-level).",
    )

# ── 4g  Per-chromosome type grouping (metacentric / submetacentric / acrocentric)
print("\n" + "─" * 80)
print("4g. CHROMOSOME-TYPE BREAKDOWN")
print("    NucFlag category composition grouped by chromosome morphology:")
print("    metacentric, submetacentric, acrocentric.")
print("─" * 80)

CHR_TYPE = {}
for c in chroms:
    CHR_TYPE[c] = _classify_chr(c)

chrtype_rows: list = []
for ds, mat in nf_matrix.items():
    mat_c = mat.copy()
    mat_c["chr_type"] = mat_c["chr_display"].map(CHR_TYPE)
    cats_ds = [c for c in all_present_cats if c in mat.columns and mat[c].sum() > 0]
    print(f"\n    {ds}")
    for ctype in sorted(mat_c["chr_type"].unique()):
        sub = mat_c.loc[mat_c["chr_type"] == ctype]
        total = sub[cats_ds].sum().sum()
        n_chr_type = sub["chr_display"].nunique()
        print(f"      {ctype} ({n_chr_type} chr, total: {total:.1f} kbp):")
        for cat in cats_ds:
            s = sub[cat].sum()
            frac = s / total if total else 0
            chrtype_rows.append({
                "dataset": ds, "chr_type": ctype,
                "category": cat,
                "total_kbp": f"{s:.3f}",
                "fraction": f"{frac:.6f}" if total else "NA",
            })
            if s > 0:
                print(f"        {cat:<20s}  {s:>8.3f} kbp  ({frac:.4f})")

if chrtype_rows:
    _write_section(
        "4g_chrtype_breakdown.tsv", pd.DataFrame(chrtype_rows),
        "NucFlag category breakdown by chromosome type (metacentric / "
        "submetacentric / acrocentric), with kbp and fraction.",
    )

# ── Final combined stats TSV (backward-compatible) ─────────────────────────
all_section_dfs = []
for fname in sorted(os.listdir(stats_dir)):
    if fname.endswith(".tsv"):
        fpath = os.path.join(stats_dir, fname)
        section_name = fname.replace(".tsv", "")
        df_sec = pd.read_csv(fpath, sep="\t")
        df_sec.insert(0, "section", section_name)
        all_section_dfs.append(df_sec)

if all_section_dfs:
    df_combined = pd.concat(all_section_dfs, ignore_index=True)
    df_combined.to_csv(stats_out, sep="\t", index=False)
    print(f"\n[SAVED] Combined stats: {stats_out}  ({len(df_combined)} rows)")

print("\n" + "=" * 80)
print("[DONE] H9 telomere NucFlag plots & statistics generated.")
print(f"       Figures: {plots_dir}")
print(f"       Stats:   {stats_dir}")
print("=" * 80)
