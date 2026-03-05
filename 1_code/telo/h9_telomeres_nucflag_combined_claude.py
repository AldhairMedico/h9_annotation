#!/usr/bin/env python3
"""
H9 terminal telomere quality assessment: combined HiFi + ONT NucFlag analysis.

Uses the base-pair-resolution combined NucFlag TSV (hifi_name × ont_name per
interval) intersected with Teloscope terminal telomere annotations to assess
assembly quality at telomeric ends using orthogonal sequencing technologies.

Generates:
  Plot 1 – Per-chromosome stacked bars (mirrored p/q), one panel per tech
           plus a combined-quality panel.
  Plot 2 – Concordance confusion matrix: HiFi × ONT categories (bp heatmap).
  Plot 3 – Combined quality summary per haplotype (both_correct / rescued /
           both_error).
  Plot 4 – a) NucFlag composition heatmap (fraction per chr × category),
           b) ONT-vs-HiFi proportion scatter,
           c) Technology rescue per telomere,
           d) Telomere quality ranking.
  Stats  – Comprehensive TSVs: per-telomere breakdown, concordance, rescue,
           enrichment, quality ranking.  All values in bp + fraction.

Input files:
  - nucflag_telo_combined.tsv
  - H9_T2T_v0.1_dip.fasta_terminal_telomeres.bed
  - nucflag_palette.tsv
"""

import os
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from typing import Dict

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcolors

# ═══════════════════════════════════════════════════════════════════════════════
# 0. Paths
# ═══════════════════════════════════════════════════════════════════════════════
script_dir  = os.path.dirname(os.path.abspath(__file__))
base_dir    = os.path.dirname(os.path.dirname(script_dir))
data_dir    = os.path.join(base_dir, "2_data", "2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures", "3.1_draft")
run_label   = "v2"

plots_dir = os.path.join(figures_dir, "26.03.01_telomeres_nucflag_combined_claude")
os.makedirs(plots_dir, exist_ok=True)
stats_dir = os.path.join(plots_dir, f"stats_{run_label}")
os.makedirs(stats_dir, exist_ok=True)

nucflag_dir   = os.path.join(data_dir, "nucflag")
teloscope_dir = os.path.join(data_dir, "25.12.10_teloscope_compiled")

COMBINED_TSV   = os.path.join(nucflag_dir, "nucflag_telo_combined.tsv")
PALETTE_TSV    = os.path.join(nucflag_dir, "nucflag_palette.tsv")
TERMINAL_BED   = os.path.join(teloscope_dir,
                              "H9_T2T_v0.1_dip.fasta_terminal_telomeres.bed")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Constants & palette
# ═══════════════════════════════════════════════════════════════════════════════
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
ACROCENTRIC    = {"chr13", "chr14", "chr15", "chr21", "chr22"}
METACENTRIC    = {"chr1", "chr3", "chr16", "chr19", "chr20"}
SUBMETACENTRIC = set(CHR_ORDER) - ACROCENTRIC - METACENTRIC

HAP_ORDER = ["hap1", "hap2"]
EXTS = ("png", "pdf", "svg")

def _chr_key(c: str) -> int:
    return CHR_ORDER.index(c) if c in CHR_ORDER else 100

def _classify_chr(c: str) -> str:
    if c in ACROCENTRIC:    return "acrocentric"
    if c in METACENTRIC:    return "metacentric"
    return "submetacentric"

# Load NucFlag palette from TSV → {name: "#RRGGBB"}
_pal = pd.read_csv(PALETTE_TSV, sep="\t")
NUCFLAG_COLORS: Dict[str, str] = {}
for _, r in _pal.iterrows():
    rgb = r["itemRgb"].strip('"').split(",")
    NUCFLAG_COLORS[r["name"]] = "#{:02x}{:02x}{:02x}".format(
        int(rgb[0]), int(rgb[1]), int(rgb[2])
    )
NUCFLAG_COLORS["uncovered"]    = "#D4D4D4"
NUCFLAG_COLORS["dinucleotide"] = NUCFLAG_COLORS.get("dinucleotide",
                                                     NUCFLAG_COLORS.get("homopolymer", "#ECEC00"))
NUCFLAG_COLORS["other_repeat"] = NUCFLAG_COLORS.get("other_repeat",
                                                     NUCFLAG_COLORS.get("simple_repeat", "#FF0080"))
NUCFLAG_COLORS["scaffold"]     = NUCFLAG_COLORS.get("scaffold", "#999999")

# Colours for the 4-tier combined quality classification
QUALITY_COLORS = {
    "both_correct":      "#2D8E47",  # deep green — confident good
    "hifi_only_correct": "#6BAED6",  # steel blue — HiFi validates
    "ont_only_correct":  "#FDAE6B",  # peach — ONT validates
    "both_error":        "#CB3335",  # crimson — confident problem
}
QUALITY_ORDER = ["both_correct", "hifi_only_correct", "ont_only_correct",
                 "both_error"]
QUALITY_LABELS = {
    "both_correct":      "Both correct",
    "hifi_only_correct": "HiFi correct only",
    "ont_only_correct":  "ONT correct only",
    "both_error":        "Both error",
}

CHRTYPE_COLORS = {
    "metacentric":    "#3C7DC4",
    "submetacentric": "#E69F00",
    "acrocentric":    "#CB3335",
}

# ── matplotlib style ────────────────────────────────────────────────────────
_font_family = "Arial"
try:
    from matplotlib.font_manager import findfont, FontProperties
    if findfont(FontProperties(family="Arial")) == findfont(FontProperties()):
        _font_family = "DejaVu Sans"
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
TICK  = 6
LABEL = 7

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Helpers
# ═══════════════════════════════════════════════════════════════════════════════
def _save(fig, stem: str) -> None:
    for ext in EXTS:
        fig.savefig(
            os.path.join(plots_dir, f"{stem}.{ext}"),
            dpi=600 if ext == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(f"  [OK] {stem}")

def _box(ax):
    for side in ("left", "bottom"):
        ax.spines[side].set_edgecolor("black")
        ax.spines[side].set_linewidth(0.4)
        ax.spines[side].set_visible(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

def _write_section(filename: str, df: pd.DataFrame, description: str) -> None:
    path = os.path.join(stats_dir, filename)
    df.to_csv(path, sep="\t", index=False)
    print(f"  [SAVED] {path}  ({len(df)} rows)")
    print(f"          {description}")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Load data
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 80)
print("Loading data...")
print("=" * 80)

# 3a. Combined NucFlag (full terminal windows, both haplotypes)
df_nf = pd.read_csv(COMBINED_TSV, sep="\t")
for c in ("chromStart", "chromEnd", "length"):
    df_nf[c] = pd.to_numeric(df_nf[c])
print(f"  Combined NucFlag: {len(df_nf)} intervals")

# 3b. Teloscope terminal telomeres (diploid BED, 10 columns)
telo_cols = ["chr", "start", "end", "length", "arm",
             "fwdCounts", "revCounts", "canCounts", "nonCanCounts",
             "chr_length"]
df_telo = pd.read_csv(TERMINAL_BED, sep="\t", header=None, names=telo_cols)
for c in ("start", "end", "length", "chr_length"):
    df_telo[c] = pd.to_numeric(df_telo[c])
df_telo["chr_display"] = df_telo["chr"].str.replace(r"_hap[12]$", "", regex=True)
df_telo["hap"] = df_telo["chr"].str.extract(r"_(hap[12])$")[0]
df_telo = df_telo.loc[df_telo["hap"].notna()].copy()
print(f"  Teloscope telomeres: {len(df_telo)} entries "
      f"({df_telo['chr_display'].nunique()} chromosomes)")

chroms    = sorted(df_telo["chr_display"].unique(), key=_chr_key)
n_chr     = len(chroms)
chr_short = [c.replace("chr", "") for c in chroms]

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Intersect combined NucFlag with telomere regions
# ═══════════════════════════════════════════════════════════════════════════════
print("\nIntersecting NucFlag with telomere regions...")

rows_intersected: list = []
for _, telo in df_telo.iterrows():
    chrom     = telo["chr"]
    t_start   = telo["start"]
    t_end     = telo["end"]
    arm       = telo["arm"]
    hap       = telo["hap"]
    chr_disp  = telo["chr_display"]
    telo_len  = telo["length"]

    # Get NucFlag intervals for this chromosome
    nf_chr = df_nf.loc[df_nf["chrom"] == chrom]

    covered_bp = 0
    for _, nf in nf_chr.iterrows():
        ov_start = max(t_start, nf["chromStart"])
        ov_end   = min(t_end, nf["chromEnd"])
        if ov_start >= ov_end:
            continue
        ov_bp = ov_end - ov_start
        # Skip the giant "uncovered" middle region (sanity check)
        if nf["hifi_name"] == "uncovered" and nf["ont_name"] == "uncovered":
            continue
        covered_bp += ov_bp
        rows_intersected.append({
            "chr_display": chr_disp,
            "arm":         arm,
            "hap":         hap,
            "chrom":       chrom,
            "ov_start":    int(ov_start),
            "ov_end":      int(ov_end),
            "overlap_bp":  int(ov_bp),
            "hifi_name":   nf["hifi_name"],
            "ont_name":    nf["ont_name"],
            "hifi_score":  nf["hifi_score"],
            "ont_score":   nf["ont_score"],
            "concordant":  nf["concordant"],
            "telo_length": int(telo_len),
        })

    # Add uncovered telomere bp (if NucFlag doesn't cover the full telomere)
    uncov = telo_len - covered_bp
    if uncov > 0:
        rows_intersected.append({
            "chr_display": chr_disp,
            "arm":         arm,
            "hap":         hap,
            "chrom":       chrom,
            "ov_start":    -1,
            "ov_end":      -1,
            "overlap_bp":  int(uncov),
            "hifi_name":   "uncovered",
            "ont_name":    "uncovered",
            "hifi_score":  0,
            "ont_score":   0,
            "concordant":  True,
            "telo_length": int(telo_len),
        })

df_ix = pd.DataFrame(rows_intersected)
print(f"  Intersection: {len(df_ix)} intervals, "
      f"{df_ix['overlap_bp'].sum():,} bp total telomere coverage")

# ── Derive quality tier for every interval ──────────────────────────────────
def _quality_tier(hifi, ont):
    if hifi == "correct" and ont == "correct":
        return "both_correct"
    if hifi == "correct":
        return "hifi_only_correct"
    if ont == "correct":
        return "ont_only_correct"
    return "both_error"

df_ix["quality"] = df_ix.apply(
    lambda r: _quality_tier(r["hifi_name"], r["ont_name"]), axis=1
)

# ── Build per-telomere aggregation matrices ─────────────────────────────────
# For HiFi categories
hifi_cats_present = sorted(
    [c for c in df_ix["hifi_name"].unique() if c != "uncovered"],
    key=lambda c: list(NUCFLAG_COLORS.keys()).index(c)
         if c in NUCFLAG_COLORS else 999
)
# For ONT categories
ont_cats_present = sorted(
    [c for c in df_ix["ont_name"].unique() if c != "uncovered"],
    key=lambda c: list(NUCFLAG_COLORS.keys()).index(c)
         if c in NUCFLAG_COLORS else 999
)
# Union of all NucFlag categories present
all_cats = sorted(
    set(hifi_cats_present) | set(ont_cats_present),
    key=lambda c: list(NUCFLAG_COLORS.keys()).index(c)
         if c in NUCFLAG_COLORS else 999
)
if "uncovered" not in all_cats:
    all_cats.append("uncovered")

print(f"  HiFi categories: {hifi_cats_present}")
print(f"  ONT categories:  {ont_cats_present}")

# Aggregation: per telomere × tech category
def _build_tech_matrix(df, tech_col, cats):
    """Pivot to one row per (chr_display, arm, hap) with bp per category."""
    agg = (df.groupby(["chr_display", "arm", "hap", tech_col])["overlap_bp"]
             .sum().reset_index())
    piv = agg.pivot_table(
        index=["chr_display", "arm", "hap"],
        columns=tech_col, values="overlap_bp", fill_value=0,
    ).reset_index()
    # Ensure all categories present
    for cat in cats:
        if cat not in piv.columns:
            piv[cat] = 0
    # Add uncovered if not present
    if "uncovered" not in piv.columns:
        piv["uncovered"] = 0
    # Telomere length
    tlen = (df.groupby(["chr_display", "arm", "hap"])["telo_length"]
              .first().reset_index())
    piv = piv.merge(tlen, on=["chr_display", "arm", "hap"], how="left")
    return piv

mat_hifi = _build_tech_matrix(df_ix, "hifi_name", hifi_cats_present)
mat_ont  = _build_tech_matrix(df_ix, "ont_name", ont_cats_present)

# Aggregation: per telomere × quality tier
mat_qual = (df_ix.groupby(["chr_display", "arm", "hap", "quality"])["overlap_bp"]
                 .sum().reset_index())
mat_qual = mat_qual.pivot_table(
    index=["chr_display", "arm", "hap"],
    columns="quality", values="overlap_bp", fill_value=0,
).reset_index()
for q in QUALITY_ORDER:
    if q not in mat_qual.columns:
        mat_qual[q] = 0
tlen = (df_ix.groupby(["chr_display", "arm", "hap"])["telo_length"]
             .first().reset_index())
mat_qual = mat_qual.merge(tlen, on=["chr_display", "arm", "hap"], how="left")

total_telo_bp = int(df_telo["length"].sum())
print(f"\n  Total telomere length (Teloscope): {total_telo_bp:,} bp")
for q in QUALITY_ORDER:
    bp = int(mat_qual[q].sum())
    frac = bp / total_telo_bp
    print(f"    {QUALITY_LABELS[q]:<25s} {bp:>10,d} bp  ({frac:.4f})")


# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
#                              P L O T S
# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("Generating plots...")
print("=" * 80)

bw      = 0.30
offsets = {"hap1": -bw / 2 - 0.02, "hap2": bw / 2 + 0.02}
pad     = 0.5


def _get_val(mat, hap, chrom, arm, cat):
    m = ((mat["hap"] == hap) & (mat["chr_display"] == chrom)
         & (mat["arm"] == arm))
    rows = mat.loc[m]
    if rows.empty or cat not in rows.columns:
        return 0.0
    return float(rows[cat].iloc[0])


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Per-chromosome stacked bars (mirrored p/q)
#          Three variants: HiFi, ONT, Combined quality
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 1: per-chromosome stacked bars ──")


def plot1_horizontal(mat, cats, color_map, ds_label, legend_labels=None):
    """Horizontal mirrored: ← p-arm | chr labels | q-arm →"""
    fig, (ax_p, ax_q) = plt.subplots(
        1, 2, figsize=(3.5, 4.5), sharey=True, constrained_layout=True,
    )
    yy = np.arange(n_chr)

    # Compute max for axis scaling (in kbp)
    raw_max = 0
    for arm_label in ("p", "q"):
        for hap in HAP_ORDER:
            for ci in range(n_chr):
                total = sum(_get_val(mat, hap, chroms[ci], arm_label, c)
                            for c in cats)
                raw_max = max(raw_max, total)
    x_max = float(np.ceil(raw_max / 1000 / 5) * 5)  # kbp, rounded to 5

    for arm_label, ax, invert in [("p", ax_p, True), ("q", ax_q, False)]:
        for hap in HAP_ORDER:
            bottoms = np.zeros(n_chr)
            for cat in cats:
                vals = np.array([_get_val(mat, hap, c, arm_label, cat)
                                 for c in chroms]) / 1000.0  # → kbp
                ax.barh(
                    yy + offsets[hap], vals, left=bottoms, height=bw,
                    color=color_map.get(cat, "#D4D4D4"),
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

    labels = legend_labels or {c: c for c in cats}
    handles = [Patch(fc=color_map.get(c, "#D4D4D4"),
                     label=labels.get(c, c), linewidth=0)
               for c in cats]
    handles += [
        Patch(fc="white", ec="black", lw=0.4, label="hap1 (left bar)"),
        Patch(fc="white", ec="black", lw=0.4, label="hap2 (right bar)"),
    ]
    fig.legend(handles=handles, fontsize=TICK - 1,
               loc="lower center", bbox_to_anchor=(0.5, -0.08),
               ncol=min(len(handles), 4), frameon=False,
               handlelength=0.8, handleheight=0.6,
               handletextpad=0.4, borderpad=0.2, labelspacing=0.25)
    _save(fig, f"fig1_{ds_label}_horizontal_{run_label}")


def plot1_vertical(mat, cats, color_map, ds_label, legend_labels=None):
    """Vertical mirrored: ↑ p-arm | chr labels | q-arm ↓"""
    fig, (ax_p, ax_q) = plt.subplots(
        2, 1, figsize=(5, 3), sharex=True, constrained_layout=True,
    )
    xx = np.arange(n_chr)

    raw_max = 0
    for arm_label in ("p", "q"):
        for hap in HAP_ORDER:
            for ci in range(n_chr):
                total = sum(_get_val(mat, hap, chroms[ci], arm_label, c)
                            for c in cats)
                raw_max = max(raw_max, total)
    y_max = float(np.ceil(raw_max / 1000 / 5) * 5)

    for arm_label, ax, invert in [("p", ax_p, False), ("q", ax_q, True)]:
        for hap in HAP_ORDER:
            bottoms = np.zeros(n_chr)
            for cat in cats:
                vals = np.array([_get_val(mat, hap, c, arm_label, cat)
                                 for c in chroms]) / 1000.0
                ax.bar(
                    xx + offsets[hap], vals, bottom=bottoms, width=bw,
                    color=color_map.get(cat, "#D4D4D4"),
                    edgecolor="none", alpha=0.9,
                )
                bottoms += vals
        if invert:
            ax.set_ylim(y_max, 0)
        else:
            ax.set_ylim(0, y_max)
        ax.set_xlim(-pad, n_chr - 1 + pad)
        ax.set_ylabel("Telomere length (kbp)", fontsize=LABEL)
        ax.set_title(f"{arm_label}-arm", fontsize=TICK)
        ax.tick_params(axis="y", labelsize=TICK)
        _box(ax)

    ax_q.set_xticks(xx)
    ax_q.set_xticklabels(chr_short, fontsize=TICK, rotation=0)

    labels = legend_labels or {c: c for c in cats}
    handles = [Patch(fc=color_map.get(c, "#D4D4D4"),
                     label=labels.get(c, c), linewidth=0)
               for c in cats]
    handles += [
        Patch(fc="white", ec="black", lw=0.4, label="hap1 (left bar)"),
        Patch(fc="white", ec="black", lw=0.4, label="hap2 (right bar)"),
    ]
    fig.legend(handles=handles, fontsize=TICK - 1,
               loc="lower center", bbox_to_anchor=(0.5, -0.12),
               ncol=min(len(handles), 6), frameon=False,
               handlelength=0.8, handleheight=0.6,
               handletextpad=0.4, borderpad=0.2, labelspacing=0.25)
    _save(fig, f"fig1_{ds_label}_vertical_{run_label}")


# 1a. HiFi per-chromosome bars
cats_hifi = [c for c in hifi_cats_present if mat_hifi[c].sum() > 0]
if "uncovered" in mat_hifi.columns and mat_hifi["uncovered"].sum() > 0:
    cats_hifi.append("uncovered")
plot1_horizontal(mat_hifi, cats_hifi, NUCFLAG_COLORS, "nucflag_hifi")
plot1_vertical(mat_hifi, cats_hifi, NUCFLAG_COLORS, "nucflag_hifi")

# 1b. ONT per-chromosome bars
cats_ont = [c for c in ont_cats_present if mat_ont[c].sum() > 0]
if "uncovered" in mat_ont.columns and mat_ont["uncovered"].sum() > 0:
    cats_ont.append("uncovered")
plot1_horizontal(mat_ont, cats_ont, NUCFLAG_COLORS, "nucflag_ont")
plot1_vertical(mat_ont, cats_ont, NUCFLAG_COLORS, "nucflag_ont")

# 1c. Combined quality per-chromosome bars
plot1_horizontal(mat_qual, QUALITY_ORDER, QUALITY_COLORS,
                 "quality_combined", QUALITY_LABELS)
plot1_vertical(mat_qual, QUALITY_ORDER, QUALITY_COLORS,
               "quality_combined", QUALITY_LABELS)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Concordance confusion matrix (HiFi × ONT, bp)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 2: concordance confusion matrix ──")


def plot2_concordance():
    """Heatmap: rows = HiFi categories, cols = ONT categories, values = bp."""
    # Build confusion matrix
    conf = (df_ix.groupby(["hifi_name", "ont_name"])["overlap_bp"]
                 .sum().reset_index())
    # All categories in order
    cats_h = [c for c in all_cats if c in df_ix["hifi_name"].unique()]
    cats_o = [c for c in all_cats if c in df_ix["ont_name"].unique()]

    mat = np.zeros((len(cats_h), len(cats_o)))
    for _, r in conf.iterrows():
        if r["hifi_name"] in cats_h and r["ont_name"] in cats_o:
            i = cats_h.index(r["hifi_name"])
            j = cats_o.index(r["ont_name"])
            mat[i, j] = r["overlap_bp"]

    # Log scale for better visibility (add 1 to avoid log(0))
    mat_log = np.log10(mat + 1)

    fig, ax = plt.subplots(figsize=(4, 3.5), constrained_layout=True)
    im = ax.imshow(mat_log, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(cats_o)))
    ax.set_xticklabels(cats_o, fontsize=TICK - 1, rotation=60, ha="right")
    ax.set_yticks(range(len(cats_h)))
    ax.set_yticklabels(cats_h, fontsize=TICK - 1)
    ax.set_xlabel("ONT category", fontsize=LABEL)
    ax.set_ylabel("HiFi category", fontsize=LABEL)

    # Annotate cells with bp values
    for i in range(len(cats_h)):
        for j in range(len(cats_o)):
            v = mat[i, j]
            if v > 0:
                txt = f"{v / 1000:.1f}" if v >= 1000 else f"{int(v)}"
                color = "white" if mat_log[i, j] > mat_log.max() * 0.6 else "black"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=TICK - 2, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.ax.tick_params(labelsize=TICK - 1)
    cbar.set_label("log10(bp + 1)", fontsize=LABEL)
    _save(fig, f"fig2_concordance_heatmap_{run_label}")

    # Also save as TSV
    df_conf = pd.DataFrame(mat, index=cats_h, columns=cats_o)
    df_conf.index.name = "hifi_category"
    df_conf.to_csv(os.path.join(stats_dir, "concordance_matrix_bp.tsv"),
                   sep="\t")
    return mat, cats_h, cats_o


conf_mat, conf_cats_h, conf_cats_o = plot2_concordance()


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — Combined quality summary per haplotype
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 3: combined quality summary ──")


def plot3_quality_summary():
    """Stacked horizontal bars: one per haplotype, coloured by quality tier."""
    fig, ax = plt.subplots(figsize=(3.5, 1.2), constrained_layout=True)
    yy = np.arange(len(HAP_ORDER))
    bar_h = 0.5

    for i, hap in enumerate(HAP_ORDER):
        sub = mat_qual.loc[mat_qual["hap"] == hap]
        left = 0.0
        for q in QUALITY_ORDER:
            v = sub[q].sum() / 1000.0  # kbp
            ax.barh(yy[i], v, left=left, height=bar_h,
                    color=QUALITY_COLORS[q], edgecolor="none", alpha=0.9)
            left += v

    ax.set_yticks(yy)
    ax.set_yticklabels(HAP_ORDER, fontsize=TICK)
    ax.set_xlabel("Total telomere length (kbp)", fontsize=LABEL)
    ax.tick_params(axis="x", labelsize=TICK)
    _box(ax)

    handles = [Patch(fc=QUALITY_COLORS[q], label=QUALITY_LABELS[q], linewidth=0)
               for q in QUALITY_ORDER]
    fig.legend(handles=handles, fontsize=TICK - 1,
               loc="lower center", bbox_to_anchor=(0.5, -0.30),
               ncol=2, frameon=False, handlelength=0.6, handleheight=0.5,
               handletextpad=0.3, borderpad=0.2, labelspacing=0.2)
    _save(fig, f"fig3_quality_summary_{run_label}")


plot3_quality_summary()


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 4a — NucFlag composition heatmap (fraction per chr, both techs)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 4a: NucFlag composition heatmaps ──")


def plot4a_heatmap(mat, cats, tech_label):
    """Heatmap: rows = chr×arm (summed across haps), columns = NucFlag categories."""
    grp = mat.groupby(["chr_display", "arm"])[cats].sum()
    grp["total"] = grp.sum(axis=1)
    frac = grp[cats].div(grp["total"], axis=0).fillna(0)

    row_order = []
    for c in chroms:
        for a in ("p", "q"):
            if (c, a) in frac.index:
                row_order.append((c, a))
    frac = frac.loc[row_order]
    row_labels = [f"{c} {a}" for c, a in row_order]

    nonzero = [c for c in cats if frac[c].sum() > 0]
    heat = frac[nonzero].values

    fig, ax = plt.subplots(figsize=(2.5, 5), constrained_layout=True)
    im = ax.imshow(heat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(len(nonzero)))
    ax.set_xticklabels(nonzero, fontsize=TICK - 1, rotation=60, ha="right")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=TICK - 1)
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=TICK - 1)
    cbar.set_label("Fraction", fontsize=LABEL)
    _box(ax)
    _save(fig, f"fig4a_heatmap_{tech_label}_{run_label}")


cats_hifi_nz = [c for c in hifi_cats_present if mat_hifi[c].sum() > 0]
cats_ont_nz  = [c for c in ont_cats_present if mat_ont[c].sum() > 0]
plot4a_heatmap(mat_hifi, cats_hifi_nz, "hifi")
plot4a_heatmap(mat_ont, cats_ont_nz, "ont")


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 4b — ONT vs HiFi proportion scatter
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 4b: ONT vs HiFi scatter ──")


def plot4b_scatter():
    """Scatter: fraction of total telomere bp per NucFlag category, ONT vs HiFi."""
    # Compute fractions from intersection data
    hifi_sums = df_ix.groupby("hifi_name")["overlap_bp"].sum()
    ont_sums  = df_ix.groupby("ont_name")["overlap_bp"].sum()
    hifi_total = hifi_sums.sum()
    ont_total  = ont_sums.sum()

    cats_union = sorted(set(hifi_sums.index) | set(ont_sums.index),
                        key=lambda c: list(NUCFLAG_COLORS.keys()).index(c)
                        if c in NUCFLAG_COLORS else 999)

    fig, ax = plt.subplots(figsize=(3, 3), constrained_layout=True)
    for cat in cats_union:
        x = hifi_sums.get(cat, 0) / hifi_total
        y = ont_sums.get(cat, 0) / ont_total
        ax.scatter(x, y, color=NUCFLAG_COLORS.get(cat, "#D4D4D4"), s=30,
                   edgecolors="black", linewidths=0.3, zorder=3)
        ax.annotate(cat, (x, y), fontsize=TICK - 1.5,
                    textcoords="offset points", xytext=(3, 3))

    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], ls="--", lw=0.5, color="grey", zorder=1)
    ax.set_xlabel("HiFi fraction", fontsize=LABEL)
    ax.set_ylabel("ONT fraction", fontsize=LABEL)
    ax.tick_params(labelsize=TICK)
    _box(ax)
    _save(fig, f"fig4b_ont_vs_hifi_{run_label}")


plot4b_scatter()


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 4c — Technology rescue per telomere
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 4c: technology rescue ──")


def plot4c_rescue():
    """
    Per-chromosome bar chart: bp where one technology calls 'correct'
    but the other calls an error category.
    """
    rescue_data = []
    for _, telo in df_telo.iterrows():
        chr_d = telo["chr_display"]
        arm   = telo["arm"]
        hap   = telo["hap"]
        tlen  = telo["length"]
        sub   = df_ix.loc[(df_ix["chr_display"] == chr_d) &
                          (df_ix["arm"] == arm) & (df_ix["hap"] == hap)]
        hifi_rescue = int(sub.loc[(sub["hifi_name"] != "correct") &
                                  (sub["ont_name"] == "correct"),
                                  "overlap_bp"].sum())
        ont_rescue  = int(sub.loc[(sub["ont_name"] != "correct") &
                                  (sub["hifi_name"] == "correct"),
                                  "overlap_bp"].sum())
        rescue_data.append({
            "chr_display": chr_d, "arm": arm, "hap": hap,
            "telo_length_bp": tlen,
            "hifi_error_ont_correct_bp": hifi_rescue,
            "ont_error_hifi_correct_bp": ont_rescue,
            "hifi_error_ont_correct_frac": hifi_rescue / tlen if tlen else 0,
            "ont_error_hifi_correct_frac": ont_rescue / tlen if tlen else 0,
        })

    df_rescue = pd.DataFrame(rescue_data)

    # Bar chart: per chromosome (summed across arms & haps)
    chr_rescue = df_rescue.groupby("chr_display").agg({
        "hifi_error_ont_correct_bp": "sum",
        "ont_error_hifi_correct_bp": "sum",
        "telo_length_bp": "sum",
    }).reindex(chroms).fillna(0)

    fig, ax = plt.subplots(figsize=(5, 2.5), constrained_layout=True)
    xx = np.arange(n_chr)
    w = 0.35
    ax.bar(xx - w / 2, chr_rescue["ont_error_hifi_correct_bp"] / 1000,
           width=w, color=QUALITY_COLORS["hifi_only_correct"],
           label="HiFi correct, ONT error", alpha=0.9)
    ax.bar(xx + w / 2, chr_rescue["hifi_error_ont_correct_bp"] / 1000,
           width=w, color=QUALITY_COLORS["ont_only_correct"],
           label="ONT correct, HiFi error", alpha=0.9)
    ax.set_xticks(xx)
    ax.set_xticklabels(chr_short, fontsize=TICK, rotation=0)
    ax.set_ylabel("Rescued telomere length (kbp)", fontsize=LABEL)
    ax.tick_params(axis="y", labelsize=TICK)
    ax.legend(fontsize=TICK - 1, frameon=False)
    _box(ax)
    _save(fig, f"fig4c_rescue_{run_label}")

    return df_rescue


df_rescue = plot4c_rescue()


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 4d — Telomere quality ranking
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Plot 4d: telomere quality ranking ──")


def plot4d_ranking():
    """Horizontal bars: all telomeres ranked by % both_correct."""
    ranking = []
    for _, telo in df_telo.iterrows():
        chr_d = telo["chr_display"]
        arm   = telo["arm"]
        hap   = telo["hap"]
        tlen  = telo["length"]
        sub = df_ix.loc[(df_ix["chr_display"] == chr_d) &
                        (df_ix["arm"] == arm) & (df_ix["hap"] == hap)]
        both_c = int(sub.loc[(sub["hifi_name"] == "correct") &
                             (sub["ont_name"] == "correct"),
                             "overlap_bp"].sum())
        ranking.append({
            "chr_display": chr_d, "arm": arm, "hap": hap,
            "telo_length_bp": tlen,
            "both_correct_bp": both_c,
            "both_correct_frac": both_c / tlen if tlen else 0,
            "label": f"{chr_d} {arm} {hap}",
        })

    df_rank = pd.DataFrame(ranking).sort_values("both_correct_frac",
                                                 ascending=True)
    n = len(df_rank)
    fig, ax = plt.subplots(figsize=(3.5, max(4, n * 0.12)),
                           constrained_layout=True)
    yy = np.arange(n)

    # Background: total length
    ax.barh(yy, df_rank["telo_length_bp"].values / 1000, height=0.7,
            color="#E0E0E0", edgecolor="none")
    # Foreground: both_correct
    ax.barh(yy, df_rank["both_correct_bp"].values / 1000, height=0.7,
            color=QUALITY_COLORS["both_correct"], edgecolor="none", alpha=0.9)

    ax.set_yticks(yy)
    ax.set_yticklabels(df_rank["label"].values, fontsize=TICK - 1.5)
    ax.set_xlabel("Telomere length (kbp)", fontsize=LABEL)
    ax.tick_params(axis="x", labelsize=TICK)
    _box(ax)

    handles = [
        Patch(fc=QUALITY_COLORS["both_correct"], label="Both correct"),
        Patch(fc="#E0E0E0", label="Total telomere"),
    ]
    ax.legend(handles=handles, fontsize=TICK - 1, frameon=False,
              loc="lower right")
    _save(fig, f"fig4d_quality_ranking_{run_label}")
    return df_rank


df_rank = plot4d_ranking()


# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
#                         S T A T I S T I C S
# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STATISTICS — H9 Telomere × NucFlag Combined Analysis")
print("=" * 80)


# ── 7a. Per-telomere breakdown (both techs, bp + fraction) ──────────────────
print("\n" + "─" * 80)
print("7a. PER-TELOMERE BREAKDOWN (both technologies)")
print("─" * 80)

rows_7a: list = []
for _, telo in df_telo.iterrows():
    chr_d = telo["chr_display"]
    arm   = telo["arm"]
    hap   = telo["hap"]
    tlen  = telo["length"]

    sub = df_ix.loc[(df_ix["chr_display"] == chr_d) &
                    (df_ix["arm"] == arm) & (df_ix["hap"] == hap)]

    entry = {
        "chromosome": chr_d, "arm": arm, "haplotype": hap,
        "telomere_length_bp": tlen,
    }

    # HiFi categories
    for cat in all_cats:
        bp = int(sub.loc[sub["hifi_name"] == cat, "overlap_bp"].sum())
        entry[f"hifi_{cat}_bp"]   = bp
        entry[f"hifi_{cat}_frac"] = f"{bp / tlen:.6f}" if tlen else "NA"

    # ONT categories
    for cat in all_cats:
        bp = int(sub.loc[sub["ont_name"] == cat, "overlap_bp"].sum())
        entry[f"ont_{cat}_bp"]   = bp
        entry[f"ont_{cat}_frac"] = f"{bp / tlen:.6f}" if tlen else "NA"

    # Quality tiers
    for q in QUALITY_ORDER:
        bp = int(sub.loc[sub["quality"] == q, "overlap_bp"].sum())
        entry[f"{q}_bp"]   = bp
        entry[f"{q}_frac"] = f"{bp / tlen:.6f}" if tlen else "NA"

    rows_7a.append(entry)

df_7a = pd.DataFrame(rows_7a)
_write_section("7a_per_telomere_breakdown.tsv", df_7a,
               "Per-telomere NucFlag breakdown (HiFi + ONT + quality tiers), "
               "absolute bp and fraction of telomere length.")

# Console summary
for hap in HAP_ORDER:
    sub = df_7a.loc[df_7a["haplotype"] == hap]
    total = sub["telomere_length_bp"].sum()
    n_telo = len(sub)
    bc = sub["both_correct_bp"].sum()
    print(f"  {hap}: {n_telo} telomeres, {total:,} bp total, "
          f"{bc:,} bp both_correct ({bc / total:.4f})")


# ── 7b. Global category summary ────────────────────────────────────────────
print("\n" + "─" * 80)
print("7b. GLOBAL CATEGORY SUMMARY")
print("─" * 80)

rows_7b: list = []
for tech, col_name in [("HiFi", "hifi_name"), ("ONT", "ont_name")]:
    sums = df_ix.groupby(col_name)["overlap_bp"].sum()
    total = sums.sum()
    print(f"\n  {tech} (total: {total:,} bp)")
    print(f"  {'Category':<20s} {'Total (bp)':>14s} {'Fraction':>10s}")
    print(f"  {'─' * 20} {'─' * 14} {'─' * 10}")
    for cat in all_cats:
        bp = int(sums.get(cat, 0))
        frac = bp / total if total else 0
        rows_7b.append({
            "technology": tech, "category": cat,
            "total_bp": bp, "fraction": f"{frac:.6f}",
            "total_bp_relative_teloscope": bp,
            "fraction_of_teloscope": f"{bp / total_telo_bp:.6f}",
        })
        if bp > 0:
            print(f"  {cat:<20s} {bp:>14,d} {frac:>10.4f}")

df_7b = pd.DataFrame(rows_7b)
_write_section("7b_global_category_summary.tsv", df_7b,
               "Global NucFlag category summary per technology, with absolute "
               "bp and fraction of total telomere length (Teloscope-defined).")


# ── 7c. Concordance statistics ──────────────────────────────────────────────
print("\n" + "─" * 80)
print("7c. CONCORDANCE STATISTICS")
print("─" * 80)

rows_7c: list = []
conf = df_ix.groupby(["hifi_name", "ont_name"])["overlap_bp"].sum().reset_index()
for _, r in conf.iterrows():
    bp = int(r["overlap_bp"])
    rows_7c.append({
        "hifi_category": r["hifi_name"],
        "ont_category":  r["ont_name"],
        "overlap_bp":    bp,
        "fraction_of_total_telomere": f"{bp / total_telo_bp:.6f}",
        "concordant":    r["hifi_name"] == r["ont_name"],
    })

df_7c = pd.DataFrame(rows_7c).sort_values("overlap_bp", ascending=False)

conc_bp   = int(df_7c.loc[df_7c["concordant"], "overlap_bp"].sum())
disc_bp   = int(df_7c.loc[~df_7c["concordant"], "overlap_bp"].sum())
print(f"  Concordant:  {conc_bp:>12,d} bp  ({conc_bp / total_telo_bp:.4f})")
print(f"  Discordant:  {disc_bp:>12,d} bp  ({disc_bp / total_telo_bp:.4f})")
print(f"  Total:       {total_telo_bp:>12,d} bp")

_write_section("7c_concordance.tsv", df_7c,
               "Full concordance matrix: every (HiFi × ONT) category pair "
               "with absolute bp and fraction of Teloscope telomere length.")


# ── 7d. Technology rescue ───────────────────────────────────────────────────
print("\n" + "─" * 80)
print("7d. TECHNOLOGY RESCUE (discordant correct calls)")
print("─" * 80)

total_hifi_resc = int(df_rescue["ont_error_hifi_correct_bp"].sum())
total_ont_resc  = int(df_rescue["hifi_error_ont_correct_bp"].sum())
print(f"  HiFi correct, ONT error:  {total_hifi_resc:>10,d} bp  "
      f"({total_hifi_resc / total_telo_bp:.4f})")
print(f"  ONT correct, HiFi error:  {total_ont_resc:>10,d} bp  "
      f"({total_ont_resc / total_telo_bp:.4f})")

_write_section("7d_rescue.tsv", df_rescue,
               "Per-telomere 'rescue' analysis: bp where one technology "
               "calls correct but the other calls an error.")

# Top rescued telomeres
print("\n  Top 10 telomeres rescued by HiFi (HiFi=correct, ONT=error):")
top_hifi = df_rescue.nlargest(10, "ont_error_hifi_correct_bp")
for _, r in top_hifi.iterrows():
    if r["ont_error_hifi_correct_bp"] > 0:
        print(f"    {r['chr_display']} {r['arm']} {r['hap']}: "
              f"{int(r['ont_error_hifi_correct_bp']):,d} bp "
              f"({r['ont_error_hifi_correct_frac']:.4f} of telomere)")

print("\n  Top 10 telomeres rescued by ONT (ONT=correct, HiFi=error):")
top_ont = df_rescue.nlargest(10, "hifi_error_ont_correct_bp")
for _, r in top_ont.iterrows():
    if r["hifi_error_ont_correct_bp"] > 0:
        print(f"    {r['chr_display']} {r['arm']} {r['hap']}: "
              f"{int(r['hifi_error_ont_correct_bp']):,d} bp "
              f"({r['hifi_error_ont_correct_frac']:.4f} of telomere)")


# ── 7e. Hypergeometric enrichment ──────────────────────────────────────────
print("\n" + "─" * 80)
print("7e. HYPERGEOMETRIC ENRICHMENT TEST")
print("    Tests each NucFlag category enrichment in telomeric ends")
print("    vs. the full terminal NucFlag window.")
print("─" * 80)

# Full window = all of df_nf (excluding the giant uncovered middle)
df_nf_terminal = df_nf.loc[~((df_nf["hifi_name"] == "uncovered") &
                              (df_nf["ont_name"] == "uncovered"))].copy()

hyper_rows: list = []
for tech, col_name in [("HiFi", "hifi_name"), ("ONT", "ont_name")]:
    # Total bp in NucFlag window (terminal regions only)
    total_window_bp = int(df_nf_terminal["length"].sum())
    # Total telomere bp
    total_tel_bp = int(df_ix.loc[~((df_ix["hifi_name"] == "uncovered") &
                                    (df_ix["ont_name"] == "uncovered")),
                                   "overlap_bp"].sum())

    cats_tech = sorted(
        df_ix[col_name].unique(),
        key=lambda c: list(NUCFLAG_COLORS.keys()).index(c)
        if c in NUCFLAG_COLORS else 999
    )

    print(f"\n  {tech} (window: {total_window_bp:,} bp, "
          f"telomere: {total_tel_bp:,} bp)")
    print(f"  {'Category':<20s} {'obs (bp)':>12s} {'exp (bp)':>12s} "
          f"{'fold':>8s} {'p-value':>12s}")
    print(f"  {'─' * 20} {'─' * 12} {'─' * 12} {'─' * 8} {'─' * 12}")

    for cat in cats_tech:
        if cat == "uncovered":
            continue
        # K = total bp of this category in full window
        K = int(df_nf_terminal.loc[df_nf_terminal[col_name] == cat,
                                    "length"].sum())
        # k = observed bp of this category in telomeres
        k = int(df_ix.loc[df_ix[col_name] == cat, "overlap_bp"].sum())
        N = total_window_bp
        n = total_tel_bp

        if N > 0 and K > 0 and n > 0:
            pval = sp_stats.hypergeom.sf(k - 1, N, K, n)
            expected = n * K / N
            fold = k / expected if expected > 0 else float("inf")
            sig = " *" if pval < 0.05 else ""
            hyper_rows.append({
                "technology": tech, "category": cat,
                "observed_bp": k, "expected_bp": f"{expected:.1f}",
                "K_category_in_window": K, "N_window_total": N,
                "n_telomere_total": n,
                "fold_enrichment": f"{fold:.3f}",
                "p_value": f"{pval:.4e}",
                "significant_0.05": "yes" if pval < 0.05 else "no",
            })
            print(f"  {cat:<20s} {k:>12,d} {expected:>12,.1f} "
                  f"{fold:>8.2f} {pval:>12.4e}{sig}")

if hyper_rows:
    _write_section("7e_hypergeometric_enrichment.tsv",
                   pd.DataFrame(hyper_rows),
                   "Hypergeometric enrichment: NucFlag categories in telomeric "
                   "ends vs full terminal window, per technology.")


# ── 7f. ONT vs HiFi Fisher exact test ──────────────────────────────────────
print("\n" + "─" * 80)
print("7f. ONT vs HiFi CATEGORY PROPORTIONS (Fisher exact)")
print("─" * 80)

hifi_sums = df_ix.groupby("hifi_name")["overlap_bp"].sum()
ont_sums  = df_ix.groupby("ont_name")["overlap_bp"].sum()
hifi_total_bp = int(hifi_sums.sum())
ont_total_bp  = int(ont_sums.sum())

fisher_rows: list = []
cats_union = sorted(set(hifi_sums.index) | set(ont_sums.index),
                    key=lambda c: list(NUCFLAG_COLORS.keys()).index(c)
                    if c in NUCFLAG_COLORS else 999)

print(f"\n  {'Category':<20s} {'HiFi frac':>10s} {'ONT frac':>10s} "
      f"{'OR':>8s} {'p-value':>12s}")
print(f"  {'─' * 20} {'─' * 10} {'─' * 10} {'─' * 8} {'─' * 12}")

for cat in cats_union:
    a = int(hifi_sums.get(cat, 0))
    b = hifi_total_bp - a
    c = int(ont_sums.get(cat, 0))
    d = ont_total_bp - c

    if min(a + c, a + b) > 0:
        odds, pval = sp_stats.fisher_exact([[a, b], [c, d]],
                                           alternative="two-sided")
        hifi_f = a / (a + b) if (a + b) else 0
        ont_f  = c / (c + d) if (c + d) else 0
        sig = " *" if pval < 0.05 else ""
        fisher_rows.append({
            "category": cat,
            "hifi_bp": a, "hifi_other_bp": b,
            "ont_bp": c, "ont_other_bp": d,
            "hifi_fraction": f"{hifi_f:.6f}",
            "ont_fraction": f"{ont_f:.6f}",
            "odds_ratio": f"{odds:.4f}",
            "p_value": f"{pval:.4e}",
            "significant_0.05": "yes" if pval < 0.05 else "no",
        })
        print(f"  {cat:<20s} {hifi_f:>10.4f} {ont_f:>10.4f} "
              f"{odds:>8.3f} {pval:>12.4e}{sig}")

if fisher_rows:
    _write_section("7f_ont_vs_hifi_fisher.tsv", pd.DataFrame(fisher_rows),
                   "Fisher exact test: NucFlag category proportions "
                   "HiFi vs ONT (2×2 contingency, bp level).")


# ── 7g. Chromosome-type breakdown ──────────────────────────────────────────
print("\n" + "─" * 80)
print("7g. CHROMOSOME-TYPE BREAKDOWN")
print("─" * 80)

chrtype_rows: list = []
for tech, col_name in [("HiFi", "hifi_name"), ("ONT", "ont_name")]:
    df_ct = df_ix.copy()
    df_ct["chr_type"] = df_ct["chr_display"].map(_classify_chr)

    print(f"\n  {tech}")
    for ctype in ["metacentric", "submetacentric", "acrocentric"]:
        sub = df_ct.loc[df_ct["chr_type"] == ctype]
        total = int(sub["overlap_bp"].sum())
        n_c = sub["chr_display"].nunique()
        print(f"    {ctype} ({n_c} chr, total: {total:,} bp):")

        cat_sums = sub.groupby(col_name)["overlap_bp"].sum().sort_values(
            ascending=False)
        for cat, bp in cat_sums.items():
            bp = int(bp)
            frac = bp / total if total else 0
            chrtype_rows.append({
                "technology": tech, "chr_type": ctype,
                "category": cat, "total_bp": bp,
                "fraction": f"{frac:.6f}",
                "fraction_of_teloscope": f"{bp / total_telo_bp:.6f}",
            })
            if bp > 0:
                print(f"      {cat:<20s} {bp:>10,d} bp  ({frac:.4f})")

if chrtype_rows:
    _write_section("7g_chrtype_breakdown.tsv", pd.DataFrame(chrtype_rows),
                   "NucFlag category breakdown by chromosome morphology "
                   "(metacentric / submetacentric / acrocentric).")


# ── 7h. Quality ranking (best / worst telomeres) ───────────────────────────
print("\n" + "─" * 80)
print("7h. TELOMERE QUALITY RANKING")
print("─" * 80)

df_rank_out = df_rank.copy()
df_rank_out["both_correct_frac_of_teloscope"] = (
    df_rank_out["both_correct_bp"] / total_telo_bp
)
df_rank_out = df_rank_out.sort_values("both_correct_frac", ascending=False)

_write_section("7h_quality_ranking.tsv",
               df_rank_out[["chr_display", "arm", "hap", "telo_length_bp",
                            "both_correct_bp", "both_correct_frac",
                            "both_correct_frac_of_teloscope"]],
               "All telomeres ranked by fraction of 'both correct' bp.")

# Print best and worst 10
print("\n  TOP 10 — highest quality telomeres (% both correct):")
for i, (_, r) in enumerate(df_rank_out.head(10).iterrows()):
    print(f"    {i+1:>2d}. {r['chr_display']} {r['arm']} {r['hap']}: "
          f"{r['both_correct_frac']:.4f}  "
          f"({int(r['both_correct_bp']):,d} / {int(r['telo_length_bp']):,d} bp)")

print("\n  BOTTOM 10 — lowest quality telomeres:")
for i, (_, r) in enumerate(df_rank_out.tail(10).iloc[::-1].iterrows()):
    print(f"    {i+1:>2d}. {r['chr_display']} {r['arm']} {r['hap']}: "
          f"{r['both_correct_frac']:.4f}  "
          f"({int(r['both_correct_bp']):,d} / {int(r['telo_length_bp']):,d} bp)")

# ── 7i. Assembly quality summary (the headline numbers) ────────────────────
print("\n" + "─" * 80)
print("7i. ASSEMBLY QUALITY SUMMARY (headline numbers for the paper)")
print("─" * 80)

bc_bp = int(mat_qual["both_correct"].sum())
ho_bp = int(mat_qual["hifi_only_correct"].sum())
oo_bp = int(mat_qual["ont_only_correct"].sum())
be_bp = int(mat_qual["both_error"].sum())
any_correct_bp = bc_bp + ho_bp + oo_bp

print(f"\n  Total telomere length (Teloscope):     {total_telo_bp:>12,d} bp")
print(f"  Both technologies correct:             {bc_bp:>12,d} bp  "
      f"({bc_bp / total_telo_bp:.4f})")
print(f"  At least one technology correct:        {any_correct_bp:>12,d} bp  "
      f"({any_correct_bp / total_telo_bp:.4f})")
print(f"    - HiFi correct only (ONT error):     {ho_bp:>12,d} bp  "
      f"({ho_bp / total_telo_bp:.4f})")
print(f"    - ONT correct only (HiFi error):     {oo_bp:>12,d} bp  "
      f"({oo_bp / total_telo_bp:.4f})")
print(f"  Both technologies error:               {be_bp:>12,d} bp  "
      f"({be_bp / total_telo_bp:.4f})")

# Per-haplotype
for hap in HAP_ORDER:
    sub = mat_qual.loc[mat_qual["hap"] == hap]
    t = int(df_telo.loc[df_telo["hap"] == hap, "length"].sum())
    bc = int(sub["both_correct"].sum())
    ac = bc + int(sub["hifi_only_correct"].sum()) + int(sub["ont_only_correct"].sum())
    print(f"\n  {hap}: total {t:,} bp, both_correct {bc:,} bp ({bc / t:.4f}), "
          f"any_correct {ac:,} bp ({ac / t:.4f})")

# Perfect telomeres (100% both_correct)
perfect = df_rank_out.loc[df_rank_out["both_correct_frac"] >= 0.999]
print(f"\n  'Perfect' telomeres (≥99.9% both correct): {len(perfect)} / {len(df_rank_out)}")
for _, r in perfect.iterrows():
    print(f"    {r['chr_display']} {r['arm']} {r['hap']}: "
          f"{int(r['telo_length_bp']):,d} bp")

# Worst telomeres (< 50% both_correct)
worst = df_rank_out.loc[df_rank_out["both_correct_frac"] < 0.5]
print(f"\n  Problematic telomeres (<50% both correct): {len(worst)} / {len(df_rank_out)}")
for _, r in worst.sort_values("both_correct_frac").iterrows():
    print(f"    {r['chr_display']} {r['arm']} {r['hap']}: "
          f"{r['both_correct_frac']:.4f} "
          f"({int(r['both_correct_bp']):,d} / {int(r['telo_length_bp']):,d} bp)")

# Save summary
summary_rows = [
    {"metric": "total_telomere_bp", "value_bp": total_telo_bp,
     "fraction": "1.000000"},
    {"metric": "both_correct_bp", "value_bp": bc_bp,
     "fraction": f"{bc_bp / total_telo_bp:.6f}"},
    {"metric": "any_correct_bp", "value_bp": any_correct_bp,
     "fraction": f"{any_correct_bp / total_telo_bp:.6f}"},
    {"metric": "hifi_only_correct_bp", "value_bp": ho_bp,
     "fraction": f"{ho_bp / total_telo_bp:.6f}"},
    {"metric": "ont_only_correct_bp", "value_bp": oo_bp,
     "fraction": f"{oo_bp / total_telo_bp:.6f}"},
    {"metric": "both_error_bp", "value_bp": be_bp,
     "fraction": f"{be_bp / total_telo_bp:.6f}"},
    {"metric": "n_telomeres_total", "value_bp": len(df_rank_out),
     "fraction": "NA"},
    {"metric": "n_perfect_telomeres", "value_bp": len(perfect),
     "fraction": f"{len(perfect) / len(df_rank_out):.6f}"},
    {"metric": "n_problematic_telomeres", "value_bp": len(worst),
     "fraction": f"{len(worst) / len(df_rank_out):.6f}"},
]
_write_section("7i_assembly_quality_summary.tsv", pd.DataFrame(summary_rows),
               "Headline assembly quality numbers for the paper.")


# ── Final combined stats TSV ────────────────────────────────────────────────
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
    stats_out = os.path.join(plots_dir,
                             f"h9_nucflag_telomere_stats_{run_label}.tsv")
    df_combined.to_csv(stats_out, sep="\t", index=False)
    print(f"\n[SAVED] Combined stats: {stats_out}  ({len(df_combined)} rows)")

print("\n" + "=" * 80)
print("[DONE] H9 telomere NucFlag combined analysis.")
print(f"       Figures: {plots_dir}")
print(f"       Stats:   {stats_dir}")
print("=" * 80)
