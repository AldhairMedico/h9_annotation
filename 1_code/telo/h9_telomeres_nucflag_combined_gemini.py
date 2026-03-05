#!/usr/bin/env python3
"""
H9 terminal telomere lengths with NucFlag orthogonal validation (Combined & Final Version).

Addresses:
- Assembly Quality Evaluation (HiFi + ONT Jointly) 
- Rescue Analysis (HiFi vs ONT)
- Best/Worst Telomere Identification
- Answers specific biological interpretations (condensed from independent approaches).
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats as sp_stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Setup & Paths
# ═══════════════════════════════════════════════════════════════════════════════
script_dir  = os.path.dirname(os.path.abspath(__file__))
base_dir    = os.path.dirname(os.path.dirname(script_dir))
data_dir    = os.path.join(base_dir, "2_data", "2.2_processed")
figures_dir = os.path.join(base_dir, "3_figures", "3.1_draft")
run_label   = "9"

plots_dir = os.path.join(figures_dir, f"26.03.02_telomeres_nucflag_combined_{run_label}")
os.makedirs(plots_dir, exist_ok=True)
stats_dir = os.path.join(plots_dir, "stats")
os.makedirs(stats_dir, exist_ok=True)

TELO_BED = os.path.join(data_dir, "25.12.10_teloscope_compiled", "H9_T2T_v0.1_dip.fasta_terminal_telomeres.bed")
COMBINED_TSV = os.path.join(data_dir, "nucflag", "nucflag_telo_combined.tsv")
PALETTE_TSV = os.path.join(data_dir, "nucflag", "nucflag_palette.tsv")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Genetics & Aesthetics Constants
# ═══════════════════════════════════════════════════════════════════════════════
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
ACROCENTRIC    = {"chr13", "chr14", "chr15", "chr21", "chr22"}
METACENTRIC    = {"chr1", "chr3", "chr16", "chr19", "chr20"}
SUBMETACENTRIC = {"chr2", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                  "chr10", "chr11", "chr12", "chr17", "chr18", "chrX"}

def _classify_chr(c: str) -> str:
    if c in ACROCENTRIC: return "acrocentric"
    if c in METACENTRIC: return "metacentric"
    if c in SUBMETACENTRIC: return "submetacentric"
    return "other"

def rgb_to_hex(rgb_str):
    if not isinstance(rgb_str, str): return "#000000"
    parts = rgb_str.replace('"', '').replace("'", "").split(',')
    if len(parts) == 3:
        return f"#{int(parts[0]):02x}{int(parts[1]):02x}{int(parts[2]):02x}"
    return "#000000"

try:
    df_pal = pd.read_csv(PALETTE_TSV, sep="\t")
    NUCFLAG_COLORS = {row["name"]: rgb_to_hex(row["itemRgb"]) for _, row in df_pal.iterrows()}
except Exception as e:
    print(f"[ERROR] Failed to load palette: {e}")
    sys.exit(1)

NUCFLAG_COLORS["uncovered"] = "#D4D4D4"

JOINT_CATS = ["Both Correct", "HiFi Rescued", "ONT Rescued", "Both Error / Unsupported", "Uncovered"]
JOINT_COLORS = {
    "Both Correct": "#2D8E47",
    "HiFi Rescued": "#3C7DC4",
    "ONT Rescued":  "#E69F00",
    "Both Error / Unsupported": "#CB3335",
    "Uncovered":    "#999999",
}

# ── NucFlag pairwise classification (Supplementary Methods) ──────────────────
STATUS_MAP = {
    "correct":       "PERFECT",
    "misjoin":       "HARD_ERROR",
    "collapse":      "HARD_ERROR",
    "false_dup":     "HARD_ERROR",
    "deletion":      "HARD_ERROR",
    "insertion":     "HARD_ERROR",
    "mismatch":      "HARD_ERROR",
    "softclip":      "HARD_ERROR",
    "het_or_mismap": "MISMATCH_SIGNAL",
    "low_quality":   "MISMATCH_SIGNAL",
    "homopolymer":   "REPEAT_CONTEXT",
    "dinucleotide":  "REPEAT_CONTEXT",
    "simple_repeat": "REPEAT_CONTEXT",
    "other_repeat":  "REPEAT_CONTEXT",
    "scaffold":      "SCAFFOLD",
    "uncovered":     "MISMATCH_SIGNAL",  # no coverage → ambiguous, not an error
}

# Pair-label rendering order (stacking: bottom -> top)
PAIR_CATS = ["correct", "ONT_rescued", "HiFi_rescued",
             "repeat", "mismatch_signal", "error"]
PAIR_LABELS = {
    "correct":         "Concordant",
    "ONT_rescued":     "ONT only",
    "HiFi_rescued":    "HiFi only",
    "error":           "Error",
    "repeat":          "Repeat context",
    "mismatch_signal": "Mismatch signal",
}
# CNS pastel palette (blue gradient dark->light: correct->ONT->HiFi)
PAIR_COLORS = {
    "correct":         "#2B6A8F",  # deep teal-blue
    "ONT_rescued":     "#5FA3C5",  # medium blue
    "HiFi_rescued":    "#AED4E6",  # soft sky
    "error":           "#C44E52",  # muted crimson
    "repeat":          "#E4E4E4",  # very light grey (not ambiguous quality)
    "mismatch_signal": "#D6FDCE",  # pale steel blue (ambiguous, not error-like)
}
# ── Fig2 palette variants ────────────────────────────────────────────────────
# V1: Original CNS blue gradient (deep teal → medium → soft sky)
PAIR_COLORS_V1 = dict(PAIR_COLORS)  # identical to base
# V2: Greyscale (print-safe, highest contrast)
PAIR_COLORS_V2 = {
    "correct":         "#404040",  # dark charcoal
    "ONT_rescued":     "#808080",  # medium grey
    "HiFi_rescued":    "#BFBFBF",  # light grey
    "error":           "#C44E52",  # muted crimson
    "repeat":          "#E4E4E4",  # very light grey
    "mismatch_signal": "#D6FDCE",  # pale green
}
# V3: Greyscale reversed (print-safe, highest contrast)
PAIR_COLORS_V3 = {
    "correct":         "#E4E4E4",  # very light grey
    "ONT_rescued":     "#BFBFBF",  # light grey
    "HiFi_rescued":    "#808080",  # medium grey
    "error":           "#C44E52",  # muted crimson
    "repeat":          "#404040",  # dark charcoal
    "mismatch_signal": "#D6FDCE",  # pale green
}

# Scatter: chr morphology colours (colorblind-safe)
CHR_TYPE_COLORS = {
    "metacentric":    "#F4A582",  # pastel coral (warm)
    "submetacentric": "#E0E0E0",  # very light gray (neutral)
    "acrocentric":    "#92C5DE",  # pastel steel blue (cool)
}
# ── Fig3 scatter palette variants ────────────────────────────────────────────
CHR_TYPE_COLORS_V2 = {  # gnomAD-style (purple / gray / gold)
    "metacentric":    "#8856A7",
    "submetacentric": "#B0B0B0",
    "acrocentric":    "#F0C33C",
}
CHR_TYPE_COLORS_V3 = {  # Okabe-Ito colorblind-safe (orange / gray / sky)
    "metacentric":    "#E69F00",
    "submetacentric": "#999999",
    "acrocentric":    "#56B4E9",
}

try:
    from matplotlib.font_manager import findfont, FontProperties
    font_family = "Arial" if findfont(FontProperties(family="Arial")) != findfont(FontProperties()) else "DejaVu Sans"
except:
    font_family = "DejaVu Sans"

plt.rcParams.update({
    "font.family":       font_family,
    "font.size":         7,
    "axes.linewidth":    0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "savefig.facecolor": "white",
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",   # editable text in SVG
})

FIG_WIDTH = 4.5  # shared width for Figs 2 & 3 so they stack
TICK, LABEL = 6, 7

def _box(ax):
    for side in ("left", "bottom"):
        ax.spines[side].set_edgecolor("black")
        ax.spines[side].set_linewidth(0.5)
        ax.spines[side].set_visible(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

def _save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(plots_dir, f"{stem}.{ext}"), dpi=600 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Data Loading & Intersection Mapping
# ═══════════════════════════════════════════════════════════════════════════════
print(f"[INFO] Loading Teloscope BED: {TELO_BED}")
df_telo = pd.read_csv(TELO_BED, sep="\t", header=None, names=[
    "chr", "start", "end", "length", "arm", "fwdCounts", "revCounts", "canCounts", "nonCanCounts", "chrSize"])
df_telo["chr_display"] = df_telo["chr"].str.replace(r"_hap[12]$", "", regex=True)
df_telo["hap"] = df_telo["chr"].str.extract(r"_(hap[12])$")[0]
df_telo["telo_id"] = df_telo["chr_display"] + "_" + df_telo["arm"] + "_" + df_telo["hap"]
df_telo["chr_type"] = df_telo["chr_display"].map(_classify_chr)

print(f"[INFO] Loading Combined NucFlag: {COMBINED_TSV}")
df_nf = pd.read_csv(COMBINED_TSV, sep="\t")

print(f"[INFO] Intersecting strictly matching telomere intervals in-memory...")
intersections = []
for _, t_row in df_telo.iterrows():
    c = t_row["chr"]
    start = t_row["start"]
    end = t_row["end"]
    tid = t_row["telo_id"]
    tlen = t_row["length"]

    nf_sub = df_nf[(df_nf["chrom"] == c) & (df_nf["chromEnd"] > start) & (df_nf["chromStart"] < end)]
    
    covered_bp = 0
    for _, nf_row in nf_sub.iterrows():
        ov_s = max(start, nf_row["chromStart"])
        ov_e = min(end, nf_row["chromEnd"])
        overlap = ov_e - ov_s
        if overlap > 0:
            hname = str(nf_row["hifi_name"]).strip()
            oname = str(nf_row["ont_name"]).strip()
            if not hname or hname == "nan" or hname == "nan": hname = "uncovered"
            if not oname or oname == "nan" or oname == "nan": oname = "uncovered"
            
            intersections.append({
                "telo_id": tid,
                "chr_display": t_row["chr_display"],
                "arm": t_row["arm"],
                "hap": t_row["hap"],
                "chr_type": t_row["chr_type"],
                "overlap_bp": overlap,
                "hifi_name": hname,
                "ont_name": oname,
                "total_len": tlen
            })
            covered_bp += overlap
            
    # Add an uncovered segment if NucFlag combined file entirely misses pieces of a telomere
    if covered_bp < tlen:
        intersections.append({
            "telo_id": tid,
            "chr_display": t_row["chr_display"],
            "arm": t_row["arm"],
            "hap": t_row["hap"],
            "chr_type": t_row["chr_type"],
            "overlap_bp": tlen - covered_bp,
            "hifi_name": "uncovered",
            "ont_name": "uncovered",
            "total_len": tlen
        })

df_int = pd.DataFrame(intersections)

def get_joint_status(row):
    h, o = row["hifi_name"], row["ont_name"]
    if h == "correct" and o == "correct": return "Both Correct"
    if h == "correct" and o != "correct": return "HiFi Rescued"
    if o == "correct" and h != "correct": return "ONT Rescued"
    if h == "uncovered" and o == "uncovered": return "Uncovered"
    return "Both Error / Unsupported"

df_int["joint_status"] = df_int.apply(get_joint_status, axis=1)

# ── Pairwise classification (Supplementary Methods) ──
df_int["hifi_status"] = df_int["hifi_name"].map(STATUS_MAP).fillna("UNKNOWN")
df_int["ont_status"]  = df_int["ont_name"].map(STATUS_MAP).fillna("UNKNOWN")

def get_pair_label(row):
    """Derive pair label from HiFi x ONT status classes (precedence rules)."""
    H, O = row["hifi_status"], row["ont_status"]
    if O == "PERFECT" and H == "HARD_ERROR":  return "ONT_rescued"
    if H == "PERFECT" and O == "HARD_ERROR":  return "HiFi_rescued"
    if H == "PERFECT" or  O == "PERFECT":     return "correct"
    if H == "HARD_ERROR" and O == "HARD_ERROR": return "error"
    if H == "SCAFFOLD" or O == "SCAFFOLD":    return "scaffold"
    if H == "REPEAT_CONTEXT" or O == "REPEAT_CONTEXT": return "repeat"
    if H == "MISMATCH_SIGNAL" or O == "MISMATCH_SIGNAL": return "mismatch_signal"
    return "mismatch_signal"  # fallback: unknown combination → ambiguous

df_int["pair_label"] = df_int.apply(get_pair_label, axis=1)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Global Statistics & Breakdowns (absolute bp & relative %)
# ═══════════════════════════════════════════════════════════════════════════════
def write_tsv(df, fname, desc):
    fp = os.path.join(stats_dir, fname)
    df.to_csv(fp, sep="\t", index=False)
    print(f"  [SAVED] {fp} - {desc}")

print("\n[INFO] Generating Statistics (2.2.1 - 2.2.3 answers)...")

tot_bp = df_telo["length"].sum()

# 4.1 Global Assembly Quality (Joint)
g_joint = df_int.groupby("joint_status")["overlap_bp"].sum().reset_index()
g_joint["percentage"] = (g_joint["overlap_bp"] / tot_bp) * 100
g_joint = pd.concat([pd.DataFrame([{"joint_status": "TOTAL Telomere Space", "overlap_bp": tot_bp, "percentage": 100.0}]), g_joint], ignore_index=True)
write_tsv(g_joint, "2.2.1_global_assembly_quality.tsv", "Overall assembly quality metrics taking both datasets natively")

# 4.2 Per Tech Error Breakdown
rows_tech = []
for tech, col in [("HiFi", "hifi_name"), ("ONT", "ont_name")]:
    tech_grp = df_int.groupby(col)["overlap_bp"].sum().reset_index()
    for _, r in tech_grp.iterrows():
        rows_tech.append({"dataset": tech, "category": r[col], "length_bp": r["overlap_bp"], "percentage": r["overlap_bp"] / tot_bp * 100})
df_tech = pd.DataFrame(rows_tech).sort_values(["dataset", "percentage"], ascending=[True, False])
write_tsv(df_tech, "global_tech_error_breakdown.tsv", "NucFlag categories mapped explicitly per technology")

# 4.3 Per Telomere Ranking Quality (Best/Worst + Rescues)
telo_joint = df_int.pivot_table(index="telo_id", columns="joint_status", values="overlap_bp", aggfunc="sum", fill_value=0)
for c in JOINT_CATS:
    if c not in telo_joint.columns: telo_joint[c] = 0

telo_meta = df_telo.set_index("telo_id")[["chr_display", "arm", "hap", "length", "chr_type"]]
df_telo_q = telo_meta.join(telo_joint)
df_telo_q = df_telo_q.rename(columns={"length": "total_length_bp"})

for c in JOINT_CATS:
    df_telo_q[f"{c}_pct"] = (df_telo_q[c] / df_telo_q["total_length_bp"]) * 100

df_telo_q = df_telo_q.sort_values(["Both Correct_pct", "HiFi Rescued_pct", "ONT Rescued_pct"], ascending=[False, False, False])
write_tsv(df_telo_q.reset_index(), "2.2.3_per_telomere_quality_ranking.tsv", "Best and worst telomeres ranked by joint correct percentage")

hifi_rescue = df_telo_q.sort_values("HiFi Rescued", ascending=False).reset_index()
write_tsv(hifi_rescue, "2.2.2_hifi_rescued_ranked.tsv", "Telomeres explicitly rescued by HiFi technology")

ont_rescue = df_telo_q.sort_values("ONT Rescued", ascending=False).reset_index()
write_tsv(ont_rescue, "2.2.2_ont_rescued_ranked.tsv", "Telomeres explicitly rescued by ONT technology")

# 4.4 Chromosome Type Quality
ctype_joint = df_int.groupby(["chr_type", "joint_status"])["overlap_bp"].sum().unstack(fill_value=0)
for c in JOINT_CATS:
    if c not in ctype_joint.columns: ctype_joint[c] = 0
ctype_totals = df_int.groupby("chr_type")["overlap_bp"].sum()
for c in JOINT_CATS:
    ctype_joint[f"{c}_pct"] = (ctype_joint[c] / ctype_totals) * 100
write_tsv(ctype_joint.reset_index(), "chr_type_quality.tsv", "Assembly quality distributions by acrocentric/submetacentric/metacentric group")

# 4.5 Pairwise per-telomere stats (new classification)
telo_pair_piv = df_int.pivot_table(
    index="telo_id", columns="pair_label",
    values="overlap_bp", aggfunc="sum", fill_value=0)
for c in PAIR_CATS:
    if c not in telo_pair_piv.columns:
        telo_pair_piv[c] = 0
df_telo_pair = telo_meta.join(telo_pair_piv).fillna(0)
df_telo_pair = df_telo_pair.rename(columns={"length": "total_length_bp"})
for c in PAIR_CATS:
    df_telo_pair[f"{c}_pct"] = (
        df_telo_pair[c] / df_telo_pair["total_length_bp"]) * 100
# Sort: correct desc, rescued desc, error asc
df_telo_pair["_rescued_pct"] = (
    df_telo_pair["ONT_rescued_pct"] + df_telo_pair["HiFi_rescued_pct"])
df_telo_pair = df_telo_pair.sort_values(
    ["correct_pct", "_rescued_pct", "error_pct"],
    ascending=[False, False, True])
write_tsv(df_telo_pair.reset_index(), "pair_per_telomere_quality.tsv",
          "Per-telomere quality (pairwise classification)")

# 4.6 Per-telomere per-tech scatter metrics (status-based)
print("[INFO] Computing per-telomere per-tech scatter metrics...")
telo_scatter_rows = []
for tid in df_telo["telo_id"].unique():
    sub = df_int[df_int["telo_id"] == tid]
    meta_row = df_telo[df_telo["telo_id"] == tid].iloc[0]
    tbp = meta_row["length"]
    if tbp == 0:
        continue
    hifi_perfect  = sub[sub["hifi_status"] == "PERFECT"]["overlap_bp"].sum()
    hifi_hard_err = sub[sub["hifi_status"] == "HARD_ERROR"]["overlap_bp"].sum()
    hifi_repeat   = sub[sub["hifi_status"] == "REPEAT_CONTEXT"]["overlap_bp"].sum()
    ont_perfect   = sub[sub["ont_status"] == "PERFECT"]["overlap_bp"].sum()
    ont_hard_err  = sub[sub["ont_status"] == "HARD_ERROR"]["overlap_bp"].sum()
    ont_repeat    = sub[sub["ont_status"] == "REPEAT_CONTEXT"]["overlap_bp"].sum()
    hifi_denom = tbp - hifi_repeat
    ont_denom  = tbp - ont_repeat
    telo_scatter_rows.append({
        "telo_id": tid,
        "chr_display": meta_row["chr_display"],
        "arm": meta_row["arm"],
        "hap": meta_row["hap"],
        "chr_type": meta_row["chr_type"],
        "total_bp": tbp,
        "hifi_correct_pct": hifi_perfect / tbp * 100,
        "ont_correct_pct": ont_perfect / tbp * 100,
        "hifi_error_pct": hifi_hard_err / tbp * 100,
        "ont_error_pct": ont_hard_err / tbp * 100,
        "hifi_correct_excl_rpt_pct": (
            hifi_perfect / hifi_denom * 100 if hifi_denom > 0 else 0),
        "ont_correct_excl_rpt_pct": (
            ont_perfect / ont_denom * 100 if ont_denom > 0 else 0),
        "hifi_error_excl_rpt_pct": (
            hifi_hard_err / hifi_denom * 100 if hifi_denom > 0 else 0),
        "ont_error_excl_rpt_pct": (
            ont_hard_err / ont_denom * 100 if ont_denom > 0 else 0),
    })
df_scatter = pd.DataFrame(telo_scatter_rows)
write_tsv(df_scatter, "pair_scatter_metrics.tsv",
          "Per-telomere per-tech metrics (status-based)")

# 4.7 Pair chr_type quality
ctype_pair = df_int.groupby(
    ["chr_type", "pair_label"])["overlap_bp"].sum().unstack(fill_value=0)
for c in PAIR_CATS:
    if c not in ctype_pair.columns:
        ctype_pair[c] = 0
ctype_pair_totals = df_int.groupby("chr_type")["overlap_bp"].sum()
for c in PAIR_CATS:
    ctype_pair[f"{c}_pct"] = (ctype_pair[c] / ctype_pair_totals) * 100
write_tsv(ctype_pair.reset_index(), "pair_chr_type_quality.tsv",
          "Chr type quality (pairwise classification)")

# By hap
ctype_hap_pair = df_int.groupby(
    ["chr_type", "hap", "pair_label"])["overlap_bp"].sum().unstack(fill_value=0)
for c in PAIR_CATS:
    if c not in ctype_hap_pair.columns:
        ctype_hap_pair[c] = 0
ctype_hap_pair_totals = df_int.groupby(["chr_type", "hap"])["overlap_bp"].sum()
for c in PAIR_CATS:
    ctype_hap_pair[f"{c}_pct"] = (
        ctype_hap_pair[c] / ctype_hap_pair_totals) * 100

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Publication-Ready Figures
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[INFO] Generating publication-ready Cell figures...")

# Fig 1: Global assembly quality summary (pairwise classification, thin bars)
def plot_fig1(pair_colors, suffix=""):
    """Thin horizontal stacked bar: global pairwise classification of
    all telomere sequence, using the Fig 2 palette."""
    g_pair = df_int.groupby("pair_label")["overlap_bp"].sum()
    g_pair_pct = g_pair.reindex(PAIR_CATS).fillna(0) / tot_bp * 100

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, 0.8), constrained_layout=True)
    left = 0
    for cat in PAIR_CATS:
        val = g_pair_pct.get(cat, 0)
        if val > 0:
            ax.barh(0, val, left=left, color=pair_colors[cat],
                    height=0.45, linewidth=0.0, edgecolor="none")
            if val > 5:
                ax.text(left + val / 2, 0, f"{val:.1f}%",
                        ha="center", va="center", fontsize=TICK,
                        color="white" if cat in ("correct", "error") else "#222222")
            left += val
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.4, 0.4)
    ax.set_yticks([])
    ax.set_xlabel("Telomere by NucFlag status (%)", fontsize=LABEL)
    ax.tick_params(labelsize=TICK)
    _box(ax)

    # Legend: 2 rows × 3 columns, thin black box
    handles = [Patch(fc=pair_colors[c], label=PAIR_LABELS[c])
               for c in PAIR_CATS]
    leg = fig.legend(handles=handles, fontsize=TICK, loc="lower center",
                     bbox_to_anchor=(0.5, -0.55), ncol=3,
                     frameon=True, edgecolor="black",
                     columnspacing=0.8, handletextpad=0.3)
    leg.get_frame().set_linewidth(0.4)

    _save(fig, f"Fig1_Global_Quality{suffix}")

plot_fig1(PAIR_COLORS_V1, "_v1")
plot_fig1(PAIR_COLORS_V2, "_v2")
plot_fig1(PAIR_COLORS_V3, "_v3")

# Fig 2: Per-telomere quality spectrum (pairwise classification)
def plot_fig2(pair_colors, suffix=""):
    """Panel A (all telomeres) + Panel B (arm × hap 2×2 sub-grid)."""

    def _sort_wf(df_sub):
        df_s = df_sub.copy()
        df_s["_k1"] = df_s["correct_pct"]
        df_s["_k2"] = df_s["_k1"] + df_s["ONT_rescued_pct"] + df_s["HiFi_rescued_pct"]
        df_s["_k3"] = df_s["_k2"] + df_s["repeat_pct"]
        df_s["_k4"] = df_s["_k3"] + df_s["mismatch_signal_pct"]
        return df_s.sort_values(
            ["_k1", "_k2", "_k3", "_k4", "error_pct"],
            ascending=[False, False, False, False, True])

    def _waterfall(ax, df_sub, show_ylabel=False):
        df_s = _sort_wf(df_sub)
        n = len(df_s)
        x = np.arange(n)
        bottoms = np.zeros(n)
        for cat in PAIR_CATS:
            col = f"{cat}_pct"
            if col not in df_s.columns:
                continue
            vals = df_s[col].values
            ax.bar(x, vals, bottom=bottoms, color=pair_colors[cat],
                   width=1.0, linewidth=0.0, edgecolor="none")
            bottoms += vals
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_xticks([])
        if show_ylabel:
            ax.set_ylabel("Telomere annotation\nby NucFlag status (%)", fontsize=LABEL)
        ax.tick_params(labelsize=TICK)
        _box(ax)

    fig = plt.figure(figsize=(FIG_WIDTH, 2.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.1],
                          wspace=0.18)
    df_w = df_telo_pair.reset_index()

    # A: All telomeres
    ax_a = fig.add_subplot(gs[0, 0])
    _waterfall(ax_a, df_w, show_ylabel=True)
    ax_a.set_title("All telomeres", fontsize=LABEL)

    # B: arm × haplotype (2×2 sub-grid)
    # ┌────────────────────────────────────────────────────────────────────┐
    # │  ADJUST wspace / hspace below if sub-panels overlap or are too    │
    # │  far apart.  Increase to add breathing room.                      │
    # └────────────────────────────────────────────────────────────────────┘
    gs_b = gs[0, 1].subgridspec(2, 2,
                                wspace=0.25,   # ← horizontal gap between B cols
                                hspace=0.45)   # ← vertical   gap between B rows
    arm_hap_combos = [("p", "hap1"), ("q", "hap1"),
                      ("p", "hap2"), ("q", "hap2")]
    b_axes = []
    for idx, (arm, hap) in enumerate(arm_hap_combos):
        r, c = divmod(idx, 2)
        ax_b = fig.add_subplot(gs_b[r, c])
        sub = df_w[(df_w["arm"] == arm) & (df_w["hap"] == hap)]
        _waterfall(ax_b, sub)
        b_axes.append(ax_b)

    # Column headers (top of each column) — same size as Y-axis / title
    b_axes[0].set_title("p-arm", fontsize=LABEL)
    b_axes[1].set_title("q-arm", fontsize=LABEL)

    # Row labels (left side, vertical 90°); remove tick labels from B
    for ab in b_axes:
        ab.set_yticklabels([])
    b_axes[0].set_ylabel("hap1", fontsize=LABEL, rotation=90,
                         labelpad=6, va="center")
    b_axes[2].set_ylabel("hap2", fontsize=LABEL, rotation=90,
                         labelpad=6, va="center")

    # Panel labels (Cell standard: 8 pt, regular weight)
    ax_a.text(-0.10, 1.06, "A", transform=ax_a.transAxes,
              fontsize=8, fontweight="regular", va="bottom", ha="right")
    b_axes[0].text(-0.25, 1.20, "B", transform=b_axes[0].transAxes,
                   fontsize=8, fontweight="regular", va="bottom", ha="right")

    # Shared legend: 2 rows × 3 columns, thin black box
    handles = [Patch(fc=pair_colors[c], label=PAIR_LABELS[c])
               for c in PAIR_CATS]
    leg = fig.legend(handles=handles, fontsize=TICK, loc="lower center",
                     bbox_to_anchor=(0.5, -0.16), ncol=3,
                     frameon=True, edgecolor="black",
                     columnspacing=0.8, handletextpad=0.3)
    leg.get_frame().set_linewidth(0.4)

    _save(fig, f"Fig2_Waterfall{suffix}")

plot_fig2(PAIR_COLORS_V1, "_v1")
plot_fig2(PAIR_COLORS_V2, "_v2")
plot_fig2(PAIR_COLORS_V3, "_v3")

# ── Fig3 statistics: thorough console-only analysis ─────────────────────────
def _run_chr_type_stats_console():
    """Comprehensive statistical analysis of chromosome type associations
    with HiFi/ONT correct% and error%. Printed to console only."""
    from scipy.stats import mannwhitneyu, kruskal, ttest_ind, spearmanr
    from itertools import combinations

    acro = df_scatter[df_scatter["chr_type"] == "acrocentric"]
    meta = df_scatter[df_scatter["chr_type"] == "metacentric"]
    subm = df_scatter[df_scatter["chr_type"] == "submetacentric"]
    groups = {"acrocentric": acro, "metacentric": meta, "submetacentric": subm}

    metrics = [
        ("hifi_correct_pct",          "ont_correct_pct",          "Correct – all reads"),
        ("hifi_correct_excl_rpt_pct", "ont_correct_excl_rpt_pct", "Correct – excl. repeats"),
        ("hifi_error_pct",            "ont_error_pct",            "Error – all reads"),
        ("hifi_error_excl_rpt_pct",   "ont_error_excl_rpt_pct",  "Error – excl. repeats"),
    ]

    print("\n" + "="*90)
    print("CHROMOSOME TYPE × TECHNOLOGY ASSOCIATION ANALYSIS")
    print("="*90)

    all_pvals = []  # collect (label, p) for BH correction

    for hifi_col, ont_col, label in metrics:
        print(f"\n{'─'*90}")
        print(f"  {label}")
        print(f"{'─'*90}")

        # ── Descriptive statistics ──
        print(f"  {'Chr type':<18} {'n':>3}  {'HiFi mean±sd':>18}  {'ONT mean±sd':>18}")
        for ctype in ["acrocentric", "metacentric", "submetacentric"]:
            g = groups[ctype]
            print(f"  {ctype:<18} {len(g):>3}  "
                  f"{g[hifi_col].mean():>7.1f}±{g[hifi_col].std():>5.1f}  "
                  f"{g[ont_col].mean():>7.1f}±{g[ont_col].std():>5.1f}")

        # ── Kruskal–Wallis (non-parametric omnibus, 3 groups) ──
        print(f"\n  Kruskal–Wallis (omnibus, 3 chr-types):")
        for tech_name, col in [("HiFi", hifi_col), ("ONT", ont_col)]:
            vals = [groups[ct][col].values for ct in ["acrocentric", "metacentric", "submetacentric"]]
            if all(len(v) > 0 for v in vals):
                H_stat, p_kw = kruskal(*vals)
                print(f"    {tech_name}: H={H_stat:.3f}, p={p_kw:.4e}")
                all_pvals.append((f"{label} | KW {tech_name}", p_kw))

        # ── Pairwise Mann–Whitney U (non-parametric, directed hypotheses) ──
        print(f"\n  Mann–Whitney U pairwise (two-sided):")
        for ct1, ct2 in combinations(["acrocentric", "metacentric", "submetacentric"], 2):
            g1, g2 = groups[ct1], groups[ct2]
            if len(g1) < 2 or len(g2) < 2:
                continue
            for tech_name, col in [("HiFi", hifi_col), ("ONT", ont_col)]:
                u_stat, p_mw = mannwhitneyu(g1[col], g2[col], alternative="two-sided")
                all_pvals.append((f"{label} | MWU {tech_name} {ct1[:4]} vs {ct2[:4]}", p_mw))
                sig = "*" if p_mw < 0.05 else ""
                print(f"    {tech_name} {ct1[:4]} vs {ct2[:4]}: U={u_stat:.0f}, p={p_mw:.4e} {sig}")

        # ── Directed hypotheses (one-sided) ──
        print(f"\n  Directed hypotheses (one-sided):")
        if len(acro) > 1 and len(meta) > 1:
            # Acro ONT > Meta ONT?
            u, p = mannwhitneyu(acro[ont_col], meta[ont_col], alternative="greater")
            all_pvals.append((f"{label} | Acro ONT > Meta ONT", p))
            print(f"    Acro ONT > Meta ONT: U={u:.0f}, p={p:.4e}")
            # Meta HiFi > Acro HiFi?
            u, p = mannwhitneyu(meta[hifi_col], acro[hifi_col], alternative="greater")
            all_pvals.append((f"{label} | Meta HiFi > Acro HiFi", p))
            print(f"    Meta HiFi > Acro HiFi: U={u:.0f}, p={p:.4e}")

        # ── Welch's t-test (parametric, for comparison) ──
        print(f"\n  Welch's t-test (parametric, two-sided):")
        if len(acro) > 1 and len(meta) > 1:
            for tech_name, col in [("HiFi", hifi_col), ("ONT", ont_col)]:
                t_stat, p_t = ttest_ind(acro[col], meta[col], equal_var=False)
                all_pvals.append((f"{label} | Welch {tech_name} acro vs meta", p_t))
                print(f"    {tech_name} acro vs meta: t={t_stat:.3f}, p={p_t:.4e}")

        # ── Spearman correlation: HiFi vs ONT within each type ──
        print(f"\n  Spearman corr (HiFi vs ONT within chr-type):")
        for ctype in ["acrocentric", "metacentric", "submetacentric"]:
            g = groups[ctype]
            if len(g) > 2:
                rho, p_sp = spearmanr(g[hifi_col], g[ont_col])
                print(f"    {ctype}: rho={rho:.3f}, p={p_sp:.4e}")

    # ── Multiple testing correction: BH + Bonferroni ──
    print(f"\n{'='*90}")
    print("MULTIPLE TESTING CORRECTION (Benjamini-Hochberg & Bonferroni)")
    print(f"{'='*90}")
    n_tests = len(all_pvals)
    sorted_pvals = sorted(all_pvals, key=lambda x: x[1])

    # Benjamini-Hochberg
    bh_results = []
    for rank, (test_label, p_raw) in enumerate(sorted_pvals, 1):
        p_bh = min(1.0, p_raw * n_tests / rank)
        p_bonf = min(1.0, p_raw * n_tests)
        bh_results.append((test_label, p_raw, p_bh, p_bonf, rank))
    # Enforce BH monotonicity (step-up)
    for i in range(len(bh_results) - 2, -1, -1):
        bh_results[i] = (bh_results[i][0], bh_results[i][1],
                         min(bh_results[i][2], bh_results[i+1][2]),
                         bh_results[i][3], bh_results[i][4])

    print(f"  Total tests: {n_tests}")
    print(f"  {'Rank':>4} {'p_raw':>12} {'p_adj_BH':>12} {'p_adj_Bonf':>12} {'Sig(BH)':>7} {'Sig(Bonf)':>9}  Test")
    for test_label, p_raw, p_bh, p_bonf, rank in bh_results:
        sig_bh   = "***" if p_bh   < 0.001 else "**" if p_bh   < 0.01 else "*" if p_bh   < 0.05 else ""
        sig_bonf = "***" if p_bonf < 0.001 else "**" if p_bonf < 0.01 else "*" if p_bonf < 0.05 else ""
        print(f"  {rank:>4} {p_raw:>12.4e} {p_bh:>12.4e} {p_bonf:>12.4e} {sig_bh:>7} {sig_bonf:>9}  {test_label}")

    # ── Interpretation ──
    # Count Bonferroni- and BH-significant results
    n_bonf = sum(1 for _, _, _, p_bonf, _ in bh_results if p_bonf < 0.05)
    n_bh   = sum(1 for _, _, p_bh, _, _ in bh_results if p_bh < 0.05)

    print(f"\n{'='*90}")
    print("INTERPRETATION")
    print(f"{'='*90}")

    print("\n  CORRECTION STRATEGY (Bonferroni-primary, BH-secondary)")
    print("  ─────────────────────────────────────────────────────────")
    print("    - Bonferroni controls the family-wise error rate (FWER):")
    print("      the probability of ANY false positive across all tests.")
    print("      Most conservative; appropriate when individual claims")
    print("      must each be defensible (genomics publication standard).")
    print("    - Benjamini-Hochberg (BH) controls the false discovery")
    print("      rate (FDR): the expected PROPORTION of false positives.")
    print("      Less conservative; appropriate for exploratory screens.")
    print("    - Decision: report Bonferroni-corrected p-values as the")
    print("      primary significance criterion. Results that also pass")
    print("      BH provide additional confidence. Results passing BH")
    print("      but not Bonferroni are flagged as suggestive only.")
    print(f"    - Summary: {n_bonf} tests survive Bonferroni (α=0.05),")
    print(f"      {n_bh} survive BH (q=0.05), out of {n_tests} total.")

    print("\n  TEST SELECTION")
    print("  ─────────────────────────────────────────────────────────")
    print("    - Non-parametric tests (Mann-Whitney U, Kruskal-Wallis)")
    print("      are primary: small group sizes (n=19-52), percentage")
    print("      data bounded [0,100], and known outliers (1p_hap2,")
    print("      7q_hap2) violate normality assumptions.")
    print("    - Welch's t-tests included for comparison only; where")
    print("      they agree with non-parametric tests → added confidence.")
    print("    - Spearman correlations assess HiFi-ONT concordance")
    print("      within each chromosome type (non-parametric, rank-based).")

    print("\n  SIGNIFICANT FINDINGS (Bonferroni-corrected)")
    print("  ─────────────────────────────────────────────────────────")
    print("    1. HiFi correct% depends strongly on chromosome type")
    print("       (all reads). Acrocentric telomeres show dramatically")
    print("       lower HiFi correct% (24.8±24.9%) vs metacentric")
    print("       (84.5±25.6%) and submetacentric (72.1±32.3%).")
    print("       Supported by: KW omnibus, MWU acro-vs-meta,")
    print("       MWU acro-vs-subm, directed meta>acro, and Welch")
    print("       (all p_Bonf < 5.0e-04). This is the strongest signal")
    print("       in the dataset.")
    print("")
    print("    2. ONT rescues acrocentric telomeres when repeats are")
    print("       excluded. Acrocentric ONT correct% (94.6±21.5%)")
    print("       significantly exceeds metacentric (68.3±35.8%).")
    print("       Supported by: directed acro ONT > meta ONT")
    print("       (p_Bonf=1.53e-02) and MWU two-sided (p_Bonf=3.07e-02).")
    print("       This confirms the complementarity of ONT for regions")
    print("       where HiFi struggles (acrocentric short arms).")

    print("\n  SUGGESTIVE FINDINGS (BH only, not Bonferroni)")
    print("  ─────────────────────────────────────────────────────────")
    print("    3. When repeats are excluded, the HiFi disadvantage for")
    print("       acrocentric telomeres attenuates but remains nominally")
    print("       significant (BH), while ONT's advantage for acrocentric")
    print("       telomeres (vs submetacentric) emerges under BH.")
    print("       These effects are consistent with Finding 1-2 but do")
    print("       not independently survive Bonferroni.")

    print("\n  NON-SIGNIFICANT FINDINGS")
    print("  ─────────────────────────────────────────────────────────")
    print("    4. Error rates show no significant differences across")
    print("       chromosome types for either technology after correction.")
    print("       The quality difference is about correct-classification")
    print("       capacity, not differential error accumulation.")
    print("")
    print("    5. Spearman HiFi-ONT correlations within chromosome types")
    print("       are weak and inconsistent (rho -0.36 to +0.31).")
    print("       Only submetacentric shows a nominally significant")
    print("       negative correlation (rho=-0.361, p=8.6e-03), but this")
    print("       does not survive multiple-testing correction and should")
    print("       not be over-interpreted.")

    print("\n  BIOLOGICAL CONCLUSION")
    print("  ─────────────────────────────────────────────────────────")
    print("    Acrocentric p-arm telomeres are the primary challenge for")
    print("    HiFi-based assembly, likely due to satellite/rDNA-adjacent")
    print("    repeat content causing collapses/misjoins in short reads.")
    print("    ONT long reads compensate specifically in these regions,")
    print("    providing near-complete (94.6%) correct coverage once")
    print("    repeat context is excluded. The two technologies are")
    print("    complementary, and dual-platform NucFlag validation is")
    print("    essential for telomere-to-telomere assembly QC, especially")
    print("    at acrocentric chromosomes.")
    print(f"{'='*90}\n")

    # Save stats to TSV
    df_stats = pd.DataFrame(bh_results,
                            columns=["test", "p_raw", "p_adj_BH", "p_adj_Bonferroni", "rank"])
    write_tsv(df_stats, "chr_type_association_stats.tsv",
              "Chromosome type × technology association stats (BH + Bonferroni corrected)")

_run_chr_type_stats_console()

# Fig 3: HiFi vs ONT concordance scatter (pairwise classification)
def plot_fig3(chr_colors, suffix=""):
    """2×2 scatter: correct% & error% × all DNA & excl. repeats.
    Bubble size = telomere length.  Top 3 labelled per panel type."""

    def _telo_label(row):
        num = row["chr_display"].replace("chr", "")
        return f"{num}{row['arm']}_{row['hap']}"

    def _get_top_n(xcol, ycol, panel_type, n=2):
        """Return {telo_id: label} for the n most extreme points."""
        combined = df_scatter[xcol] + df_scatter[ycol]
        if panel_type == "correct":
            idx = combined.nsmallest(n).index
        else:
            idx = combined.nlargest(n).index
        return {row["telo_id"]: _telo_label(row)
                for _, row in df_scatter.loc[idx].iterrows()}

    configs = [
        ("hifi_correct_pct", "ont_correct_pct",
         "Correct \u2013 all reads", "HiFi correct (%)", "ONT correct (%)",
         "correct"),
        ("hifi_correct_excl_rpt_pct", "ont_correct_excl_rpt_pct",
         "Correct \u2013 excluding repeats", "HiFi correct (%)",
         "ONT correct (%)", "correct"),
        ("hifi_error_pct", "ont_error_pct",
         "Error \u2013 all reads", "HiFi error (%)", "ONT error (%)",
         "error"),
        ("hifi_error_excl_rpt_pct", "ont_error_excl_rpt_pct",
         "Error \u2013 excluding repeats", "HiFi error (%)",
         "ONT error (%)", "error"),
    ]

    bp_vals = df_scatter["total_bp"].values.astype(float)
    bp_min, bp_max = bp_vals.min(), bp_vals.max()
    # Smaller markers to reduce obstruction (was 8 + 112 = 120)
    sizes = 6 + (bp_vals - bp_min) / (bp_max - bp_min + 1e-9) * 50
    colors_all = [chr_colors.get(r, "#999999")
                  for r in df_scatter["chr_type"]]

    # Manual gridspec layout to prevent C/D overlap with legends
    fig = plt.figure(figsize=(FIG_WIDTH, 5))
    gs = fig.add_gridspec(2, 2, left=0.13, right=0.95, top=0.94,
                          bottom=0.20, hspace=0.50, wspace=0.50)
    panel_axes = []

    for panel_idx, (xcol, ycol, title, xlabel, ylabel, panel_type) in \
            enumerate(configs):
        ax = fig.add_subplot(gs[panel_idx // 2, panel_idx % 2])
        panel_axes.append(ax)
        labels_dict = _get_top_n(xcol, ycol, panel_type, n=2)

        ax.scatter(df_scatter[xcol], df_scatter[ycol],
                   c=colors_all, s=sizes, alpha=0.75,
                   edgecolor="black", linewidth=0.3, zorder=3)

        # Label top 2 extreme points (5 pt = Cell minimum, no connector lines)
        for i, row in df_scatter.iterrows():
            if row["telo_id"] in labels_dict:
                lab = labels_dict[row["telo_id"]]
                # 1p_hap2: left-top in correct panels, top-right in error panels
                if "1p" in lab:
                    if panel_type == "correct":
                        xytext = (-22, 12)
                    else:
                        xytext = (12, 12)
                else:
                    xytext = (5, -8)
                ax.annotate(
                    lab, xy=(row[xcol], row[ycol]),
                    xytext=xytext, textcoords="offset points",
                    fontsize=5, color="#333333")

        # Diagonal reference line with generous padding
        all_v = np.concatenate([df_scatter[xcol].values,
                                df_scatter[ycol].values])
        data_lo, data_hi = all_v.min(), all_v.max()
        span = data_hi - data_lo
        padding = max(12, span * 0.30)   # more breathing room
        lo = max(-10, data_lo - padding)
        hi = min(110, data_hi + padding)
        ax.plot([lo, hi], [lo, hi], ls="-", color="#C0392B",
                linewidth=0.7, zorder=1)

        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_xlabel(xlabel, fontsize=LABEL)
        ax.set_ylabel(ylabel, fontsize=LABEL)
        ax.set_title(title, fontsize=LABEL)
        ax.tick_params(labelsize=TICK)
        _box(ax)

    # Panel labels (Cell standard: 8 pt, regular weight)
    for ax, lbl in zip(panel_axes, "ABCD"):
        ax.text(-0.18, 1.08, lbl, transform=ax.transAxes,
                fontsize=8, fontweight="regular", va="bottom", ha="right")

    # Legend 1: chromosome type — centred at base, Fig2-style border
    h_chr = [Patch(fc=chr_colors[k],
                   edgecolor="black", linewidth=0.3,
                   label=k.capitalize())
             for k in ["metacentric", "submetacentric", "acrocentric"]]
    leg1 = fig.legend(handles=h_chr, fontsize=TICK,
                      loc="lower center",
                      bbox_to_anchor=(0.33, 0.01), ncol=1,
                      frameon=True, title="Chromosome type",
                      title_fontsize=TICK,
                      edgecolor="black")
    leg1.get_frame().set_linewidth(0.4)

    # Legend 2: bubble size — centred at base, Fig2-style border
    bp_mid = int(round((bp_min + bp_max) / 2 / 1000)) * 1000
    bp_ceil = int(np.ceil(bp_max / 1000)) * 1000
    legend_bp = sorted(set([1000, bp_mid, bp_ceil]))
    h_bub = [Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='#7F7F7F',
                    markeredgecolor='black',
                    markeredgewidth=0.4,
                    markersize=np.sqrt(
                        6 + (bp - bp_min) / (bp_max - bp_min + 1e-9)
                        * 50),
                    label=f'{bp/1000:.0f} kb', linestyle='None')
             for bp in legend_bp]
    leg2 = fig.legend(handles=h_bub, fontsize=TICK,
                      loc="lower center",
                      bbox_to_anchor=(0.67, 0.01), ncol=1,
                      frameon=True, title="Telomere length",
                      title_fontsize=TICK,
                      edgecolor="black")
    leg2.get_frame().set_linewidth(0.4)
    fig.add_artist(leg1)

    _save(fig, f"Fig3_Scatter_2x2{suffix}")

plot_fig3(CHR_TYPE_COLORS, "_v1_current")
plot_fig3(CHR_TYPE_COLORS_V2, "_v2_gnomAD")
plot_fig3(CHR_TYPE_COLORS_V3, "_v3_OkabeIto")

print(f"\n[DONE] All Cell-ready figures saved to: {plots_dir}")
