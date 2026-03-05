#!/usr/bin/env python3
"""H9 telomere paired NucFlag analysis using combined HiFi+ONT segments."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parents[1]
DATA_DIR = BASE_DIR / "2_data" / "2.2_processed"
NUCFLAG_DIR = DATA_DIR / "nucflag"
TELOSCOPE_BED = DATA_DIR / "25.12.10_teloscope_compiled" / "H9_T2T_v0.1_dip.fasta_terminal_telomeres.bed"
COMBINED_TSV = NUCFLAG_DIR / "nucflag_telo_combined.tsv"
PALETTE_TSV = NUCFLAG_DIR / "nucflag_palette.tsv"

RUN_LABEL = "v2"
PLOTS_DIR = BASE_DIR / "3_figures" / "3.1_draft" / "26.03.01_telomeres_nucflag_combined_codex"
STATS_DIR = PLOTS_DIR / f"stats_{RUN_LABEL}"
FIGDATA_DIR = PLOTS_DIR / f"figure_data_{RUN_LABEL}"
COMBINED_STATS_OUT = PLOTS_DIR / f"h9_nucflag_telomere_stats_{RUN_LABEL}.tsv"
INTERPRETATION_TXT = STATS_DIR / f"interpretation_summary_{RUN_LABEL}.txt"

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
EXTS = ("png", "pdf", "svg")

FALLBACK_COLORS: Dict[str, str] = {
    "correct": "#2D8E47", "misjoin": "#CB3335", "collapse": "#7B5EA7", "false_dup": "#8B6D52",
    "het_or_mismap": "#C77CBA", "low_quality": "#E69F00", "deletion": "#3C7DC4", "insertion": "#00B6AD",
    "mismatch": "#D4A535", "softclip": "#5A6F80", "homopolymer": "#7CB950", "dinucleotide": "#B5CC3E",
    "simple_repeat": "#E07B42", "other_repeat": "#47A697", "scaffold": "#999999", "uncovered": "#D4D4D4",
}


def init_plot_style() -> None:
    font_family = "Arial"
    try:
        from matplotlib.font_manager import FontProperties, findfont
        if findfont(FontProperties(family="Arial")) == findfont(FontProperties()):
            font_family = "DejaVu Sans"
    except Exception:
        font_family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": font_family, "font.size": 7, "axes.linewidth": 0.4,
        "xtick.major.width": 0.4, "ytick.major.width": 0.4,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def chr_key(chrom: str) -> int:
    return CHR_ORDER.index(chrom) if chrom in CHR_ORDER else 100


def pct(part: float, whole: float) -> float:
    return float(part) / float(whole) if whole else 0.0


def fmt_bp(value: float) -> str:
    return f"{int(round(value)):,}"


def rgb_to_hex(rgb_text: str) -> str:
    raw = str(rgb_text).strip().replace('"', "")
    r, g, b = [int(x) for x in raw.split(",")]
    return f"#{r:02X}{g:02X}{b:02X}"


def save_figure(fig: plt.Figure, stem: str) -> None:
    for ext in EXTS:
        fig.savefig(PLOTS_DIR / f"{stem}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Figure: {stem}")


def write_tsv(df: pd.DataFrame, path: Path, note: str) -> None:
    df.to_csv(path, sep="\t", index=False)
    print(f"[OK] Table: {path.name} ({len(df)} rows) - {note}")


def load_palette(path: Path) -> Tuple[pd.DataFrame, Dict[str, str]]:
    pal = pd.read_csv(path, sep="\t")
    if not {"name", "itemRgb"}.issubset(pal.columns):
        raise ValueError(f"Palette missing columns: {path}")
    pal = pal.copy()
    pal["hex"] = pal["itemRgb"].map(rgb_to_hex)
    color_map = dict(FALLBACK_COLORS)
    color_map.update(dict(zip(pal["name"], pal["hex"])))
    color_map.setdefault("uncovered", "#D4D4D4")
    return pal, color_map


def load_teloscope(path: Path) -> pd.DataFrame:
    cols = ["chrom", "start", "end", "length", "arm", "fwdCounts", "revCounts", "canCounts", "nonCanCounts", "chr_length"]
    telo = pd.read_csv(path, sep="\t", header=None, names=cols)
    for c in ("start", "end", "length", "chr_length"):
        telo[c] = pd.to_numeric(telo[c], errors="coerce")
    telo = telo.loc[telo["chrom"].notna()].copy()
    telo["hap"] = telo["chrom"].str.extract(r"_(hap[12])$")[0]
    telo["chr_display"] = telo["chrom"].str.replace(r"_hap[12]$", "", regex=True)
    telo["telomere_id"] = telo["chrom"] + ":" + telo["arm"]
    telo = telo.sort_values(["chr_display", "hap", "arm"], key=lambda s: s.map(chr_key) if s.name == "chr_display" else s)
    print(f"[INFO] Teloscope: {len(telo)} rows, {telo['chrom'].nunique()} chromosomes, total={fmt_bp(telo['length'].sum())} bp")
    return telo


def load_combined(path: Path) -> pd.DataFrame:
    comb = pd.read_csv(path, sep="\t")
    req = {"chrom", "chromStart", "chromEnd", "length", "hifi_name", "ont_name", "concordant"}
    miss = req - set(comb.columns)
    if miss:
        raise ValueError(f"Combined TSV missing columns: {sorted(miss)}")
    for c in ("chromStart", "chromEnd", "length"):
        comb[c] = pd.to_numeric(comb[c], errors="coerce")
    for c in ("hifi_name", "ont_name"):
        comb[c] = comb[c].fillna("uncovered").astype(str)
    comb["concordant"] = comb["concordant"].astype(str).str.lower().map({"true": True, "false": False}).fillna(False)
    for c in ("hifi_score", "hifi_zscore", "hifi_af", "ont_score", "ont_zscore", "ont_af"):
        if c in comb.columns:
            comb[c] = pd.to_numeric(comb[c], errors="coerce")
    print(f"[INFO] Combined TSV: {len(comb)} rows, {comb['chrom'].nunique()} chromosomes")
    return comb

def intersect_telomeres(comb: pd.DataFrame, telo: pd.DataFrame) -> pd.DataFrame:
    rows: List[dict] = []
    comb_by_chr = {k: g for k, g in comb.groupby("chrom", sort=False)}
    for _, t in telo.iterrows():
        chrom = t["chrom"]
        bs, be = int(t["start"]), int(t["end"])
        csub = comb_by_chr.get(chrom)
        if csub is None or csub.empty:
            continue
        cand = csub.loc[(csub["chromEnd"] > bs) & (csub["chromStart"] < be)]
        if cand.empty:
            continue
        for _, c in cand.iterrows():
            s = max(bs, int(c["chromStart"]))
            e = min(be, int(c["chromEnd"]))
            ov = e - s
            if ov <= 0:
                continue
            rows.append({
                "chrom": chrom, "chr_display": t["chr_display"], "hap": t["hap"], "arm": t["arm"],
                "telomere_id": t["telomere_id"], "telomere_start": bs, "telomere_end": be,
                "seg_start": s, "seg_end": e, "overlap_bp": int(ov),
                "hifi_name": c["hifi_name"], "ont_name": c["ont_name"], "concordant": bool(c["concordant"]),
                "hifi_score": c.get("hifi_score", np.nan), "hifi_zscore": c.get("hifi_zscore", np.nan), "hifi_af": c.get("hifi_af", np.nan),
                "ont_score": c.get("ont_score", np.nan), "ont_zscore": c.get("ont_zscore", np.nan), "ont_af": c.get("ont_af", np.nan),
            })
    ov = pd.DataFrame(rows)
    if ov.empty:
        raise RuntimeError("No overlap segments found")
    exp_bp = int(telo["length"].sum())
    obs_bp = int(ov["overlap_bp"].sum())
    if exp_bp != obs_bp:
        raise RuntimeError(f"Overlap mismatch: expected {exp_bp}, observed {obs_bp}")
    print(f"[INFO] Overlap segments: {len(ov)} rows, total={fmt_bp(obs_bp)} bp")
    return ov


def ordered_categories(overlap_df: pd.DataFrame, palette_df: pd.DataFrame) -> List[str]:
    observed = set(overlap_df["hifi_name"]).union(set(overlap_df["ont_name"]))
    pal_order = [x for x in palette_df["name"].tolist() if x in observed]
    cats: List[str] = []
    if "correct" in observed:
        cats.append("correct")
    for cat in pal_order:
        if cat not in cats:
            cats.append(cat)
    cats.extend(sorted(observed - set(cats) - {"uncovered"}))
    if "uncovered" in observed:
        cats.append("uncovered")
    return cats


def build_tables(overlap_df: pd.DataFrame, telo_df: pd.DataFrame, cats: List[str]) -> Dict[str, pd.DataFrame]:
    total_bp = int(telo_df["length"].sum())

    both_correct = (overlap_df["hifi_name"] == "correct") & (overlap_df["ont_name"] == "correct")
    one_correct = (overlap_df["hifi_name"] == "correct") ^ (overlap_df["ont_name"] == "correct")
    both_non = (overlap_df["hifi_name"] != "correct") & (overlap_df["ont_name"] != "correct")

    global_q = pd.DataFrame([
        {"metric": "both_correct", "bp": int(overlap_df.loc[both_correct, "overlap_bp"].sum()), "description": "HiFi and ONT both correct"},
        {"metric": "one_correct_rescued", "bp": int(overlap_df.loc[one_correct, "overlap_bp"].sum()), "description": "Only one technology correct"},
        {"metric": "both_non_correct", "bp": int(overlap_df.loc[both_non, "overlap_bp"].sum()), "description": "Neither technology correct"},
        {"metric": "concordant", "bp": int(overlap_df.loc[overlap_df["concordant"], "overlap_bp"].sum()), "description": "HiFi and ONT agree"},
        {"metric": "discordant", "bp": int(overlap_df.loc[~overlap_df["concordant"], "overlap_bp"].sum()), "description": "HiFi and ONT disagree"},
    ])
    global_q["pct_teloscope"] = global_q["bp"].map(lambda x: pct(x, total_bp))

    tech_rows: List[dict] = []
    for tech, col in (("HiFi", "hifi_name"), ("ONT", "ont_name")):
        g = overlap_df.groupby(col)["overlap_bp"].sum()
        for cat in cats:
            bp_v = int(g.get(cat, 0))
            tech_rows.append({"technology": tech, "category": cat, "bp": bp_v, "pct_teloscope": pct(bp_v, total_bp)})
    tech_cat = pd.DataFrame(tech_rows)

    pair_long = overlap_df.groupby(["hifi_name", "ont_name"], as_index=False)["overlap_bp"].sum().rename(columns={"overlap_bp": "bp"})
    pair_long["pct_teloscope"] = pair_long["bp"].map(lambda x: pct(x, total_bp))

    pair_wide = pair_long.pivot_table(index="hifi_name", columns="ont_name", values="bp", fill_value=0)
    pair_wide = pair_wide.reindex(index=cats, columns=cats, fill_value=0).reset_index().rename(columns={"hifi_name": "hifi_category"})

    def per_tel(g: pd.DataFrame) -> pd.Series:
        total = int(g["overlap_bp"].sum())
        both = int(g.loc[(g["hifi_name"] == "correct") & (g["ont_name"] == "correct"), "overlap_bp"].sum())
        one = int(g.loc[(g["hifi_name"] == "correct") ^ (g["ont_name"] == "correct"), "overlap_bp"].sum())
        both_non_v = int(g.loc[(g["hifi_name"] != "correct") & (g["ont_name"] != "correct"), "overlap_bp"].sum())
        hifi_ok = int(g.loc[g["hifi_name"] == "correct", "overlap_bp"].sum())
        ont_ok = int(g.loc[g["ont_name"] == "correct", "overlap_bp"].sum())
        conc = int(g.loc[g["concordant"], "overlap_bp"].sum())
        disc = int(g.loc[~g["concordant"], "overlap_bp"].sum())
        hifi_rescue = hifi_ok - both
        ont_rescue = ont_ok - both
        q = np.select([
            (g["hifi_name"] == "correct") & (g["ont_name"] == "correct"),
            (g["hifi_name"] == "correct") ^ (g["ont_name"] == "correct")
        ], [1, 0], default=-1)
        qscore = float(np.sum(q * g["overlap_bp"]) / total) if total else 0.0
        win = "Tie"
        if hifi_ok > ont_ok:
            win = "HiFi"
        elif ont_ok > hifi_ok:
            win = "ONT"
        pair_mode = g.groupby(["hifi_name", "ont_name"])["overlap_bp"].sum().sort_values(ascending=False)
        dom_pair = f"{pair_mode.index[0][0]}|{pair_mode.index[0][1]}" if not pair_mode.empty else ""
        return pd.Series({
            "chrom": g["chrom"].iat[0], "chr_display": g["chr_display"].iat[0], "hap": g["hap"].iat[0], "arm": g["arm"].iat[0], "telomere_id": g["telomere_id"].iat[0],
            "total_bp": total, "pct_teloscope": pct(total, total_bp),
            "hifi_correct_bp": hifi_ok, "hifi_correct_pct_telomere": pct(hifi_ok, total),
            "ont_correct_bp": ont_ok, "ont_correct_pct_telomere": pct(ont_ok, total),
            "both_correct_bp": both, "both_correct_pct_telomere": pct(both, total), "both_correct_pct_teloscope": pct(both, total_bp),
            "one_correct_bp": one, "one_correct_pct_telomere": pct(one, total), "one_correct_pct_teloscope": pct(one, total_bp),
            "both_non_correct_bp": both_non_v, "both_non_correct_pct_telomere": pct(both_non_v, total), "both_non_correct_pct_teloscope": pct(both_non_v, total_bp),
            "concordant_bp": conc, "concordant_pct_telomere": pct(conc, total), "discordant_bp": disc, "discordant_pct_telomere": pct(disc, total),
            "hifi_rescue_bp": hifi_rescue, "hifi_rescue_pct_telomere": pct(hifi_rescue, total), "hifi_rescue_pct_teloscope": pct(hifi_rescue, total_bp),
            "ont_rescue_bp": ont_rescue, "ont_rescue_pct_telomere": pct(ont_rescue, total), "ont_rescue_pct_teloscope": pct(ont_rescue, total_bp),
            "winner": win, "quality_score": qscore, "dominant_pair": dom_pair,
        })

    per_tel_df = overlap_df.groupby("telomere_id").apply(per_tel).reset_index(drop=True)
    per_tel_df = per_tel_df.sort_values(["chr_display", "hap", "arm"], key=lambda s: s.map(chr_key) if s.name == "chr_display" else s)

    best = per_tel_df.sort_values(["both_correct_pct_telomere", "both_correct_bp"], ascending=[False, False]).head(15).copy()
    best.insert(0, "rank_group", "best"); best.insert(1, "rank", range(1, len(best) + 1))
    worst = per_tel_df.sort_values(["both_correct_pct_telomere", "both_correct_bp"], ascending=[True, True]).head(15).copy()
    worst.insert(0, "rank_group", "worst"); worst.insert(1, "rank", range(1, len(worst) + 1))
    best_worst = pd.concat([best, worst], ignore_index=True)

    hifi_res = per_tel_df.loc[per_tel_df["hifi_rescue_bp"] > 0].sort_values(["hifi_rescue_bp", "both_correct_bp"], ascending=[False, False]).copy()
    hifi_res.insert(0, "rescue_source", "HiFi"); hifi_res.insert(1, "rank", range(1, len(hifi_res) + 1))
    ont_res = per_tel_df.loc[per_tel_df["ont_rescue_bp"] > 0].sort_values(["ont_rescue_bp", "both_correct_bp"], ascending=[False, False]).copy()
    ont_res.insert(0, "rescue_source", "ONT"); ont_res.insert(1, "rank", range(1, len(ont_res) + 1))

    winners = per_tel_df.groupby("winner", as_index=False).size().rename(columns={"size": "n_telomeres"})
    winners["fraction_telomeres"] = winners["n_telomeres"].map(lambda x: pct(x, len(per_tel_df)))

    return {
        "telomere_overlap_segments": overlap_df,
        "global_quality_metrics": global_q,
        "tech_category_composition": tech_cat,
        "pair_matrix_long": pair_long,
        "pair_matrix_wide": pair_wide,
        "per_telomere_quality": per_tel_df,
        "best_worst_telomeres": best_worst,
        "rescue_ranking_hifi": hifi_res,
        "rescue_ranking_ont": ont_res,
        "winner_counts": winners,
    }

def build_interpretation(tables: Dict[str, pd.DataFrame], total_bp: int) -> pd.DataFrame:
    gq = tables["global_quality_metrics"].set_index("metric")
    per_t = tables["per_telomere_quality"]
    hifi_res = tables["rescue_ranking_hifi"]
    ont_res = tables["rescue_ranking_ont"]
    pairs = tables["pair_matrix_long"].copy()

    if per_t.empty:
        raise RuntimeError("Per-telomere table is empty; cannot build interpretation.")

    best = per_t.sort_values(["both_correct_pct_telomere", "both_correct_bp"], ascending=[False, False]).iloc[0]
    worst = per_t.sort_values(["both_correct_pct_telomere", "both_correct_bp"], ascending=[True, True]).iloc[0]
    top_hifi = hifi_res.iloc[0] if not hifi_res.empty else None
    top_ont = ont_res.iloc[0] if not ont_res.empty else None
    pairs = pairs.loc[~((pairs["hifi_name"] == "correct") & (pairs["ont_name"] == "correct"))]
    if pairs.empty:
        top_pair = pd.Series({"hifi_name": "NA", "ont_name": "NA", "bp": 0})
    else:
        top_pair = pairs.sort_values("bp", ascending=False).iloc[0]

    rows = [
        {"question": "2.2.1 assembly quality", "metric": "both_correct", "bp": int(gq.loc["both_correct", "bp"]), "pct_teloscope": pct(gq.loc["both_correct", "bp"], total_bp), "interpretation": "Bp with both technologies calling correct."},
        {"question": "2.2.1 assembly quality", "metric": "one_correct_rescued", "bp": int(gq.loc["one_correct_rescued", "bp"]), "pct_teloscope": pct(gq.loc["one_correct_rescued", "bp"], total_bp), "interpretation": "Bp rescued by only one technology."},
        {"question": "2.2.1 assembly quality", "metric": "both_non_correct", "bp": int(gq.loc["both_non_correct", "bp"]), "pct_teloscope": pct(gq.loc["both_non_correct", "bp"], total_bp), "interpretation": "Bp where neither technology is correct."},
        {"question": "2.2.1 assembly quality", "metric": "concordance", "bp": int(gq.loc["concordant", "bp"]), "pct_teloscope": pct(gq.loc["concordant", "bp"], total_bp), "interpretation": f"Concordant bp; discordant={fmt_bp(gq.loc['discordant', 'bp'])} bp ({100*pct(gq.loc['discordant', 'bp'], total_bp):.2f}%)."},
        {"question": "2.2.2 rescue", "metric": "hifi_rescue_total", "bp": int(per_t["hifi_rescue_bp"].sum()), "pct_teloscope": pct(per_t["hifi_rescue_bp"].sum(), total_bp), "interpretation": "Total HiFi-only rescue bp."},
        {"question": "2.2.2 rescue", "metric": "ont_rescue_total", "bp": int(per_t["ont_rescue_bp"].sum()), "pct_teloscope": pct(per_t["ont_rescue_bp"].sum(), total_bp), "interpretation": "Total ONT-only rescue bp."},
        {"question": "2.2.2 rescue", "metric": "top_hifi_rescued_telomere", "bp": int(top_hifi["hifi_rescue_bp"]) if top_hifi is not None else 0, "pct_teloscope": pct(top_hifi["hifi_rescue_bp"], total_bp) if top_hifi is not None else 0.0, "interpretation": f"{top_hifi['telomere_id']} has largest HiFi-only rescue." if top_hifi is not None else "No HiFi rescue detected."},
        {"question": "2.2.2 rescue", "metric": "top_ont_rescued_telomere", "bp": int(top_ont["ont_rescue_bp"]) if top_ont is not None else 0, "pct_teloscope": pct(top_ont["ont_rescue_bp"], total_bp) if top_ont is not None else 0.0, "interpretation": f"{top_ont['telomere_id']} has largest ONT-only rescue." if top_ont is not None else "No ONT rescue detected."},
        {"question": "2.2.3 best/worst", "metric": "best_telomere_by_both_correct", "bp": int(best["both_correct_bp"]), "pct_teloscope": pct(best["both_correct_bp"], total_bp), "interpretation": f"{best['telomere_id']} has {100*best['both_correct_pct_telomere']:.2f}% both-correct within telomere."},
        {"question": "2.2.3 best/worst", "metric": "worst_telomere_by_both_correct", "bp": int(worst["both_correct_bp"]), "pct_teloscope": pct(worst["both_correct_bp"], total_bp), "interpretation": f"{worst['telomere_id']} has {100*worst['both_correct_pct_telomere']:.2f}% both-correct within telomere."},
        {"question": "2.2.4 creative", "metric": "dominant_discordant_pair", "bp": int(top_pair["bp"]), "pct_teloscope": pct(top_pair["bp"], total_bp), "interpretation": f"Largest non-perfect pair is {top_pair['hifi_name']} (HiFi) vs {top_pair['ont_name']} (ONT)."},
    ]
    return pd.DataFrame(rows)


def write_interpretation_text(interp_df: pd.DataFrame, total_bp: int) -> None:
    lines = [
        "H9 paired NucFlag telomere interpretation",
        f"Total Teloscope telomere bp (100% reference): {fmt_bp(total_bp)} bp",
        "",
    ]
    for _, r in interp_df.iterrows():
        lines.append(f"- {r['question']} | {r['metric']}: {fmt_bp(r['bp'])} bp ({100*float(r['pct_teloscope']):.2f}% of Teloscope). {r['interpretation']}")
    INTERPRETATION_TXT.write_text("\n".join(lines), encoding="utf-8")
    print(f"[OK] Text: {INTERPRETATION_TXT.name}")


def plot_global_quality(global_q: pd.DataFrame) -> None:
    order = ["both_correct", "one_correct_rescued", "both_non_correct", "concordant", "discordant"]
    cmap = {"both_correct": "#2D8E47", "one_correct_rescued": "#D4A535", "both_non_correct": "#CB3335", "concordant": "#3C7DC4", "discordant": "#8B6D52"}
    p = global_q.set_index("metric").reindex(order).reset_index()
    vals = 100.0 * p["pct_teloscope"].values
    fig, ax = plt.subplots(figsize=(4.8, 2.5), constrained_layout=True)
    ax.bar(range(len(order)), vals, color=[cmap[k] for k in order], edgecolor="black", linewidth=0.25)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["Both\ncorrect", "One\ncorrect", "Both\nnon-correct", "Concordant", "Discordant"], fontsize=6)
    ax.set_ylabel("% of Teloscope telomere bp", fontsize=7)
    ax.set_title("Global assembly quality from paired HiFi+ONT labels", fontsize=7)
    ax.set_ylim(0, max(vals) * 1.15 if len(vals) else 1)
    for i, (bp_v, pct_v) in enumerate(zip(p["bp"], vals)):
        ax.text(i, pct_v + ax.get_ylim()[1] * 0.01, f"{fmt_bp(bp_v)}\n({pct_v:.1f}%)", ha="center", va="bottom", fontsize=5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    save_figure(fig, f"fig1_global_quality_{RUN_LABEL}")


def plot_tech_category(tech_cat: pd.DataFrame, cats: List[str], color_map: Dict[str, str]) -> None:
    techs = ["HiFi", "ONT"]
    fig, ax = plt.subplots(figsize=(4.4, 2.6), constrained_layout=True)
    bottoms = np.zeros(len(techs), dtype=float)
    for cat in cats:
        vals = []
        for tech in techs:
            m = tech_cat.loc[(tech_cat["technology"] == tech) & (tech_cat["category"] == cat), "pct_teloscope"]
            vals.append(100.0 * float(m.iloc[0]) if not m.empty else 0.0)
        ax.bar(techs, vals, bottom=bottoms, color=color_map.get(cat, "#BFBFBF"), edgecolor="none", label=cat)
        bottoms += np.array(vals)
    ax.set_ylabel("% of Teloscope telomere bp", fontsize=7)
    ax.set_title("Technology-wise category composition within telomeres", fontsize=7)
    ax.set_ylim(0, max(bottoms) * 1.1 if len(bottoms) else 1)
    handles = [Patch(facecolor=color_map.get(c, "#BFBFBF"), label=c) for c in cats]
    ax.legend(handles=handles, fontsize=5.5, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.23))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    save_figure(fig, f"fig2_tech_category_stacked_{RUN_LABEL}")


def plot_pair_heatmap(pair_long: pd.DataFrame, cats: List[str]) -> None:
    if not cats:
        print("[WARN] No categories available for pair heatmap; skipping.")
        return
    mat = pair_long.pivot_table(index="hifi_name", columns="ont_name", values="pct_teloscope", fill_value=0)
    mat = mat.reindex(index=cats, columns=cats, fill_value=0)
    data = mat.values * 100.0
    n = len(cats)
    fig, ax = plt.subplots(figsize=(max(5.0, 0.32 * n + 1.8), max(4.0, 0.32 * n + 1.4)), constrained_layout=True)
    im = ax.imshow(data, cmap="YlOrRd", vmin=0)
    ax.set_xticks(np.arange(n)); ax.set_yticks(np.arange(n))
    ax.set_xticklabels(cats, rotation=60, ha="right", fontsize=6); ax.set_yticklabels(cats, fontsize=6)
    ax.set_xlabel("ONT category", fontsize=7); ax.set_ylabel("HiFi category", fontsize=7)
    ax.set_title("HiFi x ONT category-pair matrix (% of Teloscope bp)", fontsize=7)
    max_v = float(np.nanmax(data)) if data.size else 0.0
    for i in range(n):
        for j in range(n):
            v = data[i, j]
            if v >= 0.5:
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5, color="black" if v < max_v * 0.65 else "white")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.set_ylabel("% of Teloscope bp", fontsize=7); cb.ax.tick_params(labelsize=6)
    save_figure(fig, f"fig3_hifi_ont_pair_heatmap_{RUN_LABEL}")

def plot_quality_scatter(per_tel: pd.DataFrame) -> None:
    if per_tel.empty:
        print("[WARN] No per-telomere data for quality scatter; skipping.")
        return
    x = per_tel["hifi_correct_pct_telomere"].values * 100.0
    y = per_tel["ont_correct_pct_telomere"].values * 100.0
    c = per_tel["both_correct_pct_telomere"].values * 100.0
    max_bp = float(per_tel["total_bp"].max()) if "total_bp" in per_tel.columns else 0.0
    if max_bp <= 0:
        max_bp = 1.0
    s = 25 + 180 * np.sqrt(per_tel["total_bp"].values / max_bp)
    fig, ax = plt.subplots(figsize=(4.0, 3.6), constrained_layout=True)
    sc = ax.scatter(x, y, c=c, s=s, cmap="viridis", edgecolors="black", linewidths=0.25, alpha=0.9)
    ax.plot([0, 100], [0, 100], linestyle="--", linewidth=0.6, color="gray")
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_xlabel("HiFi correct within telomere (%)", fontsize=7)
    ax.set_ylabel("ONT correct within telomere (%)", fontsize=7)
    ax.set_title("Per-telomere quality: HiFi vs ONT correctness", fontsize=7)
    best = per_tel.sort_values("both_correct_pct_telomere", ascending=False).head(3)
    worst = per_tel.sort_values("both_correct_pct_telomere", ascending=True).head(3)
    for _, r in pd.concat([best, worst], ignore_index=True).iterrows():
        label = f"{r['chr_display']}_{r['hap'][-1]}{r['arm']}"
        ax.annotate(label, (100 * r["hifi_correct_pct_telomere"], 100 * r["ont_correct_pct_telomere"]), textcoords="offset points", xytext=(3, 3), fontsize=5)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.set_ylabel("Both-correct within telomere (%)", fontsize=7); cb.ax.tick_params(labelsize=6)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    save_figure(fig, f"fig4_telomere_quality_scatter_{RUN_LABEL}")


def plot_rescue_diverging(per_tel: pd.DataFrame) -> None:
    rescue = per_tel.loc[(per_tel["hifi_rescue_bp"] > 0) | (per_tel["ont_rescue_bp"] > 0)].copy()
    if rescue.empty:
        print("[WARN] No rescue signal found; skipping rescue plot.")
        return
    rescue["max_rescue_bp"] = rescue[["hifi_rescue_bp", "ont_rescue_bp"]].max(axis=1)
    rescue = rescue.sort_values("max_rescue_bp", ascending=False).head(20).sort_values("max_rescue_bp", ascending=True)
    labels = rescue["telomere_id"].tolist(); y = np.arange(len(rescue))
    fig, ax = plt.subplots(figsize=(6.1, max(3.2, 0.19 * len(rescue) + 1.4)), constrained_layout=True)
    ax.barh(y, -rescue["hifi_rescue_bp"].values, color="#E69F00", label="HiFi-only rescue")
    ax.barh(y, rescue["ont_rescue_bp"].values, color="#3C7DC4", label="ONT-only rescue")
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.5)
    ax.set_xlabel("Rescued bp (negative = HiFi, positive = ONT)", fontsize=7)
    ax.set_title("Top telomeres rescued by one technology but not the other", fontsize=7)
    lim = max(rescue["hifi_rescue_bp"].max(), rescue["ont_rescue_bp"].max())
    ax.set_xlim(-1.1 * lim, 1.1 * lim)
    ax.legend(fontsize=6, frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    save_figure(fig, f"fig5_rescue_diverging_{RUN_LABEL}")


def main() -> None:
    init_plot_style()
    for d in (PLOTS_DIR, STATS_DIR, FIGDATA_DIR):
        d.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("H9 telomere paired NucFlag analysis")
    print("=" * 90)

    for p in (COMBINED_TSV, TELOSCOPE_BED, PALETTE_TSV):
        if not p.is_file():
            raise FileNotFoundError(f"Missing input file: {p}")

    palette_df, color_map = load_palette(PALETTE_TSV)
    telo_df = load_teloscope(TELOSCOPE_BED)
    comb_df = load_combined(COMBINED_TSV)
    overlap_df = intersect_telomeres(comb_df, telo_df)
    total_bp = int(telo_df["length"].sum())
    cats = ordered_categories(overlap_df, palette_df)
    print(f"[INFO] Overlap categories: {cats}")

    tables = build_tables(overlap_df, telo_df, cats)

    write_tsv(tables["telomere_overlap_segments"], STATS_DIR / f"table_telomere_overlap_segments_{RUN_LABEL}.tsv", "Paired overlap segments")
    write_tsv(tables["global_quality_metrics"], STATS_DIR / f"table_global_quality_metrics_{RUN_LABEL}.tsv", "Global quality metrics")
    write_tsv(tables["tech_category_composition"], STATS_DIR / f"table_tech_category_composition_{RUN_LABEL}.tsv", "Per-tech category composition")
    write_tsv(tables["pair_matrix_long"], STATS_DIR / f"table_hifi_ont_pair_matrix_long_{RUN_LABEL}.tsv", "HiFi x ONT pair matrix long")
    write_tsv(tables["pair_matrix_wide"], STATS_DIR / f"table_hifi_ont_pair_matrix_wide_{RUN_LABEL}.tsv", "HiFi x ONT pair matrix wide")
    write_tsv(tables["per_telomere_quality"], STATS_DIR / f"table_per_telomere_quality_{RUN_LABEL}.tsv", "Per-telomere quality metrics")
    write_tsv(tables["best_worst_telomeres"], STATS_DIR / f"table_best_worst_telomeres_{RUN_LABEL}.tsv", "Best/worst by both-correct")
    write_tsv(tables["rescue_ranking_hifi"], STATS_DIR / f"table_rescue_ranking_hifi_{RUN_LABEL}.tsv", "HiFi rescue ranking")
    write_tsv(tables["rescue_ranking_ont"], STATS_DIR / f"table_rescue_ranking_ont_{RUN_LABEL}.tsv", "ONT rescue ranking")
    write_tsv(tables["winner_counts"], STATS_DIR / f"table_technology_winner_counts_{RUN_LABEL}.tsv", "Technology winner counts")

    write_tsv(tables["global_quality_metrics"], FIGDATA_DIR / f"fig1_global_quality_data_{RUN_LABEL}.tsv", "Figure 1 data")
    write_tsv(tables["tech_category_composition"], FIGDATA_DIR / f"fig2_tech_category_data_{RUN_LABEL}.tsv", "Figure 2 data")
    write_tsv(tables["pair_matrix_long"], FIGDATA_DIR / f"fig3_pair_heatmap_data_{RUN_LABEL}.tsv", "Figure 3 data")
    write_tsv(tables["per_telomere_quality"], FIGDATA_DIR / f"fig4_quality_scatter_data_{RUN_LABEL}.tsv", "Figure 4 data")
    write_tsv(tables["per_telomere_quality"], FIGDATA_DIR / f"fig5_rescue_data_{RUN_LABEL}.tsv", "Figure 5 data")

    interp_df = build_interpretation(tables, total_bp)
    write_tsv(interp_df, STATS_DIR / f"table_interpretation_summary_{RUN_LABEL}.tsv", "Question-targeted interpretation")
    write_interpretation_text(interp_df, total_bp)

    plot_global_quality(tables["global_quality_metrics"])
    plot_tech_category(tables["tech_category_composition"], cats, color_map)
    plot_pair_heatmap(tables["pair_matrix_long"], cats)
    plot_quality_scatter(tables["per_telomere_quality"])
    plot_rescue_diverging(tables["per_telomere_quality"])

    parts = []
    for name, df in tables.items():
        x = df.copy(); x.insert(0, "section", name); parts.append(x)
    y = interp_df.copy(); y.insert(0, "section", "interpretation_summary"); parts.append(y)
    combined_df = pd.concat(parts, ignore_index=True, sort=False)
    combined_df.to_csv(COMBINED_STATS_OUT, sep="\t", index=False)
    print(f"[OK] Combined stats: {COMBINED_STATS_OUT}")

    print("-" * 90)
    print("DONE")
    print(f"Figures: {PLOTS_DIR}")
    print(f"Stats:   {STATS_DIR}")
    print("-" * 90)


if __name__ == "__main__":
    main()
