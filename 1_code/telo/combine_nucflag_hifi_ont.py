#!/usr/bin/env python3
"""
Combine NucFlag telomere annotations from HiFi and ONT into a single file.

For every chromosome, the breakpoints from both data sets are merged so that
each output interval carries both the HiFi and ONT NucFlag annotation that
covers that sub-interval.  This enables direct positional comparison of the
two technologies.

Input
-----
    nucflag_telo_hifi.tsv   BED9+2 (NucFlag annotations, HiFi reads)
    nucflag_telo_ont.tsv    BED9+2 (NucFlag annotations, ONT reads)

Output
------
    nucflag_telo_combined.tsv
        chrom, chromStart, chromEnd, length,
        hifi_name, hifi_score, hifi_zscore, hifi_af,
        ont_name,  ont_score,  ont_zscore,  ont_af,
        concordant   (bool: hifi_name == ont_name)

Usage
-----
    python 1_code/telo/combine_nucflag_hifi_ont.py
"""

import os
import sys
from typing import List, Tuple

import numpy as np
import pandas as pd

# ─── Paths ────────────────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir   = os.path.dirname(os.path.dirname(script_dir))
data_dir   = os.path.join(base_dir, "2_data", "2.2_processed", "nucflag")

HIFI_TSV = os.path.join(data_dir, "nucflag_telo_hifi.tsv")
ONT_TSV  = os.path.join(data_dir, "nucflag_telo_ont.tsv")
OUT_TSV  = os.path.join(data_dir, "nucflag_telo_combined.tsv")

# ─── Column names (BED9 + zscore + af) ───────────────────────────────────────
BED_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "zscore", "af",
]

# ─── Read a NucFlag TSV ──────────────────────────────────────────────────────
def read_nucflag(path: str) -> pd.DataFrame:
    """Read a BED9+2 NucFlag TSV, returning only the columns we need."""
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=BED_COLS)
    return df[["chrom", "chromStart", "chromEnd", "name", "score", "zscore", "af"]]


def merge_intervals(
    hifi: pd.DataFrame,
    ont: pd.DataFrame,
    chrom: str,
) -> List[dict]:
    """
    Merge HiFi and ONT intervals for one chromosome using a
    breakpoint-union (sweep-line) approach.

    Returns a list of dicts, one per sub-interval.
    """
    h = hifi[hifi["chrom"] == chrom].sort_values("chromStart")
    o = ont[ont["chrom"] == chrom].sort_values("chromStart")

    if h.empty and o.empty:
        return []

    # Collect all unique breakpoints
    bp = sorted(
        set(h["chromStart"].tolist() + h["chromEnd"].tolist() +
            o["chromStart"].tolist() + o["chromEnd"].tolist())
    )

    rows = []
    for i in range(len(bp) - 1):
        s, e = bp[i], bp[i + 1]
        if s >= e:
            continue

        # Find HiFi interval covering this sub-interval
        hm = h[(h["chromStart"] <= s) & (h["chromEnd"] >= e)]
        # Find ONT interval covering this sub-interval
        om = o[(o["chromStart"] <= s) & (o["chromEnd"] >= e)]

        hifi_name  = hm.iloc[0]["name"]   if not hm.empty else "uncovered"
        hifi_score = hm.iloc[0]["score"]  if not hm.empty else np.nan
        hifi_zsc   = hm.iloc[0]["zscore"] if not hm.empty else np.nan
        hifi_af    = hm.iloc[0]["af"]     if not hm.empty else np.nan

        ont_name   = om.iloc[0]["name"]   if not om.empty else "uncovered"
        ont_score  = om.iloc[0]["score"]  if not om.empty else np.nan
        ont_zsc    = om.iloc[0]["zscore"] if not om.empty else np.nan
        ont_af     = om.iloc[0]["af"]     if not om.empty else np.nan

        rows.append({
            "chrom":       chrom,
            "chromStart":  s,
            "chromEnd":    e,
            "length":      e - s,
            "hifi_name":   hifi_name,
            "hifi_score":  hifi_score,
            "hifi_zscore": hifi_zsc,
            "hifi_af":     hifi_af,
            "ont_name":    ont_name,
            "ont_score":   ont_score,
            "ont_zscore":  ont_zsc,
            "ont_af":      ont_af,
            "concordant":  hifi_name == ont_name,
        })

    return rows


# ─── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    for p in (HIFI_TSV, ONT_TSV):
        if not os.path.isfile(p):
            sys.exit(f"Error: {p} not found")

    print(f"Reading  {HIFI_TSV}")
    hifi = read_nucflag(HIFI_TSV)
    print(f"  rows: {len(hifi)}")

    print(f"Reading  {ONT_TSV}")
    ont = read_nucflag(ONT_TSV)
    print(f"  rows: {len(ont)}")

    chroms = sorted(
        set(hifi["chrom"].unique()) | set(ont["chrom"].unique()),
        key=lambda c: (
            int(c.split("_")[0].replace("chr", "").replace("X", "23")),
            c.split("_")[1],
        ),
    )

    print(f"\nMerging {len(chroms)} chromosomes …")
    all_rows: List[dict] = []
    for chrom in chroms:
        merged = merge_intervals(hifi, ont, chrom)
        all_rows.extend(merged)
        print(f"  {chrom:20s} → {len(merged):4d} sub-intervals")

    combined = pd.DataFrame(all_rows)

    # ── Summary statistics ────────────────────────────────────────────────
    total = len(combined)
    conc  = combined["concordant"].sum()
    print(f"\nTotal sub-intervals : {total}")
    print(f"Concordant          : {conc}  ({100 * conc / total:.1f}%)")
    print(f"Discordant          : {total - conc}  ({100 * (total - conc) / total:.1f}%)")

    total_bp   = combined["length"].sum()
    conc_bp    = combined.loc[combined["concordant"], "length"].sum()
    print(f"\nTotal bp            : {total_bp:,}")
    print(f"Concordant bp       : {conc_bp:,}  ({100 * conc_bp / total_bp:.1f}%)")
    print(f"Discordant bp       : {total_bp - conc_bp:,}  ({100 * (total_bp - conc_bp) / total_bp:.1f}%)")

    # ── Cross-tab of HiFi × ONT categories ───────────────────────────────
    xtab = pd.crosstab(
        combined["hifi_name"], combined["ont_name"],
        values=combined["length"], aggfunc="sum",
    ).fillna(0).astype(int)
    print("\nHiFi × ONT cross-tabulation (bp):")
    print(xtab.to_string())

    # ── Write output ──────────────────────────────────────────────────────
    combined.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\nWrote {OUT_TSV}")


if __name__ == "__main__":
    main()
