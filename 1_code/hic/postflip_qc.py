#!/usr/bin/env python3
"""
Post-flip QC: validate E1 eigenvector phasing after sign correction.

Reads flipped bedgraph files and coding_cov bedgraph files, computes
correlations and Delta on flipped E1 vs coding coverage, generates QC
scatter plots and summary barplots, and performs hap1-vs-hap2 concordance.

Exit code 0 = QC pass, 1 = QC fail.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr


HAPS = ["hap1", "hap2"]
RESOLUTIONS = [5000, 10000, 20000, 50000, 100000, 200000]
BASE = "H9_T2T_v0.1"

MIN_AUTOSOME_POSITIVE_FRAC = 0.8  # >=80% of autosomes must pass


def compute_delta(e1_values, coding_values):
    """Delta = median(coding | top 20% E1) - median(coding | bottom 20% E1)."""
    e1 = np.asarray(e1_values, dtype=float)
    cod = np.asarray(coding_values, dtype=float)
    if len(e1) < 10:
        return np.nan
    q80, q20 = np.percentile(e1, 80), np.percentile(e1, 20)
    top, bot = e1 >= q80, e1 <= q20
    if top.sum() < 2 or bot.sum() < 2:
        return np.nan
    return float(np.median(cod[top]) - np.median(cod[bot]))


def strip_hap_suffix(chrom_name):
    """'chr1_hap1' -> 'chr1'."""
    for suffix in ("_hap1", "_hap2"):
        if chrom_name.endswith(suffix):
            return chrom_name[: -len(suffix)]
    return chrom_name


def run_postflip_qc():
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

    data_dir = os.path.join(REPO_ROOT, "2_data", "2.2_processed")
    flipped_dir = os.path.join(data_dir, "26.02.07_cooltools_flipped_E1")
    fig_dir = os.path.join(REPO_ROOT, "3_figures", "3.1_draft", "26.02.08_hic")
    os.makedirs(fig_dir, exist_ok=True)

    all_qc_results = []
    hap_chrom_e1 = {}  # (hap, res, base_chrom) -> Series indexed by (start, end)

    for hap in HAPS:
        for res in RESOLUTIONS:
            flipped_path = os.path.join(
                flipped_dir, f"{BASE}_{hap}.{res}.AB.flipped.bedgraph"
            )
            coding_path = os.path.join(
                data_dir, f"{BASE}_{hap}.bwa.{res}.coding_cov.bedgraph"
            )

            if not os.path.exists(flipped_path) or not os.path.exists(coding_path):
                print(f"[!] Missing files for {hap} @ {res}, skipping")
                continue

            df_flipped = pd.read_csv(
                flipped_path,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "E1_flipped"],
            )

            df_coding = pd.read_csv(coding_path, sep="\t")
            value_col = df_coding.columns[4]
            if value_col != "coding":
                df_coding.rename(columns={value_col: "coding"}, inplace=True)

            df = pd.merge(
                df_flipped,
                df_coding[["chrom", "start", "end", "coding"]],
                on=["chrom", "start", "end"],
                how="inner",
            ).dropna(subset=["E1_flipped", "coding"])

            out_dir = os.path.join(fig_dir, f"{BASE}_{hap}_{res}")
            os.makedirs(out_dir, exist_ok=True)

            chrom_stats = {}
            for chrom, grp in df.groupby("chrom"):
                if len(grp) < 10:
                    continue

                rho, rho_p = spearmanr(grp["coding"], grp["E1_flipped"])
                r, r_p = pearsonr(grp["coding"], grp["E1_flipped"])
                delta = compute_delta(grp["E1_flipped"].values, grp["coding"].values)

                chrom_stats[chrom] = {
                    "spearman_rho": rho,
                    "spearman_p": rho_p,
                    "pearson_r": r,
                    "pearson_p": r_p,
                    "delta": delta,
                }

                base_chrom = strip_hap_suffix(chrom)
                hap_chrom_e1[(hap, res, base_chrom)] = grp.set_index(
                    ["start", "end"]
                )["E1_flipped"]

                all_qc_results.append(
                    {
                        "hap": hap,
                        "res": res,
                        "chrom": chrom,
                        "spearman_rho": rho,
                        "delta": delta,
                    }
                )

            chroms_sorted = sorted(chrom_stats.keys())

            # Summary TSV
            summary_rows = [
                {"chrom": c, **chrom_stats[c]} for c in chroms_sorted
            ]
            pd.DataFrame(summary_rows).to_csv(
                os.path.join(out_dir, "postflip_coding_cov_summary_stats.tsv"),
                sep="\t",
                index=False,
            )

            # Per-chromosome scatter plots
            for chrom in chroms_sorted:
                grp = df[df["chrom"] == chrom]
                stats = chrom_stats[chrom]
                fig, ax = plt.subplots(figsize=(4, 4))
                ax.scatter(grp["coding"], grp["E1_flipped"], s=5, alpha=0.5)
                ax.axhline(0, color="gray", linewidth=1)
                ax.set_xlabel("Coding gene coverage")
                ax.set_ylabel("E1 (flipped)")
                stats_text = (
                    f"\u03c1 = {stats['spearman_rho']:.3f}\n"
                    f"\u0394 = {stats['delta']:.4f}"
                )
                ax.text(
                    0.95,
                    0.95,
                    stats_text,
                    transform=ax.transAxes,
                    fontsize=7,
                    va="top",
                    ha="right",
                    bbox=dict(
                        boxstyle="round,pad=0.3", facecolor="lightgrey", alpha=0.9
                    ),
                )
                fig.tight_layout()
                fig.savefig(
                    os.path.join(out_dir, f"postflip_{chrom}_coding.png"), dpi=600
                )
                plt.close(fig)

            # Summary barplot: Spearman
            spearman_vals = [chrom_stats[c]["spearman_rho"] for c in chroms_sorted]
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.bar(chroms_sorted, spearman_vals)
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_xticks(range(len(chroms_sorted)))
            ax.set_xticklabels(chroms_sorted, rotation=90, fontsize=6)
            ax.set_ylabel("Spearman \u03c1 (flipped E1 vs coding_cov)")
            fig.tight_layout()
            fig.savefig(
                os.path.join(
                    out_dir, f"postflip_{hap}_{res}_spearman_summary.png"
                ),
                dpi=600,
            )
            plt.close(fig)

            # Summary barplot: Delta
            delta_vals = [chrom_stats[c]["delta"] for c in chroms_sorted]
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.bar(chroms_sorted, delta_vals)
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_xticks(range(len(chroms_sorted)))
            ax.set_xticklabels(chroms_sorted, rotation=90, fontsize=6)
            ax.set_ylabel("\u0394 (flipped E1 vs coding_cov)")
            fig.tight_layout()
            fig.savefig(
                os.path.join(
                    out_dir, f"postflip_{hap}_{res}_delta_summary.png"
                ),
                dpi=600,
            )
            plt.close(fig)

    # ── Hap1 vs hap2 concordance ────────────────────────────────────────
    hap_concordance = []
    for res in RESOLUTIONS:
        hap1_chroms = {k[2] for k in hap_chrom_e1 if k[0] == "hap1" and k[1] == res}
        hap2_chroms = {k[2] for k in hap_chrom_e1 if k[0] == "hap2" and k[1] == res}
        shared = sorted(hap1_chroms & hap2_chroms)

        for base_chrom in shared:
            s1 = hap_chrom_e1.get(("hap1", res, base_chrom))
            s2 = hap_chrom_e1.get(("hap2", res, base_chrom))
            if s1 is None or s2 is None:
                continue

            joined = pd.concat([s1, s2], axis=1, join="inner")
            joined.columns = ["hap1", "hap2"]
            joined = joined.dropna()
            if len(joined) < 10:
                continue

            rho, _ = spearmanr(joined["hap1"], joined["hap2"])
            hap_concordance.append(
                {
                    "res": res,
                    "chrom": base_chrom,
                    "hap1_vs_hap2_spearman": rho,
                    "n_bins": len(joined),
                }
            )

    if hap_concordance:
        df_conc = pd.DataFrame(hap_concordance)
        df_conc.to_csv(
            os.path.join(fig_dir, "hap1_vs_hap2_concordance.tsv"),
            sep="\t",
            index=False,
        )

    # ── QC pass/fail summary ────────────────────────────────────────────
    df_qc = pd.DataFrame(all_qc_results)
    print("\n" + "=" * 70)
    print("POST-FLIP QC SUMMARY")
    print("=" * 70)

    qc_pass = True

    for res in RESOLUTIONS:
        for hap in HAPS:
            subset = df_qc[(df_qc["hap"] == hap) & (df_qc["res"] == res)]
            if len(subset) == 0:
                continue

            n_pos_rho = (subset["spearman_rho"] > 0).sum()
            n_pos_delta = (subset["delta"] > 0).sum()
            n_total = len(subset)
            frac_rho = n_pos_rho / n_total
            frac_delta = n_pos_delta / n_total

            rho_ok = frac_rho >= MIN_AUTOSOME_POSITIVE_FRAC
            delta_ok = frac_delta >= MIN_AUTOSOME_POSITIVE_FRAC

            status = "PASS" if (rho_ok and delta_ok) else "FAIL"
            if status == "FAIL":
                qc_pass = False

            print(f"\n  {hap} @ {res}:")
            print(
                f"    Spearman > 0: {n_pos_rho}/{n_total} "
                f"({frac_rho:.0%}) {'OK' if rho_ok else 'FAIL'}"
            )
            print(
                f"    Delta > 0:    {n_pos_delta}/{n_total} "
                f"({frac_delta:.0%}) {'OK' if delta_ok else 'FAIL'}"
            )

            failing = subset[
                (subset["spearman_rho"] <= 0) | (subset["delta"] <= 0)
            ]
            if len(failing) > 0:
                for _, row in failing.iterrows():
                    print(
                        f"      WARNING: {row['chrom']} "
                        f"rho={row['spearman_rho']:.3f} delta={row['delta']:.4f}"
                    )

    # Hap concordance QC
    if hap_concordance:
        df_conc = pd.DataFrame(hap_concordance)
        print(f"\n  Hap1 vs Hap2 concordance:")
        for res in RESOLUTIONS:
            res_conc = df_conc[df_conc["res"] == res]
            if len(res_conc) == 0:
                continue
            mean_rho = res_conc["hap1_vs_hap2_spearman"].mean()
            n_negative = (res_conc["hap1_vs_hap2_spearman"] < 0).sum()
            print(
                f"    {res}: mean rho = {mean_rho:.3f}, "
                f"anticorrelated = {n_negative}/{len(res_conc)}"
            )
            if n_negative > len(res_conc) * 0.2:
                print(
                    f"      FAIL: >20% of chromosomes anticorrelated between haps"
                )
                qc_pass = False

    print(f"\n{'=' * 70}")
    print(f"OVERALL QC: {'PASS' if qc_pass else 'FAIL'}")
    print(f"{'=' * 70}\n")

    return 0 if qc_pass else 1


if __name__ == "__main__":
    sys.exit(run_postflip_qc())
