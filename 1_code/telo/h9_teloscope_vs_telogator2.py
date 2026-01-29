# Given the current tree:(base) aldhair@jackalienware:/mnt/d/research/h9_annotation$ tree -L 4
# .
# ├── 2_data
# │   ├── 2.1_raw
# │   │   ├── GCA_000001405.29_GRCh38.p14_genomic.chr.fna
# │   │   ├── GCA_000001405.29_GRCh38.p14_genomic.fna
# │   │   ├── GCA_009914755.4_T2T-CHM13v2.0_genomic.chr.fna
# │   │   ├── GCA_009914755.4_T2T-CHM13v2.0_genomic.fna
# │   │   ├── GCA_018852605.3_hg002v1.1.pat_genomic.fna
# │   │   ├── GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna
# │   │   ├── GCA_018852615.3_hg002v1.1.mat_genomic.fna
# │   │   ├── GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.chr.fna
# │   │   ├── GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.fna
# │   │   ├── GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.chr.fna
# │   │   ├── GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.fna
# │   │   ├── GWHDOOG00000000.genome.chr.fasta.gz
# │   │   ├── GWHDOOG00000000.genome.fasta.gz
# │   │   ├── GWHDQZJ00000000.genome.chr.fasta.gz
# │   │   ├── GWHDQZJ00000000.genome.fasta.gz
# │   │   ├── GWHGEYB00000000.1.genome.fasta.gz
# │   │   ├── GWHGEYC00000000.1.genome.fasta.gz
# │   │   ├── H9_T2T_v0.1_hap1.fasta
# │   │   ├── H9_T2T_v0.1_hap2.fasta
# │   │   ├── I002Cv0.7.hap1.fasta.gz
# │   │   └── I002Cv0.7.hap2.fasta.gz
# │   ├── 2.2_processed
# │   │   ├── 25.12.10_asms_x1_TTAGGG_v1.3.assembly_metrics.tsv
# │   │   ├── 25.12.10_asms_x1_TTAGGG_v1.3.path_summary.tsv
# │   │   ├── 25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres.bed
# │   │   ├── 25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres_extended.bed
# │   │   ├── 25.12_10_asms_x1_TTAGGG_v1.3.terminal_telomeres_filtered.bed
# │   │   ├── 26.01.07_telogator2_hifi_default
# │   │   │   ├── all_final_alleles.png
# │   │   │   ├── qc
# │   │   │   ├── temp
# │   │   │   ├── tlens_by_allele.tsv
# │   │   │   └── violin_atl.png
# │   │   ├── 26.01.08_telogator2_ont_default
# │   │   │   ├── all_final_alleles.png
# │   │   │   ├── qc
# │   │   │   ├── temp
# │   │   │   ├── tlens_by_allele.tsv
# │   │   │   └── violin_atl.png

# Compare telomere lengths from Teloscope (assembly-based) vs Telogator2 (read-based).
#
# Teloscope: measures TL directly from assembled genome
# Telogator2: measures TL from individual reads, reports sorted list in "read_TLs" column
#
# We compute from read_TLs: TL_p75, TL_p90, TL_max
# We compare: TL alone and TL+TVR (telomere variant repeats)
# We report: Pearson and Spearman correlations
#
# Output: 12 plots (3 percentiles x 2 metrics x 2 correlations)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy import stats

# Science-ready plot style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.linewidth': 1.2,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# Define paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "2_data" / "2.2_processed"
OUTPUT_DIR = BASE_DIR / "3_figures" / "3.1_draft" / "26.01.12_teloscope_vs_telogator2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Teloscope files
TELOSCOPE_HAP1 = DATA_DIR / "teloscope" / "H9_T2T_v0.1_hap1.fasta_terminal_telomeres.bed"
TELOSCOPE_HAP2 = DATA_DIR / "teloscope" / "H9_T2T_v0.1_hap2.fasta_terminal_telomeres.bed"

# Telogator2 files
TELOGATOR2_HIFI = DATA_DIR / "26.01.07_telogator2_hifi_default" / "tlens_by_allele.tsv"
TELOGATOR2_ONT = DATA_DIR / "26.01.08_telogator2_ont_default" / "tlens_by_allele.tsv"


def load_teloscope_data(hap1_path: Path, hap2_path: Path) -> pd.DataFrame:
    """Load and combine Teloscope data from both haplotypes."""
    cols = ["chr", "start", "end", "length", "label", "fwdCounts", "revCounts",
            "canCounts", "nonCanCounts", "chrSize"]

    df1 = pd.read_csv(hap1_path, sep="\t", header=None, names=cols)
    df2 = pd.read_csv(hap2_path, sep="\t", header=None, names=cols)

    df = pd.concat([df1, df2], ignore_index=True)

    # Extract chromosome number and create a key matching Telogator2 format
    # e.g., chr1_hap1 + p -> chr1p
    df["chr_num"] = df["chr"].str.extract(r"(chr\d+|chrX|chrY)")[0]
    df["teloscope_key"] = df["chr_num"] + df["label"]

    return df[["teloscope_key", "length", "chr"]].rename(columns={"length": "teloscope_length"})


def parse_read_tls(read_tls_str: str) -> list:
    """Parse the read_TLs column into a list of integers."""
    return [int(x) for x in read_tls_str.split(",")]


def load_telogator2_data(filepath: Path) -> pd.DataFrame:
    """Load Telogator2 data and compute TL statistics from read_TLs column.

    Each row has a sorted list of TL measurements in read_TLs (e.g., "-267,-13,2882,2984").
    We compute: TL_p75, TL_p90, TL_max from these values.

    Note: Telogator2 may report multiple alleles per chromosome arm (different ref_samp).
    We aggregate by chromosome arm taking median across alleles.
    """
    df = pd.read_csv(filepath, sep="\t")

    # Skip rows with multiple chr assignments (e.g., "chr5q,chr1p")
    df = df[~df["#chr"].str.contains(",", na=False)].copy()

    # Parse read_TLs and compute percentiles for each row
    df["read_TLs_parsed"] = df["read_TLs"].apply(parse_read_tls)
    df["TL_p75"] = df["read_TLs_parsed"].apply(lambda x: np.percentile(x, 75))
    df["TL_p90"] = df["read_TLs_parsed"].apply(lambda x: np.percentile(x, 90))
    df["TL_max"] = df["read_TLs_parsed"].apply(lambda x: max(x))

    # Aggregate by chromosome arm (median across multiple alleles)
    agg_df = df.groupby("#chr").agg({
        "TL_p75": "median",
        "TL_p90": "median",
        "TL_max": "median",
        "tvr_len": "median"
    }).reset_index()

    agg_df = agg_df.rename(columns={"#chr": "telogator_key"})

    # Compute TL + TVR for each percentile
    agg_df["TL_p75_plus_tvr"] = agg_df["TL_p75"] + agg_df["tvr_len"]
    agg_df["TL_p90_plus_tvr"] = agg_df["TL_p90"] + agg_df["tvr_len"]
    agg_df["TL_max_plus_tvr"] = agg_df["TL_max"] + agg_df["tvr_len"]

    return agg_df


def create_scatter_plot(teloscope_df: pd.DataFrame, telogator_df: pd.DataFrame,
                        telogator_col: str, output_path: Path, y_label: str,
                        correlation_type: str = "pearson"):
    """Create a minimalist science-ready scatter plot.

    Args:
        correlation_type: "pearson" or "spearman"
    """
    merged = teloscope_df.merge(telogator_df, left_on="teloscope_key",
                                right_on="telogator_key", how="inner")

    if merged.empty:
        print(f"Warning: No matching data for {output_path.name}")
        return None

    x = merged["teloscope_length"].values
    y = merged[telogator_col].values

    # Compute correlations
    if correlation_type == "pearson":
        r, p_value = stats.pearsonr(x, y)
        corr_label = "r"
    else:
        r, p_value = stats.spearmanr(x, y)
        corr_label = "ρ"

    # Linear regression for trend line
    slope, intercept, _, _, _ = stats.linregress(x, y)

    # Create figure
    fig, ax = plt.subplots(figsize=(4.5, 4.5))

    # Scatter points
    ax.scatter(x, y, s=40, c="#2C3E50", alpha=0.7, edgecolors="white",
               linewidth=0.5, zorder=3)

    # Regression line
    x_range = np.array([x.min(), x.max()])
    ax.plot(x_range, slope * x_range + intercept, color="#E74C3C",
            linewidth=1.5, linestyle="-", zorder=2)

    # Identity line (y=x)
    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lims, lims, color="#95A5A6", linewidth=1, linestyle="--",
            alpha=0.8, zorder=1)

    # Axis labels
    ax.set_xlabel("Teloscope (bp)", fontsize=11, fontweight="medium")
    ax.set_ylabel(y_label, fontsize=11, fontweight="medium")

    # Correlation annotation (bottom-right corner)
    corr_text = f"{corr_label} = {r:.2f}\np = {p_value:.1e}\nn = {len(merged)}"
    ax.text(0.97, 0.03, corr_text, transform=ax.transAxes, fontsize=9,
            verticalalignment="bottom", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#BDC3C7", alpha=0.9))

    # Clean up axes
    ax.tick_params(axis="both", which="major", labelsize=9)

    # Equal aspect with some padding
    ax.set_aspect("equal", adjustable="datalim")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, facecolor="white", edgecolor="none")
    plt.close()

    return {"r": r, "p": p_value, "n": len(merged), "type": correlation_type}


def main():
    print("=" * 60)
    print("Teloscope vs Telogator2 Comparison")
    print("=" * 60)

    # Load data
    print("\nLoading Teloscope data...")
    teloscope_df = load_teloscope_data(TELOSCOPE_HAP1, TELOSCOPE_HAP2)
    print(f"  Loaded {len(teloscope_df)} telomere entries")

    print("\nLoading Telogator2 HiFi data...")
    telogator_hifi_df = load_telogator2_data(TELOGATOR2_HIFI)
    print(f"  Loaded {len(telogator_hifi_df)} chromosome arm entries")

    print("\nLoading Telogator2 ONT data...")
    telogator_ont_df = load_telogator2_data(TELOGATOR2_ONT)
    print(f"  Loaded {len(telogator_ont_df)} chromosome arm entries")

    # Define plot configurations
    # 3 percentiles x 2 metrics x 2 correlations = 12 plots
    percentiles = [
        ("TL_p75", "TL_p75_plus_tvr", "p75"),
        ("TL_p90", "TL_p90_plus_tvr", "p90"),
        ("TL_max", "TL_max_plus_tvr", "max"),
    ]
    correlations = ["pearson", "spearman"]

    # Collect all results for summary
    results = []

    print("\n" + "=" * 60)
    print("Generating plots...")
    print("=" * 60)

    datasets = [
        ("hifi", "HiFi", telogator_hifi_df),
        ("ont", "ONT", telogator_ont_df),
    ]

    for dataset_prefix, dataset_label, telogator_df in datasets:
        for tl_col, tl_tvr_col, percentile_name in percentiles:
            for corr_type in correlations:
                corr_suffix = "pearson" if corr_type == "pearson" else "spearman"

                # TL only
                output_name = f"{dataset_prefix}_{percentile_name}_TL_{corr_suffix}.png"
                y_label = f"Telogator2 {dataset_label} {percentile_name.upper()} (bp)"
                result = create_scatter_plot(
                    teloscope_df, telogator_df, tl_col,
                    OUTPUT_DIR / output_name, y_label, corr_type
                )
                if result:
                    result["metric"] = f"TL_{percentile_name}"
                    result["dataset"] = dataset_label
                    results.append(result)

                # TL + TVR
                output_name = f"{dataset_prefix}_{percentile_name}_TL+TVR_{corr_suffix}.png"
                y_label = f"Telogator2 {dataset_label} {percentile_name.upper()}+TVR (bp)"
                result = create_scatter_plot(
                    teloscope_df, telogator_df, tl_tvr_col,
                    OUTPUT_DIR / output_name, y_label, corr_type
                )
                if result:
                    result["metric"] = f"TL+TVR_{percentile_name}"
                    result["dataset"] = dataset_label
                    results.append(result)

    # Print summary statistics
    print("\n" + "=" * 70)
    print("CORRELATION SUMMARY")
    print("=" * 70)
    print(f"{'Dataset':<8} {'Metric':<16} {'Type':<10} {'r/ρ':>8} {'p-value':>12} {'n':>5}")
    print("-" * 70)

    for res in results:
        print(f"{res['dataset']:<8} {res['metric']:<16} {res['type']:<10} {res['r']:>8.3f} {res['p']:>12.2e} {res['n']:>5}")

    print("-" * 70)
    print(f"\nAll plots saved to: {OUTPUT_DIR}")
    print(f"Total plots generated: {len(results)}")


if __name__ == "__main__":
    main()