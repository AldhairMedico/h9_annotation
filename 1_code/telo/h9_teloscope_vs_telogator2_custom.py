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
import matplotlib.lines as mlines

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
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "2_data" / "2.2_processed"
OUTPUT_DIR = BASE_DIR / "3_figures" / "3.1_draft" / "26.02.23_teloscope_vs_telogator2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Teloscope files
TELOSCOPE_HAP1 = DATA_DIR / "25.12.10_teloscope_asm" / "H9_T2T_v0.1_hap1.fasta_terminal_telomeres.bed"
TELOSCOPE_HAP2 = DATA_DIR / "25.12.10_teloscope_asm" / "H9_T2T_v0.1_hap2.fasta_terminal_telomeres.bed"

# Telogator2 files
TELOGATOR2_HIFI = DATA_DIR / "26.02.19_telogator2_hifi_n4_fixed" / "tlens_by_allele.tsv"
TELOGATOR2_ONT = DATA_DIR / "26.02.19_telogator2_ont_n4_fixed" / "tlens_by_allele.tsv"


def load_teloscope_data(hap1_path: Path, hap2_path: Path) -> pd.DataFrame:
    """Load and combine Teloscope data from both haplotypes."""
    cols = ["chr", "start", "end", "length", "label", "fwdCounts", "revCounts",
            "canCounts", "nonCanCounts", "chrSize"]

    df1 = pd.read_csv(hap1_path, sep="\t", header=None, names=cols)
    df1["hap"] = "hap1"
    df2 = pd.read_csv(hap2_path, sep="\t", header=None, names=cols)
    df2["hap"] = "hap2"

    df = pd.concat([df1, df2], ignore_index=True)

    # Extract chromosome number and create a key matching Telogator2 format
    # e.g., chr1p_hap1 (includes haplotype to match ref_samp in Telogator2)
    df["chr_num"] = df["chr"].str.extract(r"(chr\d+|chrX|chrY)")[0]
    df["teloscope_key"] = df["chr_num"] + df["label"] + "_" + df["hap"]

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

    # Extract haplotype from ref_samp (e.g., "H9-hap1" -> "hap1")
    df["hap"] = df["ref_samp"].str.extract(r"(hap\d+)")[0]
    df["telogator_key"] = df["#chr"] + "_" + df["hap"]

    # Parse read_TLs and compute percentiles for each row
    df["read_TLs_parsed"] = df["read_TLs"].apply(parse_read_tls)
    df["TL_p75"] = df["read_TLs_parsed"].apply(lambda x: np.percentile(x, 75))
    df["TL_p90"] = df["read_TLs_parsed"].apply(lambda x: np.percentile(x, 90))
    df["TL_max"] = df["read_TLs_parsed"].apply(lambda x: max(x))

    # Aggregate by chromosome arm + haplotype (median across multiple alleles)
    agg_df = df.groupby("telogator_key").agg({
        "TL_p75": "median",
        "TL_p90": "median",
        "TL_max": "median",
        "tvr_len": "median"
    }).reset_index()

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

    # Extract haplotype, arm, and chromosome label from teloscope_key
    # Key format: "chr1p_hap1", "chrXq_hap2", etc.
    merged = merged.copy()
    merged["chr_label"] = merged["teloscope_key"].str.extract(r"chr(\d+|X|Y)")[0]
    merged["arm"]       = merged["teloscope_key"].str.extract(r"chr(?:\d+|X|Y)([pq])")[0]
    merged["hap"]       = merged["teloscope_key"].str.extract(r"(hap\d+)")[0]

    # 1.1 Palette by haplotype
    palette = {"hap1": "#00796B", "hap2": "#80CBC4"}
    # 1.2 Marker by arm: p = square, q = circle
    markers = {"p": "s", "q": "o"}

    # Create figure
    fig, ax = plt.subplots(figsize=(4.5, 4.5))

    # Plot each (haplotype × arm) group with the correct colour and marker
    for hap, color in palette.items():
        for arm, marker in markers.items():
            mask = (merged["hap"] == hap) & (merged["arm"] == arm)
            subset = merged[mask]
            if subset.empty:
                continue
            ax.scatter(subset["teloscope_length"], subset[telogator_col],
                       s=40, c=color, marker=marker, alpha=0.85,
                       edgecolors="white", linewidth=0.4, zorder=3)

    # 1.3 Chromosome number labels — tiny, non-overlapping
    texts = []
    for _, row in merged.iterrows():
        if pd.notna(row["chr_label"]):
            txt = ax.text(row["teloscope_length"], row[telogator_col],
                          str(row["chr_label"]), fontsize=5.5,
                          ha="center", va="center", color="#222222", zorder=4)
            texts.append(txt)

    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=0.4))
    except ImportError:
        pass  # Labels may overlap slightly without adjustText

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

    # Legend: haplotype colours + arm shapes
    legend_handles = [
        mlines.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="#00796B", markersize=7, label="H9 hap1"),
        mlines.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="#80CBC4", markersize=7, label="H9 hap2"),
        mlines.Line2D([0], [0], marker="s", color="w",
                      markerfacecolor="#777777", markersize=7, label="p arm"),
        mlines.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="#777777", markersize=7, label="q arm"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="upper left",
              framealpha=0.9, edgecolor="#CCCCCC")

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


def create_hifi_vs_ont_scatter(hifi_df: pd.DataFrame, ont_df: pd.DataFrame,
                                col: str, output_path: Path,
                                axis_label: str,
                                correlation_type: str = "pearson"):
    """Scatter plot comparing Telogator2 HiFi vs ONT for a given metric.

    Points are coloured by haplotype and shaped by arm (p/q),
    with chromosome labels, same style as the Teloscope plots.
    """
    merged = hifi_df.merge(ont_df, on="telogator_key", how="inner",
                           suffixes=("_hifi", "_ont"))

    if merged.empty:
        print(f"Warning: No matching data for {output_path.name}")
        return None

    x = merged[f"{col}_hifi"].values
    y = merged[f"{col}_ont"].values

    if correlation_type == "pearson":
        r, p_value = stats.pearsonr(x, y)
        corr_label = "r"
    else:
        r, p_value = stats.spearmanr(x, y)
        corr_label = "ρ"

    slope, intercept, _, _, _ = stats.linregress(x, y)

    merged = merged.copy()
    merged["chr_label"] = merged["telogator_key"].str.extract(r"chr(\d+|X|Y)")[0]
    merged["arm"]       = merged["telogator_key"].str.extract(r"chr(?:\d+|X|Y)([pq])")[0]
    merged["hap"]       = merged["telogator_key"].str.extract(r"(hap\d+)")[0]

    palette = {"hap1": "#00796B", "hap2": "#80CBC4"}
    markers = {"p": "s", "q": "o"}

    fig, ax = plt.subplots(figsize=(4.5, 4.5))

    for hap, color in palette.items():
        for arm, marker in markers.items():
            mask = (merged["hap"] == hap) & (merged["arm"] == arm)
            subset = merged[mask]
            if subset.empty:
                continue
            ax.scatter(subset[f"{col}_hifi"], subset[f"{col}_ont"],
                       s=40, c=color, marker=marker, alpha=0.85,
                       edgecolors="white", linewidth=0.4, zorder=3)

    texts = []
    for _, row in merged.iterrows():
        if pd.notna(row["chr_label"]):
            txt = ax.text(row[f"{col}_hifi"], row[f"{col}_ont"],
                          str(row["chr_label"]), fontsize=5.5,
                          ha="center", va="center", color="#222222", zorder=4)
            texts.append(txt)

    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=0.4))
    except ImportError:
        pass

    x_range = np.array([x.min(), x.max()])
    ax.plot(x_range, slope * x_range + intercept, color="#E74C3C",
            linewidth=1.5, linestyle="-", zorder=2)

    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lims, lims, color="#95A5A6", linewidth=1, linestyle="--",
            alpha=0.8, zorder=1)

    ax.set_xlabel(f"Telogator2 HiFi {axis_label} (bp)", fontsize=11, fontweight="medium")
    ax.set_ylabel(f"Telogator2 ONT {axis_label} (bp)", fontsize=11, fontweight="medium")

    legend_handles = [
        mlines.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="#00796B", markersize=7, label="H9 hap1"),
        mlines.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="#80CBC4", markersize=7, label="H9 hap2"),
        mlines.Line2D([0], [0], marker="s", color="w",
                      markerfacecolor="#777777", markersize=7, label="p arm"),
        mlines.Line2D([0], [0], marker="o", color="w",
                      markerfacecolor="#777777", markersize=7, label="q arm"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="upper left",
              framealpha=0.9, edgecolor="#CCCCCC")

    corr_text = f"{corr_label} = {r:.2f}\np = {p_value:.1e}\nn = {len(merged)}"
    ax.text(0.97, 0.03, corr_text, transform=ax.transAxes, fontsize=9,
            verticalalignment="bottom", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#BDC3C7", alpha=0.9))

    ax.tick_params(axis="both", which="major", labelsize=9)
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

    # HiFi vs ONT comparison
    print("\n" + "=" * 60)
    print("Generating HiFi vs ONT plots...")
    print("=" * 60)

    for tl_col, tl_tvr_col, percentile_name in percentiles:
        for corr_type in correlations:
            corr_suffix = "pearson" if corr_type == "pearson" else "spearman"

            # TL only
            output_name = f"hifi_vs_ont_{percentile_name}_TL_{corr_suffix}.png"
            axis_label = percentile_name.upper()
            result = create_hifi_vs_ont_scatter(
                telogator_hifi_df, telogator_ont_df, tl_col,
                OUTPUT_DIR / output_name, axis_label, corr_type
            )
            if result:
                result["metric"] = f"TL_{percentile_name}"
                result["dataset"] = "HiFi vs ONT"
                results.append(result)

            # TL + TVR
            output_name = f"hifi_vs_ont_{percentile_name}_TL+TVR_{corr_suffix}.png"
            axis_label = f"{percentile_name.upper()}+TVR"
            result = create_hifi_vs_ont_scatter(
                telogator_hifi_df, telogator_ont_df, tl_tvr_col,
                OUTPUT_DIR / output_name, axis_label, corr_type
            )
            if result:
                result["metric"] = f"TL+TVR_{percentile_name}"
                result["dataset"] = "HiFi vs ONT"
                results.append(result)

    # Print summary statistics
    print("\n" + "=" * 70)
    print("CORRELATION SUMMARY")
    print("=" * 70)
    print(f"{'Dataset':<12} {'Metric':<16} {'Type':<10} {'r/ρ':>8} {'p-value':>12} {'n':>5}")
    print("-" * 70)

    for res in results:
        print(f"{res['dataset']:<12} {res['metric']:<16} {res['type']:<10} {res['r']:>8.3f} {res['p']:>12.2e} {res['n']:>5}")

    print("-" * 70)
    print(f"\nAll plots saved to: {OUTPUT_DIR}")
    print(f"Total plots generated: {len(results)}")


if __name__ == "__main__":
    main()