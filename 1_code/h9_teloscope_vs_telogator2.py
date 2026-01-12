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

# We have to compare telomere lengths from Teloscope and Telogator2 on a scatter plots.
# Files used:
# 1. Teloscope:
# h9_annotation/2_data/2.2_processed/teloscope/H9_T2T_v0.1_hap1.fasta_terminal_telomeres.bed
# h9_annotation/2_data/2.2_processed/teloscope/H9_T2T_v0.1_hap2.fasta_terminal_telomeres.bed

# Row format: chr, start, end, length, label, fwdCounts, revCounts, canCounts, nonCanCounts, chrSize
# Row example: chr1_hap1	5	4540	4535	p	753	0	731	22	250382078
# Note: We need to extract the chrNum from chr so by adding the label it's comparable to Telogator2 nomenclature.

# 2. Telogator2:
# h9_annotation/2_data/2.2_processed/26.01.07_telogator2_hifi_default/tlens_by_allele.tsv
# h9_annotation/2_data/2.2_processed/26.01.08_telogator2_ont_default/tlens_by_allele.tsv

# Columns of interest: #chr, TL_p75, tvr_len
# Example of #chr: chr1p (equivalent to chrNum and label in Teloscope)

# Plot 1: Teloscope vs Telogator2 HiFi
# 1A: Teloscope length vs Telogator2 HiFi "TL_p75"
# 1B: Teloscope length vs Telogator2 HiFi "TL_p75"+"tvr_len"
# Plot 2: Teloscope vs Telogator2 ONT
# 2A: Teloscope length vs Telogator2 ONT "TL_p75"
# 2B: Teloscope length vs Telogator2 ONT "TL_p75"+"tvr_len"
# Save plots at: h9_annotation/3_figures/3.1_draft/26.01.09_teloscope_vs_telogator2

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy import stats

# Define paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "2_data" / "2.2_processed"
OUTPUT_DIR = BASE_DIR / "3_figures" / "3.1_draft" / "26.01.09_teloscope_vs_telogator2"
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


def load_telogator2_data(filepath: Path) -> pd.DataFrame:
    """Load Telogator2 data and aggregate by chromosome arm."""
    df = pd.read_csv(filepath, sep="\t")

    # Group by chromosome arm and get median TL_p75 and tvr_len
    agg_df = df.groupby("#chr").agg({
        "TL_p75": "median",
        "tvr_len": "median"
    }).reset_index()

    agg_df["TL_p75_plus_tvr"] = agg_df["TL_p75"] + agg_df["tvr_len"]
    agg_df = agg_df.rename(columns={"#chr": "telogator_key"})

    return agg_df


def create_scatter_plot(teloscope_df: pd.DataFrame, telogator_df: pd.DataFrame,
                        telogator_col: str, title: str, output_path: Path,
                        y_label: str):
    """Create a scatter plot comparing Teloscope vs Telogator2."""
    # Merge data
    merged = teloscope_df.merge(telogator_df, left_on="teloscope_key",
                                 right_on="telogator_key", how="inner")

    if merged.empty:
        print(f"Warning: No matching data for {title}")
        return

    x = merged["teloscope_length"]
    y = merged[telogator_col]

    # Calculate correlation
    r, p_value = stats.pearsonr(x, y)

    # Calculate linear regression
    slope, intercept, _, _, _ = stats.linregress(x, y)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(x, y, alpha=0.7, edgecolors="black", linewidth=0.5, s=60)

    # Add regression line
    x_line = np.array([x.min(), x.max()])
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, "r--", linewidth=1.5, label=f"Linear fit")

    # Add identity line (y=x)
    max_val = max(x.max(), y.max())
    min_val = min(x.min(), y.min())
    ax.plot([min_val, max_val], [min_val, max_val], "k:", linewidth=1, alpha=0.5, label="y = x")

    # Add labels for each point
    for _, row in merged.iterrows():
        ax.annotate(row["teloscope_key"].replace("chr", ""),
                   (row["teloscope_length"], row[telogator_col]),
                   fontsize=7, alpha=0.7, ha="left", va="bottom")

    ax.set_xlabel("Teloscope Length (bp)", fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(f"{title}\n(r = {r:.3f}, p = {p_value:.2e}, n = {len(merged)})", fontsize=12)
    ax.legend(loc="upper left")

    # Equal aspect ratio
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output_path}")
    print(f"  Correlation: r = {r:.3f}, p = {p_value:.2e}")
    print(f"  N points: {len(merged)}")


def main():
    print("Loading Teloscope data...")
    teloscope_df = load_teloscope_data(TELOSCOPE_HAP1, TELOSCOPE_HAP2)
    print(f"  Loaded {len(teloscope_df)} telomere entries")

    print("\nLoading Telogator2 HiFi data...")
    telogator_hifi_df = load_telogator2_data(TELOGATOR2_HIFI)
    print(f"  Loaded {len(telogator_hifi_df)} chromosome arm entries")

    print("\nLoading Telogator2 ONT data...")
    telogator_ont_df = load_telogator2_data(TELOGATOR2_ONT)
    print(f"  Loaded {len(telogator_ont_df)} chromosome arm entries")

    # Plot 1A: Teloscope vs Telogator2 HiFi TL_p75
    print("\nCreating Plot 1A: Teloscope vs Telogator2 HiFi TL_p75...")
    create_scatter_plot(
        teloscope_df, telogator_hifi_df, "TL_p75",
        "Teloscope vs Telogator2 HiFi (TL_p75)",
        OUTPUT_DIR / "1A_teloscope_vs_telogator2_hifi_TL_p75.png",
        "Telogator2 HiFi TL_p75 (bp)"
    )

    # Plot 1B: Teloscope vs Telogator2 HiFi TL_p75 + tvr_len
    print("\nCreating Plot 1B: Teloscope vs Telogator2 HiFi TL_p75 + tvr_len...")
    create_scatter_plot(
        teloscope_df, telogator_hifi_df, "TL_p75_plus_tvr",
        "Teloscope vs Telogator2 HiFi (TL_p75 + TVR)",
        OUTPUT_DIR / "1B_teloscope_vs_telogator2_hifi_TL_p75_plus_tvr.png",
        "Telogator2 HiFi TL_p75 + TVR (bp)"
    )

    # Plot 2A: Teloscope vs Telogator2 ONT TL_p75
    print("\nCreating Plot 2A: Teloscope vs Telogator2 ONT TL_p75...")
    create_scatter_plot(
        teloscope_df, telogator_ont_df, "TL_p75",
        "Teloscope vs Telogator2 ONT (TL_p75)",
        OUTPUT_DIR / "2A_teloscope_vs_telogator2_ont_TL_p75.png",
        "Telogator2 ONT TL_p75 (bp)"
    )

    # Plot 2B: Teloscope vs Telogator2 ONT TL_p75 + tvr_len
    print("\nCreating Plot 2B: Teloscope vs Telogator2 ONT TL_p75 + tvr_len...")
    create_scatter_plot(
        teloscope_df, telogator_ont_df, "TL_p75_plus_tvr",
        "Teloscope vs Telogator2 ONT (TL_p75 + TVR)",
        OUTPUT_DIR / "2B_teloscope_vs_telogator2_ont_TL_p75_plus_tvr.png",
        "Telogator2 ONT TL_p75 + TVR (bp)"
    )

    print(f"\nAll plots saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()