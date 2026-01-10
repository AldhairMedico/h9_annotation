#!/usr/bin/env python3
"""
Generate Contig Nx and Scaffold Nx plots from gfastats output.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from typing import Dict, List, Tuple

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
DATA_DIR = os.path.join(REPO_ROOT, "2_data", "2.2_processed", "gfastats")
OUT_DIR = os.path.join(REPO_ROOT, "3_figures", "3.1_draft")

# Palette (user-updated; do not change)
PALETTE: Dict[str, str] = {
    "GRCh38"    : "#999999",  # Neutral Gray
    "CHM13"     : "#F0E442",  # Yellow
    "YAO mat"   : "#FF5353",  # Red
    "YAO pat"   : "#FFA5A5",  # Light Red
    "HG002 mat" : "#0072B2",  # Okabe-Ito Blue
    "HG002 pat" : "#56B4E9",  # Okabe-Ito Sky Blue
    "RPE1 hap1" : "#984EA3",  # Purple
    "RPE1 hap2" : "#CC79A7",  # Magenta
    "I002C hap1": "#D55E00",  # Vermillion
    "I002C hap2": "#E69F00",  # Orange
    "H9 hap1"   : "#00796B",  # Teal 700
    "H9 hap2"   : "#80CBC4",  # Teal 200
}

# Map raw assembly names to display names for palette lookup
ASSEMBLY_MAP: Dict[str, str] = {
    "GCA_000001405.29_GRCh38.p14_genomic.chr.fna": "GRCh38",
    "GCA_009914755.4_T2T-CHM13v2.0_genomic.chr.fna": "CHM13",
    "GCA_018852605.3_hg002v1.1.pat_genomic.fna": "HG002 pat",
    "GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna": "HG002 mat",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.chr.fna": "RPE1 hap2",
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.chr.fna": "RPE1 hap1",
    "GWHDOOG00000000.genome.chr.fasta.gz": "YAO pat",
    "GWHDQZJ00000000.genome.chr.fasta.gz": "YAO mat",
    "GWHGEYB00000000.1.genome.fasta.gz": "YAO pat",
    "GWHGEYC00000000.1.genome.fasta.gz": "YAO mat",
    "H9_T2T_v0.1_hap1.fasta": "H9 hap1",
    "H9_T2T_v0.1_hap2.fasta": "H9 hap2",
    "I002Cv0.7.hap1.fasta.gz": "I002C hap1",
    "I002Cv0.7.hap2.fasta.gz": "I002C hap2",
    # Non-.chr versions (for quality assessment including unplaced contigs)
    "GCA_000001405.29_GRCh38.p14_genomic": "GRCh38",
    "GCA_009914755.4_T2T-CHM13v2.0_genomic": "CHM13",
    "GCA_018852605.3_hg002v1.1.pat_genomic": "HG002 pat",
    "GCA_018852615.3_hg002v1.1.mat_genomic": "HG002 mat",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic": "RPE1 hap2",
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic": "RPE1 hap1",
    "GWHDOOG00000000.genome": "YAO pat",
    "GWHDQZJ00000000.genome": "YAO mat",
    "GWHGEYB00000000.1.genome": "YAO pat",
    "GWHGEYC00000000.1.genome": "YAO mat",
    "H9_T2T_v0.1_hap1": "H9 hap1",
    "H9_T2T_v0.1_hap2": "H9 hap2",
    "I002Cv0.7.hap1": "I002C hap1",
    "I002Cv0.7.hap2": "I002C hap2",
}


def parse_nx_file(filepath: str) -> Tuple[str, List[float], List[float]]:
    """
    Parse an Nx TSV file.
    Format: assembly_name<TAB>size1,cumfrac1<TAB>size2,cumfrac2<TAB>...

    Returns: (assembly_name, nx_percentages, lengths)
    """
    with open(filepath, 'r') as f:
        line = f.readline().strip()

    parts = line.split('\t')
    asm_name = parts[0]

    nx_percentages = []
    lengths = []

    for entry in parts[1:]:
        if ',' in entry:
            size_str, cumfrac_str = entry.split(',')
            length = float(size_str)
            cumfrac = float(cumfrac_str) * 100  # Convert to percentage
            lengths.append(length)
            nx_percentages.append(cumfrac)

    return asm_name, nx_percentages, lengths


def load_nx_data(nx_type: str) -> Dict[str, Tuple[List[float], List[float]]]:
    """
    Load all Nx data from individual assembly directories.
    nx_type: 'Contig' or 'Scaffold'

    Returns: {display_name: (nx_percentages, lengths)}
    """
    data = {}

    if not os.path.exists(DATA_DIR):
        print(f"Data directory not found: {DATA_DIR}")
        return data

    for asm_dir in os.listdir(DATA_DIR):
        asm_path = os.path.join(DATA_DIR, asm_dir)
        if not os.path.isdir(asm_path):
            continue

        nx_file = os.path.join(asm_path, f"gfastatsNx{nx_type}_{asm_dir}.tsv")
        if not os.path.exists(nx_file):
            continue

        try:
            asm_name, nx_pct, lengths = parse_nx_file(nx_file)

            # Map to display name
            display_name = ASSEMBLY_MAP.get(asm_dir, asm_dir)

            # Skip if we already have this display name (prefer non-.chr versions)
            if display_name in data:
                # If current is non-.chr version, prefer it
                if ".chr" not in asm_dir:
                    data[display_name] = (nx_pct, lengths)
            else:
                data[display_name] = (nx_pct, lengths)

        except Exception as e:
            print(f"Error parsing {nx_file}: {e}")

    return data


def plot_nx(data: Dict[str, Tuple[List[float], List[float]]],
            nx_type: str,
            output_path: str) -> None:
    """
    Create an Nx plot.

    Args:
        data: {display_name: (nx_percentages, lengths)}
        nx_type: 'Contig' or 'Scaffold'
        output_path: Path to save the figure
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    # Sort by palette order for consistent legend
    palette_order = list(PALETTE.keys())
    sorted_names = sorted(data.keys(),
                          key=lambda x: palette_order.index(x) if x in palette_order else len(palette_order))

    for name in sorted_names:
        nx_pct, lengths = data[name]
        color = PALETTE.get(name, "#8D8D8D")
        ax.plot(nx_pct, lengths, label=name, color=color, linewidth=2, rasterized=True)

    ax.set_xlabel("Nx (%)", fontsize=12)
    ax.set_ylabel(f"{nx_type} length (bp)", fontsize=12)
    ax.set_yscale('log')
    ax.set_xlim(0, 100)
    ax.set_title(f"{nx_type} Nx Plot", fontsize=14)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load and plot Contig Nx
    contig_data = load_nx_data("Contig")
    if contig_data:
        plot_nx(contig_data, "Contig", os.path.join(OUT_DIR, "gfastats_Nx_Contig.png"))
    else:
        print("No Contig Nx data found")

    # Load and plot Scaffold Nx
    scaffold_data = load_nx_data("Scaffold")
    if scaffold_data:
        plot_nx(scaffold_data, "Scaffold", os.path.join(OUT_DIR, "gfastats_Nx_Scaffold.png"))
    else:
        print("No Scaffold Nx data found")


if __name__ == "__main__":
    main()
