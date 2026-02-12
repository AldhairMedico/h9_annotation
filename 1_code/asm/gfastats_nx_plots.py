#!/usr/bin/env python3
"""
Generate Contig Nx and Scaffold Nx plots from gfastats output.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

sns.set_style("ticks")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
import os
from typing import Dict, List, Tuple

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "2_data", "2.2_processed", "26.01.09_gfastats")
OUT_DIR = os.path.join(REPO_ROOT, "3_figures", "3.1_draft", "26.02.11_nx_plots")

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
    "GCA_000001405.29_GRCh38.p14_genomic": "GRCh38",
    "GCA_009914755.4_T2T-CHM13v2.0_genomic": "CHM13",
    "GCA_018852605.3_hg002v1.1.pat_genomic": "HG002 pat",
    "GCA_018852615.3_hg002v1.1.mat_genomic": "HG002 mat",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic": "RPE1 hap2",
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic": "RPE1 hap1",
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

    for asm_dir, display_name in ASSEMBLY_MAP.items():
        nx_file = os.path.join(DATA_DIR, asm_dir, f"gfastatsNx{nx_type}_{asm_dir}.tsv")
        if not os.path.exists(nx_file):
            continue

        try:
            _, nx_pct, lengths = parse_nx_file(nx_file)
            data[display_name] = (nx_pct, lengths)
        except Exception as e:
            print(f"Error parsing {nx_file}: {e}")

    return data


def plot_nx(data: Dict[str, Tuple[List[float], List[float]]],
            nx_type: str,
            output_path: str,
            log_scale: bool = True) -> None:
    """
    Create an Nx plot.

    Args:
        data: {display_name: (nx_percentages, lengths)}
        nx_type: 'Contig' or 'Scaffold'
        output_path: Path to save the figure
        log_scale: If True, use log Y axis; otherwise linear with Mbp ticks
    """
    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    # Sort by palette order for consistent legend
    palette_order = list(PALETTE.keys())
    sorted_names = sorted(data.keys(),
                          key=lambda x: palette_order.index(x) if x in palette_order else len(palette_order))

    for name in sorted_names:
        nx_pct, lengths = data[name]
        color = PALETTE.get(name, "#8D8D8D")
        ax.plot(nx_pct, lengths, label=name, color=color, linewidth=1, rasterized=True)

    ax.set_xlabel("Nx (%)")
    ax.set_xlim(0, 105)

    if log_scale:
        ax.set_yscale('log')
        ax.set_ylabel(f"{nx_type} length (bp)")
    else:
        mbp = 1e6
        ax.set_yticks([0, 50*mbp, 100*mbp, 150*mbp, 200*mbp, 250*mbp])
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v/mbp)}"))
        ax.set_ylim(0, 260*mbp)
        ax.set_ylabel(f"{nx_type} length (Mbp)")

    # Legend: two-column, compact font; bottom-left for log, top-right for linear
    legend_loc = 'lower left' if log_scale else 'upper right'
    # Legend handles: small circles instead of long lines
    handles, labels = ax.get_legend_handles_labels()
    legend_handles = [plt.Line2D([0], [0], marker='o', color=h.get_color(),
                                 linestyle='', markersize=4) for h in handles]
    ax.legend(legend_handles, labels, loc=legend_loc, fontsize=6, ncol=2,
              columnspacing=0.25, handletextpad=0.3)

    sns.despine(ax=ax)

    # Axes padding: Increasing left/bottom or decreasing right/top shrinks the plot area
    fig.subplots_adjust(left=0.16, right=0.95, top=0.90, bottom=0.16)
    base, _ = os.path.splitext(output_path)
    for ext in (".png", ".pdf", ".svg"):
        path = base + ext
        plt.savefig(path, dpi=600, bbox_inches='tight')
        print(f"Saved: {path}")
    plt.close()


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    for nx_type in ("Contig", "Scaffold"):
        data = load_nx_data(nx_type)
        if not data:
            print(f"No {nx_type} Nx data found")
            continue
        plot_nx(data, nx_type, os.path.join(OUT_DIR, f"gfastats_Nx_{nx_type}_log.png"), log_scale=True)
        plot_nx(data, nx_type, os.path.join(OUT_DIR, f"gfastats_Nx_{nx_type}_linear.png"), log_scale=False)


if __name__ == "__main__":
    main()
