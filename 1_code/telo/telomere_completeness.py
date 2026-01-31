# telomere_completeness.py
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

# ──────────────────────────────────────────────────────────────────────────────
# Function to annotate generic stacked horizontal bars with their widths (integers)
def annotate_bars(ax):
    for patch in ax.patches:
        w = patch.get_width()
        if w > 0:
            ax.text(
                patch.get_x() + w/2,
                patch.get_y() + patch.get_height()/2,
                f"{int(round(w))}",
                ha='center', va='center',
                fontsize=8, color='black'
            )

# Annotate stacked horizontal bars with custom (absolute) labels for each patch
def annotate_bars_with_values(ax, values):
    """
    values: list/iterable of numbers (absolute counts) in the same order as ax.patches
    """
    for patch, val in zip(ax.patches, values):
        w = patch.get_width()
        if w > 0:
            ax.text(
                patch.get_x() + w/2,
                patch.get_y() + patch.get_height()/2,
                f"{int(round(val))}",
                ha='center', va='center',
                fontsize=8, color='black'
            )

# ──────────────────────────────────────────────────────────────────────────────
# Paths - relative to repo root
# Get the repo root directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

# Use assembly metrics TSV (rows=metric, cols=assembly)
reports_file = os.path.join(REPO_ROOT, "2_data", "2.2_processed", "25.12.10_asms_x1_TTAGGG_v1.3.assembly_metrics.tsv")

# Output dir and formats
plot_dir = os.path.join(REPO_ROOT, "3_figures", "3.1_draft" ,"26.01.29_telomeres")
os.makedirs(plot_dir, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Global style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.weight': 'regular',
    'font.size': 9,
    'axes.linewidth': 0.5,
    'axes.titlesize': 9,
    'axes.labelsize': 9,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'legend.frameon': False,
    'figure.dpi': 600,
    'savefig.dpi': 600
})

# ──────────────────────────────────────────────────────────────────────────────
# Load data
# TSV → rows per metric, columns per assembly. We transpose and
#        derive the same columns expected downstream.
raw = pd.read_csv(reports_file, sep="\t").fillna(0)
raw.set_index('metric', inplace=True)
df = raw.T

# Map metric rows → columns used downstream
df['total_paths']          = df['Total paths']
df['total_telomeres']      = df['Total telomeres']
df['two_telomeres']        = df['Two telomeres']
df['one_telomere']         = df['One telomere']
df['zero_telomeres']       = df['Zero telomeres']

df['t2t']                  = df['T2T']
df['gapped_t2t']           = df['Gapped T2T']
df['missassembled']        = df['Missassembled']
df['gapped_missassembled'] = df['Gapped missassembled']
df['incomplete']           = df['Incomplete']
df['gapped_incomplete']    = df['Gapped incomplete']
df['discordant']           = df['Discordant']
df['gapped_discordant']    = df['Gapped discordant']
df['no_telomeres']         = df['No telomeres']
df['gapped_no_telomeres']  = df['Gapped no telomeres']

# Optional manual list of accessions that INCLUDE mitochondrial DNA.
# If provided, we subtract 1 from: total_paths, zero_telomeres, no_telomeres
# for those accessions. If empty, auto-correct any row with no_telomeres>0.
mito_accessions = []  # e.g., ["GCA_048771995.1", "GCA_051427915.1"]

def _safe_sub_one(s: pd.Series) -> pd.Series:
    s = s.astype(int)
    return (s - 1).clip(lower=0)

if mito_accessions:
    mask = df.index.isin(mito_accessions)
    for col in ['total_paths', 'zero_telomeres', 'no_telomeres']:
        if col in df.columns:
            df.loc[mask, col] = _safe_sub_one(df.loc[mask, col])
else:
    # Auto rule: assembled CM → only mitochondrion lacks telomeres/gaps.
    # Subtract 1 where there is at least one "no_telomeres".
    if 'no_telomeres' in df.columns:
        mask = df['no_telomeres'].astype(int) > 0
        for col in ['total_paths', 'zero_telomeres', 'no_telomeres']:
            if col in df.columns:
                df.loc[mask, col] = _safe_sub_one(df.loc[mask, col])

# ──────────────────────────────────────────────────────────────────────────────
# Display order + labels for y-axis
# Based on assemblies in data/2_processed/25.12.10_asms_x1_TTAGGG_v1.3.assembly_metrics.tsv
DISPLAY_LABELS = {
    "GCA_009914755.4_T2T-CHM13v2.0_genomic.chr.fna": "CHM13",
    "GCA_018852605.3_hg002v1.1.pat_genomic.fna": "HG002 pat",
    "GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna": "HG002 mat",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.chr.fna": "RPE1 hap2",
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.chr.fna": "RPE1 hap1",
    "GWHGEYB00000000.1.genome.fasta.gz": "YAO pat",
    "GWHGEYC00000000.1.genome.fasta.gz": "YAO mat",
    "H9_T2T_v0.1_hap1.fasta": "H9 hap1",
    "H9_T2T_v0.1_hap2.fasta": "H9 hap2",
    "I002Cv0.7.hap1.fasta.gz": "I002C hap1",
    "I002Cv0.7.hap2.fasta.gz": "I002C hap2",
}

# Split into two groups: human assemblies and human cell lines
human_assemblies = [
    "GCA_009914755.4_T2T-CHM13v2.0_genomic.chr.fna",
    "GWHGEYB00000000.1.genome.fasta.gz",
    "GWHGEYC00000000.1.genome.fasta.gz",
    "GCA_018852605.3_hg002v1.1.pat_genomic.fna",
    "GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna",
    "I002Cv0.7.hap1.fasta.gz",
    "I002Cv0.7.hap2.fasta.gz",
]

cell_lines = [
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.chr.fna",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.chr.fna",
    "H9_T2T_v0.1_hap1.fasta",
    "H9_T2T_v0.1_hap2.fasta",
]

# Combined order for full dataframe
order = human_assemblies + cell_lines

# ──────────────────────────────────────────────────────────────────────────────
# GRCh38 special handling: none of the telomeres found by Teloscope were
# actually terminal (>10kbp from ends), so we force zeros for GRCh38
GRCH38_ID = "GCA_000001405.29_GRCh38.p14_genomic.chr.fna"
# Keep only those present and reindex in desired order
order = [acc for acc in order if acc in df.index]
df = df.reindex(order)

# Override GRCh38 values: 0 found telomeres, 48 missing, 48 gapped_no_telomeres
if GRCH38_ID in df.index:
    df.loc[GRCH38_ID, 'total_telomeres'] = 0
    df.loc[GRCH38_ID, 'two_telomeres'] = 0
    df.loc[GRCH38_ID, 'one_telomere'] = 0
    df.loc[GRCH38_ID, 'zero_telomeres'] = 48
    # Reset all completeness categories to 0 except gapped_no_telomeres
    for col in ['t2t', 'gapped_t2t', 'incomplete', 'gapped_incomplete',
                'missassembled', 'gapped_missassembled', 'discordant',
                'gapped_discordant', 'no_telomeres']:
        df.loc[GRCH38_ID, col] = 0
    df.loc[GRCH38_ID, 'gapped_no_telomeres'] = 24

# Pretty y tick labels (no years, just assembly / group labels)
yticklabels = [DISPLAY_LABELS.get(acc, acc) for acc in df.index]

# Create separate dataframes for each group
human_asm_order = [acc for acc in human_assemblies if acc in df.index]
cell_line_order = [acc for acc in cell_lines if acc in df.index]
df_human = df.loc[human_asm_order]
df_cells = df.loc[cell_line_order]
yticklabels_human = [DISPLAY_LABELS.get(acc, acc) for acc in human_asm_order]
yticklabels_cells = [DISPLAY_LABELS.get(acc, acc) for acc in cell_line_order]

# ──────────────────────────────────────────────────────────────────────────────
# Color map (A uses Found + Missing; C uses full palette; "Missing" = #F0F0F0)
colors = {
    'Found':                   '#BDA7CA',  # Blue - contrasts with green completeness palette
    'Missing':                 '#F0F0F0',

    't2t':                     '#3C995B',
    'gapped_t2t':              '#9CCF60',
    'incomplete':              '#FFC754',
    'gapped_incomplete':       '#FFE885',
    'missassembled':           '#D6594C',
    'gapped_missassembled':    '#F58B6D',
    'discordant':              '#8278F4',
    'gapped_discordant':       '#B395EB',
    'no_telomeres':            '#C8C8C8',
    'gapped_no_telomeres':     '#F0F0F0',
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper to place legend at bottom with 2 rows
# LAYOUT: bbox_to_anchor=(0.5, Y) controls vertical position; more negative = lower
def bottom_legend(ax, handles, title=None, ncol=None):
    if ncol is None:
        ncol = (len(handles) + 1) // 2  # 2 rows
    ax.legend(handles=handles, title=title, frameon=False,
              loc='upper center', bbox_to_anchor=(0.5, -0.6),  # LAYOUT: Y = legend below x-axis label
              ncol=ncol, borderaxespad=0.)

# ──────────────────────────────────────────────────────────────────────────────
# MODE 1 — Absolute counts vs expected from assembled chromosomes
#   • Plot A: Found telomeres vs Expected–Found (expected = 2×paths)
#   • Plot C: Absolute category counts (no "Missing" category)
#   Layout: 2×2 grid (top=human assemblies, bottom=cell lines)

# Prepare data for both groups
def prepare_data_A(df_subset):
    p = df_subset['total_paths'].astype(int)
    tt = df_subset['total_telomeres'].astype(int)
    baseline = (2 * p).astype(int)
    missing = (baseline - tt).clip(lower=0).astype(int)
    return pd.DataFrame({'Found': tt, 'Missing': missing}, index=df_subset.index), tt, missing

def prepare_data_C(df_subset):
    return pd.DataFrame({
        't2t':                  df_subset['t2t'].astype(int),
        'gapped_t2t':           df_subset['gapped_t2t'].astype(int),
        'incomplete':           df_subset['incomplete'].astype(int),
        'gapped_incomplete':    df_subset['gapped_incomplete'].astype(int),
        'missassembled':        df_subset['missassembled'].astype(int),
        'gapped_missassembled': df_subset['gapped_missassembled'].astype(int),
        'discordant':           df_subset['discordant'].astype(int),
        'gapped_discordant':    df_subset['gapped_discordant'].astype(int),
        'no_telomeres':         df_subset['no_telomeres'].astype(int),
        'gapped_no_telomeres':  df_subset['gapped_no_telomeres'].astype(int),
    }, index=df_subset.index)

# Prepare data
dataA_human, tt_human, miss_human = prepare_data_A(df_human)
dataA_cells, tt_cells, miss_cells = prepare_data_A(df_cells)
dataC_human = prepare_data_C(df_human)
dataC_cells = prepare_data_C(df_cells)

# Combine for legend filtering (only show categories with data)
dataC_all = pd.concat([dataC_human, dataC_cells])

# Calculate height ratios based on number of assemblies in each group
n_human = len(human_asm_order)
n_cells = len(cell_line_order)
height_ratios = [n_human, n_cells]  # LAYOUT: proportional heights (7:4 for human:cells)

# ──────────────────────────────────────────────────────────────────────────────
# LAYOUT PARAMETERS (adjust these to fine-tune spacing):
#   - figsize=(W, H): overall figure dimensions in inches
#   - hspace: vertical gap between top and bottom rows (fraction of avg subplot height)
#   - wspace: horizontal gap between left and right columns
#   - bottom: space reserved at bottom for x-labels and legend (fraction of figure)
#   - width_ratios=[L, R]: relative widths of left vs right columns (e.g., [1, 2] = right is 2x wider)
# ──────────────────────────────────────────────────────────────────────────────
abs_fig = plt.figure(figsize=(6, 4))  # LAYOUT: (width, height) in inches
abs_fig.subplots_adjust(bottom=0.25, hspace=0.10)  # LAYOUT: bottom=space for legend, hspace=gap between rows
abs_gs = gridspec.GridSpec(2, 2, wspace=0.20, width_ratios=[1, 1.2], height_ratios=height_ratios)  # LAYOUT: width_ratios=[1,1.5], wspace=0.20

# Top row: human assemblies
axA_top = plt.subplot(abs_gs[0, 0])
axC_top = plt.subplot(abs_gs[0, 1])
# Bottom row: cell lines
axA_bot = plt.subplot(abs_gs[1, 0])
axC_bot = plt.subplot(abs_gs[1, 1])

# Set x-axis locators
for ax in [axA_top, axA_bot]:
    ax.xaxis.set_major_locator(MultipleLocator(20))
for ax in [axC_top, axC_bot]:
    ax.xaxis.set_major_locator(MultipleLocator(10))

# --- Panel A: Assembly telomeres ---
# Top (human assemblies)
dataA_human.plot.barh(stacked=True, ax=axA_top,
                      color=[colors['Found'], colors['Missing']],
                      width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axA_top.set_title('Assembly telomeres')
axA_top.set_yticklabels(yticklabels_human)
axA_top.set_ylabel('')
axA_top.set_xlabel('')
axA_top.set_xticklabels([])
axA_top.invert_yaxis()
abs_vals_A_human = list(tt_human.values) + list(miss_human.values)
annotate_bars_with_values(axA_top, abs_vals_A_human)

# Bottom (cell lines)
dataA_cells.plot.barh(stacked=True, ax=axA_bot,
                      color=[colors['Found'], colors['Missing']],
                      width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axA_bot.set_yticklabels(yticklabels_cells)
axA_bot.set_ylabel('')
axA_bot.set_xlabel('Total telomeres')
axA_bot.invert_yaxis()
abs_vals_A_cells = list(tt_cells.values) + list(miss_cells.values)
annotate_bars_with_values(axA_bot, abs_vals_A_cells)

# Legend for Panel A (below bottom row)
bottom_legend(axA_bot, [
    Patch(facecolor=colors['Found'],   label='Found'),
    Patch(facecolor=colors['Missing'], label='Missing'),
], ncol=1)  # LAYOUT: ncol=1 for 2 rows, 1 column

# --- Panel C: Chrs by telomere completeness ---
# Top (human assemblies)
dataC_human.plot.barh(stacked=True, ax=axC_top,
                      color=[colors[k] for k in dataC_human.columns],
                      width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axC_top.set_title('Chrs. by telomere completeness')
axC_top.set_yticklabels([])
axC_top.set_ylabel('')
axC_top.set_xlabel('')
axC_top.set_xticklabels([])
axC_top.invert_yaxis()
annotate_bars(axC_top)

# Bottom (cell lines)
dataC_cells.plot.barh(stacked=True, ax=axC_bot,
                      color=[colors[k] for k in dataC_cells.columns],
                      width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axC_bot.set_yticklabels([])
axC_bot.set_ylabel('')
axC_bot.set_xlabel('Assembled chromosomes')
axC_bot.invert_yaxis()
annotate_bars(axC_bot)

# Legend for Panel C (below bottom row, only categories with data)
bottom_legend(axC_bot, [Patch(facecolor=colors[k], label=k.replace('_',' '))
                        for k in dataC_all.columns if dataC_all[k].sum() > 0], ncol=2)

# Save MODE 1
for ext, dpi in [('pdf', None), ('svg', None), ('png', 600)]:
    save_kwargs = {'format': ext, 'bbox_inches': 'tight'}
    if dpi is not None:
        save_kwargs['dpi'] = dpi
    abs_fig.savefig(os.path.join(plot_dir, f"telomere_completeness_absolute.{ext}"), **save_kwargs)

# ──────────────────────────────────────────────────────────────────────────────
# MODE 2 — Percent of expected from assembled chromosomes
#   • Plot A: Observed% (= found / (2*paths) × 100), remainder = Expected − observed
#   • Plot C: Each category / paths × 100
#   Layout: 2×2 grid (top=human assemblies, bottom=cell lines)

def prepare_data_A_pct(df_subset):
    p = df_subset['total_paths'].astype(int)
    tt = df_subset['total_telomeres'].astype(int)
    baseline = (2 * p).astype(int)
    missing_abs = (baseline - tt).clip(lower=0).astype(int)
    with np.errstate(divide='ignore', invalid='ignore'):
        baseline_nonzero = baseline.replace(0, np.nan)
        found_pct = (tt / baseline_nonzero * 100).fillna(0).astype(float)
    missing_pct = (100 - found_pct).clip(lower=0)
    return pd.DataFrame({'Found': found_pct, 'Missing': missing_pct}, index=df_subset.index), tt, missing_abs

def prepare_data_C_pct(df_subset):
    p = df_subset['total_paths'].astype(int)
    with np.errstate(divide='ignore', invalid='ignore'):
        denom = p.replace(0, np.nan).astype(float)
        return pd.DataFrame({
            't2t':                  (df_subset['t2t'].astype(int) / denom * 100),
            'gapped_t2t':           (df_subset['gapped_t2t'].astype(int) / denom * 100),
            'incomplete':           (df_subset['incomplete'].astype(int) / denom * 100),
            'gapped_incomplete':    (df_subset['gapped_incomplete'].astype(int) / denom * 100),
            'missassembled':        (df_subset['missassembled'].astype(int) / denom * 100),
            'gapped_missassembled': (df_subset['gapped_missassembled'].astype(int) / denom * 100),
            'discordant':           (df_subset['discordant'].astype(int) / denom * 100),
            'gapped_discordant':    (df_subset['gapped_discordant'].astype(int) / denom * 100),
            'no_telomeres':         (df_subset['no_telomeres'].astype(int) / denom * 100),
            'gapped_no_telomeres':  (df_subset['gapped_no_telomeres'].astype(int) / denom * 100),
        }, index=df_subset.index).fillna(0.0)

# Prepare percent data
dataA_human_pct, tt_human_pct, miss_human_pct = prepare_data_A_pct(df_human)
dataA_cells_pct, tt_cells_pct, miss_cells_pct = prepare_data_A_pct(df_cells)
dataC_human_pct = prepare_data_C_pct(df_human)
dataC_cells_pct = prepare_data_C_pct(df_cells)

# Combine for legend filtering
dataC_all_pct = pd.concat([dataC_human_pct, dataC_cells_pct])

# ──────────────────────────────────────────────────────────────────────────────
# LAYOUT PARAMETERS (same as MODE 1):
#   - width_ratios=[1, 1.2]: left panel is 2/5, right panel is 3/5 of total width
# ──────────────────────────────────────────────────────────────────────────────
pct_fig = plt.figure(figsize=(6, 4))  # LAYOUT: (width, height) in inches
pct_fig.subplots_adjust(bottom=0.25, hspace=0.10)  # LAYOUT: bottom=space for legend, hspace=gap between rows
pct_gs = gridspec.GridSpec(2, 2, wspace=0.20, width_ratios=[1, 1.2], height_ratios=height_ratios)  # LAYOUT: width_ratios=[1,1.5], wspace=0.20

# Top row: human assemblies
axA_top_pct = plt.subplot(pct_gs[0, 0])
axC_top_pct = plt.subplot(pct_gs[0, 1])
# Bottom row: cell lines
axA_bot_pct = plt.subplot(pct_gs[1, 0])
axC_bot_pct = plt.subplot(pct_gs[1, 1])

# Set x-axis locators
for ax in [axA_top_pct, axA_bot_pct, axC_top_pct, axC_bot_pct]:
    ax.xaxis.set_major_locator(MultipleLocator(20))

# --- Panel A: Assembly telomeres (percent) ---
# Top (human assemblies)
dataA_human_pct.plot.barh(stacked=True, ax=axA_top_pct,
                          color=[colors['Found'], colors['Missing']],
                          width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axA_top_pct.set_title('Assembly telomeres')
axA_top_pct.set_yticklabels(yticklabels_human)
axA_top_pct.set_ylabel('')
axA_top_pct.set_xlabel('')
axA_top_pct.set_xticklabels([])
axA_top_pct.invert_yaxis()
abs_vals_A_human_pct = list(tt_human_pct.values) + list(miss_human_pct.values)
annotate_bars_with_values(axA_top_pct, abs_vals_A_human_pct)

# Bottom (cell lines)
dataA_cells_pct.plot.barh(stacked=True, ax=axA_bot_pct,
                          color=[colors['Found'], colors['Missing']],
                          width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axA_bot_pct.set_yticklabels(yticklabels_cells)
axA_bot_pct.set_ylabel('')
axA_bot_pct.set_xlabel('Total telomeres (%)')
axA_bot_pct.invert_yaxis()
abs_vals_A_cells_pct = list(tt_cells_pct.values) + list(miss_cells_pct.values)
annotate_bars_with_values(axA_bot_pct, abs_vals_A_cells_pct)

# Legend for Panel A (below bottom row)
bottom_legend(axA_bot_pct, [
    Patch(facecolor=colors['Found'],   label='Found'),
    Patch(facecolor=colors['Missing'], label='Missing'),
], ncol=1)  # LAYOUT: ncol=1 for 2 rows, 1 column

# --- Panel C: Chrs by telomere completeness (percent) ---
# Top (human assemblies)
dataC_human_pct.plot.barh(stacked=True, ax=axC_top_pct,
                          color=[colors[k] for k in dataC_human_pct.columns],
                          width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axC_top_pct.set_title('Chrs. by telomere completeness')
axC_top_pct.set_yticklabels([])
axC_top_pct.set_ylabel('')
axC_top_pct.set_xlabel('')
axC_top_pct.set_xticklabels([])
axC_top_pct.invert_yaxis()
# Annotate with absolute counts for human assemblies
abs_vals_C_human = []
for col in dataC_human_pct.columns:
    abs_vals_C_human.extend(list(df_human[col].astype(int).values))
annotate_bars_with_values(axC_top_pct, abs_vals_C_human)

# Bottom (cell lines)
dataC_cells_pct.plot.barh(stacked=True, ax=axC_bot_pct,
                          color=[colors[k] for k in dataC_cells_pct.columns],
                          width=0.7, edgecolor='white', linewidth=0.5, legend=False)
axC_bot_pct.set_yticklabels([])
axC_bot_pct.set_ylabel('')
axC_bot_pct.set_xlabel('Assembled chromosomes (%)')
axC_bot_pct.invert_yaxis()
# Annotate with absolute counts for cell lines
abs_vals_C_cells = []
for col in dataC_cells_pct.columns:
    abs_vals_C_cells.extend(list(df_cells[col].astype(int).values))
annotate_bars_with_values(axC_bot_pct, abs_vals_C_cells)

# Legend for Panel C (below bottom row, only categories with data)
bottom_legend(axC_bot_pct, [Patch(facecolor=colors[k], label=k.replace('_',' '))
                            for k in dataC_all_pct.columns if dataC_all_pct[k].sum() > 0], ncol=2)

# Save MODE 2
for ext, dpi in [('pdf', None), ('svg', None), ('png', 600)]:
    save_kwargs = {'format': ext, 'bbox_inches': 'tight'}
    if dpi is not None:
        save_kwargs['dpi'] = dpi
    pct_fig.savefig(os.path.join(plot_dir, f"telomere_completeness_percent.{ext}"), **save_kwargs)

print("Saved figures to", os.path.abspath(plot_dir))
