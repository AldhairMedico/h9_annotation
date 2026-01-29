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
# Display order + labels for x-axis
# Based on assemblies in data/2_processed/25.12.10_asms_x1_TTAGGG_v1.3.assembly_metrics.tsv
DISPLAY_LABELS = {
    # "GCA_000001405.29_GRCh38.p14_genomic.chr.fna": "GRCh38",
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

order = [
    # "GCA_000001405.29_GRCh38.p14_genomic.chr.fna",
    "GCA_009914755.4_T2T-CHM13v2.0_genomic.chr.fna",
    "GWHGEYB00000000.1.genome.fasta.gz",
    "GWHGEYC00000000.1.genome.fasta.gz",
    "GCA_018852605.3_hg002v1.1.pat_genomic.fna",
    "GCA_018852615.3_hg002v1.1.mat_genomic.chr.fna",
    "I002Cv0.7.hap1.fasta.gz",
    "I002Cv0.7.hap2.fasta.gz",
    "GCA_050656345.1_RPE1V1.1_Haplotype_1_genomic.chr.fna",
    "GCA_050656315.1_RPE1V1.1_Haplotype_2_genomic.chr.fna",
    "H9_T2T_v0.1_hap1.fasta",
    "H9_T2T_v0.1_hap2.fasta",
]

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

# ──────────────────────────────────────────────────────────────────────────────
# Series per accession (ensure integer dtype where applicable)
paths                 = df['total_paths'].astype(int)
total_telomeres       = df['total_telomeres'].astype(int)
two_telomeres         = df['two_telomeres'].astype(int)
one_telomere          = df['one_telomere'].astype(int)
zero_telomeres        = df['zero_telomeres'].astype(int)

# Completeness classes (full palette)
t2t                   = df['t2t'].astype(int)
gapped_t2t            = df['gapped_t2t'].astype(int)
missassembled         = df['missassembled'].astype(int)
gapped_missassembled  = df['gapped_missassembled'].astype(int)
incomplete            = df['incomplete'].astype(int)
gapped_incomplete     = df['gapped_incomplete'].astype(int)
discordant            = df.get('discordant', 0).astype(int)
gapped_discordant     = df.get('gapped_discordant', 0).astype(int)
no_telomeres          = df.get('no_telomeres', 0).astype(int)
gapped_no_telomeres   = df.get('gapped_no_telomeres', 0).astype(int)

# ──────────────────────────────────────────────────────────────────────────────
# Color map (A uses Found + Missing; C uses full palette; "Missing" = #F0F0F0)
colors = {
    'Found':                   '#4DCCBD',
    'Missing':                 '#F0F0F0',

    't2t':                     '#1A9641',
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
def bottom_legend(ax, handles, title=None, ncol=None):
    if ncol is None:
        ncol = (len(handles) + 1) // 2  # 2 rows
    ax.legend(handles=handles, title=title, frameon=False,
              loc='upper center', bbox_to_anchor=(0.5, -0.2),
              ncol=ncol, borderaxespad=0.)

# ──────────────────────────────────────────────────────────────────────────────
# MODE 1 — Absolute counts vs expected from assembled chromosomes
#   • Plot A: Found telomeres vs Expected–Found (expected = 2×paths)
#   • Plot C: Absolute category counts (no "Missing" category)
abs_fig = plt.figure(figsize=(8, 3))
abs_fig.subplots_adjust(bottom=0.05)
abs_gs  = gridspec.GridSpec(1, 2, wspace=0.15, width_ratios=[1, 1])
axA_abs = plt.subplot(abs_gs[0, 0])
axC_abs = plt.subplot(abs_gs[0, 1])

axA_abs.xaxis.set_major_locator(MultipleLocator(20))
axC_abs.xaxis.set_major_locator(MultipleLocator(10))

baseline_tel = (2 * paths).astype(int)
missing_tel_abs = (baseline_tel - total_telomeres).clip(lower=0).astype(int)

dataA_abs = pd.DataFrame({
    'Found':   total_telomeres,
    'Missing': missing_tel_abs
}, index=df.index)

dataA_abs.plot.barh(stacked=True, ax=axA_abs,
                    color=[colors['Found'], colors['Missing']],
                    width=0.7, edgecolor='white', linewidth=0.5)
axA_abs.set_title('Assembly telomeres')
axA_abs.set_xlabel('Total telomeres')
axA_abs.set_yticklabels(yticklabels)
axA_abs.set_ylabel('')
axA_abs.invert_yaxis()
bottom_legend(axA_abs, [
    Patch(facecolor=colors['Found'],   label='Found'),
    Patch(facecolor=colors['Missing'], label='Missing'),
], ncol=2)
# annotate with absolute numbers per segment
abs_vals_A = list(total_telomeres.values) + list(missing_tel_abs.values)
annotate_bars_with_values(axA_abs, abs_vals_A)

# Panel C (absolute)
dataC_abs = pd.DataFrame({
    't2t':                     t2t,
    'gapped_t2t':              gapped_t2t,
    'incomplete':              incomplete,
    'gapped_incomplete':       gapped_incomplete,
    'missassembled':           missassembled,
    'gapped_missassembled':    gapped_missassembled,
    'discordant':              discordant,
    'gapped_discordant':       gapped_discordant,
    'no_telomeres':            no_telomeres,
    'gapped_no_telomeres':     gapped_no_telomeres,
}, index=df.index)

dataC_abs.plot.barh(stacked=True, ax=axC_abs,
                    color=[colors[k] for k in dataC_abs.columns],
                    width=0.7, edgecolor='white', linewidth=0.5)
axC_abs.set_title('Chrs. by telomere completeness')
axC_abs.set_xlabel('Assembled chromosomes')
axC_abs.set_yticklabels([])  # shared y-axis labels with panel A
axC_abs.set_ylabel('')
axC_abs.invert_yaxis()
bottom_legend(axC_abs, [Patch(facecolor=colors[k], label=k.replace('_',' '))
                        for k in dataC_abs.columns if dataC_abs[k].sum() > 0], ncol=2)
annotate_bars(axC_abs)

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
pct_fig = plt.figure(figsize=(8, 3))
pct_fig.subplots_adjust(bottom=0.05)
pct_gs  = gridspec.GridSpec(1, 2, wspace=0.15, width_ratios=[1, 1])
axA_pct = plt.subplot(pct_gs[0, 0])
axC_pct = plt.subplot(pct_gs[0, 1])

# MODE 2 (percent)
for ax in (axA_pct, axC_pct):
    ax.xaxis.set_major_locator(MultipleLocator(20))

# A (percent)
with np.errstate(divide='ignore', invalid='ignore'):
    baseline_tel_nonzero = baseline_tel.replace(0, np.nan)
    found_pct = (total_telomeres / baseline_tel_nonzero * 100).fillna(0)
found_pct = found_pct.astype(float)
missing_pct = (100 - found_pct).clip(lower=0)

dataA_pct = pd.DataFrame({
    'Found':   found_pct,
    'Missing': missing_pct
}, index=df.index)

dataA_pct.plot.barh(stacked=True, ax=axA_pct,
                    color=[colors['Found'], colors['Missing']],
                    width=0.7, edgecolor='white', linewidth=0.5)
axA_pct.set_title('Assembly telomeres')
axA_pct.set_xlabel('Total telomeres %')
axA_pct.set_yticklabels(yticklabels)
axA_pct.set_ylabel('')
axA_pct.invert_yaxis()
# show absolute values in labels
abs_vals_Apct = list(total_telomeres.values) + list(missing_tel_abs.values)
annotate_bars_with_values(axA_pct, abs_vals_Apct)
bottom_legend(axA_pct, [
    Patch(facecolor=colors['Found'],   label='Found'),
    Patch(facecolor=colors['Missing'], label='Missing'),
], ncol=2)

# C (percent)
with np.errstate(divide='ignore', invalid='ignore'):
    denom = paths.replace(0, np.nan).astype(float)
    dataC_pct = pd.DataFrame({
        't2t':                   (t2t / denom * 100),
        'gapped_t2t':            (gapped_t2t / denom * 100),
        'incomplete':            (incomplete / denom * 100),
        'gapped_incomplete':     (gapped_incomplete / denom * 100),
        'missassembled':         (missassembled / denom * 100),
        'gapped_missassembled':  (gapped_missassembled / denom * 100),
        'discordant':            (discordant / denom * 100),
        'gapped_discordant':     (gapped_discordant / denom * 100),
        'no_telomeres':          (no_telomeres / denom * 100),
        'gapped_no_telomeres':   (gapped_no_telomeres / denom * 100),
    }, index=df.index).fillna(0.0)

dataC_pct.plot.barh(stacked=True, ax=axC_pct,
                    color=[colors[k] for k in dataC_pct.columns],
                    width=0.7, edgecolor='white', linewidth=0.5)
axC_pct.set_title('Chrs. by telomere completeness')
axC_pct.set_xlabel('Assembled chromosomes %')
axC_pct.set_yticklabels([])  # shared y-axis labels with panel A
axC_pct.set_ylabel('')
axC_pct.invert_yaxis()
# annotate with ABSOLUTE counts, not percents
abs_vals_Cpct = []
for col in ['t2t','gapped_t2t','incomplete','gapped_incomplete',
            'missassembled','gapped_missassembled','discordant','gapped_discordant',
            'no_telomeres','gapped_no_telomeres']:
    abs_vals_Cpct.extend(list(df[col].astype(int).values))
annotate_bars_with_values(axC_pct, abs_vals_Cpct)
bottom_legend(axC_pct, [Patch(facecolor=colors[k], label=k.replace('_',' '))
                        for k in dataC_pct.columns if dataC_pct[k].sum() > 0], ncol=2)

# Save MODE 2
for ext, dpi in [('pdf', None), ('svg', None), ('png', 600)]:
    save_kwargs = {'format': ext, 'bbox_inches': 'tight'}
    if dpi is not None:
        save_kwargs['dpi'] = dpi
    pct_fig.savefig(os.path.join(plot_dir, f"telomere_completeness_percent.{ext}"), **save_kwargs)

print("Saved figures to", os.path.abspath(plot_dir))
