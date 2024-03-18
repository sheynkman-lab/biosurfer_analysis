# %%
from pathlib import Path
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib as mpl

from plot_config import NTERM_CLASSES, NTERM_COLORS, pblocks

# %% Output paths
output = Path('../D_nterm_summary_plots')
output.mkdir(exist_ok=True)

# %%
nterm_pblocks = pblocks[~pblocks['nterm'].isna() & (pblocks['nterm'] != 'ALTERNATIVE_ORF') & (pblocks['cterm'].isna())].copy()
nterm_pblocks['nterm'].replace(NTERM_CLASSES, inplace=True)
nterm_pblocks['altTSS'] = nterm_pblocks['events'].apply(lambda x: eval(x).intersection('BbPp')).astype(bool)

# %% Fig3 panel A (both Alt TSS and 5' UTR AS)
tss_fig = plt.figure(figsize=(5, 4))
ax = sns.countplot(
    data = nterm_pblocks,
    y = 'nterm',
    order = NTERM_COLORS.keys(),
    palette = NTERM_COLORS,
    edgecolor = 'k',
    saturation = 1,
)
sns.countplot(
    ax = ax,
    data = nterm_pblocks[nterm_pblocks['altTSS']],
    y = 'nterm',
    order = NTERM_COLORS.keys(),
    palette = NTERM_COLORS,
    edgecolor = 'k',
    fill = False,
    hatch = '//',
)
ax.legend(
    loc = (0, 1),
    frameon = False,
    handles = [Patch(facecolor='w', edgecolor='k', hatch='///'), Patch(facecolor='w', edgecolor='k')],
    labels = ['Alternative transcription start site', '5\' UTR alternative splicing'],
)
ax.set_xlabel('Number of alternative isoforms')
ax.set_ylabel(None)
plt.savefig(output/'nterm-counts-all_mechanism.png', dpi=500, facecolor=None, bbox_inches='tight')

#Output source data
nterm_pblocks.query("nterm in ['Mutually exclusive starts (MSX)', 'Shared downstream start (SDS)']")[['anchor','other','nterm','altTSS']].to_csv(output/'nterm-counts-all_mechanism.tsv', sep='\t')

# %% Fig3 panel C: MXS vs SDS scatterplot 
font = {
    'family': 'sans-serif',
    'sans-serif': ['Arial'],
    'weight': 'normal',
    'size': 10
}
mpl.rc('font', **font)

# Filter the dataframe for 'Mutually exclusive starts (MXS)' and 'Shared downstream start (SDS)'
msx_data = nterm_pblocks[nterm_pblocks['nterm'] == 'Mutually exclusive starts (MSX)']
sds_data = nterm_pblocks[nterm_pblocks['nterm'] == 'Shared downstream start (SDS)']
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(5.5, 5.5))
msx_color = (0.565498, 0.84243, 0.262877)
sds_color = (0.20803, 0.718701, 0.472873)

sns.scatterplot(data=msx_data, x='aa_loss', y='aa_gain', marker='.', ax=axes[0], alpha=0.2,
                color=msx_color)
axes[0].set_title('Mutually exclusive starts (MXS)', fontsize=11)  
axes[0].set_xlabel('Reference \n(amino acids)', fontsize=10) 
axes[0].set_ylabel('Alternative \n(amino acids)', fontsize=10) 
axes[0].set_xlim(0, 2000)
axes[0].set_ylim(0, 2000)
axes[0].set_aspect('equal') 
axes[0].grid(True, linestyle='--', linewidth=0.5)

sns.scatterplot(data=sds_data, x='aa_loss', y='aa_gain', marker='.', ax=axes[1], alpha=0.2, color=sds_color)
axes[1].set_title('Shared downstream start (SDS)', fontsize=11)  
axes[1].set_xlabel('Reference \n(amino acids)', fontsize=10) 
axes[1].set_ylabel('Alternative \n(amino acids)', fontsize=10) 
axes[1].set_xlim(0, 2000)
axes[1].set_ylim(0, 2000)
axes[1].set_aspect('equal')  
axes[1].grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()

# Save plot
plt.savefig(output/'nterm-rel-length-change_scatterplot.png', dpi=600, facecolor=None, bbox_inches='tight')
plt.show()

#Output source data
nterm_pblocks.query("nterm in ['Mutually exclusive starts (MSX)', 'Shared downstream start (SDS)']")[['anchor','other','aa_loss','aa_gain']].to_csv(output/'nterm_mechanism_affected_len.tsv', sep='\t')



