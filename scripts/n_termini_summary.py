# %%
from pathlib import Path
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from plot_config import NTERM_CLASSES, NTERM_COLORS, pblocks

# %% Output paths
output = Path('../D_nterm_summary_plots')
output.mkdir(exist_ok=True)

# %%
nterm_pblocks = pblocks[~pblocks['nterm'].isna() & (pblocks['nterm'] != 'ALTERNATIVE_ORF') & (pblocks['cterm'].isna())].copy()
nterm_pblocks['nterm'].replace(NTERM_CLASSES, inplace=True)
nterm_pblocks['altTSS'] = nterm_pblocks['events'].apply(lambda x: eval(x).intersection('BbPp')).astype(bool)

# %%
fig = plt.figure(figsize=(5, 4))
ax = sns.countplot(
    data = nterm_pblocks,
    y = 'nterm',
    order = NTERM_COLORS.keys(),
    palette = NTERM_COLORS,
    edgecolor = 'k',
    saturation = 1,
)
ax.set_xlabel('Number of altered N-terminal regions')
ax.set_ylabel(None)

fig.savefig(output/'nterm-mechanism-counts.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
nterm_pblocks[['anchor','other','nterm']].to_csv(output/'nterm-mechanism-counts-table.tsv', sep='\t')

# %%
tss_fig = plt.figure(figsize=(5, 2))
ax = sns.countplot(
    data = nterm_pblocks,
    y = 'nterm',
    order = ('Mutually exclusive starts', 'Shared downstream start'),
    palette = NTERM_COLORS,
    edgecolor = 'k',
    saturation = 1,
)
sns.countplot(
    ax = ax,
    data = nterm_pblocks[nterm_pblocks['altTSS']],
    y = 'nterm',
    order = ('Mutually exclusive starts', 'Shared downstream start'),
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
plt.savefig(output/'nterm-altTSS-counts.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
nterm_pblocks.query("nterm in ['Mutually exclusive starts', 'Shared downstream start']")[['anchor','other','nterm','altTSS']].to_csv(output/'nterm-altTSS-counts-table.tsv', sep='\t')

# %%
nterm_length_fig = plt.figure(figsize=(5, 2))
ax = sns.boxenplot(
    data = nterm_pblocks,
    x = nterm_pblocks['anchor_relative_length_change'].abs(),
    y = 'nterm',
    order = ('Mutually exclusive starts', 'Shared downstream start'),
    k_depth = 'trustworthy',
    trust_alpha = 0.01,
    palette = NTERM_COLORS,
    saturation = 1,
    linewidth = 1,
    box_kws = {'edgecolor': 'k'},
    line_kws = {'color': 'k', 'alpha': 1},
)
xmax = max(ax.get_xlim())
ymin, ymax = ax.get_ylim()
ax.vlines(x=0, ymin=ymin, ymax=ymax, color='#808080', linewidth=1, linestyle='--')
ax.set_xlim(0, 1)
ax.set_xlabel('Change in N-terminal length\n(fraction of reference isoform length)')
ax.set_ylabel(None)

nterm_length_fig.savefig(output/'nterm-rel-length-change.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
nterm_pblocks.query("nterm in ['Mutually exclusive starts', 'Shared downstream start']")[['anchor','other','nterm','altTSS']].to_csv(output/'nterm-rel-length-change-table.tsv', sep='\t')

# %%
mxs_rel_lengths = nterm_pblocks[nterm_pblocks['nterm'] == 'Mutually exclusive starts']['anchor_relative_length_change'].abs()
sds_rel_lengths = nterm_pblocks[nterm_pblocks['nterm'] == 'Shared downstream start']['anchor_relative_length_change'].abs()
mannwhitneyu(mxs_rel_lengths, sds_rel_lengths)

# %%
