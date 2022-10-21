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
nterm_pblocks = pblocks[~pblocks['nterm'].isna() & (pblocks['nterm'] != 'ALTERNATIVE_ORF') & (pblocks['cterm'] != 'ALTERNATIVE_ORF')].copy()
nterm_pblocks['nterm'].replace(NTERM_CLASSES, inplace=True)
nterm_pblocks['altTSS'] = nterm_pblocks['events'].apply(lambda x: set(x).intersection('BbPp')).astype(bool)

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

# %%
tss_fig = plt.figure(figsize=(5, 2))
ax = sns.countplot(
    data = nterm_pblocks,
    y = 'nterm',
    order = ('Alternative initiation exon', 'Hybrid initiation exon'),
    palette = NTERM_COLORS,
    edgecolor = 'k',
    saturation = 1,
)
sns.countplot(
    ax = ax,
    data = nterm_pblocks[nterm_pblocks['altTSS']],
    y = 'nterm',
    order = ('Alternative initiation exon', 'Hybrid initiation exon'),
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

# %%
nterm_length_fig = plt.figure(figsize=(5, 2))
ax = sns.boxplot(
    data = nterm_pblocks,
    x = 'anchor_relative_length_change',
    y = 'nterm',
    order = ('Alternative initiation exon', 'Hybrid initiation exon'),
    palette = NTERM_COLORS,
    saturation = 1,
    fliersize = 2,
    flierprops = {'marker': 'x'}
)
xmax = max(ax.get_xlim())
ymin, ymax = ax.get_ylim()
ax.vlines(x=0, ymin=ymin, ymax=ymax, color='#808080', linewidth=1, linestyle='--')
ax.set_xlim(-1, 1)
ax.set_xlabel('Change in N-terminal length\n(fraction of reference isoform length)')
ax.set_ylabel(None)

nterm_length_fig.savefig(output/'nterm-rel-length-change.png', dpi=200, facecolor=None, bbox_inches='tight')

# %%
aie_rel_lengths = nterm_pblocks[nterm_pblocks['nterm'] == 'Alternative initiation exon']['anchor_relative_length_change']
hie_rel_lengths = nterm_pblocks[nterm_pblocks['nterm'] == 'Hybrid initiation exon']['anchor_relative_length_change']
mannwhitneyu(aie_rel_lengths, hie_rel_lengths)

# %%
