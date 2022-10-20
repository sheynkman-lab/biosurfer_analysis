# %%
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from plot_config import NTERM_CLASSES, NTERM_COLORS, pblocks

# %% Output paths
output = Path('../D_nterm_summary_plots')
output.mkdir(exist_ok=True)

# %%
nterm_pblocks = pblocks[~pblocks['nterm'].isnull()].copy()
nterm_pblocks['nterm'].replace(NTERM_CLASSES, inplace=True)
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
nterm_length_fig = plt.figure(figsize=(6, 4))
ax = sns.violinplot(
    data = nterm_pblocks,
    x = 'anchor_relative_length_change',
    y = 'nterm',
    order = ('Alternative initiation exon', 'Hybrid initiation exon'),
    gridsize = 200,
    palette = NTERM_COLORS,
    saturation = 1,
    scale = 'area',
)
xmax = max(ax.get_xlim())
ymin, ymax = ax.get_ylim()
ax.vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linewidth=1, linestyle=':')
ax.set_xlim(-1, 1)
ax.set_ylim(ymin, ymax)
ax.set_xlabel('Change in N-terminal length\n(fraction of reference isoform length)')
ax.set_ylabel(None)

# %%
