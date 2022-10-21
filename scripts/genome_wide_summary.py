# %% Importing libraries
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from plot_config import PBLOCK_COLORS, pblocks

# %% Output paths
output = Path('../C_altered_region_summary_plots')
output.mkdir(exist_ok=True)

# %% Plot 3A: Number of observed pblocks per alternative protein isoforms
fig = plt.figure(figsize=(4, 2.4))
ax = sns.histplot(
    x = pd.cut(
        pblocks.groupby(['anchor', 'other']).size(),
        bins = [1, 2, 3, 4, 5, 14],
        right = False,
        labels = ['1', '2', '3', '4', '5+']
    ),
    color = '#888888',
    edgecolor = 'k',
    alpha = 1,
)
ax.set_xlabel('Number of altered regions\nper isoform')
ax.set_ylabel('Number of alternative\nprotein isoforms')
ax.set_ylim(0, 30000)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

fig.savefig(output/'altered-regions-per-isoform.png', dpi=200, facecolor=None, bbox_inches='tight')

# %% Plot 3B: Distribution of lengths of the insertion, deletion and substituion affected regions for proteins 
aa_loss = pblocks[pblocks['pblock_category'].isin({'DELETION', 'SUBSTITUTION'})].reset_index()[['pblock_category', 'aa_loss']]
aa_loss['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (reference)', inplace=True)
aa_loss.rename(columns={'aa_loss': 'length'}, inplace=True)
aa_gain = pblocks[pblocks['pblock_category'].isin({'INSERTION', 'SUBSTITUTION'})].reset_index()[['pblock_category', 'aa_gain']]
aa_gain['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (alternative)', inplace=True)
aa_gain.rename(columns={'aa_gain': 'length'}, inplace=True)
affected_lengths = pd.concat([aa_loss, aa_gain])

binwidth = 50
xmax = 600

fig = plt.figure(figsize=(4.5, 2))
data = affected_lengths[affected_lengths['pblock_category'] != 'SUBSTITUTION (alternative)']
ax = sns.histplot(
    data = data,
    x = 'length',
    binwidth = binwidth,
    binrange = (0, xmax),
    stat = 'count',
    color = '#808080',
    alpha = 1,
)
ax.set_xlabel('Length of altered region (amino acids)')
ax.set_ylabel('Number of\naltered regions')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
ax.vlines(data['length'].median(), *ax.get_ylim(), color='#b0b0b0', linestyle='-', linewidth=1)

fig.savefig(output/'altered-region-affected-lengths.png', dpi=200, facecolor=None, bbox_inches='tight')

# %%
facets = sns.displot(
    data = affected_lengths,
    x = 'length',
    binwidth = binwidth,
    binrange = (0, xmax),
    stat = 'count',
    row = 'pblock_category',
    hue = 'pblock_category',
    palette = PBLOCK_COLORS,
    row_order = ('DELETION', 'INSERTION', 'SUBSTITUTION (reference)', 'SUBSTITUTION (alternative)'),
    legend = False,
    alpha = 1,
    height = 2,
    aspect = 2.5
)
facets.set_xlabels('Length of altered region (amino acids)')
facets.set_ylabels('Number of\naltered regions')
for category, ax in facets.axes_dict.items():
    ax.set_title(category.capitalize())
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax.vlines(affected_lengths[affected_lengths['pblock_category'] == category]['length'].median(), *ax.get_ylim(), color='#808080', linestyle='-', linewidth=1)

facets.fig.savefig(output/'altered-region-affected-lengths-categories.png', dpi=200, facecolor=None, bbox_inches='tight')

# %% Plot 3C: Substitution scatter plot 
plt.figure(figsize=(4.8, 3.6))
ax = sns.histplot(
    data = pblocks[pblocks['pblock_category'] == 'SUBSTITUTION'],
    x = 'aa_loss',
    y = 'aa_gain',
    binwidth = binwidth/2,
    stat = 'count',
    color = PBLOCK_COLORS['SUBSTITUTION'],
    legend = False,
    cbar = True,
    cbar_kws = {
        'label': 'Number of regions',
    },
    alpha = 1,
)
# ax.spines.right.set_visible(False)
# ax.spines.top.set_visible(False)
ax.set_xlim(0, xmax)
ax.set_ylim(0, xmax)
ax.set_xlabel('Length of substitution region \nin reference isoform (AA)')
ax.set_ylabel('Length of substitution region \nin alternative isoform (AA)')
plt.savefig(output/'substitution-reference-alternative-lengths.png', dpi=200, facecolor=None, bbox_inches='tight')

# %% Pie chart
category_counts = pblocks['pblock_category'].value_counts()
fig = plt.figure()
wedges, texts = plt.pie(
    category_counts,
    colors = category_counts.index.map(PBLOCK_COLORS),
    wedgeprops = {'width': 0.4},
    startangle = 0,
)
for wedge in wedges:
    wedge.set_edgecolor('k')

fig.savefig(output/'altered-region-category-donut.png', dpi=200, facecolor=None, bbox_inches='tight')

# %%
