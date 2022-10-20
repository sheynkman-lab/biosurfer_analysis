# %%
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plot_config import PBLOCK_COLORS

# %%
for font_path in fm.findSystemFonts():
    fm.fontManager.addfont(font_path)

font = {
    'family': 'sans-serif',
    'sans-serif': ['Arial'],
    'weight': 'normal',
    'size': 16
}
mpl.rc('font', **font)

# %%
source = Path('../B_hybrid_aln_gencode_v42')
output = Path('../D_altered_region_summary_plots')
output.mkdir(exist_ok=True)

pblocks = pd.read_csv(source/'pblocks.tsv', sep='\t', index_col=['other', 'pblock_number'])
pblocks

# %%
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

# %%
aa_loss = pblocks[pblocks['pblock_category'].isin({'DELETION', 'SUBSTITUTION'})].reset_index()[['pblock_category', 'aa_loss']]
aa_loss['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (reference)', inplace=True)
aa_loss.rename(columns={'aa_loss': 'length'}, inplace=True)
aa_gain = pblocks[pblocks['pblock_category'].isin({'INSERTION', 'SUBSTITUTION'})].reset_index()[['pblock_category', 'aa_gain']]
aa_gain['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (alternative)', inplace=True)
aa_gain.rename(columns={'aa_gain': 'length'}, inplace=True)
affected_lengths = pd.concat([aa_loss, aa_gain])

xmax = 500

# fig = plt.figure()
facets = sns.displot(
    data = affected_lengths[affected_lengths['length'] <= xmax],
    x = 'length',
    binwidth = 20,
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
# facets.set_titles('{row_name}')
facets.set_xlabels('Length of altered region (amino acids)')
facets.set_ylabels('Number of\naltered regions')
for category, ax in facets.axes_dict.items():
    ax.set_title(category.capitalize())
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax.vlines(affected_lengths[affected_lengths['pblock_category'] == category]['length'].median(), *ax.get_ylim(), color='#808080', linestyle='-', linewidth=1)

facets.fig.savefig(output/'altered-region-affected-lengths.png', dpi=200, facecolor=None, bbox_inches='tight')

# %%
pblock_per_isoform=pblocks.groupby(['anchor', 'other'])
sns.scatterplot(data=pblocks[pblocks["pblock_category"]=="SUBSTITUTION"], x="aa_loss", y='aa_gain')

# %%
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
nterm_pblocks = pblocks[~pblocks['nterm'].isnull()]

