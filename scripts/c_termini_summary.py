#!/usr/bin/env python

#%%Importing libraries
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import scipy.stats as stats
import matplotlib as mpl

from plot_config import CTERM_CLASSES, CTERM_PALETTE, cterm_splice_palette, cterm_frameshift_palette, pblocks

# %% Output paths
output = Path('../E_cterm_summary_plots')
output.mkdir(exist_ok=True)

#%%
cterm_pblocks = pblocks[~pblocks['cterm'].isna() & (pblocks['nterm'].isna()) & (pblocks['cterm'] != "ALTERNATIVE_ORF") & (pblocks['cterm'] != "UNKNOWN")].copy()
cterm_pblocks['cterm'] = cterm_pblocks['cterm'].map(CTERM_CLASSES).astype('category')
# Changed string to set for intersection
cterm_pblocks['APA'] = cterm_pblocks['events'].apply(lambda x: set(x).intersection('BbPp')).astype(bool)

#%% Fig5 panel A: Frequency of splice-driven and frameshift-driven C-terminal events

cterm_fig = plt.figure(figsize=(3.8, 2))
ax = sns.countplot(
    data = cterm_pblocks,
    y = 'cterm',
    order = CTERM_CLASSES.values(),
    palette = CTERM_PALETTE,
    saturation = 1,
    linewidth = 1,
    edgecolor = 'k',
)
ax.set_xlabel('Number of alternative isoforms')
ax.set_ylabel('')
plt.savefig(output/'cterm-class-counts.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
cterm_pblocks.query("cterm in ['Splice-driven', 'Frameshift-driven']")[['anchor','other','cterm']].to_csv(output/'cterm-class-counts-table.tsv', sep='\t')

# %% Fig5 panel B: Frequency of splice-driven patterns
cterm_pblock_events = cterm_pblocks['up_stop_events'].combine(cterm_pblocks['down_stop_events'], lambda x, y: (x, y))
single_ATE = (cterm_pblocks['cterm'] == 'Splice-driven') & cterm_pblocks['tblock_events'].isin({('B', 'b'), ('b', 'B')})
cterm_splice_subcats = pd.DataFrame(
    {
        'Exon extension introduces termination': cterm_pblocks['up_stop_events'].isin({'P', 'I', 'D'}),
        'Alternative terminal exon(s)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}),
        'Poison exon inclusion': cterm_pblocks['up_stop_events'] == 'E',
        'Other': [True for _ in cterm_pblocks.index]
        #'Alternative last exon in UTR': cterm_pblocks['cblocks'].apply(lambda x: 'TRANSLATED' in x and 'DELETION' in x and 'UNTRANSLATED' not in x)
    }
)
cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

cterm_splice_palette_dict = dict(zip(
    cterm_splice_subcats.columns,
    cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + ['#bbbbbb']
))
splice_subcat_order = tuple(cterm_splice_subcats.keys())

cterm_pblock_events = cterm_pblocks['up_stop_events'].combine(cterm_pblocks['down_stop_events'], lambda x, y: (x, y))
single_ATE = (cterm_pblocks['cterm'] == 'Splice-driven') & cterm_pblocks['tblock_events'].isin({('B', 'b'), ('b', 'B')})


cterm_splice_subcats = pd.DataFrame(
    {
        'Exon extension introduces \n termination (EXIT)': cterm_pblocks['up_stop_events'].isin({'P', 'I', 'D'}),
        'Alternative terminal \n exon(s) (ATE)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}),
        'Alternative last exon \n in UTR (ALE in UTR)': cterm_pblocks.apply(lambda row: 'TRANSLATED' in row['cblocks'] and 'DELETION' in row['cblocks'] and 'UNTRANSLATED' not in row['cblocks'] if row['cterm'] == 'Splice-driven' and row['splice_subcat'] == 'Other' else False, axis=1),
        'Poison exon inclusion': cterm_pblocks['up_stop_events'] == 'E',
        'Cut-out splice terminal \n exon (COSTE)': cterm_pblocks.apply(lambda row: 'DELETION' in row['cblocks'] and 'INSERTION' in row['cblocks'] and 'TRANSLATED' not in row['cblocks'] and 'UNTRANSLATED' not in row['cblocks'] and 'FRAME' not in row['cblocks'] and 'p' in row['tblock_events'] and row['tblock_events'].count('B') == 1 if row['cterm'] == 'Splice-driven' and row['splice_subcat'] == 'Other' else False, axis=1),
        'Other': [True for _ in cterm_pblocks.index]
    }
)
cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

cterm_splice_palette_dict = dict(zip(
    cterm_splice_subcats.columns,
    cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + cterm_splice_palette[3:4] + cterm_splice_palette[4:5] +  ['#bbbbbb']
))
splice_subcat_order = tuple(cterm_splice_subcats.keys())

cterm_splice_fig, axs = plt.subplots(1, 2, figsize=(9, 4))
sns.countplot(
    ax = axs[0],
    data = cterm_pblocks[cterm_pblocks['cterm'] == 'Splice-driven'],
    y = 'splice_subcat',
    order = splice_subcat_order,
    palette = cterm_splice_palette_dict,
    saturation = 1,
    edgecolor = 'k',
    linewidth = 1,
)
axs[0].set_xlabel('Number of alternative isoforms')
axs[0].set_ylabel(None)

plt.savefig(output/'cterm-splicing-subcats.png', dpi=700, facecolor=None, bbox_inches='tight')

#Output source data
cterm_pblocks.assign(anchor_relative_length_change = cterm_pblocks['anchor_relative_length_change'].abs())[['anchor','other', 'splice_subcat','anchor_relative_length_change']].to_csv(output/'cterm-splicing-subcats-table.tsv', sep='\t')

#%% 
#TODO: Mann-Whitney U Test signed ranked test here between SE, alt Acc .. vs Intron
cterm_frameshift=cterm_pblocks[cterm_pblocks['cterm'] == 'Frameshift-driven']
cterm_intron = cterm_pblocks[cterm_pblocks['frame_subcat'] == 'Intron']
cterm_se = cterm_pblocks[cterm_pblocks['frame_subcat'] == 'Single exon']
cterm_altacc = cterm_pblocks[cterm_pblocks['frame_subcat'] == 'Alt. acceptor']
cterm_altdonor = cterm_pblocks[cterm_pblocks['frame_subcat'] == 'Alt. donor']

data = [[cterm_intron],[cterm_frameshift],[cterm_se],[cterm_altacc],[cterm_altdonor]]
data = [[4890, 3499, 3301, 551], [6819, 2247, 1185, 1096, 457, 437]]
stat, p, dof, expected = chi2_contingency(data)


# %% Alternative Last Exon in 3' UTR case from Splice-driven 'Other' category.

cterm_pblock_splice = cterm_pblocks[cterm_pblocks['cterm'] == 'Splice-driven']
cterm_splice_other = cterm_pblock_splice[cterm_pblock_splice['splice_subcat']=='Other']
condition1 = cterm_splice_other['cblocks'].apply(lambda x: 'DELETION' in x and 'TRANSLATED' in x)
condition2 = cterm_splice_other['cblocks'].apply(lambda x: 'UNTRANSLATED' not in x)
cterm_aleutr = cterm_splice_other[condition1 & condition2].copy()
cterm_aleutr.to_csv(output/'cterm-splice-driven-ALEinUTR.tsv', sep='\t')

# %% Cut-out splice terminal exon case from Splice-driven 'Other' category.

condition3 = cterm_splice_other['cblocks'].apply(lambda x: 'DELETION' in x and 'INSERTION' in x)
condition4 = cterm_splice_other['cblocks'].apply(lambda x: 'TRANSLATED' not in x and 'UNTRANSLATED' not in x and 'FRAME' not in x)
condition5 = cterm_splice_other['tblock_events'].apply(lambda x: x.count('B') == 1 and 'p' in x)
cterm_other_new = cterm_splice_other[condition3 & condition4 & condition5].copy()
cterm_other_new.to_csv(output/'cterm-splice-driven-other-NEW.tsv', sep='\t')

# %%  Fig5 panel C & D: 2D scatter plot v2 splice-driven vs frameshift-driven

font = {
    'family': 'sans-serif',
    'sans-serif': ['Arial'],
    'weight': 'normal',
    'size': 10
}
mpl.rc('font', **font)
msx_data = cterm_pblocks[cterm_pblocks['cterm'] == 'Splice-driven']
sds_data = cterm_pblocks[cterm_pblocks['cterm'] == 'Frameshift-driven']
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6, 6))
msx_color = (0.5048212226066897, 0.00392156862745098, 0.47021914648212226)
sds_color = (0.7885121107266436, 0.03238754325259515, 0.13656286043829297)

sns.scatterplot(data=msx_data, x='aa_loss', y='aa_gain', marker='o', ax=axes[0], alpha=0.2,
                color=msx_color)
axes[0].set_title('Splice-driven', fontsize=13)  
axes[0].set_xlabel('Reference \n(amino acids)', fontsize=12)  
axes[0].set_ylabel('Alternative \n(amino acids)', fontsize=12)  
axes[0].set_xlim(0, 3000)
axes[0].set_ylim(0, 3000)
axes[0].set_aspect('equal')  
axes[0].grid(True, linestyle='--', linewidth=0.5)

sns.scatterplot(data=sds_data, x='aa_loss', y='aa_gain', marker='o', ax=axes[1], alpha=0.2, color=sds_color)
axes[1].set_title('Frameshift-driven', fontsize=13)  
axes[1].set_xlabel('Reference \n(amino acids)', fontsize=12)  
axes[1].set_ylabel('Alternative \n(amino acids)', fontsize=12) 
axes[1].set_xlim(0, 3000)
axes[1].set_ylim(0, 3000)
axes[1].set_aspect('equal')  
axes[1].grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.savefig(output/'cterm-rel-length-change_scatterplot.png', dpi=800, facecolor=None, bbox_inches='tight')
plt.show()

#Output source data
cterm_pblocks.query("cterm in ['Splice-driven', 'Frameshift-driven']")[['anchor','other','aa_loss','aa_gain']].to_csv(output/'cterm_mechanism_affected_len.tsv', sep='\t')

# %% Supplementary Figure S5: 2D scatter plot v2 frameshift-driven subcats

d1 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Exon extension introduces \n termination (EXIT)']
d2 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Alternative terminal \n exon(s) (ATE)']
d3 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Alternative last exon \n in UTR (ALE in UTR)']
d4 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Poison exon inclusion']
d5 = cterm_pblocks[cterm_pblocks['splice_subcat'] == 'Cut-out splice terminal \n exon (COSTE)']

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(6, 15))
colors = [(0.5048212226066897, 0.00392156862745098, 0.47021914648212226),
          (0.735840061514802, 0.061960784313725495, 0.5225682429834679), (0.9094502114571319, 0.2894886582083814, 0.6086120722798923), (0.9754555940023067, 0.5330257593233372, 0.6768935024990388), (0.9859592464436755, 0.7293041138023837, 0.7404229142637447)]

sns.scatterplot(data=d1, x='aa_loss', y='aa_gain', marker='o', ax=axes[0], alpha=0.2,
                color=colors[0])
axes[0].set_title('Exon extension introduces termination', fontsize=30, pad=20)
axes[0].set_xlabel('Reference \n(amino acids)', fontsize=25)
axes[0].set_ylabel('Alternative \n(amino acids)', fontsize=25)
axes[0].set_xlim(0, 3000)
axes[0].set_ylim(0, 3000)
axes[0].grid(True, linestyle='--', linewidth=0.5)

sns.scatterplot(data=d2, x='aa_loss', y='aa_gain', marker='o', ax=axes[1], alpha=0.2,
                color=colors[1])
axes[1].set_title('Alternative terminal exon(s)', fontsize=30, pad=20)
axes[1].set_xlabel('Reference \n(amino acids)', fontsize=25)
axes[1].set_ylabel('Alternative \n(amino acids)', fontsize=25)
axes[1].set_xlim(0, 3000)
axes[1].set_ylim(0, 3000)
axes[1].grid(True, linestyle='--', linewidth=0.5)

sns.scatterplot(data=d3, x='aa_loss', y='aa_gain', marker='o', ax=axes[2], alpha=0.2,
                color=colors[2])
axes[2].set_title('Alternative last exon in UTR', fontsize=30, pad=20)
axes[2].set_xlabel('Reference \n(amino acids)', fontsize=25)
axes[2].set_ylabel('Alternative \n(amino acids)', fontsize=25)
axes[2].set_xlim(0, 3000)
axes[2].set_ylim(0, 3000)
axes[2].grid(True, linestyle='--', linewidth=0.5)

sns.scatterplot(data=d4, x='aa_loss', y='aa_gain', marker='o', ax=axes[3], alpha=0.2,
                color=colors[3])
axes[3].set_title('Poison exon inclusion', fontsize=30, pad=20)
axes[3].set_xlabel('Reference \n(amino acids)', fontsize=25)
axes[3].set_ylabel('Alternative \n(amino acids)', fontsize=25)
axes[3].set_xlim(0, 3000)
axes[3].set_ylim(0, 3000)
axes[3].grid(True, linestyle='--', linewidth=0.5)

sns.scatterplot(data=d5, x='aa_loss', y='aa_gain', marker='o', ax=axes[4], alpha=0.2,
                color=colors[4])
axes[4].set_title('Cut-out splice terminal exon', fontsize=30, pad=20)
axes[4].set_xlabel('Reference \n(amino acids)', fontsize=25)
axes[4].set_ylabel('Alternative \n(amino acids)', fontsize=25)
axes[4].set_xlim(0, 3000)
axes[4].set_ylim(0, 3000)
axes[4].grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()
plt.savefig(output/'cterm-rel-splice-driven-subcat-length-change_scatterplot.png', dpi=900, facecolor=None, bbox_inches='tight')
plt.show()

# Output source data
#cterm_pblocks.query("splice_subcat in ['Exon extension introduces \n termination (EXIT)', 'Alternative terminal \n exon(s) (ATE)', 'Alternative last exon \n in UTR (ALE in UTR)', 'Poison exon inclusion', 'Cut-out splice terminal \n exon (COSTE)']")[['anchor','other','aa_loss','aa_gain']].to_csv(output/'cterm_mechanism_affected_len.tsv', sep='\t')


