#!/usr/bin/env python
# coding: utf-8

#%%Importing libraries

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from plot_config import NTERM_CLASSES, NTERM_COLORS, CTERM_CLASSES, CTERM_PALETTE, cterm_splice_palette, cterm_frameshift_palette, pblocks

# %% Output paths
output = Path('../E_cterm_summary_plots')
output.mkdir(exist_ok=True)

#%%

cterm_pblocks = pblocks[~pblocks['cterm'].isna() & (pblocks['nterm'] != "ALTERNATIVE_ORF") & (pblocks['cterm'] != "ALTERNATIVE_ORF") & (pblocks['cterm'] != "UNKNOWN")].copy()
cterm_pblocks['cterm'] = cterm_pblocks['cterm'].astype('category')
# Changed string to set for intersection
cterm_pblocks['APA'] = cterm_pblocks['events'].apply(lambda x: set(x).intersection('BbPp')).astype(bool)

display(pd.crosstab(cterm_pblocks['up_stop_cblock_category'], cterm_pblocks['down_stop_cblock_category'], margins=True))


#%% Plot 1: Distribution of C-termini categories for alternative isoforms

cterm_fig = plt.figure(figsize=(3,2))
ax = sns.countplot(
    data = cterm_pblocks,
    y = 'cterm',
    order = CTERM_CLASSES.keys(),
    palette = CTERM_PALETTE,
    saturation = 1,
    linewidth = 1,
    edgecolor = 'k'
)
ax.set_xlabel('Number of alternative isoforms')
ax.set_ylabel('')
ax.set_yticklabels(['',''])
plt.savefig(output/'E_cterm_plots/cterm-class-counts.svg', dpi=200, facecolor=None)


# %% Plot 2: 

cterm_pblock_events = cterm_pblocks['up_stop_events'].combine(cterm_pblocks['down_stop_events'], lambda x, y: (x, y))
single_ATE = (cterm_pblocks['cterm'] == "SPLICING") & cterm_pblocks['tblock_events'].isin({('B', 'b'), ('b', 'B')})
cterm_splice_subcats = pd.DataFrame(
    {
        'exon extension introduces termination': cterm_pblocks['up_stop_events'].isin({'P', 'I', 'D'}),
        'alternative terminal exon(s)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}),
        'poison exon': cterm_pblocks['up_stop_events'] == 'E',
        'other': [True for _ in cterm_pblocks.index]
    }
)
cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

cterm_splice_palette_dict = dict(zip(
    cterm_splice_subcats.columns,
    cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + ['#bbbbbb']
))
splice_subcat_order = tuple(cterm_splice_subcats.keys())

cterm_splice_fig, axs = plt.subplots(1, 2, figsize=(10, 5))

sns.countplot(
    ax = axs[0],
    data = cterm_pblocks[cterm_pblocks['cterm'] == "SPLICING"],
    y = 'splice_subcat',
    order = splice_subcat_order,
    palette = cterm_splice_palette_dict,
    saturation = 1,
    edgecolor = 'k',
)
axs[0].set_xlabel('Number of alternative isoforms')
axs[0].set_ylabel(None)
axs[0].set_yticklabels(('Exon extension introduces termination','Alternative terminal exon(s)','Poison exon', 'Other'))

sns.boxplot(
    ax = axs[1],
    data = cterm_pblocks[cterm_pblocks['cterm'] == "SPLICING"],
    x = 'anchor_relative_length_change',
    y = 'splice_subcat',
    order = splice_subcat_order,
    palette = cterm_splice_palette_dict,
    saturation = 1,
    linewidth = 1,
)
xmax = max(axs[1].get_xlim())
ymin, ymax = axs[1].get_ylim()
axs[1].vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linewidth=1, linestyle=':')
axs[1].set(yticklabels=[])
axs[1].set_xlim([-1,1])
axs[1].set_ylim([ymin, ymax])
axs[1].set_xlabel('Change in C-terminal length\n (fraction of anchor isoform length)')
axs[1].set_ylabel('')
cterm_splice_fig.subplots_adjust( wspace=0.5)
plt.savefig(output/'E_cterm_plots/cterm-splicing-subcats.svg', dpi=200, facecolor=None, bbox_inches='tight')


#%% Plot 3: 
cterm_frame_subcats = pd.DataFrame(
    {
        'exon': cterm_pblocks['up_stop_cblock_events'].isin({'E', 'e'}),
        'acceptor': cterm_pblocks['up_stop_cblock_events'].isin({'A', 'a'}),
        'donor': cterm_pblocks['up_stop_cblock_events'].isin({'D', 'd'}),
        'intron': cterm_pblocks['up_stop_cblock_events'].isin({'I', 'i'}),
        'other': [True for _ in cterm_pblocks.index]
    }
)
cterm_pblocks['frame_subcat'] = cterm_frame_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_frame_subcats.columns, ordered=True))

cterm_frameshift_palette_dict = dict(zip(
    cterm_frame_subcats.columns,
    cterm_frameshift_palette + ['#bbbbbb']
))

frame_subcat_order = cterm_pblocks[cterm_pblocks['cterm'] == "FRAMESHIFT"]['frame_subcat'].value_counts().index

cterm_frameshift_fig, axs = plt.subplots(1, 2, figsize=(10, 5))
sns.countplot(
    ax = axs[0],
    data = cterm_pblocks[cterm_pblocks['cterm'] == "FRAMESHIFT"],
    y = 'frame_subcat',
    order = frame_subcat_order,
    palette = cterm_frameshift_palette_dict,
    saturation = 1,
    linewidth = 1,
    edgecolor='k',
)
axs[0].set_xlabel('Number of alternative isoforms')
axs[0].set_ylabel(None)
axs[0].set_yticklabels(('Exon','Acceptor','Donor', 'Intron', 'Other'))

sns.boxplot(
    ax = axs[1],
    data = cterm_pblocks[cterm_pblocks['cterm'] == "FRAMESHIFT"],
    x = 'anchor_relative_length_change',
    y = 'frame_subcat',
    order = frame_subcat_order,
    palette = cterm_frameshift_palette_dict,
    saturation = 1,
    linewidth = 1,
)
xmax = max(axs[1].get_xlim())
ymin, ymax = axs[1].get_ylim()
axs[1].vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linewidth=1, linestyle=':')
axs[1].set(yticklabels=[])
axs[1].set_xlim([-1,1])
axs[1].set_ylim([ymin, ymax])
axs[1].set_xlabel('Change in C-terminal length\n (fraction of anchor isoform length)')
axs[1].set_ylabel('')
cterm_frameshift_fig.subplots_adjust(wspace=0.5)
plt.savefig(output/'E_cterm_plots/cterm-frameshift-subcats.svg', dpi=200, facecolor=None, bbox_inches='tight')




# %%
