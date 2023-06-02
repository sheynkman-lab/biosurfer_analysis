#!/usr/bin/env python
# coding: utf-8

#%%Importing libraries

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

from plot_config import CTERM_CLASSES, CTERM_PALETTE, cterm_splice_palette, cterm_frameshift_palette, pblocks

# %% Output paths
output = Path('../E_cterm_summary_plots')
output.mkdir(exist_ok=True)

#%%
from plot_config import CTERM_CLASSES, CTERM_PALETTE, cterm_splice_palette, cterm_frameshift_palette, pblocks

cterm_pblocks = pblocks[~pblocks['cterm'].isna() & (pblocks['nterm'].isna()) & (pblocks['cterm'] != "ALTERNATIVE_ORF") & (pblocks['cterm'] != "UNKNOWN")].copy()
cterm_pblocks['cterm'] = cterm_pblocks['cterm'].map(CTERM_CLASSES).astype('category')
# Changed string to set for intersection
cterm_pblocks['APA'] = cterm_pblocks['events'].apply(lambda x: set(x).intersection('BbPp')).astype(bool)

#%% Plot 1: Distribution of C-termini categories for alternative isoforms

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
### output
plt.savefig(output/'cterm-class-counts.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
cterm_pblocks.query("cterm in ['Splice-driven', 'Frameshift-driven']")[['anchor','other','cterm']].to_csv(output/'cterm-class-counts-table.tsv', sep='\t')

# %% Plot 2: 
cterm_pblock_events = cterm_pblocks['up_stop_events'].combine(cterm_pblocks['down_stop_events'], lambda x, y: (x, y))
single_ATE = (cterm_pblocks['cterm'] == 'Splice-driven') & cterm_pblocks['tblock_events'].isin({('B', 'b'), ('b', 'B')})
cterm_splice_subcats = pd.DataFrame(
    {
        'Exon extension introduces termination': cterm_pblocks['up_stop_events'].isin({'P', 'I', 'D'}),
        'Alternative terminal exon(s)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}),
        'Poison exon inclusion': cterm_pblocks['up_stop_events'] == 'E',
        'Other': [True for _ in cterm_pblocks.index]
    }
)
cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

cterm_splice_palette_dict = dict(zip(
    cterm_splice_subcats.columns,
    cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + ['#bbbbbb']
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

sns.boxenplot(
    ax = axs[1],
    data = cterm_pblocks[cterm_pblocks['cterm'] == 'Splice-driven'],
    x = cterm_pblocks['anchor_relative_length_change'].abs(),
    y = 'splice_subcat',
    k_depth = 'trustworthy',
    trust_alpha = 0.01,
    order = splice_subcat_order,
    palette = cterm_splice_palette_dict,
    saturation = 1,
    linewidth = 1,
    box_kws = {'edgecolor': 'k'},
    line_kws = {'color': 'k', 'alpha': 1},
)
xmax = max(axs[1].get_xlim())
ymin, ymax = axs[1].get_ylim()
axs[1].set(yticklabels=[])
axs[1].set_xlim([0, 1])
axs[1].set_ylim([ymin, ymax])
axs[1].set_xlabel('Change in C-terminal length\n (fraction of anchor isoform length)')
axs[1].set_ylabel('')
### output
plt.savefig(output/'cterm-splicing-subcats.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
cterm_pblocks.assign(anchor_relative_length_change = cterm_pblocks['anchor_relative_length_change'].abs())[['anchor','other', 'splice_subcat','anchor_relative_length_change']].to_csv(output/'cterm-splicing-subcats-table.tsv', sep='\t')

#%% Plot 3: 
cterm_frame_subcats = pd.DataFrame(
    {
        'Single exon': cterm_pblocks['up_stop_cblock_events'].isin({'E', 'e'}),
        'Alt. acceptor': cterm_pblocks['up_stop_cblock_events'].isin({'A', 'a'}),
        'Multi-exon skipping': cterm_pblocks['up_stop_cblock_events'].apply(lambda x: set(x) if x else set()) == {'e'},
        'Alt. donor': cterm_pblocks['up_stop_cblock_events'].isin({'D', 'd'}),
        'Intron': cterm_pblocks['up_stop_cblock_events'].isin({'I', 'i'}),
        'Other': [True for _ in cterm_pblocks.index]
    }
)
cterm_pblocks['frame_subcat'] = cterm_frame_subcats.idxmax(axis=1)#.astype(pd.CategoricalDtype(cterm_frame_subcats.columns, ordered=True))

cterm_frameshift_palette_dict = dict(zip(
    cterm_frame_subcats.columns,
    cterm_frameshift_palette + ['#bbbbbb']
))

frame_subcat_order = cterm_pblocks[cterm_pblocks['cterm'] == 'Frameshift-driven']['frame_subcat'].value_counts().index

cterm_frameshift_fig, axs = plt.subplots(1, 2, figsize=(9, 4))
sns.countplot(
    ax = axs[0],
    data = cterm_pblocks[cterm_pblocks['cterm'] == 'Frameshift-driven'],
    y = 'frame_subcat',
    order = frame_subcat_order,
    palette = cterm_frameshift_palette_dict,
    saturation = 1,
    edgecolor = 'k',
    linewidth = 1,
)
axs[0].set_xlabel('Number of alternative isoforms')
axs[0].set_ylabel(None)

sns.boxenplot(
    ax = axs[1],
    data = cterm_pblocks[cterm_pblocks['cterm'] == 'Frameshift-driven'],
    x = cterm_pblocks['anchor_relative_length_change'].abs(),
    y = 'frame_subcat',
    k_depth = 'trustworthy',
    trust_alpha = 0.01,
    order = frame_subcat_order,
    palette = cterm_frameshift_palette_dict,
    saturation = 1,
    linewidth = 1,
    box_kws = {'edgecolor': 'k'},
    line_kws = {'color': 'k', 'alpha': 1},
)
xmax = max(axs[1].get_xlim())
ymin, ymax = axs[1].get_ylim()
axs[1].set(yticklabels=[])
axs[1].set_xlim([0, 1])
axs[1].set_ylim([ymin, ymax])
axs[1].set_xlabel('Change in C-terminal length\n (fraction of reference isoform length)')
axs[1].set_ylabel('')
### output
plt.savefig(output/'cterm-frameshift-subcats.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
cterm_pblocks.assign(anchor_relative_length_change = cterm_pblocks['anchor_relative_length_change'].abs()).query("cterm in 'Frameshift-driven'")[['anchor','other', 'frame_subcat','anchor_relative_length_change']].to_csv(output/'cterm-frameshift-subcats-table.tsv', sep='\t')

# %%
