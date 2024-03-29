# %% Importing libraries
from pathlib import Path
from re import M
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import csv
from plot_config import PBLOCK_COLORS, SECTION_COLORS, pblocks

# %% Output paths
output = Path('../C_altered_region_summary_plots')
output.mkdir(exist_ok=True)

# %% Fig2 panel A: Number of altered isoforms per gene vs number of genes
fig = plt.figure(figsize=(4, 2.4))
bins = list(range(1, 11)) + [100]
ax = sns.histplot(
    x = pd.cut(
        pblocks.groupby('anchor')['other'].nunique(),
        bins = bins,
        right = False,
        labels = [str(x) for x in bins[:-2]] + [f'{bins[-2]}+'],
    ),
    shrink = 0.75,
    color = '#888888',
    edgecolor = 'k',
    alpha = 1,
)
ax.set_xlabel('Number of alternative isoforms\nper gene')
ax.set_ylabel('Number of genes')
ax.set_ylim(0, 5000)
### output
fig.savefig(output/'alternative-isoforms-per-gene.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
###
pblocks.groupby('anchor')['other'].nunique().to_frame(name='count').to_csv(output/'alternative-isoforms-per-gene-table.tsv', sep='\t')

# %% Fig2 panel B: Number of observed pblocks per alternative protein isoforms
fig = plt.figure(figsize=(4, 2.4))
ax = sns.histplot(
    x = pd.cut(
        pblocks.groupby(['anchor', 'other']).size(),
        bins = [1, 2, 3, 4, 5, 14],
        right = False,
        labels = ['1', '2', '3', '4', '5+']
    ),
    shrink = 0.75,
    color = '#888888',
    edgecolor = 'k',
    alpha = 1,
)
ax.set_xlabel('Number of altered regions\nper isoform')
ax.set_ylabel('Number of alternative\nprotein isoforms')
# ax.set_ylim(0, 30000)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

### output
fig.savefig(output/'altered-regions-per-isoform.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
pblocks.groupby(['anchor', 'other']).size().to_frame(name='num_alt_regions').to_csv(output/'altered-regions-per-isoform-table.tsv', sep='\t')

# %% Fig2 panel C: Distribution of lengths of the insertion, deletion and substituion affected regions for proteins 
aa_loss = pblocks[pblocks['pblock_category'].isin({'DELETION', 'SUBSTITUTION'})].reset_index()[['anchor','other','pblock_category', 'aa_loss']]
aa_loss['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (reference)', inplace=True)
aa_loss.rename(columns={'aa_loss': 'length'}, inplace=True)
aa_gain = pblocks[pblocks['pblock_category'].isin({'INSERTION', 'SUBSTITUTION'})].reset_index()[['anchor','other','pblock_category', 'aa_gain']]
aa_gain['pblock_category'].replace('SUBSTITUTION', 'SUBSTITUTION (alternative)', inplace=True)
aa_gain.rename(columns={'aa_gain': 'length'}, inplace=True)
affected_lengths = pd.concat([aa_loss, aa_gain])

binwidth = 50
xmax = 600
xtick = 200

fig = plt.figure(figsize=(5, 2))
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

### output
fig.savefig(output/'altered-region-affected-lengths.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
affected_lengths[affected_lengths['pblock_category'] != 'SUBSTITUTION (alternative)'].to_csv(output/'altered-region-affected-lengths-table.tsv', sep='\t')

# %% Fig2 panel D: Distribution of the length of altered protein regions across the annotated proteome
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
    ax.set_xticks(range(0, xmax+1, xtick))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax.vlines(affected_lengths[affected_lengths['pblock_category'] == category]['length'].median(), *ax.get_ylim(), color='#808080', linestyle='-', linewidth=1)

### output
facets.fig.savefig(output/'altered-region-affected-lengths-categories.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
affected_lengths.to_csv(output/'altered-region-affected-lengths-categories-table.tsv', sep='\t')

# %%
#  T-test for the means of two independent samples of scores.
insertion = affected_lengths[affected_lengths['pblock_category'] == 'INSERTION']['length']
deletion = affected_lengths[affected_lengths['pblock_category'] == 'DELETION']['length']
substitution_ref = affected_lengths[affected_lengths['pblock_category'] == 'SUBSTITUTION (reference)']['length']
substitution_alt = affected_lengths[affected_lengths['pblock_category'] == 'SUBSTITUTION (alternative)']['length']

stats.mannwhitneyu(insertion, deletion) #MannwhitneyuResult(statistic=37940405.5, pvalue=0.0)
stats.mannwhitneyu(insertion, substitution_ref) #MannwhitneyuResult(statistic=20566509.0, pvalue=0.0)
stats.mannwhitneyu(insertion, substitution_alt) #MannwhitneyuResult(statistic=51275501.0, pvalue=1.1174708861070221e-07)

# Also alternative package
import pingouin as pg
pg.mwu(insertion, deletion, alternative='two-sided')
pg.mwu(insertion, substitution_ref, alternative='two-sided')

# Mann-Whitney U Test in R using stats library
# wilcox.test(insertion_data$V1, deletion_data$V1, alternative = "two.sided")
# data:  insertion_data$V1 and deletion_data$V1
# W = 37940406, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


# %% Fig2 panel I =: Substitution scatter plot 
plt.figure(figsize=(4.8, 3.6))
ax = sns.histplot(
    data = pblocks[pblocks['pblock_category'] == 'SUBSTITUTION'],
    x = 'aa_gain',
    y = 'aa_loss',
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
ax.set_xlim(0, xmax)
ax.set_ylim(0, xmax)
ax.set_xticks(range(0, xmax+1, xtick))
ax.set_yticks(range(0, xmax+1, xtick))
ax.set_xlabel('Length of substitution region \nin alternative isoform (AA)')
ax.set_ylabel('Length of substitution region \nin reference isoform (AA)')
### output
plt.savefig(output/'substitution-reference-alternative-lengths.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
pblocks.query("pblock_category == 'SUBSTITUTION'")[['anchor', 'other','pblock_category','aa_gain','aa_loss']].to_csv(output/'substitution-reference-alternative-lengths-table.tsv', sep='\t')

# %% Fig2 panel D: Pie chart
category_counts = pblocks['pblock_category'].value_counts()
total_pblocks = category_counts.sum()
fig, ax = plt.subplots()
wedges, texts, autotexts = plt.pie(
    category_counts,
    colors = category_counts.index.map(PBLOCK_COLORS),
    wedgeprops = {'width': 0.4},
    startangle = 180,
    counterclock = False,
    autopct = lambda x: f'{np.round(total_pblocks*x/100):.0f}\n({x:.0f}%)',
    pctdistance = 1.3,
)
for i, wedge in enumerate(wedges):
    wedge.set_edgecolor('k')
### output
fig.savefig(output/'altered-region-category-donut.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
pblocks['pblock_category'].value_counts().to_csv(output/'altered-region-category-donut-table.tsv', sep='\t')

# %%
def get_section(nterm, cterm):
    if nterm and cterm:
        return 'Full-length'
    elif nterm:
        return 'N-terminal'
    elif cterm:
        return 'C-terminal'
    else:
        return 'Internal'

pblocks['protein_section'] = list(map(get_section, ~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()))
pblock_sections = pblocks['protein_section'].value_counts()

fig, ax = plt.subplots(figsize=(6, 1))
left = 0
for section, color in SECTION_COLORS.items():
    val = pblock_sections[section]
    label = f'{val:g}\n({100*val/pblock_sections.sum():0.1f}%)'
    if section == 'Full-length':
        left += 5000
        label_type = 'edge'
        padding = 5
    else:
        label_type = 'center'
        padding = 0
    bar = plt.barh(
        [0],
        val,
        left = left,
        color = color,
        edgecolor = 'k',
        label = section,
    )
    plt.bar_label(bar, labels=[label], label_type=label_type, padding=padding)
    left = left + pblock_sections[section]
ax.legend(loc='upper left', bbox_to_anchor=(0, 0, 1, -0.1), ncols=2, frameon=False)
plt.axis('off')
### output
fig.savefig(output/'protein-section-counts.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
with open(output/'protein-section-counts-table.tsv', 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=SECTION_COLORS.keys(), delimiter='\t')
    writer.writeheader()
    writer.writerow(SECTION_COLORS)

# %% Identifying frameshift cases 
pblocks[pblocks['cblocks'].str.count('FRAME') > 2]
