# %%
from pathlib import Path
from matplotlib.patches import Patch
from scipy.stats.contingency import chi2_contingency
from itertools import combinations

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from plot_config import PBLOCK_COLORS, SPLICE_EVENT_COLORS, pblocks

# %% Output paths
output = Path('../F_internal_summary_plots')
output.mkdir(exist_ok=True)

# %%
internal_pblocks = (
    pblocks[pblocks['internal']].
    drop(columns=[col for col in pblocks.columns if 'start' in col or 'stop' in col]).
    copy()
)
# convert string repr back to Python object
internal_pblocks['tblock_events'] = internal_pblocks['tblock_events'].map(eval)
internal_pblocks['events'] = internal_pblocks['events'].map(eval)  

internal_subcats = pd.DataFrame(
    {
        'Snapback': internal_pblocks['frameshift'],
        'Intron': internal_pblocks['tblock_events'].isin({('I',), ('i',)}),
        'Alt. donor': internal_pblocks['tblock_events'].isin({('D',), ('d',)}),
        'Alt. acceptor': internal_pblocks['tblock_events'].isin({('A',), ('a',)}),
        'Single exon': internal_pblocks['tblock_events'].isin({('E',), ('e',)}),
        'Compound': [True for _ in internal_pblocks.index]
    }
)
subcat_order = ('Intron', 'Alt. donor', 'Alt. acceptor', 'Single exon', 'Compound', 'Snapback')
internal_pblocks['splice_event'] = internal_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(subcat_order, ordered=True))

# %%
internal_pblocks_fig = plt.figure(figsize=(4.6, 3.8))
ax = sns.countplot(
    data = internal_pblocks.sort_values('pblock_category', ascending=True),
    y = 'splice_event',
    dodge = True,
    hue = 'pblock_category',
    palette = PBLOCK_COLORS,
    saturation = 1,
    edgecolor = 'k',
)
plt.legend(loc='upper right', labels=['Deletions', 'Insertions', 'Substitutions'])
ax.set_xlabel('Number of altered internal regions')
ax.set_ylabel(None)
### output
internal_pblocks_fig.savefig(output/'internal-events.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
### output
internal_pblocks[['anchor','other','splice_event','pblock_category']].to_csv(output/'internal-events-table.tsv', sep='\t')

# %%
internal_pblocks_ragged_fig = plt.figure(figsize=(4.6, 3.8))
ax = sns.countplot(
    data = internal_pblocks.sort_values('pblock_category', ascending=True),
    y = 'splice_event',
    palette = SPLICE_EVENT_COLORS,
    saturation = 1,
    edgecolor = 'k',
)
sns.countplot(
    ax = ax,
    data = internal_pblocks[internal_pblocks['split_codons']].sort_values('pblock_category', ascending=True),
    y = 'splice_event',
    fill = False,
    edgecolor = 'k',
    hatch = '///',
)
plt.legend(loc='upper right', handles=[Patch(facecolor='w', edgecolor='k', hatch='///')], labels=['Contains \nragged codons'])
ax.set_xlabel('Number of altered internal regions')
ax.set_ylabel(None)
internal_pblocks_ragged_fig.savefig(output/'internal-events-ragged.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
internal_pblocks[['splice_event','split_codons']].to_csv(output/'internal-events-ragged-table.tsv', sep='\t')
# %%
alpha = 0.01
ragged_contingency = pd.crosstab(internal_pblocks['split_codons'], internal_pblocks['splice_event'])
chi2, p_all, dof, expected = chi2_contingency(ragged_contingency)

ps = dict()
for event1, event2 in combinations(internal_subcats.columns, 2):
    sub_contingency = ragged_contingency[[event1, event2]]
    _, ps[event1, event2], _, _ = chi2_contingency(sub_contingency)

ps_sig = {k: p for k, p in ps.items() if p < alpha/len(ps)}
ps_insig = {k: p for k, p in ps.items() if k not in ps_sig}

# %%
internal_rel_length_change = plt.figure(figsize=(4.6, 3.8))
ax = sns.boxenplot(
    data = internal_pblocks,
    x = internal_pblocks['anchor_relative_length_change'].abs(),
    y = 'splice_event',
    k_depth = 'trustworthy',
    trust_alpha = 0.01,
    palette = SPLICE_EVENT_COLORS,
    saturation = 1,
    linewidth = 1,
    box_kws = {'edgecolor': 'k'},
    line_kws = {'color': 'k', 'alpha': 1},
)
ax.set_xlim(0, 1)
ax.set_xlabel('Change in protein length\n(fraction of reference isoform length)')
ax.set_ylabel(None)
internal_rel_length_change.savefig(output/'internal-rel-length-change.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
df = internal_pblocks[['splice_event','anchor_relative_length_change']]
df['anchor_relative_length_change'] = df['anchor_relative_length_change'].abs()
df.to_csv(output/'internal-rel-length-change-table.tsv', sep='\t')

# %%
nagnag_pblocks = internal_pblocks[(internal_pblocks['splice_event'] == 'Alt. acceptor') & (internal_pblocks['length_change'].abs() == 1)]

# %%
internal_compound_pblocks = internal_pblocks[internal_pblocks['splice_event'] == 'Compound'].copy()

internal_compound_subcats = pd.DataFrame(
    {
        'Multi-exon skipping': internal_compound_pblocks['events'] == frozenset('e'),
        'Exon skipping + \nalt. donor/acceptor': internal_compound_pblocks['events'].isin({
            frozenset(sorted('de')),
            frozenset(sorted('De')),
            frozenset(sorted('ea')),
            frozenset(sorted('eA')),
            frozenset(sorted('dea')),
            frozenset(sorted('Dea')),
            frozenset(sorted('deA')),
            frozenset(sorted('DeA')),
        }),
        'Mutually exclusive exons': internal_compound_pblocks['tblock_events'].isin({('E', 'e'), ('e', 'E')}),
        'Multi-exon inclusion': internal_compound_pblocks['events'] == frozenset('E'),
        'Alt. donor + alt. acceptor': internal_compound_pblocks['events'].isin({
            frozenset(sorted('ad')),
            frozenset(sorted('Ad')),
            frozenset(sorted('aD')),
            frozenset(sorted('AD')),
        }),
        'Exon inclusion + \nalt. donor/acceptor': internal_compound_pblocks['events'].isin({
            frozenset(sorted('dE')),
            frozenset(sorted('DE')),
            frozenset(sorted('Ea')),
            frozenset(sorted('EA')),
            frozenset(sorted('dEa')),
            frozenset(sorted('DEa')),
            frozenset(sorted('dEA')),
            frozenset(sorted('DEA')),
        }),
        'Other': [True for _ in internal_compound_pblocks.index]
    }
)
internal_compound_pblocks['compound_subcat'] = internal_compound_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(internal_compound_subcats.columns, ordered=True))

internal_pblocks_compound_fig = plt.figure(figsize=(3, 3))
ax = sns.countplot(
        data = internal_compound_pblocks,
        y = 'compound_subcat',
        palette = 'Greys_r',
        saturation = 1,
        edgecolor = 'k',
)
ax.set_xlabel('Number of altered\ninternal regions'),
ax.set_ylabel(None)
internal_pblocks_compound_fig.savefig(output/'internal-compound-events.png', dpi=200, facecolor=None, bbox_inches='tight')

#Output source data
internal_compound_pblocks[['anchor','other','compound_subcat']].to_csv(output/'internal-compound-events-table.tsv', sep='\t')
# %%
