import colorsys

import matplotlib as mpl
import matplotlib.colors as mc
import matplotlib.font_manager as fm
import pandas as pd
from seaborn import color_palette

# Setting configurations for plotting
for font_path in fm.findSystemFonts():
    fm.fontManager.addfont(font_path)

font = {
    'family': 'sans-serif',
    'sans-serif': ['Arial'],
    'weight': 'normal',
    'size': 16
}
mpl.rc('font', **font)

# from https://stackoverflow.com/a/49601444
def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    cnew = colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])
    return mc.to_hex(cnew)

PBLOCK_COLORS = {
    'DELETION': '#ff0082',
    'INSERTION': '#05e0ff',
    'SUBSTITUTION': '#ffd700',
}

PBLOCK_COLORS['SUBSTITUTION (reference)'] = adjust_lightness(PBLOCK_COLORS['SUBSTITUTION'], 1)
PBLOCK_COLORS['SUBSTITUTION (alternative)'] = adjust_lightness(PBLOCK_COLORS['SUBSTITUTION'], 1)

NTERM_CLASSES = {
    'MUTUALLY_EXCLUSIVE': 'Alternative initiation exon',
    'DOWNSTREAM_SHARED': 'Hybrid initiation exon',
    'UPSTREAM_SHARED': 'Splicing-dependent initiation',
    'MUTUALLY_SHARED': 'Alternative translation initiation'
}
NTERM_COLORS = dict(zip(
    NTERM_CLASSES.values(),
    color_palette('viridis_r', n_colors=len(NTERM_CLASSES)+1)[:-1]
))

SPLICE_EVENT_COLORS = {
    'Intron': '#EBA85F',
    'Alt. donor': '#86BB6F',
    'Alt. acceptor': '#A26FBB',
    'Single exon': '#649FD2',
    'Mutually exclusive exons': '#3674AF',
    'Compound': '#888888'
}

# # GENCODE v42
pblocks = pd.read_csv('../B_hybrid_aln_gencode_v42/pblocks.tsv', sep='\t')
# # WTC11
# pblocks = pd.read_csv('../B_hybrid_aln_wtc11_v42/pblocks.tsv', sep='\t')
