import matplotlib.colors as mc
import colorsys

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

NTERM_COLORS = {
    
}