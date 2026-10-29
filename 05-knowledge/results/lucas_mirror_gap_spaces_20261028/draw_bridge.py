"""Illustration only; check_bridges.py supplies the exact algebraic checks."""
from pathlib import Path
import cmath
from itertools import combinations
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'lucas-mirror-20261028'
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from check_bridges import words, reader

root = Path(__file__).parent
W = words(5)
zeta = cmath.exp(2j*cmath.pi/5)
z = {b: sum((b >> i & 1)*zeta**i for i in range(5)) for b in W}
labels = {b: reader(b, 5, 11, 4) for b in W}
fig, ax = plt.subplots(figsize=(7, 7))
fig.patch.set_facecolor('#fcfaf5')
ax.set_facecolor('#fcfaf5')
for a, b in combinations(W, 2):
    if (a ^ b).bit_count() != 1:
        continue
    if a.bit_count() > b.bit_count():
        a, b = b, a
    ax.annotate('', (z[b].real, z[b].imag), (z[a].real, z[a].imag),
                arrowprops={'arrowstyle': '-|>', 'color': '#a3aaa7',
                            'lw': 1.4, 'shrinkA': 13, 'shrinkB': 13}, zorder=1)
palette = {0: '#66727b', 1: '#b08224', 2: '#1a7c7c'}
for b in W:
    ax.scatter([z[b].real], [z[b].imag], s=600, color=palette[b.bit_count()],
               edgecolor='white', linewidth=1.5, zorder=3)
    ax.text(z[b].real, z[b].imag, str(labels[b]), color='white',
            ha='center', va='center', fontsize=13, weight='bold', zorder=4)
bylabel = {x: b for b, x in labels.items()}
za, zb = z[bylabel[1]], z[bylabel[2]]
ax.annotate('', (zb.real, zb.imag), (za.real, za.imag),
            arrowprops={'arrowstyle': '-|>', 'color': '#ad4143', 'lw': 1.8,
                        'linestyle': '--', 'shrinkA': 14, 'shrinkB': 14,
                        'connectionstyle': 'arc3,rad=-0.19'}, zorder=2)
ax.set_title('Eleven Lucas states in the golden plane', fontsize=17, pad=17)
ax.text(0.5, -0.07, 'Gold: one bit  •  Teal: two bits  •  Gray: empty word\n'
        'Solid arrows: the 15 unit digit moves\n'
        'Dashed 1 → 2: the extra frame selector (not unit length)',
        transform=ax.transAxes, ha='center', va='top', fontsize=11, linespacing=1.6)
ax.set_aspect('equal')
ax.set_xlim(-1.28, 1.28)
ax.set_ylim(-1.14, 1.14)
ax.axis('off')
fig.subplots_adjust(bottom=0.20, top=0.9, left=0.05, right=0.95)
for ext in ('png', 'svg'):
    fig.savefig(root / ('golden_lucas_bridge.'+ext), dpi=180,
                facecolor=fig.get_facecolor(), bbox_inches='tight',
                metadata={'Date': None} if ext == 'svg' else None)
svg = root / 'golden_lucas_bridge.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
plt.close(fig)
