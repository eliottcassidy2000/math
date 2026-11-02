"""Illustration only; all proof certificates are in verify.py."""
from pathlib import Path
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify import vertices
from fractions import Fraction

HERE = Path(__file__).resolve().parent
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
for ax, b in zip(axs[:2], (2, 3)):
    pts = vertices(Fraction(b), 4)
    for w, (xf, yf) in pts.items():
        if not w:
            continue
        x, y = float(xf), float(yf)
        px, py = map(float, pts[w[:-1]])
        side = 1 if w[-1] == '1' else -1
        centre = px-side*(3*b+5)/(4*b)*py
        radius = math.hypot(px-centre, py)
        xx = np.linspace(min(x, px), max(x, px), 50)
        yy = np.sqrt(np.maximum(0, radius*radius-(xx-centre)**2))
        ax.plot(xx, yy, color='#286891', linewidth=1.2)
    ax.scatter([float(p[0]) for p in pts.values()], [float(p[1]) for p in pts.values()],
               s=15, color='#102f44', zorder=3)
    if b == 3:
        ax.axvspan(1/3, 2/3, color='#edb66b', alpha=.22)
    ax.set(xlim=(-.03, 1.03), ylim=(0, 1.06), xlabel='Horizontal coordinate', ylabel='Height')
    ax.set_title('b=2: adjoining child regions' if b == 2 else 'b=3: a gap retains branch separation')
    ax.spines[['top','right']].set_visible(False)
k = np.arange(1, 13)
for b, color in ((2, '#9e4268'), (3, '#286891')):
    distance = [2*math.asinh(((b-2)*b**(int(t)-1)+1)/2) for t in k]
    axs[2].plot(k, distance, 'o-', label=f'b={b}', color=color, markersize=4)
axs[2].plot(k, 2*k*math.log(3)-2*math.log(3), '--', color='#ca8b37', label='Proved b=3 lower bound')
axs[2].set(xlabel='Common depth k; tree distance is 2k', ylabel='Hyperbolic distance',
           title='Two addresses across the root split')
axs[2].legend(frameon=False, fontsize=9)
axs[2].spines[['top','right']].set_visible(False)
fig.suptitle('The metric stays hyperbolic; a boundary gap repairs the lost tree distance', fontsize=14)
fig.savefig(HERE/'gap_geometry.png', dpi=160)
fig.savefig(HERE/'gap_geometry.svg')
print(HERE/'gap_geometry.png')
