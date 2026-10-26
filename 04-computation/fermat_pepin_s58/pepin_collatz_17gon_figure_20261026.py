#!/usr/bin/env python3
"""Figure (SVG, stdlib only): Pepin's test as the valuation-two Collatz map on the regular 17-gon, and the squaring
tree of F_17^* (opus-2026-10-26-S58).
Left: vertices = Z/17; A(x) = (3x+1)/4 = 5x + 13 (mod 17) fixes vertex 1 (the trivial Collatz cycle) and runs one
cycle through the other 16 vertices (Pepin); the integer Fermat line 17 -> 13 -> 10 is its first arc; A^8 is the
mirror x -> 2 - x, i.e. the reflection of the polygon through vertex 1.
Right: x -> x^2 on F_17^* is one binary tree of depth 4 rooted at 1 (Somer-Krizek 2004, Thm 4.1/4.4); 3 is a deepest
leaf and its squaring path 3 -> 9 -> 13 -> 16 -> 1 passes through -1 = 16 (Pepin's congruence).
Palette: reference slots 1-3 (blue, orange, aqua; documented all-pairs-validated), neutral ink for text.
"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', '05-knowledge', 'results', 's58_fermat_pepin_collatz_20261026',
                   'pepin_collatz_17gon.svg')
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#c9c8c2'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
p = 17
W, H = 1180, 660
els = []


def txt(x, y, s, size=13, color=INK, anchor='middle', weight='normal'):
    els.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" text-anchor="{anchor}" '
               f'font-family="Helvetica, Arial, sans-serif" font-weight="{weight}" dominant-baseline="middle">{s}</text>')


def line(x0, y0, x1, y1, color, w=1.2, dash=None, marker=None, shrink=0.0):
    if shrink:
        L = math.hypot(x1 - x0, y1 - y0)
        ux, uy = (x1 - x0) / L, (y1 - y0) / L
        x0, y0, x1, y1 = x0 + ux * shrink, y0 + uy * shrink, x1 - ux * shrink, y1 - uy * shrink
    d = f' stroke-dasharray="{dash}"' if dash else ''
    m = f' marker-end="url(#{marker})"' if marker else ''
    els.append(f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" stroke="{color}" stroke-width="{w}"'
               f' stroke-linecap="round"{d}{m}/>')


def circ(x, y, r, fill, stroke, w=1.6):
    els.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{w}"/>')


# ---- left panel: the 17-gon
cx, cy, R = 300, 378, 195


def pos(k, r=R):
    th = math.pi / 2 - 2 * math.pi * k / p
    return cx + r * math.cos(th), cy - r * math.sin(th)


A = lambda x: (5 * x + 13) % p
assert A(1) == 1
cyc = [0]
while A(cyc[-1]) != 0:
    cyc.append(A(cyc[-1]))
assert len(cyc) == 16 and 1 not in cyc
assert all(pow(5, 8, p) == p - 1 for _ in [0])           # A^8 = x -> 2 - x
for k in range(p):
    line(*pos(k), *pos((k + 1) % p), GRID, 1.0)
u = pos(1, R)
v = pos(1, -R * 1.08)
line(*u, *v, AQUA, 2.0, dash='7 5')
for i in range(16):
    line(*pos(cyc[i]), *pos(cyc[(i + 1) % 16]), INK2, 1.0, marker='arrG', shrink=13)
for a_, b_ in ((0, 13), (13, 10)):
    line(*pos(a_), *pos(b_), ORANGE, 2.6, marker='arrO', shrink=13)
for k in range(p):
    x, y = pos(k)
    if k == 1:
        circ(x, y, 11, BLUE, BLUE)
    else:
        circ(x, y, 8, SURF, INK2, 1.4)
    lx, ly = pos(k, R + 27)
    txt(lx, ly, str(k), 14)
txt(300, 34, 'Pépin on the 17-gon: A(x) = (3x+1)/4 = 5x + 13 (mod 17)', 17, weight='bold')
txt(300, 636, 'one 16-cycle through all vertices but 1;  A⁸ = the mirror through vertex 1', 12.5, INK2)
# legend
lgx, lgy = 22, 66
circ(lgx + 8, lgy, 7, BLUE, BLUE)
txt(lgx + 24, lgy, 'vertex 1 = the trivial Collatz cycle (fixed by A)', 12, anchor='start')
line(lgx, lgy + 24, lgx + 18, lgy + 24, ORANGE, 2.6)
txt(lgx + 24, lgy + 24, 'integer Fermat line 17 → 13 → 10', 12, anchor='start')
line(lgx + 330, lgy, lgx + 348, lgy, INK2, 1.0)
txt(lgx + 354, lgy, 'the Pépin cycle of A (16 steps)', 12, anchor='start')
line(lgx + 330, lgy + 24, lgx + 348, lgy + 24, AQUA, 2.0, dash='7 5')
txt(lgx + 354, lgy + 24, 'mirror axis: A⁸(x) = 2 − x', 12, anchor='start')

# ---- right panel: the squaring tree
levels = {}
for x in range(1, p):
    y, d = x, 0
    while y != 1:
        y = y * y % p
        d += 1
    levels[x] = d
order = {0: [1]}
for d in range(1, 5):
    order[d] = []
    for par in order[d - 1]:
        order[d] += sorted(x for x in range(1, p) if levels[x] == d and x * x % p == par)
left, right, top = 690, 1150, 150
coords = {}
for d in range(5):
    row = order[d]
    for i, x in enumerate(row):
        coords[x] = (left + (i + 0.5) / len(row) * (right - left), top + d * 105)
path = [3, 9, 13, 16, 1]
for x in range(2, p):
    (x0, y0), (x1, y1) = coords[x], coords[x * x % p]
    line(x0, y0, x1, y1, GRID, 1.3)
for a_, b_ in zip(path, path[1:]):
    (x0, y0), (x1, y1) = coords[a_], coords[b_]
    line(x0, y0, x1, y1, ORANGE, 2.8)
for x, (px, py) in coords.items():
    if x == 1:
        circ(px, py, 15, BLUE, BLUE)
        txt(px, py, '1', 13, SURF, weight='bold')
    else:
        circ(px, py, 15, SURF, ORANGE if x in path else INK2, 1.8)
        txt(px, py, str(x), 13)
for d in range(5):
    txt(left - 26, top + d * 105, f'depth {d}', 12, INK2, anchor='end')
txt(905, 34, 'Collatz for squaring: x → x² on F₁₇* is one binary tree', 17, weight='bold')
txt(905, 612, 'leaves = the 8 non-residues = the primitive roots; Pépin: 3 → 9 → 13 → −1 → 1', 12.5, INK2)
txt(905, 634, 'tree depths at F₀..F₄ are 1, 2, 4, 8, 16 (Somer–Křížek 2004, Thm 4.4)', 12.5, INK2)

defs = (f'<defs><marker id="arrG" markerWidth="9" markerHeight="9" refX="7" refY="4.5" orient="auto">'
        f'<path d="M0,0 L9,4.5 L0,9 z" fill="{INK2}"/></marker>'
        f'<marker id="arrO" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{ORANGE}"/></marker></defs>')
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">'
       f'<title>Pépin on the 17-gon and the squaring tree of F17*</title>'
       f'<rect width="{W}" height="{H}" fill="{SURF}"/>{defs}' + ''.join(els) + '</svg>')
open(OUT, 'w', encoding='utf-8').write(svg)
print('wrote', os.path.normpath(OUT), '; Pepin cycle:', cyc)
