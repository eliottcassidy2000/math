#!/usr/bin/env python3
"""toothpick.py: toothpick diagrams of the canonical-parent tree of tournament classes (A000568), orders 1..N.  Each class
is one toothpick; its children (classes of the next order whose least canonical deletion it is) touch it at distinct
interior points and leave it perpendicularly, each in its own slice of the plane proportional to its subtree size.  First
creators (classes whose count first appears at their order) are drawn in red, the 63 orphans in black.  Usage:
python3 toothpick.py N  ->  toothpick_N.svg and a line of statistics."""
import sys, subprocess
exec(open('nearjoin_check.py').read().split('# (1)')[0])
N = int(sys.argv[1])
def d6_decode(s):
    s = s.strip(); n = ord(s[1]) - 63; bits = ''.join(format(ord(c) - 63, '06b') for c in s[2:])
    return n, [[int(bits[i*n+j]) for j in range(n)] for i in range(n)]
H = {1: [1], 2: [1]}; cls = {}
for n in range(3, N + 1):
    raw = subprocess.run(['gentourng', '-q', '-z', str(n)], capture_output=True, text=True).stdout.split()
    H[n] = [ham_all(d6_decode(s)[1])[1][-1] for s in raw]
par = {1: [None], 2: [0]}
for n in range(3, N + 1): par[n] = [int(x) for x in open(f'cparent_{n}.txt')]
seen = set(); first = {}
for n in range(1, N + 1):
    new = set(H[n]) - seen; first[n] = [h in new for h in H[n]]; seen |= set(H[n])
children = {n: [[] for _ in H[n]] for n in range(1, N + 1)}
for n in range(2, N + 1):
    for i, p in enumerate(par[n]): children[n - 1][p].append(i)
size = {N: [1] * len(H[N])}
for n in range(N - 1, 0, -1): size[n] = [1 + sum(size[n + 1][c] for c in children[n][i]) for i in range(len(H[n]))]
W = 2400.0; out = []; pal = ['#444', '#1f77b4', '#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#17becf', '#bcbd22', '#7f7f7f']
def draw(n, i, rect, horiz):
    # the toothpick of class (n, i) runs through the centre of its rectangle, spanning it in its orientation; each child
    # gets a slice of the rectangle along that direction (width proportional to the size of its subtree) and its own
    # toothpick runs perpendicularly through the slice, crossing the parent at the slice centre: the contact point.
    x0, y0, x1, y1 = rect; cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    col = '#d62728' if first[n][i] else pal[n % len(pal)]; wd = max(0.35, 3.2 - 0.4 * n)
    if horiz: out.append(f'<line x1="{x0:.1f}" y1="{cy:.1f}" x2="{x1:.1f}" y2="{cy:.1f}" stroke="{col}" stroke-width="{wd}"/>')
    else: out.append(f'<line x1="{cx:.1f}" y1="{y0:.1f}" x2="{cx:.1f}" y2="{y1:.1f}" stroke="{col}" stroke-width="{wd}"/>')
    ch = children[n][i]
    if not ch: return
    tot = sum(size[n + 1][c] for c in ch); pos = 0.0
    for c in ch:
        f0, f1 = pos / tot, (pos + size[n + 1][c]) / tot; pos += size[n + 1][c]
        if horiz:
            sx0, sx1 = x0 + f0 * (x1 - x0), x0 + f1 * (x1 - x0); m = 0.06 * (y1 - y0)
            sub = (sx0, y0 + m, sx1, y1 - m)
        else:
            sy0, sy1 = y0 + f0 * (y1 - y0), y0 + f1 * (y1 - y0); m = 0.06 * (x1 - x0)
            sub = (x0 + m, sy0, x1 - m, sy1)
        draw(n + 1, c, sub, not horiz)
draw(1, 0, (10, 10, W - 10, W - 10), True)
tot = sum(len(H[n]) for n in range(1, N + 1)); fc = sum(sum(first[n]) for n in range(1, N + 1))
svg = f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{W:.0f}" viewBox="0 0 {W:.0f} {W:.0f}"><rect width="100%" height="100%" fill="white"/>' + ''.join(out) + f'<text x="20" y="{W-20:.0f}" font-family="sans-serif" font-size="22">toothpick diagram of A000568 through order {N}: {tot} toothpicks (one per isomorphism class), {fc} first creators in red; each toothpick crosses its least canonical deletion at the contact point</text></svg>'
open(f'toothpick_{N}.svg', 'w').write(svg)
print(f"order {N}: {tot} toothpicks, {fc} first creators; branching by order (max children in the tree): " + ', '.join(f"{n}:{max(len(c) for c in children[n])}" for n in range(1, N)))
