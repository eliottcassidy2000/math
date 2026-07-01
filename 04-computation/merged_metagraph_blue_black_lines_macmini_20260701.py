#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S83 -- THE BLUE/BLACK LINE STRUCTURE of the MERGED METAGRAPH (tiling explorer).

Ground truth (tournament-tiling-explorer.html):
- TILES: for y=1..n-2, x=n..y+2 -> tile (x,y);  m = C(n-1,2) tiles; base path fixed.
- isGridSym(bits): invariant under the anti-diagonal reflection TRANS_MAP: (x,y) -> (n-y+1, n-x+1).
- LINE = {tiling t, flip(t)} where flip = complement-tiling (flip ALL tile bits). 2^m tilings pair into
  2^{m-1} lines. A line is BLUE iff isGridSym(t) (== isGridSym(flip(t))), BLACK otherwise.
- MERGE = by transpose (classTransposeTarget); SC class = transpose-self; NS classes pair up.
- Per CLAUDE.md: SC(transpose-self) are PURE BLUE or MIXED; NS are PURE BLACK.

We build the merged metagraph, classify each node PURE-BLUE / PURE-BLACK / MIXED, and count per node:
  blue/black lines to OTHER merged nodes, and blue/black SELF-loops (both endpoints same merged node).
Then verify the owner's parity claims and the self-loop-only-on-mixed conjecture.
"""
import itertools
from math import comb
from collections import defaultdict

def build(n):
    VERTS = list(range(n, 0, -1))
    TILES = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            TILES.append((x, y))
    m = len(TILES); assert m == comb(n-1, 2)
    tileIdx = {t: i for i, t in enumerate(TILES)}
    TRANS = [tileIdx[(n-y+1, n-x+1)] for (x, y) in TILES]
    vpos = {v: i for i, v in enumerate(VERTS)}
    perms = list(itertools.permutations(range(n)))
    def bits_to_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1] = 1
        for i, (xL, yL) in enumerate(TILES):
            xi, yi = vpos[xL], vpos[yL]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A
    def canon(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best
    tilings = []
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        gs = all(TRANS[i] == i or bits[i] == bits[TRANS[i]] for i in range(m))
        flipmask = mask ^ ((1 << m) - 1)
        # transpose tiling mask
        tb = [0]*m
        for i in range(m): tb[TRANS[i]] = bits[i]
        transmask = sum(b << k for k, b in enumerate(tb))
        tilings.append(dict(mask=mask, canon=canon(A), gs=gs, flip=flipmask, trans=transmask))
    # class indices
    sigs = sorted(set(t['canon'] for t in tilings))
    cidx = {s: i for i, s in enumerate(sigs)}
    for t in tilings: t['ci'] = cidx[t['canon']]
    mask2t = {t['mask']: t for t in tilings}
    # transpose target per class -> merge
    trans_target = {}
    for t in tilings:
        trans_target.setdefault(t['ci'], mask2t[t['trans']]['ci'])
    # merged node id: SC (self) or min of pair
    def merged(ci): return min(ci, trans_target[ci])
    return n, m, tilings, mask2t, merged, len(sigs), trans_target

def analyze(n):
    n, m, tilings, mask2t, merged, nclasses, trans_target = build(n)
    # node = merged id; collect member classes & tiling counts & gridSym content
    node_tilings = defaultdict(list)
    for t in tilings: node_tilings[merged(t['ci'])].append(t)
    # classify node: pure blue (all gs), pure black (none gs), mixed
    def classify(node):
        gsvals = set(t['gs'] for t in node_tilings[node])
        if gsvals == {True}: return 'PURE_BLUE'
        if gsvals == {False}: return 'PURE_BLACK'
        return 'MIXED'
    # lines: dedupe {mask, flipmask}
    seen = set(); lines = []
    for t in tilings:
        key = min(t['mask'], t['flip'])
        if key in seen: continue
        seen.add(key)
        a = merged(t['ci']); b = merged(mask2t[t['flip']]['ci'])
        lines.append((a, b, 'blue' if t['gs'] else 'black'))
    assert len(lines) == (1 << (m-1))
    # per-node degree counts
    stat = defaultdict(lambda: dict(blue_other=0, black_other=0, blue_self=0, black_self=0, count=0))
    for node, ts in node_tilings.items(): stat[node]['count'] = len(ts)
    for a, b, col in lines:
        if a == b:
            stat[a][f'{col}_self'] += 1
        else:
            stat[a][f'{col}_other'] += 1; stat[b][f'{col}_other'] += 1
    return n, m, node_tilings, classify, stat, nclasses, lines

print("Verifying the owner's structure on the merged metagraph, n=4,5,6:\n")
for n in [4, 5, 6]:
    n, m, node_tilings, classify, stat, nclasses, lines = analyze(n)
    cats = defaultdict(list)
    for node in node_tilings: cats[classify(node)].append(node)
    nlines = len(lines)
    nblue = sum(1 for *_, c in lines if c == 'blue'); nblack = nlines - nblue
    print(f"=== n={n}: m={m}, classes={nclasses}, merged nodes={len(node_tilings)}, "
          f"lines=2^(m-1)={nlines} (blue={nblue}, black={nblack}) ===")
    for cat in ['PURE_BLUE', 'PURE_BLACK', 'MIXED']:
        nodes = cats.get(cat, [])
        print(f"  {cat}: {len(nodes)} nodes")
        for node in sorted(nodes):
            s = stat[node]
            print(f"     node {node:2d}: tilings={s['count']:3d} | blue_other={s['blue_other']} "
                  f"black_other={s['black_other']} | blue_self={s['blue_self']} black_self={s['black_self']}")
    # parity checks
    def par(x): return 'even' if x % 2 == 0 else 'ODD'
    pb_black_other_even = all(stat[v]['black_other'] % 2 == 0 for v in cats.get('PURE_BLACK', []))
    mixed_black_even = all(stat[v]['black_other'] % 2 == 0 for v in cats.get('MIXED', []))
    mixed_blue_odd = all(stat[v]['blue_other'] % 2 == 1 for v in cats.get('MIXED', []))
    pblue_blue_odd = all(stat[v]['blue_other'] % 2 == 1 for v in cats.get('PURE_BLUE', []))
    selfloops_only_mixed = all(stat[v]['blue_self']+stat[v]['black_self'] == 0
                               for cat in ['PURE_BLUE','PURE_BLACK'] for v in cats.get(cat, []))
    print(f"  CLAIMS: pure-black black_other even={pb_black_other_even}; "
          f"mixed black_other even={mixed_black_even}, blue_other odd={mixed_blue_odd}; "
          f"pure-blue blue_other odd={pblue_blue_odd}; self-loops only on mixed={selfloops_only_mixed}\n")
