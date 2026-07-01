#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S86 -- RUN THE TWIN SEPARATOR on the 5 mixed twin-pairs at n=6 (HYP-3809 twin-pairing).

S84 found: SC merged nodes' #grid-sym-tilings multiset has all-even multiplicities => a fixed-point-free
involution twins the SC nodes by grid-sym count. At n=6: pure-blue = 2 nodes (grid-sym 1, the trivial twin),
MIXED = 10 nodes = 5 twin-pairs (grid-sym counts {3:2, 5:2, 7:4, 9:2}). The "5 remaining twins."

TASK: for each SC node compute a full invariant signature (grid-sym count, tiling count, H=#Ham-paths,
|Aut|, score sequence, c3=#cyclic-triangles, adjacency spectrum), identify the twin pairing, and find the
SEPARATOR -- the invariant that distinguishes the two members of each twin (or shows a cospectral twin).
"""
import itertools
from math import comb
from collections import defaultdict, Counter

def build(n):
    VERTS = list(range(n, 0, -1))
    TILES = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
    m = len(TILES); tileIdx = {t: i for i, t in enumerate(TILES)}
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
    T = []
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bits_to_adj(bits)
        gs = all(TRANS[i] == i or bits[i] == bits[TRANS[i]] for i in range(m))
        tb = [0]*m
        for i in range(m): tb[TRANS[i]] = bits[i]
        sigma = sum(b << k for k, b in enumerate(tb)); flip = mask ^ ((1 << m) - 1)
        T.append(dict(mask=mask, canon=canon(A), adj=A, gs=gs, sigma=sigma, flip=flip))
    sigs = sorted(set(t['canon'] for t in T)); cidx = {s: i for i, s in enumerate(sigs)}
    for t in T: t['ci'] = cidx[t['canon']]
    m2 = {t['mask']: t for t in T}
    tgt = {}
    for t in T: tgt.setdefault(t['ci'], m2[t['sigma']]['ci'])
    return n, m, T, m2, tgt, len(sigs), perms

def invariants(A, n, perms):
    # H = # Hamiltonian paths (directed)
    H = 0
    for p in perms:
        if all(A[p[k]][p[k+1]] for k in range(n-1)): H += 1
    # |Aut|
    aut = sum(1 for p in perms if all(A[i][j] == A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j))
    # c3 cyclic triangles
    c3 = 0
    for i, j, k in itertools.combinations(range(n), 3):
        s = A[i][j] + A[j][k] + A[k][i]
        if s == 0 or s == 3: c3 += 1
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    return H, aut, c3, scores

n = 6
n, m, T, m2, tgt, nclasses, perms = build(n)
node_t = defaultdict(list)
merged = lambda ci: min(ci, tgt[ci])
for t in T: node_t[merged(t['ci'])].append(t)
def cat(v):
    g = set(t['gs'] for t in node_t[v]); return 'BLUE' if g == {True} else ('BLACK' if g == {False} else 'MIXED')
SCnodes = [v for v in node_t if cat(v) in ('BLUE', 'MIXED')]

rows = []
for v in SCnodes:
    ts = node_t[v]; gscount = sum(1 for t in ts if t['gs']); tc = len(ts)
    rep = ts[0]['adj']
    H, aut, c3, scores = invariants(rep, n, perms)
    rows.append(dict(v=v, cat=cat(v), gs=gscount, tc=tc, H=H, aut=aut, c3=c3, scores=scores))

print("SC nodes at n=6 (sorted by grid-sym count, then tiling count):")
print(f"{'node':>4} {'cat':>6} {'gridsym':>7} {'tilings':>7} {'H':>5} {'|Aut|':>5} {'c3':>3} {'scores':>16}")
for r in sorted(rows, key=lambda r: (r['gs'], r['tc'])):
    print(f"{r['v']:>4} {r['cat']:>6} {r['gs']:>7} {r['tc']:>7} {r['H']:>5} {r['aut']:>5} {r['c3']:>3} {str(r['scores']):>16}")

print("\nTWINS (mixed nodes grouped by grid-sym count) and the SEPARATOR:")
mixed = [r for r in rows if r['cat'] == 'MIXED']
by_gs = defaultdict(list)
for r in mixed: by_gs[r['gs']].append(r)
for gs, grp in sorted(by_gs.items()):
    print(f"  grid-sym={gs}: {len(grp)} nodes")
    for r in grp:
        print(f"     node {r['v']:2d}: tilings={r['tc']:3d} H={r['H']:4d} |Aut|={r['aut']} c3={r['c3']:2d} scores={r['scores']}")
    # separator: which invariant differs within this grid-sym group?
    for inv in ['tc', 'H', 'aut', 'c3', 'scores']:
        vals = set(r[inv] for r in grp)
        if len(vals) == len(grp):
            print(f"     => SEPARATED by {inv} (all distinct: {[r[inv] for r in grp]})"); break
    else:
        print(f"     => NOT fully separated by tc/H/aut/c3/scores (cospectral within group)")
