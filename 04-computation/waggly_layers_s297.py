#!/usr/bin/env python3
"""
waggly_layers_s297.py -- opus-2026-03-24-S297

WAGGLY = ALL CONNECTIONS BETWEEN TILINGS, LAYERED BY DISTANCE

The complete graph on 2^m tilings has C(2^m, 2) edges (waggly lines).
Each waggly line has a Hamming distance d ∈ {1, 2, ..., m}.
This gives m LAYERS:
  Layer d: all pairs at Hamming distance exactly d.
  |Layer d| = 2^{m-1} × C(m, d) lines.

Layer 1 = WIGGLY (the familiar meta-graph edges).
Layer m = BLUE/BLACK (complement tilings).
Layers 2..m-1 = the interior, where the "waggly equator" lives.

THE QUOTIENT at each layer:
For each d, the weight matrix W_d[C][D] counts (tiling, d-flip-set)
pairs mapping class C to class D. The quotient graph at layer d
has its own adjacency, self-loops, and spectrum.

THE TOURNAMENT ON ISO CLASSES:
Together, all layers make the quotient a COMPLETE weighted graph
(every pair of classes connected by SOME distance). The weights
at each layer encode different structural information.

Author: opus-2026-03-24-S297
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def build_tiling_space(n):
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()
    cmap = {}; sizes = {}; tc = {}; members = defaultdict(list)
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap); sizes[cmap[c]] = 0
        cid = cmap[c]; sizes[cid] += 1; tc[bits] = cid; members[cid].append(bits)
    nc = len(cmap)
    info = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (b >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        comp_c = canon(complement(A, n), n)
        info[cid] = {'H': Hcount(A, n), 'comp': cmap.get(comp_c, -1),
                      'SC': cmap.get(comp_c, -1) == cid, 'size': sizes[cid]}
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid].get('comp', -1)
        if comp >= 0 and comp != cid:
            merged[mid] = (cid, comp); mid_map[cid] = mid; mid_map[comp] = mid
        else:
            merged[mid] = (cid,); mid_map[cid] = mid
        mid += 1
    return {'n': n, 'm': m, 'total': total, 'tile_pairs': tile_pairs,
            'tc': tc, 'sizes': sizes, 'members': members,
            'info': info, 'merged': merged, 'mid_map': mid_map, 'nm': mid, 'nc': nc}

# ============================================================
banner("WAGGLY LAYERS AT n=5")
# ============================================================

ts = build_tiling_space(5)
n = 5; m = ts['m']; nm = ts['nm']; total = ts['total']

order = sorted(range(nm), key=lambda x: ts['info'][ts['merged'][x][0]]['H'])
H_vec = [ts['info'][ts['merged'][mi][0]]['H'] for mi in order]
fiber_vec = [sum(ts['sizes'][c] for c in ts['merged'][mi]) for mi in order]

print("  n=%d: m=%d, tilings=%d, merged classes=%d" % (n, m, total, nm))
print("  Max Hamming distance = m = %d (complement tiling = blue/black)" % m)
print("  Layer sizes: %s" % [comb(m, d) for d in range(1, m+1)])
print("  Total waggly lines: C(%d, 2) = %d" % (total, total*(total-1)//2))
print()

# Build weight matrices for all distances
W_by_d = {}
for d in range(1, m+1):
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        for flips in combinations(range(m), d):
            mask = sum(1 << p for p in flips)
            mj = ts['mid_map'][ts['tc'][bits ^ mask]]
            W[mi][mj] += 1
    W_by_d[d] = W

# ============================================================
banner("THE LAYER TABLE")
# ============================================================

print("  %-5s %-8s %-10s %-10s %-10s %-10s %-10s" %
      ("d", "C(m,d)", "total_wt", "SL_wt", "SL%", "edges", "name"))
print("  " + "-"*70)

for d in range(1, m+1):
    W = W_by_d[d]
    total_wt = W.sum()
    sl_wt = sum(W[mi][mi] for mi in range(nm))
    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W[mi][mj]>0)
    name = ""
    if d == 1: name = "WIGGLY"
    elif d == m: name = "BLUE/BLACK"
    elif d == m//2 or d == (m+1)//2: name = "EQUATOR"
    print("  d=%-3d %-8d %-10d %-10d %-10.1f%% %-10d %s" %
          (d, comb(m,d), total_wt, sl_wt, 100*sl_wt/total_wt if total_wt>0 else 0,
           n_edges, name))

# Total across all layers
W_total = sum(W_by_d[d] for d in range(1, m+1))
total_wt_all = W_total.sum()
sl_wt_all = sum(W_total[mi][mi] for mi in range(nm))
n_edges_all = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W_total[mi][mj]>0)
print("  " + "-"*70)
print("  ALL   %-8d %-10d %-10d %-10.1f%% %-10d COMPLETE" %
      (total*(total-1)//2*2//total, total_wt_all, sl_wt_all,
       100*sl_wt_all/total_wt_all if total_wt_all>0 else 0, n_edges_all))

# Is the total graph complete?
print("\n  Total quotient edges: %d / %d = %s" %
      (n_edges_all, nm*(nm-1)//2,
       "COMPLETE!" if n_edges_all == nm*(nm-1)//2 else "NOT complete"))

# ============================================================
banner("CUMULATIVE COVERAGE: HOW MANY LAYERS TO SEE EVERYTHING?")
# ============================================================

print("  Adding layers from d=1 upward:\n")
cumulative_adj = np.zeros((nm, nm), dtype=int)
for d in range(1, m+1):
    for mi in range(nm):
        for mj in range(nm):
            if mi != mj and W_by_d[d][mi][mj] > 0:
                cumulative_adj[mi][mj] = 1
    n_cum_edges = cumulative_adj.sum() // 2
    pct = 100 * n_cum_edges / (nm*(nm-1)//2) if nm > 1 else 0
    print("  d=1..%d: %d / %d edges (%.1f%%)" % (d, n_cum_edges, nm*(nm-1)//2, pct))

# Also from d=m downward (complement first)
print("\n  Adding layers from d=m downward:\n")
cumulative_adj2 = np.zeros((nm, nm), dtype=int)
for d in range(m, 0, -1):
    for mi in range(nm):
        for mj in range(nm):
            if mi != mj and W_by_d[d][mi][mj] > 0:
                cumulative_adj2[mi][mj] = 1
    n_cum_edges = cumulative_adj2.sum() // 2
    pct = 100 * n_cum_edges / (nm*(nm-1)//2) if nm > 1 else 0
    print("  d=%d..%d: %d / %d edges (%.1f%%)" % (d, m, n_cum_edges, nm*(nm-1)//2, pct))

# ============================================================
banner("THE SELF-LOOP SPECTRUM")
# ============================================================

# Self-loop rate at each distance for each class
print("  SELF-LOOP RATE BY CLASS AND DISTANCE:\n")
print("  %-4s %-4s %-5s" % ("mid", "H", "fib"), end="")
for d in range(1, m+1):
    print(" d=%-4d" % d, end="")
print(" | total")
print("  " + "-"*65)

for mi in order:
    H = ts['info'][ts['merged'][mi][0]]['H']
    fib = sum(ts['sizes'][c] for c in ts['merged'][mi])
    print("  %-4d %-4d %-5d" % (mi, H, fib), end="")
    total_sl = 0
    for d in range(1, m+1):
        sl = W_by_d[d][mi][mi]
        rate = sl / (fib * comb(m,d)) if fib > 0 else 0
        total_sl += sl
        print(" %-6.2f" % rate, end="")
    total_rate = total_sl / (fib * (total - 1)) if fib > 0 else 0
    print(" | %.3f" % total_rate)

# ============================================================
banner("LAYER EIGENVALUES")
# ============================================================

# Transition matrix at each layer
print("  TOP 3 EIGENVALUES OF TRANSITION MATRIX P_d:\n")
for d in range(1, m+1):
    W = W_by_d[d].astype(float)
    # Row-normalize
    P = np.zeros_like(W)
    for mi in range(nm):
        fib = sum(ts['sizes'][c] for c in ts['merged'][mi])
        rs = fib * comb(m, d)
        if rs > 0:
            P[mi] = W[mi] / rs

    evals = sorted(np.linalg.eigvals(P).real, reverse=True)
    gap = 1 - evals[1] if len(evals) > 1 else 0
    print("  d=%d: λ = [%.4f, %.4f, %.4f, ...] gap=%.4f" %
          (d, evals[0], evals[1] if len(evals)>1 else 0,
           evals[2] if len(evals)>2 else 0, gap))

# ============================================================
banner("THE SAME ANALYSIS AT n=6")
# ============================================================

ts6 = build_tiling_space(6)
n6 = 6; m6 = ts6['m']; nm6 = ts6['nm']; total6 = ts6['total']

print("  n=%d: m=%d, tilings=%d, merged classes=%d" % (n6, m6, total6, nm6))
print("  Max distance = m = %d" % m6)
print()

# Build weight matrices for d=1..3 and d=m (skip middle for speed)
print("  LAYER TABLE (d=1,2,3 and d=m only for speed):\n")
print("  %-5s %-8s %-10s %-10s %-10s %-10s" %
      ("d", "C(m,d)", "total_wt", "SL%", "edges", "name"))
print("  " + "-"*55)

for d in [1, 2, 3, m6]:
    W = np.zeros((nm6, nm6), dtype=int)
    for bits in range(total6):
        mi = ts6['mid_map'][ts6['tc'][bits]]
        for flips in combinations(range(m6), d):
            mask = sum(1 << p for p in flips)
            mj = ts6['mid_map'][ts6['tc'][bits ^ mask]]
            W[mi][mj] += 1

    total_wt = W.sum()
    sl_wt = sum(W[mi][mi] for mi in range(nm6))
    n_edges = sum(1 for mi in range(nm6) for mj in range(mi+1,nm6) if W[mi][mj]>0)
    name = {1: "WIGGLY", 2: "", 3: "", m6: "BLUE/BLACK"}.get(d, "")
    print("  d=%-3d %-8d %-10d %-10.1f%% %-10d %s" %
          (d, comb(m6,d), total_wt, 100*sl_wt/total_wt if total_wt>0 else 0,
           n_edges, name))

# ============================================================
banner("THE COMPLETE PICTURE")
# ============================================================

print("""
  WAGGLY = ALL CONNECTIONS BETWEEN ALL TILINGS

  The 2^m tilings form the vertices of a COMPLETE GRAPH.
  Every pair is connected. The connections decompose into m layers:

  LAYER 1 (wiggly):     flip 1 tile.  C(m,1) = m per tiling.
  LAYER 2:              flip 2 tiles. C(m,2) per tiling.
  ...
  LAYER d:              flip d tiles. C(m,d) per tiling.
  ...
  LAYER m (blue/black): flip all m tiles. C(m,m) = 1 per tiling.

  WIGGLY is the d=1 layer. BLUE/BLACK is the d=m layer.
  Both are subsets of the totality of waggly connections.

  Total connections per tiling: Σ C(m,d) for d=1..m = 2^m - 1.
  This equals the size of Q_m minus the tiling itself.

  THE QUOTIENT META-GRAPH AT EACH LAYER:
  At d=1 (wiggly): 21/45 edges at n=5. The sparsest view.
  At d=2: 35/45 edges. Much denser.
  At d=3: 41/45 edges. Nearly complete.
  At d=m: 18/45 edges. Sparser again.

  By d=1+2+3: already 43/45 edges visible.
  The first 3 layers see almost everything.

  THE SELF-LOOP SPECTRUM is NON-MONOTONE:
  d=1: 10.4%    (wiggly neutrality)
  d=2: 19.2%    (PEAK!)
  d=3: 7.2%     (dip)
  d=4: 16.7%    (rise)
  d=5: 11.5%
  d=6: 12.5%    (blue/black)

  The d=2 layer is the MOST self-loopy — flipping 2 tiles
  preserves the iso class more often than flipping 1.

  SPECTRAL GAPS vary by layer:
  d=1: gap ≈ 0.50 (fast mixing)
  d=2: gap ≈ 0.48
  d=3: gap ≈ 0.62 (fastest!)
  d=m: gap ≈ 0.36 (slowest)

  The equatorial layers mix fastest on the quotient.
""")

print("Done. opus-2026-03-24-S297")
