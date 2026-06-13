#!/usr/bin/env python3
"""
wiggly_waggly_definitive_s296.py -- opus-2026-03-24-S296

DEFINITIVE WIGGLY AND WAGGLY LINE COUNTS AND DEFINITIONS

DEFINITIONS (to be added to CLAUDE.md):

  WIGGLY LINE: connects two tilings at Hamming distance exactly 1.
    = one edge of Q_m. Flip exactly 1 tile.
    Count per tiling: m
    Total wiggly lines: 2^m × m / 2 = 2^{m-1} × m
    (Each line counted twice, once from each endpoint.)

  WAGGLY LINE: connects two tilings at Hamming distance ≥ 2.
    = non-edge of Q_m. Flip 2 or more tiles.
    Count per tiling: 2^m - 1 - m
    Total waggly lines: 2^m × (2^m - 1 - m) / 2

  BLUE/BLACK LINE: the special waggly line at Hamming distance m.
    = flip ALL m tiles. The complement tiling.
    Count per tiling: 1
    Total blue/black lines: 2^m / 2 = 2^{m-1}
    Blue/black ⊂ waggly (since m ≥ 2 for n ≥ 4).

  MEANINGFUL WAGGLY SUBSETS by Hamming distance d:
    Waggly-d: tilings at distance exactly d (for d = 2, 3, ..., m).
    Count per tiling at distance d: C(m, d)
    Total waggly-d lines: 2^m × C(m, d) / 2 = 2^{m-1} × C(m, d)

Author: opus-2026-03-24-S296
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
        if c not in cmap:
            cmap[c] = len(cmap); sizes[cmap[c]] = 0
        cid = cmap[c]; sizes[cid] += 1; tc[bits] = cid
        members[cid].append(bits)
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
            'info': info, 'merged': merged, 'mid_map': mid_map,
            'nm': mid, 'nc': nc}

# ============================================================
banner("LINE COUNTS: VERIFYING THE FORMULAS")
# ============================================================

for n in [3, 4, 5, 6]:
    m = comb(n-1, 2)
    total = 2**m
    wiggly_per = m
    waggly_per = total - 1 - m
    total_pairs = total * (total - 1) // 2
    total_wiggly = total * m // 2
    total_waggly = total_pairs - total_wiggly

    print("  n=%d: m=%d, tilings=2^%d=%d" % (n, m, m, total))
    print("    Per tiling: %d wiggly + %d waggly = %d total neighbors" %
          (wiggly_per, waggly_per, wiggly_per + waggly_per))
    print("    Total lines: %d wiggly + %d waggly = %d total pairs" %
          (total_wiggly, total_waggly, total_pairs))
    print("    Waggly breakdown by distance:")
    for d in range(2, m+1):
        count_per = comb(m, d)
        total_d = total * count_per // 2
        label = ""
        if d == m: label = " ← BLUE/BLACK"
        elif d == 2: label = " ← waggly-2"
        print("      d=%d: %d per tiling, %d total%s" % (d, count_per, total_d, label))
    # Verify: sum of all distances = total pairs
    check = sum(total * comb(m, d) // 2 for d in range(1, m+1))
    print("    Check: Σ = %d = %d ✓" % (check, total_pairs))
    print()

# ============================================================
banner("QUOTIENT WEIGHT MATRICES FOR ALL DISTANCES (n=5)")
# ============================================================

ts = build_tiling_space(5)
n = 5; m = ts['m']; nm = ts['nm']; total = ts['total']

order = sorted(range(nm), key=lambda x: ts['info'][ts['merged'][x][0]]['H'])
H_vec = [ts['info'][ts['merged'][mi][0]]['H'] for mi in order]
fiber_vec = [sum(ts['sizes'][c] for c in ts['merged'][mi]) for mi in order]

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

# Also build TOTAL waggly (d≥2) by summing
W_waggly_total = sum(W_by_d[d] for d in range(2, m+1))

print("  n=5: QUOTIENT WEIGHT MATRICES\n")
for d_label, W in [("WIGGLY (d=1)", W_by_d[1]),
                    ("ALL WAGGLY (d≥2)", W_waggly_total),
                    ("BLUE/BLACK (d=%d)" % m, W_by_d[m])]:
    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1, nm) if W[mi][mj] > 0)
    n_self = sum(1 for mi in range(nm) if W[mi][mi] > 0)
    total_w = W.sum()
    sl_w = sum(W[mi][mi] for mi in range(nm))
    print("  %s:" % d_label)
    print("    Edges in quotient: %d, Self-loops: %d/%d" % (n_edges, n_self, nm))
    print("    Total weight: %d, Self-loop weight: %d (%.1f%%)" %
          (total_w, sl_w, 100*sl_w/total_w if total_w > 0 else 0))
    print("    H:", H_vec)
    for mi in order:
        row = " ".join("%5d" % W[mi][mj] for mj in order)
        print("    H=%2d: %s  |fib=%d" % (H_vec[order.index(mi)], row, fiber_vec[order.index(mi)]))
    print()

# ============================================================
banner("WIGGLY vs WAGGLY: THE ADJACENCY COMPARISON (n=5)")
# ============================================================

adj_wig = np.array([[1 if mi!=mj and W_by_d[1][mi][mj]>0 else 0 for mj in range(nm)] for mi in range(nm)])
adj_wag = np.array([[1 if mi!=mj and W_waggly_total[mi][mj]>0 else 0 for mj in range(nm)] for mi in range(nm)])
adj_bb = np.array([[1 if mi!=mj and W_by_d[m][mi][mj]>0 else 0 for mj in range(nm)] for mi in range(nm)])

e_wig = adj_wig.sum()//2
e_wag = adj_wag.sum()//2
e_bb = adj_bb.sum()//2

print("  ADJACENCY:")
print("    Wiggly:     %d / %d edges" % (e_wig, nm*(nm-1)//2))
print("    Waggly:     %d / %d edges" % (e_wag, nm*(nm-1)//2))
print("    Blue/black: %d / %d edges" % (e_bb, nm*(nm-1)//2))
print()

# What does wiggly see that waggly doesn't? And vice versa?
wig_only = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_wig[mi][mj] and not adj_wag[mi][mj])
wag_only = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_wag[mi][mj] and not adj_wig[mi][mj])
both = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_wig[mi][mj] and adj_wag[mi][mj])
neither = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if not adj_wig[mi][mj] and not adj_wag[mi][mj])

print("  OVERLAP:")
print("    Both wiggly AND waggly: %d" % both)
print("    Only wiggly:            %d" % wig_only)
print("    Only waggly:            %d" % wag_only)
print("    Neither:                %d (disconnected pairs)" % neither)
print()

if wig_only > 0:
    print("  WIGGLY-ONLY edges (connected by d=1 but NOT by d≥2):")
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if adj_wig[mi][mj] and not adj_wag[mi][mj]:
                Hi = ts['info'][ts['merged'][mi][0]]['H']
                Hj = ts['info'][ts['merged'][mj][0]]['H']
                print("    {%d(H=%d) ↔ %d(H=%d)}" % (mi, Hi, mj, Hj))

if wag_only > 0:
    print("  WAGGLY-ONLY edges (connected by d≥2 but NOT by d=1):")
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if adj_wag[mi][mj] and not adj_wig[mi][mj]:
                Hi = ts['info'][ts['merged'][mi][0]]['H']
                Hj = ts['info'][ts['merged'][mj][0]]['H']
                # Which distances connect them?
                dists = [d for d in range(2, m+1) if W_by_d[d][mi][mj] > 0]
                print("    {%d(H=%d) ↔ %d(H=%d)}: at distances %s" % (mi, Hi, mj, Hj, dists))

if neither > 0:
    print("  DISCONNECTED pairs (not in wiggly OR waggly):")
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if not adj_wig[mi][mj] and not adj_wag[mi][mj]:
                Hi = ts['info'][ts['merged'][mi][0]]['H']
                Hj = ts['info'][ts['merged'][mj][0]]['H']
                print("    {%d(H=%d) ↔ %d(H=%d)}" % (mi, Hi, mj, Hj))

# ============================================================
banner("MEANINGFUL WAGGLY SUBSETS")
# ============================================================

print("""
  CREATIVE WAGGLY SUBSETS — beyond just Hamming distance:

  1. STRIDE SUBSETS: flip tiles in the same ROW of the staircase.
     Row r has (n-2-r) tiles. Flipping all tiles in one row is a
     "row mutation" — it reverses all non-base arcs involving one vertex.
     This is meaningful because it corresponds to REVERSING one vertex's
     role in the tournament (making all its tile arcs go the other way).

  2. DIAGONAL SUBSETS: flip tiles on the same anti-diagonal.
     The anti-diagonal of the staircase connects tile (x,y) with x+y = const.
     This is the "wave" direction — it corresponds to Walsh degree.

  3. RANGE SUBSETS: flip all tiles of the same range r = x - y.
     Range-r tiles connect vertices at distance r in the base path.
     Flipping all range-r tiles is a "harmonic" of the staircase.

  4. VERTEX-STAR SUBSETS: flip all tiles incident to vertex v.
     This reverses all of v's non-base arcs. It's a local perturbation
     centered on one vertex.

  5. TRIANGLE SUBSETS: flip tiles forming a triangle in the staircase.
     Three tiles (a,b), (b,c), (a,c) form a 3-cycle. Flipping all three
     is a "3-cycle rotation" — it cyclically permutes the relative order
     of three vertices while leaving others unchanged.
""")

# Compute the RANGE subsets
tile_pairs = ts['tile_pairs']
range_groups = defaultdict(list)
for k, (lo, hi) in enumerate(tile_pairs):
    r = hi - lo
    range_groups[r].append(k)

print("  RANGE SUBSETS at n=5:")
for r in sorted(range_groups.keys()):
    positions = range_groups[r]
    tiles = [tile_pairs[k] for k in positions]
    print("    Range %d: %d tiles at positions %s = %s" %
          (r, len(positions), positions, tiles))

# For each range subset: compute the waggly weight matrix
print("\n  RANGE-SUBSET WAGGLY MATRICES:\n")
for r in sorted(range_groups.keys()):
    positions = range_groups[r]
    d = len(positions)
    if d < 2: continue  # Need at least 2 to be waggly

    mask = sum(1 << k for k in positions)
    W_range = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        mj = ts['mid_map'][ts['tc'][bits ^ mask]]
        W_range[mi][mj] += 1

    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W_range[mi][mj]>0)
    n_self = sum(1 for mi in range(nm) if W_range[mi][mi]>0)
    sl_rate = sum(W_range[mi][mi] for mi in range(nm)) / total

    print("  RANGE-%d FLIP (flip %d tiles): %d edges, %d self-loops, neutrality=%.3f" %
          (r, d, n_edges, n_self, sl_rate))

# Compute VERTEX-STAR subsets
print("\n  VERTEX-STAR SUBSETS at n=5:")
for v in range(n):
    # Tiles incident to vertex v (where v is one of the pair)
    positions = [k for k, (lo, hi) in enumerate(tile_pairs) if lo == v or hi == v]
    d = len(positions)
    tiles = [tile_pairs[k] for k in positions]
    print("    Vertex %d: %d tiles at %s" % (v, d, tiles))

    if d < 2: continue
    mask = sum(1 << k for k in positions)
    W_star = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        mj = ts['mid_map'][ts['tc'][bits ^ mask]]
        W_star[mi][mj] += 1

    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W_star[mi][mj]>0)
    sl_rate = sum(W_star[mi][mi] for mi in range(nm)) / total
    print("      → %d edges, neutrality=%.3f" % (n_edges, sl_rate))

# Compute TRIANGLE subsets
print("\n  TRIANGLE SUBSETS at n=5:")
print("  (3-cycles in the tile pair graph)\n")
triangles = []
for i in range(n):
    for j in range(i+1, n):
        for k in range(j+1, n):
            # Check if (i,j), (j,k), (i,k) are all tile pairs
            pairs_needed = [(min(i,j),max(i,j)), (min(j,k),max(j,k)), (min(i,k),max(i,k))]
            positions = []
            for p in pairs_needed:
                if p in tile_pairs:
                    positions.append(tile_pairs.index(p))
                else:
                    break
            if len(positions) == 3:
                triangles.append((i,j,k, positions))

for i, j, k, positions in triangles:
    mask = sum(1 << p for p in positions)
    W_tri = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        mj = ts['mid_map'][ts['tc'][bits ^ mask]]
        W_tri[mi][mj] += 1

    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W_tri[mi][mj]>0)
    sl_rate = sum(W_tri[mi][mi] for mi in range(nm)) / total
    print("  Triangle {%d,%d,%d} (tiles %s): %d edges, neutrality=%.3f" %
          (i, j, k, [tile_pairs[p] for p in positions], n_edges, sl_rate))

# ============================================================
banner("THE COMPLETE n=5 WAGGLY DECOMPOSITION")
# ============================================================

print("  SELF-LOOP RATES FOR ALL STRUCTURED SUBSETS:\n")
print("  %-30s %-8s %-8s %-8s" % ("Subset", "#tiles", "edges", "neutral"))
print("  " + "-"*55)

# Individual tiles (wiggly)
for k, (lo, hi) in enumerate(tile_pairs):
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        mj = ts['mid_map'][ts['tc'][bits ^ (1<<k)]]
        W[mi][mj] += 1
    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W[mi][mj]>0)
    sl = sum(W[mi][mi] for mi in range(nm)) / total
    print("  %-30s %-8d %-8d %-8.3f" % ("Tile (%d,%d)" % (lo,hi), 1, n_edges, sl))

print()

# Range flips
for r in sorted(range_groups.keys()):
    positions = range_groups[r]
    d = len(positions)
    mask = sum(1 << k for k in positions)
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        mj = ts['mid_map'][ts['tc'][bits ^ mask]]
        W[mi][mj] += 1
    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W[mi][mj]>0)
    sl = sum(W[mi][mi] for mi in range(nm)) / total
    print("  %-30s %-8d %-8d %-8.3f" % ("Range-%d flip" % r, d, n_edges, sl))

print()

# Vertex star flips
for v in range(n):
    positions = [k for k, (lo, hi) in enumerate(tile_pairs) if lo == v or hi == v]
    d = len(positions)
    if d == 0: continue
    mask = sum(1 << k for k in positions)
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        mj = ts['mid_map'][ts['tc'][bits ^ mask]]
        W[mi][mj] += 1
    n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W[mi][mj]>0)
    sl = sum(W[mi][mi] for mi in range(nm)) / total
    print("  %-30s %-8d %-8d %-8.3f" % ("Star(v=%d)" % v, d, n_edges, sl))

print()

# Triangle flips
for i, j, k, positions in triangles:
    mask = sum(1 << p for p in positions)
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        mj = ts['mid_map'][ts['tc'][bits ^ mask]]
        W[mi][mj] += 1
    n_edges_t = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W[mi][mj]>0)
    sl = sum(W[mi][mi] for mi in range(nm)) / total
    print("  %-30s %-8d %-8d %-8.3f" % ("Tri{%d,%d,%d}" % (i,j,k), 3, n_edges_t, sl))

print()

# Blue/black (complement)
mask = (1 << m) - 1
W = np.zeros((nm, nm), dtype=int)
for bits in range(total):
    mi = ts['mid_map'][ts['tc'][bits]]
    mj = ts['mid_map'][ts['tc'][bits ^ mask]]
    W[mi][mj] += 1
n_edges = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if W[mi][mj]>0)
sl = sum(W[mi][mi] for mi in range(nm)) / total
print("  %-30s %-8d %-8d %-8.3f" % ("COMPLEMENT (all tiles)", m, n_edges, sl))

# ============================================================
banner("WHICH SUBSETS ARE MOST DISTINCTIVE?")
# ============================================================

print("""
  NEUTRALITY RANKINGS (fraction of tilings where the flip is silent):

  The MOST NEUTRAL subsets are those where flipping preserves the iso class
  most often. These correspond to "hidden symmetries" of tournaments.

  The LEAST NEUTRAL subsets are the most "disruptive" mutations —
  they almost always change the iso class.

  MEANINGFUL PATTERNS TO LOOK FOR:
  - Do range-2 flips have higher neutrality than range-3?
    (= are "closer" arcs more interchangeable?)
  - Do vertex-star flips have uniform neutrality across vertices?
    (= are all vertices equally "important" in the tiling model?)
  - Do triangle flips involving base-path-adjacent vertices differ
    from those involving far-apart vertices?
""")

print("Done. opus-2026-03-24-S296")
