#!/usr/bin/env python3
"""
three_graphs_corrected_s295b.py -- opus-2026-03-24-S295b

CORRECTED THREE-GRAPH ANALYSIS (MISTAKE-033 fix)

Blue/black lines connect tiling t to its COMPLEMENT TILING:
  t_comp = t XOR (2^m - 1)

This flips all m tile bits. It stays INSIDE Q_m.
It does NOT flip base-path arcs, so it is NOT T^op.

The complement tiling gives a tournament where:
- Base-path arcs: PRESERVED (still n-1→n-2→...→1→0)
- Tile arcs: ALL REVERSED

This is a different tournament, potentially in a different iso class.

Author: opus-2026-03-24-S295b
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
    # Tile pairs: non-consecutive vertices (skip ≥ 2)
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    cmap = {}; sizes = {}; tc = {}; members = defaultdict(list)
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1  # base path: higher→lower
        for k, (lo, hi) in enumerate(tile_pairs):
            # pair (lo, hi) with lo < hi
            if (bits >> k) & 1:
                A[lo][hi] = 1  # lo beats hi (flipped)
            else:
                A[hi][lo] = 1  # hi beats lo (natural)
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
            merged[mid] = (cid, comp)
            mid_map[cid] = mid; mid_map[comp] = mid
        else:
            merged[mid] = (cid,)
            mid_map[cid] = mid
        mid += 1
    return {'n': n, 'm': m, 'total': total, 'tile_pairs': tile_pairs,
            'tc': tc, 'sizes': sizes, 'members': members,
            'info': info, 'merged': merged, 'mid_map': mid_map,
            'nm': mid, 'nc': nc}

# ============================================================
banner("CORRECTED THREE-GRAPH ANALYSIS")
# ============================================================

for n in [4, 5, 6]:
    ts = build_tiling_space(n)
    m = ts['m']; nm = ts['nm']; total = ts['total']
    mask = (1 << m) - 1  # all-ones mask for m bits

    banner("n=%d: m=%d tiles, %d tilings, %d merged classes" % (n, m, total, nm))

    order = sorted(range(nm), key=lambda x: ts['info'][ts['merged'][x][0]]['H'])
    H_vec = [ts['info'][ts['merged'][mi][0]]['H'] for mi in order]
    fiber_vec = [sum(ts['sizes'][c] for c in ts['merged'][mi]) for mi in order]

    # GRAPH 1: WIGGLY (d=1)
    W_wig = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        for k in range(m):
            mj = ts['mid_map'][ts['tc'][bits ^ (1 << k)]]
            W_wig[mi][mj] += 1

    # GRAPH 2: WAGGLY-2 (d=2)
    W_wag2 = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        for j in range(m):
            for k in range(j+1, m):
                mj = ts['mid_map'][ts['tc'][bits ^ (1<<j) ^ (1<<k)]]
                W_wag2[mi][mj] += 1

    # BLUE/BLACK (CORRECT): complement TILING = XOR all m bits
    W_bb = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        comp_bits = bits ^ mask  # flip ALL tiles
        mj = ts['mid_map'][ts['tc'][comp_bits]]
        W_bb[mi][mj] += 1

    # Print all three
    for name, W in [("WIGGLY (d=1)", W_wig),
                     ("WAGGLY-2 (d=2)", W_wag2),
                     ("BLUE/BLACK (tiling complement, d=%d)" % m, W_bb)]:
        n_edges = sum(1 for mi in range(nm) for mj in range(mi+1, nm) if W[mi][mj] > 0)
        n_self = sum(1 for mi in range(nm) if W[mi][mi] > 0)
        total_w = W.sum()
        sl_w = sum(W[mi][mi] for mi in range(nm))

        print("  %s:" % name)
        print("    Edges: %d, Self-loops: %d/%d, Total: %d, SL%%: %.1f%%" %
              (n_edges, n_self, nm, total_w, 100*sl_w/total_w if total_w > 0 else 0))
        print("    H:", H_vec)
        for mi in order:
            row = " ".join("%4d" % W[mi][mj] for mj in order)
            print("    H=%2d: %s  |fib=%d" % (H_vec[order.index(mi)], row, fiber_vec[order.index(mi)]))
        print()

    # VENN DIAGRAM
    adj_wig = np.array([[1 if mi != mj and W_wig[mi][mj] > 0 else 0
                         for mj in range(nm)] for mi in range(nm)])
    adj_wag = np.array([[1 if mi != mj and W_wag2[mi][mj] > 0 else 0
                         for mj in range(nm)] for mi in range(nm)])
    adj_bb = np.array([[1 if mi != mj and W_bb[mi][mj] > 0 else 0
                        for mj in range(nm)] for mi in range(nm)])

    e_wig = adj_wig.sum()//2; e_wag = adj_wag.sum()//2; e_bb = adj_bb.sum()//2

    both_ww = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_wig[mi][mj] and adj_wag[mi][mj])
    both_wb = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_wig[mi][mj] and adj_bb[mi][mj])
    both_bw = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_bb[mi][mj] and adj_wag[mi][mj])
    all3 = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_wig[mi][mj] and adj_wag[mi][mj] and adj_bb[mi][mj])
    union = sum(1 for mi in range(nm) for mj in range(mi+1,nm) if adj_wig[mi][mj] or adj_wag[mi][mj] or adj_bb[mi][mj])

    print("  VENN DIAGRAM:")
    print("    Wiggly edges:       %d" % e_wig)
    print("    Waggly-2 edges:     %d" % e_wag)
    print("    Blue/black edges:   %d" % e_bb)
    print("    Wig ∩ Wag:          %d" % both_ww)
    print("    Wig ∩ BB:           %d" % both_wb)
    print("    Wag ∩ BB:           %d" % both_bw)
    print("    All three:          %d" % all3)
    print("    Union:              %d / %d possible" % (union, nm*(nm-1)//2))
    print()

    # Unique edges
    only_wig = sum(1 for mi in range(nm) for mj in range(mi+1,nm)
                   if adj_wig[mi][mj] and not adj_wag[mi][mj] and not adj_bb[mi][mj])
    only_wag = sum(1 for mi in range(nm) for mj in range(mi+1,nm)
                   if adj_wag[mi][mj] and not adj_wig[mi][mj] and not adj_bb[mi][mj])
    only_bb = sum(1 for mi in range(nm) for mj in range(mi+1,nm)
                  if adj_bb[mi][mj] and not adj_wig[mi][mj] and not adj_wag[mi][mj])

    print("  UNIQUE EDGES:")
    print("    Only wiggly:     %d" % only_wig)
    print("    Only waggly-2:   %d" % only_wag)
    print("    Only blue/black: %d" % only_bb)

    # List blue/black cross-class edges
    bb_cross = [(mi, mj) for mi in range(nm) for mj in range(mi+1, nm) if adj_bb[mi][mj]]
    if bb_cross:
        print("\n  BLUE/BLACK CROSS-CLASS EDGES:")
        for mi, mj in bb_cross:
            Hi = ts['info'][ts['merged'][mi][0]]['H']
            Hj = ts['info'][ts['merged'][mj][0]]['H']
            w = W_bb[mi][mj] + W_bb[mj][mi]
            in_wig = "also wiggly" if adj_wig[mi][mj] else ""
            print("    {%d(H=%d) ↔ %d(H=%d)}: weight=%d %s" % (mi, Hi, mj, Hj, w, in_wig))

    # Blue/black self-loop analysis
    print("\n  BLUE/BLACK SELF-LOOPS (SC classes in tiling model):")
    for mi in order:
        if W_bb[mi][mi] > 0:
            H = ts['info'][ts['merged'][mi][0]]['H']
            fib = fiber_vec[order.index(mi)]
            print("    mid=%d (H=%d): %d self-loops out of %d tilings (%.1f%%)" %
                  (mi, H, W_bb[mi][mi], fib, 100*W_bb[mi][mi]/fib))

    # How does blue/black relate to wiggly and waggly?
    print("\n  BLUE/BLACK vs WIGGLY OVERLAP:")
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if adj_bb[mi][mj]:
                Hi = ts['info'][ts['merged'][mi][0]]['H']
                Hj = ts['info'][ts['merged'][mj][0]]['H']
                in_wig = adj_wig[mi][mj]
                in_wag = adj_wag[mi][mj]
                print("    {%d(H=%d)↔%d(H=%d)}: wiggly=%s, waggly-2=%s" %
                      (mi, Hi, mj, Hj, bool(in_wig), bool(in_wag)))

# ============================================================
banner("KEY IDENTITIES AND RELATIONSHIPS")
# ============================================================

print("""
  CORRECTED DEFINITIONS:

  1. WIGGLY LINE: flip 1 tile. Hamming distance 1 in Q_m.
     = single-tile mutation. Each tiling has exactly m neighbors.

  2. WAGGLY LINE: any non-wiggly connection in Q_m. Distance ≥ 2.
     Includes the complement tiling at distance m.

  3. BLUE/BLACK LINE: connect tiling t to its complement tiling
     t_comp = t XOR (2^m - 1). This is a SPECIFIC waggly line at
     distance m. It flips all TILES but preserves the BASE PATH.

     This is NOT the tournament complement T^op!
     T^op would require also flipping the base-path arcs.
     The complement tiling gives a DIFFERENT tournament than T^op.

  CRITICAL: blue/black IS a waggly line (at max distance m).
  It IS inside Q_m. It DOES connect different iso classes.

  The previous analysis (S295) was WRONG because it computed T^op
  (flipping all arcs) instead of the complement tiling (flipping
  only non-base-path arcs).

  CORRECT HIERARCHY:
  wiggly ⊂ Q_m edges (d=1)
  waggly ⊃ blue/black (d=2,...,m including d=m)
  blue/black = the d=m waggly line = complement tiling

  Blue/black has cross-class edges. They are the connections between
  classes that are NOT self-complementary in the TILING sense
  (which is different from SC in the TOURNAMENT sense).
""")

print("Done. opus-2026-03-24-S295b")
