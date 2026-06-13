#!/usr/bin/env python3
"""
wiggly_vs_lines_s275.py -- opus-2026-03-24-S275

THE THREE TYPES OF TILING CONNECTIONS (CORRECTED)

In the tournament-tiling-explorer, tilings are connected by:

1. BLUE/BLACK LINES: Connect each tiling to its COMPLEMENT.
   = flip ALL bits simultaneously. Each tiling has exactly 1 line partner.
   BLUE if tiling is grid-symmetric, BLACK if not.
   These form a PERFECT MATCHING on the tiling space.

2. WIGGLY LINES: Connect each tiling to its single-tile-flip neighbors.
   = flip ONE bit. Each tiling has m = C(n,2) wiggly neighbors.
   Wiggly class A,B,C,... = which tile was flipped.
   These form the m-CUBE skeleton Q_m.

3. The MERGED META-GRAPH: Merge complement pairs (blue/black lines) into
   single nodes. Then the wiggly lines BETWEEN merged nodes give the
   meta-graph edges.

KEY INSIGHT: Blue/black lines and wiggly lines are DIFFERENT operations!
  Blue/black: complement (flip all bits) — a single partner per tiling
  Wiggly: single-tile flip — m partners per tiling

A wiggly line from tiling T to T' (single flip) may connect:
  (a) T and T' in the SAME merged node → self-loop (neutral arc)
  (b) T and T' in DIFFERENT merged nodes → meta-graph edge

The meta-graph edges come ENTIRELY from wiggly lines.
Blue/black lines contribute ONLY to the merging structure (Z_2 quotient).

QUESTION: Which merged-node pairs are connected by wiggly lines but NOT
by blue/black lines? Answer: ALMOST ALL of them, because blue/black only
connect complement pairs (within or between merged nodes).

Author: opus-2026-03-24-S275
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
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

def is_grid_sym(bits, m):
    """Check if tiling is grid-symmetric (complement = same tiling under vertex reversal)."""
    # Grid-symmetric: the tiling is invariant under (x,y) -> (n-y+1, n-x+1)
    # In bit representation: this reverses and complements
    return bits == ((2**m - 1) ^ bits)  # This is complement = flip all bits
    # Wait: grid-sym means the tiling EQUALS its complement under the grid reflection.
    # This is NOT the same as self-complement.
    # Need to implement properly. For now, skip.

for n in [4, 5]:
    banner("WIGGLY vs BLUE/BLACK LINES AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build iso classes
    cmap = {}; sizes = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        sizes[cmap[c]] += 1; tc[bits] = cmap[c]
    nc = len(cmap)

    info = {}
    for cid in range(nc):
        bits_rep = next(b for b in range(2**m) if tc[b] == cid)
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits_rep>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        comp_c = canon(complement(A, n), n)
        info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                     'SC': cmap[comp_c] == cid}

    # Merged nodes
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    print("  V_unmerged = %d, V_merged = %d, m = %d\n" % (nc, nm, m))

    # ============================================================
    # BLUE/BLACK LINES: complement pairing
    # ============================================================
    # Each tiling bits has complement (2^m - 1) ^ bits = all-flip.
    # The complement tiling maps to a (possibly different) iso class.
    # The LINE connects bits to its complement.

    # At the MERGED level: a blue/black line connects merged nodes mi and mj
    # where mi contains bits and mj contains complement(bits).
    # Since complement maps cid to info[cid]['comp'], and the merged node
    # ALREADY merges cid with comp: the blue/black line connects mi to ITSELF
    # (within the same merged node) for EVERY tiling!

    # Wait — that's the POINT of merging. Complement pairs are merged into
    # one node. So blue/black lines are ALL self-loops in the merged graph.

    # But some complement pairs map to the SAME unmerged class (SC classes).
    # For SC classes: complement = same class, so the blue/black line
    # stays within the same unmerged class too.

    # For non-SC classes: complement maps to a different unmerged class,
    # but both are in the SAME merged node. So still a self-loop.

    # CONCLUSION: In the merged meta-graph, blue/black lines are ALL
    # internal to nodes. They contribute ZERO inter-node edges.

    print("  BLUE/BLACK LINES (complement pairing):")
    print("  Every line connects a tiling to its complement.")
    print("  In the merged graph: ALL blue/black lines are SELF-LOOPS")
    print("  (because complement pairs are merged into the same node).\n")

    # ============================================================
    # WIGGLY LINES: single-tile flips
    # ============================================================
    # Each tiling bits connects to m neighbors bits ^ (1<<k) for k=0..m-1.
    # Each neighbor maps to some merged class.
    # If neighbor's merged class ≠ self's merged class → meta-graph edge.

    wiggly_edges = set()  # (mi, mj) edges from wiggly
    wiggly_self = 0
    wiggly_cross = 0

    for bits in range(2**m):
        mi = mid_map[tc[bits]]
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            if mi == mj:
                wiggly_self += 1
            else:
                wiggly_cross += 1
                wiggly_edges.add((min(mi,mj), max(mi,mj)))

    E_wiggly = len(wiggly_edges)

    print("  WIGGLY LINES (single-tile flips):")
    print("  Total wiggly lines: %d (= 2^m × m / 2 = %d)" %
          ((wiggly_self + wiggly_cross) // 2, 2**m * m // 2))
    print("  Self-loops (same merged node): %d" % (wiggly_self // 2))
    print("  Cross-edges (different merged nodes): %d" % (wiggly_cross // 2))
    print("  Distinct merged-node pairs connected: %d\n" % E_wiggly)

    # ============================================================
    # COMPARISON: what does each type connect?
    # ============================================================

    print("  COMPARISON:")
    print("  Blue/black lines: 0 inter-node edges (all self-loops)")
    print("  Wiggly lines: %d inter-node edges" % E_wiggly)
    print("  The meta-graph has E = %d edges (ALL from wiggly lines)" % E_wiggly)
    print()
    print("  Blue/black lines connect TILINGS (one-to-one complement pairing)")
    print("  Wiggly lines connect TILINGS (one-to-m single-flip)")
    print("  Meta-graph edges = wiggly lines projected to merged nodes")
    print()

    # ============================================================
    # THE KEY: wiggly-only edges (not reachable by blue/black)
    # ============================================================

    # Blue/black never creates inter-node edges in the merged graph.
    # So ALL meta-graph edges are wiggly-only.
    # There are NO edges that are BOTH blue/black and wiggly.

    # But in the UNMERGED graph, we can ask:
    # Which unmerged class pairs are connected by wiggly but not by blue/black?

    # Blue/black connects each unmerged class to its complement class.
    bb_edges = set()
    for cid in range(nc):
        comp_cid = info[cid]['comp']
        if cid != comp_cid:
            bb_edges.add((min(cid, comp_cid), max(cid, comp_cid)))

    # Wiggly connects unmerged classes that differ by one arc flip
    wig_unmerged = set()
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            cj = tc[bits ^ (1 << k)]
            if ci != cj:
                wig_unmerged.add((min(ci, cj), max(ci, cj)))

    both = bb_edges & wig_unmerged
    bb_only = bb_edges - wig_unmerged
    wig_only = wig_unmerged - bb_edges

    print("  UNMERGED GRAPH:")
    print("  Blue/black edges (complement pairs): %d" % len(bb_edges))
    print("  Wiggly edges (single-flip pairs): %d" % len(wig_unmerged))
    print("  Both: %d" % len(both))
    print("  Blue/black only: %d" % len(bb_only))
    print("  Wiggly only: %d" % len(wig_only))

    if both:
        print("\n  Edges that are BOTH blue/black AND wiggly:")
        for (ci, cj) in sorted(both):
            Hi = info[ci]['H']; Hj = info[cj]['H']
            print("    {%d(H=%d), %d(H=%d)}: complement pair AND single-flip neighbors" %
                  (ci, Hi, cj, Hj))
        print("  These are class pairs where the complement is ALSO reachable by 1 flip!")
    else:
        print("\n  NO edges are both blue/black and wiggly.")
        print("  The complement is NEVER reachable by a single arc flip.")

banner("DEFINITIVE TERMINOLOGY UPDATE")

print("""
  THE THREE TYPES OF TILING CONNECTIONS:

  1. BLUE LINE: tiling → its complement, where the tiling is grid-symmetric.
     = flip ALL bits simultaneously.
     Connects each tiling to exactly 1 partner.
     In the merged meta-graph: ALL blue lines are internal (self-loops).

  2. BLACK LINE: tiling → its complement, where the tiling is NOT grid-symmetric.
     Same as blue but for non-grid-symmetric tilings.
     Also internal to merged nodes.

  3. WIGGLY LINE (class X): tiling → its single-tile-flip neighbor at cell X.
     = flip ONE bit.
     Each tiling has m = C(n,2) wiggly lines (one per cell).
     Wiggly class A,B,C,... labels which cell.
     In the merged meta-graph: wiggly lines create ALL the edges.

  THE RELATIONSHIP:
  - Blue/black lines create the MERGING (Z_2 quotient)
  - Wiggly lines create the EDGES (meta-graph structure)
  - These are DIFFERENT operations on different scales
  - A complement is almost NEVER reachable by a single flip
    (verified: 0 cases at n=4, need to check n=5,6)

  PREVIOUS ERRORS TO CORRECT:
  - "translucent" → DEPRECATED (was confusing)
  - The meta-graph edges were sometimes called "blue/black" → WRONG
    Meta-graph edges come from WIGGLY lines, not from complement pairing
  - "opaque" → redundant, since ALL inter-node connections are wiggly
""")

print("Done. opus-2026-03-24-S275")
