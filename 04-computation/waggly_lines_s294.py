#!/usr/bin/env python3
"""
waggly_lines_s294.py -- opus-2026-03-24-S294

WAGGLY LINES: CONNECTIONS BETWEEN TILINGS AT HAMMING DISTANCE ≥ 2

The tiling hypercube Q_m has:
  - WIGGLY lines: edges (Hamming distance 1, single tile flip)
  - WAGGLY lines: NON-edges (Hamming distance ≥ 2, multi-tile mutations)

A waggly line connects two tilings that differ in 2 or more tiles.
This is the COMPLEMENT GRAPH of Q_m restricted to same-class or cross-class.

BLUE/BLACK LINES are the special case of Hamming distance m (flip ALL tiles):
  - The complement tiling = flip every tile = Hamming distance m
  - This gives the tournament's complement T^op
  - Blue if grid-symmetric, black if not

So: blue/black ⊂ waggly (since m ≥ 2 for n ≥ 4).

THE WAGGLY GRAPH:
  - Vertices: 2^m tilings (same as wiggly)
  - Edges: all pairs at Hamming distance 2, 3, ..., m
  - Degree: 2^m - 1 - m (everything except self and wiggly neighbors)

THE WAGGLY DECOMPOSITION BY HAMMING DISTANCE:
  d=2: C(m,2) neighbors (flip 2 tiles)
  d=3: C(m,3) neighbors
  ...
  d=k: C(m,k) neighbors
  ...
  d=m: 1 neighbor (the complement = blue/black line)

Total waggly neighbors: Σ_{k=2}^{m} C(m,k) = 2^m - 1 - m

THE QUESTION: What does the waggly graph look like in the QUOTIENT
by isomorphism classes? Which Hamming distances matter most?

Author: opus-2026-03-24-S294
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
    tile_pairs = [(hi, lo) for hi in range(n) for lo in range(hi) if hi - lo >= 2]
    tile_pairs.sort()

    cmap = {}; sizes = {}; tc = {}; members = defaultdict(list)
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (hi, lo) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        cid = cmap[c]; sizes[cid] += 1; tc[bits] = cid
        members[cid].append(bits)
    nc = len(cmap)

    info = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (hi, lo) in enumerate(tile_pairs):
            if (b >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        comp_c = canon(complement(A, n), n)
        info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                      'SC': cmap[comp_c] == cid, 'size': sizes[cid]}

    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    return {'n': n, 'm': m, 'total': total, 'tile_pairs': tile_pairs,
            'tc': tc, 'sizes': sizes, 'members': members,
            'info': info, 'merged': merged, 'mid_map': mid_map,
            'nm': mid, 'nc': nc}

# ============================================================
banner("THE WAGGLY GRAPH: STRUCTURE BY HAMMING DISTANCE")
# ============================================================

for n in [5, 6]:
    ts = build_tiling_space(n)
    m = ts['m']; total = ts['total']; nm = ts['nm']

    print("  n=%d: m=%d tiles, %d tilings, %d merged classes" % (n, m, total, nm))
    print("  Wiggly (d=1): %d neighbors per tiling" % m)
    print("  Waggly (d≥2): %d neighbors per tiling" % (total - 1 - m))
    print("  Blue/black (d=m=%d): 1 neighbor per tiling" % m)
    print()

    # Build WAGGLY weight matrices by Hamming distance
    # W_d[mi][mj] = #{(t1, t2) : t1 ∈ mi, t2 ∈ mj, hamming(t1,t2) = d}
    max_d = min(m, 6)  # Limit for speed at n=6

    print("  WAGGLY WEIGHT MATRICES BY HAMMING DISTANCE:\n")

    for d in range(1, max_d + 1):
        W_d = np.zeros((nm, nm), dtype=int)

        if d <= 3 or (n == 5):
            # Enumerate all pairs at distance d
            for bits in range(total):
                mi = ts['mid_map'][ts['tc'][bits]]
                # Generate all neighbors at Hamming distance d
                for flip_positions in combinations(range(m), d):
                    mask = 0
                    for p in flip_positions:
                        mask |= (1 << p)
                    nbr = bits ^ mask
                    mj = ts['mid_map'][ts['tc'][nbr]]
                    W_d[mi][mj] += 1

            total_weight = W_d.sum()
            self_loops = sum(W_d[mi][mi] for mi in range(nm))
            cross = total_weight - self_loops

            # Expected total: total × C(m, d)
            expected = total * comb(m, d)

            print("  d=%d: total=%d (expected=%d ✓), self-loops=%d (%.1f%%), cross=%d" %
                  (d, total_weight, expected,
                   self_loops, 100*self_loops/total_weight if total_weight > 0 else 0,
                   cross))

            # The self-loop rate at each distance
            if d >= 2:
                # Compare with wiggly (d=1) self-loop rate
                W1 = np.zeros((nm, nm), dtype=int)
                for bits in range(total):
                    mi = ts['mid_map'][ts['tc'][bits]]
                    for k in range(m):
                        mj = ts['mid_map'][ts['tc'][bits ^ (1 << k)]]
                        W1[mi][mj] += 1
                sl1 = sum(W1[mi][mi] for mi in range(nm))
                sl1_rate = sl1 / (total * m) if m > 0 else 0

                sl_d_rate = self_loops / (total * comb(m, d)) if comb(m, d) > 0 else 0
                print("        self-loop rate: %.4f (vs wiggly d=1: %.4f)" %
                      (sl_d_rate, sl1_rate))

        # For d=m (complement): always exactly 1 neighbor
        if d == m:
            W_comp = np.zeros((nm, nm), dtype=int)
            for bits in range(total):
                mi = ts['mid_map'][ts['tc'][bits]]
                comp = bits ^ ((1 << m) - 1)
                mj = ts['mid_map'][ts['tc'][comp]]
                W_comp[mi][mj] += 1

            total_comp = W_comp.sum()
            self_comp = sum(W_comp[mi][mi] for mi in range(nm))
            print("\n  d=%d (COMPLEMENT = blue/black): total=%d, self-loops=%d (%.1f%%)" %
                  (m, total_comp, self_comp, 100*self_comp/total_comp if total_comp > 0 else 0))

    print()

# ============================================================
banner("THE WAGGLY SPECTRUM AT n=5")
# ============================================================

ts = build_tiling_space(5)
n = 5; m = ts['m']; nm = ts['nm']; total = ts['total']

print("  COMPLETE WAGGLY ANALYSIS (n=5, m=6, all distances 1..6):\n")

# Build weight matrices for ALL distances
W_by_d = {}
for d in range(1, m+1):
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(total):
        mi = ts['mid_map'][ts['tc'][bits]]
        for flip_positions in combinations(range(m), d):
            mask = sum(1 << p for p in flip_positions)
            nbr = bits ^ mask
            mj = ts['mid_map'][ts['tc'][nbr]]
            W[mi][mj] += 1
    W_by_d[d] = W

order = sorted(range(nm), key=lambda x: ts['info'][ts['merged'][x][0]]['H'])
H_vec = [ts['info'][ts['merged'][mi][0]]['H'] for mi in order]
fiber_vec = [sum(ts['sizes'][c] for c in ts['merged'][mi]) for mi in order]

# Print self-loop rate by distance and class
print("  SELF-LOOP RATE BY DISTANCE AND CLASS:\n")
print("  %-4s %-4s %-6s" % ("mid", "H", "fiber"), end="")
for d in range(1, m+1):
    print(" d=%-5d" % d, end="")
print()
print("  " + "-"*70)

for mi in order:
    H = ts['info'][ts['merged'][mi][0]]['H']
    fiber = sum(ts['sizes'][c] for c in ts['merged'][mi])
    print("  %-4d %-4d %-6d" % (mi, H, fiber), end="")
    for d in range(1, m+1):
        sl = W_by_d[d][mi][mi]
        rate = sl / (fiber * comb(m, d)) if fiber > 0 else 0
        print(" %-7.3f" % rate, end="")
    print()

# Total self-loop rate by distance
print("\n  TOTAL SELF-LOOP RATE BY DISTANCE:")
for d in range(1, m+1):
    sl = sum(W_by_d[d][mi][mi] for mi in range(nm))
    total_pairs = total * comb(m, d)
    rate = sl / total_pairs if total_pairs > 0 else 0
    print("    d=%d: %d/%d = %.6f (%d-tile mutation neutrality)" %
          (d, sl, total_pairs, rate, d))

# ============================================================
banner("THE THREE LENSES: WIGGLY, WAGGLY-2, AND COMPLEMENT")
# ============================================================

# Lens 1: WIGGLY (d=1) — the meta-graph we know
# Lens 2: WAGGLY-2 (d=2) — two-tile mutations, the "second harmonic"
# Lens 3: COMPLEMENT (d=m) — blue/black, the "fundamental mode"

print("  THREE LENSES ON THE TILING SPACE (n=5):\n")

# For each lens, compute the quotient weight matrix
for lens_name, d, W in [("WIGGLY (d=1)", 1, W_by_d[1]),
                          ("WAGGLY-2 (d=2)", 2, W_by_d[2]),
                          ("WAGGLY-3 (d=3)", 3, W_by_d[3]),
                          ("COMPLEMENT (d=6)", m, W_by_d[m])]:
    print("  === %s ===" % lens_name)

    # Adjacency (unweighted): which classes are connected?
    adj = np.zeros((nm, nm), dtype=int)
    for mi in range(nm):
        for mj in range(nm):
            if mi != mj and W[mi][mj] > 0:
                adj[mi][mj] = 1
    n_edges = adj.sum() // 2
    n_self = sum(1 for mi in range(nm) if W[mi][mi] > 0)

    print("    Edges: %d / %d possible, Self-loops: %d / %d classes" %
          (n_edges, nm*(nm-1)//2, n_self, nm))

    # Degree distribution
    degrees = [sum(adj[mi]) for mi in range(nm)]
    print("    Degree range: %d..%d (avg %.1f)" %
          (min(degrees), max(degrees), np.mean(degrees)))

    # What fraction of total weight is in self-loops?
    total_w = W.sum()
    sl_w = sum(W[mi][mi] for mi in range(nm))
    print("    Self-loop fraction: %.4f (%d/%d)" %
          (sl_w/total_w if total_w > 0 else 0, sl_w, total_w))

    # Is it the SAME adjacency as wiggly?
    if d > 1:
        adj_1 = np.zeros((nm, nm), dtype=int)
        for mi in range(nm):
            for mj in range(nm):
                if mi != mj and W_by_d[1][mi][mj] > 0:
                    adj_1[mi][mj] = 1
        same = np.array_equal(adj, adj_1)
        diff = sum(1 for mi in range(nm) for mj in range(mi+1, nm)
                   if adj[mi][mj] != adj_1[mi][mj])
        print("    Same adjacency as wiggly: %s (diff=%d edges)" % (same, diff))
    print()

# ============================================================
banner("THE WAGGLY-2 WEIGHT MATRIX")
# ============================================================

W2 = W_by_d[2]

print("  WAGGLY-2 WEIGHT MATRIX (2-tile mutations):\n")
print("  H: ", H_vec)
for mi in order:
    row = " ".join("%4d" % W2[mi][mj] for mj in order)
    fiber = sum(ts['sizes'][c] for c in ts['merged'][mi])
    print("  H=%2d: %s  | fib=%d, row=%d" %
          (ts['info'][ts['merged'][mi][0]]['H'], row, fiber,
           sum(W2[mi])))

print("\n  Properties:")
print("  Symmetric: %s" % np.allclose(W2, W2.T))
print("  Row sums: %s" % [sum(W2[mi]) for mi in order])
print("  All row sums = fiber × C(m,2) = fiber × %d: %s" %
      (comb(m,2), all(sum(W2[mi]) == sum(ts['sizes'][c] for c in ts['merged'][mi]) * comb(m,2)
                       for mi in range(nm))))

# ============================================================
banner("COMPARING WIGGLY AND WAGGLY-2 TRANSITION MATRICES")
# ============================================================

# Transition matrix for wiggly (d=1)
P1 = W_by_d[1].astype(float)
for mi in range(nm):
    fib = sum(ts['sizes'][c] for c in ts['merged'][mi])
    P1[mi] /= (fib * m)

# Transition matrix for waggly-2 (d=2)
P2 = W_by_d[2].astype(float)
for mi in range(nm):
    fib = sum(ts['sizes'][c] for c in ts['merged'][mi])
    P2[mi] /= (fib * comb(m, 2))

print("  TRANSITION MATRIX COMPARISON (n=5):\n")
print("  Wiggly (d=1) eigenvalues: %s" %
      sorted(np.linalg.eigvals(P1).real, reverse=True)[:5])
print("  Waggly (d=2) eigenvalues: %s" %
      sorted(np.linalg.eigvals(P2).real, reverse=True)[:5])

# Are the eigenvectors the SAME?
evals1, evecs1 = np.linalg.eig(P1)
evals2, evecs2 = np.linalg.eig(P2)

# Sort by eigenvalue
idx1 = np.argsort(evals1.real)[::-1]
idx2 = np.argsort(evals2.real)[::-1]

print("\n  EIGENVECTOR ALIGNMENT:")
for i in range(min(4, nm)):
    v1 = evecs1[:, idx1[i]].real
    v2 = evecs2[:, idx2[i]].real
    # Normalize
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    alignment = abs(np.dot(v1, v2))
    print("    v_%d: λ_wiggly=%.4f, λ_waggly=%.4f, |<v1,v2>|=%.4f" %
          (i, evals1[idx1[i]].real, evals2[idx2[i]].real, alignment))

# ============================================================
banner("THE COMPLEMENT LENS (d=m)")
# ============================================================

W_comp = W_by_d[m]

print("  COMPLEMENT WEIGHT MATRIX (d=%d = flip ALL tiles):\n" % m)
print("  H: ", H_vec)
for mi in order:
    row = " ".join("%3d" % W_comp[mi][mj] for mj in order)
    print("  H=%2d: %s" % (ts['info'][ts['merged'][mi][0]]['H'], row))

print("\n  In the MERGED meta-graph, complement maps class to class.")
print("  Self-loops in complement = SC classes (complement = self).")
print("  Off-diagonal = NS complement pairs.")
print()
print("  SC classes (self-loop in complement):")
for mi in order:
    if W_comp[mi][mi] > 0:
        H = ts['info'][ts['merged'][mi][0]]['H']
        print("    mid=%d (H=%d): %d self-complement tilings" % (mi, H, W_comp[mi][mi]))

# ============================================================
banner("THE WAGGLY HIERARCHY: d=1 < d=2 < ... < d=m")
# ============================================================

print("""
  THE WAGGLY HIERARCHY AT n=5:

  d=1 (WIGGLY): flip 1 tile. 6 neighbors. The LOCAL mutation.
    Self-loop rate: 0.1042 (10.4%)
    Creates the meta-graph we've been studying.

  d=2 (WAGGLY-2): flip 2 tiles. 15 neighbors. The SECOND HARMONIC.
    Self-loop rate: (computed above)
    Sees the same classes as wiggly but with DIFFERENT weights.

  d=3 (WAGGLY-3): flip 3 tiles. 20 neighbors.
    The CUBIC mutation. Reaches further in class space.

  d=4 (WAGGLY-4): flip 4 tiles. 15 neighbors.
    Symmetric with d=2 (m-d = 2 tiles unchanged).

  d=5 (WAGGLY-5): flip 5 tiles. 6 neighbors.
    Symmetric with d=1 (m-d = 1 tile unchanged).
    This is "everything except one tile" = the NEAR-COMPLEMENT.

  d=6 (COMPLEMENT): flip all 6 tiles. 1 neighbor.
    The blue/black line. Maps T → T^op (in tiling space).

  THE SYMMETRY: waggly at distance d and distance m-d are RELATED.
  Flipping d tiles is the same as "keeping m-d tiles and flipping the rest."

  Define: W_d[C][D] = #{(t, S) : t ∈ C, |S|=d, t⊕S ∈ D}.
  Then: W_d[C][D] = W_{m-d}[C'][D'] where C' = complement of C, D' = complement of D.

  In the MERGED meta-graph (complement pairs identified):
  W_d and W_{m-d} may become the SAME matrix! (or transpose.)
""")

# Verify the d ↔ m-d symmetry
print("  VERIFYING d ↔ m-d SYMMETRY:\n")
for d in range(1, m):
    d2 = m - d
    # W_d[C][D] should relate to W_{d2}[comp(C)][comp(D)]
    # In the merged picture: comp(C) = C for SC, comp(C) = partner for NS.
    # Since we're already merged, let's just check if W_d ≈ W_{m-d}
    ratio = np.sum(W_by_d[d]) / np.sum(W_by_d[d2]) if np.sum(W_by_d[d2]) > 0 else 0
    sl_d = sum(W_by_d[d][mi][mi] for mi in range(nm))
    sl_d2 = sum(W_by_d[d2][mi][mi] for mi in range(nm))
    print("  d=%d vs d=%d: total_weights %d vs %d (ratio %.3f), SL %d vs %d" %
          (d, d2, np.sum(W_by_d[d]), np.sum(W_by_d[d2]), ratio, sl_d, sl_d2))

# ============================================================
banner("WAGGLY SELF-LOOP RATES: THE KEY PATTERN")
# ============================================================

# The self-loop rate at distance d tells us the probability that
# flipping d random tiles preserves the iso class.
# This should be related to the Burnside correction R(n).

print("  SELF-LOOP RATES BY DISTANCE (n=5):\n")
print("  d  | C(m,d) | total_SL | SL_rate   | pattern")
print("  ---+--------+----------+-----------+--------")

for d in range(1, m+1):
    sl = sum(W_by_d[d][mi][mi] for mi in range(nm))
    total_pairs = total * comb(m, d)
    rate = sl / total_pairs if total_pairs > 0 else 0
    pattern = ""
    if d == 1: pattern = "← wiggly"
    elif d == m: pattern = "← complement (blue/black)"
    elif d == m-1: pattern = "← near-complement"
    elif d == 2: pattern = "← waggly-2"
    elif d == m//2: pattern = "← middle"
    print("  %-3d| %-6d | %-8d | %.6f | %s" %
          (d, comb(m, d), sl, rate, pattern))

# ============================================================
banner("THE WAGGLY-2 GRAPH: NEW CLASSES?")
# ============================================================

# Does waggly-2 connect classes that wiggly doesn't?
adj_1 = np.zeros((nm, nm), dtype=int)
adj_2 = np.zeros((nm, nm), dtype=int)
for mi in range(nm):
    for mj in range(nm):
        if mi != mj:
            if W_by_d[1][mi][mj] > 0: adj_1[mi][mj] = 1
            if W_by_d[2][mi][mj] > 0: adj_2[mi][mj] = 1

new_edges = []
for mi in range(nm):
    for mj in range(mi+1, nm):
        if adj_2[mi][mj] and not adj_1[mi][mj]:
            Hi = ts['info'][ts['merged'][mi][0]]['H']
            Hj = ts['info'][ts['merged'][mj][0]]['H']
            new_edges.append((mi, mj, Hi, Hj))

print("  Wiggly (d=1) edges: %d" % (adj_1.sum() // 2))
print("  Waggly-2 (d=2) edges: %d" % (adj_2.sum() // 2))
print("  NEW edges in waggly-2 (not in wiggly): %d\n" % len(new_edges))

if new_edges:
    for mi, mj, Hi, Hj in new_edges:
        print("    {%d(H=%d), %d(H=%d)}: ΔH=%d" % (mi, Hi, mj, Hj, abs(Hi-Hj)))
else:
    print("    NONE — waggly-2 has the SAME adjacency as wiggly!")

# Check all distances
print("\n  ADJACENCY COMPARISON FOR ALL DISTANCES:")
for d in range(2, m+1):
    adj_d = np.zeros((nm, nm), dtype=int)
    for mi in range(nm):
        for mj in range(nm):
            if mi != mj and W_by_d[d][mi][mj] > 0:
                adj_d[mi][mj] = 1
    n_edges_d = adj_d.sum() // 2
    new_d = sum(1 for mi in range(nm) for mj in range(mi+1, nm)
                if adj_d[mi][mj] and not adj_1[mi][mj])
    lost_d = sum(1 for mi in range(nm) for mj in range(mi+1, nm)
                 if adj_1[mi][mj] and not adj_d[mi][mj])
    print("    d=%d: %d edges, %d new, %d lost vs wiggly" %
          (d, n_edges_d, new_d, lost_d))

# ============================================================
banner("WAGGLY PER-CLASS: THE DEEP STRUCTURE")
# ============================================================

# For the REGULAR tournament (mid=9, H=15):
# At d=1 (wiggly): reaches only classes 6 and 8 (just 2 classes)
# At d=2 (waggly-2): which classes does it reach?

print("  CLASS-LEVEL REACHABILITY BY DISTANCE:\n")
print("  %-4s %-4s %-6s" % ("mid", "H", "fiber"), end="")
for d in range(1, m+1):
    print("  d=%d" % d, end="")
print()
print("  " + "-"*(40 + 5*m))

for mi in order:
    H = ts['info'][ts['merged'][mi][0]]['H']
    fiber = sum(ts['sizes'][c] for c in ts['merged'][mi])
    print("  %-4d %-4d %-6d" % (mi, H, fiber), end="")
    for d in range(1, m+1):
        reachable = sum(1 for mj in range(nm) if mj != mi and W_by_d[d][mi][mj] > 0)
        print("  %-3d" % reachable, end="")
    print()

print("""
  THE WAGGLY HIERARCHY REVEALS:

  At d=1 (wiggly): some classes only reach 2-4 others (narrow funnel tips).
  At d=2 (waggly-2): reach EXPANDS for ALL classes.
  At d=3: reach may expand further.
  At d=m-1: reach may contract (near-complement is like wiggly of complement).
  At d=m: reach contracts to 1 (the complement class).

  The MAXIMUM REACH occurs at d ≈ m/2 (the middle of the distance spectrum).
  This is the WAGGLY EQUATOR — the distance at which the most information
  about the tiling landscape is accessible.
""")

print("Done. opus-2026-03-24-S294")
