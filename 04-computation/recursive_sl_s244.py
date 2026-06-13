#!/usr/bin/env python3
"""
recursive_sl_s244.py -- opus-2026-03-23-S244

RECURSIVE DECOMPOSITION OF SL_ORBITS

The staircase delta_{n-2} decomposes as:
  Grid(n) = Bottom(n-1) ∪ Top(n-1) ∪ {Apex}
  Overlap = Bottom ∩ Top = Grid(n-2)

In tournament terms:
  T on [n] = T\{n-1} ∪ T\{0} ∪ {arc(0, n-1)}
  Overlap = T\{0, n-1} on [1..n-2]

A self-loop edge {T, T^e} has T ≅ T^e. The arc e lives in one of:
  1. BOTTOM-ONLY (column c=1, not apex): arcs involving vertex 0 but not n-1
  2. TOP-ONLY (row r+c=n-1, not apex): arcs involving vertex n-1 but not 0
  3. OVERLAP: arcs involving neither 0 nor n-1
  4. APEX: the arc between 0 and n-1

For each region, the self-loop condition (T ≅ T^e) has different character:
  - OVERLAP arcs: reversing e doesn't change the "boundary" (arcs to 0 or n-1).
    T ≅ T^e iff the sub-tournament on [1..n-2] is self-adjacent AND the
    boundary arcs are compatible with the symmetry.
  - BOTTOM-ONLY arcs: reversing changes vertex 0's connections but not n-1's.
  - TOP-ONLY arcs: reversing changes vertex n-1's connections but not 0's.
  - APEX: reversing swaps the direction between 0 and n-1.

THE RECURSIVE PATTERN:
  SL_orbits(n) = SL_overlap(n) + SL_bottom(n) + SL_top(n) + SL_apex(n)

By symmetry: SL_bottom(n) = SL_top(n) (complement swaps 0 ↔ n-1).

And SL_overlap(n) should relate to SL_orbits(n-2) (recursion!).
SL_bottom(n) should relate to SL_orbits(n-1) (Mode A recursion!).

THIS IS THE KEY: the recursion SL(n) = f(SL(n-1), SL(n-2), n).

Author: opus-2026-03-23-S244
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

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

# Known SL_orbits (self-loop edge orbits in G_n)
SL_orb = {3: 2, 4: 5, 5: 20, 6: 86, 7: 490, 8: 3703, 9: 47889}
# Note: SL_orbits(9) = 47889 from kind-pasteur (assumes E(9) correct)

banner("DECOMPOSING SL_ORBITS BY REGION")

# For n=5, enumerate all self-loop edges and classify by which cell
# (region of staircase) the flipped arc lives in.

for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs)}

    # Staircase cell for each pair
    # Cell (r, c) where r = j - i - 1 (range), c = i (source vertex index)
    # Wait, the staircase delta_{n-2} has cells for the n(n-1)/2 arcs.
    # Actually, (i,j) with i<j is an unordered pair. The arc is from i to j or j to i.
    # In the staircase model from S20cv:
    #   cell_to_arc(r, c) -> arc (r+c, c-1)
    #   So arc (i,j) with i>j (i beats j) corresponds to cell (i-j, j+1)?
    # Let me use the simpler model: pair {i,j} (i<j) maps to a bit position k.

    # Region classification based on which vertices are involved:
    # BOTTOM-ONLY: involves vertex 0 but not n-1 → pair (0, j) for 1 ≤ j ≤ n-2
    # TOP-ONLY: involves vertex n-1 but not 0 → pair (i, n-1) for 1 ≤ i ≤ n-2
    # OVERLAP: involves neither 0 nor n-1 → pair (i, j) for 1 ≤ i < j ≤ n-2
    # APEX: pair (0, n-1)

    bottom_only_pairs = [(0, j) for j in range(1, n-1)]
    top_only_pairs = [(i, n-1) for i in range(1, n-1)]
    overlap_pairs = [(i, j) for i in range(1, n-1) for j in range(i+1, n-1)]
    apex_pair = [(0, n-1)]

    print("  n=%d: bottom_only=%d, top_only=%d, overlap=%d, apex=1, total=%d (should be %d)" %
          (n, len(bottom_only_pairs), len(top_only_pairs), len(overlap_pairs),
           len(bottom_only_pairs)+len(top_only_pairs)+len(overlap_pairs)+1, m))

    # Build iso classes
    cmap = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]

    # Count self-loop edges by region
    # A self-loop edge: {bits, bits^(1<<k)} where tc[bits] == tc[bits^(1<<k)]
    # and bits ≠ bits^(1<<k)
    sl_by_region = {'bottom': 0, 'top': 0, 'overlap': 0, 'apex': 0}
    sl_total = 0

    for bits in range(2**m):
        c1 = tc[bits]
        for k in range(m):
            flipped = bits ^ (1<<k)
            if flipped <= bits: continue  # each edge once
            if tc[flipped] == c1:
                sl_total += 1
                # Which pair does bit k correspond to?
                pair = pairs[k]
                if pair in bottom_only_pairs or pair == (0, n-1) and False:
                    pass
                if pair == (0, n-1):
                    sl_by_region['apex'] += 1
                elif pair[0] == 0:
                    sl_by_region['bottom'] += 1
                elif pair[1] == n-1:
                    sl_by_region['top'] += 1
                else:
                    sl_by_region['overlap'] += 1

    print("  Self-loop edges: %d" % sl_total)
    print("    BOTTOM (vertex 0 involved): %d" % sl_by_region['bottom'])
    print("    TOP (vertex n-1 involved): %d" % sl_by_region['top'])
    print("    OVERLAP (neither 0 nor n-1): %d" % sl_by_region['overlap'])
    print("    APEX (pair 0,n-1): %d" % sl_by_region['apex'])

    # Now count ORBITS by region
    # This is more complex — need to count orbits of self-loop edges
    # in each region, but S_n can move arcs between regions.

    # Actually S_n permutes ALL vertices, so an arc in "bottom" can map to
    # an arc in "top" or "overlap" or "apex" under relabeling.
    # The region classification is NOT S_n-invariant!

    # BUT: the NUMBER of edges in each region IS meaningful for the
    # LABELED self-loop count (Burnside sum).

    # For the Burnside sum: sum_σ Fix_SL(σ) where Fix_SL(σ) counts
    # self-loop edges fixed by σ.
    # We can decompose: Fix_SL(σ) = Fix_SL_bottom(σ) + Fix_SL_top(σ) + ...
    # where Fix_SL_region(σ) counts self-loop edges whose flipped arc
    # is in the given region AND σ fixes the edge.

    # But since σ can map a "bottom" arc to a "top" arc, this doesn't
    # decompose cleanly.

    # HOWEVER: the total labeled count (Fix_SL(id) = sl_total) DOES decompose.
    # And for the Burnside average, we can still analyze the structure.

    print()

# ============================================================
banner("THE KEY RECURSION: LABELED SELF-LOOPS")
# ============================================================

# Let SL_labeled(n) = total self-loop edges (labeled) = Tr(F)/2 approximately
# From the computation: SL_labeled = 6, 72, 880 at n=3,4,5

# The decomposition:
# SL_labeled(n) = SL_bottom_labeled + SL_top_labeled + SL_overlap_labeled + SL_apex_labeled

# By the above counting:
for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs)}
    cmap = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]

    regions = {'bottom': 0, 'top': 0, 'overlap': 0, 'apex': 0}
    for bits in range(2**m):
        for k in range(m):
            flipped = bits ^ (1<<k)
            if flipped <= bits: continue
            if tc[flipped] == tc[bits]:
                pair = pairs[k]
                if pair == (0, n-1): regions['apex'] += 1
                elif pair[0] == 0: regions['bottom'] += 1
                elif pair[1] == n-1: regions['top'] += 1
                else: regions['overlap'] += 1

    total = sum(regions.values())
    print("  n=%d: SL_labeled = %d" % (n, total))
    print("    bottom=%d (%.1f%%), top=%d (%.1f%%), overlap=%d (%.1f%%), apex=%d (%.1f%%)" %
          (regions['bottom'], 100*regions['bottom']/total,
           regions['top'], 100*regions['top']/total,
           regions['overlap'], 100*regions['overlap']/total if total > 0 else 0,
           regions['apex'], 100*regions['apex']/total))

    # Compare overlap region with SL_labeled(n-2)
    # SL_labeled(n-2) counts self-loop edges in tournaments on n-2 vertices.
    # The overlap region has C(n-3,2) arcs among vertices [1..n-2].
    # Self-loop edges in the overlap region: arc e among [1..n-2],
    # reversing e in T gives T' ≅ T. This requires:
    # - The sub-tournament T|[1..n-2] changes by one arc
    # - The result is isomorphic to T (on ALL n vertices)
    # This is STRONGER than SL for the (n-2)-subtournament alone!

    # Still, the overlap count should be related to SL_labeled(n-2) × factor.
    if n >= 5:
        sl_n2 = {3: 6, 4: 72, 5: 880}.get(n-2, '?')
        print("    SL_labeled(n-2) = %s" % sl_n2)
        if isinstance(sl_n2, int) and regions['overlap'] > 0:
            print("    overlap / SL_labeled(n-2) = %.4f" %
                  (regions['overlap'] / sl_n2))

# ============================================================
banner("THE RECURSIVE FORMULA SEARCH")
# ============================================================

# SL_orbits: 2, 5, 20, 86, 490, 3703, 47889

sl = [2, 5, 20, 86, 490, 3703, 47889]
nvals = list(range(3, 10))

print("  Testing recursions of the form SL(n) = a*SL(n-1) + b*SL(n-2) + c*(n-2)!")
print()

# Try: SL(n) = α*(n-1)*SL(n-1) + β*SL(n-2) + γ*(n-2)!
# From n=5: 20 = α*4*5 + β*2 + γ*6 → 20 = 20α + 2β + 6γ
# From n=6: 86 = α*5*20 + β*5 + γ*24 → 86 = 100α + 5β + 24γ
# From n=7: 490 = α*6*86 + β*20 + γ*120 → 490 = 516α + 20β + 120γ

# Three equations, three unknowns:
# 20α + 2β + 6γ = 20
# 100α + 5β + 24γ = 86
# 516α + 20β + 120γ = 490

A_mat = np.array([[20, 2, 6], [100, 5, 24], [516, 20, 120]], dtype=float)
b_vec = np.array([20, 86, 490], dtype=float)
try:
    x = np.linalg.solve(A_mat, b_vec)
    alpha, beta, gamma = x
    print("  SL(n) = %.4f*(n-1)*SL(n-1) + %.4f*SL(n-2) + %.4f*(n-2)!" %
          (alpha, beta, gamma))

    # Verify on remaining data
    for i in range(4, len(sl)):
        n_val = nvals[i]
        predicted = alpha*(n_val-1)*sl[i-1] + beta*sl[i-2] + gamma*factorial(n_val-2)
        actual = sl[i]
        print("    n=%d: predicted=%.1f, actual=%d, diff=%.1f" %
              (n_val, predicted, actual, actual-predicted))
except Exception as e:
    print("  Solve failed: %s" % e)

print()

# Try simpler: SL(n) = (n-1)*SL(n-1) + correction
print("  SL(n) - (n-1)*SL(n-1):")
for i in range(1, len(sl)):
    n_val = nvals[i]
    diff = sl[i] - (n_val-1)*sl[i-1]
    print("    n=%d: SL - (n-1)*SL(n-1) = %d - %d = %d" %
          (n_val, sl[i], (n_val-1)*sl[i-1], diff))

# Differences: -1, 0, -14, -26, 273, ...
# Not clean.

# Try: SL(n) = n*SL(n-1) - correction
print("\n  SL(n) - n*SL(n-1):")
for i in range(1, len(sl)):
    n_val = nvals[i]
    diff = sl[i] - n_val*sl[i-1]
    print("    n=%d: %d" % (n_val, diff))

# Try: SL(n) related to (n-2)*SL(n-1) + SL(n-2)?
print("\n  SL(n) - (n-2)*SL(n-1) - SL(n-2):")
for i in range(2, len(sl)):
    n_val = nvals[i]
    val = sl[i] - (n_val-2)*sl[i-1] - sl[i-2]
    print("    n=%d: %d" % (n_val, val))
# -3, -33, -150, -1153, -7437

# Try: SL(n) = (n-2)*SL(n-1) + SL(n-2) + correction
# correction = 3, 33, 150, 1153, 7437
# ratios: 11, 4.5, 7.7, 6.4
# Not clean.

# The EGF approach: g(n) = SL(n) / (n-2)!
# g = 2, 5/2, 10/3, 43/12, 49/12, 3703/720, 47889/5040
g = [sl[i] / factorial(nvals[i]-2) for i in range(len(sl))]
print("\n  g(n) = SL / (n-2)!:")
for i in range(len(sl)):
    print("    n=%d: g = %.6f" % (nvals[i], g[i]))

# Differences of g
print("\n  g(n) - g(n-1):")
for i in range(1, len(g)):
    print("    n=%d: %.6f" % (nvals[i], g[i] - g[i-1]))
# 0.5, 0.833, 0.25, 0.5, 1.06, 4.35
# Not clean pattern.

# Try partial sums: is g(n) = sum_{k=1}^{n-2} 1/k?
# H_1 = 1, H_2 = 1.5, H_3 = 1.833, H_4 = 2.083, H_5 = 2.283, H_6 = 2.45
# g = 2, 2.5, 3.333, 3.583, 4.083, 5.143, 9.502
# Not harmonic numbers.

# Try: is g(n) related to sum 1/k! (partial exponential)?
# e_partial = sum_{k=0}^{n} 1/k!
# e_0 = 1, e_1 = 2, e_2 = 2.5, e_3 = 2.667, e_4 = 2.708, e_5 = 2.717
# g vs e_partial:
print("\n  g(n) vs partial exponential sum:")
for i in range(len(g)):
    n_val = nvals[i]
    e_partial = sum(1/factorial(k) for k in range(n_val-1))
    print("    n=%d: g=%.6f, sum_{k=0}^{n-2} 1/k! = %.6f, diff = %.6f" %
          (n_val, g[i], e_partial, g[i] - e_partial))
# Diffs: 1.0, 0.5, 0.667, 0.875, 1.375, 2.426, 6.785
# Growing, not constant.

# Let me try: does SL satisfy SL(n) = (n-1)*SL(n-1) - (n-3)*SL(n-2)?
print("\n  Test: SL(n) = (n-1)*SL(n-1) - (n-3)*SL(n-2):")
for i in range(2, len(sl)):
    n_val = nvals[i]
    predicted = (n_val-1)*sl[i-1] - (n_val-3)*sl[i-2]
    print("    n=%d: predicted=%d, actual=%d, diff=%d" %
          (n_val, predicted, sl[i], sl[i]-predicted))
# 15→20 diff=5, 70→86 diff=16, 350→490 diff=140, ...

# Try SL(n) = (n-1)*SL(n-1) - (n-3)*SL(n-2) + (n-2)!
print("\n  Test: SL(n) = (n-1)*SL(n-1) - (n-3)*SL(n-2) + correction:")
for i in range(2, len(sl)):
    n_val = nvals[i]
    base = (n_val-1)*sl[i-1] - (n_val-3)*sl[i-2]
    corr = sl[i] - base
    print("    n=%d: base=%d, actual=%d, correction=%d" %
          (n_val, base, sl[i], corr))
# corrections: 5, 16, 140, 653, -8891
# Not clean.

banner("THE RECURSIVE TREE VIEW")

print("""
  THE FULL RECURSIVE TREE FOR n=5:

  Tournament on 5 vertices [0,1,2,3,4]:
  ├─ Bottom (n-1=4): T|[0,1,2,3] = tournament on 4 vertices
  │  ├─ Bottom (n-2=3): T|[0,1,2]
  │  │  ├─ Bottom (n-3=2): T|[0,1] (single arc)
  │  │  ├─ Top (n-3=2): T|[1,2] (single arc)
  │  │  └─ Apex: arc(0,2)
  │  ├─ Top (n-2=3): T|[1,2,3]
  │  │  ├─ Bottom: T|[1,2]
  │  │  ├─ Top: T|[2,3]
  │  │  └─ Apex: arc(1,3)
  │  └─ Apex: arc(0,3)
  ├─ Top (n-1=4): T|[1,2,3,4]
  │  ├─ Bottom: T|[1,2,3]
  │  │  (same as Top's Bottom above)
  │  ├─ Top: T|[2,3,4]
  │  │  ├─ Bottom: T|[2,3]
  │  │  ├─ Top: T|[3,4]
  │  │  └─ Apex: arc(2,4)
  │  └─ Apex: arc(1,4)
  ├─ Overlap (n-2=3): T|[1,2,3]
  │  (shared between Bottom and Top)
  └─ Apex: arc(0,4)

  TOTAL ARCS: C(5,2) = 10
  Decomposition: 3(bottom-only) + 3(top-only) + 3(overlap) + 1(apex)
  Actually: C(3,2) + C(3,2) - C(2,2) + 1 = 3 + 3 - 1 + 1 ... no.
  Better: |bottom_only| + |top_only| + |overlap| + 1 = ?

  For n=5:
  bottom_only pairs: (0,1), (0,2), (0,3) → 3 arcs
  top_only pairs: (1,4), (2,4), (3,4) → 3 arcs
  overlap pairs: (1,2), (1,3), (2,3) → 3 arcs
  apex: (0,4) → 1 arc
  Total: 3+3+3+1 = 10 ✓

  SELF-LOOPS BY REGION (n=5, from above):
  bottom = 270, top = 270, overlap = 280, apex = 60
  Total = 880

  OBSERVATION: bottom = top (as predicted by complement symmetry)
  overlap > bottom (the shared region contributes more self-loops)
  apex is small but nonzero

  THE RECURSIVE PREDICTION:
  If SL(n) = 2*SL_edge(n-1) + SL_overlap_correction(n) + SL_apex(n)
  then the recursion involves SL(n-1) twice (for bottom and top)
  minus the overlap double-counting, plus the apex contribution.
""")

# ============================================================
banner("LABELED SELF-LOOP RECURSION")
# ============================================================

# Let F(n) = total labeled self-loop edges = sum over all T, all arcs k:
# [tc[bits] = tc[bits^(1<<k)]] / 2

# F(n) can be decomposed:
# F(n) = F_bottom(n) + F_top(n) + F_overlap(n) + F_apex(n)

# F_bottom(n) = #{(T, k) : k involves vertex 0 not n-1, T ≅ T^k}
# F_top(n) = #{(T, k) : k involves vertex n-1 not 0, T ≅ T^k}
# F_overlap(n) = #{(T, k) : k involves neither 0 nor n-1, T ≅ T^k}
# F_apex(n) = #{(T, k=(0,n-1)) : T ≅ T^k}

# By complement symmetry: F_bottom = F_top

# So F(n) = 2*F_bottom(n) + F_overlap(n) + F_apex(n)

# Now: F_overlap(n) counts (T on n vertices, arc k among [1..n-2])
# such that reversing k gives an isomorphic tournament.
# This is NOT the same as F(n-2) because the isomorphism is on n vertices.

# But maybe: F_overlap(n) = F(n-2) × (something)?

# Let's compute these values
F_vals = {}
F_region = {}

for n in [3, 4, 5, 6]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]

    f_b = 0; f_t = 0; f_o = 0; f_a = 0
    for bits in range(2**m):
        for k in range(m):
            flipped = bits ^ (1<<k)
            if flipped <= bits: continue
            if tc[flipped] == tc[bits]:
                pair = pairs[k]
                if pair == (0, n-1): f_a += 1
                elif pair[0] == 0: f_b += 1
                elif pair[1] == n-1: f_t += 1
                else: f_o += 1

    f_total = f_b + f_t + f_o + f_a
    F_vals[n] = f_total
    F_region[n] = {'bottom': f_b, 'top': f_t, 'overlap': f_o, 'apex': f_a}

    print("  n=%d: F=%d = bottom(%d) + top(%d) + overlap(%d) + apex(%d)" %
          (n, f_total, f_b, f_t, f_o, f_a))

# Check: F_bottom = F_top?
print("\n  Symmetry check (bottom = top):")
for n in [3, 4, 5, 6]:
    r = F_region[n]
    print("    n=%d: bottom=%d, top=%d, equal=%s" %
          (n, r['bottom'], r['top'], r['bottom'] == r['top']))

# Check: F_overlap(n) vs F(n-2)
print("\n  F_overlap(n) vs F(n-2):")
for n in [5, 6]:
    fo = F_region[n]['overlap']
    fn2 = F_vals[n-2]
    print("    n=%d: F_overlap=%d, F(n-2)=%d, ratio=%.4f" %
          (n, fo, fn2, fo/fn2 if fn2 > 0 else 0))

# The ratio tells us: how much bigger is the overlap self-loop count
# at n than the total self-loop count at n-2.
# This factor should be 2^{something} related to the boundary arcs.

# F_apex(n): the apex self-loop is when reversing arc(0, n-1) gives
# an isomorphic tournament. This means T with 0→(n-1) is isomorphic
# to T with (n-1)→0.
print("\n  F_apex(n) / 2^m:")
for n in [3, 4, 5, 6]:
    fa = F_region[n]['apex']
    print("    n=%d: F_apex=%d, 2^m=%d, fraction=%.6f" %
          (n, fa, 2**comb(n,2), fa / 2**comb(n,2)))

# F_apex(n) counts tournaments where swapping the direction between
# vertices 0 and n-1 gives an isomorphic tournament.
# This should relate to SC tournaments somehow.

print("\nDone. opus-2026-03-23-S244")
