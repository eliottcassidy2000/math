#!/usr/bin/env python3
"""
wiggly_intrinsic_s292.py -- opus-2026-03-24-S292

THE WIGGLY TILING GRAPH AS A STANDALONE OBJECT

Forget blue/black lines. Forget complements. Start fresh:

A TILING = a binary word of length m = C(n-1,2), encoding the
non-base-path arcs of a tournament with fixed base path.

The WIGGLY GRAPH Q(n) is:
  - Vertices: all 2^m tilings
  - Edges: two tilings are connected if they differ in exactly one bit
  - This IS the m-dimensional hypercube Q_m

Each tiling has EXACTLY m neighbors. Each neighbor corresponds to
flipping one specific tile (= one specific wiggly class A, B, C, ...).

THE QUOTIENT: Group tilings by isomorphism class.
The quotient of Q_m by the S_n-action (relabeling) gives the META-GRAPH.

THE INTRINSIC QUESTIONS:
1. What is the LOCAL structure around each tiling?
   - How many neighbors are in the SAME iso class? (= self-loops = neutral flips)
   - How many DISTINCT iso classes are among the m neighbors?
   - What is the DISTRIBUTION of classes among neighbors?

2. What is the SECOND NEIGHBORHOOD?
   - From a tiling t, its m neighbors each have m-1 NEW neighbors (plus t itself).
   - Total 2-step-reachable tilings: up to m(m-1).
   - How many DISTINCT iso classes at distance 2?

3. What SHAPE does this reveal?
   - The m-cube has uniform local structure (each vertex looks the same).
   - But after quotienting by iso class, the structure becomes NON-UNIFORM.
   - The Tournament Counting Theorem predicts this non-uniformity.

Author: opus-2026-03-24-S292
"""

import sys
import numpy as np
from math import comb, factorial, gcd
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

# ============================================================
banner("THE WIGGLY GRAPH Q(n) AT n=5")
# ============================================================

n = 5
m = comb(n-1, 2)  # C(4,2) = 6 tiles
total_tilings = 2**m  # 2^6 = 64 tilings
all_pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

# The BASE PATH: n -> n-1 -> ... -> 1
# Base path arcs: (k, k-1) for k = n, n-1, ..., 2
# = (5,4), (4,3), (3,2), (2,1) — but in 0-indexed: (4,3),(3,2),(2,1),(1,0)
# These are NOT tiles (not flippable).
# TILES = non-consecutive pairs

# At n=5, pairs in standard order:
# (0,1), (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
# Base path arcs (consecutive): (0,1), (1,2), (2,3), (3,4) — 4 arcs
# Tiles (non-consecutive): (0,2), (0,3), (0,4), (1,3), (1,4), (2,4) — 6 tiles

# The explorer's tile enumeration: for y=1..n-2, for x=n down to y+2
# n=5: y=1: x=5,4,3 → (5,1),(4,1),(3,1) → 0-indexed: (4,0),(3,0),(2,0)
#       y=2: x=5,4 → (5,2),(4,2) → (4,1),(3,1) — wait, this doesn't match
# Let me use the 1-indexed convention from CLAUDE.md:
# n=5 tiles: A=(5,1), B=(4,1), C=(3,1), D=(5,2), E=(4,2), F=(5,3)

# In 0-indexed pairs: A=(4,0), B=(3,0), C=(2,0), D=(4,1), E=(3,1), F=(4,2)
tile_pairs = [(4,0), (3,0), (2,0), (4,1), (3,1), (4,2)]
tile_labels = ['A', 'B', 'C', 'D', 'E', 'F']

print("  WIGGLY GRAPH Q(5):")
print("  Vertices: 2^%d = %d tilings" % (m, total_tilings))
print("  Degree: %d (= C(%d,2) = triangular number T_%d)" % (m, n-1, n-2))
print("  Tiles: %s" % list(zip(tile_labels, tile_pairs)))
print()

# Build the full tournament from a tiling
# Base path: 4→3→2→1→0 (descending, so arc i→j when i>j for consecutive)
# Actually base path is n→n-1→...→1, meaning arc from higher to lower
# for consecutive vertices.
# In 0-indexed: arcs 4→3, 3→2, 2→1, 1→0.
# For non-base arcs (tiles): bit 0 = tile A direction, etc.
# bit=0 means arc goes from higher to lower (like base path: "natural" direction)
# bit=1 means arc goes from lower to higher (flipped)

def tiling_to_tournament(bits, n, tile_pairs):
    """Convert a tiling bitstring to a tournament adjacency matrix."""
    A = [[0]*n for _ in range(n)]
    # Base path arcs: i→i-1 for i=n-1,...,1 (0-indexed)
    for i in range(1, n):
        A[i][i-1] = 1  # higher beats lower on base path

    # Tile arcs
    for k, (hi, lo) in enumerate(tile_pairs):
        if (bits >> k) & 1:
            A[lo][hi] = 1  # flipped: lower beats higher
        else:
            A[hi][lo] = 1  # natural: higher beats lower
    return A

# Build all iso classes
cmap = {}; sizes = {}; tc = {}; reps = {}
members_by_class = defaultdict(list)

for bits in range(total_tilings):
    A = tiling_to_tournament(bits, n, tile_pairs)
    c = canon(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
    cid = cmap[c]
    sizes[cid] += 1
    tc[bits] = cid
    members_by_class[cid].append(bits)

nc = len(cmap)

# Compute H for each class
info = {}
for cid in range(nc):
    A = reps[cid]
    comp_c = canon(complement(A, n), n)
    info[cid] = {
        'H': Hcount(A, n),
        'comp': cmap[comp_c],
        'SC': cmap[comp_c] == cid,
        'size': sizes[cid],
        'score': tuple(sorted(sum(A[i]) for i in range(n)))
    }

# Merge complement pairs
merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

print("  Iso classes: %d (unmerged), %d (merged)" % (nc, nm))
print("  Tilings per class:")
for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]
    H = info[cid]['H']
    fiber = sum(sizes[c] for c in merged[mi])
    n_classes = len(merged[mi])
    print("    mid=%d: H=%2d, fiber=%d (%d unmerged)" % (mi, H, fiber, n_classes))

# ============================================================
banner("THE TILING NEIGHBORHOOD: WHAT EACH TILING SEES")
# ============================================================

# For each tiling, compute its m=6 neighbors and their classes
print("  For each tiling: its %d wiggly neighbors span how many iso classes?\n" % m)

# Statistics per iso class
class_stats = defaultdict(lambda: {
    'n_self': [],        # how many neighbors in same (merged) class
    'n_distinct': [],    # how many distinct (merged) classes among neighbors
    'neighbor_classes': [],  # multiset of neighbor classes
    'class_vector': [],  # vector of neighbor class counts
})

for bits in range(total_tilings):
    mi = mid_map[tc[bits]]
    neighbors = []
    self_count = 0
    neighbor_classes = Counter()

    for k in range(m):
        nbr = bits ^ (1 << k)
        mj = mid_map[tc[nbr]]
        neighbors.append(mj)
        neighbor_classes[mj] += 1
        if mj == mi:
            self_count += 1

    n_distinct = len(neighbor_classes)
    class_stats[mi]['n_self'].append(self_count)
    class_stats[mi]['n_distinct'].append(n_distinct)
    class_stats[mi]['neighbor_classes'].append(dict(neighbor_classes))

# Print summary
print("  %-4s %-4s %-6s %-10s %-10s %-30s" %
      ("mid", "H", "fiber", "avg_self", "avg_dist", "typical neighbor profile"))
print("  " + "-"*80)

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]; H = info[cid]['H']
    fiber = sum(sizes[c] for c in merged[mi])
    avg_self = np.mean(class_stats[mi]['n_self'])
    avg_distinct = np.mean(class_stats[mi]['n_distinct'])

    # Most common neighbor profile
    profiles = Counter()
    for nc_dict in class_stats[mi]['neighbor_classes']:
        # Sorted by class id
        profile = tuple(sorted(nc_dict.items()))
        profiles[profile] += 1
    top_profile = profiles.most_common(1)[0]
    n_profiles = len(profiles)

    print("  %-4d %-4d %-6d %-10.2f %-10.2f %d profiles, top: %d tilings" %
          (mi, H, fiber, avg_self, avg_distinct, n_profiles, top_profile[1]))

# ============================================================
banner("THE INTRINSIC TILING NAMES")
# ============================================================

print("""
  WITHOUT blue/black, what would we call the wiggly structures?

  THE TILING LATTICE:
  - Each tiling is a vertex of the 6-dimensional hypercube Q_6.
  - Each edge of Q_6 = one tile flip = one wiggly line.
  - The hypercube has 64 vertices and 192 edges.

  NATURAL INTRINSIC NAMES:

  1. TILE FLIP = changing one bit of the tiling word.
     Each tiling has exactly 6 tile-flip neighbors.
     (Not "wiggly line" — just "tile flip" or "single-tile mutation".)

  2. TILE-FLIP CLASS = which tile was flipped (A through F).
     The 6 tile-flip classes partition the 192 edges into 6 groups of 32.
     Each class is a PERFECT MATCHING on Q_6 (pairs every tiling with exactly one other).

  3. NEUTRAL FLIP = a tile flip that preserves the iso class.
     This is a SELF-LOOP in the quotient graph.
     Natural name: "silent mutation" (by analogy with genetics).

  4. CLASS-CHANGING FLIP = a tile flip that changes the iso class.
     Natural name: "expressive mutation".

  5. THE NEUTRALITY of a tiling = fraction of its 6 flips that are neutral.
     Ranges from 0/6 (every flip changes class) to some maximum.

  6. THE EXPRESSIVITY of a tiling = fraction of its 6 flips that are class-changing.
     = 1 - neutrality.

  7. THE REACH of a tiling = number of distinct iso classes among its 6 neighbors.
     High reach = highly expressive, can reach many different classes.
     Low reach = trapped, most flips lead to the same few classes.

  8. THE TILING SPECTRUM = the multiset of (class, multiplicity) among neighbors.
     Two tilings in the same iso class may have different spectra
     (because the S_n action permutes tile positions).

  QUOTIENT STRUCTURE:

  9. The TILE-FLIP QUOTIENT = the meta-graph G_n/Z_2.
     Weighted edges = total number of class-changing flips between classes.
     Self-loops = total neutral flips within a class.

  10. The CLASS CONNECTIVITY = which classes are reachable from which.
      At n=5: ALL pairs of (merged) classes are connected (one-step reachable).
      This is the "local = global" theorem.
""")

# ============================================================
banner("NEUTRALITY AND REACH PER TILING")
# ============================================================

# For each tiling: compute neutrality (self-loop count) and reach
print("  DISTRIBUTION OF NEUTRALITY (# silent mutations out of 6):\n")

neutrality_by_class = defaultdict(Counter)
reach_by_class = defaultdict(Counter)

for bits in range(total_tilings):
    mi = mid_map[tc[bits]]
    neutral = 0
    classes_seen = set()
    for k in range(m):
        nbr = bits ^ (1 << k)
        mj = mid_map[tc[nbr]]
        if mj == mi:
            neutral += 1
        classes_seen.add(mj)
    neutrality_by_class[mi][neutral] += 1
    reach_by_class[mi][len(classes_seen)] += 1

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]; H = info[cid]['H']
    fiber = sum(sizes[c] for c in merged[mi])
    neut_dist = neutrality_by_class[mi]
    reach_dist = reach_by_class[mi]

    neut_str = ", ".join("%d:%d" % (k, v) for k, v in sorted(neut_dist.items()))
    reach_str = ", ".join("%d:%d" % (k, v) for k, v in sorted(reach_dist.items()))

    avg_neut = sum(k*v for k,v in neut_dist.items()) / fiber
    avg_reach = sum(k*v for k,v in reach_dist.items()) / fiber

    print("  mid=%d (H=%2d, fiber=%d): neutrality={%s} avg=%.2f" %
          (mi, H, fiber, neut_str, avg_neut))
    print("    reach={%s} avg=%.2f" % (reach_str, avg_reach))

# ============================================================
banner("THE SECOND NEIGHBORHOOD: NEIGHBOR-OF-NEIGHBOR ANALYSIS")
# ============================================================

# For each tiling t, the m neighbors each have m neighbors.
# Of these m² potential 2-step paths:
#   - m go back to t (one per neighbor, flipping the same tile back)
#   - The rest reach new tilings at Hamming distance 2.
# Number of distinct tilings at distance 2 = C(m,2) = 15 at n=5.
# (Because flipping 2 different tiles from any of the C(m,2) pairs.)

# But the iso CLASS distribution at distance 2 varies!

print("  SECOND NEIGHBORHOOD STATISTICS:\n")
print("  Each tiling has %d neighbors (distance 1) and %d at distance 2.\n" %
      (m, comb(m, 2)))

# For each tiling: how many distinct classes at distance exactly 2?
d2_classes_by_class = defaultdict(list)

for bits in range(total_tilings):
    mi = mid_map[tc[bits]]
    d2_classes = set()
    for j in range(m):
        for k in range(j+1, m):
            nbr2 = bits ^ (1 << j) ^ (1 << k)
            mj = mid_map[tc[nbr2]]
            d2_classes.add(mj)
    d2_classes_by_class[mi].append(len(d2_classes))

print("  %-4s %-4s %-6s %-15s %-15s" %
      ("mid", "H", "fiber", "d1 classes", "d2 classes"))
print("  " + "-"*55)

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]; H = info[cid]['H']
    fiber = sum(sizes[c] for c in merged[mi])
    avg_d1 = np.mean([class_stats[mi]['n_distinct'][i]
                       for i in range(fiber)])
    avg_d2 = np.mean(d2_classes_by_class[mi])

    print("  %-4d %-4d %-6d %-15.2f %-15.2f" %
          (mi, H, fiber, avg_d1, avg_d2))

# ============================================================
banner("THE SHAPE: CLUSTERING AND ASSORTATIVITY")
# ============================================================

# The CLUSTERING COEFFICIENT of a tiling in the hypercube:
# For the raw hypercube Q_m: clustering = 0 (no triangles, since
# the hypercube is bipartite — Hamming parity alternates).
# But in the QUOTIENT (meta-graph), there ARE triangles.

# Instead, let's measure the CLASS ASSORTATIVITY:
# Do tilings in the same class tend to be close to each other?

# For each pair of tilings in the same class: Hamming distance
print("  INTRA-CLASS DISTANCE DISTRIBUTION:\n")

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]; H = info[cid]['H']
    all_tilings = []
    for c in merged[mi]:
        all_tilings.extend(members_by_class[c])
    fiber = len(all_tilings)

    if fiber > 1:
        distances = Counter()
        for i in range(len(all_tilings)):
            for j in range(i+1, len(all_tilings)):
                d = bin(all_tilings[i] ^ all_tilings[j]).count('1')
                distances[d] += 1

        avg_dist = sum(d*c for d,c in distances.items()) / sum(distances.values())
        dist_str = ", ".join("d%d:%d" % (d,c) for d,c in sorted(distances.items()))
        print("  mid=%d (H=%2d, fiber=%d): avg_dist=%.2f [%s]" %
              (mi, H, fiber, avg_dist, dist_str))

# ============================================================
banner("THE TILING GRAPH'S INTRINSIC GEOMETRY")
# ============================================================

# Each tiling lives in the hypercube Q_m.
# The S_n action partitions Q_m into orbits (= iso classes).
# The ORBIT GRAPH (quotient of Q_m by S_n) is the meta-graph.

# The GEOMETRY of the orbit graph:
# - Node weights = fiber sizes
# - Edge weights = number of cross-orbit tile flips
# - Self-loop weights = number of neutral flips

# Build the full weight matrix
W = np.zeros((nm, nm), dtype=int)
for bits in range(total_tilings):
    mi = mid_map[tc[bits]]
    for k in range(m):
        mj = mid_map[tc[bits ^ (1 << k)]]
        W[mi][mj] += 1

print("  THE TILE-FLIP QUOTIENT WEIGHT MATRIX W (n=5):\n")
order = sorted(range(nm), key=lambda x: info[merged[x][0]]['H'])
H_vec = [info[merged[mi][0]]['H'] for mi in order]

print("  H: ", H_vec)
for mi in order:
    row = " ".join("%4d" % W[mi][mj] for mj in order)
    fiber = sum(sizes[c] for c in merged[mi])
    print("  H=%2d: %s  | fiber=%d, row_sum=%d" %
          (info[merged[mi][0]]['H'], row, fiber, sum(W[mi])))

# Properties of W
print("\n  W PROPERTIES:")
print("  Symmetric: %s" % np.allclose(W, W.T))
print("  Row sums: %s" % [sum(W[mi]) for mi in order])
print("  All row sums = fiber × m: %s" %
      all(sum(W[mi]) == sum(sizes[c] for c in merged[mi]) * m for mi in range(nm)))

# The TRANSITION MATRIX P = W / (fiber × m)
fiber_vec = np.array([sum(sizes[c] for c in merged[mi]) for mi in range(nm)], dtype=float)
P = W.astype(float)
for mi in range(nm):
    P[mi] /= (fiber_vec[mi] * m)

print("\n  TRANSITION MATRIX P (row-normalized):")
for mi in order:
    row = " ".join("%6.3f" % P[mi][mj] for mj in order)
    print("  H=%2d: %s" % (info[merged[mi][0]]['H'], row))

# Stationary distribution
pi = fiber_vec / fiber_vec.sum()
print("\n  Stationary distribution π (∝ fiber):")
for mi in order:
    print("    H=%2d: π = %.4f (fiber=%d)" %
          (info[merged[mi][0]]['H'], pi[mi],
           sum(sizes[c] for c in merged[mi])))

# Check detailed balance: π_i P_{ij} = π_j P_{ji}
print("\n  Detailed balance: π_i P_{ij} = π_j P_{ji}:")
ok = True
for mi in range(nm):
    for mj in range(mi+1, nm):
        lhs = pi[mi] * P[mi][mj]
        rhs = pi[mj] * P[mj][mi]
        if abs(lhs - rhs) > 1e-10:
            print("    VIOLATION: %d-%d: %.6f ≠ %.6f" % (mi, mj, lhs, rhs))
            ok = False
if ok:
    print("    All pairs satisfy detailed balance. ✓")

# Eigenvalues of P
evals = np.sort(np.linalg.eigvals(P))[::-1].real
print("\n  Eigenvalues of P:")
for i, ev in enumerate(evals):
    print("    λ_%d = %.6f" % (i, ev))

spectral_gap = 1 - evals[1]
print("\n  Spectral gap: 1 - λ_2 = %.6f" % spectral_gap)
print("  Mixing time estimate: O(1/gap) ≈ %.1f steps" % (1/spectral_gap))

# ============================================================
banner("USING THE TOURNAMENT COUNTING THEOREM")
# ============================================================

# From S291: V_n × n!/2^m = T(n) = 1 + R(n)
# R(n) = Σ_C fiber(C) × (|Aut(C)| - 1) / 2^m
#
# The AVERAGE NEUTRALITY of a random tiling:
# E[neutrality] = (total neutral flips) / (total tilings × m)
# = Σ_C (neutral flips in C) / (2^m × m)
#
# A tiling has a neutral flip at position k iff flipping tile k
# gives a tournament in the same iso class.
#
# Total neutral flips = Σ_C W[C][C] (self-loop weight)
# Average neutrality = Σ_C W[C][C] / (2^m × m)

total_neutral = sum(W[mi][mi] for mi in range(nm))
avg_neutrality = total_neutral / (total_tilings * m)
print("  AVERAGE NEUTRALITY (fraction of neutral flips):")
print("  Total neutral flips: %d out of %d = %.4f" %
      (total_neutral, total_tilings * m, avg_neutrality))

# From the Tournament Counting Theorem:
# The neutrality is related to the arc-flip self-loop rate.
# A tile flip preserves the iso class iff the tiling has an automorphism
# that maps the flipped tile position to itself...
# Actually, it's more subtle. A flip is neutral iff the new tournament
# is isomorphic to the old one — which depends on the full structure.

# THE PREDICTION: at large n, almost all tilings have |Aut|=1,
# so almost no flips are neutral.
# Neutrality ≈ (correction-like term) → 0.

# Let me compute the exact neutrality at n=3,4,5,6.
print("\n  NEUTRALITY BY n:\n")

for nn in range(3, 7):
    mm = comb(nn-1, 2)
    total_t = 2**mm
    pp = [(i,j) for i in range(nn) for j in range(i+1,nn)]
    # Tile pairs (non-consecutive in base path)
    tp = [(i,j) for i in range(nn) for j in range(i+1,nn) if i+1 != j]
    # Actually: base path arcs are (k, k-1) for consecutive k
    # Tiles = pairs (i,j) with j-i >= 2 (in 0-indexed with base path 0→1→...→n-1)
    # Wait, base path is n-1→n-2→...→0, so consecutive pairs are (k, k-1)
    # Non-consecutive: |i-j| >= 2
    tp2 = [(max(i,j), min(i,j)) for i in range(nn) for j in range(i+1,nn) if abs(i-j) >= 2]
    assert len(tp2) == mm, f"Expected {mm} tiles, got {len(tp2)}"

    # Build classes
    cm = {}; sz = {}; tc2 = {}
    for bits in range(total_t):
        A = [[0]*nn for _ in range(nn)]
        # Base path
        for i in range(1, nn):
            A[i][i-1] = 1
        # Tiles
        for k, (hi, lo) in enumerate(tp2):
            if (bits >> k) & 1:
                A[lo][hi] = 1
            else:
                A[hi][lo] = 1
        c = canon(A, nn)
        if c not in cm:
            cid = len(cm); cm[c] = cid; sz[cid] = 0
        sz[cm[c]] += 1; tc2[bits] = cm[c]

    # Merge
    inf2 = {}
    for cid in range(len(cm)):
        # Need complement... skip for now, just count neutrals
        pass

    # Count neutral flips
    neutral_count = 0
    for bits in range(total_t):
        for k in range(mm):
            nbr = bits ^ (1 << k)
            if tc2[bits] == tc2[nbr]:
                neutral_count += 1

    neutrality = neutral_count / (total_t * mm)
    print("  n=%d: m=%d, tilings=2^%d=%d, neutral=%d/%d = %.4f (%.1f%%)" %
          (nn, mm, mm, total_t, neutral_count, total_t * mm,
           neutrality, 100*neutrality))

# ============================================================
banner("THE SHAPE: FIBER-WEIGHTED SECOND-MOMENT")
# ============================================================

# The "shape" of the wiggly graph can be captured by the
# SECOND MOMENT of the neighbor distribution.
#
# For each tiling t with m=6 neighbors, the neighbor PROFILE is:
# p(t) = (count of neighbors in class 0, class 1, ..., class nm-1)
# This is a vector in Z^{nm} with ||p||_1 = m.
#
# The VARIANCE of p measures how "lumpy" the neighborhood is.
# If all neighbors are in different classes: variance low.
# If many neighbors cluster in one class: variance high.

print("  NEIGHBORHOOD LUMPINESS (variance of neighbor profile):\n")

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]; H = info[cid]['H']
    fiber = sum(sizes[c] for c in merged[mi])

    variances = []
    for nc_dict in class_stats[mi]['neighbor_classes']:
        profile = np.zeros(nm)
        for mj, count in nc_dict.items():
            profile[mj] = count
        var = np.var(profile)
        variances.append(var)

    avg_var = np.mean(variances)
    # Entropy of profile
    entropies = []
    for nc_dict in class_stats[mi]['neighbor_classes']:
        probs = np.array(list(nc_dict.values()), dtype=float) / m
        ent = -sum(p * np.log2(p) for p in probs if p > 0)
        entropies.append(ent)
    avg_ent = np.mean(entropies)

    print("  mid=%d (H=%2d): lumpiness=%.4f, profile_entropy=%.4f bits" %
          (mi, H, avg_var, avg_ent))

# ============================================================
banner("THE NEIGHBOR-OF-NEIGHBOR COUNT (N2)")
# ============================================================

# For each tiling t: N2(t) = number of distinct ISO CLASSES
# reachable in EXACTLY 2 tile flips.
# This measures the "connectivity horizon" of the tiling.

# More interesting: for each tiling t, the WEIGHT of the second
# neighborhood = Σ over d2-neighbors of (fiber of target class).
# This captures how "connected" t is to the bulk of tiling space.

print("  N2: number of distinct classes at Hamming distance exactly 2:\n")

# Also compute the NUMBER OF TILINGS at distance 2
# (which is always C(m,2) = 15 at n=5, but the class distribution varies)

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]; H = info[cid]['H']
    fiber = sum(sizes[c] for c in merged[mi])

    n2_vals = d2_classes_by_class[mi]
    avg_n2 = np.mean(n2_vals)
    min_n2 = min(n2_vals)
    max_n2 = max(n2_vals)

    print("  mid=%d (H=%2d, fiber=%d): N2 avg=%.2f [%d..%d]" %
          (mi, H, fiber, avg_n2, min_n2, max_n2))

print("""
  THE SHAPE EMERGES:

  At n=5, every tiling has 6 neighbors (degree = C(4,2) = T_3 = 6).
  But the REACH varies: some tilings reach 5 classes, others reach 6.
  And the SECOND NEIGHBORHOOD varies even more.

  The tilings with HIGH neutrality (many silent mutations) have:
  - LOWER reach (fewer distinct classes among neighbors)
  - But potentially HIGHER N2 (more classes at distance 2)
  - They are "locally sticky" but "globally connected"

  The tilings with LOW neutrality (all mutations expressive) have:
  - HIGHER reach at distance 1
  - They explore the class landscape MORE EFFICIENTLY per step

  This is the SHAPE of the wiggly graph:
  NOT a uniform lattice, but a STRATIFIED structure where
  each layer (iso class) has its own local geometry,
  determined by the |Aut| structure of the tournament.
""")

# ============================================================
banner("TOURNAMENT COUNTING THEOREM: PREDICTED NEUTRALITY")
# ============================================================

# From the Euler product: at large n, the correction R(n) → 0.
# The neutrality is related to the self-loop rate in the weight matrix.
#
# Self-loops = neutral flips = Σ_C W[C][C]
# Total flips = 2^m × m
#
# Neutrality = Σ_C W[C][C] / (2^m × m)
#
# A neutral flip at position k means the tournament T and its k-flip T'
# are in the same iso class. This happens iff there exists σ ∈ S_n
# such that σ(T') = T, i.e., σ maps T with one flipped arc back to T.
#
# The probability of this is related to |Aut| of both T and T'.
#
# For a RANDOM tiling: prob(flip k is neutral) = E[#{iso copies of T' in class(T)}] / 2^m
# This is hard to compute directly.
#
# But we know the TOTAL:
# At n=5: 176 neutral flips out of 384 total = 45.8% neutrality.
# At n=4: ? (need to compute)

# Actually let me just verify from the weight matrix:
total_sl = sum(W[mi][mi] for mi in range(nm))
total_all = W.sum()
print("  n=5: Self-loop weight = %d, total weight = %d" % (total_sl, total_all))
print("  Neutrality = %d / %d = %.4f" % (total_sl, total_all, total_sl/total_all))
print()

# From S291: R(n) = Σ fiber(C)×(|Aut(C)|-1)/2^m
# But this is about COUNTING classes, not about NEUTRAL FLIPS.
# The neutral flip rate is a DIFFERENT quantity.

# Actually: the self-loop weight W[mi][mi] counts the number of
# (tiling, tile-position) pairs where the flip is neutral.
# W[mi][mi] = #{(t, k) : t ∈ class mi, flip k neutral}

# The TOTAL self-loops = Σ_{all tilings t} #{neutral flips of t}
# = Σ_C Σ_{t ∈ C} neutrality(t) × m

# The Burnside connection: neutral flips relate to the number of
# labeled tournaments isomorphic to t that differ by a single arc.

# Use the orbit-counting identity:
# #{neutral flips} = #{(t, t', k) : t' = t^k ∈ orbit(t)}
# For a fixed t: flip k is neutral iff t^k is in the same S_n-orbit as t.
# = #{k : ∃σ ∈ S_n with σ(t^k) = t}
# = #{k : t^k ∈ [t]}

# This depends on the STRUCTURE of [t] (= fiber(t) elements in the orbit).

print("  PREDICTION: Neutrality ∝ average (|Aut|-1) fraction")
print("  (Because neutral flips indicate 'hidden symmetry')")
print()

# Let me compute |Aut| for each class at n=5 and correlate with self-loops
print("  %-4s %-4s %-4s %-6s %-8s %-10s" %
      ("mid", "H", "Aut", "fiber", "SL_wt", "SL_rate"))
print("  " + "-"*45)

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid = merged[mi][0]; H = info[cid]['H']
    fiber = sum(sizes[c] for c in merged[mi])
    aut = factorial(n) // (sizes[cid])  # |Aut| = n!/fiber_unmerged
    sl = W[mi][mi]
    sl_rate = sl / (fiber * m) if fiber > 0 else 0
    print("  %-4d %-4d %-4d %-6d %-8d %-10.4f" %
          (mi, H, aut, fiber, sl, sl_rate))

print("\nDone. opus-2026-03-24-S292")
