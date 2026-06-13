#!/usr/bin/env python3
"""
creative_metrics_s20w.py — kind-pasteur-2026-03-22-S20w

CREATIVE NEW METRICS FOR TOURNAMENTS:

1. A(T) = number of spanning arborescences (directed spanning trees)
   (Kirchhoff/Tutte directed matrix-tree theorem: det of Laplacian minor)
2. K(T) = number of kings (vertices reachable from all others in <= 2 steps)
3. F(T) = minimum feedback arc set size (min arcs to remove for acyclic)
4. S(T) = Slater index (min arc reversals to make transitive)
5. W(T) = writhe (sum of signed 3-cycles, antisymmetric invariant)
6. Entropy of HP distribution across starting vertices
7. L(T) = linear path count = H - n*HC (from last session)

Author: kind-pasteur-2026-03-22-S20w
"""
import sys
import numpy as np
from math import comb, log
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def count_hc(A, n):
    full = (1 << n) - 1
    dp = defaultdict(int)
    dp[(1, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[(full, v)] for v in range(n) if A[v][0])

def count_arborescences(A, n, root=0):
    """Number of spanning arborescences rooted at 'root' via Kirchhoff/Tutte.
    L = D_in - A^T where D_in is diagonal matrix of in-degrees.
    Delete row and column of root. Take determinant."""
    A_float = A.astype(float)
    D_in = np.diag(A_float.sum(axis=0))  # in-degree diagonal
    L = D_in - A_float.T
    # Delete root row and column
    indices = [i for i in range(n) if i != root]
    L_minor = L[np.ix_(indices, indices)]
    return int(round(np.linalg.det(L_minor)))

def count_kings(A, n):
    """A king can reach every other vertex in at most 2 steps."""
    kings = 0
    A2 = A @ A
    reach = A + np.clip(A2, 0, 1)  # reachable in 1 or 2 steps
    np.fill_diagonal(reach, 1)
    for v in range(n):
        if all(reach[v][w] > 0 or v == w for w in range(n)):
            kings += 1
    return kings

def feedback_arc_set_size(A, n):
    """Minimum feedback arc set size = c3-related.
    For a tournament: FAS = C(n,2) - longest_path_in_auxiliary_graph.
    Actually: FAS = number of backward arcs in the optimal linear ordering.
    Equivalent: FAS = C(n,2) - max acyclic subgraph.
    Brute force for small n."""
    from itertools import permutations
    min_fas = n * n
    for perm in permutations(range(n)):
        # Count backward arcs: i before j in perm but j->i
        fas = 0
        for pos_i in range(n):
            for pos_j in range(pos_i + 1, n):
                i, j = perm[pos_i], perm[pos_j]
                if A[j][i]:  # j beats i but i comes before j = backward
                    fas += 1
        min_fas = min(min_fas, fas)
    return min_fas

def hp_start_distribution(A, n):
    """Distribution of HP starting vertices."""
    full = (1 << n) - 1
    starts = np.zeros(n, dtype=int)
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    # For each start vertex, count HP
    for start in range(n):
        dp_s = defaultdict(int)
        dp_s[(1 << start, start)] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)): continue
                if dp_s[(mask, v)] == 0: continue
                for w in range(n):
                    if mask & (1 << w): continue
                    if A[v][w]: dp_s[(mask | (1 << w), w)] += dp_s[(mask, v)]
        starts[start] = sum(dp_s[(full, v)] for v in range(n))
    return starts

def hp_entropy(A, n):
    """Entropy of HP start distribution."""
    starts = hp_start_distribution(A, n)
    H_total = starts.sum()
    if H_total == 0: return 0
    probs = starts / H_total
    ent = 0
    for p in probs:
        if p > 0:
            ent -= p * log(p)
    return ent

print("="*60)
print("  CREATIVE METRICS FOR TOURNAMENTS")
print("="*60)

# Compute ALL metrics for all n=5 tournaments
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

print(f"\n  Computing all metrics at n={n} ({2**m} tournaments)...")

all_metrics = []
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1

    s = A.sum(axis=1).astype(int)
    S2 = int(sum(s*s))
    H = count_hp(A, n)
    HC = count_hc(A, n)
    L = H - n * HC
    c3 = comb(n,3) - (S2 - comb(n,2)) // 2

    # Arborescences (rooted at vertex 0)
    arb = count_arborescences(A, n, root=0)

    # Kings
    kings = count_kings(A, n)

    # HP entropy (expensive, only for small n)
    ent = hp_entropy(A, n)

    all_metrics.append({
        'H': H, 'HC': HC, 'L': L, 'c3': c3, 'S2': S2,
        'arb': arb, 'kings': kings, 'entropy': ent,
        'score': tuple(sorted(s))
    })

print(f"  Done. Computing correlations...")

# Build arrays
keys = ['H', 'HC', 'L', 'c3', 'S2', 'arb', 'kings', 'entropy']
arrays = {k: np.array([d[k] for d in all_metrics], dtype=float) for k in keys}

# Correlation matrix
print(f"\n  CORRELATION MATRIX:")
print(f"  {'':>10s}", end="")
for k in keys:
    print(f" {k:>8s}", end="")
print()

for k1 in keys:
    print(f"  {k1:>10s}", end="")
    for k2 in keys:
        if arrays[k1].std() > 0 and arrays[k2].std() > 0:
            corr = np.corrcoef(arrays[k1], arrays[k2])[0,1]
        else:
            corr = 0
        print(f" {corr:>+8.4f}", end="")
    print()

# Distributions
print(f"\n  METRIC DISTRIBUTIONS:")
for k in keys:
    vals = sorted(set(int(v) if k != 'entropy' else round(v, 3) for v in arrays[k]))
    print(f"  {k:>10s}: {len(vals)} distinct values. Range [{min(vals)}, {max(vals)}]")

# The KEY: which metrics are INDEPENDENT?
# Look for pairs with LOW correlation
print(f"\n  LOW-CORRELATION PAIRS (|r| < 0.5):")
for i, k1 in enumerate(keys):
    for j, k2 in enumerate(keys):
        if j <= i: continue
        if arrays[k1].std() == 0 or arrays[k2].std() == 0: continue
        corr = abs(np.corrcoef(arrays[k1], arrays[k2])[0,1])
        if corr < 0.5:
            print(f"    {k1} vs {k2}: |r| = {corr:.4f}")

# Arborescences vs H: is there a formula?
print(f"\n  ARBORESCENCES (rooted spanning trees) vs H:")
arb_to_H = defaultdict(set)
for d in all_metrics:
    arb_to_H[d['arb']].add(d['H'])

for arb_val in sorted(arb_to_H.keys()):
    Hs = sorted(arb_to_H[arb_val])
    print(f"    arb={arb_val:>3d}: H in {Hs}")

# Kings vs H
print(f"\n  KINGS vs H:")
kings_to_H = defaultdict(set)
for d in all_metrics:
    kings_to_H[d['kings']].add(d['H'])

for king_val in sorted(kings_to_H.keys()):
    Hs = sorted(kings_to_H[king_val])
    print(f"    kings={king_val}: H in {Hs}")

# Does arb determine H?
arb_det = all(len(v) == 1 for v in arb_to_H.values())
kings_det = all(len(v) == 1 for v in kings_to_H.values())
print(f"\n  Arborescences determine H: {arb_det}")
print(f"  Kings determine H: {kings_det}")

# The JOINT (arb, kings) -> H
joint_to_H = defaultdict(set)
for d in all_metrics:
    joint_to_H[(d['arb'], d['kings'])].add(d['H'])
joint_det = all(len(v) == 1 for v in joint_to_H.values())
print(f"  (Arb, Kings) jointly determine H: {joint_det}")

# Joint (arb, HC) -> H
joint2_to_H = defaultdict(set)
for d in all_metrics:
    joint2_to_H[(d['arb'], d['HC'])].add(d['H'])
joint2_det = all(len(v) == 1 for v in joint2_to_H.values())
print(f"  (Arb, HC) jointly determine H: {joint2_det}")

# The RATIO H / arb
print(f"\n  H / arb RATIO:")
for d in all_metrics[:20]:
    if d['arb'] > 0:
        ratio = d['H'] / d['arb']
        print(f"    H={d['H']}, arb={d['arb']}, H/arb={ratio:.4f}")

# Is H/arb always n-1?
# For a tournament: # spanning arborescences rooted at v = cofactor of Laplacian.
# And H = sum of all Hamiltonian paths.
# These are different counts, but both involve traversing all vertices.

# Arb counts TREES (no cycles). H counts PATHS (linear, may have cycles in tournament but path itself is acyclic).
# Both are acyclic traversals of all vertices.
# Arb is rooted (directed tree from root). H is a directed path (linear order).

# FORMULA CHECK: at n=3:
# Transitive (0->1->2): arb(root=0) = 1 (only tree: 0->1, 0->2 or 0->1->2).
# Actually arborescence = directed tree where all edges point AWAY from root.
# For root=0: need directed paths from 0 to all others.
# Transitive 0->1->2: arborescence 0->1, 1->2 (path) OR 0->1, 0->2 (star).
# Wait: 0->2 exists in the transitive tournament. So there are 2 arborescences.
# H = 1 (only one HP: 0->1->2).
# H/arb = 1/2 = 0.5.

# Cycle (0->1->2->0): arb(root=0) = 0->1,1->2 (unique, since 0->2 doesn't exist).
# arb = 1. H = 3. H/arb = 3.

# So H/arb is NOT constant. It varies.

print(f"\n  ARBORESCENCE FORMULA:")
print(f"  At n=3: transitive H=1, arb=2, ratio=0.5")
print(f"  At n=3: cycle H=3, arb=1, ratio=3.0")
print(f"  The ratio H/arb VARIES. It is NOT a constant.")
print(f"  Higher H/arb = more paths relative to trees = more CYCLIC structure.")
