#!/usr/bin/env python3
"""
analytic_edges_s20bw.py -- kind-pasteur-2026-03-22-S20bw

ANALYTICAL EDGE THEORY: Deriving everything about G_n edges.

KNOWN:
  edges(G_n) = (1/2) * sum_C degree(C)
  degree(C) = #{distinct neighbor classes from arc flips}
  For |Aut|=1: degree = m - self_loops (all arcs are their own orbit)
  For |Aut|>1: degree = #{Aut-orbits on arcs} - self_orbit_count

KEY IDENTITY: sum_C degree(C) * |C| = total_cross_flips
  where total_cross_flips = 2^m * m * (1 - f(n))
  and f(n) = (1/2)_{n-2}/(n-2)! is the self-loop fraction.

So: sum_C degree(C) = total_cross_flips / avg_orbit_size
    = 2^m * m * (1-f(n)) / ???

The missing piece: the WEIGHTED average of degree over classes.

APPROACH: Instead of edges, compute the SELF-LOOP COUNT analytically.
self_loops(C) = #{arcs whose flip stays in C} = #{arcs (a,b) such that
  the tournament obtained by flipping (a,b) is isomorphic to T}.

For a tournament T with |Aut(T)| = 1:
  Flipping arc (a,b) gives T'. T' is isomorphic to T iff
  there exists a permutation sigma such that sigma(T') = T.
  This sigma sends the flipped arc back, effectively "undoing" the flip
  via a relabeling. sigma is NOT an automorphism of T (since |Aut|=1),
  but it IS a "near-automorphism" that works modulo the single flip.

Author: kind-pasteur-2026-03-22-S20bw
"""
import sys
import numpy as np
from math import comb, factorial
from fractions import Fraction
from collections import defaultdict, Counter
from itertools import permutations
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

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

print("=" * 70)
print("  ANALYTICAL EDGE THEORY FOR G_n")
print("=" * 70)

# ================================================================
# 1. THE SELF-LOOP DECOMPOSITION
# ================================================================
print(f"\n{'='*70}")
print(f"  1. SELF-LOOP ANALYSIS: WHAT DETERMINES SELF-LOOPS?")
print(f"{'='*70}\n")

for n in [4, 5]:
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    canon_map = defaultdict(list)
    H_map = {}
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        H = count_hp(A, n)
        H_map[bits] = H
        cf = canonical(A, n)
        canon_map[cf].append(bits)

    classes = []
    cf_to_id = {}
    for cf, members in sorted(canon_map.items(), key=lambda x: H_map[x[1][0]]):
        cid = len(classes)
        cf_to_id[cf] = cid
        classes.append({'id': cid, 'H': H_map[members[0]], 'size': len(members), 'members': set(members)})

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    print(f"  n={n}:")
    print(f"  {'Class':>6s} {'H':>4s} {'|Aut|':>6s} {'Self':>5s} {'Cross':>6s} {'Degree':>7s} {'Self/m':>7s}")

    total_self = 0
    total_cross = 0
    for c in classes:
        T = list(c['members'])[0]
        aut_size = factorial(n) // c['size']

        # Count self-loops for this representative
        self_count = 0
        neighbor_classes = set()
        for k in range(m):
            nb = T ^ (1 << k)
            if bits_to_class[nb] == c['id']:
                self_count += 1
            else:
                neighbor_classes.add(bits_to_class[nb])

        degree = len(neighbor_classes)
        cross = m - self_count
        total_self += self_count * c['size']
        total_cross += cross * c['size']

        print(f"  {c['id']:>6d} {c['H']:>4d} {aut_size:>6d} {self_count:>5d} {cross:>6d} {degree:>7d} {self_count/m:>7.3f}")

    total_flips = 2**m * m
    self_frac = Fraction(total_self, total_flips)
    print(f"\n  Total self-loops: {total_self}/{total_flips} = {self_frac}")
    print(f"  Expected (1/2)_k/k!: {Fraction(comb(2*(n-2), n-2), 4**(n-2))}")
    print(f"  Match: {self_frac == Fraction(comb(2*(n-2), n-2), 4**(n-2))}")

# ================================================================
# 2. SELF-LOOPS BY SCORE CLASS
# ================================================================
print(f"\n{'='*70}")
print(f"  2. SELF-LOOPS CORRELATE WITH SCORE SEQUENCE")
print(f"{'='*70}\n")

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Rebuild
canon_map = defaultdict(list)
H_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H_map[bits] = count_hp(A, n)
    s = tuple(sorted(A.sum(axis=1).astype(int)))
    cf = canonical(A, n)
    canon_map[cf].append((bits, s))

classes = []
for cf, members in sorted(canon_map.items(), key=lambda x: H_map[x[0][0]]):
    bits0, score = members[0]
    cid = len(classes)
    classes.append({'id': cid, 'H': H_map[bits0], 'size': len(members),
                   'score': score, 'bits': bits0, 'members': set(b for b, _ in members)})

bits_to_class = {}
for c in classes:
    for b in c['members']:
        bits_to_class[b] = c['id']

print(f"  n={n}: Self-loops by score class")
print(f"  {'Score':>20s} {'H':>4s} {'|Aut|':>6s} {'Self':>5s} {'S/m':>6s}")

score_self = defaultdict(list)
for c in classes:
    T = c['bits']
    aut_size = factorial(n) // c['size']
    self_count = sum(1 for k in range(m) if bits_to_class[T ^ (1 << k)] == c['id'])
    score_self[c['score']].append((self_count, c['H'], aut_size))
    print(f"  {str(list(c['score'])):>20s} {c['H']:>4d} {aut_size:>6d} {self_count:>5d} {self_count/m:>6.3f}")

# ================================================================
# 3. THE ANALYTICAL FORMULA ATTEMPT
# ================================================================
print(f"\n{'='*70}")
print(f"  3. CAN WE PREDICT SELF-LOOPS FROM SCORE + |Aut|?")
print(f"{'='*70}\n")

# Self-loop count for |Aut|=1 classes:
aut1_selfs = [self_count for c in classes
              if factorial(n)//c['size'] == 1
              for self_count in [sum(1 for k in range(m) if bits_to_class[c['bits'] ^ (1 << k)] == c['id'])]]

print(f"  Self-loop counts for |Aut|=1 classes: {sorted(set(aut1_selfs))}")
print(f"  Values: {Counter(aut1_selfs)}")
print(f"  Range: [{min(aut1_selfs)}, {max(aut1_selfs)}]")

# For |Aut|=3:
aut3_selfs = [sum(1 for k in range(m) if bits_to_class[c['bits'] ^ (1 << k)] == c['id'])
              for c in classes if factorial(n)//c['size'] == 3]
print(f"  Self-loop counts for |Aut|=3 classes: {aut3_selfs}")

# For |Aut|=5:
aut5_selfs = [sum(1 for k in range(m) if bits_to_class[c['bits'] ^ (1 << k)] == c['id'])
              for c in classes if factorial(n)//c['size'] == 5]
print(f"  Self-loop counts for |Aut|=5 classes: {aut5_selfs}")

# ================================================================
# 4. THE EDGE FORMULA DECOMPOSITION
# ================================================================
print(f"\n{'='*70}")
print(f"  4. EDGE COUNT DECOMPOSITION")
print(f"{'='*70}\n")

# edges = (1/2) * sum degree(C)
# degree(C) = m - self_loops(C) - duplicate_neighbors(C)
# where duplicate_neighbors = #{arcs going to the same neighbor as another arc}

# For |Aut|=1: no arc orbit grouping, so duplicates come from
# distinct arcs whose flips produce ISOMORPHIC tournaments in DIFFERENT classes.
# Wait: distinct arcs CAN go to the same class even without Aut grouping.

# Let me compute: for each class, how many DISTINCT classes are reached?
print(f"  {'Class':>6s} {'H':>4s} {'|Aut|':>6s} {'Deg':>4s} {'Self':>5s} {'Cross':>6s} {'Dups':>5s}")

for c in classes:
    T = c['bits']
    aut_size = factorial(n) // c['size']
    self_count = 0
    neighbor_counts = defaultdict(int)

    for k in range(m):
        nb = T ^ (1 << k)
        nb_class = bits_to_class[nb]
        if nb_class == c['id']:
            self_count += 1
        else:
            neighbor_counts[nb_class] += 1

    degree = len(neighbor_counts)
    cross = m - self_count
    dups = cross - degree  # arcs going to same class as another arc

    print(f"  {c['id']:>6d} {c['H']:>4d} {aut_size:>6d} {degree:>4d} {self_count:>5d} {cross:>6d} {dups:>5d}")

# ================================================================
# 5. THE KEY ANALYTICAL RELATIONSHIPS
# ================================================================
print(f"\n{'='*70}")
print(f"  5. ANALYTICAL RELATIONSHIPS")
print(f"{'='*70}\n")

# Relationship 1: total self-loops = 2^m * m * f(n)
total_self = sum(
    sum(1 for k in range(m) if bits_to_class[list(c['members'])[0] ^ (1 << k)] == c['id']) * c['size']
    for c in classes
)
predicted_self = 2**m * m * float(Fraction(comb(2*(n-2), n-2), 4**(n-2)))
print(f"  Total self-loops: {total_self}")
print(f"  Predicted (2^m * m * f(n)): {predicted_self:.0f}")
print(f"  Match: {total_self == int(predicted_self)}")

# Relationship 2: average self-loops per class (weighted by size)
avg_self_weighted = total_self / (2**m * m) * m
print(f"  Average self-loops per tournament: {avg_self_weighted:.4f}")
print(f"  = m * f(n) = {m * float(Fraction(comb(2*(n-2), n-2), 4**(n-2))):.4f}")

# Relationship 3: for |Aut|=1 classes, average degree
aut1_degrees = []
for c in classes:
    if factorial(n) // c['size'] != 1: continue
    T = c['bits']
    nbs = set()
    for k in range(m):
        nb_class = bits_to_class[T ^ (1 << k)]
        if nb_class != c['id']:
            nbs.add(nb_class)
    aut1_degrees.append(len(nbs))

print(f"\n  |Aut|=1 classes: {len(aut1_degrees)} classes")
print(f"    Degree range: [{min(aut1_degrees)}, {max(aut1_degrees)}]")
print(f"    Average degree: {sum(aut1_degrees)/len(aut1_degrees):.2f}")
print(f"    Sum of degrees: {sum(aut1_degrees)}")

# Total edge contribution from |Aut|=1 classes
aut1_edge_contrib = sum(aut1_degrees)
# Total edge contribution from |Aut|>1 classes
other_degrees = []
for c in classes:
    if factorial(n) // c['size'] == 1: continue
    T = c['bits']
    nbs = set()
    for k in range(m):
        nb_class = bits_to_class[T ^ (1 << k)]
        if nb_class != c['id']:
            nbs.add(nb_class)
    other_degrees.append(len(nbs))

other_edge_contrib = sum(other_degrees)
print(f"\n  |Aut|>1 classes: {len(other_degrees)} classes")
print(f"    Sum of degrees: {other_edge_contrib}")
print(f"\n  Total degree: {aut1_edge_contrib + other_edge_contrib}")
print(f"  Edges: {(aut1_edge_contrib + other_edge_contrib) // 2}")

# The proportion from |Aut|=1
print(f"\n  Degree from |Aut|=1: {aut1_edge_contrib}/{aut1_edge_contrib + other_edge_contrib} = {aut1_edge_contrib/(aut1_edge_contrib + other_edge_contrib):.3f}")
print(f"  At large n, almost all classes have |Aut|=1.")
print(f"  So: edges ~ (1/2) * A000568(n) * (m - avg_self_loops)")
print(f"  where avg_self_loops ~ m * f(n) = m * (1/2)_{{n-2}}/(n-2)!")
print(f"  giving: edges ~ (1/2) * A000568(n) * m * (1 - f(n))")
print(f"  = (1/2) * A000568(n) * C(n,2) * (1 - (1/2)_{{n-2}}/(n-2)!)")

# Test this approximation
for nn in [3, 4, 5, 6]:
    V = {3:2, 4:4, 5:12, 6:56}[nn]
    mm = comb(nn, 2)
    fn = float(Fraction(comb(2*(nn-2), nn-2), 4**(nn-2)))
    E_pred = 0.5 * V * mm * (1 - fn)
    E_actual = {3:1, 4:5, 5:30, 6:290}[nn]
    print(f"  n={nn}: predicted={E_pred:.0f}, actual={E_actual}, ratio={E_actual/E_pred:.3f}")
