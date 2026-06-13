#!/usr/bin/env python3
"""
edge_burnside_s20bv.py -- kind-pasteur-2026-03-22-S20bv

EDGE COUNT VIA BURNSIDE ON EDGES.

The BYPASS: instead of finding a formula for the edge count,
use Burnside on the EDGE SET of the hypercube Q_m.

Q_m has 2^m vertices (tournaments) and m*2^{m-1} edges (arc flips).
S_n acts on both vertices and edges.

Edge orbits under S_n = number of "distinct flip types" at the iso class level.
But we want QUOTIENT EDGES (distinct pairs of iso classes connected by a flip),
not edge orbits.

The BYPASS idea: the quotient graph Q_m / S_n has:
  - Vertices = vertex orbits of S_n on Q_m = A000568(n) (iso classes)
  - Edges = pairs of vertex orbits connected by an edge in Q_m

Two vertex orbits C1, C2 are connected iff there exists a tournament T in C1
and an arc flip e such that flip(T,e) is in C2. This is NOT the same as
edge orbits (which may split the connection between C1 and C2 into multiple orbits).

BUT: we can compute the degree of each iso class directly.
degree(C) = number of DISTINCT iso classes reachable by flipping one arc of any T in C.

And: degree(C) depends on the AUTOMORPHISM GROUP Aut(T) and the
action of Aut(T) on the m arcs of T.

Specifically: degree(C) = m - |self-loops| - |duplicates|
where self-loops = arcs whose flip stays in C,
and duplicates = arcs whose flips land in the same neighbor class.

Author: kind-pasteur-2026-03-22-S20bv
"""
import sys
import numpy as np
from math import comb, factorial
from collections import defaultdict
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
print("  EDGE COUNT VIA DEGREE DECOMPOSITION")
print("=" * 70)

# The idea: for each iso class C, compute:
#   m arcs total
#   s(C) = arcs whose flip stays in C (self-loops)
#   The remaining m - s(C) arcs go to other classes
#   But multiple arcs may go to the SAME other class
#   degree(C) = number of DISTINCT neighbor classes

# The degree depends on the Aut(T) action on arcs.
# Aut(T) permutes the m arcs. The arcs partition into Aut-orbits.
# Within each Aut-orbit, ALL arcs go to the SAME neighbor class
# (because automorphisms map flips to flips within the same class pair).

# So: degree(C) = number of Aut-ORBITS on arcs that go to DIFFERENT classes.

# Total edges = (1/2) * sum over classes of degree(C)
# (divide by 2 because each edge is counted from both endpoints)

for n in [3, 4, 5]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    # Build iso classes
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
        classes.append({'id': cid, 'H': H_map[members[0]], 'size': len(members),
                       'members': set(members)})

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    N = len(classes)

    # For each class: compute degree, self-loops, Aut-orbit structure
    print(f"  {'Class':>6s} {'H':>4s} {'Size':>6s} {'|Aut|':>6s} {'ArcOrb':>7s} {'SelfOrb':>8s} {'CrossOrb':>9s} {'Degree':>7s}")

    total_degree = 0
    for c in classes:
        aut_size = factorial(n) // c['size']
        T = list(c['members'])[0]  # representative

        # Compute Aut(T) explicitly
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (T >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1

        aut_perms = []
        for perm in permutations(range(n)):
            if all(A[perm[i]][perm[j]] == A[i][j] for i in range(n) for j in range(n)):
                aut_perms.append(perm)

        # Aut-orbits on arcs: group arcs by Aut action
        arc_orbits = []
        arc_visited = [False] * m
        for k in range(m):
            if arc_visited[k]: continue
            orbit = set()
            for perm in aut_perms:
                # perm maps arc k=(i,j) to the arc (perm[i], perm[j])
                i, j = pairs[k]
                pi, pj = perm[i], perm[j]
                # Find which arc index this corresponds to
                if pi < pj:
                    new_arc = pairs.index((pi, pj))
                else:
                    new_arc = pairs.index((pj, pi))
                orbit.add(new_arc)
            for a in orbit:
                arc_visited[a] = True
            arc_orbits.append(orbit)

        # For each arc orbit: where does flipping go?
        self_orbs = 0
        cross_orbs = 0
        neighbor_classes = set()

        for orbit in arc_orbits:
            arc = list(orbit)[0]
            nb_bits = T ^ (1 << arc)
            nb_class = bits_to_class[nb_bits]
            if nb_class == c['id']:
                self_orbs += 1
            else:
                cross_orbs += 1
                neighbor_classes.add(nb_class)

        degree = len(neighbor_classes)
        total_degree += degree

        print(f"  {c['id']:>6d} {c['H']:>4d} {c['size']:>6d} {aut_size:>6d} {len(arc_orbits):>7d} {self_orbs:>8d} {cross_orbs:>9d} {degree:>7d}")

    edges = total_degree // 2
    print(f"\n  Total degree: {total_degree}")
    print(f"  Edges = total_degree / 2 = {edges}")
    print(f"  Expected: {[1, 5, 30, 290][[3,4,5,6].index(n)] if n in [3,4,5,6] else '?'}")
    print(f"  Match: {edges == [1, 5, 30, 290][[3,4,5,6].index(n)] if n in [3,4,5,6] else '?'}")

    # THE KEY FORMULA:
    # degree(C) = cross_orbits = m / |Aut| - self_loops / |Aut|
    # Wait: arc_orbits = m / |Aut| (by orbit-counting).
    # No: number of arc orbits = m only if |Aut| = 1.
    # By Burnside: #arc_orbits = (1/|Aut|) * sum_{g in Aut} Fix_arcs(g)
    # For |Aut|=1: #arc_orbits = m. For |Aut|=5 (regular): #arc_orbits = m/5 = 2.
    print(f"\n  ARC ORBIT ANALYSIS:")
    for c in classes:
        aut_size = factorial(n) // c['size']
        predicted_orbits = m // aut_size if m % aut_size == 0 else f"~{m/aut_size:.1f}"
        print(f"    Class {c['id']}: |Aut|={aut_size}, arc_orbits=?, m/|Aut|={m/aut_size:.1f}")

print(f"\n{'='*70}")
print(f"  THE DEGREE FORMULA")
print(f"{'='*70}\n")

print(f"""  For each iso class C with automorphism group Aut(C):

  1. The m arcs partition into Aut-ORBITS.
     #orbits = (1/|Aut|) * sum_{{g in Aut}} Fix_arcs(g)
     For |Aut|=1: #orbits = m (every arc is its own orbit)
     For |Aut|=k: #orbits < m (arcs are grouped)

  2. Each arc orbit maps to a SINGLE neighbor class (all arcs in
     the orbit go to the same class because Aut maps them to each other).

  3. Some orbits are SELF-LOOPS (flip stays in C).
     Others are CROSS-ORBITS (flip goes to a different class).

  4. degree(C) = number of DISTINCT classes among the cross-orbits.
     Multiple cross-orbits can go to the SAME neighbor, reducing degree.

  5. edges(G_n) = (1/2) * sum_C degree(C).

  THE BYPASS: We don't need a closed-form formula for edges.
  We can compute degree(C) from Aut(C) and the orbit structure.
  This is O(|S_n| * |Aut|) per class, much faster than brute force.

  For n=7 (456 classes): feasible with this approach.
  For n=8+ (6880+ classes): might need Aut-based shortcuts.
""")
