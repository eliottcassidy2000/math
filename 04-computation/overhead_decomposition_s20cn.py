#!/usr/bin/env python3
"""
overhead_decomposition_s20cn.py -- kind-pasteur-2026-03-22-S20cn

EXACT DECOMPOSITION of the overhead T_n - 2|E_n|.

The overhead = SL_orbits + MW where:
  SL_orbits = self-loop arc orbits (flipping arc stays in same class)
  MW = multi-weight correction (multiple arc orbits connecting same pair of classes)

For each iso class C:
  Arc orbits under Aut(C) partition into:
    - self-loop orbits (go back to C)
    - cross-orbits to distinct neighbors (contribute to degree)
    - cross-orbits duplicating a neighbor (contribute to MW)

  So per class: arc_orbits = SL_orbs + degree + MW_orbs

  Summing: T = Σ arc_orbits = Σ SL_orbs + Σ degree + Σ MW_orbs
           T = SL_total + 2|E| + MW_total
  => T - 2|E| = SL_total + MW_total

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys
import numpy as np
from math import comb, factorial
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
print("  OVERHEAD DECOMPOSITION: T_n - 2|E_n| = SL + MW")
print("  kind-pasteur-2026-03-22-S20cn")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    # Build iso classes
    canon_map = defaultdict(list)
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        cf = canonical(A, n)
        canon_map[cf].append(bits)

    classes = []
    for cf, members in canon_map.items():
        cid = len(classes)
        aut_size = factorial(n) // len(members)
        classes.append({'id': cid, 'aut': aut_size, 'size': len(members),
                       'members': set(members)})

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    N = len(classes)

    # For each class: compute arc orbits under Aut, classify each
    total_SL = 0
    total_MW = 0
    total_degree = 0
    total_arc_orbits = 0

    SL_by_aut = Counter()  # SL orbits by |Aut| value
    MW_by_aut = Counter()

    for c in classes:
        T = list(c['members'])[0]  # representative

        # Build adjacency
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (T >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1

        # Find Aut(T)
        aut_perms = []
        for perm in permutations(range(n)):
            if all(A[perm[i]][perm[j]] == A[i][j] for i in range(n) for j in range(n)):
                aut_perms.append(perm)

        # Compute arc orbits under Aut
        # Map each arc index to its canonical representative under Aut
        arc_orbit_map = {}  # arc_index -> orbit representative
        for k in range(m):
            if k in arc_orbit_map:
                continue
            orbit = set()
            for perm in aut_perms:
                i, j = pairs[k]
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    new_arc = pairs.index((pi, pj))
                else:
                    new_arc = pairs.index((pj, pi))
                orbit.add(new_arc)
            rep = min(orbit)
            for a in orbit:
                arc_orbit_map[a] = rep

        # Get distinct orbits
        orbit_reps = sorted(set(arc_orbit_map.values()))
        n_orbits = len(orbit_reps)

        # For each orbit: which class does flipping go to?
        orbit_targets = {}
        for rep in orbit_reps:
            nb_bits = T ^ (1 << rep)
            nb_class = bits_to_class[nb_bits]
            orbit_targets[rep] = nb_class

        # Count self-loop orbits, distinct neighbors, multi-weight
        self_loop_orbits = sum(1 for r in orbit_reps if orbit_targets[r] == c['id'])
        neighbor_classes = set(orbit_targets[r] for r in orbit_reps if orbit_targets[r] != c['id'])
        degree = len(neighbor_classes)
        cross_orbits = n_orbits - self_loop_orbits
        multi_weight = cross_orbits - degree

        total_SL += self_loop_orbits
        total_MW += multi_weight
        total_degree += degree
        total_arc_orbits += n_orbits

        SL_by_aut[c['aut']] += self_loop_orbits
        MW_by_aut[c['aut']] += multi_weight

    edges = total_degree // 2
    overhead = total_SL + total_MW
    T_n = total_arc_orbits

    print(f"  N = {N} iso classes")
    print(f"  T_n = {T_n} transition orbits")
    print(f"  |E| = {edges} edges")
    print(f"  2|E| = {total_degree}")
    print(f"  overhead = T - 2|E| = {overhead}")
    print(f"    SL (self-loop orbits) = {total_SL}")
    print(f"    MW (multi-weight)     = {total_MW}")
    print(f"    SL + MW = {total_SL + total_MW}")
    print(f"  Check: T = SL + 2E + MW = {total_SL + total_degree + total_MW} (should be {T_n})")

    # By |Aut| value
    print(f"\n  Breakdown by |Aut|:")
    all_auts = sorted(set(c['aut'] for c in classes))
    for a in all_auts:
        count = sum(1 for c in classes if c['aut'] == a)
        sl = SL_by_aut.get(a, 0)
        mw = MW_by_aut.get(a, 0)
        print(f"    |Aut|={a}: {count} classes, SL={sl}, MW={mw}")

    # Self-loop analysis
    print(f"\n  Self-loop analysis:")
    print(f"    SL_n / V_n = {total_SL / N:.4f}")
    print(f"    MW_n / V_n = {total_MW / N:.4f}")
    print(f"    SL_n / T_n = {total_SL / T_n:.4f}")
    print(f"    MW_n / T_n = {total_MW / T_n:.4f}")

# Summary table
print(f"\n{'='*70}")
print(f"  SUMMARY TABLE")
print(f"{'='*70}\n")
print(f"  {'n':>3s} {'V':>6s} {'m':>3s} {'T':>6s} {'|E|':>6s} {'OH':>5s} {'SL':>4s} {'MW':>5s} {'SL/V':>7s} {'MW/V':>7s}")
