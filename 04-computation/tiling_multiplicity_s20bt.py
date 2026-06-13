#!/usr/bin/env python3
"""
tiling_multiplicity_s20bt.py -- kind-pasteur-2026-03-22-S20bt

H = TILING MULTIPLICITY: Extending and understanding.

From S20bs: for single-tile blue neighbor classes at range d >= 2,
the number of tilings = H = 1 + 2^d.

NOW: Does this generalize beyond single-tile classes?
For EVERY iso class, is H = number of tilings mapping to it?

If yes: H has a direct combinatorial interpretation as a FIBER COUNT
in the tiling fibration (staircase tilings -> iso classes).

Author: kind-pasteur-2026-03-22-S20bt
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
print("  H = TILING MULTIPLICITY? TESTING FOR ALL ISO CLASSES")
print("=" * 70)

for n in [4, 5]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    non_base = [(i,j) for (i,j) in pairs if j > i+1]
    n_tiles = len(non_base)

    # Build iso classes from ALL tournaments
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

    # For TILINGS: base path is 0->1->...->n-1
    # A tiling is a binary vector on non-base arcs
    # Each tiling produces ONE labeled tournament
    # Map each tiling to its iso class

    tiling_to_class = {}
    class_tiling_count = defaultdict(int)

    for tbits in range(2**n_tiles):
        # Build tournament from tiling
        A = np.zeros((n,n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                A[i][j] = 1  # start transitive
        for t in range(n_tiles):
            if (tbits >> t) & 1:
                a, b = non_base[t]
                A[a][b] = 0; A[b][a] = 1

        cf = canonical(A, n)
        cid = cf_to_id[cf]
        tiling_to_class[tbits] = cid
        class_tiling_count[cid] += 1

    # Compare: H vs tiling count for each class
    print(f"  {'Class':>6s} {'H':>5s} {'Tilings':>8s} {'ClassSize':>10s} {'H=Tilings?':>11s} {'Ratio':>8s}")

    all_match = True
    for c in classes:
        tc = class_tiling_count.get(c['id'], 0)
        match = tc == c['H']
        ratio = tc / c['H'] if c['H'] > 0 else 0
        if not match: all_match = False
        print(f"  {c['id']:>6d} {c['H']:>5d} {tc:>8d} {c['size']:>10d} {str(match):>11s} {ratio:>8.4f}")

    print(f"\n  H = tiling count for ALL classes: {all_match}")

    if not all_match:
        print(f"\n  MISMATCHES (H != tiling count):")
        for c in classes:
            tc = class_tiling_count.get(c['id'], 0)
            if tc != c['H']:
                print(f"    Class {c['id']}: H={c['H']}, tilings={tc}, diff={tc - c['H']}")

    # Summary statistics
    total_tilings = 2**n_tiles
    total_H = sum(c['H'] * c['size'] for c in classes) // factorial(n)
    print(f"\n  Total tilings: {total_tilings}")
    print(f"  Sum of H/class: {sum(class_tiling_count.values())}")
    print(f"  These should be equal: {total_tilings == sum(class_tiling_count.values())}")

    # Relationship between tiling count and class size
    print(f"\n  RELATIONSHIP: tiling_count vs H vs class_size")
    print(f"  {'Class':>6s} {'H':>5s} {'Tilings':>8s} {'Size':>6s} {'|Aut|':>6s} {'T/H':>6s} {'T*|Aut|':>8s}")
    for c in classes:
        tc = class_tiling_count.get(c['id'], 0)
        aut = factorial(n) // c['size']
        ratio = tc / c['H'] if c['H'] > 0 else 0
        t_aut = tc * aut
        print(f"  {c['id']:>6d} {c['H']:>5d} {tc:>8d} {c['size']:>6d} {aut:>6d} {ratio:>6.3f} {t_aut:>8d}")

print(f"\n{'='*70}")
print(f"  SYNTHESIS")
print(f"{'='*70}\n")

print(f"""  THE QUESTION: Does H(T) = number of tilings mapping to iso class of T?

  AT n=4: Check each class.
  AT n=5: Check each class.

  If H = tiling multiplicity for ALL classes:
    Then H has a direct combinatorial interpretation as a FIBER COUNT.
    The map (tilings) -> (iso classes) is a fibration, and H(T) = |fiber(T)|.

  If H != tiling multiplicity for some classes:
    Then the relationship is more subtle -- H = tiling count only
    for SPECIFIC classes (like the single-tile blue neighbors).
    The general relationship involves the base path choice.

  THE KEY: each tiling uses a SPECIFIC base path (0->1->...->n-1).
  Different base paths give different tiling-to-class maps.
  The TOTAL over all n! base paths should give: n! * H / something.

  Specifically: each labeled tournament appears as a tiling for exactly
  (n-1)! base paths (since the base path is a Hamiltonian path of the
  transitive tournament, and each starts at a different vertex).
  Wait: the base path is fixed as 0->1->...->n-1.
  A different base path would be a different Hamiltonian path.
  The transitive tournament has exactly 1 Hamiltonian path: 0->1->...->n-1.
  So the tiling model with THIS base path gives each tournament exactly ONCE
  if the tournament contains the path 0->1->...->n-1 as a Hamiltonian subpath.

  But NOT every tournament contains this specific HP!
  Tournaments that DON'T have 0->1->...->n-1 as an HP can't be reached
  from ANY tiling with this base path.
""")
