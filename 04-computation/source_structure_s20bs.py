#!/usr/bin/env python3
"""
source_structure_s20bs.py -- kind-pasteur-2026-03-22-S20bs

THE SOURCE'S BLUE LINE: Deep investigation.

From S20ak: flipping the anti-diagonal tile (range n-2) of the transitive
tournament gives a tournament with H = 1 + 2^(n-2) via a BLUE (SC-preserving) flip.

The user observes: this blue neighbor class has "1 + 2^(n-2) tilings."

Let me trace EXACTLY:
- How many LABELED tournaments are in this blue neighbor class?
- How many TILINGS (base-path representations) map to this class?
- What is the automorphism group?
- Can similar statements be made for OTHER blue neighbors?

Author: kind-pasteur-2026-03-22-S20bs
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

def is_sc(A, n):
    A_comp = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            if i != j: A_comp[i][j] = 1 - A[i][j]
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[j]] == A_comp[i][j] for i in range(n) for j in range(n) if i != j):
            return True
    return False

print("=" * 70)
print("  THE SOURCE'S BLUE LINE: DEEP STRUCTURE")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    non_base = [(i,j) for (i,j) in pairs if j > i+1]
    base_arcs = [(i,j) for (i,j) in pairs if j == i+1]

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

    classes = {}
    for cf, members in canon_map.items():
        H = H_map[members[0]]
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (members[0] >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        sc = is_sc(A, n) if n <= 6 else None
        cid = len(classes)
        classes[cf] = {'id': cid, 'H': H, 'size': len(members), 'sc': sc, 'members': set(members)}

    bits_to_class = {}
    for cf, info in classes.items():
        for b in info['members']:
            bits_to_class[b] = info

    # The transitive tournament: bits = 0
    trans_class = bits_to_class[0]
    print(f"  Transitive class: H={trans_class['H']}, size={trans_class['size']}, SC={trans_class['sc']}")

    # Find ALL blue neighbors of the transitive
    print(f"\n  BLUE NEIGHBORS of transitive (single arc flip, SC-preserving):")
    print(f"  {'Arc':>10s} {'Type':>6s} {'Range':>6s} {'H':>5s} {'ClassSize':>10s} {'SC':>4s} {'Color':>6s}")

    blue_classes = {}
    for k in range(m):
        arc = pairs[k]
        nb_bits = 0 ^ (1 << k)  # flip one bit of transitive (bits=0)
        nb_class = bits_to_class[nb_bits]
        is_base = arc in [(i, i+1) for i in range(n-1)]
        arc_type = "BASE" if is_base else "TILE"
        d = arc[1] - arc[0] - 1 if not is_base else 0
        color = "BLUE" if nb_class['sc'] == trans_class['sc'] else "BLACK"

        # Track unique blue class neighbors
        if color == "BLUE" and nb_class['id'] != trans_class['id']:
            blue_classes[nb_class['id']] = {
                'H': nb_class['H'], 'size': nb_class['size'],
                'sc': nb_class['sc'], 'range': d, 'arc': arc
            }

        if not is_base:  # only show tiles
            print(f"  ({arc[0]},{arc[1]}): {arc_type:>6s} {d:>6d} {nb_class['H']:>5d} {nb_class['size']:>10d} {str(nb_class['sc']):>4s} {color:>6s}")

    # The anti-diagonal blue neighbor
    print(f"\n  BLUE NEIGHBOR CLASSES (distinct, excluding self):")
    for cid, info in sorted(blue_classes.items(), key=lambda x: x[1]['range']):
        aut_size = factorial(n) // info['size']
        print(f"    Range {info['range']}: H={info['H']}, class_size={info['size']}, |Aut|={aut_size}, arc={info['arc']}")

    # Count tilings that map to each blue neighbor class
    # A tiling = a binary vector on non-base arcs (the cycle space)
    # Each tiling gives a specific labeled tournament (with base path 0->1->...->n-1)
    n_tiles = len(non_base)
    print(f"\n  TILINGS (base-path {list(range(n))}) mapping to blue neighbor classes:")

    for cid, info in sorted(blue_classes.items(), key=lambda x: x[1]['range']):
        # Count tilings whose tournament is in this class
        tiling_count = 0
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
            if cf in classes and classes[cf]['id'] == cid:
                tiling_count += 1

        formula = 1 + 2**(info['range']) if info['range'] > 0 else 1
        print(f"    Range {info['range']}: H={info['H']}, tilings={tiling_count}, 1+2^range={formula}, match={tiling_count==formula}")

    # THE KEY: the anti-diagonal blue neighbor
    if blue_classes:
        max_range_class = max(blue_classes.values(), key=lambda x: x['range'])
        d_max = max_range_class['range']
        formula = 1 + 2**d_max
        print(f"\n  ANTI-DIAGONAL (max range d={d_max}):")
        print(f"    H = 1 + 2^{d_max} = {formula}")
        print(f"    Class size = {max_range_class['size']}")
        print(f"    Expected 'tilings' = 1 + 2^(n-2) = {1 + 2**(n-2)}")

print(f"\n{'='*70}")
print(f"  SYNTHESIS")
print(f"{'='*70}\n")

print(f"""  THE SOURCE'S BLUE LINE STRUCTURE:

  From the transitive tournament, each non-base arc flip gives a
  tournament with H = 1 + 2^d (where d = range of the flipped arc).

  The BLUE flips (SC-preserving) go to specific classes.
  The number of TILINGS mapping to each blue neighbor class
  should tell us about the class's internal structure.

  The key observation: the anti-diagonal flip (maximum range d = n-2)
  always gives a BLUE neighbor with H = 1 + 2^(n-2).

  The number of tilings mapping to this class from the fixed base path
  reveals the TILING MULTIPLICITY of each iso class -- how many
  different base-path representations lead to the same iso class.
""")
