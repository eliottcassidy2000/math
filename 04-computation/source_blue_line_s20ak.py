#!/usr/bin/env python3
"""
source_blue_line_s20ak.py -- kind-pasteur-2026-03-22-S20ak

THE SOURCE'S BLUE LINE: How does the transitive tournament connect
to the iso class graph, and what is the 1+2^(n-2) pattern?

The user observes: the transitive class (unique source) connects via
a single blue (SC-preserving) edge to a class containing 1+2^(n-2) tilings.
  n=3: 1+2=3
  n=4: 1+4=5
  n=5: 1+8=9
  n=6: 1+16=17

This session:
1. Map every arc flip from the transitive tournament to its iso class
2. Color each flip blue (SC-preserving) or black (SC-changing)
3. Trace the blue line from the source
4. Understand the tiling-to-iso-class map through edge colors

Author: kind-pasteur-2026-03-22-S20ak
"""
import sys
import numpy as np
from math import comb
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

def canonical_form(A, n):
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
print("  THE SOURCE'S BLUE LINE: TILING-TO-ISO-CLASS MAP")
print("=" * 70)

for n in [3, 4, 5, 6]:
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
        cf = canonical_form(A, n)
        canon_map[cf].append(bits)

    classes = []
    cf_to_id = {}
    for cf, members in sorted(canon_map.items(), key=lambda x: (H_map[x[1][0]], len(x[1]))):
        cid = len(classes)
        cf_to_id[cf] = cid
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (members[0] >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        score = tuple(sorted(A.sum(axis=1).astype(int)))
        sc = is_sc(A, n) if n <= 6 else None
        classes.append({
            'id': cid, 'H': H_map[members[0]], 'score': score,
            'size': len(members), 'sc': sc, 'members': set(members)
        })

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    N = len(classes)

    # Find the transitive tournament (all bits = 0 means i->j for all i<j)
    # With our encoding: bit k for pair (i,j) with i<j, bit=1 means j->i
    # So bits=0 is the transitive tournament i->j for i<j
    trans_bits = 0
    trans_class = bits_to_class[trans_bits]
    trans_H = H_map[trans_bits]
    trans_sc = classes[trans_class]['sc']

    print(f"  Transitive tournament: class {trans_class}, H={trans_H}, SC={trans_sc}, size={classes[trans_class]['size']}")

    # Find ALL single-arc-flip neighbors of the transitive tournament
    print(f"\n  Single arc flips from transitive (bits=0):")
    print(f"  {'arc':>10s} {'new_class':>10s} {'new_H':>6s} {'new_SC':>6s} {'class_size':>10s} {'color':>8s}")

    flip_to_class = {}
    flip_colors = defaultdict(int)  # (trans_class, neighbor_class) -> count of blue/black

    for k in range(m):
        nb_bits = trans_bits ^ (1 << k)
        nb_class = bits_to_class[nb_bits]
        nb_H = H_map[nb_bits]
        nb_sc = classes[nb_class]['sc']
        nb_size = classes[nb_class]['size']

        # Color: blue if both SC or both non-SC
        color = "BLUE" if (trans_sc == nb_sc) else "BLACK"

        arc = pairs[k]
        print(f"  ({arc[0]}->{arc[1]}): {nb_class:>10d} {nb_H:>6d} {str(nb_sc):>6s} {nb_size:>10d} {color:>8s}")

        flip_to_class[k] = nb_class

    # Count distinct neighbor classes
    neighbor_classes = set(flip_to_class.values())
    print(f"\n  Distinct neighbor classes: {len(neighbor_classes)}")

    # Group by color
    blue_neighbors = set()
    black_neighbors = set()
    for k in range(m):
        nb_class = flip_to_class[k]
        nb_sc = classes[nb_class]['sc']
        if trans_sc == nb_sc:
            blue_neighbors.add(nb_class)
        else:
            black_neighbors.add(nb_class)

    print(f"  Blue (SC-preserving) neighbor classes: {sorted(blue_neighbors)}")
    print(f"  Black (SC-changing) neighbor classes: {sorted(black_neighbors)}")
    print(f"  # Blue neighbors: {len(blue_neighbors)}")
    print(f"  # Black neighbors: {len(black_neighbors)}")

    # For each neighbor class, how many flips lead to it?
    print(f"\n  Flips per neighbor class:")
    class_flip_count = defaultdict(int)
    for k in range(m):
        class_flip_count[flip_to_class[k]] += 1

    for cid in sorted(neighbor_classes):
        c = classes[cid]
        color = "BLUE" if c['sc'] == trans_sc else "BLACK"
        print(f"    Class {cid}: {class_flip_count[cid]} flips, H={c['H']}, size={c['size']}, SC={c['sc']}, score={list(c['score'])}, color={color}")

    # Check the 1 + 2^(n-2) pattern
    target = 1 + 2**(n-2)
    print(f"\n  1 + 2^(n-2) = {target}")
    for cid in sorted(neighbor_classes):
        c = classes[cid]
        if c['size'] == target:
            color = "BLUE" if c['sc'] == trans_sc else "BLACK"
            print(f"  CLASS WITH SIZE {target} FOUND: class {cid}, H={c['H']}, SC={c['sc']}, color={color}")

    # Also check if any neighbor has size matching other patterns
    print(f"\n  All neighbor class sizes: {sorted(set(classes[cid]['size'] for cid in neighbor_classes))}")
    print(f"  n! = {[6,24,120,720][n-3]}")
    print(f"  n!/2 = {[6,24,120,720][n-3]//2}")

    # HOW MANY LABELED TOURNAMENTS come from flipping one arc of ANY transitive?
    # There are n! transitive tournaments (one per permutation).
    # Each has m arc-flip neighbors.
    # Each flip produces one labeled tournament.
    # So: n! * m labeled tournaments reachable in one flip from transitive.
    # But many are the same labeled tournament. Let's count distinct ones per class.

    # Actually: from the SPECIFIC transitive tournament (bits=0),
    # each of the m flips gives a distinct labeled tournament (all different bits).
    # These m tournaments distribute among the neighbor classes.
    # class_flip_count[cid] = how many of the m flips land in class cid.
    # class_size[cid] = total labeled tournaments in that class.
    # So: class_flip_count[cid] / class_size[cid] = fraction of class reachable from THIS transitive.

    print(f"\n  Reachability fractions (flips from trans / class size):")
    for cid in sorted(neighbor_classes):
        c = classes[cid]
        frac = class_flip_count[cid] / c['size']
        color = "BLUE" if c['sc'] == trans_sc else "BLACK"
        print(f"    Class {cid}: {class_flip_count[cid]}/{c['size']} = {frac:.4f}, color={color}")

    # The STAIRCASE TILING perspective:
    # In the staircase model, the transitive tournament is the "empty tiling"
    # (all tiles = 0). Flipping one tile changes one non-base-path arc.
    # But the base path arcs (i -> i+1 for i=0..n-2) are FIXED.
    # So flipping tile k corresponds to flipping arc (i,j) where (i,j) is
    # a non-consecutive pair (j > i+1).

    base_path_arcs = set()
    for i in range(n-1):
        base_path_arcs.add((i, i+1))

    non_base = [(i,j) for (i,j) in pairs if (i,j) not in base_path_arcs]
    base = [(i,j) for (i,j) in pairs if (i,j) in base_path_arcs]

    print(f"\n  STAIRCASE TILING VIEW:")
    print(f"  Base path arcs: {base} ({len(base)} arcs, fixed)")
    print(f"  Non-base arcs (tiles): {non_base} ({len(non_base)} = C(n-1,2) tiles)")
    print(f"  C(n-1,2) = {comb(n-1,2)}")

    # Which of the m flips are base-path flips vs tile flips?
    for k in range(m):
        arc = pairs[k]
        is_base = arc in base_path_arcs
        nb_class = flip_to_class[k]
        c = classes[nb_class]
        color = "BLUE" if c['sc'] == trans_sc else "BLACK"
        arc_type = "BASE" if is_base else "TILE"
        if k < 20 or is_base:  # print all base, some tiles
            print(f"    arc {arc}: {arc_type}, -> class {nb_class} (H={c['H']}, SC={c['sc']}, color={color})")

print(f"\n{'='*70}")
print(f"  SUMMARY: THE BLUE LINE FROM THE SOURCE")
print(f"{'='*70}\n")

for n in [3, 4, 5, 6]:
    target = 1 + 2**(n-2)
    print(f"  n={n}: 1+2^(n-2) = {target}")
