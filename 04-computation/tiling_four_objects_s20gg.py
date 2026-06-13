#!/usr/bin/env python3
"""
tiling_four_objects_s20gg.py — Each tiling corresponds to four objects
kind-pasteur-2026-03-25-S20gg

One tiling (C(n-1,2) bits on delta_{n-2}) gives:

  PATH-FIXED (the tiling + base path together):
    1. TOURNAMENT: directed complete graph with base path forward
    2. EULER GRAPH: XOR of fundamental cycles (via cycle space)

  PATH-FREE (forget the base path):
    3. TOURNAMENT ISO CLASS: the tournament's isomorphism class
    4. EULER GRAPH ISO CLASS: the Euler graph's isomorphism class

The COMPLEMENT tiling (flip all bits) gives:
    1'. COMPLEMENT TOURNAMENT (same base path, all tiles flipped)
    2'. COMPLEMENT EULER GRAPH (K_n XOR original in cycle space)
    3'. COMPLEMENT TOURNAMENT CLASS
    4'. COMPLEMENT EULER GRAPH CLASS

So each tiling really gives TWO tournaments and TWO Euler graphs
(original + complement), connected by the Z_2 bit-flip involution.

QUESTION: How do the four objects (tournament class, complement class,
Euler class, complement Euler class) relate across all tilings?
Is there a 4-way correspondence?
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  ONE TILING -> FOUR OBJECTS")
print("  kind-pasteur-2026-03-25-S20gg")
print("=" * 80)

for n in [4, 5]:
    t0 = time.time()
    N = n
    edges = [(i,j) for i in range(N) for j in range(i+1,N)]
    m_full = len(edges)
    base_edges = [(i,i+1) for i in range(N-1)]
    non_base = [e for e in edges if e not in base_edges]
    m = len(non_base)
    all_perms = list(permutations(range(N)))
    edge_idx = {e: k for k, e in enumerate(edges)}

    # Fundamental cycles
    fund_cycles = {}
    for e in non_base:
        i, j = e
        cycle_bits = 0
        for k in range(i, j):
            cycle_bits |= (1 << edge_idx[(k, k+1)])
        cycle_bits |= (1 << edge_idx[(i, j)])
        fund_cycles[e] = cycle_bits

    def tiling_to_tournament(tile_bits):
        """Tournament adjacency from tiling."""
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for k, (i,j) in enumerate(non_base):
            if tile_bits & (1 << k):
                A[j][i] = 1
            else:
                A[i][j] = 1
        for k in range(N-1):
            A[k+1][k] = 0
        return A

    def tiling_to_euler(tile_bits):
        """Euler graph from tiling via cycle space."""
        result = 0
        for k, e in enumerate(non_base):
            if tile_bits & (1 << k):
                result ^= fund_cycles[e]
        return result

    def tournament_canon(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def graph_canon(bits):
        best = None
        for p in all_perms:
            nb = 0
            for k, (i,j) in enumerate(edges):
                pi, pj = min(p[i],p[j]), max(p[i],p[j])
                nk = edge_idx[(pi, pj)]
                if bits & (1 << k): nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    print(f"\n{'#'*60}")
    print(f"  n = {n}, tiles = {m}")
    print(f"{'#'*60}")

    # For each tiling: compute all four objects + complements
    quadruples = []  # (tourn_cn, tourn_comp_cn, euler_cn, euler_comp_cn)

    for tile_bits in range(1 << m):
        comp_bits = tile_bits ^ ((1 << m) - 1)

        # Tournaments
        T = tiling_to_tournament(tile_bits)
        T_comp = tiling_to_tournament(comp_bits)
        t_cn = tournament_canon(T)
        tc_cn = tournament_canon(T_comp)

        # Euler graphs
        E = tiling_to_euler(tile_bits)
        E_comp = tiling_to_euler(comp_bits)
        e_cn = graph_canon(E)
        ec_cn = graph_canon(E_comp)

        quadruples.append((t_cn, tc_cn, e_cn, ec_cn))

    # How many distinct quadruples?
    distinct_quads = len(set(quadruples))
    print(f"\n  Distinct (T, T_comp, E, E_comp) quadruples: {distinct_quads}")

    # Relationship between tournament class and Euler class
    tourn_to_euler = defaultdict(set)
    euler_to_tourn = defaultdict(set)
    for t_cn, tc_cn, e_cn, ec_cn in quadruples:
        tourn_to_euler[t_cn].add(e_cn)
        euler_to_tourn[e_cn].add(t_cn)

    tourn_classes = set(t for t, _, _, _ in quadruples)
    euler_classes = set(e for _, _, e, _ in quadruples)

    print(f"\n  Tournament classes: {len(tourn_classes)}")
    print(f"  Euler graph classes: {len(euler_classes)}")
    print(f"\n  TOURNAMENT -> EULER MAP:")
    for t_cn in sorted(tourn_classes):
        euler_set = tourn_to_euler[t_cn]
        print(f"    T class -> {len(euler_set)} Euler classes")

    print(f"\n  EULER -> TOURNAMENT MAP:")
    for e_cn in sorted(euler_classes):
        tourn_set = euler_to_tourn[e_cn]
        print(f"    E class -> {len(tourn_set)} Tournament classes")

    # Complement pairing
    tourn_comp_pairs = set()
    euler_comp_pairs = set()
    for t_cn, tc_cn, e_cn, ec_cn in quadruples:
        tourn_comp_pairs.add((min(t_cn, tc_cn), max(t_cn, tc_cn)))
        euler_comp_pairs.add((min(e_cn, ec_cn), max(e_cn, ec_cn)))

    tourn_sc = sum(1 for a,b in tourn_comp_pairs if a == b)
    euler_sc = sum(1 for a,b in euler_comp_pairs if a == b)

    print(f"\n  COMPLEMENT STRUCTURE:")
    print(f"    Tournament SC pairs: {tourn_sc}")
    print(f"    Euler SC pairs: {euler_sc}")
    print(f"    Tournament complement pairs: {len(tourn_comp_pairs)}")
    print(f"    Euler complement pairs: {len(euler_comp_pairs)}")

    # The CROSS-TABLE: (tournament class, euler class) pairs
    cross_table = Counter()
    for t_cn, tc_cn, e_cn, ec_cn in quadruples:
        cross_table[(t_cn, e_cn)] += 1

    print(f"\n  CROSS-TABLE (tournament x euler):")
    print(f"    Distinct (T-class, E-class) pairs: {len(cross_table)}")
    print(f"    If independent: {len(tourn_classes)} x {len(euler_classes)} = {len(tourn_classes) * len(euler_classes)}")
    print(f"    Observed / independent: {len(cross_table) / (len(tourn_classes) * len(euler_classes)):.3f}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
