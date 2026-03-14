#!/usr/bin/env python3
"""
h21_contradiction.py — INVESTIGATE: if (α₁,α₂)=(10,0) exists at n=7,
then H should equal 21. But exhaustive says H≠21. Find the bug.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def fast_hp(A, n):
    """DP Hamiltonian path count."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if not c or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def get_cycle_sets(A, n):
    """Get odd-cycle vertex sets."""
    cs = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    cs.add(frozenset(verts))
                    break
    return cs

print("INVESTIGATING (α₁,α₂) = (10,0) CONTRADICTION")
print("=" * 60)

n = 7
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

found = 0
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Quick filter: count 3-cycles
    dc3 = 0
    for v0, v1, v2 in combinations(range(n), 3):
        if A[v0][v1] and A[v1][v2] and A[v2][v0]:
            dc3 += 1
        elif A[v0][v2] and A[v2][v1] and A[v1][v0]:
            dc3 += 1

    # No filter on dc3 — check all tournaments
    if dc3 > 15:  # Quick skip for very cyclic tournaments
        continue

    cs = get_cycle_sets(A, n)
    a1 = len(cs)

    if a1 != 10:
        continue

    # Count disjoint pairs
    cs_list = list(cs)
    a2 = sum(1 for i in range(len(cs_list)) for j in range(i+1, len(cs_list))
             if not (cs_list[i] & cs_list[j]))

    if a2 != 0:
        continue

    # FOUND (10,0) tournament! Check H.
    H = fast_hp(A, n)

    # Compute I(Ω,2) from the independence polynomial directly
    # Since all pairs intersect, the intersection graph is K₁₀
    # I(K₁₀, 2) = 1 + 10·2 = 21
    # But wait — are there independent TRIPLES? With α₂=0 (no disjoint pair),
    # certainly α₃=0. So I = 1 + 10x, I(2) = 21.

    # BUT: is the OCF using the FULL independence polynomial or just
    # the vertex sets? Let me check.

    # Actually, maybe the issue is that the OCF uses DIRECTED cycles
    # as separate elements, not vertex sets?

    # Count directed cycles (not vertex sets)
    dir_cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    dir_cycles.append(frozenset(verts))

    # Each directed cycle on a vertex set: for 3-vertices, exactly 1 direction
    # For 5-vertices: could be multiple directed Hamiltonian cycles

    # Count unique directed cycles (not just vertex sets)
    dir_cycle_count = len(dir_cycles)  # includes duplicates for same vertex set
    unique_vsets = set(dir_cycles)

    found += 1
    if found <= 5:
        print(f"\nTournament bits={bits}")
        print(f"  dc3={dc3}")
        print(f"  Cycle vertex sets (α₁={a1}): {[sorted(s) for s in cs_list]}")
        print(f"  Directed cycle count (with multiplicity): {dir_cycle_count}")
        print(f"  Unique vertex sets: {len(unique_vsets)}")
        print(f"  α₂ (disjoint pairs) = {a2}")
        print(f"  H (from DP) = {H}")
        print(f"  Expected H from I(K₁₀, 2) = 21")
        print(f"  MATCH: {'YES' if H == 21 else 'NO — DISCREPANCY!'}")

        if H != 21:
            print(f"  *** BUG FOUND! H={H} ≠ 21 ***")
            # Check if there are higher-order independent sets
            # Actually, I was wrong: α₂=0 means all pairs intersect.
            # Then the intersection graph IS K_{α₁}.
            # I(K_n, x) = 1 + n*x. I(K_10, 2) = 21.
            # So if H ≠ 21 with (α₁,α₂)=(10,0), then the OCF is wrong?!
            # NO — the OCF is proved. So either my (10,0) identification is wrong.

            # Let me double-check disjoint pairs
            print("  Checking ALL pairs for intersection:")
            for i in range(len(cs_list)):
                for j in range(i+1, len(cs_list)):
                    inter = cs_list[i] & cs_list[j]
                    if not inter:
                        print(f"    DISJOINT: {sorted(cs_list[i])} ∩ {sorted(cs_list[j])} = ∅")

            # Check cycle lengths
            print("  Cycle lengths:")
            for s in cs_list:
                print(f"    {sorted(s)} (length {len(s)})")

    if found >= 10:
        break

print(f"\nTotal (10,0) tournaments found in first {bits+1}: {found}")
