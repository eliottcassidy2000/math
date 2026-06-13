#!/usr/bin/env python3
"""
overlap_weight_H_formula.py -- kind-pasteur-2026-03-13-S60

Orbit-based alpha decomposition: express H via orbit pairwise-disjointness.

For circulant tournaments on Z_p, cycles are organized by Z_p orbits.
Each orbit O has p member vertex sets, all with the same n(O) directed cycles.

alpha_j = sum over j-tuples of pairwise-disjoint cycle instances,
where two cycles are disjoint iff their vertex sets don't overlap.

The orbit structure simplifies counting dramatically.
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total

def compute_H_heldkarp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


for p in [7, 11]:
    m = (p - 1) // 2
    S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
    A = build_adj(p, S_qr)
    H = compute_H_heldkarp(A, p)

    print(f"\n{'='*60}")
    print(f"  ORBIT-BASED ALPHA: Paley p={p}, H={H}")
    print(f"{'='*60}")

    # Get active vertex sets with their orbit and cycle count
    active_vsets = []  # (frozenset, k, nV)
    for k in range(3, p+1, 2):
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                active_vsets.append((fs, k, nc))

    N_total = sum(nc for _, _, nc in active_vsets)
    print(f"  {len(active_vsets)} active vertex sets, N_total = {N_total}")

    # Build conflict graph on DIRECTED CYCLES (not vertex sets)
    # Each active vertex set V with n(V) directed cycles contributes n(V) nodes.
    # Two nodes conflict iff their vertex sets overlap.
    # Since all n(V) cycles on V have the same vertex set, they all conflict
    # with each other AND with any cycle on an overlapping set.

    # So the conflict graph is: complete biclique between V and W whenever V & W != 0,
    # plus complete graph within each V (all n(V) cycles conflict pairwise).

    # Independent sets: pick at most one cycle from each vertex set,
    # and the vertex sets must be pairwise disjoint.

    # alpha_j = sum over j-element collections of pairwise-disjoint active V-sets
    #           of product n(V_i)

    # For j=0: alpha_0 = 1
    # For j=1: alpha_1 = sum n(V) = N_total

    # For j=2: alpha_2 = sum_{V1<V2, V1 cap V2 = 0} n(V1) * n(V2)
    alpha_2 = 0
    for i in range(len(active_vsets)):
        for j in range(i+1, len(active_vsets)):
            V1, _, n1 = active_vsets[i]
            V2, _, n2 = active_vsets[j]
            if not (V1 & V2):
                alpha_2 += n1 * n2

    print(f"  alpha_1 = {N_total}")
    print(f"  alpha_2 = {alpha_2}")

    H2 = 1 + 2 * N_total + 4 * alpha_2
    print(f"  H(j<=2) = 1 + 2*{N_total} + 4*{alpha_2} = {H2}")

    if p <= 11:
        # Compute alpha_3
        # Need triples of pairwise disjoint active V-sets
        alpha_3 = 0

        # Build list of disjoint pairs first for efficiency
        disjoint_pairs = []
        for i in range(len(active_vsets)):
            for j in range(i+1, len(active_vsets)):
                V1, _, n1 = active_vsets[i]
                V2, _, n2 = active_vsets[j]
                if not (V1 & V2):
                    disjoint_pairs.append((i, j))

        print(f"  Disjoint pairs: {len(disjoint_pairs)}")

        # For alpha_3: extend each disjoint pair with a third disjoint set
        for idx_pair, (i, j) in enumerate(disjoint_pairs):
            V1, k1, n1 = active_vsets[i]
            V2, k2, n2 = active_vsets[j]
            used = V1 | V2
            for l in range(j+1, len(active_vsets)):
                V3, k3, n3 = active_vsets[l]
                if not (used & V3):
                    alpha_3 += n1 * n2 * n3

        print(f"  alpha_3 = {alpha_3}")
        H3 = H2 + 8 * alpha_3
        print(f"  H(j<=3) = {H3}")

        if H3 != H:
            # Need alpha_4 at p=11? Max disjoint 3-subsets: floor(11/3)=3
            # So max j with 3-cycles only = 3.
            # But with 5-cycles: 3+5=8, 3+3+5=11 (tight packing).
            # Four pairwise disjoint: 3+3+3+3=12 > 11, impossible.
            # 3+3+5 = 11: uses all vertices, no room for 4th.
            # So alpha_4 should be 0 at p=11 for odd-k cycles.

            # Wait: 3+3+3=9 < 11. Can we fit a 4th? Need 2 vertices left.
            # No odd cycle on 2 vertices. But a 1-vertex "cycle"? No.
            # So alpha_4 = 0 when cycles use >= 3 vertices each.

            # Actually, this is wrong for small enough cycles.
            # We need j pairwise disjoint odd cycles, total vertices <= p.
            # With 3-vertex cycles: max j = floor(p/3) = floor(11/3) = 3.

            print(f"  Deficit: {H - H3}")
            print(f"  Max j at p={p}: floor({p}/3) = {p // 3}")
            print(f"  So alpha_{p//3 + 1} should be 0")
        else:
            print(f"  EXACT MATCH at j=3!")

    print(f"  H from Held-Karp: {H}")
    print(f"  Match: {H2 == H if p == 7 else H3 == H}")

# Now cross-orientation comparison at p=7
print(f"\n\n{'='*60}")
print(f"  ALL ORIENTATIONS at p=7")
print(f"{'='*60}")

p = 7
m = 3
for bits in range(1 << m):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)
    H = compute_H_heldkarp(A, p)

    active = []
    for k in range(3, p+1, 2):
        for subset in combinations(range(p), k):
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                active.append((frozenset(subset), k, nc))

    N_total = sum(nc for _, _, nc in active)

    alpha_2 = 0
    for i in range(len(active)):
        for j in range(i+1, len(active)):
            if not (active[i][0] & active[j][0]):
                alpha_2 += active[i][2] * active[j][2]

    H_check = 1 + 2 * N_total + 4 * alpha_2
    match = H_check == H

    # Breakdown
    c_k = defaultdict(int)
    for _, k, nc in active:
        c_k[k] += nc

    print(f"  bits={bits}: S={S}, H={H}, N={N_total}, "
          f"a2={alpha_2}, c_k={dict(sorted(c_k.items()))}, "
          f"H_check={H_check}, {'OK' if match else 'FAIL'}")

print("\nDONE.")
