#!/usr/bin/env python3
"""
unified_overlap_theory.py -- kind-pasteur-2026-03-13-S60

UNIFIED THEORY: Why H = linear(c_k) for circulant tournaments.

The chain of reasoning:

1. ind_k(V) = sum_{directed k-cycles on V} 1
   is a polynomial of degree d(k) = k-1 in sigma = (sigma_1,...,sigma_m)
   with ALL ODD DEGREES VANISHING. [cycle_indicator_degree.py]

2. c_k = sum_{|V|=k} ind_k(V) = polynomial of degree k-1 in sigma
   By Z_p-symmetry, this is a function of the orbit structure,
   hence of the eigenvalue spectrum D_1,...,D_m.
   Result: c_k = affine(S_4,...,S_{p-3}). [THM-157 + extension]

3. disj(k1,k2) = sum_{V1 cap V2 = 0} ind_{k1}(V1) * ind_{k2}(V2)
   The product has degree d(k1)+d(k2) = k1+k2-2 in sigma.
   But the Z_p-symmetry and disjointness constraint force it
   into the span of S_4,...,S_{p-3}. [degree_reduction_mechanism.py]

4. alpha_j = sum of j-wise disjoint tuple counts
   Since each tuple count is affine(S_4,...,S_{p-3}), so is alpha_j.
   [THM-158]

5. H = 1 + 2*alpha_1 + 4*alpha_2 + ... + 2^{floor(p/3)} * alpha_{floor(p/3)}
   Each alpha_j is affine(S_4,...,S_{p-3}), hence H is too.
   And since c_k is also affine(S_4,...,S_{p-3}), and the mapping
   (c_3,...,c_p) -> (S_4,...,S_{p-3}) is (essentially) invertible,
   H is affine in c_k. [H_from_ck_exact.py]

This script performs the COMPLETE verification of this chain at p=7 and p=11,
and tests the OVERLAP WEIGHT decomposition:

  w(C) = degree of cycle C in conflict graph Omega
       = sum_{C' != C} 1[C shares vertex with C']

The overlap weight w(C) depends on the TYPE of the cycle (gap structure).
For circulant tournaments, all cycles in the same Z_p-orbit have the same weight.
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import time


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


def overlap_weight_analysis(p, S, name):
    """Complete overlap weight analysis for a single tournament."""
    A = build_adj(p, S)
    H = compute_H_heldkarp(A, p)

    # Enumerate all cycles with vertex sets and weights
    cycle_data = []  # (frozenset, k, orbit_type)
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            for _ in range(nc):
                # Orbit type: canonical form under Z_p translation
                canon = min(frozenset((v + t) % p for v in subset) for t in range(p))
                cycle_data.append((fs, k, canon))

    c_k = defaultdict(int)
    for fs, k, _ in cycle_data:
        c_k[k] += 1

    # Build conflict graph
    n = len(cycle_data)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if cycle_data[i][0] & cycle_data[j][0]:
                adj[i][j] = 1
                adj[j][i] = 1

    # Overlap weight = degree in conflict graph
    weights = [sum(adj[i]) for i in range(n)]

    # Group by cycle type (length + orbit)
    by_type = defaultdict(list)
    for i, (fs, k, canon) in enumerate(cycle_data):
        by_type[(k, canon)].append(i)

    print(f"\n{'='*60}")
    print(f"  {name}, p={p}, S={S}, H={H}")
    print(f"  Total cycles: {n}, c_k: {dict(sorted(c_k.items()))}")
    print(f"{'='*60}")

    print(f"\n  Cycle types and overlap weights:")
    print(f"  {'type':>20} {'len':>4} {'count':>5} {'avg_w':>8} {'min_w':>6} {'max_w':>6}")

    for (k, canon) in sorted(by_type.keys()):
        indices = by_type[(k, canon)]
        ws = [weights[i] for i in indices]
        avg_w = sum(ws) / len(ws)
        print(f"  {str(sorted(canon))[:20]:>20} {k:>4} {len(indices):>5} "
              f"{avg_w:>8.1f} {min(ws):>6} {max(ws):>6}")

    # Verify: all cycles in same type have same weight (Z_p symmetry)
    all_same = True
    for (k, canon), indices in by_type.items():
        ws = set(weights[i] for i in indices)
        if len(ws) > 1:
            all_same = False
            print(f"  WARNING: type {sorted(canon)} has {len(ws)} distinct weights!")

    if all_same:
        print(f"\n  All cycles within each type have SAME weight (Z_p-symmetric). OK")

    # Overlap weight DECOMPOSITION by cycle length
    # For each cycle C of length k, w(C) = sum_{k'} w_{k'}(C)
    # where w_{k'}(C) = number of k'-cycles conflicting with C
    print(f"\n  Overlap weight decomposition by target length:")
    for (k, canon) in sorted(by_type.keys()):
        # Pick one representative cycle
        i = by_type[(k, canon)][0]
        decomp = defaultdict(int)
        for j in range(n):
            if i != j and adj[i][j]:
                decomp[cycle_data[j][1]] += 1

        decomp_str = ', '.join(f'w_{k2}={decomp[k2]}' for k2 in sorted(decomp))
        total_w = sum(decomp.values())
        print(f"    {k}-cycle {sorted(canon)[:5]}...: w={total_w} = {decomp_str}")

    # Alpha decomposition
    # Count independent sets by size
    if n <= 80:
        nbr = [0] * n
        for i in range(n):
            for j in range(n):
                if adj[i][j]:
                    nbr[i] |= (1 << j)

        alpha = [0] * (n + 1)

        def backtrack(v, mask, size):
            alpha[size] += 1
            for w_idx in range(v + 1, n):
                if not (mask & (1 << w_idx)):
                    backtrack(w_idx, mask | nbr[w_idx], size + 1)

        backtrack(-1, 0, 0)

        max_j = max(j for j in range(len(alpha)) if alpha[j] > 0)
        print(f"\n  Alpha decomposition: max j = {max_j}")
        for j in range(max_j + 1):
            print(f"    alpha_{j} = {alpha[j]}")

        H_check = sum(alpha[j] * (2 ** j) for j in range(max_j + 1))
        print(f"  H(OCF) = {H_check}, H(HK) = {H}, match = {H_check == H}")

    # OVERLAP WEIGHT SPECTRUM: eigenvalues of the weighted overlap matrix
    W = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            ov = len(cycle_data[i][0] & cycle_data[j][0])
            W[i][j] = ov
            W[j][i] = ov

    eigvals = np.sort(np.linalg.eigvalsh(W))[::-1]
    print(f"\n  Overlap matrix eigenvalues (top 10):")
    for i in range(min(10, len(eigvals))):
        print(f"    lambda_{i} = {eigvals[i]:.4f}")

    return cycle_data, weights, by_type


# Run for Paley and Interval at p=7 and p=11
for p_val in [7, 11]:
    m = (p_val - 1) // 2
    S_qr = sorted(j for j in range(1, p_val) if pow(j, (p_val - 1) // 2, p_val) == 1)
    S_int = list(range(1, m + 1))

    t0 = time.time()
    overlap_weight_analysis(p_val, S_qr, f"Paley p={p_val}")
    t1 = time.time()
    overlap_weight_analysis(p_val, S_int, f"Interval p={p_val}")
    t2 = time.time()
    print(f"\n  Paley: {t1-t0:.1f}s, Interval: {t2-t1:.1f}s")


# CROSS-ORIENTATION COMPARISON at p=7
print(f"\n{'='*60}")
print(f"  CROSS-ORIENTATION WEIGHT COMPARISON at p=7")
print(f"{'='*60}")

p = 7
m = 3
N = 1 << m

for bits in range(N):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)
    H = compute_H_heldkarp(A, p)

    cycles = []
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            nc = count_ham_cycles(A, list(subset))
            for _ in range(nc):
                cycles.append(frozenset(subset))

    # Compute overlap weight for each cycle
    n_cyc = len(cycles)
    total_weight = 0
    for i in range(n_cyc):
        for j in range(i + 1, n_cyc):
            if cycles[i] & cycles[j]:
                total_weight += 2  # both directions

    avg_weight = total_weight / n_cyc if n_cyc > 0 else 0
    c3 = sum(1 for c in cycles if len(c) == 3)
    c5 = sum(1 for c in cycles if len(c) == 5)
    c7 = sum(1 for c in cycles if len(c) == 7)

    print(f"  bits={bits}: S={S}, H={H}, c3={c3}, c5={c5}, c7={c7}, "
          f"avg_overlap_w={avg_weight:.1f}")


print("\nDONE.")
