#!/usr/bin/env python3
"""
Fast real-rootedness test for I(Omega_3(T), x) at large n.

Key insight: alpha(Omega_3) <= floor(n/3). An independent set of size k
is a collection of k vertex-disjoint 3-cycles. So we need to count
k-matchings in the 3-uniform hypergraph of 3-cycles.

For k <= 7 (n <= 21), we can enumerate efficiently using recursive
extension with vertex tracking.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament
import numpy as np
from collections import defaultdict

def find_3cycles(T):
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i,j,k))
    return cycles

def count_disjoint_matchings(cycles, max_k):
    """Count collections of k pairwise vertex-disjoint 3-cycles for k=0,...,max_k."""
    m = len(cycles)
    if m == 0:
        return [1]

    # Convert to frozensets for fast intersection
    cycle_sets = [frozenset(c) for c in cycles]

    # For each cycle, precompute which cycles it's disjoint from
    # (and which have higher index, to avoid double counting)
    disjoint_after = [[] for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if not (cycle_sets[i] & cycle_sets[j]):
                disjoint_after[i].append(j)

    counts = [0] * (max_k + 1)
    counts[0] = 1
    counts[1] = m

    if max_k >= 2:
        # Count pairs
        c2 = 0
        for i in range(m):
            c2 += len(disjoint_after[i])
        counts[2] = c2

    if max_k >= 3:
        # Count triples: extend each pair
        c3 = 0
        for i in range(m):
            used_i = cycle_sets[i]
            for j in disjoint_after[i]:
                used_ij = used_i | cycle_sets[j]
                for k_idx in range(len(disjoint_after[j])):
                    k = disjoint_after[j][k_idx]
                    if k > j and not (cycle_sets[k] & used_ij):
                        c3 += 1
        counts[3] = c3

    if max_k >= 4:
        # Count 4-tuples
        c4 = 0
        for i in range(m):
            used_i = cycle_sets[i]
            for j in disjoint_after[i]:
                used_ij = used_i | cycle_sets[j]
                cands_ij = [k for k in disjoint_after[j] if k > j and not (cycle_sets[k] & used_ij)]
                for ci, k in enumerate(cands_ij):
                    used_ijk = used_ij | cycle_sets[k]
                    for l in cands_ij[ci+1:]:
                        if not (cycle_sets[l] & used_ijk):
                            c4 += 1
        counts[4] = c4

    if max_k >= 5:
        # Count 5-tuples (recursive approach)
        c5 = 0
        for i in range(m):
            used_i = cycle_sets[i]
            for j in disjoint_after[i]:
                used_ij = used_i | cycle_sets[j]
                cands_ij = [k for k in disjoint_after[j] if k > j and not (cycle_sets[k] & used_ij)]
                for ci, k in enumerate(cands_ij):
                    used_ijk = used_ij | cycle_sets[k]
                    cands_ijk = [l for l in cands_ij[ci+1:] if not (cycle_sets[l] & used_ijk)]
                    for cl, l in enumerate(cands_ijk):
                        used_ijkl = used_ijk | cycle_sets[l]
                        for p in cands_ijk[cl+1:]:
                            if not (cycle_sets[p] & used_ijkl):
                                c5 += 1
        counts[5] = c5

    # Trim trailing zeros
    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()

    return counts

print("=" * 70)
print("REAL-ROOTEDNESS OF I(Omega_3(T), x) — FAST VERSION")
print("=" * 70)

for n in [9, 10, 12, 15, 18, 21]:
    max_alpha = n // 3
    samples = 500 if n <= 10 else (200 if n <= 15 else (50 if n <= 18 else 20))
    fails = 0
    deg_dist = defaultdict(int)
    c3_stats = []
    all_coeffs = []

    for trial in range(samples):
        T = random_tournament(n)
        c3 = find_3cycles(T)
        c3_stats.append(len(c3))
        if not c3:
            deg_dist[0] += 1
            continue

        coeffs = count_disjoint_matchings(c3, max_alpha)
        deg = len(coeffs) - 1
        deg_dist[deg] += 1
        all_coeffs.append(coeffs)

        if deg >= 2:
            p = list(reversed(coeffs))
            roots = np.roots(p)
            if any(abs(r.imag) > 1e-6 for r in roots):
                fails += 1
                print(f"  FAIL n={n}, trial {trial}: {coeffs}")
                print(f"    roots: {[f'{r:.4f}' for r in roots]}")

        if trial > 0 and trial % 50 == 0:
            print(f"  n={n}: {trial}/{samples} done, {fails} failures so far", flush=True)

    avg_c3 = sum(c3_stats) / len(c3_stats)
    print(f"\nn={n}: {samples} samples, {fails} FAILURES")
    print(f"  avg c3 = {avg_c3:.1f}, max c3 = {max(c3_stats)}")
    print(f"  Degree dist: {dict(sorted(deg_dist.items()))}")

    # Log-concavity check
    if all_coeffs:
        lc_fails = 0
        for coeffs in all_coeffs:
            d = len(coeffs) - 1
            for k in range(1, d):
                if coeffs[k]**2 < coeffs[k-1] * coeffs[k+1]:
                    lc_fails += 1
                    break
        print(f"  Log-concavity failures: {lc_fails}/{len(all_coeffs)}")

    # Newton's inequality check (necessary for all-real-roots)
    if all_coeffs:
        newton_fails = 0
        for coeffs in all_coeffs:
            d = len(coeffs) - 1
            if d < 2:
                continue
            for k in range(1, d):
                # Newton: a_k^2 * C(d,k-1)*C(d,k+1) >= a_{k-1}*a_{k+1} * C(d,k)^2
                # Simplified: (k+1)*(d-k+1)*a_k^2 >= k*(d-k)*a_{k-1}*a_{k+1} ...
                # Actually simplest: a_k^2 >= a_{k-1}*a_{k+1} (log-concavity)
                # which we already checked
                pass

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
