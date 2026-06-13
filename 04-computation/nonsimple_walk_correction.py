#!/usr/bin/env python3
"""
nonsimple_walk_correction.py — The non-simple walk correction that saves Paley

KEY FINDING: The trace-based H approximation H_approx = 1 + sum 2^{(k-1)/2} * tr(A^k)/k
ALWAYS favors the interval, even at p=7 and p=11 where Paley actually wins H.

This means: Paley's actual H advantage at small p comes entirely from the
CORRECTION term c_k - tr(A^k)/k for k >= 7. This correction accounts for
non-simple walks (walks that revisit vertices) and must favor Paley.

Goal: Compute actual cycle counts c_k for both Paley and Interval at p=7, 11
and compare with tr(A^k)/k to find the correction.

Author: kind-pasteur-2026-03-12-S56c
"""

import cmath
import math
import numpy as np
from itertools import permutations


def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p - 1) // 2, p) == 1


def adjacency_matrix(S, p):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S:
                A[i][j] = 1
    return A


def count_directed_cycles(A, k):
    """Count directed k-cycles using Held-Karp DP on each k-subset."""
    n = len(A)
    total = 0
    from itertools import combinations

    for subset in combinations(range(n), k):
        # Count Hamiltonian cycles in the subgraph induced by subset
        vertices = list(subset)
        m = len(vertices)
        if m != k:
            continue

        # Build subgraph adjacency
        sub_A = np.zeros((m, m), dtype=np.int8)
        for i in range(m):
            for j in range(m):
                sub_A[i][j] = A[vertices[i]][vertices[j]]

        # Count Hamiltonian cycles starting from vertex 0
        # Using DP: dp[mask][v] = # paths from 0 visiting exactly mask, ending at v
        full = (1 << m) - 1
        dp = np.zeros((1 << m, m), dtype=np.int64)
        dp[1 << 0][0] = 1

        for mask in range(1, 1 << m):
            for v in range(m):
                if not (mask & (1 << v)):
                    continue
                c = dp[mask][v]
                if c == 0:
                    continue
                for u in range(m):
                    if mask & (1 << u):
                        continue
                    if sub_A[v][u]:
                        dp[mask | (1 << u)][u] += c

        # Count cycles: paths visiting all vertices, ending at v with edge v->0
        cycles = 0
        for v in range(1, m):  # Don't count v=0 (that would be length-0 cycle)
            if sub_A[v][0]:
                cycles += int(dp[full][v])

        total += cycles

    # Each k-cycle is counted k times (once per starting vertex in the cycle)
    # But we fixed vertex 0 as start, so actually each cycle on a fixed k-subset
    # is counted once per rotation... wait.
    # We fix vertex 0 as start. So each directed cycle through vertex 0 is
    # counted exactly once. But the cycle may not go through vertex 0.
    # Actually we enumerate over all k-subsets, and for each we count
    # Hamiltonian cycles starting from the first vertex (index 0 in the subset).
    # A directed k-cycle on this subset will be counted once: when vertex 0
    # in the subset is the lexicographically first vertex of the cycle.
    # Wait no - we fix subset[0] as start, so we count directed cycles
    # that start at subset[0]. Each directed cycle passes through subset[0]
    # exactly once, so it's counted once.
    # But we enumerate ALL k-subsets, and the cycle is in exactly one subset.
    # So total counts each directed k-cycle ONCE (started from its lowest vertex).
    # Wait, no. The cycle has k vertices. We enumerate the k-subset containing
    # all of them. In that subset, we start from index 0 = the first element.
    # The directed cycle visits all k vertices. How many times is it counted?
    # Once: when we pick this specific subset, and the cycle starts at subset[0].
    # Actually no - the cycle doesn't have to start at subset[0].
    # Let me reconsider. dp counts paths starting at local index 0.
    # The cycle completes when dp[full][v] and sub_A[v][0] = 1.
    # So we count directed Hamiltonian cycles that include vertex 0
    # (in local indexing). Since vertex 0 is subset[0], we count each
    # directed cycle that passes through the smallest-indexed vertex.
    # But every Hamiltonian cycle on this subset passes through ALL vertices
    # including the smallest. So we count it once.
    # Since there are k directed rotations of each undirected cycle,
    # and we only count cycles starting at the first vertex... hmm.
    #
    # Actually: a directed cycle on vertices {v1,...,vk} can be written
    # as (v_{sigma(1)}, ..., v_{sigma(k)}) for k rotations. When we
    # fix start = smallest vertex, we get exactly ONE rotation.
    # So total = number of DIRECTED cycles (each counted once).
    # That's correct.

    return total


def count_hp(A):
    """Count Hamiltonian paths using DP."""
    n = len(A)
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = int(dp[mask][v])
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += c
    return int(np.sum(dp[full]))


def main():
    print("=" * 70)
    print("NON-SIMPLE WALK CORRECTION ANALYSIS")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2

        S_paley = frozenset(j for j in range(1, p) if is_qr(j, p))
        S_interval = frozenset(range(1, m + 1))

        A_P = adjacency_matrix(S_paley, p)
        A_I = adjacency_matrix(S_interval, p)

        # Compute actual H
        H_P = count_hp(A_P)
        H_I = count_hp(A_I)

        print(f"\n  p={p}: H(Paley)={H_P}, H(Interval)={H_I}, diff={H_P-H_I}")

        omega = cmath.exp(2j * cmath.pi / p)
        eigs_p = [sum(omega ** (k * s) for s in S_paley) for k in range(p)]
        eigs_i = [sum(omega ** (k * s) for s in S_interval) for k in range(p)]

        print(f"\n    Cycle comparison (c_k actual vs tr(A^k)/k approx):")
        print(f"    {'k':>4} {'c_k(P)':>10} {'tr/k(P)':>10} {'corr(P)':>10} "
              f"{'c_k(I)':>10} {'tr/k(I)':>10} {'corr(I)':>10} "
              f"{'Dc_k':>10} {'Dtr/k':>10} {'Dcorr':>10}")

        # OCF: H = 1 + sum 2^{(k-1)/2} * c_k
        H_P_ocf = 1
        H_I_ocf = 1
        H_P_approx = 1
        H_I_approx = 1

        for k in range(3, p + 1, 2):
            tr_p = sum(e ** k for e in eigs_p).real
            tr_i = sum(e ** k for e in eigs_i).real
            tr_div_k_p = tr_p / k
            tr_div_k_i = tr_i / k

            if k <= p:  # can compute actual cycle count
                c_k_p = count_directed_cycles(A_P, k)
                c_k_i = count_directed_cycles(A_I, k)
            else:
                c_k_p = c_k_i = None

            if c_k_p is not None:
                corr_p = c_k_p - tr_div_k_p
                corr_i = c_k_i - tr_div_k_i
                Dc_k = c_k_p - c_k_i
                Dtr = tr_div_k_p - tr_div_k_i
                Dcorr = corr_p - corr_i

                weight = 2 ** ((k - 1) // 2)
                H_P_ocf += weight * c_k_p
                H_I_ocf += weight * c_k_i
                H_P_approx += weight * tr_div_k_p
                H_I_approx += weight * tr_div_k_i

                print(f"    {k:>4} {c_k_p:>10} {tr_div_k_p:>10.1f} {corr_p:>10.1f} "
                      f"{c_k_i:>10} {tr_div_k_i:>10.1f} {corr_i:>10.1f} "
                      f"{Dc_k:>10} {Dtr:>10.1f} {Dcorr:>10.1f}")

        print(f"\n    H(Paley) OCF = {H_P_ocf}, H(Interval) OCF = {H_I_ocf}")
        print(f"    H(Paley) approx = {H_P_approx:.1f}, H(Interval) approx = {H_I_approx:.1f}")
        print(f"    OCF match: Paley {H_P_ocf == H_P}, Interval {H_I_ocf == H_I}")
        print(f"    Correction swings H by: {(H_P_ocf - H_I_ocf) - (H_P_approx - H_I_approx):.1f}")


if __name__ == '__main__':
    main()
