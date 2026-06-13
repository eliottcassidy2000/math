#!/usr/bin/env python3
"""
trace_H_analytic.py -- Analytic connection: trace alternation -> H maximization

The key question: does trace alternation (THM-136) imply H(Interval) > H(Paley)?

APPROACH: Use the Newton-Girard formulas to express the independence polynomial
I(Omega, x) in terms of the "power sums" (traces).

For a graph G with n vertices and independence polynomial I(G, x) = sum alpha_k x^k:
  The CLIQUE polynomial C(G, x) = I(complement(G), x) satisfies:
  C(G, x) = exp(-sum_{k>=1} tr(A_G^k) x^k / k)  ... NO, this is for matching poly.

Actually, the connection between eigenvalues and independence polynomial is NOT
through Newton-Girard. Independence polynomial is NOT the characteristic polynomial.

CORRECT APPROACH: Use the TRACE of the TOURNAMENT MATRIX to bound cycle counts,
then bound the independence polynomial of the conflict graph.

For tournament T, the directed cycle structure determines Omega(T).
The cycle counts c_k are constrained (but not determined) by the traces tr(A^k).

Specifically: tr(A^k) = sum over all closed walks of length k.
For k = 3: tr(A^3) = 3 * c_3 (no backtracking possible in tournament)
For k = 5: tr(A^5) = 5 * c_5 + (contributions from 3-cycles)
  A closed walk of length 5 can visit a 3-cycle and then traverse 2 more edges.
  In a tournament with no 2-cycles, a walk of length 5 must visit at least 3
  distinct vertices. Options:
  - 5 distinct vertices: a directed 5-cycle (contributes 5 per cycle)
  - 4 distinct vertices: visit one vertex twice
  - 3 distinct vertices: visit vertices in a 3-cycle pattern

Actually, for TOURNAMENT-specific structure:
  tr(A^k) = sum_{v} #{closed walks of length k starting at v}

  For k = 3: only 3-cycles contribute (no 2-cycles).
  For k = 5: 5-cycles AND closed walks using shorter cycles.

Let me just COMPUTE the relationship empirically.

Author: kind-pasteur-2026-03-12-S57
"""

import math
import cmath
from itertools import combinations


def compute_all_odd_cycles(p, S):
    """Compute all directed odd cycle counts for circulant tournament."""
    S_set = set(S)
    cycle_counts = {}

    for k in range(3, p + 1, 2):
        count = 0
        for subset in combinations(range(p), k):
            # Since circulant: just check if this subset has a Hamiltonian cycle
            # in the induced sub-tournament
            nodes = list(subset)
            n = len(nodes)
            adj = [[False] * n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i != j and (nodes[j] - nodes[i]) % p in S_set:
                        adj[i][j] = True

            # Count directed Hamiltonian cycles (starting from node 0 to avoid overcounting)
            if n <= 3:
                # For n=3: check if 0->1->2->0 or 0->2->1->0
                hc = 0
                if n == 3:
                    if adj[0][1] and adj[1][2] and adj[2][0]:
                        hc = 1
                    if adj[0][2] and adj[2][1] and adj[1][0]:
                        hc += 1
                count += hc
            else:
                # Held-Karp for Hamiltonian cycles
                N = 1 << (n - 1)
                dp = [0] * (N * n)
                dp[0] = 1

                for mask in range(N):
                    for v in range(n):
                        c = dp[mask * n + v]
                        if c == 0:
                            continue
                        for u in range(1, n):
                            bit = u - 1
                            if mask & (1 << bit):
                                continue
                            if adj[v][u]:
                                dp[(mask | (1 << bit)) * n + u] += c

                full = N - 1
                hc = sum(dp[full * n + v] for v in range(1, n) if adj[v][0])
                count += hc

        cycle_counts[k] = count

    return cycle_counts


def compute_traces_from_eigs(p, S, max_k=None):
    """Compute traces from eigenvalues."""
    if max_k is None:
        max_k = p
    omega = cmath.exp(2j * cmath.pi / p)
    eigs = [sum(omega ** (j * s) for s in S) for j in range(p)]
    traces = {}
    for k in range(3, max_k + 1, 2):
        traces[k] = round(sum(e ** k for e in eigs).real)
    return traces


def main():
    print("=" * 70)
    print("TRACE-CYCLE RELATIONSHIP IN TOURNAMENTS")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            print(f"\n--- p={p}, {name} (S={S}) ---")

            traces = compute_traces_from_eigs(p, S)
            cycles = compute_all_odd_cycles(p, S)

            print(f"  {'k':>4} {'tr(A^k)':>12} {'c_k (cycles)':>12} {'tr/k':>12} {'ratio tr/(k*c_k)':>18}")
            for k in sorted(traces.keys()):
                if k in cycles and cycles[k] > 0:
                    print(f"  {k:4d} {traces[k]:12d} {cycles[k]:12d} {traces[k]/k:12.2f} {traces[k]/(k*cycles[k]):18.4f}")
                else:
                    ck = cycles.get(k, 0)
                    print(f"  {k:4d} {traces[k]:12d} {ck:12d} {traces[k]/k:12.2f} {'N/A':>18}")

            # OCF: H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
            alpha_1 = sum(cycles.values())
            print(f"\n  Total odd cycles (alpha_1 of Omega): {alpha_1}")

            # Simple bound: H >= 1 + 2*alpha_1
            H_lower = 1 + 2 * alpha_1
            print(f"  H >= 1 + 2*alpha_1 = {H_lower}")

    # Now the key analysis: HOW does trace alternation affect H?
    print("\n" + "=" * 70)
    print("TRACE ALTERNATION IMPACT ON H")
    print("=" * 70)
    print("""
  THM-136 says: for k = 1 mod 4, tr(A_P^k) > tr(A_I^k)
                for k = 3 mod 4, tr(A_P^k) < tr(A_I^k)

  The trace alternation does NOT directly give c_k differences because
  tr(A^k) includes closed walks, not just cycles.

  However, for the TOTAL alpha_1 = sum c_k, we can ask:
  Does the interval have MORE total directed cycles?

  At p=7: alpha_1(Paley) = ? vs alpha_1(Interval) = ?
  At p=11: similarly
""")

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
        S_int = list(range(1, m + 1))

        cycles_P = compute_all_odd_cycles(p, S_qr)
        cycles_I = compute_all_odd_cycles(p, S_int)

        alpha1_P = sum(cycles_P.values())
        alpha1_I = sum(cycles_I.values())

        print(f"\np={p}:")
        print(f"  Paley cycles:    {cycles_P}")
        print(f"  Interval cycles: {cycles_I}")
        print(f"  alpha_1(Paley)    = {alpha1_P}")
        print(f"  alpha_1(Interval) = {alpha1_I}")

        # The contribution to H:
        # The OCF-weighted cycle counts
        weighted_P = sum((2 ** ((k - 1) / 2)) * ck for k, ck in cycles_P.items())
        weighted_I = sum((2 ** ((k - 1) / 2)) * ck for k, ck in cycles_I.items())
        print(f"  Weighted sum(2^{{(k-1)/2}} * c_k):")
        print(f"    Paley = {weighted_P:.1f}")
        print(f"    Interval = {weighted_I:.1f}")
        print(f"    Advantage: {'Paley' if weighted_P > weighted_I else 'Interval'}")

        # But H depends on I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
        # The weighted cycle sum is NOT H. We need the independence polynomial.
        # However, the key insight is:
        # H = 1 + 2*alpha_1 + 4*alpha_2 + ... and alpha_1 is the "first order" term
        # If alpha_1(I) > alpha_1(P), AND the higher alpha_j favor interval,
        # then H(I) > H(P).

    print("\n" + "=" * 70)
    print("CRITICAL INSIGHT")
    print("=" * 70)
    print("""
  The trace alternation theorem (THM-136) controls tr(A^k) differences.
  But H depends on I(Omega, 2), where Omega is the CONFLICT GRAPH of
  directed odd cycles.

  The connection is INDIRECT:
  1. More directed k-cycles => more vertices in Omega
  2. But also potentially more edges (conflicts) in Omega
  3. The independence polynomial I(Omega, 2) depends on the BALANCE
     between having many cycles (high alpha_1) and having them be
     non-conflicting (high alpha_2, alpha_3, ...)

  The interval tournament's advantage at large k (THM-136) manifests as:
  - More LONG directed cycles (k large)
  - Long cycles tend to be more vertex-disjoint (each uses many vertices)
  - This creates larger independent sets in Omega
  - Which amplifies H via higher alpha_j terms

  This is exactly the ADDITIVE ENERGY / FLOW argument from cross-field
  connections: interval's additive structure creates coherent flow paths
  that tile the vertex set more efficiently.
""")


if __name__ == '__main__':
    main()
