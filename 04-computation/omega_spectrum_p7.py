#!/usr/bin/env python3
"""
omega_spectrum_p7.py -- Spectrum of the conflict graph Omega at p=7

For each of the 2 orbit classes at p=7, compute the adjacency spectrum
of Omega(T) and see how it relates to the Walsh structure.

At p=7, all cycles are 3-cycles, 5-cycles, and 7-cycles.
The conflict graph has: alpha_1 vertices, with edges between overlapping cycles.

Since H = I(Omega, 2) and we want to understand how H varies with orientation,
the spectrum of Omega is the key.

For claw-free graphs (our case at p<=8), I(G,x) has all real roots.
The roots determine the alpha_j via Vieta's formulas.

Author: kind-pasteur-2026-03-12-S60
"""

from itertools import combinations
from collections import defaultdict
import math


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
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


def main():
    p = 7
    m = (p - 1) // 2
    half = 1 << m
    pairs = [(s, p - s) for s in range(1, m + 1)]

    print("=" * 70)
    print(f"CONFLICT GRAPH OMEGA SPECTRUM at p={p}")
    print("=" * 70)

    for bits in [0b000, 0b011]:  # One from each orbit: Interval and Paley-orbit
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                   for i in range(m))

        A = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in S:
                A[v][(v + s) % p] = 1

        # Enumerate all directed odd cycles
        all_cycles = []
        by_len = defaultdict(list)
        for k in range(3, p + 1, 2):
            for subset in combinations(range(p), k):
                verts = list(subset)
                n_cyc = count_ham_cycles(A, verts)
                for _ in range(n_cyc):
                    fs = frozenset(subset)
                    all_cycles.append(fs)
                    by_len[k].append(fs)

        n_cycles = len(all_cycles)

        # Build conflict graph adjacency matrix
        adj = [[0]*n_cycles for _ in range(n_cycles)]
        for i in range(n_cycles):
            for j in range(i+1, n_cycles):
                if all_cycles[i] & all_cycles[j]:  # overlapping vertex sets
                    adj[i][j] = 1
                    adj[j][i] = 1

        # Compute degree sequence
        degrees = [sum(row) for row in adj]
        deg_dist = defaultdict(int)
        for d in degrees:
            deg_dist[d] += 1

        # Compute clique number (max clique in Omega = max overlap set)
        # Alpha_j = independent set count
        alpha_1 = n_cycles
        alpha_2 = sum(1 for i in range(n_cycles)
                       for j in range(i+1, n_cycles)
                       if not (all_cycles[i] & all_cycles[j]))
        H = 1 + 2*alpha_1 + 4*alpha_2

        # Edge count
        n_edges = sum(degrees) // 2

        # For the spectrum: need numpy
        try:
            import numpy as np
            A_np = np.array(adj, dtype=float)
            eigenvalues = np.linalg.eigvalsh(A_np)
            eigenvalues = sorted(eigenvalues, reverse=True)

            # Independence polynomial roots from eigenvalues?
            # I(G,x) = sum alpha_j x^j
            # For the complement G_bar: clique polynomial C(G_bar, x)
            # Lovasz theta function: theta(G) = max(sum x_i : sum x_i*x_j*adj_ij = 0, sum x_i^2 = 1)

            label = "Interval" if bits == 0 else "Paley-orbit"
            print(f"\n  {label} (bits={bits:03b}, S={S}):")
            print(f"    |V(Omega)| = {n_cycles}")
            print(f"    |E(Omega)| = {n_edges}")
            print(f"    by_len: {dict((k, len(v)) for k, v in sorted(by_len.items()))}")
            print(f"    alpha_1={alpha_1}, alpha_2={alpha_2}, H={H}")
            print(f"    Degree dist: {dict(sorted(deg_dist.items()))}")
            print(f"    Mean degree: {sum(degrees)/n_cycles:.4f}")

            print(f"    Top 5 eigenvalues: {[f'{e:.4f}' for e in eigenvalues[:5]]}")
            print(f"    Bottom 5 eigenvalues: {[f'{e:.4f}' for e in eigenvalues[-5:]]}")
            print(f"    Spectral gap: {eigenvalues[0] - eigenvalues[1]:.4f}")

            # Independence polynomial: I(x) = det(I + x*A_complement)^{-1}... no
            # Actually I(G,x) = sum_S x^|S| where S is independent set
            # We already have alpha_0=1, alpha_1, alpha_2
            # For p=7, alpha_j = 0 for j >= 3 (max independent set size 2)
            # So I(G,x) = 1 + alpha_1*x + alpha_2*x^2

            disc = alpha_1**2 - 4*alpha_2
            if disc >= 0:
                r1 = (-alpha_1 + math.sqrt(disc)) / (2*alpha_2) if alpha_2 > 0 else None
                r2 = (-alpha_1 - math.sqrt(disc)) / (2*alpha_2) if alpha_2 > 0 else None
                print(f"    I(x) = 1 + {alpha_1}x + {alpha_2}x^2")
                print(f"    Discriminant: {disc}")
                if r1 is not None:
                    print(f"    Roots: {r1:.6f}, {r2:.6f}")
                    print(f"    H = I(2) = 1 + {2*alpha_1} + {4*alpha_2} = {H}")
            else:
                print(f"    I(x) = 1 + {alpha_1}x + {alpha_2}x^2")
                print(f"    Discriminant: {disc} (complex roots)")

            # Chromatic-like invariant: what's the complement graph?
            n_comp_edges = n_cycles * (n_cycles - 1) // 2 - n_edges
            print(f"    |E(Omega_complement)| = {n_comp_edges}")
            print(f"    Edge density Omega: {2*n_edges/(n_cycles*(n_cycles-1)):.4f}")

            # Compare the two orbits' spectral gap
        except ImportError:
            print("  numpy not available for eigenvalue computation")

    # Now compare the TWO orbits' Omega graphs
    print(f"\n\n{'='*70}")
    print("COMPARISON: INTERVAL vs PALEY Omega")
    print("=" * 70)

    print("""
  Key differences:
  - Paley: alpha_1=80, alpha_2=7  -> fewer cycles but denser conflict graph
  - Interval: alpha_1=59, alpha_2=14 -> fewer cycles but more independent pairs

  The Paley Omega has MORE vertices (80 vs 59) but FEWER independent pairs (7 vs 14).
  This means Paley's conflict graph is DENSER: more overlap between cycles.

  This is the BIBD property in action: Paley arranges 3-cycles so they overlap
  more uniformly, leaving fewer opportunities for disjoint packing.
  But it has MORE cycles overall.

  At H = 1 + 2*alpha_1 + 4*alpha_2:
  Paley: 1 + 160 + 28 = 189
  Interval: 1 + 118 + 56 = 175
  """)

    print("DONE.")


if __name__ == '__main__':
    main()
