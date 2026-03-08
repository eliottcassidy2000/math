#!/usr/bin/env python3
"""
Ihara Zeta Function of Tournaments (kind-pasteur-S34).

The Ihara-Bass formula for a directed graph D on n vertices:
  zeta_D(u)^{-1} = det(I - A*u + (D_out - I)*u^2)
where A = adjacency matrix, D_out = diag(out-degrees).

For a tournament T on n vertices:
  - A is the {0,1} adjacency matrix with A[i][j] = 1 iff i->j
  - D_out = diag(s_0, s_1, ..., s_{n-1}) where s_i = out-degree of vertex i
  - Note: A + A^T = J - I (all-ones minus identity)

The Ihara zeta encodes ALL directed cycle counts.

We test: does zeta_T(u) at specific u values relate to I(Omega(T), x) or H(T)?

Key evaluations to try:
  - zeta_T(1/2): since OCF uses x=2 and 2^{-1} = 1/2
  - zeta_T(1/3): since I(Omega, 3) is computed
  - log derivative: -zeta'/zeta gives cycle counts
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

import numpy as np
from itertools import combinations
from collections import defaultdict


def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def count_ham_paths(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])


def ihara_bass_det(adj, n, u):
    """Compute det(I - A*u + (D_out - I)*u^2) for the Ihara-Bass formula."""
    A = np.array(adj, dtype=float)
    D_out = np.diag(A.sum(axis=1))
    I = np.eye(n)
    M = I - A * u + (D_out - I) * (u ** 2)
    return np.linalg.det(M)


def count_3cycles(adj, n):
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or \
           (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1
    return c3


def main():
    print("=== Ihara Zeta Function of Tournaments ===\n")

    for n in [4, 5]:
        print(f"--- n={n} ---")
        num_edges = n * (n - 1) // 2
        max_bits = 1 << num_edges

        h_to_zeta = defaultdict(list)
        results = []

        for bits in range(max_bits):
            adj = tournament_adj(n, bits)
            H = count_ham_paths(adj, n)
            c3 = count_3cycles(adj, n)

            # Compute Ihara-Bass determinant at various u values
            z_half = ihara_bass_det(adj, n, 0.5)
            z_third = ihara_bass_det(adj, n, 1/3)
            z_quarter = ihara_bass_det(adj, n, 0.25)

            # zeta(u)^{-1} = det(I - Au + (D-I)u^2)
            # So zeta(u) = 1 / det(...)
            # We want zeta_T(1/2) etc.

            results.append((H, c3, z_half, z_third, z_quarter, bits))
            h_to_zeta[H].append((z_half, z_third))

        # Print sample results
        print(f"\n  {'H':>4} {'c3':>3} {'z_inv(1/2)':>12} {'z_inv(1/3)':>12} {'z_inv(1/4)':>12}")
        seen = set()
        for H, c3, zh, zt, zq, bits in sorted(results):
            if (H, c3) not in seen:
                seen.add((H, c3))
                print(f"  {H:4d} {c3:3d} {zh:12.4f} {zt:12.4f} {zq:12.4f}")

        # Check: is z_inv(1/2) determined by H?
        print(f"\n  Is z_inv(1/2) determined by H?")
        for H in sorted(h_to_zeta):
            vals = set(round(v[0], 6) for v in h_to_zeta[H])
            det = "YES" if len(vals) == 1 else f"NO ({len(vals)} values)"
            print(f"    H={H:3d}: {det} {sorted(vals)[:5]}")

        # Check: is z_inv(1/3) determined by H?
        print(f"\n  Is z_inv(1/3) determined by H?")
        for H in sorted(h_to_zeta):
            vals = set(round(v[1], 6) for v in h_to_zeta[H])
            det = "YES" if len(vals) == 1 else f"NO ({len(vals)} values)"
            print(f"    H={H:3d}: {det} {sorted(vals)[:5]}")

    # Part 2: Relationship between zeta and H
    print(f"\n\n=== Part 2: Regression analysis ===")
    n = 5
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    H_vals = []
    z_vals = []
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        z = ihara_bass_det(adj, n, 0.5)
        H_vals.append(H)
        z_vals.append(z)

    H_arr = np.array(H_vals, dtype=float)
    z_arr = np.array(z_vals, dtype=float)

    # Check correlation
    corr = np.corrcoef(H_arr, z_arr)[0, 1]
    print(f"  n=5: Pearson correlation H vs z_inv(1/2): {corr:.6f}")

    # Check if log(z_inv(1/2)) ~ a * log(H) + b
    log_H = np.log(H_arr)
    log_z = np.log(np.abs(z_arr))
    corr_log = np.corrcoef(log_H, log_z)[0, 1]
    print(f"  n=5: Pearson correlation log|H| vs log|z_inv(1/2)|: {corr_log:.6f}")

    # Part 3: Characteristic polynomial of A at special values
    print(f"\n\n=== Part 3: Adjacency eigenvalues ===")
    n = 5
    num_edges = n * (n - 1) // 2
    max_bits = 1 << num_edges

    h_to_spectrum = defaultdict(list)
    for bits in range(max_bits):
        adj = tournament_adj(n, bits)
        H = count_ham_paths(adj, n)
        A = np.array(adj, dtype=complex)
        eigs = np.sort(np.linalg.eigvals(A))
        # Tournament eigenvalues: one real = (n-1)/2, rest in conjugate pairs
        real_eig = np.real(eigs)
        h_to_spectrum[H].append(tuple(np.round(real_eig, 4)))

    for H in sorted(h_to_spectrum):
        specs = set(h_to_spectrum[H])
        print(f"  H={H:3d}: {len(specs)} distinct real-part spectra")


if __name__ == "__main__":
    main()
