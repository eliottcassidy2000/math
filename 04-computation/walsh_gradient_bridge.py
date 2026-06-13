#!/usr/bin/env python3
"""
walsh_gradient_bridge.py -- Connecting Walsh decomposition to co-occurrence gradient

CONJECTURE: Walsh degree-2d in the orientation cube H(sigma) corresponds to
alpha_2 contributions from cycles of length ~(2d+1).

At Walsh degree 2: H depends on pairs of chords (sigma_i * sigma_j).
  These correspond to 3-cycle interactions (each 3-cycle involves 2 chords).
  The co-occurrence gradient at k=3 has slope b_3=1.

At Walsh degree 4: H depends on 4-tuples of chords.
  These correspond to 5-cycle interactions.
  The co-occurrence gradient at k=5 has slope b_5=C(m-2,2).

The SCALING MATCH:
  Walsh deg-2 energy ~ O(m^2) per HYP-529
  Co-occ gradient effect at k=3 ~ O(p^2) ~ O(m^2)

  Walsh deg-4 energy ~ O(m^{2+8.98}) per HYP-496
  Co-occ gradient effect at k=5 ~ O(p^6) ~ O(m^6)

This script tests the correspondence numerically at p=7,11,13.

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
from math import comb, factorial
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def all_circulant_H(p):
    """Compute H for all 2^m circulant tournaments."""
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]

    results = {}
    for bits in range(1 << m):
        sigma = []
        S = []
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.append(a)
                sigma.append(+1)
            else:
                S.append(b)
                sigma.append(-1)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        results[tuple(sigma)] = H

    return results


def walsh_decomposition(H_dict, m):
    """Compute Walsh-Hadamard transform of H on {+1,-1}^m.

    H(sigma) = sum_S h_hat[S] * prod_{i in S} sigma_i
    h_hat[S] = (1/2^m) * sum_sigma H(sigma) * prod_{i in S} sigma_i
    """
    n = 1 << m
    h_hat = {}

    for bits in range(n):
        S = [i for i in range(m) if bits & (1 << i)]
        # Compute h_hat[S]
        coeff = 0
        for sigma_bits in range(n):
            sigma = [(1 if sigma_bits & (1 << i) else -1) for i in range(m)]
            H = H_dict[tuple(sigma)]
            prod = 1
            for i in S:
                prod *= sigma[i]
            coeff += H * prod
        h_hat[tuple(S)] = coeff / n
        # Only keep degree = |S|

    return h_hat


def walsh_by_degree(h_hat, m):
    """Group Walsh coefficients by degree."""
    by_degree = defaultdict(list)
    for S, coeff in h_hat.items():
        deg = len(S)
        by_degree[deg].append((S, coeff))

    # Energy at each degree = sum of coeff^2
    energy_by_deg = {}
    for deg, items in by_degree.items():
        energy = sum(c**2 for _, c in items)
        energy_by_deg[deg] = energy

    return by_degree, energy_by_deg


def co_occ_formula(p, k, d):
    m = (p - 1) // 2
    if d > m:
        d = p - d
    return comb(2*m - 1, k - 2) - comb(m - 1, k - 2) - (m - d) * comb(m - 2, k - 3)


def main():
    print("=" * 70)
    print("WALSH-GRADIENT BRIDGE")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        # Compute all H values
        H_dict = all_circulant_H(p)
        n_tours = len(H_dict)
        print(f"  {n_tours} circulant tournaments")

        # Walsh decomposition
        h_hat = walsh_decomposition(H_dict, m)

        # Group by degree
        by_deg, energy_by_deg = walsh_by_degree(h_hat, m)

        total_energy = sum(energy_by_deg.values())
        print(f"\n  Walsh decomposition:")
        print(f"    {'deg':>4} {'#coeffs':>8} {'energy':>18} {'pct':>8}")
        for deg in sorted(energy_by_deg):
            e = energy_by_deg[deg]
            pct = 100 * e / total_energy if total_energy > 0 else 0
            n_coeffs = len(by_deg[deg])
            print(f"    {deg:>4} {n_coeffs:>8} {e:>18.1f} {pct:>8.2f}%")

        # Show h_hat[0] = mean H
        mean_H = h_hat[()]
        print(f"\n    h_hat[{{}}] = mean H = {mean_H:.1f}")

        # Degree-0 energy = mean_H^2
        print(f"    Degree-0 energy = {mean_H**2:.1f}")

        # Show all non-zero coefficients
        print(f"\n    Non-zero Walsh coefficients:")
        for deg in sorted(by_deg):
            items = sorted(by_deg[deg], key=lambda x: -abs(x[1]))
            for S, c in items:
                if abs(c) > 0.01:
                    print(f"      deg={deg}, S={set(S)}: h_hat = {c:.2f}")

        # ====== INTERVAL VALUE ======
        # Interval = all sigma_i = +1
        sigma_int = tuple([1] * m)
        H_int = H_dict[sigma_int]

        # Interval Walsh: H(1,...,1) = sum_S h_hat[S] (all products = 1)
        H_walsh_check = sum(h_hat.values())
        print(f"\n    H(Interval) = {H_int}")
        print(f"    Walsh sum check = {H_walsh_check:.1f}")

        # How much does each degree contribute to Interval's H?
        print(f"\n    Interval H by Walsh degree:")
        for deg in sorted(by_deg):
            contrib = sum(c for _, c in by_deg[deg])
            print(f"      deg {deg}: {contrib:.1f} ({100*contrib/H_int:.2f}%)")

        # ====== CO-OCCURRENCE GRADIENT ======
        print(f"\n    Co-occurrence gradient (THM-143):")
        for k in range(3, min(p + 1, 14), 2):
            b_k = comb(m - 2, k - 3)
            var_proxy = p * b_k**2 * (m**2 - 1) / 12
            print(f"      k={k}: b_k={b_k}, var_proxy={var_proxy:.0f}")

        # ====== THE BRIDGE ======
        # Walsh degree 2 involves C(m,2) pairwise products sigma_i*sigma_j
        # Each (sigma_i, sigma_j) pair corresponds to chord pair (i,j)
        # A 3-cycle uses 3 chords (gaps g1,g2,g3 with g1+g2+g3=p)
        # The Walsh degree-2 coefficient h_hat[{i,j}] depends on
        # how many 3-cycles use both chords i and j.

        # Walsh degree 4 involves C(m,4) quadruple products
        # A 5-cycle uses 5 chords, and the degree-4 coefficient
        # depends on how many 5-cycles use a specific 4-chord subset.

        if p <= 11:
            print(f"\n    BRIDGE: Walsh coeff -> cycle content")

            # For degree 2: each h_hat[{i,j}] measures the effect of
            # flipping chord i and chord j simultaneously.
            # In terms of 3-cycles: flipping chords i,j affects all
            # 3-cycles that use BOTH gaps corresponding to i and j.

            # The "gap" corresponding to chord i is the i-th pair
            # where pair i = (i+1, p-i-1).
            print(f"\n    Chord-pair gap mapping:")
            for i in range(m):
                a, b = i + 1, p - i - 1
                print(f"      chord {i}: gap {a} or {b}")

    # ====== RATIO COMPARISON ======
    print(f"\n{'='*70}")
    print("RATIO COMPARISON: Walsh energy vs co-occ gradient")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        H_dict = all_circulant_H(p)
        h_hat = walsh_decomposition(H_dict, m)
        _, energy = walsh_by_degree(h_hat, m)

        # Walsh energy ratios
        e2 = energy.get(2, 0)
        e4 = energy.get(4, 0)
        e6 = energy.get(6, 0)
        ratio_42 = e4 / e2 if e2 > 0 else float('inf')
        ratio_62 = e6 / e2 if e2 > 0 else float('inf')

        # Co-occ gradient ratios
        b3 = comb(m - 2, 0)  # =1
        b5 = comb(m - 2, 2)
        b7 = comb(m - 2, 4)
        # Variance proxy ratio
        vr_53 = (b5 / b3)**2 if b3 > 0 else 0
        vr_73 = (b7 / b3)**2 if b3 > 0 else 0

        print(f"\n  p={p}, m={m}:")
        print(f"    Walsh: E_4/E_2 = {ratio_42:.4f}, E_6/E_2 = {ratio_62:.4f}")
        print(f"    Co-occ: (b_5/b_3)^2 = {vr_53:.4f}, (b_7/b_3)^2 = {vr_73:.4f}")
        print(f"    Match quality: E_4/E_2 vs b_5^2 = {ratio_42:.4f} vs {b5**2}")


if __name__ == '__main__':
    main()
