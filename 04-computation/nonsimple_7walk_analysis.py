#!/usr/bin/env python3
"""
nonsimple_7walk_analysis.py -- Exact non-simple 7-walk correction and its Walsh structure

At k=7: tr(A^7) = 7*c_7 + N_7 where N_7 counts non-simple closed 7-walks.
These decompose into 3-cycle + 4-cycle sharing a vertex.

Key questions:
1. Is N_7 a linear function of c_4? (since c_3 is constant)
2. What is the Walsh structure of N_7?
3. Does the correction explain why c_7 has anti-product law for q=3?

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def held_karp(A, verts):
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
    print("=" * 70)
    print("NON-SIMPLE 7-WALK CORRECTION")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(j, p - j) for j in range(1, m + 1)]
        n_orient = 1 << m
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n{'='*60}")
        print(f"p={p}, m={m}")
        print("=" * 60)

        data = []  # list of dicts with all quantities per orientation

        for bits in range(n_orient):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))
            A = [[0]*p for _ in range(p)]
            for v in range(p):
                for s in S:
                    A[v][(v + s) % p] = 1

            c3 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 3))
            c4 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 4))
            c5 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 5))
            c7 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 7))

            tr7 = sum(sum(omega**(s*t) for s in S)**7 for t in range(p)).real
            N7 = round(tr7 - 7 * c7)

            data.append({'bits': bits, 'c3': c3, 'c4': c4, 'c5': c5,
                         'c7': c7, 'tr7': tr7, 'N7': N7})

        # Display
        c3_const = data[0]['c3']
        print(f"  c3 = {c3_const} (constant)")
        print(f"  c4 values: {sorted(set(d['c4'] for d in data))}")
        print(f"  c5 values: {sorted(set(d['c5'] for d in data))}")
        print(f"  c7 values: {sorted(set(d['c7'] for d in data))}")
        print(f"  N7 values: {sorted(set(d['N7'] for d in data))}")

        # Test: N7 = a * c4 + b * c5 + c?
        # First check if N7 is affine in c4 alone
        c4_vals = [d['c4'] for d in data]
        N7_vals = [d['N7'] for d in data]
        n = n_orient
        mean_c4 = sum(c4_vals) / n
        mean_N7 = sum(N7_vals) / n
        cov_c4_N7 = sum((c4_vals[i]-mean_c4)*(N7_vals[i]-mean_N7) for i in range(n)) / n
        var_c4 = sum((c4_vals[i]-mean_c4)**2 for i in range(n)) / n

        if var_c4 > 0:
            a_c4 = cov_c4_N7 / var_c4
            b_c4 = mean_N7 - a_c4 * mean_c4
            resid_c4 = [N7_vals[i] - a_c4*c4_vals[i] - b_c4 for i in range(n)]
            max_r_c4 = max(abs(r) for r in resid_c4)
            print(f"\n  Fit N7 = {a_c4:.4f}*c4 + {b_c4:.4f}, max_residual={max_r_c4:.4f}")

            # Express coefficient as fraction of c3
            ratio_a = a_c4 / c3_const if c3_const > 0 else 0
            print(f"    a / c3 = {ratio_a:.6f}")
            # Check: is ratio_a = 12/p?
            print(f"    12/p = {12/p:.6f}")
            # Or: a = 12*c3/p?
            pred_a = 12 * c3_const / p
            print(f"    12*c3/p = {pred_a:.4f} (actual a = {a_c4:.4f})")

        # Bivariate: N7 = a*c4 + b*c5 + c?
        c5_vals = [d['c5'] for d in data]
        # Use least squares: N7 = a*c4 + b*c5 + c
        # Normal equations:
        # [sum c4^2, sum c4*c5, sum c4] [a]   [sum c4*N7]
        # [sum c4*c5, sum c5^2, sum c5] [b] = [sum c5*N7]
        # [sum c4, sum c5, n]           [c]   [sum N7]
        S_c4c4 = sum(c4_vals[i]**2 for i in range(n))
        S_c5c5 = sum(c5_vals[i]**2 for i in range(n))
        S_c4c5 = sum(c4_vals[i]*c5_vals[i] for i in range(n))
        S_c4 = sum(c4_vals)
        S_c5 = sum(c5_vals)
        S_c4N7 = sum(c4_vals[i]*N7_vals[i] for i in range(n))
        S_c5N7 = sum(c5_vals[i]*N7_vals[i] for i in range(n))
        S_N7 = sum(N7_vals)

        # Solve 3x3 system
        import numpy as np
        M = np.array([[S_c4c4, S_c4c5, S_c4],
                       [S_c4c5, S_c5c5, S_c5],
                       [S_c4, S_c5, n]], dtype=float)
        rhs = np.array([S_c4N7, S_c5N7, S_N7], dtype=float)
        try:
            coeffs = np.linalg.solve(M, rhs)
            a_bi, b_bi, c_bi = coeffs
            resid_bi = [N7_vals[i] - a_bi*c4_vals[i] - b_bi*c5_vals[i] - c_bi
                        for i in range(n)]
            max_r_bi = max(abs(r) for r in resid_bi)
            print(f"\n  Bivariate fit: N7 = {a_bi:.4f}*c4 + {b_bi:.4f}*c5 + {c_bi:.4f}")
            print(f"    max_residual = {max_r_bi:.4f}")
            if max_r_bi < 0.01:
                print(f"    *** EXACT FIT! N7 is a LINEAR function of c4 and c5 ***")
                # Express as fractions
                for name, val in [("a", a_bi), ("b", b_bi), ("c", c_bi)]:
                    for denom in range(1, 33):
                        num = val * denom
                        if abs(num - round(num)) < 0.01:
                            print(f"    {name} = {round(num)}/{denom}")
                            break
        except np.linalg.LinAlgError:
            print("  Singular matrix in bivariate fit")

        # Walsh decomposition
        sigma = {}
        for bits in range(n_orient):
            sigma[bits] = tuple(1 if bits & (1 << i) else -1 for i in range(m))
        chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

        print(f"\n  Degree-2 Walsh coefficients:")
        print(f"  {'pair':>8} {'q':>3} {'h_c4':>10} {'h_c5':>10} {'h_c7':>10} {'h_N7/7':>10} {'h_tr/7':>10}")

        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)
            chi_ab = legendre(ga * gb, p)

            h_c4 = sum(data[bits]['c4'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_c5 = sum(data[bits]['c5'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_c7 = sum(data[bits]['c7'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_N7 = sum(data[bits]['N7'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_tr = sum(data[bits]['tr7']/7 * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient

            if abs(h_c7) > 0.01 or abs(h_c4) > 0.01:
                print(f"  ({a},{b}) q={q:>3}: {h_c4:>10.4f} {h_c5:>10.4f} "
                      f"{h_c7:>10.4f} {h_N7/7:>10.4f} {h_tr:>10.4f}")

        # Check relationship: h_c7 = h_tr/7 - h_N7/7
        # h_N7 should relate to h_c4 and h_c5 via the bivariate formula
        if max_r_bi < 0.1 if 'max_r_bi' in dir() else False:
            print(f"\n  Predicted Walsh from bivariate fit:")
            for a, b in chord_pairs:
                h_c4 = sum(data[bits]['c4'] * sigma[bits][a] * sigma[bits][b]
                           for bits in range(n_orient)) / n_orient
                h_c5 = sum(data[bits]['c5'] * sigma[bits][a] * sigma[bits][b]
                           for bits in range(n_orient)) / n_orient
                h_N7 = sum(data[bits]['N7'] * sigma[bits][a] * sigma[bits][b]
                           for bits in range(n_orient)) / n_orient

                h_N7_pred = a_bi * h_c4 + b_bi * h_c5
                if abs(h_N7) > 0.01:
                    q = resonance_level(a+1, b+1, p)
                    print(f"    ({a},{b}) q={q}: h_N7={h_N7:.4f}, "
                          f"pred={h_N7_pred:.4f}, match={abs(h_N7-h_N7_pred)<0.01}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
