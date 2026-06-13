#!/usr/bin/env python3
"""
H_linear_ck_structure.py -- kind-pasteur-2026-03-13-S60

Analyze the structure of H = sum w_k * c_k + C for circulant tournaments.

Key questions:
1. What are the EXACT rational coefficients w_k?
2. Do they have a pattern across primes?
3. Can we verify the formula at p=17 (where H is computable but slow)?
4. What is the relationship between w_k and the OCF structure?

The OCF says: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
             = 1 + 2*(c3+...+cp) + 4*alpha_2 + 8*alpha_3 + ...

So the "alpha_1 contribution" gives each c_k a weight of 2.
The residual (4*alpha_2 + ...) contributes ADDITIONAL weights.
Therefore: w_k = 2 + (additional weight from disjoint sets containing k-cycles)

If alpha_j = linear(c_k), then we can derive the full w_k.

At p=7: alpha_2 = disj33 = f(c3) only. So alpha_2 doesn't depend on c5, c7.
  H = 1 + 2*(c3+c5+c7) + 4*alpha_2(c3)
  = 2*c5 + 2*c7 + [1 + 2*c3 + 4*alpha_2(c3)]
  So w_5 = w_7 = 2 and the constant absorbs c3 and alpha_2(c3).
"""

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


def solve_linear_system(A_mat, b_vec):
    """Solve A*x = b using exact Fraction arithmetic."""
    n_eq = len(A_mat)
    n_var = len(A_mat[0])
    aug = [[Fraction(a) for a in row] + [Fraction(b)] for row, b in zip(A_mat, b_vec)]

    for col in range(min(n_eq, n_var)):
        pivot = None
        for row in range(col, n_eq):
            if aug[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(n_eq):
            if row == col:
                continue
            factor = aug[row][col] / aug[col][col]
            for j in range(n_var + 1):
                aug[row][j] -= factor * aug[col][j]

    return [aug[i][-1] / aug[i][i] for i in range(min(n_eq, n_var))]


print("=" * 70)
print("H = linear(c_k) COEFFICIENT STRUCTURE")
print("=" * 70)

all_results = {}

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"p = {p}")
    print(f"{'='*70}")

    t0 = time.time()

    # Collect data
    data = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        A = build_adj(p, S)
        H = compute_H_heldkarp(A, p)
        c_k = {}
        for k in range(3, p + 1, 2):
            ck = 0
            for subset in combinations(range(p), k):
                ck += count_ham_cycles(A, list(subset))
            c_k[k] = ck
        data.append({'bits': bits, 'H': H, 'c_k': c_k})

    # Get unique types
    by_H = defaultdict(list)
    for d in data:
        by_H[d['H']].append(d)
    types = [by_H[H][0] for H in sorted(by_H.keys())]

    cycle_lengths = sorted(data[0]['c_k'].keys())
    constant_k = [k for k in cycle_lengths
                  if len(set(d['c_k'][k] for d in data)) == 1]
    varying_k = [k for k in cycle_lengths
                 if len(set(d['c_k'][k] for d in data)) > 1]

    t1 = time.time()
    print(f"  {N} orientations, {len(types)} types, {t1-t0:.1f}s")
    print(f"  Constant c_k: {[(k, data[0]['c_k'][k]) for k in constant_k]}")
    print(f"  Varying c_k: {varying_k}")

    # Test H = 2*alpha_1 + 1 + residual
    print(f"\n  OCF contribution analysis:")
    for d in types:
        alpha_1 = sum(d['c_k'].values())
        R = d['H'] - 1 - 2 * alpha_1
        print(f"    H={d['H']}, alpha_1={alpha_1}, R={R}")

    # Solve for w_k using exact arithmetic
    n_types = len(types)
    n_params = len(varying_k) + 1

    if n_types >= n_params:
        A_mat = [[t['c_k'][k] for k in varying_k] + [1] for t in types[:n_params]]
        b_vec = [t['H'] for t in types[:n_params]]
        coeffs = solve_linear_system(A_mat, b_vec)

        print(f"\n  Exact coefficients (w_k for varying c_k, then constant):")
        for k, c in zip(varying_k + ['const'], coeffs):
            print(f"    w_{k} = {c} = {float(c):.10f}")

        # Verify all types
        all_ok = True
        for d in types:
            pred = sum(coeffs[i] * d['c_k'][varying_k[i]]
                       for i in range(len(varying_k))) + coeffs[-1]
            if pred != d['H']:
                all_ok = False
                print(f"    FAIL at H={d['H']}: pred={pred}")

        if all_ok:
            print(f"  All {n_types} types verified EXACT.")

        # Analyze: w_k - 2 = contribution from higher alpha_j
        print(f"\n  Delta w_k = w_k - 2 (excess beyond alpha_1 contribution):")
        for k, c in zip(varying_k, coeffs[:-1]):
            delta = c - 2
            print(f"    delta_w_{k} = {delta} = {float(delta):.10f}")

        # Common denominator
        denoms = [c.denominator for c in coeffs]
        from math import gcd
        from functools import reduce
        lcm_denom = reduce(lambda a, b: a * b // gcd(a, b), denoms)
        print(f"\n  LCM of denominators: {lcm_denom}")
        print(f"  Integer formula: {lcm_denom} * H = ", end="")
        parts = []
        for k, c in zip(varying_k, coeffs[:-1]):
            int_coeff = int(c * lcm_denom)
            parts.append(f"{int_coeff}*c{k}")
        int_const = int(coeffs[-1] * lcm_denom)
        parts.append(str(int_const))
        print(" + ".join(parts))

        all_results[p] = {
            'varying_k': varying_k,
            'coeffs': coeffs,
            'lcm_denom': lcm_denom
        }

    else:
        print(f"  Underdetermined ({n_types} types < {n_params} params)")
        # For p=7 and p=11, use least squares
        import numpy as np
        y = np.array([d['H'] for d in data], dtype=float)
        X_cols = [np.array([d['c_k'][k] for d in data], dtype=float) for k in varying_k]
        X = np.column_stack(X_cols + [np.ones(N)])
        coeffs_f, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
        print(f"  Floating-point fit (rank={rank}):")
        for k, c in zip(varying_k + ['const'], coeffs_f):
            print(f"    w_{k} = {c:.10f}")

        # Check if w_k = 2 for all varying k
        all_two = all(abs(c - 2.0) < 0.01 for c in coeffs_f[:-1])
        if all_two:
            print(f"  ALL w_k = 2 (H = 2*alpha_1 + constant!)")
            # Compute the constant
            d = types[0]
            const = d['H'] - 2 * sum(d['c_k'].values())
            print(f"  H = 2*alpha_1 + {const + 1}")

        all_results[p] = {
            'varying_k': varying_k,
            'coeffs_float': coeffs_f
        }


# Cross-prime comparison
print(f"\n{'='*70}")
print(f"CROSS-PRIME ANALYSIS")
print(f"{'='*70}")

for p in [7, 11, 13]:
    if 'coeffs' in all_results[p]:
        coeffs = all_results[p]['coeffs']
        varying_k = all_results[p]['varying_k']
        print(f"\n  p={p}:")
        for k, c in zip(varying_k, coeffs[:-1]):
            print(f"    w_{k} = {c} ({float(c):.6f})")
    elif 'coeffs_float' in all_results[p]:
        coeffs = all_results[p]['coeffs_float']
        varying_k = all_results[p]['varying_k']
        print(f"\n  p={p}:")
        for k, c in zip(varying_k, coeffs[:-1]):
            print(f"    w_{k} = {c:.6f}")

print("\nDONE.")
