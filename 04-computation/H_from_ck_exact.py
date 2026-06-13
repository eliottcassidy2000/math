#!/usr/bin/env python3
"""
H_from_ck_exact.py -- kind-pasteur-2026-03-13-S60

Key question: Is H an EXACT linear function of the individual cycle counts
c_3, c_5, ..., c_p for circulant tournaments?

H and all c_k are INTEGERS, so this test involves NO floating point.
If H = sum w_k * c_k + C, we can solve for w_k exactly using Fraction arithmetic.

At p=13: 6 orbit types, 7 cycle types (c3,...,c13). If c3 is constant,
then 6 types and 6 params (c5,...,c13 + const) => exactly determined.
If NOT constant, 6 types and 7 params => underdetermined.

This also reveals WHICH cycle lengths contribute to the higher alpha_j.
"""

from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import numpy as np
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


for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"H = linear(c_k) TEST at p={p}, m={m}, N={N}")
    print(f"{'='*70}")

    t0 = time.time()
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
                verts = list(subset)
                ck += count_ham_cycles(A, verts)
            c_k[k] = ck

        data.append({'bits': bits, 'S': S, 'H': H, 'c_k': c_k})

    t1 = time.time()
    print(f"  {N} orientations in {t1-t0:.1f}s")

    # Show unique types
    by_H = defaultdict(list)
    for d in data:
        by_H[d['H']].append(d)

    cycle_lengths = sorted(data[0]['c_k'].keys())
    print(f"\n  {len(by_H)} unique H values (orbit types)")
    print(f"  Cycle lengths: {cycle_lengths}")

    for H_val in sorted(by_H.keys()):
        group = by_H[H_val]
        d = group[0]
        ck_str = ', '.join(f'c{k}={d["c_k"][k]}' for k in cycle_lengths)
        print(f"  H={H_val} (count={len(group)}): {ck_str}")

    # Check which c_k are constant
    constant_ck = []
    varying_ck = []
    for k in cycle_lengths:
        vals = set(d['c_k'][k] for d in data)
        if len(vals) == 1:
            constant_ck.append((k, list(vals)[0]))
        else:
            varying_ck.append(k)

    print(f"\n  Constant c_k: {[(k, v) for k, v in constant_ck]}")
    print(f"  Varying c_k: {varying_ck}")

    # Test H = linear(varying c_k) using Fraction arithmetic
    n_types = len(by_H)
    n_params = len(varying_ck) + 1  # +1 for constant

    print(f"  {n_types} types, {n_params} parameters "
          f"({'overconstrained' if n_types > n_params else 'underdetermined' if n_types < n_params else 'exactly determined'})")

    # Build the linear system using representatives of each type
    type_reps = [by_H[H_val][0] for H_val in sorted(by_H.keys())]

    if n_types >= n_params:
        # Use first n_params types to solve, verify on rest
        A_mat = []
        b_vec = []
        for d in type_reps[:n_params]:
            row = [Fraction(d['c_k'][k]) for k in varying_ck] + [Fraction(1)]
            A_mat.append(row)
            b_vec.append(Fraction(d['H']))

        # Solve A_mat * x = b_vec using Gaussian elimination
        n_eq = len(A_mat)
        n_var = len(A_mat[0])
        aug = [A_mat[i] + [b_vec[i]] for i in range(n_eq)]

        for col in range(n_var):
            # Find pivot
            pivot = None
            for row in range(col, n_eq):
                if aug[row][col] != 0:
                    pivot = row
                    break
            if pivot is None:
                print(f"  *** No pivot at column {col}! System is singular. ***")
                break
            aug[col], aug[pivot] = aug[pivot], aug[col]
            for row in range(n_eq):
                if row == col:
                    continue
                factor = aug[row][col] / aug[col][col]
                for j in range(n_var + 1):
                    aug[row][j] -= factor * aug[col][j]

        coeffs = [aug[i][-1] / aug[i][i] for i in range(n_var)]
        coeff_names = [f'c{k}' for k in varying_ck] + ['const']

        print(f"\n  Solution (exact rational):")
        for name, coeff in zip(coeff_names, coeffs):
            print(f"    w_{name} = {coeff} = {float(coeff):.8f}")

        # Add constant c_k contributions
        # H = sum w_k * c_k(varying) + w_const + sum w_k_const * c_k(constant)
        # But the constant c_k are absorbed into the constant term
        # Let's separate: total_constant = w_const + sum(known w_k * c_k_const)
        print(f"\n  Formula: H = ", end="")
        parts = []
        for name, coeff in zip(coeff_names[:-1], coeffs[:-1]):
            parts.append(f"({coeff})*{name}")
        parts.append(f"({coeffs[-1]})")
        print(" + ".join(parts))

        # Verify on ALL types
        print(f"\n  Verification on all {n_types} types:")
        all_exact = True
        for d in type_reps:
            pred = sum(coeffs[i] * d['c_k'][varying_ck[i]]
                      for i in range(len(varying_ck))) + coeffs[-1]
            error = d['H'] - pred
            status = "OK" if error == 0 else f"ERROR={error}"
            if error != 0:
                all_exact = False
            print(f"    H={d['H']}: pred={pred}, {status}")

        if all_exact:
            print(f"\n  *** H = EXACT linear function of ({', '.join(f'c{k}' for k in varying_ck)}) ***")

            # Expand to include constant c_k
            total_const = coeffs[-1]
            for k, v in constant_ck:
                # We need to know the coefficient of the constant c_k
                # Since they're constant, they're absorbed into the intercept
                pass

            # Check if coefficients are "nice" (simple fractions)
            print(f"\n  Coefficient denominators: {[c.denominator for c in coeffs]}")

    else:
        print(f"  System underdetermined, testing floating point fit...")
        y = np.array([d['H'] for d in data], dtype=float)
        X_cols = [np.array([d['c_k'][k] for d in data], dtype=float) for k in varying_ck]
        X = np.column_stack(X_cols + [np.ones(N)])
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(y - pred))
        print(f"  max_err = {max_err:.6f}")

    # Also test: R = H - 1 - 2*alpha_1 as function of c_k
    print(f"\n  --- RESIDUAL R = H - 1 - 2*alpha_1 = f(c_k) ---")
    for d in type_reps:
        alpha_1 = sum(d['c_k'].values())
        R = d['H'] - 1 - 2 * alpha_1
        R_over_4 = Fraction(R, 4)
        print(f"    H={d['H']}, a1={alpha_1}, R={R}, R/4={R_over_4}")

print("\nDONE.")
