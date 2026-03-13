#!/usr/bin/env python3
"""
alpha2_moment_analysis.py -- kind-pasteur-2026-03-13-S60

With corrected backtracking (MISTAKE-019 fix), analyze alpha_2's moment dependence.
Key question: is alpha_2 = linear(moments) for all primes?

At p=7: trivially yes (only 2 orbit types).
At p=11: 4 orbit types with 4 moment parameters. Test if alpha_2 = f(S4).
At p=13: 6 orbit types, needs S4+S6 or more.

OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
If alpha_j = linear(moments) for all j, then H = linear(moments).

This script enumerates ALL directed odd cycles and counts disjoint pairs DIRECTLY
(no backtracking) to avoid any bugs.
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import time
import sys


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
    """Held-Karp for total directed Hamiltonian paths."""
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


def enumerate_cycles(A, p):
    """Return list of (frozenset_of_vertices, multiplicity) for all directed odd cycles."""
    cycles = []
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_ham_cycles(A, verts)
            for _ in range(n_cyc):
                cycles.append(frozenset(subset))
    return cycles


def direct_alpha_decomposition(cycles, max_j=3):
    """Count independent sets directly (no backtracking). Only up to max_j."""
    n = len(cycles)
    alpha = [0] * (max_j + 1)
    alpha[0] = 1
    alpha[1] = n

    if max_j >= 2:
        # Count disjoint pairs
        for i in range(n):
            for j in range(i + 1, n):
                if not (cycles[i] & cycles[j]):
                    alpha[2] += 1

    if max_j >= 3:
        # Count disjoint triples
        for i in range(n):
            for j in range(i + 1, n):
                if cycles[i] & cycles[j]:
                    continue
                for k in range(j + 1, n):
                    if not (cycles[i] & cycles[k]) and not (cycles[j] & cycles[k]):
                        alpha[3] += 1

    return alpha


def compute_moments(p, S):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    moments = {}
    for k in range(4, p, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments


# ================================================================
# MAIN
# ================================================================

for p in [7, 11]:
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"ALPHA_2 MOMENT ANALYSIS at p={p}, m={m}, N={N} orientations")
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
        H_hk = compute_H_heldkarp(A, p)
        cycles = enumerate_cycles(A, p)

        # For p=7: max_j=2 suffices (alpha_3=0 proved)
        # For p=11: try max_j=3 but alpha_3 counting might be slow
        max_j = 2 if p == 7 else 3
        alpha = direct_alpha_decomposition(cycles, max_j=max_j)

        moments = compute_moments(p, S)

        # Verify H reconstruction
        H_ocf = sum(alpha[j] * (2**j) for j in range(len(alpha)))
        if max_j < p:
            H_residual = H_hk - H_ocf
        else:
            H_residual = 0

        data.append({
            'bits': bits, 'S': S, 'H': H_hk,
            'alpha': alpha, 'moments': moments,
            'n_cycles': len(cycles), 'H_residual': H_residual
        })

        if bits % 8 == 0:
            elapsed = time.time() - t0
            print(f"  bits={bits}/{N}, elapsed={elapsed:.1f}s, "
                  f"cycles={len(cycles)}, alpha={alpha}, H={H_hk}", flush=True)

    t1 = time.time()
    print(f"\n  Done: {N} orientations in {t1-t0:.1f}s")

    # Group by orbit type (same alpha profile)
    by_alpha = defaultdict(list)
    for d in data:
        key = tuple(d['alpha'])
        by_alpha[key].append(d)

    print(f"\n  Unique alpha profiles: {len(by_alpha)}")
    for prof in sorted(by_alpha.keys()):
        group = by_alpha[prof]
        d = group[0]
        H_vals = sorted(set(d2['H'] for d2 in group))
        print(f"    alpha={list(prof)}, count={len(group)}, H={H_vals}, "
              f"cycles={d['n_cycles']}, residual={d['H_residual']}")

    # Test alpha_2 = f(moments)
    print(f"\n  --- ALPHA_2 MOMENT FIT ---")
    moment_names = [f'S{k}' for k in range(4, p, 2)]

    a2_arr = np.array([d['alpha'][2] for d in data], dtype=float)
    print(f"  alpha_2 values: {sorted(set(int(x) for x in a2_arr))}")

    for n_mom in range(1, len(moment_names) + 1):
        cols = [np.array([d['moments'][4 + 2*i] for d in data], dtype=float)
                for i in range(n_mom)]
        X = np.column_stack(cols + [np.ones(N)])
        coeffs, _, _, _ = np.linalg.lstsq(X, a2_arr, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(a2_arr - pred))
        mom_str = ','.join(moment_names[:n_mom])
        exact = max_err < 0.01
        print(f"    alpha_2 = f({mom_str}): max_err = {max_err:.6f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            coeff_str = ' + '.join(f'{c:.6f}*{name}' for c, name in
                                   zip(coeffs[:-1], moment_names[:n_mom]))
            coeff_str += f' + {coeffs[-1]:.2f}'
            print(f"    = {coeff_str}")
            break

    # Test alpha_1 = f(moments) for comparison
    print(f"\n  --- ALPHA_1 MOMENT FIT ---")
    a1_arr = np.array([d['alpha'][1] for d in data], dtype=float)
    print(f"  alpha_1 values: {sorted(set(int(x) for x in a1_arr))}")

    for n_mom in range(1, len(moment_names) + 1):
        cols = [np.array([d['moments'][4 + 2*i] for d in data], dtype=float)
                for i in range(n_mom)]
        X = np.column_stack(cols + [np.ones(N)])
        coeffs, _, _, _ = np.linalg.lstsq(X, a1_arr, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(a1_arr - pred))
        mom_str = ','.join(moment_names[:n_mom])
        exact = max_err < 0.01
        print(f"    alpha_1 = f({mom_str}): max_err = {max_err:.6f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            coeff_str = ' + '.join(f'{c:.6f}*{name}' for c, name in
                                   zip(coeffs[:-1], moment_names[:n_mom]))
            coeff_str += f' + {coeffs[-1]:.2f}'
            print(f"    = {coeff_str}")
            break

    # If we have alpha_3 data, test it too
    if len(data[0]['alpha']) > 3:
        print(f"\n  --- ALPHA_3 MOMENT FIT ---")
        a3_arr = np.array([d['alpha'][3] if len(d['alpha']) > 3 else 0
                          for d in data], dtype=float)
        a3_vals = sorted(set(int(x) for x in a3_arr))
        print(f"  alpha_3 values: {a3_vals}")

        if max(a3_vals) > 0:
            for n_mom in range(1, len(moment_names) + 1):
                cols = [np.array([d['moments'][4 + 2*i] for d in data], dtype=float)
                        for i in range(n_mom)]
                X = np.column_stack(cols + [np.ones(N)])
                coeffs, _, _, _ = np.linalg.lstsq(X, a3_arr, rcond=None)
                pred = X @ coeffs
                max_err = np.max(np.abs(a3_arr - pred))
                mom_str = ','.join(moment_names[:n_mom])
                exact = max_err < 0.01
                print(f"    alpha_3 = f({mom_str}): max_err = {max_err:.6f} "
                      f"{'EXACT' if exact else ''}")
                if exact:
                    coeff_str = ' + '.join(f'{c:.6f}*{name}' for c, name in
                                           zip(coeffs[:-1], moment_names[:n_mom]))
                    coeff_str += f' + {coeffs[-1]:.2f}'
                    print(f"    = {coeff_str}")
                    break
        else:
            print(f"  alpha_3 = 0 for all orientations (trivial)")

    # H = f(moments) directly
    print(f"\n  --- H MOMENT FIT ---")
    H_arr = np.array([d['H'] for d in data], dtype=float)
    print(f"  H values: {sorted(set(int(x) for x in H_arr))}")

    for n_mom in range(1, len(moment_names) + 1):
        cols = [np.array([d['moments'][4 + 2*i] for d in data], dtype=float)
                for i in range(n_mom)]
        X = np.column_stack(cols + [np.ones(N)])
        coeffs, _, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(H_arr - pred))
        mom_str = ','.join(moment_names[:n_mom])
        exact = max_err < 0.1
        print(f"    H = f({mom_str}): max_err = {max_err:.6f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            coeff_str = ' + '.join(f'{c:.6f}*{name}' for c, name in
                                   zip(coeffs[:-1], moment_names[:n_mom]))
            coeff_str += f' + {coeffs[-1]:.2f}'
            print(f"    = {coeff_str}")
            break

    # Verify: H_residual = 8*alpha_3 + 16*alpha_4 + ...
    if p == 11:
        print(f"\n  --- H RESIDUAL ANALYSIS ---")
        for d in data[:5]:
            residual = d['H_residual']
            print(f"    bits={d['bits']}: H={d['H']}, alpha_0-3={d['alpha']}, "
                  f"residual={residual}")
            H_partial = sum(d['alpha'][j] * (2**j) for j in range(len(d['alpha'])))
            print(f"      H_partial = {H_partial}, residual = H - H_partial = {d['H'] - H_partial}")

print("\nDONE.")
