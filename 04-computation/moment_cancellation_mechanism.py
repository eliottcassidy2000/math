#!/usr/bin/env python3
"""
moment_cancellation_mechanism.py -- kind-pasteur-2026-03-13-S60

At p=11, H = f(S4, S6, S8) exactly -- S10 cancels from the formula.
At p=13, H needs all 5 moments (S4,...,S12) -- no cancellation.

The S10 coeff in alpha_1 = sum c_k is exactly 1 (not 0).
So cancellation must come from the full OCF: H = sum alpha_j * 2^j.

This script:
1. Computes the full alpha_j decomposition at p=7 and p=11
2. For each alpha_j, tests its dependence on eigenvalue moments
3. Determines which alpha_j contributes the cancellation
4. Explores whether this is related to p mod 4

The key insight: alpha_j for j >= 3 involves overlapping cycle
interactions, which may create moment dependencies beyond c_k.
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def compute_eigenvalue_moments(p, S, max_moment=12):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments


def enumerate_all_directed_odd_cycles(A, p, max_k=None):
    """Enumerate all directed odd cycles, returning frozenset of vertex sets."""
    if max_k is None:
        max_k = p
    cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_ham_cycles(A, verts)
            for _ in range(n_cyc):
                cycles.append(frozenset(subset))
    return cycles


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


def full_alpha_decomposition(A, p, cycles):
    """Compute all alpha_j for the conflict graph of cycles."""
    n = len(cycles)
    # Build conflict graph
    nbr = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if cycles[i] & cycles[j]:
                nbr[i] |= (1 << j)
                nbr[j] |= (1 << i)

    alpha = [0] * (n + 1)

    def backtrack(v, mask, size):
        alpha[size] += 1
        for w in range(v + 1, n):
            if not (mask & (1 << w)):
                backtrack(w + 1, mask | nbr[w], size + 1)

    if n <= 25:
        backtrack(0, 0, 0)
    else:
        # Fallback: only compute alpha_0, alpha_1, alpha_2, alpha_3
        alpha[0] = 1
        alpha[1] = n
        alpha[2] = 0
        alpha[3] = 0
        adj = [[False]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if cycles[i] & cycles[j]:
                    adj[i][j] = True
                    adj[j][i] = True

        for i in range(n):
            for j in range(i+1, n):
                if not adj[i][j]:
                    alpha[2] += 1
        for i in range(n):
            for j in range(i+1, n):
                if adj[i][j]:
                    continue
                for k in range(j+1, n):
                    if not adj[i][k] and not adj[j][k]:
                        alpha[3] += 1

    return alpha


def analyze_alpha_moment_dependence(p):
    """For each alpha_j, test its dependence on eigenvalue moments."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"ALPHA-MOMENT DEPENDENCE at p={p}")
    print(f"{'='*70}")

    # Collect data for all orientations
    data = []
    import time
    t0 = time.time()

    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        moments = compute_eigenvalue_moments(p, S, max_moment=2*m)
        cycles = enumerate_all_directed_odd_cycles(A, p)
        alpha = full_alpha_decomposition(A, p, cycles)

        # Trim trailing zeros
        max_j = max(j for j in range(len(alpha)) if alpha[j] > 0) if any(alpha) else 0

        H = sum(alpha[j] * (2**j) for j in range(len(alpha)))

        data.append({
            'bits': bits, 'S': S, 'H': H,
            'moments': moments, 'alpha': alpha[:max_j+1]
        })

    t1 = time.time()
    print(f"  Computed {N} orientations in {t1-t0:.1f}s")

    # Find max alpha index
    max_alpha = max(len(d['alpha']) for d in data)
    print(f"  Max alpha index: {max_alpha - 1}")

    # For each alpha_j, test moment dependence
    moment_names = [f'S{2*j}' for j in range(2, m + 1)]
    max_moments = len(moment_names)

    for j in range(max_alpha):
        aj_arr = np.array([d['alpha'][j] if j < len(d['alpha']) else 0
                          for d in data], dtype=float)

        if np.std(aj_arr) < 1e-6:
            print(f"\n  alpha_{j}: CONSTANT = {int(aj_arr[0])}")
            continue

        print(f"\n  alpha_{j} dependence:")
        unique_vals = sorted(set(int(x) for x in aj_arr))
        print(f"    Values: {unique_vals}")

        for n_mom in range(1, max_moments + 1):
            cols = [np.array([d['moments'][2*(i+2)] for d in data], dtype=float)
                    for i in range(n_mom)]
            X = np.column_stack(cols + [np.ones(N)])
            coeffs, _, _, _ = np.linalg.lstsq(X, aj_arr, rcond=None)
            pred = X @ coeffs
            max_err = np.max(np.abs(aj_arr - pred))
            mom_str = ','.join(moment_names[:n_mom])
            exact = max_err < 0.01
            print(f"    alpha_{j} = f({mom_str}): max_err = {max_err:.4f} "
                  f"{'EXACT' if exact else ''}")
            if exact:
                # Show coefficients
                coeff_str = ' + '.join(f'{c:.6f}*{name}' for c, name in
                                       zip(coeffs[:-1], moment_names[:n_mom]))
                coeff_str += f' + {coeffs[-1]:.2f}'
                print(f"    = {coeff_str}")
                break

    # Special: test S_{p-1} coefficient in each alpha_j
    Sp1 = 2 * m  # S_{p-1} = S_{2m}
    sp1_name = f'S{Sp1}'
    print(f"\n  --- {sp1_name} (highest moment) dependence ---")

    for j in range(max_alpha):
        aj_arr = np.array([d['alpha'][j] if j < len(d['alpha']) else 0
                          for d in data], dtype=float)
        if np.std(aj_arr) < 1e-6:
            continue

        Sp1_arr = np.array([d['moments'][Sp1] for d in data], dtype=float)
        # Partial correlation: alpha_j vs S_{p-1} controlling for S4,...,S_{p-3}
        # Just compute raw correlation
        r = np.corrcoef(aj_arr, Sp1_arr)[0, 1]
        print(f"    Corr(alpha_{j}, {sp1_name}) = {r:.6f}")

    # Test: does removing S_{p-1} from the alpha_j formula break exactness for any j?
    print(f"\n  --- Can we remove {sp1_name} from alpha_j formulas? ---")

    # First check: H = f(S4,...,S_{p-3}) -- does the SECOND-highest work?
    H_arr = np.array([d['H'] for d in data], dtype=float)
    reduced_moments = max_moments - 1  # all except S_{p-1}
    if reduced_moments > 0:
        cols = [np.array([d['moments'][2*(i+2)] for d in data], dtype=float)
                for i in range(reduced_moments)]
        X = np.column_stack(cols + [np.ones(N)])
        coeffs, _, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(H_arr - pred))
        mom_str = ','.join(moment_names[:reduced_moments])
        print(f"    H = f({mom_str}): max_err = {max_err:.4f} "
              f"{'EXACT' if max_err < 0.01 else ''}")

    # Full moments
    cols_full = [np.array([d['moments'][2*(i+2)] for d in data], dtype=float)
                 for i in range(max_moments)]
    X_full = np.column_stack(cols_full + [np.ones(N)])
    coeffs_full, _, _, _ = np.linalg.lstsq(X_full, H_arr, rcond=None)
    pred_full = X_full @ coeffs_full
    max_err_full = np.max(np.abs(H_arr - pred_full))
    mom_str_full = ','.join(moment_names)
    print(f"    H = f({mom_str_full}): max_err = {max_err_full:.4f} "
          f"{'EXACT' if max_err_full < 0.01 else ''}")

    # If reduced works: S_{p-1} cancellation confirmed!
    if reduced_moments > 0 and max_err < 0.01:
        print(f"    *** {sp1_name} CANCELLATION CONFIRMED at p={p} ***")
        # Show the coefficient that would be S_{p-1}'s if we include it
        print(f"    {sp1_name} coefficient in full fit = {coeffs_full[-2]:.8f}")

    return data


# ================================================================
# MAIN
# ================================================================
print("=" * 70)
print("MOMENT CANCELLATION MECHANISM")
print("=" * 70)

for p in [7, 11]:
    analyze_alpha_moment_dependence(p)

print("\nDONE.")
