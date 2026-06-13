#!/usr/bin/env python3
"""
alpha3_moment_analysis.py -- kind-pasteur-2026-03-13-S60

With THM-157 (alpha_1 = linear(moments)) and THM-156 (alpha_2 = linear(S4)),
the question of whether H = linear(moments) reduces to whether
alpha_3 (and higher) are linear functions of eigenvalue moments.

H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

At p=7: alpha_3 = 0 for all orientations (only alpha_0=1, alpha_1, alpha_2)
So H = 1 + 2*alpha_1 + 4*alpha_2 automatically linear!

At p=11: alpha_3, alpha_4 may be nonzero.
The question: is alpha_3 = linear(S4,...,S10)?
This would confirm H = linear(moments) at p=11 (which we know is true trivially).

More importantly: at p=17 (if we could compute H), would 8*alpha_3 + ... be
a linear function of moments? If not, H is NOT linear at p>=17.

This script:
1. Computes full alpha decomposition at p=7 (quick)
2. Tests alpha_j moment dependence
3. Derives the theoretical structure of alpha_3
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
from math import comb
import time


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def enumerate_directed_cycles(A, p, max_k=None):
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


def full_alpha(cycles):
    n = len(cycles)
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
                backtrack(w, mask | nbr[w], size + 1)

    if n <= 25:
        backtrack(-1, 0, 0)
    return alpha


def compute_moments(p, S, max_moment=None):
    m = (p - 1) // 2
    if max_moment is None:
        max_moment = 2 * m
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments


def analyze_alpha_j(p):
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"FULL ALPHA DECOMPOSITION at p={p}, m={m}")
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
        cycles = enumerate_directed_cycles(A, p)
        alpha = full_alpha(cycles)
        moments = compute_moments(p, S)

        H = sum(alpha[j] * (2**j) for j in range(len(alpha)))

        nonzero = [j for j in range(len(alpha)) if alpha[j] > 0]
        max_j = max(nonzero) if nonzero else 0
        data.append({
            'bits': bits, 'S': S, 'H': H,
            'alpha': alpha[:max_j+1], 'moments': moments
        })

    t1 = time.time()
    max_alpha_idx = max(len(d['alpha']) for d in data) - 1
    print(f"  {N} orientations in {t1-t0:.1f}s")
    print(f"  Max alpha index: {max_alpha_idx}")

    # Show all unique alpha profiles
    by_alpha = defaultdict(list)
    for d in data:
        key = tuple(d['alpha'])
        by_alpha[key].append(d)

    print(f"\n  Unique alpha profiles: {len(by_alpha)}")
    for alpha_prof in sorted(by_alpha.keys()):
        d = by_alpha[alpha_prof][0]
        H = d['H']
        count = len(by_alpha[alpha_prof])
        alpha_str = ', '.join(f'a{j}={alpha_prof[j]}' for j in range(len(alpha_prof)))
        H_check = sum(alpha_prof[j] * (2**j) for j in range(len(alpha_prof)))
        print(f"    H={H:>8} (count={count:>3}): {alpha_str}")

    # Test each alpha_j for moment dependence
    moment_names = [f'S{2*j}' for j in range(2, m+1)]

    for j in range(max_alpha_idx + 1):
        aj_arr = np.array([d['alpha'][j] if j < len(d['alpha']) else 0
                          for d in data], dtype=float)
        if np.std(aj_arr) < 1e-6:
            print(f"\n  alpha_{j}: CONSTANT = {int(aj_arr[0])}")
            continue

        print(f"\n  alpha_{j} values: {sorted(set(int(x) for x in aj_arr))}")
        for n_mom in range(1, len(moment_names) + 1):
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
                coeff_str = ' + '.join(f'{c:.6f}*{name}' for c, name in
                                       zip(coeffs[:-1], moment_names[:n_mom]))
                coeff_str += f' + {coeffs[-1]:.2f}'
                print(f"    = {coeff_str}")
                break

    # Verify: H = 1 + 2*alpha_1 + 4*alpha_2 + ...
    print(f"\n  --- H DECOMPOSITION VERIFICATION ---")
    for d in data[:3]:
        alpha = d['alpha']
        parts = [f"{alpha[j]}*{2**j}" for j in range(len(alpha)) if alpha[j] > 0]
        total = sum(alpha[j] * (2**j) for j in range(len(alpha)))
        print(f"    bits={d['bits']}: H={d['H']} = {' + '.join(parts)} = {total}")

    # Derive: 8*alpha_3 + 16*alpha_4 + ... as residual
    print(f"\n  --- HIGHER ALPHA RESIDUAL ---")
    for d in data:
        alpha = d['alpha']
        alpha_1 = alpha[1] if len(alpha) > 1 else 0
        alpha_2 = alpha[2] if len(alpha) > 2 else 0
        residual = d['H'] - 1 - 2*alpha_1 - 4*alpha_2
        high_alpha = sum(alpha[j]*(2**j) for j in range(3, len(alpha)))
        if residual != high_alpha:
            print(f"    MISMATCH at bits={d['bits']}: residual={residual}, high={high_alpha}")

    # If alpha_3+ are all zero, H = 1 + 2*a1 + 4*a2 trivially
    all_zero_3plus = all(
        all(d['alpha'][j] == 0 for j in range(3, len(d['alpha'])))
        for d in data
    )
    if all_zero_3plus:
        print(f"  All alpha_j = 0 for j >= 3!")
        print(f"  H = 1 + 2*alpha_1 + 4*alpha_2")
        print(f"  This is AUTOMATICALLY linear in moments since alpha_1 and alpha_2 are linear.")
    else:
        # Show non-zero higher alpha
        for d in data:
            for j in range(3, len(d['alpha'])):
                if d['alpha'][j] > 0:
                    print(f"    alpha_{j} = {d['alpha'][j]} at bits={d['bits']}")

    return data


# ================================================================
# MAIN
# ================================================================
print("=" * 70)
print("ALPHA_3+ ANALYSIS: Do higher alphas affect H linearity?")
print("=" * 70)

# p=7: should be quick
data_7 = analyze_alpha_j(7)

# p=11: might be slow due to large number of cycles
# Let's try it but with a timeout warning
print(f"\n*** Attempting p=11 (may be slow due to ~21000 cycles) ***")
try:
    data_11 = analyze_alpha_j(11)
except Exception as e:
    print(f"  p=11 failed: {e}")

print("\nDONE.")
