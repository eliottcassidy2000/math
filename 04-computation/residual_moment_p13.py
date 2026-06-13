#!/usr/bin/env python3
"""
residual_moment_p13.py -- kind-pasteur-2026-03-13-S60

Fast approach to test alpha_j linearity at p=13 WITHOUT full cycle enumeration.

R = H - 1 - 2*alpha_1 = 4*alpha_2 + 8*alpha_3 + 16*alpha_4

H is fast (Held-Karp O(2^p * p^2)).
alpha_1 = sum c_k is fast (cycle enumeration for each k, or matrix trace).

If R is linear in moments, and alpha_1 is linear (THM-157), then
4*alpha_2 + 8*alpha_3 + 16*alpha_4 is linear. Combined with the
decomposition data from a few orientations, this strongly constrains
individual alpha_j.

Also compute: disj_{3,3} from THM-156 (function of S4 only).
Then R2 = R - 4*disj_{3,3} = 4*(alpha_2 - disj_{3,3}) + 8*alpha_3 + 16*alpha_4
       = 4*(higher pair types) + 8*alpha_3 + 16*alpha_4

Test each of these for moment linearity.
"""

import cmath
import numpy as np
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


def compute_alpha1(A, p):
    """Count total directed odd cycles by enumerating each odd-length subset."""
    total = 0
    by_k = {}
    for k in range(3, p + 1, 2):
        ck = 0
        for subset in combinations(range(p), k):
            verts = list(subset)
            ck += count_ham_cycles(A, verts)
        by_k[k] = ck
        total += ck
    return total, by_k


def compute_disj33(cycles_3):
    """Count disjoint (3,3) pairs."""
    n = len(cycles_3)
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if not (cycles_3[i] & cycles_3[j]):
                count += 1
    return count


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


p = 13
m = (p - 1) // 2
N = 1 << m

print(f"{'='*70}")
print(f"RESIDUAL MOMENT ANALYSIS at p={p}, m={m}, N={N}")
print(f"R = H - 1 - 2*alpha_1 = 4*a2 + 8*a3 + 16*a4")
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

    # Fast computations
    H = compute_H_heldkarp(A, p)
    alpha_1, c_k = compute_alpha1(A, p)

    # Enumerate 3-cycles for disj_{3,3}
    c3_list = []
    for subset in combinations(range(p), 3):
        verts = list(subset)
        nc = count_ham_cycles(A, verts)
        for _ in range(nc):
            c3_list.append(frozenset(subset))
    disj33 = compute_disj33(c3_list)

    R = H - 1 - 2 * alpha_1  # = 4*a2 + 8*a3 + 16*a4
    R2 = R - 4 * disj33  # = 4*(a2 - disj33) + 8*a3 + 16*a4 = 4*(higher pair types) + ...

    moments = compute_moments(p, S)

    data.append({
        'bits': bits, 'S': S, 'H': H, 'alpha_1': alpha_1,
        'c_k': c_k, 'disj33': disj33, 'R': R, 'R2': R2,
        'moments': moments
    })

    elapsed = time.time() - t0
    if bits % 8 == 0:
        print(f"  bits={bits}/{N}: H={H}, a1={alpha_1}, R={R}, "
              f"disj33={disj33}, R2={R2}, t={elapsed:.0f}s", flush=True)

t1 = time.time()
print(f"\nCompleted {N} orientations in {t1-t0:.0f}s")

# Group by orbit type
by_type = defaultdict(list)
for d in data:
    key = d['H']  # H determines the type at p=13 (verified earlier)
    by_type[key].append(d)

print(f"\n{len(by_type)} orbit types:")
for H_val in sorted(by_type.keys()):
    group = by_type[H_val]
    d = group[0]
    print(f"  H={H_val}: a1={d['alpha_1']}, R={d['R']}, "
          f"disj33={d['disj33']}, R2={d['R2']}, count={len(group)}")

# Moment linearity tests
moment_names = [f'S{k}' for k in range(4, p, 2)]

for label, arr_fn in [
    ('alpha_1', lambda d: d['alpha_1']),
    ('R = 4a2+8a3+16a4', lambda d: d['R']),
    ('disj_{3,3}', lambda d: d['disj33']),
    ('R2 = R-4*disj33', lambda d: d['R2']),
    ('H', lambda d: d['H']),
]:
    print(f"\n--- {label} MOMENT FIT ---")
    y = np.array([arr_fn(d) for d in data], dtype=float)
    vals = sorted(set(int(x) for x in y))
    print(f"  values: {vals}")

    if len(vals) <= 1:
        print(f"  CONSTANT = {vals[0]}")
        continue

    for n_mom in range(1, len(moment_names) + 1):
        cols = [np.array([d['moments'][4 + 2*i] for d in data], dtype=float)
                for i in range(n_mom)]
        X = np.column_stack(cols + [np.ones(N)])
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(y - pred))
        mom_str = ','.join(moment_names[:n_mom])
        exact = max_err < 1.0  # tolerance for large numbers
        print(f"  {label} = f({mom_str}): max_err = {max_err:.4f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            coeff_str = ' + '.join(f'{c:.8f}*{name}' for c, name in
                                   zip(coeffs[:-1], moment_names[:n_mom]))
            coeff_str += f' + {coeffs[-1]:.2f}'
            print(f"  = {coeff_str}")
            break

# KEY TEST: Can we separate alpha_2 from alpha_3 from alpha_4 using
# the known data points from the slow computation?
# From the slow script: bits=0 gives a2=436748, a3=87568, a4=3224
# R(bits=0) = 4*436748 + 8*87568 + 16*3224 = 1746992 + 700544 + 51584 = 2499120
# R(bits=0) should also = H - 1 - 2*a1 = 3711175 - 1 - 2*606027 = 3711175 - 1212055 = 2499120 ✓

print(f"\n--- CROSS-CHECK WITH KNOWN alpha_j ---")
# Known from slow computation:
known = [
    (0, 606027, 436748, 87568, 3224, 3711175),
    (8, 645937, 439452, 77428, 2119, 3703011),
]
for bits, a1, a2, a3, a4, H in known:
    d = next(d for d in data if d['bits'] == bits)
    R_computed = d['R']
    R_expected = 4*a2 + 8*a3 + 16*a4
    R_from_Ha1 = H - 1 - 2*a1
    print(f"  bits={bits}: R={R_computed}, expected={R_expected}, "
          f"H-1-2a1={R_from_Ha1}, match={R_computed==R_expected==R_from_Ha1}")

print("\nDONE.")
