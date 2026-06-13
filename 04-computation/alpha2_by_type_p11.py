#!/usr/bin/env python3
"""
alpha2_by_type_p11.py -- kind-pasteur-2026-03-13-S60

At p=11, alpha_2 = disj_{3,3} + disj_{3,5} + disj_{3,7} + disj_{5,5}.
Each component should individually be a linear function of moments.

This script computes alpha_2 decomposed by cycle-length pair type,
and tests each component for moment linearity.

Also computes alpha_3 = disj_{3,3,3} + disj_{3,3,5}.
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


p = 11
m = (p - 1) // 2
N = 1 << m

print(f"{'='*70}")
print(f"ALPHA_2 BY TYPE at p={p}, m={m}, N={N}")
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

    # Enumerate cycles by length
    cycles_by_k = {}
    for k in range(3, p + 1, 2):
        cyc_k = []
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = count_ham_cycles(A, verts)
            for _ in range(n_cyc):
                cyc_k.append(frozenset(subset))
        if cyc_k:
            cycles_by_k[k] = cyc_k

    c_k = {k: len(v) for k, v in cycles_by_k.items()}
    alpha_1 = sum(c_k.values())

    # Count disjoint pairs by type
    disj = {}
    for k1 in sorted(cycles_by_k.keys()):
        for k2 in sorted(cycles_by_k.keys()):
            if k1 > k2:
                continue
            if k1 + k2 > p:
                continue
            count = 0
            cyc1 = cycles_by_k[k1]
            cyc2 = cycles_by_k[k2]
            if k1 == k2:
                for i in range(len(cyc1)):
                    for j in range(i + 1, len(cyc1)):
                        if not (cyc1[i] & cyc1[j]):
                            count += 1
            else:
                for c1 in cyc1:
                    for c2 in cyc2:
                        if not (c1 & c2):
                            count += 1
            if count > 0:
                disj[(k1, k2)] = count

    alpha_2 = sum(disj.values())

    # Count disjoint triples by type
    # At p=11: only (3,3,3) and (3,3,5) possible
    disj3 = {}
    c3 = cycles_by_k.get(3, [])
    c5 = cycles_by_k.get(5, [])

    # (3,3,3) triples
    count_333 = 0
    for i in range(len(c3)):
        for j in range(i + 1, len(c3)):
            if c3[i] & c3[j]:
                continue
            for k in range(j + 1, len(c3)):
                if not (c3[i] & c3[k]) and not (c3[j] & c3[k]):
                    count_333 += 1
    disj3[(3, 3, 3)] = count_333

    # (3,3,5) triples: choose 2 disjoint 3-cycles and 1 5-cycle disjoint from both
    count_335 = 0
    for i in range(len(c3)):
        for j in range(i + 1, len(c3)):
            if c3[i] & c3[j]:
                continue
            union_ij = c3[i] | c3[j]
            for c5_cyc in c5:
                if not (union_ij & c5_cyc):
                    count_335 += 1
    disj3[(3, 3, 5)] = count_335

    alpha_3 = count_333 + count_335

    # Verify H
    H_ocf = 1 + 2 * alpha_1 + 4 * alpha_2 + 8 * alpha_3
    match = (H_ocf == H_hk)

    moments = compute_moments(p, S)

    data.append({
        'bits': bits, 'S': S, 'H': H_hk, 'c_k': c_k,
        'alpha_1': alpha_1, 'alpha_2': alpha_2, 'alpha_3': alpha_3,
        'disj': disj, 'disj3': disj3, 'moments': moments,
        'match': match
    })

    elapsed = time.time() - t0
    if bits % 4 == 0:
        print(f"  bits={bits}/{N}: c_k={c_k}, a1={alpha_1}, a2={alpha_2}, "
              f"a3={alpha_3}, H={H_hk}, match={match}, t={elapsed:.0f}s", flush=True)

t1 = time.time()
print(f"\nCompleted {N} orientations in {t1-t0:.0f}s")

# Check all H matches
all_match = all(d['match'] for d in data)
print(f"\nH(OCF) = H(HK) for ALL orientations: {all_match}")
if not all_match:
    for d in data:
        if not d['match']:
            print(f"  MISMATCH at bits={d['bits']}: OCF={1+2*d['alpha_1']+4*d['alpha_2']+8*d['alpha_3']}, HK={d['H']}")

# Show unique profiles
by_profile = defaultdict(list)
for d in data:
    key = (d['alpha_1'], d['alpha_2'], d['alpha_3'])
    by_profile[key].append(d)

print(f"\nUnique alpha profiles: {len(by_profile)}")
for prof in sorted(by_profile.keys()):
    group = by_profile[prof]
    d = group[0]
    H = d['H']
    disj_str = ', '.join(f'd{k}={v}' for k, v in sorted(d['disj'].items()))
    disj3_str = ', '.join(f'd{k}={v}' for k, v in sorted(d['disj3'].items()))
    print(f"  a1={prof[0]}, a2={prof[1]}, a3={prof[2]}, H={H}, count={len(group)}")
    print(f"    c_k = {d['c_k']}")
    print(f"    disj2: {disj_str}")
    print(f"    disj3: {disj3_str}")

# Test moment linearity
moment_names = [f'S{k}' for k in range(4, p, 2)]

for label, arr_fn in [
    ('alpha_1', lambda d: d['alpha_1']),
    ('alpha_2', lambda d: d['alpha_2']),
    ('alpha_3', lambda d: d['alpha_3']),
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
        exact = max_err < 0.1
        print(f"  {label} = f({mom_str}): max_err = {max_err:.4f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            coeff_str = ' + '.join(f'{c:.8f}*{name}' for c, name in
                                   zip(coeffs[:-1], moment_names[:n_mom]))
            coeff_str += f' + {coeffs[-1]:.2f}'
            print(f"  = {coeff_str}")
            break

# Also test disj_{k1,k2} components individually
print(f"\n--- DISJOINT PAIR COMPONENT FITS ---")
disj_types = set()
for d in data:
    disj_types.update(d['disj'].keys())

for dt in sorted(disj_types):
    label = f'disj_{dt}'
    y = np.array([d['disj'].get(dt, 0) for d in data], dtype=float)
    vals = sorted(set(int(x) for x in y))
    print(f"\n  {label} values: {vals}")

    if len(vals) <= 1:
        print(f"    CONSTANT = {vals[0]}")
        continue

    for n_mom in range(1, len(moment_names) + 1):
        cols = [np.array([d['moments'][4 + 2*i] for d in data], dtype=float)
                for i in range(n_mom)]
        X = np.column_stack(cols + [np.ones(N)])
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(y - pred))
        mom_str = ','.join(moment_names[:n_mom])
        exact = max_err < 0.1
        print(f"    {label} = f({mom_str}): max_err = {max_err:.4f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            coeff_str = ' + '.join(f'{c:.8f}*{name}' for c, name in
                                   zip(coeffs[:-1], moment_names[:n_mom]))
            coeff_str += f' + {coeffs[-1]:.2f}'
            print(f"    = {coeff_str}")
            break

print("\nDONE.")
