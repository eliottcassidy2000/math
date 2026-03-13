#!/usr/bin/env python3
"""
alpha_ocf_p13.py -- kind-pasteur-2026-03-13-S60

At p=13, m=6, N=64. alpha_j = 0 for j >= 5 (since 5*3=15>13).
H = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4.

Test: is each alpha_j = f(S4,...,S10) exactly? (m-2 = 4 moments)
S12 should cancel.

Challenges: cycle counts are large (c3=C(12,2)/3*2=52, c5~O(1000), etc.)
alpha_2 requires O(n^2) pairs, alpha_3 O(n^3) triples. Might be slow.

Strategy: compute alpha_j directly from cycle lists (no backtracking).
For alpha_4: only (3,3,3,3) possible (3+3+3+3=12<13).
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict
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


p = 13
m = (p - 1) // 2
N = 1 << m

print(f"{'='*70}")
print(f"ALPHA OCF at p={p}, m={m}, N={N}")
print(f"Max alpha index: {p//3} (since {(p//3)+1}*3={((p//3)+1)*3} > {p})")
print(f"{'='*70}")

# Possible disjoint pair types: k1+k2 <= 13
# (3,3)=6, (3,5)=8, (3,7)=10, (5,5)=10, (3,9)=12, (5,7)=12, (3,11)=14>13, etc.
print(f"Possible pair types: (3,3), (3,5), (3,7), (5,5), (3,9), (5,7)")
# Triple types: sum <= 13
# (3,3,3)=9, (3,3,5)=11, (3,3,7)=13, (3,5,5)=13
print(f"Possible triple types: (3,3,3), (3,3,5), (3,3,7), (3,5,5)")
# Quadruple types: sum <= 13
# (3,3,3,3)=12
print(f"Possible quadruple types: (3,3,3,3)")

t_start = time.time()
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
    disj2 = {}
    pair_types = [(3,3), (3,5), (3,7), (5,5), (3,9), (5,7)]
    for k1, k2 in pair_types:
        if k1 not in cycles_by_k or k2 not in cycles_by_k:
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
        disj2[(k1, k2)] = count

    alpha_2 = sum(disj2.values())

    # Count disjoint triples
    c3 = cycles_by_k.get(3, [])
    c5 = cycles_by_k.get(5, [])
    c7 = cycles_by_k.get(7, [])

    # (3,3,3): three mutually disjoint 3-cycles covering 9 vertices
    count_333 = 0
    for i in range(len(c3)):
        for j in range(i + 1, len(c3)):
            if c3[i] & c3[j]:
                continue
            for k_ in range(j + 1, len(c3)):
                if not (c3[i] & c3[k_]) and not (c3[j] & c3[k_]):
                    count_333 += 1

    # (3,3,5): two disjoint 3-cycles + disjoint 5-cycle
    count_335 = 0
    for i in range(len(c3)):
        for j in range(i + 1, len(c3)):
            if c3[i] & c3[j]:
                continue
            union_ij = c3[i] | c3[j]
            for c5_cyc in c5:
                if not (union_ij & c5_cyc):
                    count_335 += 1

    # (3,3,7): two disjoint 3-cycles + disjoint 7-cycle (3+3+7=13=p, covers all!)
    count_337 = 0
    for i in range(len(c3)):
        for j in range(i + 1, len(c3)):
            if c3[i] & c3[j]:
                continue
            union_ij = c3[i] | c3[j]
            for c7_cyc in c7:
                if not (union_ij & c7_cyc):
                    count_337 += 1

    # (3,5,5): one 3-cycle + two disjoint 5-cycles (3+5+5=13=p, covers all!)
    count_355 = 0
    for i in range(len(c5)):
        for j in range(i + 1, len(c5)):
            if c5[i] & c5[j]:
                continue
            union_ij = c5[i] | c5[j]
            for c3_cyc in c3:
                if not (union_ij & c3_cyc):
                    count_355 += 1

    alpha_3 = count_333 + count_335 + count_337 + count_355

    # Count disjoint quadruples: only (3,3,3,3) with 3*4=12 < 13
    count_3333 = 0
    # Enumerate 4 mutually disjoint 3-cycles
    for i in range(len(c3)):
        for j in range(i + 1, len(c3)):
            if c3[i] & c3[j]:
                continue
            for k_ in range(j + 1, len(c3)):
                if (c3[i] & c3[k_]) or (c3[j] & c3[k_]):
                    continue
                for l in range(k_ + 1, len(c3)):
                    if not (c3[i] & c3[l]) and not (c3[j] & c3[l]) and not (c3[k_] & c3[l]):
                        count_3333 += 1

    alpha_4 = count_3333

    H_ocf = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + 16*alpha_4
    match = (H_ocf == H_hk)

    moments = compute_moments(p, S)

    data.append({
        'bits': bits, 'S': S, 'H': H_hk, 'c_k': c_k,
        'alpha_1': alpha_1, 'alpha_2': alpha_2, 'alpha_3': alpha_3, 'alpha_4': alpha_4,
        'disj2': disj2, 'moments': moments, 'match': match,
        'disj3_333': count_333, 'disj3_335': count_335,
        'disj3_337': count_337, 'disj3_355': count_355,
        'disj4_3333': count_3333
    })

    elapsed = time.time() - t_start
    if bits % 8 == 0:
        print(f"  bits={bits}/{N}: a1={alpha_1}, a2={alpha_2}, "
              f"a3={alpha_3}, a4={alpha_4}, H={H_hk}, match={match}, "
              f"t={elapsed:.0f}s", flush=True)

t1 = time.time()
print(f"\nCompleted {N} orientations in {t1-t_start:.0f}s")

all_match = all(d['match'] for d in data)
print(f"H(OCF) = H(HK) for ALL: {all_match}")
if not all_match:
    for d in data:
        if not d['match']:
            print(f"  MISMATCH bits={d['bits']}: OCF={1+2*d['alpha_1']+4*d['alpha_2']+8*d['alpha_3']+16*d['alpha_4']}, HK={d['H']}")

# Unique profiles
by_profile = defaultdict(list)
for d in data:
    key = (d['alpha_1'], d['alpha_2'], d['alpha_3'], d['alpha_4'])
    by_profile[key].append(d)

print(f"\nUnique alpha profiles: {len(by_profile)}")
for prof in sorted(by_profile.keys(), key=lambda x: -x[0]):
    group = by_profile[prof]
    d = group[0]
    print(f"  a1={prof[0]}, a2={prof[1]}, a3={prof[2]}, a4={prof[3]}, "
          f"H={d['H']}, count={len(group)}")

# Moment linearity tests
moment_names = [f'S{k}' for k in range(4, p, 2)]

for label, arr_fn in [
    ('alpha_1', lambda d: d['alpha_1']),
    ('alpha_2', lambda d: d['alpha_2']),
    ('alpha_3', lambda d: d['alpha_3']),
    ('alpha_4', lambda d: d['alpha_4']),
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
        exact = max_err < 0.5  # slightly looser tolerance for larger numbers
        print(f"  {label} = f({mom_str}): max_err = {max_err:.4f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            coeff_str = ' + '.join(f'{c:.8f}*{name}' for c, name in
                                   zip(coeffs[:-1], moment_names[:n_mom]))
            coeff_str += f' + {coeffs[-1]:.2f}'
            print(f"  = {coeff_str}")
            break

print("\nDONE.")
