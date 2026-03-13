#!/usr/bin/env python3
"""
alpha3_p7_only.py -- kind-pasteur-2026-03-13-S60

At p=7: only 80 directed cycles, full alpha decomposition feasible.
Key question: are alpha_3+ zero (making H trivially linear in moments)?
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


def enumerate_directed_cycles(A, p):
    cycles = []
    for k in range(3, p + 1, 2):
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

    # Start from -1 so vertex 0 is included in the enumeration
    backtrack(-1, 0, 0)
    return alpha


p = 7
m = 3
N = 8

print(f"{'='*70}")
print(f"FULL ALPHA DECOMPOSITION at p={p}")
print(f"{'='*70}")

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

    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    S4 = sum(d**4 for d in D_vals)
    S6 = sum(d**6 for d in D_vals)

    H = sum(alpha[j] * (2**j) for j in range(len(alpha)))
    max_j = max((j for j in range(len(alpha)) if alpha[j] > 0), default=0)

    print(f"\n  bits={bits}, S={S}, H={H}")
    print(f"    Cycles: {len(cycles)}, Max alpha idx: {max_j}")
    for j in range(max_j + 1):
        print(f"    alpha_{j} = {alpha[j]}", end="")
        if j > 0:
            print(f"  (contribution: {alpha[j]}*{2**j} = {alpha[j]*(2**j)})", end="")
        print()
    print(f"    S4 = {S4:.4f}, S6 = {S6:.4f}")
    print(f"    H = {' + '.join(f'{alpha[j]}*{2**j}' for j in range(max_j+1) if alpha[j]>0)}")

    data.append({
        'bits': bits, 'H': H, 'alpha': alpha[:max_j+1],
        'S4': S4, 'S6': S6
    })

# Summary
print(f"\n{'='*70}")
print(f"SUMMARY")
print(f"{'='*70}")

# Check if alpha_3+ are all zero
all_zero_3plus = True
for d in data:
    for j in range(3, len(d['alpha'])):
        if d['alpha'][j] > 0:
            all_zero_3plus = False
            print(f"  alpha_{j} = {d['alpha'][j]} at bits={d['bits']}")

if all_zero_3plus:
    print(f"\n  *** All alpha_j = 0 for j >= 3 ***")
    print(f"  H = 1 + 2*alpha_1 + 4*alpha_2 for ALL orientations at p=7")
    print(f"  Since alpha_1 = linear(S4,S6) [THM-157] and alpha_2 = linear(S4) [THM-156],")
    print(f"  H is AUTOMATICALLY a linear function of (S4, S6).")
    print(f"  (In fact, S6 is redundant since there are only 2 orbit types)")
else:
    print(f"\n  Higher alphas are present. H linearity depends on their moment structure.")

# Show the two orbit types
print(f"\n  Orbit type summary:")
by_H = defaultdict(list)
for d in data:
    by_H[d['H']].append(d)
for H_val in sorted(by_H.keys()):
    group = by_H[H_val]
    d = group[0]
    print(f"    H={H_val}: alpha={d['alpha']}, S4={d['S4']:.4f}, count={len(group)}")

print("\nDONE.")
