#!/usr/bin/env python3
"""
WHY does eigenvalue 1 appear in almost every M at n=5?

Observations from exhaustive n=5:
  H=1:  evals = {2, sqrt(2), 1, -sqrt(2), -2}  -- eigenvalue 1 present
  H=3:  evals include 1 (multiplicity 1 or 2)
  H=5:  evals include 1 or don't (two patterns)
  H=9:  two patterns, both include 1
  H=11: evals include 1 (multiplicity? no, 3.0 appears twice)
  H=13: no eigenvalue 1
  H=15: all eigenvalues = 3

Actually looking more carefully: for n odd, does det(M-I) = 0?
At n=5: det(M-I) for transitive = det(diag(0,-2,0,-2,0)-offdiag) = ?

Let me check this systematically.

Also: trace(M) = H. For n=5, tr(M) = H.
If M has eigenvalue 1, then tr(M-I) = H-5...
That doesn't directly constrain eigenvalue 1.

Maybe the question is: is (M-I)v = 0 for some natural vector v?

For the transitive tournament: M has eigenvalue 1 with eigenvector...
Let me find it.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def compute_M_diagonal(T, n, a):
    val = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            pos = list(perm).index(a)
            val += (-1)**pos
    return val

def compute_M_offdiag(T, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    val = 0
    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1)**len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})
        ea = 0
        if len(S_set) == 1:
            ea = 1
        else:
            for p in permutations(S_set):
                if p[-1] != a: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                ea += prod
        bb2 = 0
        if len(R_set) == 1:
            bb2 = 1
        else:
            for p in permutations(R_set):
                if p[0] != b: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                bb2 += prod
        val += sign * ea * bb2
    return val

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        M[a, a] = compute_M_diagonal(T, n, a)
        for b in range(a+1, n):
            M[a, b] = compute_M_offdiag(T, n, a, b)
            M[b, a] = M[a, b]
    return M

def tournament_from_bits(n, bits):
    T = {}
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]:
                T[(i,j)] = 1; T[(j,i)] = 0
            else:
                T[(i,j)] = 0; T[(j,i)] = 1
            idx += 1
    return T


# ============================================================
# n=5: eigenvalue analysis
# ============================================================
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

print("=" * 70)
print("n=5: Eigenvalue analysis of M")
print("=" * 70)

# Group by isomorphism class
seen_patterns = {}

for bits in range(1 << len(pairs)):
    b_list = [(bits >> k) & 1 for k in range(len(pairs))]
    T = tournament_from_bits(n, b_list)
    H = count_H(T, n)
    M = compute_M(T, n)
    evals = sorted(np.linalg.eigvalsh(M.astype(float)))[::-1]
    evals_key = tuple(round(e, 3) for e in evals)

    if evals_key not in seen_patterns:
        # Find eigenvector for eigenvalue closest to 1
        has_one = any(abs(e - 1.0) < 0.01 for e in evals)

        seen_patterns[evals_key] = {
            'H': H, 'M': M.copy(), 'T': T.copy(),
            'evals': evals, 'has_one': has_one
        }

print(f"\n  {len(seen_patterns)} distinct eigenvalue patterns\n")

for evals_key, info in sorted(seen_patterns.items(), key=lambda x: x[1]['H']):
    H = info['H']
    M = info['M']
    evals = info['evals']
    has_one = info['has_one']

    print(f"  H={H}: evals = {[round(e, 4) for e in evals]}")
    print(f"    has eigenvalue 1: {has_one}")

    if has_one:
        # Find eigenvector
        evals_full, evecs = np.linalg.eigh(M.astype(float))
        for i, e in enumerate(evals_full):
            if abs(e - 1.0) < 0.01:
                v = evecs[:, i]
                # Normalize to make largest component 1
                v = v / v[np.argmax(np.abs(v))]
                print(f"    eigenvector for lambda=1: {[round(x, 4) for x in v]}")

    # det(M-I)
    det_MI = round(np.linalg.det(M.astype(float) - np.eye(n)))
    print(f"    det(M-I) = {det_MI}")

    # What is (M-I) * [1,1,...,1]?
    ones = np.ones(n)
    Mv = M @ ones
    print(f"    M * [1,...,1] = {[int(round(x)) for x in Mv]}")
    print(f"    (M-I) * [1,...,1] = {[int(round(x)) for x in Mv - ones]}")

    # Score sequence
    scores = tuple(sorted(sum(info['T'].get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    print(f"    scores = {scores}")
    print()


# ============================================================
# DEEPER: What vector is the eigenvalue-1 eigenvector?
# ============================================================
print("\n" + "=" * 70)
print("Eigenvalue 1 eigenvector structure")
print("=" * 70)

# For the transitive tournament
T = {}
for i in range(n):
    for j in range(n):
        if i != j:
            T[(i,j)] = 1 if i < j else 0
M = compute_M(T, n)

print(f"\n  Transitive T_5:")
print(f"    M =")
for row in M:
    print(f"      {list(row)}")

evals_full, evecs = np.linalg.eigh(M.astype(float))
print(f"    Eigenvalues: {[round(e, 4) for e in sorted(evals_full)[::-1]]}")
print(f"    Eigenvectors (columns):")
for i in range(n):
    v = evecs[:, i]
    e = evals_full[i]
    print(f"      lambda={round(e, 4)}: {[round(x, 4) for x in v]}")


# ============================================================
# Does eigenvalue 1 persist at n=7?
# ============================================================
print("\n" + "=" * 70)
print("n=7: Does eigenvalue 1 persist?")
print("=" * 70)

import random
random.seed(777)

n = 7
pairs_7 = [(i,j) for i in range(n) for j in range(i+1, n)]

has_one_count = 0
no_one_count = 0

for trial in range(500):
    T = {}
    for (i,j) in pairs_7:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    M = compute_M(T, n)
    evals = sorted(np.linalg.eigvalsh(M.astype(float)))[::-1]

    has_one = any(abs(e - 1.0) < 0.1 for e in evals)
    if has_one:
        has_one_count += 1
    else:
        no_one_count += 1
        if no_one_count <= 3:
            H = count_H(T, n)
            scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
            print(f"  No eigenvalue 1: H={H}, scores={scores}")
            print(f"    evals = {[round(e, 2) for e in evals]}")

    if trial >= 50 and has_one_count + no_one_count >= 50:
        break

print(f"\n  {has_one_count + no_one_count} tested: {has_one_count} have eigenvalue 1, {no_one_count} don't")
print(f"  Fraction with eigenvalue 1: {has_one_count/(has_one_count+no_one_count):.3f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
