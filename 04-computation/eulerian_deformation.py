#!/usr/bin/env python3
"""
The forward-edge distribution a_k(T) as a deformation of Eulerian numbers.

For the transitive tournament: a_k = A(n, k).
For a general tournament T: a_k(T) = A(n, k) + correction(T, k).

Question: can we express correction(T, k) in terms of OCF invariants
and some "deformation Eulerian numbers"?

Since W(r) = sum_k a_k(T) * (r+1/2)^k * (r-1/2)^{n-1-k},
and W(r) = F_{n-1}(r) + sum_I 2^{parts} F_{f_I}(r) I(T),
we have:

sum_k a_k(T) * p^k * q^{n-1-k} = sum_k A(n,k) * p^k * q^{n-1-k}
                                  + sum_I 2^{parts} * [sum_j A(f_I+1, j) p^j q^{f_I-j}] * I(T)

So the correction to a_k from invariant I with free degree f_I is:
  delta_a_k(I) = 2^{parts(I)} * sum over ways to distribute k ascents among
                 the f_I free positions and the constrained positions.

Actually, more precisely: the I-term contributes 2^{parts} * F_{f_I}(r) * I(T)
to W(r). We need to expand F_{f_I}(r) in the p^k q^{n-1-k} basis (where p=r+1/2, q=r-1/2).

F_{f_I}(r) = sum_j A(f_I+1, j) p^j q^{f_I-j}

But this has degree f_I, while a_k involves degree n-1. So we need to express
F_{f_I}(r) * 1 in terms of p^k q^{n-1-k}... no, it's just that the
coefficients of W(r) in the p^k q^{n-1-k} basis give us a_k directly.

Let's think of it differently. W(r) = sum_k a_k p^k q^{n-1-k} where p=r+1/2, q=r-1/2.

Also W(r) = sum_{m=0}^{n-1} w_m r^m (the monomial expansion).

The connection between the two: p^k q^{n-1-k} = (r+1/2)^k (r-1/2)^{n-1-k}.

Expanding: p^k q^{n-1-k} = sum_{i+j = some} C(k,i) C(n-1-k,j) r^{i+j} (1/2)^{k-i} (-1/2)^{n-1-k-j}

This is a change-of-basis matrix. Let's compute it explicitly for small n.

opus-2026-03-07-S32
"""
from itertools import permutations, combinations
from collections import defaultdict
from fractions import Fraction
from math import comb, factorial

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def random_tournament(n, seed=42):
    import random
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def forward_edge_dist(A, n):
    dist = defaultdict(int)
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        dist[fwd] += 1
    return dict(dist)

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cycle_len):
    if n < cycle_len: return 0
    total = 0
    for verts in combinations(range(n), cycle_len):
        sub = [[A[verts[i]][verts[j]] for j in range(cycle_len)] for i in range(cycle_len)]
        dp = [[0]*cycle_len for _ in range(1 << cycle_len)]
        dp[1][0] = 1
        for m in range(1, 1 << cycle_len):
            for v in range(cycle_len):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cycle_len):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cycle_len) - 1
        total += sum(dp[full][v] for v in range(1, cycle_len) if sub[v][0])
    return total

# ====================================================================
# For each invariant I with free degree f, the contribution to a_k is:
# delta_a_k = 2^{parts} * [coefficient of p^k q^{n-1-k} in F_{f}(r)] * I(T)
#
# But F_f(r) = sum_j A(f+1, j) p^j q^{f-j}, which has degree f < n-1.
# We need to express p^j q^{f-j} in terms of p^k q^{n-1-k}.
# Since p - q = 1, we have p = q + 1.
# p^j q^{f-j} = sum_m C(j,m) q^m * q^{f-j} = sum_m C(j,m) q^{f-j+m}
# Hmm, that's the q-basis. We need the p^k q^{n-1-k} basis.
#
# Actually: F_f(r) lives in the polynomial ring in r of degree f.
# The a_k are the coefficients in the p^k q^{n-1-k} basis of degree n-1.
# So F_f(r) is a degree-f polynomial being expressed in a degree-(n-1) basis.
# This means the "coefficient of p^k q^{n-1-k} in F_f(r)" involves
# expressing r^m as a combination of p^k q^{n-1-k}.
#
# Let's do this numerically.
# ====================================================================

print("DEFORMATION OF EULERIAN NUMBERS BY INVARIANTS")
print("=" * 70)

for n in [5, 7]:
    print(f"\nn={n}:")

    # Compute for several tournaments, subtract Eulerian baseline
    baseline = {k: eulerian_number(n, k) for k in range(n)}
    print(f"  Baseline A({n},k): {baseline}")

    samples = []
    for trial in range(20):
        A = random_tournament(n, n*100 + trial)
        dist = forward_edge_dist(A, n)
        t3 = count_t3(A, n)
        t5 = count_directed_cycles(A, n, 5) if n >= 5 else 0

        delta = {k: dist.get(k, 0) - baseline[k] for k in range(n)}
        samples.append({'t3': t3, 't5': t5, 'delta': delta, 'dist': dist})

        if trial < 5:
            print(f"  T{trial}: t3={t3}, t5={t5}, delta={delta}")

    # For each k, regress delta_k on invariants
    if n == 5:
        print(f"\n  Regression: delta_k = c0 + c1*t3 + c2*t5")
        import numpy as np
        for k in range(n):
            X = np.array([[1, s['t3'], s['t5']] for s in samples])
            y = np.array([s['delta'][k] for s in samples])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            err = np.max(np.abs(y - X @ coeffs))
            c = [round(c) for c in coeffs]
            print(f"    k={k}: delta = {c[0]} + {c[1]}*t3 + {c[2]}*t5  (err={err:.6f})")

    if n == 7:
        def count_bc(A, n):
            cyc3 = [set(t) for t in combinations(range(n), 3)
                    if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
                       A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
            return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
                       if cyc3[i].isdisjoint(cyc3[j]))

        # Recompute with bc and t7
        for trial in range(20):
            A = random_tournament(n, n*100 + trial)
            t7 = count_directed_cycles(A, n, 7)
            bc = count_bc(A, n)
            samples[trial]['t7'] = t7
            samples[trial]['bc'] = bc

        print(f"\n  Regression: delta_k = c0 + c1*t3 + c2*t5 + c3*t7 + c4*bc")
        import numpy as np
        for k in range(n):
            X = np.array([[1, s['t3'], s['t5'], s['t7'], s['bc']] for s in samples])
            y = np.array([s['delta'][k] for s in samples])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            err = np.max(np.abs(y - X @ coeffs))
            c = [round(c) for c in coeffs]
            print(f"    k={k}: delta = {c[0]:5d} + {c[1]:4d}*t3 + {c[2]:4d}*t5 + {c[3]:4d}*t7 + {c[4]:4d}*bc  (err={err:.4f})")

        # Check palindromy of deltas: delta_k should equal delta_{n-1-k}
        print(f"\n  Palindromy of deltas:")
        for k in range(n//2 + 1):
            match = all(s['delta'][k] == s['delta'][n-1-k] for s in samples)
            print(f"    delta_{k} = delta_{n-1-k}: {match}")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
