#!/usr/bin/env python3
"""
Real-rootedness of the tournament Eulerian polynomial.

E_T(x) = sum_k a_k(T) x^k

where a_k = #{perms with k forward edges in T}.

For transitive T: E_T(x) = Eulerian polynomial A_n(x), which is real-rooted.

Question: Is E_T(x) real-rooted for ALL tournaments T?

Since a_k = a_{n-1-k} (palindromic), E_T(x) = x^{n-1} E_T(1/x).
So roots come in pairs (r, 1/r). All roots are negative for A_n(x).

Also: W(r) = E_T((r+1/2)/(r-1/2)) * (r-1/2)^{n-1} (up to scaling)
Wait, more precisely:
W(r) = sum_k a_k (r+1/2)^k (r-1/2)^{n-1-k}
     = (r-1/2)^{n-1} sum_k a_k ((r+1/2)/(r-1/2))^k
     = (r-1/2)^{n-1} E_T(p/q)  where p=r+1/2, q=r-1/2

So real roots of E_T give real roots of W(r) (via the Mobius transform r -> p/q).

opus-2026-03-07-S32
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
import random

def random_tournament(n, seed=42):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def forward_edge_dist_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

print("REAL-ROOTEDNESS OF TOURNAMENT EULERIAN POLYNOMIAL")
print("=" * 70)

for n in [5, 7, 9]:
    print(f"\nn={n}:")
    num_real = 0
    num_complex = 0
    num_tested = 0

    for trial in range(50 if n <= 7 else 20):
        A = random_tournament(n, n * 333 + trial)
        dist = forward_edge_dist_dp(A, n)
        coeffs = [dist.get(k, 0) for k in range(n)]

        # numpy polynomial: coeffs[0] + coeffs[1]*x + ... + coeffs[n-1]*x^{n-1}
        roots = np.roots(list(reversed(coeffs)))  # numpy wants highest degree first

        all_real = all(abs(r.imag) < 1e-6 for r in roots)
        all_negative = all(r.real < 0 for r in roots) if all_real else False

        if all_real:
            num_real += 1
        else:
            num_complex += 1
            if num_complex <= 3:
                print(f"  COMPLEX ROOTS at trial {trial}: roots = {sorted(roots, key=lambda r: r.real)}")
                print(f"    coeffs = {coeffs}")

        num_tested += 1

    print(f"  {num_real}/{num_tested} real-rooted, {num_complex} with complex roots")

    # Check: are all real roots negative?
    if n <= 7 and num_real > 0:
        all_neg = True
        for trial in range(50):
            A = random_tournament(n, n * 333 + trial)
            dist = forward_edge_dist_dp(A, n)
            coeffs = [dist.get(k, 0) for k in range(n)]
            roots = np.roots(list(reversed(coeffs)))
            if all(abs(r.imag) < 1e-6 for r in roots):
                if not all(r.real < 1e-6 for r in roots):
                    all_neg = False
                    break
        print(f"  All real roots negative: {all_neg}")

# Also check: log-concavity of a_k
print(f"\n{'=' * 70}")
print("LOG-CONCAVITY: a_k^2 >= a_{k-1}*a_{k+1}")
print("=" * 70)

for n in [5, 7, 9]:
    lc_pass = 0
    lc_fail = 0
    for trial in range(50 if n <= 7 else 20):
        A = random_tournament(n, n * 333 + trial)
        dist = forward_edge_dist_dp(A, n)
        coeffs = [dist.get(k, 0) for k in range(n)]

        lc = all(coeffs[k]**2 >= coeffs[k-1]*coeffs[k+1] for k in range(1, n-1))
        if lc:
            lc_pass += 1
        else:
            lc_fail += 1
            if lc_fail <= 2:
                print(f"  n={n}, trial {trial}: LC FAIL. coeffs={coeffs}")
                for k in range(1, n-1):
                    if coeffs[k]**2 < coeffs[k-1]*coeffs[k+1]:
                        print(f"    k={k}: {coeffs[k]}^2 = {coeffs[k]**2} < {coeffs[k-1]}*{coeffs[k+1]} = {coeffs[k-1]*coeffs[k+1]}")

    tested = 50 if n <= 7 else 20
    print(f"  n={n}: {lc_pass}/{tested} log-concave, {lc_fail} failures")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
