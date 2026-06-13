#!/usr/bin/env python3
"""
Does the characteristic polynomial of A determine H(T)?

Verified EXACTLY at n=5 (0 ambiguous).
Test at n=6 (partial) and n=7 (sampled).

Also test: does the MULTISET of eigenvalues determine H?
(Same as char poly, since char poly = prod(x - lambda_i))
"""
import numpy as np
import random
from collections import defaultdict
from itertools import combinations

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

random.seed(42)

for n in [6, 7]:
    N = 2000 if n == 6 else 500
    cp_to_H = defaultdict(set)

    for _ in range(N):
        A = random_tournament(n)
        H = ham_path_count_dp(A, n)
        A_np = np.array(A, dtype=float)
        # Use rounded char poly coefficients as key
        coeffs = np.round(np.real(np.poly(A_np)), 4)
        key = tuple(coeffs)
        cp_to_H[key].add(H)

    ambiguous = {k: v for k, v in cp_to_H.items() if len(v) > 1}
    print(f"n={n}: {len(cp_to_H)} distinct char polys, {len(ambiguous)} ambiguous ({N} samples)")
    for k, v in sorted(ambiguous.items(), key=lambda x: min(x[1]))[:10]:
        print(f"  char_poly coeffs={k} -> H = {sorted(v)}")

# ========== ALSO: What about tr(A^k) for k=1,...,n? ==========
# The traces of A^k determine the char poly (Newton's identities).
# tr(A) = 0 always.
# tr(A^2) = sum_{i,j} A[i][j]*A[j][i] = 0 (tournament: A[i][j]*A[j][i]=0)
# tr(A^3) = 6*t3 (counts directed 3-cycles, each counted 3*2=6 times...
#   actually each directed 3-cycle contributes 3 to the trace, and each
#   3-vertex cycle comes in 2 orientations, so tr(A^3) = 6*t3? Let's verify.)

print(f"\n\n=== tr(A^k) and cycle counts ===")
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

for bits in range(0, 1 << m, 128):  # sample
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    A_np = np.array(A, dtype=float)

    # Count 3-cycles
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[i][k] and A[k][j] and A[j][i]: t3 += 1

    tr3 = np.trace(A_np @ A_np @ A_np)
    tr4 = np.trace(np.linalg.matrix_power(A_np, 4))

    print(f"  bits={bits}: t3={t3}, tr(A^3)={tr3:.0f}, tr(A^3)/3={tr3/3:.0f}, "
          f"tr(A^4)={tr4:.0f}")

# The key insight: tr(A^k) counts DIRECTED walks of length k.
# For tournaments: tr(A^2) = 0 (no 2-cycles since A[i][j]*A[j][i]=0)
# tr(A^3) = number of directed 3-cycles = 3 * t3
# (Each 3-cycle i->j->k->i contributes 1 to sum_{all orderings} A[i][j]*A[j][k]*A[k][i])
# Wait: tr(A^3) = sum_{i,j,k} A[i][j]*A[j][k]*A[k][i] = sum of products over all
# directed walks i->j->k->i. Each UNDIRECTED 3-cycle {i,j,k} can have at most
# 2 directed orientations, each contributing to the trace.
# With 3 starting points each: 2*3 = 6 walks per undirected 3-cycle?
# But tr counts over all starting points, so tr(A^3) = 2*3*t3? No...
# Actually there are 2 directed 3-cycles per undirected triple, each with
# 3 rotations, giving 6 per triple. But only directed 3-cycles contribute!
# A triple with a 3-cycle contributes either 3 (one direction) or 6 (both).
# But a tournament triple has either 0 or 1 directed 3-cycle (in one direction).
# Wait: a tournament on 3 vertices has exactly 1 directed 3-cycle iff it's a
# 3-cycle (t3 counts these). Each such cycle has 3 rotations.
# So tr(A^3) = 3 * (number of directed 3-cycles) = 3 * 2 * t3? No.
# Let me just check: for t3=1 at n=5, what is tr(A^3)?

# ========== KEY QUESTION: Does {tr(A^k) : k=1,...,n} determine H? ==========
print(f"\n\n=== Do traces tr(A^k) determine H? (n=5 exact) ===")
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

traces_to_H = defaultdict(set)
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = ham_path_count_dp(A, n)
    A_np = np.array(A, dtype=float)

    traces = tuple(round(np.trace(np.linalg.matrix_power(A_np, k)), 4) for k in range(1, n+1))
    traces_to_H[traces].add(H)

ambiguous = {k: v for k, v in traces_to_H.items() if len(v) > 1}
print(f"  {len(traces_to_H)} distinct trace vectors, {len(ambiguous)} ambiguous")
# Since traces determine char poly (Newton), this should give same result: 0 ambiguous

# ========== At n=5: tr(A^3) = 3*t3 (verify) and higher ===
print(f"\n  Traces at n=5:")
for bits in [0, 100, 500, 1023]:
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = ham_path_count_dp(A, n)
    A_np = np.array(A, dtype=float)
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[i][k] and A[k][j] and A[j][i]: t3 += 1

    traces = [round(np.trace(np.linalg.matrix_power(A_np, k)), 1) for k in range(1, 6)]
    print(f"  bits={bits}: H={H}, t3={t3}, traces={traces}")
