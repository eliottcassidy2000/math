#!/usr/bin/env python3
"""
E_T(-1) as a deformation of Euler/zigzag numbers.

A_n(-1) = sum_k (-1)^k A(n,k) = E_n (signed Euler number):
  E_1=1, E_2=0, E_3=-2, E_4=0, E_5=16, E_6=0, E_7=-272, E_8=0, E_9=7936

(Even n: E_n=0 from palindromy. Odd n: |E_n| = tangent/secant numbers.)

E_T(-1) = A_n(-1) + sum_I 2^{parts} I(T) A_{f+1}(-1) (-2)^{d-f}

For odd n, this is nonzero and gives:
  n=3: E_T(-1) = -2 + 8*t3
  n=5: E_T(-1) = 16 - 16*t3 + 32*t5
  n=7: E_T(-1) = -272 + 128*t3 - 64*t5 + 128*t7 - 128*bc

QUESTION: Does E_T(-1) have a combinatorial interpretation?
  E_T(-1) = sum_P (-1)^{fwd(P)} where fwd(P) = #forward edges
  This is a SIGNED count of permutations by parity of forward edges.

Since fwd + bwd = n-1, and (-1)^{fwd} = (-1)^{n-1} (-1)^{bwd},
  for even n: E_T(-1) = -sum_P (-1)^{bwd(P)} = -E_T(-1), so = 0.
  for odd n: E_T(-1) = sum_P (-1)^{bwd(P)} = E_T(-1), consistent.

For odd n: E_T(-1) = sum_P (-1)^{fwd(P)}.

ALSO: E_T(i) where i = sqrt(-1)?
  E_T(i) = sum_k i^k a_k = sum of {1, i, -1, -i, ...} pattern.
  For odd n (n-1 even), a_k = a_{n-1-k}, and n-1 is even...

opus-2026-03-07-S33
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def eulerian_poly_eval(n, t):
    return sum(eulerian_number(n, k) * t**k for k in range(n))

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

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cl):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

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

# ====================================================================
# Part 1: Verify E_T(-1) formula
# ====================================================================
print("E_T(-1) = DEFORMED EULER/ZIGZAG NUMBER")
print("=" * 70)

# n=3: E_T(-1) = -2 + 8*t3
n = 3
print(f"\nn={n}: E_T(-1) = -2 + 8*t3")
for seed in range(4):
    A = random_tournament(n, seed)
    dist = forward_edge_dist_dp(A, n)
    actual = sum((-1)**k * dist.get(k, 0) for k in range(n))
    t3 = count_t3(A, n)
    pred = -2 + 8*t3
    print(f"  seed={seed}: t3={t3}, actual={actual}, pred={pred}, {'OK' if actual==pred else 'FAIL'}")

# n=5: E_T(-1) = 16 - 16*t3 + 32*t5
n = 5
print(f"\nn={n}: E_T(-1) = 16 - 16*t3 + 32*t5")
for seed in range(10):
    A = random_tournament(n, seed)
    dist = forward_edge_dist_dp(A, n)
    actual = sum((-1)**k * dist.get(k, 0) for k in range(n))
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    pred = 16 - 16*t3 + 32*t5
    ok = actual == pred
    if seed < 5 or not ok:
        print(f"  seed={seed}: t3={t3}, t5={t5}, actual={actual}, pred={pred}, {'OK' if ok else 'FAIL'}")

# n=7: E_T(-1) = -272 + 128*t3 - 64*t5 + 128*t7 - 128*bc
n = 7
print(f"\nn={n}: E_T(-1) = -272 + 128*t3 - 64*t5 + 128*t7 - 128*bc")
for seed in range(15):
    A = random_tournament(n, n*100 + seed)
    dist = forward_edge_dist_dp(A, n)
    actual = sum((-1)**k * dist.get(k, 0) for k in range(n))
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)
    pred = -272 + 128*t3 - 64*t5 + 128*t7 - 128*bc
    ok = actual == pred
    if seed < 5 or not ok:
        print(f"  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc={bc}, actual={actual}, pred={pred}, {'OK' if ok else 'FAIL'}")

# ====================================================================
# Part 2: Factor structure of E_T(-1) coefficients
# ====================================================================
print(f"\n{'=' * 70}")
print("COEFFICIENT STRUCTURE OF E_T(-1)")
print("=" * 70)

# E_T(-1) = A_n(-1) + sum_I 2^parts * I(T) * A_{f+1}(-1) * (-2)^{d-f}
# Let's compute the coefficient of each invariant:

for n in [3, 5, 7, 9]:
    d = n - 1
    print(f"\nn={n}: A_n(-1) = {eulerian_poly_eval(n, -1)}")

    if n == 3:
        invs = [('t3', 0, 1)]
    elif n == 5:
        invs = [('t3', 2, 1), ('t5', 0, 1)]
    elif n == 7:
        invs = [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]
    elif n == 9:
        invs = [('t3', 6, 1), ('t5', 4, 1), ('t7', 2, 1), ('t9', 0, 1),
                ('bc', 4, 2), ('bc35', 2, 2), ('bc37', 0, 2), ('a3', 2, 3)]

    for name, f, parts in invs:
        coeff = 2**parts * eulerian_poly_eval(f+1, -1) * (-2)**(d-f)
        print(f"  {name}: 2^{parts} * A_{f+1}(-1) * (-2)^{d-f} = {2**parts} * {eulerian_poly_eval(f+1,-1)} * {(-2)**(d-f)} = {coeff}")

# ====================================================================
# Part 3: Pattern recognition
# ====================================================================
print(f"\n{'=' * 70}")
print("PATTERN: Coefficients as powers of 2 * Euler numbers")
print("=" * 70)

# The coefficients are:
# n=3: const=-2, t3: 8
# n=5: const=16, t3: -16, t5: 32
# n=7: const=-272, t3: 128, t5: -64, t7: 128, bc: -128

# Let me check if there's a simpler expression.
# 2^parts * A_{f+1}(-1) * (-2)^{d-f}
# = 2^parts * E_{f+1} * (-1)^{d-f} * 2^{d-f}
# = (-1)^{d-f} * 2^{parts + d - f} * E_{f+1}
# where E_m = A_m(-1) is the signed Euler number.

# For single cycles (parts=1): (-1)^{d-f} * 2^{1+d-f} * E_{f+1}
#   d-f = n-1-f. For t_{2m+1}: f=d-2m, so d-f=2m.
#   Coeff = (-1)^{2m} * 2^{1+2m} * E_{d-2m+1} = 2^{1+2m} * E_{n-2m}

# Check: t3 at n=7: m=1, d=6, coeff = 2^3 * E_5 = 8 * 16 = 128. ✓
# t5 at n=7: m=2, d=6, coeff = 2^5 * E_3 = 32 * (-2) = -64. ✓
# t7 at n=7: m=3, d=6, coeff = 2^7 * E_1 = 128 * 1 = 128. ✓

# For pairs (parts=2): (-1)^{d-f} * 2^{2+d-f} * E_{f+1}
# bc at n=7: f=2, d-f=4. Coeff = 2^6 * E_3 = 64 * (-2) = -128. ✓

print("For single odd (2m+1)-cycle at vertex count n:")
print("  Coeff of t_{2m+1} in E_T(-1) = 2^{1+2m} * E_{n-2m}")
print("  where E_j = A_j(-1) is the signed Euler number.")
print()
print("For pairs of disjoint cycles with total size 2S:")
print("  Coeff = 2^{2+2S} * E_{n-2S}")
print()
print("General: for invariant I with parts(I) parts and total cycle vertices 2S:")
print("  Coeff = 2^{parts + 2S} * E_{n-2S}")
print()
print("And E_T(-1) = sum over all independent sets in Omega:")
print("  E_T(-1) = sum_{S independent} 2^{2|V(S)|} * E_{n-2|V(S)|} * product of cycle-count terms")
print()
print("Or more precisely:")
print("  E_T(-1) = E_n + sum_I 2^{parts(I)+2S_I} * E_{n-2S_I} * I(T)")
print("  where S_I = total cycle size (in pairs of vertices) for invariant I")

# ====================================================================
# Part 4: Can we express E_T(-1) as a deformed I(Omega, z)?
# ====================================================================
print(f"\n{'=' * 70}")
print("E_T(-1) AS DEFORMED INDEPENDENCE POLYNOMIAL?")
print("=" * 70)

# E_T(-1) = G_T(-1, 2) where G_T(t, x) is the trivariate GF.
# Can we write E_T(-1) = some_function(I(Omega, z)) for some z?

# The structure is: each invariant at level j gets coeff 2^{parts+2S} E_{n-2S}.
# For alpha_1 (single cycles): different cycles get different E factors.
# So E_T(-1) is NOT simply I(Omega, z) for any z.

# But: E_T(-1) / E_n = 1 + (coefficients / E_n) * invariants
# = 1 + sum_I [2^{parts+2S_I} E_{n-2S_I} / E_n] * I(T)

# This is like I(Omega, x) but with x_I = 2^{parts+2S_I} E_{n-2S_I} / E_n varying by level.

n = 7
E_n = eulerian_poly_eval(n, -1)
print(f"\nn={n}: E_n = A_{n}(-1) = {E_n}")

for name, total_cyc_verts, parts in [('t3', 2, 1), ('t5', 4, 1), ('t7', 6, 1), ('bc', 4, 2)]:
    S = total_cyc_verts // 2  # half of total vertices
    # Wait, S should be the total cycle vertex count / 2? Let me recheck.
    # Actually S is just the total size parameter. Let me use d-f directly.
    # d-f = n-1-f. For t_{2m+1}: d-f = 2m. Total cycle vertices = 2m+1.

    # The coeff is 2^{parts} * E_{f+1} * (-2)^{d-f}
    f = {'t3': 4, 't5': 2, 't7': 0, 'bc': 2}[name]
    coeff = 2**parts * eulerian_poly_eval(f+1, -1) * (-2)**(n-1-f)
    ratio = coeff / E_n if E_n != 0 else float('inf')
    print(f"  {name}: coeff={coeff}, ratio coeff/E_n = {coeff}/{E_n} = {Fraction(coeff, E_n)}")

# Check: is there a pattern in the ratios?
# t3: 128/(-272) = -8/17
# t5: -64/(-272) = 4/17
# t7: 128/(-272) = -8/17
# bc: -128/(-272) = 8/17
# Hmm, all are multiples of 1/17. And 17 = |E_7| / 16 = 272/16.
# Actually 272 = 16*17. So the ratios are ±8/17, ±4/17, ±8/17.
# Not particularly clean.

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
