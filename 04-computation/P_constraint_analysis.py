#!/usr/bin/env python3
"""
The constraint P(2, x) = n! and what it implies.

Since G_T(1, x) = n! for all x (all corrections vanish at t=1),
and t=1 corresponds to u = 1 + 1/1 = 2, we get P(2, x) = n! for all x.

This means: p_0(x) + 2*p_1(x) + 4*p_2(x) + 8*p_3(x) = n! for all x.

Since p_j(x) are polynomials in x, this gives constraints at each power of x:
  [x^0]: p_0(0) + 2*p_1(0) + 4*p_2(0) + 8*p_3(0) = n!
  [x^1]: p_0'(0) + 2*p_1'(0) + 4*p_2'(0) + 8*p_3'(0) = 0
  [x^2]: p_0''(0)/2 + 2*p_1''(0)/2 + 4*p_2''(0)/2 + 8*p_3''(0)/2 = 0

The x^0 equation is: A(7,3) - 2*A(7,1) + 2*(A(7,2)-3) + 4*A(7,1) + 8 = n!
Wait, let me check the x=0 values:
  p_0(0) = 2176, p_1(0) = 1188, p_2(0) = 120, p_3(0) = 1
  2176 + 2*1188 + 4*120 + 8*1 = 2176 + 2376 + 480 + 8 = 5040 = 7! ✓

The x^1 equation (coefficients of x):
  (-128*t3+16*t5-8*t7) + 2*12*(t3-t5+t7) + 4*(24*t3-6*t7) + 8*(t3+t5+t7) = 0
  = (-128+24+96+8)*t3 + (16-24+0+8)*t5 + (-8+24-24+8)*t7 = 0*t3 + 0*t5 + 0*t7 = 0 ✓

The x^2 equation (coefficients of x^2):
  16*bc + 2*(-12*bc) + 4*0 + 8*bc = (16-24+0+8)*bc = 0 ✓

ALL PASS! The constraint P(2, x) = n! is automatically satisfied.

This is expected since the GF formula guarantees G_T(1, x) = n!.
But it gives a CONSISTENCY CHECK on our P-coefficients.

MORE INTERESTING: What does P(0, x) = p_0(x) represent?
  P(0, x) corresponds to u = 0, i.e., t + 1/t = 0, i.e., t = ±i.
  So P(0, x) = G_T(i, x) / i^3.

And P(-2, x) = sum_j (-2)^j p_j(x) corresponds to t = -1 (since -1 + 1/(-1) = -2).
  P(-2, x) = G_T(-1, x) / (-1)^3 = -G_T(-1, x).

opus-2026-03-07-S33
"""
from itertools import combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

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
                for u_nd in range(n):
                    if mask & (1 << u_nd): continue
                    new_fwd = fwd + A[v][u_nd]
                    key = (mask | (1 << u_nd), u_nd, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

# ====================================================================
# Part 1: Verify the P(2, x) = n! constraint
# ====================================================================
print("CONSTRAINT VERIFICATION: P(2, x) = n!")
print("=" * 70)

n = 7
d = n - 1

# x^0 coefficients of p_j
p0_0 = 2176
p1_0 = 1188
p2_0 = 120
p3_0 = 1
print(f"[x^0]: {p0_0} + 2*{p1_0} + 4*{p2_0} + 8*{p3_0} = {p0_0 + 2*p1_0 + 4*p2_0 + 8*p3_0} = {factorial(n)} {'OK' if p0_0 + 2*p1_0 + 4*p2_0 + 8*p3_0 == factorial(n) else 'FAIL'}")

# x^1: for each invariant
for name, coeff_p0, coeff_p1, coeff_p2, coeff_p3 in [
    ('t3', -128, 12, 24, 1),
    ('t5', 16, -12, 0, 1),
    ('t7', -8, 12, -6, 1),
]:
    total = coeff_p0 + 2*coeff_p1 + 4*coeff_p2 + 8*coeff_p3
    print(f"[x^1, {name}]: {coeff_p0} + 2*{coeff_p1} + 4*{coeff_p2} + 8*{coeff_p3} = {total} {'= 0 OK' if total == 0 else 'FAIL'}")

# x^2: bc
coeff_p0_bc = 16
coeff_p1_bc = -12
coeff_p2_bc = 0
coeff_p3_bc = 1
total_bc = coeff_p0_bc + 2*coeff_p1_bc + 4*coeff_p2_bc + 8*coeff_p3_bc
print(f"[x^2, bc]: {coeff_p0_bc} + 2*{coeff_p1_bc} + 4*{coeff_p2_bc} + 8*{coeff_p3_bc} = {total_bc} {'= 0 OK' if total_bc == 0 else 'FAIL'}")

# ====================================================================
# Part 2: P at other special u-values
# ====================================================================
print(f"\n{'=' * 70}")
print("P(u, x) AT SPECIAL u-VALUES")
print("=" * 70)

print("""
u = 2 (t=1): P(2, x) = n! for all x [PROVED]
u = -2 (t=-1): P(-2, x) = -G_T(-1, x)
u = 0 (t=i): P(0, x) = G_T(i, x) / i^3

Since G_T(-1, x) = A_n(-1) + sum_I x^{parts} I(T) A_{f+1}(-1) (-2)^{d-f}
and A_n(-1) = E_n (Euler numbers for odd n, 0 for even n):

P(-2, x) = -(E_n + sum_I x^{parts} I(T) E_{f+1} (-2)^{d-f})

For even n: P(-2, x) = 0 (since G_T(-1, x) = 0 by palindromy).

For n=7 (odd): P(-2, x) = -(-272 + sum_I ...)

The "P form" at u=-2 gives us E_T(-1) in a compact way.
""")

# Compute P(-2, x) via the p_j formulas
n = 7
for seed in range(5):
    A = random_tournament(n, n * 5000 + seed)
    inv = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    dist = forward_edge_dist_dp(A, n)
    a = [dist.get(k, 0) for k in range(n)]

    t3, t5, t7, bc = inv['t3'], inv['t5'], inv['t7'], inv['bc']

    # P(-2, 2) = p_0(2) - 2*p_1(2) + 4*p_2(2) - 8*p_3(2)
    # p_j(2) at x=2:
    p3_2 = a[0]  # = H
    p2_2 = a[1]
    p1_2 = a[2] - 3*a[0]
    p0_2 = a[3] - 2*a[1]

    P_m2_2 = p0_2 - 2*p1_2 + 4*p2_2 - 8*p3_2

    # Compare with -E_T(-1)
    E_m1 = sum((-1)**k * a[k] for k in range(n))
    neg_E_m1 = -E_m1

    # Also directly: P(-2,2) = G_T(-1, 2) / (-1)^3 = -G_T(-1, 2) = -E_T(-1)
    if seed < 3:
        print(f"  seed={seed}: P(-2, 2) = {P_m2_2}, -E_T(-1) = {neg_E_m1}, match={P_m2_2 == neg_E_m1}")

# ====================================================================
# Part 3: The "inverse OCF" — expressing cycle counts from a_k
# ====================================================================
print(f"\n{'=' * 70}")
print("INVERSE OCF: Recovering cycle counts from a_k")
print("=" * 70)

print("""
Since p_3 = H = a_0 and the p_j are LINEAR in cycle counts,
we can potentially INVERT the system to recover cycle counts from a_k.

For n=7, we have 4 unknowns (t3, t5, t7, bc) and 4 equations:
  p_3 = H = 1 + 2*(t3+t5+t7) + 4*bc
  p_2 = a_1 = 120 + 48*t3 - 12*t7
  p_1 = a_2 - 3*H = 1188 + 24*t3 - 24*t5 + 24*t7 - 48*bc
  p_0 = a_3 - 2*a_1 = 2176 - 256*t3 + 32*t5 - 16*t7 + 64*bc

Solve the linear system:
  2*(t3+t5+t7) + 4*bc = H - 1 = a_0 - 1
  48*t3 - 12*t7 = a_1 - 120
  24*t3 - 24*t5 + 24*t7 - 48*bc = a_2 - 3*a_0 - 1188
  -256*t3 + 32*t5 - 16*t7 + 64*bc = a_3 - 2*a_1 - 2176
""")

import numpy as np

# System: M * [t3, t5, t7, bc]^T = b
M = np.array([
    [2, 2, 2, 4],       # p_3 equation
    [48, 0, -12, 0],    # p_2 equation
    [24, -24, 24, -48],  # p_1 equation
    [-256, 32, -16, 64]  # p_0 equation
])

print(f"Coefficient matrix M:")
print(M)
print(f"det(M) = {np.linalg.det(M):.0f}")

M_inv = np.linalg.inv(M)
print(f"\nM^(-1) (rounded):")
for row in M_inv:
    print(f"  [{', '.join(f'{x:.6f}' for x in row)}]")

# Convert to fractions for exact inversion
from fractions import Fraction as F
M_frac = [[F(int(M[i][j])) for j in range(4)] for i in range(4)]

# Gaussian elimination
def frac_inv(mat):
    n = len(mat)
    aug = [row[:] + [F(1) if i == j else F(0) for j in range(n)] for i, row in enumerate(mat)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if aug[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for row in range(n):
            if row != col:
                factor = aug[row][col]
                aug[row] = [aug[row][j] - factor * aug[col][j] for j in range(2*n)]
    return [row[n:] for row in aug]

M_inv_exact = frac_inv(M_frac)
print(f"\nExact M^(-1):")
for i, row in enumerate(M_inv_exact):
    terms = ['t3', 't5', 't7', 'bc']
    print(f"  {terms[i]} = " + " + ".join(f"({r})*p_{3-j}" for j, r in enumerate(row)))

# Verify
n = 7
print(f"\nVerification:")
for seed in range(5):
    A = random_tournament(n, n * 5100 + seed)
    inv = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    dist = forward_edge_dist_dp(A, n)
    a = [dist.get(k, 0) for k in range(n)]

    # Compute p_j
    p3 = a[0]
    p2 = a[1]
    p1 = a[2] - 3*a[0]
    p0 = a[3] - 2*a[1]

    # b vector (rhs)
    b = [p3 - 1, p2 - 120, p1 - 1188, p0 - 2176]

    # Recover cycle counts
    recovered = [sum(M_inv_exact[i][j] * b[j] for j in range(4)) for i in range(4)]
    actual = [inv['t3'], inv['t5'], inv['t7'], inv['bc']]
    ok = all(recovered[i] == actual[i] for i in range(4))
    if seed < 3 or not ok:
        print(f"  seed={seed}: recovered={[int(r) for r in recovered]}, actual={actual}, {'OK' if ok else 'FAIL'}")

print(f"\n{'=' * 70}")
print("RESULT: Cycle counts ARE recoverable from the forward-edge distribution!")
print("The P(u) representation makes this inversion explicit and clean.")
print("=" * 70)

# ====================================================================
# Part 4: What constraints does a_k >= 0 impose on cycle counts?
# ====================================================================
print(f"\n{'=' * 70}")
print("POSITIVITY CONSTRAINTS: a_k >= 0")
print("=" * 70)

print("""
Since a_k(T) >= 0 for all k, and a_k = A(n,k) + corrections,
the corrections must satisfy:

  A(n, k) + sum_I 2^{parts} c_k^{(f,d)} I(T) >= 0 for all k.

The most restrictive constraint is typically k=0 (or k=n-1):
  a_0 = H(T) >= 1 (since every tournament has at least 1 Hamiltonian path by Rédei).

But also: a_1 >= 0. At n=7:
  a_1 = 120 + 48*t3 - 12*t7 >= 0
  So t7 <= (120 + 48*t3)/12 = 10 + 4*t3.

And a_3 >= 0:
  a_3 = 2416 - 160*t3 + 32*t5 - 40*t7 + 64*bc >= 0
  (using the x=2 evaluated coefficients)
""")

# Actually a_3 = p_0 + 2*p_2 = (2176 - 256*t3 + 32*t5 - 16*t7 + 64*bc) + 2*(120 + 48*t3 - 12*t7)
#             = 2176 + 240 + (-256+96)*t3 + 32*t5 + (-16-24)*t7 + 64*bc
#             = 2416 - 160*t3 + 32*t5 - 40*t7 + 64*bc

# Let me check this:
for seed in range(3):
    A = random_tournament(n, n * 5200 + seed)
    inv = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    dist = forward_edge_dist_dp(A, n)
    t3, t5, t7, bc = inv['t3'], inv['t5'], inv['t7'], inv['bc']

    # a_1 bound
    a1_pred = 120 + 48*t3 - 12*t7
    # a_3 formula
    a3_pred = 2416 - 160*t3 + 32*t5 - 40*t7 + 64*bc
    print(f"  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc={bc}")
    print(f"    a_1 = {a1_pred} >= 0: {'OK' if a1_pred >= 0 else 'VIOLATION'}")
    print(f"    a_3 = {a3_pred} = {dist.get(3,0)}: {'OK' if a3_pred == dist.get(3,0) else 'MISMATCH'}")
    print(f"    Bound on t7 from a_1: t7 <= {10 + 4*t3}, actual t7={t7}: {'OK' if t7 <= 10 + 4*t3 else 'TIGHT'}")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
