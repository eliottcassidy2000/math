#!/usr/bin/env python3
"""
Palindromic symmetry of G_T(t, x) and its consequences.

THEOREM (discovered this session):
  G_T(t, x) = t^{n-1} * G_T(1/t, x) for all x.

This means: for EACH x, the polynomial G_T(t, x) in t is palindromic
of degree n-1. So it can be written as t^{(n-1)/2} * P(t + 1/t, x)
for a polynomial P in (t + 1/t) and x.

PROOF: A_n(t) is palindromic: A_n(t) = t^{n-1} A_n(1/t).
A_{f+1}(t) is palindromic: A_{f+1}(t) = t^f A_{f+1}(1/t).
(t-1)^{d-f} under t -> 1/t: (1/t - 1)^{d-f} = (-1/t)^{d-f} (t-1)^{d-f}/t^{d-f}.
Since d-f is even: (-1)^{d-f} = 1, so (1/t-1)^{d-f} = (t-1)^{d-f}/t^{d-f}.

So A_{f+1}(1/t)(1/t-1)^{d-f} = t^f A_{f+1}(t)/t^f * (t-1)^{d-f}/t^{d-f}
  = A_{f+1}(t)(t-1)^{d-f} / t^{d}
  = A_{f+1}(t)(t-1)^{d-f} / t^{n-1}.

Therefore G_T(1/t, x) = A_n(1/t) + sum_I x^{parts} I(T) A_{f+1}(1/t)(1/t-1)^{d-f}
  = A_n(t)/t^{n-1} + sum_I x^{parts} I(T) A_{f+1}(t)(t-1)^{d-f}/t^{n-1}
  = G_T(t, x) / t^{n-1}.

QED.

CONSEQUENCE 1: The "k-colored independence polynomial" I_k(Omega, x) = [t^k] G_T(t, x)
satisfies I_k(Omega, x) = I_{n-1-k}(Omega, x). So the weight is "palindromic in k."

CONSEQUENCE 2: The t=i evaluation:
  G_T(i, x) = i^{n-1} G_T(1/i, x) = i^{n-1} G_T(-i, x).
  For n odd (n-1 even): G_T(i, x) = G_T(-i, x), meaning G_T(t, x) is EVEN in t
  when restricted to {i, -i}. This is consistent with E_T(i) being purely real/imaginary.

CONSEQUENCE 3: The roots of G_T(t, x) in t come in pairs (r, 1/r).
  For x=2: the roots of E_T(t) have this structure.
  For x=0: the roots of A_n(t) have this structure (well known).

This script explores these consequences and looks for new identities.

opus-2026-03-07-S33
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def inflated_eulerian(f, d, k):
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total

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
# Part 1: Verify palindromic symmetry of I_k(Omega, x)
# ====================================================================
print("PALINDROMIC SYMMETRY OF I_k(Omega, x)")
print("=" * 70)

print("""
THEOREM: I_k(Omega, x) = I_{n-1-k}(Omega, x) for all k and all x.

This follows from G_T(t, x) = t^{n-1} G_T(1/t, x), extracting the t^k coefficient.
""")

n = 7
d = n - 1
invariants = [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]

for seed in range(3):
    A = random_tournament(n, n * 3000 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }

    print(f"\n  seed={seed}: inv={inv_vals}")
    for k in range(n):
        # Compute I_k(Omega, x) as polynomial in x
        # Level 0: A(n,k)
        # Level 1: sum of 2*c_k * cycle_count
        # Level 2: sum of 4*c_k * pair_count
        x0 = eulerian_number(n, k)  # coefficient of x^0
        x1 = sum(2 * inflated_eulerian(f, d, k) * inv_vals[name]
                 for name, f, parts in invariants if parts == 1)
        x2 = sum(4 * inflated_eulerian(f, d, k) * inv_vals[name]
                 for name, f, parts in invariants if parts == 2)

        # I_k as (x0, x1, x2) polynomial in x
        I_k_poly = (x0, x1, x2)
        k_mirror = n - 1 - k
        if k < k_mirror:
            # Compute the mirror
            x0m = eulerian_number(n, k_mirror)
            x1m = sum(2 * inflated_eulerian(f, d, k_mirror) * inv_vals[name]
                      for name, f, parts in invariants if parts == 1)
            x2m = sum(4 * inflated_eulerian(f, d, k_mirror) * inv_vals[name]
                      for name, f, parts in invariants if parts == 2)
            I_km_poly = (x0m, x1m, x2m)
            match = I_k_poly == I_km_poly
            print(f"    I_{k} = ({x0}, {x1}x, {x2}x^2)  vs  I_{k_mirror} = ({x0m}, {x1m}x, {x2m}x^2)  {'MATCH' if match else 'FAIL'}")

# ====================================================================
# Part 2: The substitution u = t + 1/t
# ====================================================================
print(f"\n{'=' * 70}")
print("SUBSTITUTION: u = t + 1/t")
print("=" * 70)

print("""
Since G_T(t, x) = t^{(n-1)/2} * P(u, x) where u = t + 1/t (for odd n),
the polynomial P is of degree (n-1)/2 in u.

For n=7 (degree 6 in t, palindromic => degree 3 in u):
  G_T(t, x) = t^3 * P(u, x)
  P(u, x) = p_0(x) + p_1(x)*u + p_2(x)*u^2 + p_3(x)*u^3

The t^3 factor means: G_T(t, x) = p_3(x)*t^6 + ... + p_3(x)*t^0
which is palindromic as required.

Let me compute P(u, x) coefficients.
""")

n = 7
d = n - 1
m = d // 2  # = 3 for n=7

# For a palindromic polynomial sum_{k=0}^{d} a_k t^k with a_k = a_{d-k},
# write t^{d/2} P(u) where u = t + 1/t.
# a_k = sum_j p_j * [coefficient of t^{k-d/2} in u^j]
# u = t + 1/t, u^2 = t^2 + 2 + 1/t^2, u^3 = t^3 + 3t + 3/t + 1/t^3
# In t^{d/2} * u^j, the t^k coefficient is [t^{k-d/2} in u^j].

# For u^0 = 1: contributes to t^{d/2} only
# For u^1 = t + 1/t: contributes to t^{d/2+1} and t^{d/2-1}
# For u^2 = t^2 + 2 + 1/t^2: contributes to t^{d/2+2}, 2*t^{d/2}, t^{d/2-2}
# For u^3 = t^3 + 3t + 3/t + 1/t^3: contributes to t^{d/2+3}, 3*t^{d/2+1}, 3*t^{d/2-1}, t^{d/2-3}

# So for d=6, m=3:
# a_6 = p_3, a_5 = 3*p_3 + p_2 + ... wait, let me be more systematic.

# The Chebyshev-like expansion: for a palindromic polynomial of degree 2m,
# t^m * P(t + 1/t) = sum_k a_k t^k where P(u) = sum_j p_j u^j.
# The relation between a_k and p_j involves Chebyshev-type coefficients.

# Let's just compute directly for each tournament.

for seed in range(3):
    A = random_tournament(n, n * 3100 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    alpha_1 = inv_vals['t3'] + inv_vals['t5'] + inv_vals['t7']
    alpha_2 = inv_vals['bc']
    H = 1 + 2*alpha_1 + 4*alpha_2

    # Compute I_k for each k at x=2
    dist = forward_edge_dist_dp(A, n)

    # For palindromic (a_0, a_1, a_2, a_3, a_3, a_2, a_1, a_0) of degree 6:
    # a_0 = p_3 - 3*p_2 + ... hmm this is getting complicated.
    # Let me just substitute u-values.

    # P(u) = G_T(t, 2) / t^3 where u = t + 1/t.
    # At t=1: u=2, P(2) = n!/1 = 5040.
    # At t=-1: u=-2, P(-2) = G_T(-1, 2) / (-1)^3 = -E_T(-1)
    # At t=2: u=5/2, P(5/2) = G_T(2, 2) / 8

    E_at_2 = sum(dist.get(k, 0) * 2**k for k in range(n))
    E_at_m1 = sum(dist.get(k, 0) * (-1)**k for k in range(n))

    # Compute at several t values to determine P(u)
    # t = exp(i*theta): u = 2*cos(theta), |P| = |G| / |t|^3 = |G|
    if seed == 0:
        print(f"\n  seed={seed}: H={H}")
        print(f"    P(2) = n! = {factorial(n)} [t=1]")
        print(f"    P(-2) = -E_T(-1) = {-E_at_m1} [t=-1]")
        print(f"    P(5/2) = E_T(2)/8 = {Fraction(E_at_2, 8)} [t=2]")

        # Also compute P(0) = G_T(i, 2) / i^3 = G_T(i, 2) * i (since 1/i^3 = i)
        # But G_T(i, 2) involves complex numbers...
        # For u=0: t = i (since i + 1/i = i - i = 0)
        # P(0) = G_T(i, 2) / i^3 = G_T(i, 2) * i
        # Since for n=7 (n≡3 mod 4), E_T(i) is purely imaginary: E_T(i) = b*i
        # P(0) = b*i * i = -b

        # Compute E_T(i) at x=2
        E_at_i_real = sum(dist.get(k, 0) * [1, 0, -1, 0][k % 4] for k in range(n))
        E_at_i_imag = sum(dist.get(k, 0) * [0, 1, 0, -1][k % 4] for k in range(n))
        print(f"    E_T(i) = {E_at_i_real} + {E_at_i_imag}i")
        P_0 = -E_at_i_imag  # P(0) = -b where E_T(i) = bi
        print(f"    P(0) = {P_0} [t=i, u=0]")

        # Now we have P at u = -2, 0, 2, 5/2.
        # P is cubic in u: P(u) = p0 + p1*u + p2*u^2 + p3*u^3
        # Solve for coefficients:
        # P(-2) = p0 - 2*p1 + 4*p2 - 8*p3
        # P(0)  = p0
        # P(2)  = p0 + 2*p1 + 4*p2 + 8*p3
        # P(5/2)= p0 + 5/2*p1 + 25/4*p2 + 125/8*p3

        p0 = P_0
        # P(2) + P(-2) = 2*(p0 + 4*p2) => p2 = (P(2)+P(-2)-2*p0) / 8
        P2 = factorial(n)
        Pm2 = -E_at_m1
        p2 = Fraction(P2 + Pm2 - 2*p0, 8)
        # P(2) - P(-2) = 2*(2*p1 + 8*p3) = 4*p1 + 16*p3
        diff = P2 - Pm2
        # P(5/2) = p0 + 5/2*p1 + 25/4*p2 + 125/8*p3
        # 4*p1 + 16*p3 = diff
        # 5/2*p1 + 125/8*p3 = P(5/2) - p0 - 25/4*p2
        P52 = Fraction(E_at_2, 8)
        rhs2 = P52 - p0 - Fraction(25, 4) * p2
        # System: 4*p1 + 16*p3 = diff, 5/2*p1 + 125/8*p3 = rhs2
        # 4*p1 = diff - 16*p3, p1 = (diff - 16*p3)/4
        # 5/2 * (diff - 16*p3)/4 + 125/8*p3 = rhs2
        # 5*diff/8 - 20*p3 + 125/8*p3 = rhs2
        # (125/8 - 20)*p3 = rhs2 - 5*diff/8
        # (125/8 - 160/8)*p3 = rhs2 - 5*diff/8
        # -35/8 * p3 = rhs2 - 5*diff/8
        p3 = Fraction(rhs2 - 5*Fraction(diff, 8), Fraction(-35, 8))
        p1 = Fraction(diff - 16*p3, 4)

        print(f"\n    P(u) = {p0} + {p1}*u + {p2}*u^2 + {p3}*u^3")
        print(f"    Verify: P(2) = {p0 + 2*p1 + 4*p2 + 8*p3} = {P2}")
        print(f"    Verify: P(-2) = {p0 - 2*p1 + 4*p2 - 8*p3} = {Pm2}")
        print(f"    Verify: P(0) = {p0} = {P_0}")

# ====================================================================
# Part 3: What do the P-coefficients depend on?
# ====================================================================
print(f"\n{'=' * 70}")
print("P-COEFFICIENTS AS FUNCTIONS OF CYCLE COUNTS")
print("=" * 70)

n = 7
d = n - 1
for seed in range(10):
    A = random_tournament(n, n * 3200 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    dist = forward_edge_dist_dp(A, n)

    # Compute P(0), P(2), P(-2)
    E_m1 = sum(dist.get(k, 0) * (-1)**k for k in range(n))
    E_i_imag = sum(dist.get(k, 0) * [0, 1, 0, -1][k % 4] for k in range(n))
    E_2 = sum(dist.get(k, 0) * 2**k for k in range(n))

    P0 = -E_i_imag
    P2 = factorial(n)
    Pm2 = -E_m1
    P52 = Fraction(E_2, 8)

    p0 = Fraction(P0)
    p2 = Fraction(P2 + Pm2 - 2*P0, 8)
    rhs2 = P52 - p0 - Fraction(25, 4) * p2
    diff = P2 - Pm2
    p3 = (rhs2 - 5*Fraction(diff, 8)) / Fraction(-35, 8)
    p1 = Fraction(diff - 16*p3, 4)

    t3, t5, t7, bc = inv_vals['t3'], inv_vals['t5'], inv_vals['t7'], inv_vals['bc']
    if seed < 5:
        print(f"  seed={seed}: t3={t3},t5={t5},t7={t7},bc={bc} | p0={p0}, p1={p1}, p2={p2}, p3={p3}")

# ====================================================================
# Part 4: Roots of G_T(t, 2) and reciprocal pairing
# ====================================================================
print(f"\n{'=' * 70}")
print("ROOTS OF E_T(t) AND RECIPROCAL PAIRING")
print("=" * 70)

import numpy as np

n = 7
for seed in range(5):
    A = random_tournament(n, n * 3300 + seed)
    dist = forward_edge_dist_dp(A, n)
    coeffs = [dist.get(k, 0) for k in range(n)]

    # Reverse for numpy (highest power first)
    roots = np.roots(coeffs[::-1])

    if seed < 3:
        print(f"\n  seed={seed}: H={coeffs[0]}")
        print(f"    Roots of E_T(t):")
        for r in sorted(roots, key=lambda z: (abs(z.imag), z.real)):
            recip = 1/r if abs(r) > 1e-10 else float('inf')
            # Check if 1/r is also a root
            closest_dist = min(abs(r2 - recip) for r2 in roots)
            print(f"      {r.real:+.4f} {r.imag:+.4f}i | 1/r = {recip.real:+.4f} {recip.imag:+.4f}i | paired: {closest_dist < 0.001}")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
