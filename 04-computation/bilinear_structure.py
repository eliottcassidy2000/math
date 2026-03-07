#!/usr/bin/env python3
"""
Bilinear structure of G_T(t, x).

G_T(t, x) has clean evaluations on both axes:
  G_T(t, 0) = A_n(t)    [Eulerian polynomial]
  G_T(0, x) = I(Omega, x) [independence polynomial]

QUESTION 1: Is G_T(t, x) the PRODUCT A_n(t) * I(Omega, x)?
  At (0,0): A_n(0) * I(Omega, 0) = 1 * 1 = 1 = G_T(0,0). ✓
  At (0,2): A_n(0) * I(Omega, 2) = 1 * H = H = G_T(0,2). ✓
  At (1,0): A_n(1) * I(Omega, 0) = n! * 1 = n! = G_T(1,0). ✓
  At (1,2): A_n(1) * I(Omega, 2) = n! * H. But G_T(1,2) = n!.
  So A_n(1) * I(Omega, 2) = n! * H ≠ n! unless H = 1.
  FAILS! G_T is NOT a simple product.

QUESTION 2: Is G_T(t, x) = A_n(t) * f(x) + g(t) * I(Omega, x) - h(t,x)?
  A separation of variables approach.

QUESTION 3: What is the "cross term" G_T(t,x) - A_n(t) * I(Omega, x)?
  Let Delta(t, x) = G_T(t, x) - A_n(t) * I(Omega, x).
  Delta vanishes on both axes: Delta(t, 0) = 0, Delta(0, x) = 0.
  So Delta = t * x * R(t, x) for some polynomial R.

QUESTION 4: Write G_T(t, x) / A_n(t) when A_n(t) ≠ 0.
  This gives a "normalized" function that equals 1 at x=0 and I(Omega, x)/1 at t=0.

QUESTION 5: The mixed partial d^2G/(dt dx) at (0,0).
  At (0,0): d/dt G = A'_n(0) at x=0 and something else at x≠0.
  The mixed partial measures the "interaction" between t and x.

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

def G_T_formula(n, inv_vals, t, x):
    """Compute G_T(t, x) via the trivariate GF formula."""
    d = n - 1
    result = sum(eulerian_number(n, k) * t**k for k in range(n))

    if n == 5:
        invariants = [('t3', 2, 1), ('t5', 0, 1)]
    elif n == 7:
        invariants = [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]
    else:
        return result

    for name, f, parts in invariants:
        val = inv_vals.get(name, 0)
        if val == 0: continue
        A_f1_t = sum(eulerian_number(f+1, j) * t**j for j in range(f+1))
        result += x**parts * val * A_f1_t * (t-1)**(d-f)

    return result

# ====================================================================
# Part 1: Is G_T = product?
# ====================================================================
print("IS G_T(t, x) = A_n(t) * I(Omega, x)?")
print("=" * 70)

n = 7
d = n - 1
for seed in range(3):
    A = random_tournament(n, n * 2000 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    alpha_1 = inv_vals['t3'] + inv_vals['t5'] + inv_vals['t7']
    alpha_2 = inv_vals['bc']
    H = 1 + 2*alpha_1 + 4*alpha_2

    print(f"\n  seed={seed}: alpha=({alpha_1},{alpha_2}), H={H}")

    for t, x in [(Fraction(1,2), 1), (2, 1), (Fraction(1,2), Fraction(1,2))]:
        An_t = sum(eulerian_number(n, k) * t**k for k in range(n))
        I_x = 1 + alpha_1 * x + alpha_2 * x**2
        product = An_t * I_x
        G_val = G_T_formula(n, inv_vals, t, x)
        print(f"    (t={t}, x={x}): G_T={G_val}, A_n*I={product}, ratio={Fraction(G_val, product) if product != 0 else 'div0'}")

# ====================================================================
# Part 2: Cross term Delta = G_T - A_n * I
# ====================================================================
print(f"\n{'=' * 70}")
print("CROSS TERM: Delta(t, x) = G_T(t, x) - A_n(t) * I(Omega, x)")
print("=" * 70)

print("""
Since Delta(t, 0) = 0 and Delta(0, x) = 0, we can write:
  Delta(t, x) = sum of terms involving t * x * ...

From the GF: G_T(t, x) = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t) (t-1)^{d-f}
And: A_n(t) * I(Omega, x) = A_n(t) * (1 + sum_I x^{parts} I(T))
   = A_n(t) + A_n(t) * sum_I x^{parts} I(T)

So: Delta = sum_I x^{parts} I(T) [A_{f+1}(t) (t-1)^{d-f} - A_n(t)]

Each invariant I contributes: x^{parts} I(T) * [A_{f+1}(t)(t-1)^{d-f} - A_n(t)]

The factor [A_{f+1}(t)(t-1)^{d-f} - A_n(t)] vanishes at t=0:
  A_{f+1}(0)(0-1)^{d-f} - A_n(0) = 1 * (-1)^{d-f} - 1 = 0 (since d-f is even)
And at t=1:
  A_{f+1}(1)(1-1)^{d-f} - A_n(1) = (f+1)! * 0 - n! = -n! (for d-f > 0)
  So Delta(1, x) = -n! * sum_I x^{parts} I(T) = -n! * (I(Omega, x) - 1)
  Thus G_T(1, x) = A_n(1) * I(Omega, x) - n! * (I(Omega, x) - 1)
                  = n! * I(Omega, x) - n! * I(Omega, x) + n! = n!  ✓
""")

n = 7
d = n - 1
for seed in range(3):
    A = random_tournament(n, n * 2100 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    alpha_1 = inv_vals['t3'] + inv_vals['t5'] + inv_vals['t7']
    alpha_2 = inv_vals['bc']

    # Compute the "correction factor" for each invariant type
    for t_val in [Fraction(1,3), Fraction(1,2), Fraction(2,1)]:
        An_t = sum(eulerian_number(n, k) * t_val**k for k in range(n))
        I_omega_x2 = 1 + 2*alpha_1 + 4*alpha_2
        G_val = G_T_formula(n, inv_vals, t_val, 2)

        # Delta = G - An * I(Omega, 2)
        Delta = G_val - An_t * I_omega_x2
        if seed == 0:
            print(f"  t={t_val}: G_T(t,2)={G_val}, A_n(t)*H={An_t*I_omega_x2}, Delta={Delta}")

# ====================================================================
# Part 3: The "interaction polynomial" per invariant
# ====================================================================
print(f"\n{'=' * 70}")
print("INTERACTION POLYNOMIAL: R_f(t) = A_{f+1}(t)(t-1)^{d-f} - A_n(t)")
print("=" * 70)

print("""
For each invariant I at level f, the contribution to Delta is:
  x^{parts} I(T) * R_f(t) where R_f(t) = A_{f+1}(t)(t-1)^{d-f} - A_n(t)

R_f(t) vanishes at t=0 (since d-f even) and equals -n! at t=1.

Let's compute R_f(t) / t for the first few values:
""")

n = 7
d = n - 1
for f, label in [(4, 't3'), (2, 't5'), (0, 't7')]:
    print(f"\n  f={f} ({label}): R_f(t) = A_{f+1}(t)(t-1)^{d-f} - A_7(t)")
    for t_val in [Fraction(0), Fraction(1,4), Fraction(1,2), Fraction(3,4), 1, 2, -1]:
        A_f1 = sum(eulerian_number(f+1, j) * t_val**j for j in range(f+1))
        An = sum(eulerian_number(n, k) * t_val**k for k in range(n))
        R_f = A_f1 * (t_val - 1)**(d-f) - An
        R_f_over_t = Fraction(R_f, t_val) if t_val != 0 else "limit"
        print(f"    t={t_val}: R_f={R_f}, R_f/t={R_f_over_t}")

# ====================================================================
# Part 4: Normalized form G_T / A_n (when A_n != 0)
# ====================================================================
print(f"\n{'=' * 70}")
print("NORMALIZED: G_T(t, x) / A_n(t)")
print("=" * 70)

print("""
Define N_T(t, x) = G_T(t, x) / A_n(t) (well-defined when A_n(t) != 0).

N_T(t, 0) = 1 for all t.
N_T(0, x) = I(Omega, x) / 1 = I(Omega, x).
N_T(1, x) = n! / n! = 1. (Remarkable: N at t=1 is T-independent!)

So N_T(t, x) goes from I(Omega, x) at t=0 to 1 at t=1.
The t-direction "flattens" the independence polynomial to a constant!

N_T(t, 2) = E_T(t) / A_n(t) goes from H at t=0 to 1 at t=1.
""")

n = 7
d = n - 1
for seed in range(3):
    A = random_tournament(n, n * 2200 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    alpha_1 = inv_vals['t3'] + inv_vals['t5'] + inv_vals['t7']
    alpha_2 = inv_vals['bc']
    H = 1 + 2*alpha_1 + 4*alpha_2

    print(f"\n  seed={seed}: H={H}")
    for t_val in [Fraction(0), Fraction(1,4), Fraction(1,2), Fraction(3,4)]:
        An_t = sum(eulerian_number(n, k) * t_val**k for k in range(n))
        if An_t == 0:
            print(f"    t={t_val}: A_n(t) = 0!")
            continue
        G_2 = G_T_formula(n, inv_vals, t_val, 2)
        N = Fraction(G_2, An_t)
        print(f"    N(t={t_val}, 2) = {N} = {float(N):.4f}")

# ====================================================================
# Part 5: The "smoothing" interpretation
# ====================================================================
print(f"\n{'=' * 70}")
print("SMOOTHING: t moves H toward n!/n! = 1")
print("=" * 70)

print("""
N_T(t, 2) = E_T(t) / A_n(t) interpolates from H(T) to 1.

This means: as t -> 1, the forward-edge distribution a_k(T) approaches
A(n,k) (the Eulerian numbers), up to a constant factor.

More precisely: a_k(T) = A(n,k) + corrections, and the corrections
are weighted by (t-1)^{d-f} in the GF, so they vanish at t=1.

The "speed" at which N_T(t, 2) approaches 1 near t=1 depends on the
smallest d-f among active invariants. For n=7:
  - t7: d-f = 6 (vanishes fastest)
  - t5, bc: d-f = 4
  - t3: d-f = 2 (vanishes slowest)

So near t=1, the 3-cycle correction dominates!
This means: the tournament "looks most different" from transitive
in its 3-cycle structure.
""")

# Compute the derivative dN/dt at t=1
# N(t, 2) = G_T(t, 2) / A_n(t)
# dN/dt|_{t=1} = [G'(1,2)*A(1) - G(1,2)*A'(1)] / A(1)^2
#              = [G'(1,2)*n! - n!*A'_n(1)] / (n!)^2
#              = [G'(1,2) - A'_n(1)] / n!

# G_T(t, 2) = A_n(t) + sum_I 2^parts I(T) A_{f+1}(t) (t-1)^{d-f}
# G'(t, 2) = A'_n(t) + sum_I 2^parts I(T) [A'_{f+1}(t)(t-1)^{d-f} + A_{f+1}(t)(d-f)(t-1)^{d-f-1}]
# At t=1: (t-1)^{d-f} = 0 for d-f >= 1, and (d-f)(t-1)^{d-f-1} = 0 for d-f >= 2.
# So G'(1, 2) = A'_n(1) + sum_{I with d-f=1} 2^parts I(T) A'_{f+1}(1) * 1 * 1
#             + sum_{I with d-f=1} 2^parts I(T) A_{f+1}(1) * 1

# But d-f = 1 means the cycle uses only 1 extra position... for odd cycles, d-f is always even.
# So NO invariant has d-f = 1! Therefore G'(1, 2) = A'_n(1), and dN/dt|_{t=1} = 0.

# What about the second derivative?
# For d-f = 2 (which is the 3-cycle's d-f at n=5, 7, ...):
# (t-1)^{d-f} at d-f=2: second derivative is 2
# So d^2 N/dt^2 at t=1 involves the d-f=2 invariants.

print("At t=1: dN/dt = 0 (all corrections have d-f >= 2, which is always even)")
print("d^2N/dt^2 at t=1 is dominated by d-f=2 terms (3-cycles for n>=5)")

n = 7
d = n - 1
nfact = factorial(n)
for seed in range(3):
    A = random_tournament(n, n * 2300 + seed)
    inv_vals = {
        't3': count_t3(A, n),
        't5': count_directed_cycles(A, n, 5),
        't7': count_directed_cycles(A, n, 7),
        'bc': count_bc(A, n)
    }
    # d-f=2 terms: t3 (f=4, d-f=2)
    # Contribution to G''(1, 2) from t3:
    # 2 * 2 * t3 * [A''_5(1) * 0 + 2 * A'_5(1) * 1 + A_5(1) * 2]
    # Wait, need to be more careful with the product rule.
    # Let h(t) = A_{f+1}(t) * (t-1)^{d-f}
    # h''(t) = A''(t)(t-1)^{d-f} + 2A'(t)(d-f)(t-1)^{d-f-1} + A(t)(d-f)(d-f-1)(t-1)^{d-f-2}
    # At t=1 with d-f=2:
    # h''(1) = A''(1)*0 + 2*A'(1)*2*0 + A(1)*2*1*1 = 2*(f+1)!

    # So the d-f=2 contribution to G''(1,2):
    # sum_{I with d-f=2} 2^parts * I(T) * 2 * (f+1)!
    f = 4  # for t3
    contrib_t3 = 2 * inv_vals['t3'] * 2 * factorial(f+1)
    # A'_n(1) = ? We need this too.
    # Actually, A_n(t) = sum_k A(n,k) t^k, so A'_n(1) = sum_k k*A(n,k)
    An_prime_1 = sum(k * eulerian_number(n, k) for k in range(n))
    An_second_1 = sum(k * (k-1) * eulerian_number(n, k) for k in range(n))
    # A'_n(1) = sum k*A(n,k) = (sum of descent counts) = n! * (n-1)/2
    # because the average number of descents in S_n is (n-1)/2.

    if seed == 0:
        print(f"\n  A'_n(1) = {An_prime_1}, expected n!*(n-1)/2 = {nfact * (n-1) // 2}")
        print(f"  Contribution of t3 to G''(1,2): {contrib_t3}")
        print(f"  t3 = {inv_vals['t3']}")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
