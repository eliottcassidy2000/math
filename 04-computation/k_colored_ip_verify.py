#!/usr/bin/env python3
"""
The k-colored independence polynomial I_k(Omega, x).

DEFINITION:
  I_k(Omega(T), x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * x^{parts(I)} * I(T)

where c_k^{(f,d)} is the inflated Eulerian coefficient.

KEY PROPERTIES:
  (1) I_k(Omega, 2) = a_k(T) for all k (forward-edge distribution)
  (2) I_0(Omega, x) = I(Omega, x) (standard independence polynomial)
  (3) I_{n-1}(Omega, x) = I(Omega, x) (palindromy)
  (4) For the transitive tournament: I_k(Omega, x) = A(n,k) for all x

The GENERATING FUNCTION over k:
  sum_k I_k(Omega, x) t^k = A_n(t) + sum_I x^{parts} I(T) A_{f+1}(t) (t-1)^{d-f}

SPECIAL x VALUES:
  x = 2: a_k(T) (known)
  x = 1: "reduced" forward-edge distribution
  x = 0: A(n,k) = Eulerian numbers (invariant-independent)
  x = -1: signed contribution

THE a_0 = EMPTY SET CONNECTION:
  I(Omega, x) = sum_{S independent} x^{|S|}
  The empty set S=emptyset contributes x^0 = 1.
  In the k-colored version, the "empty set" contributes A(n,k).
  So A(n,k) is the "k-colored empty-set weight."

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
# Part 1: Compute I_k(Omega, x) as polynomial in x
# ====================================================================
print("k-COLORED INDEPENDENCE POLYNOMIAL I_k(Omega, x)")
print("=" * 70)

for n in [5, 7]:
    print(f"\nn={n}:")
    d = n - 1

    if n == 5:
        invariants = [('t3', 2, 1), ('t5', 0, 1)]
    else:
        invariants = [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]

    # For a specific tournament, compute I_k as polynomial in x
    for seed in range(3):
        A = random_tournament(n, n * 500 + seed)
        dist = forward_edge_dist_dp(A, n) if n >= 7 else {k: v for k, v in
               ((k, sum(1 for p in permutations(range(n))
                       if sum(1 for i in range(n-1) if A[p[i]][p[i+1]]) == k))
               for k in range(n))}

        # Actually, let me use the DP version for both
        dist = forward_edge_dist_dp(A, n)

        inv_vals = {}
        inv_vals['t3'] = count_t3(A, n)
        if n >= 5:
            inv_vals['t5'] = count_directed_cycles(A, n, 5)
        if n >= 7:
            inv_vals['t7'] = count_directed_cycles(A, n, 7)
            inv_vals['bc'] = count_bc(A, n)

        H = dist.get(n-1, 0)
        alpha_1 = sum(inv_vals.get(f't{c}', 0) for c in [3, 5, 7])
        alpha_2 = inv_vals.get('bc', 0)

        print(f"\n  seed={seed}: inv={inv_vals}, H={H}, alpha_1={alpha_1}, alpha_2={alpha_2}")

        # Compute I_k(Omega, x) for each k
        # I_k(x) = A(n,k) + sum_I c_k^{(f,d)} * x^{parts} * I(T)
        # Group by x-power:
        #   coeff of x^0 = A(n,k)
        #   coeff of x^1 = sum_{I with parts=1} c_k * I(T) = sum_{odd cycle types} c_k * count
        #   coeff of x^2 = sum_{I with parts=2} c_k * I(T)

        for k in range(n):
            # x^0 coefficient: A(n,k)
            x0 = eulerian_number(n, k)
            # x^1 coefficient: sum of c_k * cycle_count for single cycles
            x1 = 0
            for name, f, parts in invariants:
                if parts == 1:
                    x1 += inflated_eulerian(f, d, k) * inv_vals[name]
            # x^2 coefficient: sum of c_k * cycle_count for pairs
            x2 = 0
            for name, f, parts in invariants:
                if parts == 2:
                    x2 += inflated_eulerian(f, d, k) * inv_vals[name]

            # Verify: I_k(2) should equal a_k
            I_k_at_2 = x0 + 2*x1 + 4*x2
            actual = dist.get(k, 0)

            poly_str = f"{x0} + {x1}x + {x2}x^2" if x2 != 0 else f"{x0} + {x1}x"

            if k == 0 or k == n-1 or seed == 0:
                print(f"    I_{k}(x) = {poly_str}  | I_k(2)={I_k_at_2}, a_k={actual}, {'OK' if I_k_at_2==actual else 'FAIL'}")

        # Verify I_0 = I_{n-1} = I(Omega, x)
        # I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2
        for k_check in [0, n-1]:
            x0 = eulerian_number(n, k_check)
            x1 = sum(inflated_eulerian(f, d, k_check) * inv_vals[name]
                     for name, f, parts in invariants if parts == 1)
            x2 = sum(inflated_eulerian(f, d, k_check) * inv_vals[name]
                     for name, f, parts in invariants if parts == 2)
            I_k_poly = (x0, x1, x2)
            I_omega = (1, alpha_1, alpha_2)
            print(f"    I_{k_check}(x) = ({x0}, {x1}, {x2}) vs I(Omega,x) = {I_omega}: {'MATCH' if I_k_poly == I_omega else 'DIFFER'}")

# ====================================================================
# Part 2: I_k(Omega, 1) — what does this count?
# ====================================================================
print(f"\n{'=' * 70}")
print("I_k(Omega, 1) — the 'unweighted' count")
print("=" * 70)

# I_k(Omega, 1) = A(n,k) + sum_I c_k * I(T)
# I(Omega, 1) = 1 + alpha_1 + alpha_2 = number of independent sets including empty set
# This is M(Omega) = matching number in some contexts

n = 7
for seed in range(10):
    A = random_tournament(n, n * 600 + seed)
    dist = forward_edge_dist_dp(A, n)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)

    alpha_1 = t3 + t5 + t7
    alpha_2 = bc
    I_omega_1 = 1 + alpha_1 + alpha_2  # total independent sets

    # I_k(1) for each k
    d = n - 1
    Ik_vals = []
    for k in range(n):
        val = eulerian_number(n, k)
        for name, f, parts in [('t3', 4, 1), ('t5', 2, 1), ('t7', 0, 1), ('bc', 2, 2)]:
            val += inflated_eulerian(f, d, k) * {'t3': t3, 't5': t5, 't7': t7, 'bc': bc}[name]
        Ik_vals.append(val)

    if seed < 5:
        print(f"\n  seed={seed}: I(Omega,1)={I_omega_1}")
        print(f"    I_k(1): {Ik_vals}")
        print(f"    Sum: {sum(Ik_vals)} (should be n!={factorial(n)})")
        # Check: is I_k(1) palindromic?
        pal = all(Ik_vals[k] == Ik_vals[n-1-k] for k in range(n))
        print(f"    Palindromic: {pal}")
        # Check I_0(1) = I_{n-1}(1) = I(Omega, 1)
        print(f"    I_0(1) = {Ik_vals[0]}, I_{n-1}(1) = {Ik_vals[n-1]}, I(Omega,1) = {I_omega_1}, match: {Ik_vals[0]==I_omega_1}")

# ====================================================================
# Part 3: Ratio structure — I_k(Omega, x) / I(Omega, x)
# ====================================================================
print(f"\n{'=' * 70}")
print("RATIO: I_k(Omega, x) / I(Omega, x) at x=2")
print("=" * 70)

n = 7
print(f"\nn={n}: ratio a_k / H(T):")
for seed in range(5):
    A = random_tournament(n, n * 700 + seed)
    dist = forward_edge_dist_dp(A, n)
    H = dist.get(n-1, 0)
    ratios = [dist.get(k, 0) / H if H > 0 else 0 for k in range(n)]
    print(f"  seed={seed}: H={H}, ratios={[f'{r:.3f}' for r in ratios]}")

# For the transitive tournament: a_k = A(n,k), H = 1
# ratios = A(n,k) / 1 = A(n,k) — very spiky
# For random: ratios should be smoother

print(f"\n  Transitive: ratios = {[eulerian_number(n, k) for k in range(n)]}")

# ====================================================================
# Part 4: The TRIVARIATE generating function
# ====================================================================
print(f"\n{'=' * 70}")
print("TRIVARIATE GF: G(t, x) = sum_k I_k(Omega, x) t^k")
print("=" * 70)

print("""
G_T(t, x) = A_n(t) + sum_I x^{parts(I)} * I(T) * A_{f_I+1}(t) * (t-1)^{n-1-f_I}

Special evaluations:
  G_T(t, 0) = A_n(t) for all T  [Eulerian poly, T-independent]
  G_T(t, 2) = E_T(t) = sum_k a_k t^k  [tournament Eulerian poly]
  G_T(0, x) = I_0(Omega, x) = I(Omega, x) = H at x=2
  G_T(1, x) = n! for all x  [since (t-1)^{d-f} = 0 at t=1]

The point (t, x) = (0, 2) gives H(T).
The point (t, x) = (1, any) gives n!.
The line x = 0 gives the Eulerian polynomial (T-independent).
The line x = 2 gives the tournament Eulerian polynomial.
""")

print(f"{'=' * 70}")
print("DONE")
print("=" * 70)
