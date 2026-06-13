#!/usr/bin/env python3
"""
ROOTS OF F(T,x) AND LOG-CONCAVITY

The Eulerian polynomial A_n(x) has only real, negative roots (Frobenius).
This implies the coefficients A(n,k) are log-concave and unimodal.

QUESTION: Does F(T,x) also have real roots? Negative roots?

If F has only real negative roots, then F_k would be log-concave:
  F_k^2 >= F_{k-1} * F_{k+1}

This would be a strong structural constraint on HP ascent distributions.

Also explore: the "gamma-expansion" of F(T,x).
For palindromic polynomials with nonneg coefficients and only real roots,
there's a gamma-expansion: p(x) = sum gamma_k x^k (1+x)^{d-2k}.
Since F is NOT palindromic, this doesn't directly apply.

But F(T,x) + x^{n-1}F(T,1/x) IS palindromic (it's F(T,x) + F(T^op,x)).
Interesting angle to explore.
"""
from itertools import permutations
from math import comb, factorial
import numpy as np
from collections import defaultdict, Counter

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def forward_edge_poly(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        asc = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[asc] += 1
    return F

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

# ===== SECTION 1: Roots of F(T,x) =====
print("=" * 70)
print("ROOTS OF F(T,x)")
print("=" * 70)

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}, degree d={d}:")

    all_real_count = 0
    all_neg_count = 0
    has_zero_root = 0
    total = 0
    root_data = []

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)

        # Polynomial coefficients (highest degree first for np.roots)
        coeffs = list(reversed(F))
        roots = np.roots(coeffs)

        all_real = all(abs(r.imag) < 1e-8 for r in roots)
        all_neg = all_real and all(r.real < 1e-8 for r in roots)
        has_zero = F[0] == 0  # constant term is 0 => x=0 is a root

        if all_real: all_real_count += 1
        if all_neg: all_neg_count += 1
        if has_zero: has_zero_root += 1
        total += 1

        root_data.append((H, F, roots, all_real, all_neg))

        if total <= 5:
            real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-8])
            complex_roots = [(r.real, r.imag) for r in roots if abs(r.imag) > 1e-8]
            print(f"  F={F}, H={H}")
            print(f"    real roots: {[f'{r:.4f}' for r in real_roots]}")
            if complex_roots:
                print(f"    complex roots: {[(f'{r:.4f}', f'{i:.4f}') for r,i in complex_roots]}")
            print(f"    all real={all_real}, all neg={all_neg}")

    print(f"\n  SUMMARY: all real roots: {all_real_count}/{total}")
    print(f"  all negative real: {all_neg_count}/{total}")
    print(f"  has x=0 root (F_0=0): {has_zero_root}/{total}")

    # For those with all real roots, are they all negative?
    if all_real_count > all_neg_count:
        print(f"  {all_real_count - all_neg_count} have all real but some positive roots")

# ===== SECTION 2: Log-concavity =====
print("\n\n" + "=" * 70)
print("LOG-CONCAVITY OF F_k")
print("=" * 70)
print("Check F_k^2 >= F_{k-1} * F_{k+1} for all k")

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    lc_count = 0
    total = 0

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)

        is_lc = True
        for k in range(1, d):
            if F[k]**2 < F[k-1] * F[k+1]:
                is_lc = False
                break

        if is_lc: lc_count += 1
        total += 1

    print(f"  Log-concave: {lc_count}/{total}")

# ===== SECTION 3: Unimodality =====
print("\n\n" + "=" * 70)
print("UNIMODALITY OF F_k")
print("=" * 70)

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    unimodal_count = 0
    total = 0

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)

        # Unimodal: F increases then decreases
        peak = max(range(n), key=lambda k: F[k])
        is_uni = all(F[k] <= F[k+1] for k in range(peak)) and \
                 all(F[k] >= F[k+1] for k in range(peak, d))
        if is_uni: unimodal_count += 1
        total += 1

    print(f"  Unimodal: {unimodal_count}/{total}")

# ===== SECTION 4: The palindromic combination F + F^op =====
print("\n\n" + "=" * 70)
print("G(T,x) = F(T,x) + x^{n-1}F(T,1/x) = F(T,x) + F(T^op,x)")
print("=" * 70)
print("G is palindromic by construction. Study its properties.")

for n in [4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    all_real_G = 0
    gamma_nonneg = 0
    total = 0

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        F_op = list(reversed(F))  # F_{n-1-k}(T) = F_k(T^op)
        G = [F[k] + F_op[k] for k in range(n)]  # G_k = F_k + F_{d-k}

        # G is palindromic: G_k = G_{d-k}
        assert all(G[k] == G[d-k] for k in range(n))

        # Roots of G
        coeffs_G = list(reversed(G))
        roots_G = np.roots(coeffs_G)
        all_real = all(abs(r.imag) < 1e-8 for r in roots_G)
        if all_real: all_real_G += 1

        # Gamma expansion: G(x) = sum gamma_k x^k (1+x)^{d-2k}
        # For d=3 (n=4): G = gamma_0(1+x)^3 + gamma_1 x(1+x)
        # Expanding: gamma_0 + 3*gamma_0*x + 3*gamma_0*x^2 + gamma_0*x^3
        #          + gamma_1*x + gamma_1*x^2
        # G_0 = gamma_0
        # G_1 = 3*gamma_0 + gamma_1
        # By palindromicity, G_2 = G_1, G_3 = G_0 (automatically)
        if d == 3:
            gamma_0 = G[0]
            gamma_1 = G[1] - 3*G[0]
            gammas = [gamma_0, gamma_1]
        elif d == 4:
            # G(x) = gamma_0(1+x)^4 + gamma_1 x(1+x)^2 + gamma_2 x^2
            # G_0 = gamma_0
            # G_1 = 4*gamma_0 + gamma_1
            # G_2 = 6*gamma_0 + 2*gamma_1 + gamma_2
            gamma_0 = G[0]
            gamma_1 = G[1] - 4*G[0]
            gamma_2 = G[2] - 6*G[0] - 2*gamma_1
            gammas = [gamma_0, gamma_1, gamma_2]
        else:
            gammas = None

        if gammas and all(g >= 0 for g in gammas):
            gamma_nonneg += 1

        total += 1

        if total <= 5:
            print(f"  G={G}, gammas={gammas}")
            if gammas:
                print(f"    gamma nonneg: {all(g >= 0 for g in gammas)}")

    print(f"\n  G all real roots: {all_real_G}/{total}")
    print(f"  G gamma-nonneg: {gamma_nonneg}/{total}")

# ===== SECTION 5: The anti-palindromic part D(T,x) = F(T,x) - F(T^op,x) =====
print("\n\n" + "=" * 70)
print("D(T,x) = F(T,x) - F(T^op,x) [ANTI-PALINDROMIC PART]")
print("=" * 70)
print("D_k = F_k - F_{d-k}, so D_k + D_{d-k} = 0 (anti-palindromic)")
print("F = (G + D)/2 where G is palindromic and D is anti-palindromic")

for n in [4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    # What determines D?
    all_D = []
    all_data = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        D = [F[k] - F[d-k] for k in range(n)]
        H = sum(F)

        t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

        all_D.append(D)
        all_data.append((H, t3, F, D))

        if len(all_data) <= 5:
            print(f"  F={F}, D={D}, H={H}")

    # D(1) = F(1) - F^op(1) = H - H = 0. Good.
    # D has n independent components but D_k = -D_{d-k}, so floor(n/2) independent.

    # For n=4 (d=3): D_0 = -D_3, D_1 = -D_2. So 2 independent: D_0, D_1.
    # For n=5 (d=4): D_0 = -D_4, D_1 = -D_3, D_2 = 0. So 2 independent: D_0, D_1.

    print(f"\n  D at x=1: {[sum(D) for D in all_D[:3]]} (should all be 0)")

    # D_0 = F_0 - F_{d} = #{decreasing HPs} - #{increasing HPs}
    D0_vals = Counter([D[0] for D in all_D])
    print(f"\n  D_0 = F_0 - F_d distribution: {dict(sorted(D0_vals.items()))}")

    if n == 5:
        # D_0 and D_1 characterize the anti-palindromic part
        D01_counter = Counter([(D[0], D[1]) for D in all_D])
        print(f"  (D_0, D_1) has {len(D01_counter)} distinct values")

        # Is D determined by some tournament invariant?
        by_D01 = defaultdict(list)
        for H, t3, F, D in all_data:
            by_D01[(D[0], D[1])].append((H, t3))

        # For each (D_0, D_1), what H values appear?
        print(f"\n  Does (D_0, D_1) determine H?")
        for d01 in sorted(by_D01.keys()):
            Hs = set(h for h, _ in by_D01[d01])
            if len(Hs) > 1:
                print(f"    D_0={d01[0]}, D_1={d01[1]}: H in {sorted(Hs)}")
                break
        else:
            print(f"    Yes! (D_0, D_1) determines H for all values")

# ===== SECTION 6: F(T,x) for specific tournament families =====
print("\n\n" + "=" * 70)
print("F(T,x) FOR SPECIFIC TOURNAMENT FAMILIES")
print("=" * 70)

# 1. Transitive tournaments
print("\n--- Transitive tournaments ---")
for n in [3, 4, 5, 6, 7]:
    # Transitive: A[i][j] = 1 iff i < j
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1

    if n <= 5:
        F = forward_edge_poly(A, n)
    else:
        # DP version
        F = [0] * n
        dp = {}
        for v in range(n):
            dp[(1 << v, v, 0)] = 1  # (mask, last_vertex, ascent_count)

        for mask in range(1, 1 << n):
            for v in range(n):
                for a in range(n):
                    if (mask, v, a) not in dp: continue
                    if not (mask & (1 << v)): continue
                    for u in range(n):
                        if mask & (1 << u): continue
                        if A[v][u]:
                            new_a = a + (1 if v < u else 0)
                            key = (mask | (1 << u), u, new_a)
                            dp[key] = dp.get(key, 0) + dp[(mask, v, a)]

        full = (1 << n) - 1
        for v in range(n):
            for a in range(n):
                if (full, v, a) in dp:
                    F[a] += dp[(full, v, a)]

    H = sum(F)
    print(f"  n={n}: F={F}, H={H}")
    # For transitive tournament, ALL edges go from lower to higher label.
    # Every HP is a permutation where each step goes along a tournament edge (i<j => i->j).
    # So ALL steps are "forward" (v->u with A[v][u]=1, and since it's transitive, v<u).
    # Wait no: the HP visits vertices in some order, and v<u iff there's an edge v->u.
    # But the HP visits v then u where v->u, so v < u. So ALL edges are ascents.
    # Hence F = [0, 0, ..., 0, 1] (one HP with n-1 ascents).
    # But the only HP of a transitive tournament is the identity permutation (0,1,...,n-1).
    # Wait, is that the only one? In a transitive tournament on {0,...,n-1} with A[i][j]=1 for i<j,
    # a HP requires A[P_i][P_{i+1}] = 1, i.e., P_i < P_{i+1} for all i.
    # So the only HP is (0, 1, 2, ..., n-1), which has n-1 ascents.
    # Hence F = [0, ..., 0, 1]. Correct!

# 2. "Doubly regular" tournaments (exist for n = 4k+3)
print("\n--- Near-regular tournaments (n=5) ---")
for A in all_tournaments(5):
    scores = sorted([sum(A[i]) for i in range(5)])
    if scores == [2, 2, 2, 2, 2]:
        F = forward_edge_poly(A, 5)
        H = sum(F)
        print(f"  scores=(2,2,2,2,2): F={F}, H={H}")
        break

# 3. Tournaments with maximum H
print("\n--- Maximum H tournaments ---")
for n in [4, 5]:
    max_H = 0
    max_F = None
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        if H > max_H:
            max_H = H
            max_F = F
    print(f"  n={n}: max H={max_H}, F={max_F}")

# ===== SECTION 7: Can we express F(T,x) via a transfer matrix? =====
print("\n\n" + "=" * 70)
print("F(T,x) VIA TRANSFER MATRIX")
print("=" * 70)
print("""
Define the x-weighted adjacency matrix M_x:
  M_x[v,u] = A[v][u] * x^{[v < u]}

Then F(T,x) = sum of all products along HPs = "permanent-like" computation.

More precisely, F(T,x) can be computed via the DP:
  dp[mask, v] = polynomial in x
  dp[{v}, v] = 1
  dp[mask|{u}, u] += dp[mask, v] * M_x[v,u]

F(T,x) = sum_v dp[all, v]
""")

for n in [4]:
    print(f"\nn={n}:")
    for idx, A in enumerate(list(all_tournaments(n))[:3]):
        F = forward_edge_poly(A, n)
        H = sum(F)

        # Build M_x matrix (symbolic, but evaluate at x=1,2,3 to verify)
        for x_val in [1, 2]:
            Mx = [[A[v][u] * x_val**(1 if v < u else 0) for u in range(n)] for v in range(n)]

            # DP to compute sum of HP weights
            dp = {}
            for v in range(n):
                dp[(1 << v, v)] = 1

            for mask in range(1, 1 << n):
                for v in range(n):
                    if not (mask & (1 << v)) or (mask, v) not in dp: continue
                    for u in range(n):
                        if mask & (1 << u): continue
                        if Mx[v][u] > 0:
                            key = (mask | (1 << u), u)
                            dp[key] = dp.get(key, 0) + dp[(mask, v)] * Mx[v][u]

            full = (1 << n) - 1
            Fx = sum(dp.get((full, v), 0) for v in range(n))
            Fx_poly = sum(F[k] * x_val**k for k in range(n))
            print(f"  T#{idx}: F={F}, F({x_val})={Fx_poly}, DP={Fx}, match={Fx==Fx_poly}")

print("\n\nDone.")
