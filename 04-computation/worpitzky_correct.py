#!/usr/bin/env python3
"""
CORRECTED Worpitzky expansion of W(T,r) using the CORRECT definition:
  W(T,r) = sum_P prod_{i} (r + s_e(P_i, P_{i+1}))
  where s_e = +1/2 if P_i < P_{i+1} (ascent), -1/2 if P_i > P_{i+1} (descent)

The previous code used s_e = A[P_i][P_{i+1}] - 1/2 which is always +1/2
since all edges in HPs exist (A[P_i][P_{i+1}]=1). That was trivially
W = H*(r+1/2)^{n-1}, explaining the "universal" type B Eulerian numbers.

Now let's see what the REAL W-polynomial looks like in Worpitzky basis.

Connection to F(T,x): W(T,r) = sum_k F_k (r+1/2)^k (r-1/2)^{n-1-k}
where k = number of ascents = forward edges.
"""
from itertools import permutations
from math import comb, factorial
import random
import numpy as np

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def W_poly_correct(A, n):
    """W(T,r) using ascent/descent s_e. Return polynomial coefficients in r."""
    W = np.zeros(n)
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        poly = np.array([1.0])
        for i in range(n-1):
            s_e = 0.5 if P[i] < P[i+1] else -0.5
            factor = np.array([s_e, 1.0])
            poly = np.convolve(poly, factor)
        for k in range(min(n, len(poly))):
            W[k] += poly[k]
    return W

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

def worpitzky_basis_matrix(d):
    n = d + 1
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i][j] = comb(i + j, d)
    return M

def expand_in_worpitzky(poly_coeffs):
    d = len(poly_coeffs) - 1
    n = d + 1
    vals = [sum(poly_coeffs[k] * x**k for k in range(n)) for x in range(n)]
    M = worpitzky_basis_matrix(d)
    w = np.linalg.solve(M, vals)
    return w

random.seed(42)

print("=" * 70)
print("CORRECT W(T,r) IN WORPITZKY BASIS")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n{'='*50}")
    print(f"EXHAUSTIVE n={n}")
    print(f"{'='*50}")

    palindromic_count = 0
    nonneg_count = 0
    integer_count = 0
    total = 0

    data = []
    for A in all_tournaments(n):
        total += 1
        W = W_poly_correct(A, n)
        w = expand_in_worpitzky(W)
        H = ham_path_count_dp(A, n)

        is_pal = all(abs(w[i] - w[n-1-i]) < 1e-8 for i in range(n))
        is_nonneg = all(x >= -1e-8 for x in w)
        is_int = all(abs(x - round(x)) < 1e-8 for x in w)

        if is_pal: palindromic_count += 1
        if is_nonneg: nonneg_count += 1
        if is_int: integer_count += 1

        data.append((H, w, W))

        if total <= 5:
            print(f"  T#{total}: H={H}")
            print(f"    W(r) coeffs: {[round(x,6) for x in W]}")
            print(f"    Worpitzky: {[round(x,6) for x in w]}")
            print(f"    palindromic={is_pal}, nonneg={is_nonneg}, integer={is_int}")

    print(f"\n  SUMMARY: palindromic={palindromic_count}/{total}, nonneg={nonneg_count}/{total}, integer={integer_count}/{total}")

    # Check if proportional to H
    scale = 2**(n-1)
    print(f"\n  Checking if w_k = H * universal / 2^{n-1}:")
    all_proportional = True
    for H, w, W in data:
        ratios = [w[k] * scale / H if H != 0 else None for k in range(n)]
        if total <= 10 or not all_proportional:
            pass  # print for debugging
        # Check constancy
        break  # just check the structure

    # Group by H and see if w/H is constant
    from collections import defaultdict
    by_H = defaultdict(list)
    for H, w, W in data:
        by_H[H].append([round(x, 8) for x in w])

    for H in sorted(by_H.keys()):
        ws = by_H[H]
        # Normalize by H
        normalized = [tuple(round(x/H, 8) for x in w) for w in ws]
        unique = set(normalized)
        if len(unique) == 1:
            print(f"  H={H}: w/H = {list(unique)[0]} (CONSTANT)")
        else:
            print(f"  H={H}: {len(unique)} distinct w/H values")
            for u in sorted(unique)[:3]:
                count = sum(1 for x in normalized if x == u)
                print(f"    w/H = {list(u)}, count={count}")

# ===== Now with n=6,7 sampling =====
for n in [6, 7]:
    print(f"\n{'='*50}")
    print(f"SAMPLED n={n}")
    print(f"{'='*50}")

    palindromic_count = 0
    nonneg_count = 0
    total = 100

    by_H = defaultdict(list)
    for trial in range(total):
        A = random_tournament(n)
        W = W_poly_correct(A, n)
        w = expand_in_worpitzky(W)
        H = ham_path_count_dp(A, n)

        is_pal = all(abs(w[i] - w[n-1-i]) < 1e-6 for i in range(n))
        is_nonneg = all(x >= -1e-6 for x in w)

        if is_pal: palindromic_count += 1
        if is_nonneg: nonneg_count += 1

        by_H[H].append([round(x, 6) for x in w])

        if trial < 3:
            print(f"  Trial {trial}: H={H}")
            print(f"    W(r) coeffs: {[round(x,4) for x in W]}")
            print(f"    Worpitzky: {[round(x,4) for x in w]}")
            print(f"    palindromic={is_pal}, nonneg={is_nonneg}")

    print(f"\n  SUMMARY: palindromic={palindromic_count}/{total}, nonneg={nonneg_count}/{total}")

    # Check proportionality to H
    print(f"\n  Proportionality check:")
    for H in sorted(by_H.keys()):
        ws = by_H[H]
        if len(ws) < 2: continue
        normalized = [tuple(round(x/H, 6) for x in w) for w in ws]
        unique = set(normalized)
        if len(unique) == 1:
            print(f"  H={H} ({len(ws)} samples): w/H CONSTANT = {list(unique)[0][:3]}...")
        else:
            print(f"  H={H} ({len(ws)} samples): {len(unique)} distinct w/H values")
            break

# ===== Check at r=1/2: should give H since W(T,1/2) = H/2^{n-1} =====
print("\n\n" + "=" * 70)
print("W(T,r) PROPERTIES AT ODD n")
print("=" * 70)

for n in [3, 5, 7]:
    print(f"\nn={n}:")
    for trial in range(3):
        A = random_tournament(n)
        W = W_poly_correct(A, n)
        H = ham_path_count_dp(A, n)

        # Check only even powers
        even_only = all(abs(W[k]) < 1e-8 for k in range(n) if k % 2 == 1)
        print(f"  Trial {trial}: H={H}, W_coeffs = {[round(x,6) for x in W]}")
        print(f"    Even powers only? {even_only}")

        # Evaluate at r = 1/2
        W_at_half = sum(W[k] * 0.5**k for k in range(n))
        print(f"    W(1/2) = {W_at_half:.6f}, expected = 0 (since F is palindromic and eval at x=1 gives H but W(1/2) = ...)")

        # Actually W(T, 1/2) = sum_P prod(1/2 + s_e) where s_e = ±1/2
        # So each factor is either 1 (if ascent) or 0 (if descent)
        # So W(T, 1/2) = #{HPs with no descents} = #{monotone increasing HPs} = 0 or 1

        # W(T, -1/2) = sum_P prod(-1/2 + s_e) = sum_P prod(0 or -1)
        # = sum_P (-1)^{#descents} * 0^{#ascents} = #{HPs with no ascents} * (-1)^{n-1}

        # H = F(T,1) = sum F_k. But W(T,r) = sum_k F_k (r+1/2)^k (r-1/2)^{n-1-k}
        # At r=1/2: W = sum F_k * 1^k * 0^{n-1-k} = F_{n-1} (HPs that are fully ascending)
        print(f"    W(1/2) should be F_{{n-1}} = #{'{'}HPs where all P_i < P_{{i+1}}{'}'}.")
