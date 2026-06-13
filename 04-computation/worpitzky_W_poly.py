#!/usr/bin/env python3
"""
W-POLYNOMIAL IN WORPITZKY BASIS — Deep investigation

Key finding: W(T,r) expanded in Worpitzky basis C(r+k, n-1) gives
PALINDROMIC coefficients at odd n.

This could be a new structural result connecting:
- Palindromicity of F(T,x) (known: F_k = F_{n-1-k})
- Even-r property of W at odd n (known: THM-L)
- Worpitzky / Eulerian number structure

Also investigate:
1. The exact relationship between palindromicity of F and palindromicity
   of W's Worpitzky expansion
2. Whether the W Worpitzky coefficients are always non-negative
3. Connection to the h-vector of some polytope
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

def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
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

def W_poly_coeffs(A, n):
    """Compute W(T,r) as polynomial in r."""
    W = np.zeros(n)
    for perm in permutations(range(n)):
        if not all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            continue
        poly = np.array([1.0])
        for i in range(n-1):
            s_e = A[perm[i]][perm[i+1]] - 0.5
            factor = np.array([s_e, 1.0])
            poly = np.convolve(poly, factor)
        for k in range(len(poly)):
            if k < n:
                W[k] += poly[k]
    return W

def worpitzky_basis_matrix(d):
    n = d + 1
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i][j] = comb(i + j, d)
    return M

def expand_in_worpitzky_rational(poly_coeffs):
    """Expand polynomial in Worpitzky basis, return exact rationals if possible."""
    d = len(poly_coeffs) - 1
    n = d + 1
    vals = [sum(poly_coeffs[k] * x**k for k in range(n)) for x in range(n)]
    M = worpitzky_basis_matrix(d)
    w = np.linalg.solve(M, vals)
    return w

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

print("=" * 70)
print("W(T,r) IN WORPITZKY BASIS — PALINDROMICITY TEST")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n{'='*60}")
    print(f"n = {n}, degree = {n-1}")
    print(f"{'='*60}")

    palindromic_count = 0
    nonneg_count = 0
    total = 0

    for A in all_tournaments(n):
        total += 1
        W = W_poly_coeffs(A, n)
        w = expand_in_worpitzky_rational(W)

        is_pal = all(abs(w[i] - w[n-1-i]) < 1e-8 for i in range(n))
        is_nonneg = all(x >= -1e-8 for x in w)

        if is_pal:
            palindromic_count += 1
        if is_nonneg:
            nonneg_count += 1

        if total <= 5 or (total <= 20 and not is_pal):
            H = ham_path_count_dp(A, n)
            print(f"  T#{total}: H={H}, W_worpitzky = {[round(x,6) for x in w]}")
            print(f"    palindromic={is_pal}, nonneg={is_nonneg}")

    print(f"\n  SUMMARY: palindromic={palindromic_count}/{total}, nonneg={nonneg_count}/{total}")

# At n=6,7 do sampling
for n in [6, 7]:
    print(f"\n{'='*60}")
    print(f"n = {n} (SAMPLED)")
    print(f"{'='*60}")

    palindromic_count = 0
    nonneg_count = 0
    total = 200

    for trial in range(total):
        A = random_tournament(n)
        W = W_poly_coeffs(A, n)
        w = expand_in_worpitzky_rational(W)

        is_pal = all(abs(w[i] - w[n-1-i]) < 1e-6 for i in range(n))
        is_nonneg = all(x >= -1e-6 for x in w)

        if is_pal:
            palindromic_count += 1
        if is_nonneg:
            nonneg_count += 1

        if trial < 5:
            H = ham_path_count_dp(A, n)
            print(f"  Trial {trial}: H={H}, W_worpitzky = {[round(x,4) for x in w]}")
            print(f"    palindromic={is_pal}, nonneg={is_nonneg}")

    print(f"\n  SUMMARY: palindromic={palindromic_count}/{total}, nonneg={nonneg_count}/{total}")

# ===== Now the KEY question: denominator structure =====
print("\n\n" + "=" * 70)
print("DENOMINATOR STRUCTURE OF W WORPITZKY COEFFICIENTS")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\nn={n}:")
    # Multiply by 2^{n-1} to get integer coefficients
    scale = 2**(n-1)

    for trial in range(10):
        A = random_tournament(n)
        W = W_poly_coeffs(A, n)
        w = expand_in_worpitzky_rational(W)
        w_scaled = [x * scale for x in w]

        # Check if scaled values are integers
        all_int = all(abs(x - round(x)) < 1e-6 for x in w_scaled)

        if trial < 3:
            print(f"  Trial {trial}: w = {[round(x,6) for x in w]}")
            print(f"    2^{n-1} * w = {[round(x,2) for x in w_scaled]}, all_int={all_int}")

# ===== CONNECTION: W Worpitzky vs F Worpitzky =====
print("\n\n" + "=" * 70)
print("RELATIONSHIP: W WORPITZKY <-> F WORPITZKY")
print("=" * 70)
print("W(T,r) = sum_P prod(r + s_e)")
print("F(T,x) = sum_P x^{fwd(P)}")
print("Relationship: W(T,r) = sum_k F_k (r+1/2)^k (r-1/2)^{n-1-k}")
print("")

for n in [3, 4, 5]:
    print(f"\nn={n}:")
    from math import comb as C

    for trial in range(3):
        A = random_tournament(n)

        # Compute F_k
        F = [0] * n
        for P in permutations(range(n)):
            if not all(A[P[i]][P[i+1]] for i in range(n-1)):
                continue
            fwd = sum(1 for i in range(n-1) if P[i] < P[i+1])
            F[fwd] += 1

        # Compute F Worpitzky expansion
        d = n - 1
        F_vals = [sum(F[k] * x**k for k in range(n)) for x in range(n)]
        M = worpitzky_basis_matrix(d)
        wF = np.linalg.solve(M, F_vals)

        # Compute W Worpitzky expansion
        W = W_poly_coeffs(A, n)
        wW = expand_in_worpitzky_rational(W)

        H = sum(F)
        print(f"  Trial {trial}: H={H}")
        print(f"    F_worpitzky  = {[int(round(x)) for x in wF]}")
        print(f"    W_worpitzky  = {[round(x,6) for x in wW]}")

        # The relationship between wF and wW should involve the Mobius transform
        # W(r) = (r-1/2)^{n-1} F((2r+1)/(2r-1))
        # In Worpitzky basis, this becomes...
        # Let's check: wW_k = wF_k / 2^{n-1} * something?
        ratios = []
        for k in range(n):
            if abs(wF[k]) > 0.01:
                ratios.append(wW[k] / wF[k])
            else:
                ratios.append(None)
        print(f"    wW/wF ratios = {[f'{r:.6f}' if r is not None else 'N/A' for r in ratios]}")

# ===== MULTIPLICATIVE STRUCTURE =====
print("\n\n" + "=" * 70)
print("WORPITZKY COEFFICIENTS AND CYCLE COUNTS")
print("=" * 70)

for n in [5]:
    print(f"\nn={n}: exhaustive")

    data = []
    for A in all_tournaments(n):
        F = [0] * n
        for P in permutations(range(n)):
            if not all(A[P[i]][P[i+1]] for i in range(n-1)):
                continue
            fwd = sum(1 for i in range(n-1) if P[i] < P[i+1])
            F[fwd] += 1

        H = sum(F)
        t3 = sum(1 for i in range(n) for j in range(n) for k in range(n)
                 if i<j<k and A[i][j] and A[j][k] and A[k][i]) + \
             sum(1 for i in range(n) for j in range(n) for k in range(n)
                 if i<j<k and A[j][i] and A[k][j] and A[i][k])

        d = n - 1
        F_vals = [sum(F[k] * x**k for k in range(n)) for x in range(n)]
        M_mat = worpitzky_basis_matrix(d)
        wF = np.linalg.solve(M_mat, F_vals)
        wF_int = tuple(int(round(x)) for x in wF)

        data.append((H, t3, wF_int, F))

    # Group by (H, t3) and see if Worpitzky coefficients are determined
    from collections import defaultdict
    groups = defaultdict(list)
    for H, t3, wF, F in data:
        groups[(H, t3)].append((wF, F))

    print("\n  By (H, t3):")
    for (H, t3), entries in sorted(groups.items()):
        unique_w = set(e[0] for e in entries)
        if len(unique_w) > 1:
            print(f"  H={H:3d}, t3={t3}: {len(unique_w)} distinct w vectors")
            for w in sorted(unique_w):
                count = sum(1 for e in entries if e[0] == w)
                print(f"    w={list(w)}, count={count}")

    # What determines w?
    print("\n  What invariants determine w?")
    print("  Testing if w is determined by (H, t3, score_sequence)...")

    groups2 = defaultdict(list)
    for A_idx, A in enumerate(all_tournaments(n)):
        if A_idx >= len(data): break
        H, t3, wF, F = data[A_idx]
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        groups2[(H, t3, scores)].append(wF)

    ambiguous = 0
    for key, ws in groups2.items():
        unique = set(ws)
        if len(unique) > 1:
            ambiguous += 1
    print(f"  Ambiguous (H, t3, scores): {ambiguous}/{len(groups2)}")
