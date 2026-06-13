#!/usr/bin/env python3
"""
STRUCTURE OF W(T,r) WORPITZKY COEFFICIENTS

CONFIRMED: W(T,r) = sum_k w_k C(r+k, n-1) where:
1. w_k >= 0 always (non-negative)
2. w_k = w_{n-1-k} (palindromic)
3. 2^{n-1} * w_k is always a positive integer

QUESTIONS:
A. What is the formula for 2^{n-1} * w_k?
B. Is 2^{n-1} * w_k = H * something_universal?
C. What is the relationship between w_k and cycle counts?
D. Can we prove palindromicity + non-negativity algebraically?
E. Does this give a new proof of OCF?
"""
from itertools import permutations
from math import comb, factorial
import random
import numpy as np
from collections import defaultdict

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def W_poly_coeffs(A, n):
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

def expand_in_worpitzky(poly_coeffs):
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

# ===== QUESTION A: Is 2^{n-1} * w_k = H * universal_k? =====
print("=" * 70)
print("QUESTION A: Is 2^{n-1} * w_k proportional to H?")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\nn={n}:")
    scale = 2**(n-1)

    ratios_all = defaultdict(list)
    for A in all_tournaments(n):
        W = W_poly_coeffs(A, n)
        w = expand_in_worpitzky(W)
        H = ham_path_count_dp(A, n)

        scaled = [int(round(x * scale)) for x in w]
        for k in range(n):
            if H > 0:
                ratios_all[k].append(scaled[k] / H)

    for k in range(n):
        vals = set(round(r, 8) for r in ratios_all[k])
        if len(vals) == 1:
            print(f"  k={k}: 2^{n-1}*w_k/H = {vals.pop()} (CONSTANT!)")
        else:
            print(f"  k={k}: 2^{n-1}*w_k/H takes {len(vals)} values: {sorted(vals)[:5]}...")

# ===== QUESTION B: Exact formula for scaled coefficients =====
print("\n\n" + "=" * 70)
print("QUESTION B: What is 2^{n-1} * w_k exactly?")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\nn={n}:")
    scale = 2**(n-1)

    # Collect all (H, scaled_w) pairs
    hw_pairs = defaultdict(int)
    for A in all_tournaments(n):
        W = W_poly_coeffs(A, n)
        w = expand_in_worpitzky(W)
        H = ham_path_count_dp(A, n)
        scaled = tuple(int(round(x * scale)) for x in w)
        hw_pairs[(H, scaled)] += 1

    for (H, scaled), count in sorted(hw_pairs.items()):
        print(f"  H={H:3d}: 2^{n-1}*w = {list(scaled)}, count={count}")

# ===== QUESTION C: Eulerian number connection =====
print("\n\n" + "=" * 70)
print("QUESTION C: Connection to Eulerian numbers")
print("=" * 70)

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

for n in [3, 4, 5]:
    d = n - 1
    E = [eulerian_number(d, k) for k in range(d+1)]
    print(f"\nn={n}: Eulerian numbers A({d},k) = {E}")
    print(f"  sum = {sum(E)} = {factorial(d)}")

    scale = 2**(n-1)
    # For transitive tournament: H=1, F = [0,...,0,1]
    # W(transitive, r) = (r+1/2)^{n-1} since all edges are forward
    # Worpitzky of (r+1/2)^{n-1} = sum_k A(n-1,k) C(r+1/2+k, n-1)... hmm
    # Actually (r+1/2)^{n-1} is not the same as r^{n-1}.
    # Let me compute it directly.
    A_trans = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            A_trans[i][j] = 1
    W_trans = W_poly_coeffs(A_trans, n)
    w_trans = expand_in_worpitzky(W_trans)
    scaled_trans = [int(round(x * scale)) for x in w_trans]
    print(f"  Transitive: 2^{d}*w = {scaled_trans}")
    print(f"  Ratios to Eulerian: {[f'{s/e:.4f}' if e != 0 else 'N/A' for s, e in zip(scaled_trans, E + [0])]}")

# ===== QUESTION D: Can we PROVE non-negativity? =====
print("\n\n" + "=" * 70)
print("QUESTION D: Algebraic structure of W Worpitzky coefficients")
print("=" * 70)

# W(T,r) = sum_P prod_{edges e} (r + s_e)
# where s_e = ±1/2.
#
# For a single path P with f forward edges:
# prod_e (r + s_e) = (r+1/2)^f * (r-1/2)^{n-1-f}
#
# So W(T,r) = sum_f F_f * (r+1/2)^f * (r-1/2)^{n-1-f}
#
# Now (r+1/2)^f * (r-1/2)^{n-1-f} in Worpitzky basis...
# Let's compute this.

print("\n(r+1/2)^f * (r-1/2)^{n-1-f} in Worpitzky basis:")
for n in [4, 5]:
    d = n - 1
    print(f"\nn={n} (degree {d}):")
    scale = 2**d
    for f in range(n):
        # Polynomial (r+1/2)^f * (r-1/2)^{n-1-f}
        # Expand as polynomial in r
        poly = np.array([1.0])
        for i in range(f):
            poly = np.convolve(poly, [0.5, 1.0])  # (r+1/2)
        for i in range(d - f):
            poly = np.convolve(poly, [-0.5, 1.0])  # (r-1/2)
        w = expand_in_worpitzky(poly)
        scaled = [round(x * scale) for x in w]
        all_nonneg = all(x >= -0.001 for x in w)
        print(f"  f={f}: 2^{d}*w = {scaled}, nonneg={all_nonneg}")

# So the key question is: is each (r+1/2)^f (r-1/2)^{d-f} non-negative in Worpitzky basis?
# If YES, then W = sum F_f * [non-negative] is automatically non-negative.

print("\n\n" + "=" * 70)
print("KEY TEST: Is (r+1/2)^f (r-1/2)^{d-f} non-neg in Worpitzky?")
print("=" * 70)

for d in range(2, 10):
    all_nonneg = True
    for f in range(d+1):
        poly = np.array([1.0])
        for i in range(f):
            poly = np.convolve(poly, [0.5, 1.0])
        for i in range(d - f):
            poly = np.convolve(poly, [-0.5, 1.0])
        w = expand_in_worpitzky(poly)
        if any(x < -1e-10 for x in w):
            all_nonneg = False
            break
    print(f"  d={d}: all (r+1/2)^f(r-1/2)^{{d-f}} non-neg in Worpitzky? {all_nonneg}")

# ===== QUESTION E: What IS C(r+k, d) combinatorially? =====
print("\n\n" + "=" * 70)
print("C(r+k, d) AT HALF-INTEGER POINTS")
print("=" * 70)
# At r = 1/2: W(T, 1/2) = H/2^{n-1}
# So sum_k w_k C(1/2 + k, d) = H/2^{n-1}
# C(1/2 + k, d) = product_{j=0}^{d-1} (1/2 + k - j) / d!
#               = product_{j=0}^{d-1} (2k+1-2j) / (2^d * d!)

for d in range(2, 7):
    print(f"\nd={d}:")
    for k in range(d+1):
        val = 1.0
        for j in range(d):
            val *= (0.5 + k - j)
        val /= factorial(d)
        print(f"  C(1/2+{k}, {d}) = {val:.6f} = {val * 2**d * factorial(d):.0f} / (2^{d} * {d}!)")

# ===== VERIFY: sum_k w_k C(1/2+k, d) = H / 2^d =====
print("\n\n" + "=" * 70)
print("VERIFY: sum_k w_k C(1/2+k, d) = H / 2^d")
print("=" * 70)

for n in [5, 7]:
    d = n - 1
    print(f"\nn={n} (d={d}):")

    # Compute C(1/2+k, d) values
    C_vals = []
    for k in range(d+1):
        val = 1.0
        for j in range(d):
            val *= (0.5 + k - j)
        val /= factorial(d)
        C_vals.append(val)

    for trial in range(5):
        A = random_tournament(n)
        W = W_poly_coeffs(A, n)
        w = expand_in_worpitzky(W)
        H = ham_path_count_dp(A, n)

        check = sum(w[k] * C_vals[k] for k in range(d+1))
        expected = H / 2**d
        print(f"  Trial {trial}: H={H}, sum w_k C(1/2+k,{d}) = {check:.8f}, H/2^{d} = {expected:.8f}, match={abs(check-expected)<1e-8}")
