#!/usr/bin/env python3
"""
Paley H/E[H] Ratio — Seeking a Closed-Form Formula

Observation: H(P_p) / E[H_random] appears to approach a constant ~2.4.

E[H_random] = p! / 2^(p-1) for tournaments on p vertices.

If H(P_p) / E[H] → c, then H(P_p) ~ c · p! / 2^(p-1).

This script investigates:
1. The exact ratio H(P_p) / (p!/2^{p-1}) for small p
2. Whether the ratio has a closed form (e.g., involving p, sqrt(p), etc.)
3. Connection to the Paley determinant formula det(I+A) = (p+1)^{(p+1)/2}/2^p
4. Whether H(P_p) can be expressed in terms of Gauss sums
5. Asymptotic behavior as p → ∞

Author: opus-2026-03-13-S67k
"""

import numpy as np
from math import factorial, log, log2, sqrt, comb
from collections import defaultdict

def paley_adjacency(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A

def count_hp(A):
    n = A.shape[0]
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for j in range(n):
            if not (mask & (1 << j)):
                continue
            prev_mask = mask ^ (1 << j)
            if prev_mask == 0:
                continue
            for i in range(n):
                if (prev_mask & (1 << i)) and A[i][j]:
                    dp[mask][j] += dp[prev_mask][i]
    full = (1 << n) - 1
    return sum(dp[full][j] for j in range(n))

print("=" * 70)
print("PALEY TOURNAMENT H / E[H] RATIO ANALYSIS")
print("=" * 70)

# Compute H for Paley tournaments
# p ≡ 3 mod 4: p = 3, 7, 11
# p = 3: n=3, feasible
# p = 7: n=7, feasible (7! = 5040)
# p = 11: n=11, HP computation takes O(2^11 * 11) = O(22528) steps, feasible
# p = 19: too expensive for HP (2^19 * 19 ≈ 10M) — might work with patience

print("\n--- Exact H(P_p) for small primes ---\n")

results = []
for p in [3, 7, 11]:
    A = paley_adjacency(p)
    H = count_hp(A)
    E_H = factorial(p) / (2 ** (p - 1))
    ratio = H / E_H
    m = (p - 1) // 2

    # Also compute det(I+A)
    det_IpA = np.linalg.det(np.eye(p) + A.astype(float))
    det_formula = (p + 1) ** ((p + 1) / 2) / (2 ** p)

    results.append({
        'p': p, 'H': H, 'E_H': E_H, 'ratio': ratio,
        'det': det_IpA, 'det_formula': det_formula, 'm': m
    })

    print(f"p = {p}:")
    print(f"  H(P_{p}) = {H}")
    print(f"  E[H_random] = p!/2^(p-1) = {E_H:.4f}")
    print(f"  H/E[H] = {ratio:.6f}")
    print(f"  det(I+A) = {det_IpA:.2f} (formula: {det_formula:.2f})")
    print(f"  H/det(I+A) = {H/det_IpA:.6f}")
    print(f"  sqrt(p) = {sqrt(p):.6f}")
    print(f"  H/(p-1)! = {H/factorial(p-1):.6f}")
    print(f"  H·2^m/p! = {H * 2**m / factorial(p):.6f}")
    print()

# Check if H(P_p) has a nice factorization
print("--- Factorization of H(P_p) ---\n")
for r in results:
    p, H = r['p'], r['H']
    # Factor H
    h = H
    factors = []
    for f in range(2, h + 1):
        while h % f == 0:
            factors.append(f)
            h //= f
        if h == 1:
            break
    print(f"H(P_{p}) = {H} = {'·'.join(map(str, factors))}")

# Check specific formulas
print("\n--- Testing specific closed-form candidates ---\n")
for r in results:
    p, H = r['p'], r['H']

    # Candidate 1: H = p!! * something (double factorial)
    # Candidate 2: H = (p-1)!! * p
    # Candidate 3: H involves Catalan-like numbers
    # Candidate 4: H = (2m+1)!/(2^m) * correction

    # p=3: H=3. 3! = 6, 3!/2 = 3. So H = p!/2^1 = p!/2^m where m=1.
    # p=7: H=189. 7!/2^3 = 5040/8 = 630. H/630 = 0.3.
    #   189 = 7·27 = 7·3^3. Also 189 = 7·27 = (p)·(p-1)^m... 7·3^3... 3^3 = (p-1)/2)^m = 3^3
    # p=11: H=95095. 95095 = 5·7·11·13·19 ... let me compute

    # Check: H = p · ((p-1)/2)^((p-1)/2)?
    m = (p - 1) // 2
    candidate = p * (m ** m)
    print(f"p={p}: H={H}, p·m^m = {candidate}, ratio = {H/candidate:.4f}")

print("\n--- More formula candidates ---\n")
for r in results:
    p, H = r['p'], r['H']
    m = (p - 1) // 2

    # Check: H = C(p, m) * something?
    binom = comb(p, m)
    print(f"p={p}: H={H}, C(p,m)={binom}, H/C(p,m)={H/binom:.4f}")

    # Check: H = (2m)! / m! * something?
    ratio1 = H * factorial(m) / factorial(2 * m)
    print(f"  H·m!/(2m)! = {ratio1:.6f}")

    # Check: H = m! · (m+1)?
    candidate2 = factorial(m) * (m + 1)
    print(f"  m!·(m+1) = {candidate2}, H/{candidate2} = {H/candidate2:.4f}")

    # p=3: 3 / (1!·2) = 1.5
    # p=7: 189 / (3!·4) = 189/24 = 7.875
    # p=11: 95095 / (5!·6) = 95095/720 = 132.076

    # Check H against permanent of A
    # permanent = sum over permutations of product of A[i,sigma(i)]
    # For small p, this is feasible
    from itertools import permutations
    if p <= 11:
        perm_sum = 0
        for sigma in permutations(range(p)):
            prod = 1
            for i in range(p):
                prod *= A[i][sigma[i]]
            perm_sum += prod
        print(f"  perm(A) = {perm_sum}")
        print(f"  H/perm(A) = {H/perm_sum:.6f}" if perm_sum > 0 else "  perm(A) = 0")

    print()

# Look at the H/E[H] ratio growth
print("--- H/E[H] ratio pattern ---\n")
for r in results:
    p = r['p']
    ratio = r['ratio']
    m = (p - 1) // 2
    # Is ratio ~ e^{something}?
    print(f"p={p}: ratio={ratio:.6f}, ln(ratio)={log(ratio):.6f}, ratio/sqrt(p)={ratio/sqrt(p):.6f}")
    print(f"  ratio/p^(1/4) = {ratio/p**0.25:.6f}")
    print(f"  ratio/ln(p) = {ratio/log(p):.6f}")

# Try to compute H for P_19 (2^19 · 19 ≈ 10M, might take a few minutes)
# Actually Held-Karp is O(2^n · n²) = O(2^19 · 361) ≈ 190M ops. Might work.
print("\n--- Attempting P_19 H computation ---")

try:
    import sys
    p = 19
    A = paley_adjacency(p)
    print(f"Starting Held-Karp for P_{p} (2^{p} = {2**p})...")
    sys.stdout.flush()

    # Optimized Held-Karp using numpy
    n = p
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for i in range(n):
        dp[1 << i][i] = 1

    for mask in range(1, 1 << n):
        for j in range(n):
            if not (mask & (1 << j)):
                continue
            prev_mask = mask ^ (1 << j)
            if prev_mask == 0:
                continue
            val = 0
            for i in range(n):
                if (prev_mask & (1 << i)) and A[i][j]:
                    val += dp[prev_mask][i]
            dp[mask][j] = val

    full = (1 << n) - 1
    H_19 = int(dp[full].sum())
    E_H_19 = factorial(p) / (2 ** (p - 1))
    ratio_19 = H_19 / E_H_19
    print(f"H(P_19) = {H_19}")
    print(f"E[H] = {E_H_19:.2f}")
    print(f"H/E[H] = {ratio_19:.6f}")

    # Factorize
    h = H_19
    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    for f in small_primes:
        count = 0
        while h % f == 0:
            count += 1
            h //= f
        if count > 0:
            print(f"  {f}^{count}", end="")
    if h > 1:
        print(f" · {h}", end="")
    print()

except MemoryError:
    print("MemoryError: P_19 Held-Karp needs ~5GB, skipping")
except Exception as e:
    print(f"Error: {e}")

print("\nDone.")
