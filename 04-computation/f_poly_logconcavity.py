#!/usr/bin/env python3
"""
f_poly_logconcavity.py — Is F(T,x) log-concave? Ultra-log-concave?

F(T,x) = sum_k F_k x^k is palindromic: F_k = F_{n-1-k}.

LOG-CONCAVITY: F_k^2 >= F_{k-1} * F_{k+1} for all k.
This is a STRONG structural property. The Eulerian polynomial A_n(x) IS
log-concave (proved by Brenti). But F(T,x) for a fixed tournament may
or may not be.

ULTRA-LOG-CONCAVITY: F_k^2 / C(n-1,k)^2 >= F_{k-1}*F_{k+1} / (C(n-1,k-1)*C(n-1,k+1))
Equivalent to F_k/C(n-1,k) being log-concave.

REAL-ROOTEDNESS: If F(T,x) has all real roots, then it's log-concave.
From S44: F(T,x) zeros at n=5 H=9 are ALL on unit circle. But are they real?
No — palindromic polynomials of even degree have roots in reciprocal pairs,
which for unit circle means conjugate pairs. Not necessarily real.

GAMMA-POSITIVITY: A palindromic polynomial of degree d is gamma-positive if
F(x) = sum_{k=0}^{d/2} gamma_k x^k (1+x)^{d-2k} with gamma_k >= 0.
Gamma-positivity implies log-concavity and unimodality.

Author: opus-2026-03-07-S46
"""
from itertools import permutations, combinations
from math import comb, factorial

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def is_log_concave(F):
    n = len(F)
    for k in range(1, n-1):
        if F[k]**2 < F[k-1] * F[k+1]:
            return False
    return True

def gamma_coefficients(F):
    """Compute gamma-expansion: F(x) = sum gamma_k x^k (1+x)^{d-2k}."""
    d = len(F) - 1  # degree
    # F(x) = sum_{k=0}^{d} F_k x^k
    # gamma expansion: F(x) = sum_{j=0}^{d//2} gamma_j x^j (1+x)^{d-2j}
    # Expanding (1+x)^{d-2j} = sum_i C(d-2j,i) x^i
    # So F_k = sum_{j=0}^{min(k,d//2)} gamma_j * C(d-2j, k-j)

    gammas = [0] * (d//2 + 1)
    remaining = list(F)

    for j in range(d//2 + 1):
        # gamma_j is determined by remaining[j]
        gammas[j] = remaining[j]
        # Subtract gamma_j * x^j * (1+x)^{d-2j}
        for i in range(d - 2*j + 1):
            remaining[j + i] -= gammas[j] * comb(d - 2*j, i)

    # Check: remaining should be all zeros
    if any(abs(r) > 0.001 for r in remaining):
        return gammas, False  # decomposition didn't work cleanly
    return gammas, True

# ============================================================
# LOG-CONCAVITY TEST
# ============================================================
print("=" * 60)
print("LOG-CONCAVITY OF F(T,x)")
print("=" * 60)

for n in [4, 5, 6, 7]:
    m = n*(n-1)//2
    seen = set()
    lc_pass = 0
    lc_fail = 0

    import random
    random.seed(42)
    num = min(1 << m, 50000)

    for trial in range(num):
        if n <= 5:
            bits = trial
        else:
            bits = random.getrandbits(m)

        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        if is_log_concave(F):
            lc_pass += 1
        else:
            lc_fail += 1
            if lc_fail <= 3:
                print(f"  n={n} FAIL: F={F}")
                for k in range(1, n-1):
                    ratio = F[k]**2 / (F[k-1]*F[k+1]) if F[k-1]*F[k+1] > 0 else float('inf')
                    fail_mark = " ***" if F[k]**2 < F[k-1]*F[k+1] else ""
                    print(f"    k={k}: F_k^2={F[k]**2}, F_{k-1}*F_{k+1}={F[k-1]*F[k+1]}, ratio={ratio:.4f}{fail_mark}")

    print(f"\nn={n}: {len(seen)} distinct, log-concave={lc_pass}, FAILS={lc_fail}")

# ============================================================
# GAMMA-POSITIVITY TEST
# ============================================================
print("\n" + "=" * 60)
print("GAMMA-POSITIVITY OF F(T,x)")
print("=" * 60)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    seen = set()
    gp_pass = 0
    gp_fail = 0

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        gammas, valid = gamma_coefficients(F)
        if valid and all(g >= 0 for g in gammas):
            gp_pass += 1
        else:
            gp_fail += 1
            if gp_fail <= 3:
                print(f"  n={n} FAIL: F={F}, gammas={gammas}")

    print(f"\nn={n}: {len(seen)} distinct, gamma-positive={gp_pass}, FAILS={gp_fail}")

# ============================================================
# NORMALIZED F: is F_k / C(n-1,k) log-concave?
# ============================================================
print("\n" + "=" * 60)
print("ULTRA-LOG-CONCAVITY: F_k/C(n-1,k)")
print("=" * 60)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    seen = set()
    ulc_pass = 0
    ulc_fail = 0

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        # Normalized: g_k = F_k / C(n-1,k)
        g = [F[k] / comb(n-1, k) for k in range(n)]
        ulc = True
        for k in range(1, n-1):
            if g[k]**2 < g[k-1] * g[k+1] - 1e-10:
                ulc = False
                break
        if ulc:
            ulc_pass += 1
        else:
            ulc_fail += 1
            if ulc_fail <= 3:
                print(f"  n={n} FAIL: F={F}, g={[f'{x:.2f}' for x in g]}")

    print(f"\nn={n}: {len(seen)} distinct, ultra-log-concave={ulc_pass}, FAILS={ulc_fail}")

# ============================================================
# F(T,x) / (n-1)! — is this a probability distribution?
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) AS DISTRIBUTION: p_k = F_k / n!")
print("=" * 60)

# p_k = F_k / n! is the probability that a random permutation of T has
# exactly k forward edges. This is a probability distribution on {0,...,n-1}.
# Variance, skewness, kurtosis?

for n in [5, 7]:
    m = n*(n-1)//2
    seen = set()

    import random
    random.seed(42)
    num = min(1 << m, 50000)

    variances = []
    for trial in range(num):
        if n <= 5:
            bits = trial
        else:
            bits = random.getrandbits(m)

        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        total = sum(F)
        mean = sum(k * F[k] for k in range(n)) / total
        var = sum((k - mean)**2 * F[k] for k in range(n)) / total
        variances.append(var)

    # By palindrome, mean = (n-1)/2 always
    print(f"\nn={n}: mean = (n-1)/2 = {(n-1)/2}")
    print(f"  variance range: [{min(variances):.4f}, {max(variances):.4f}]")
    print(f"  variance at transitive: {variances[0]:.4f}")
    # Variance of uniform on {0,...,n-1} would be (n^2-1)/12
    print(f"  max possible variance (Bernoulli at endpoints): {((n-1)/2)**2:.4f}")
    print(f"  variance of A_n(x): {(n-1)/12:.4f}")
    # Number of distinct variances
    rounded = set(round(v, 4) for v in variances)
    print(f"  distinct variances: {len(rounded)}")

# ============================================================
# ENTROPY of F(T,x) distribution
# ============================================================
print("\n" + "=" * 60)
print("ENTROPY H(F) = -sum p_k log p_k")
print("=" * 60)

import math

for n in [5]:
    m = n*(n-1)//2
    seen = set()

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        total = sum(F)
        entropy = -sum((F[k]/total) * math.log2(F[k]/total) for k in range(n) if F[k] > 0)

        if len(seen) <= 8:
            H = sum(F)
            t3 = sum(1 for triple in combinations(range(n), 3)
                     for perm in [(triple[0],triple[1],triple[2]), (triple[0],triple[2],triple[1])]
                     if adj[perm[0]][perm[1]] and adj[perm[1]][perm[2]] and adj[perm[2]][perm[0]])
            print(f"  F={F}, t3={t3}, entropy={entropy:.4f}")

    # Max entropy = log2(n) for uniform distribution
    print(f"  Max entropy (uniform) = log2({n}) = {math.log2(n):.4f}")
