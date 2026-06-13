#!/usr/bin/env python3
"""
chi_sign_pattern.py — opus-2026-03-13-S71

Investigate:
1. chi/n for Interval tournaments: is it (-1)^{(n-1)/2} for n≥7?
2. prod(Q_k, k≠0) = 1 for Interval: prove this algebraically
3. The alternating sum of Omega: alt_Omega = chi for all n
"""

import numpy as np
from math import comb

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def get_regular_paths(A, m):
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

# ============================================================
print("="*70)
print("PRODUCT IDENTITY: prod(Q_k, k≠0)")
print("="*70)

for n in [3, 5, 7, 9, 11, 13, 15]:
    m = (n-1)//2
    S = set(range(1, m+1))
    omega = np.exp(2j * np.pi / n)

    Q = []
    for k in range(n):
        S_hat = sum(omega**(k*s) for s in S)
        Q.append(abs(S_hat)**2)

    prod_Q = np.prod([Q[k] for k in range(1, n)])
    sum_Q = sum(Q[k] for k in range(1, n))

    print(f"  n={n:2d}: prod(Q_k,k≠0)={prod_Q:12.4f}, "
          f"sum(Q_k,k≠0)={(n*n-1)/4:.1f} (pred), {sum_Q:.4f} (actual)")

# ============================================================
print(f"\n{'='*70}")
print("CHI / n FOR INTERVAL TOURNAMENTS")
print("="*70)

# The Euler characteristic is the alternating sum of Omega/Betti
# For per-eigenspace: chi_per = sum(-1)^m * omega_m_per_eig
# Since all eigenspaces are uniform: chi = n * chi_per

for n in [3, 5, 7, 9]:
    m_half = (n-1)//2
    S = set(range(1, m_half+1))
    A = circulant_tournament(n, S)

    omegas = []
    for m in range(n):
        paths = get_regular_paths(A, m)
        omegas.append(len(paths))

    alt_omega = sum((-1)**m * omegas[m] for m in range(n))
    per_eig_omega = [o // n for o in omegas]
    per_eig_chi = sum((-1)**m * per_eig_omega[m] for m in range(n))

    predicted_sign = (-1)**m_half if n >= 7 else 0
    print(f"  n={n}: Omega={omegas}")
    print(f"        per_eig_Omega={per_eig_omega}")
    print(f"        alt_Omega={alt_omega}, chi/n={alt_omega//n}")
    print(f"        per_eig_chi={per_eig_chi}")
    print(f"        predicted (-1)^{m_half}={predicted_sign}")
    print()

# ============================================================
print("="*70)
print("ALGEBRAIC ANALYSIS: prod(Q_k) for Interval")
print("="*70)

# S = {1,2,...,m} for n=2m+1
# S_hat(k) = sum_{s=1}^{m} omega^{ks}
# = omega^k * (1 - omega^{mk}) / (1 - omega^k)
# = omega^k * (omega^{mk} - 1) / (omega^k - 1)  [careful with sign]

# |S_hat(k)|^2 = |sum_{s=1}^m omega^{ks}|^2
# = |omega^k * (1 - omega^{mk}) / (1 - omega^k)|^2
# = |(1 - omega^{mk}) / (1 - omega^k)|^2
# = sin^2(m*pi*k/n) / sin^2(pi*k/n)

# Using the Chebyshev identity: this is U_{m-1}^2(cos(pi*k/n))
# where U is a Chebyshev polynomial of the second kind.

# prod_{k=1}^{n-1} sin^2(m*pi*k/n) / sin^2(pi*k/n)

# Known: prod_{k=1}^{n-1} sin(pi*k/n) = n / 2^{n-1}
# And:   prod_{k=1}^{n-1} sin(m*pi*k/n) = ?

# Actually for the Interval tournament with S={1,...,m}, n=2m+1:
# Q_k = sin^2(m*pi*k/n) / sin^2(pi*k/n)

# Let's verify numerically:
for n in [5, 7, 9, 11]:
    m = (n-1)//2
    omega = np.exp(2j * np.pi / n)
    print(f"\n  n={n}, m={m}:")

    prod_sin_ratio = 1.0
    for k in range(1, n):
        sin_m = np.sin(m * np.pi * k / n)
        sin_1 = np.sin(np.pi * k / n)
        Q_k_trig = sin_m**2 / sin_1**2

        S_hat = sum(omega**(k*s) for s in range(1, m+1))
        Q_k_direct = abs(S_hat)**2

        prod_sin_ratio *= sin_m**2 / sin_1**2

    print(f"    prod Q_k = {prod_sin_ratio:.6f}")

    # Theoretical: prod_{k=1}^{n-1} sin(m*pi*k/n)
    prod_sin_m = np.prod([np.sin(m*np.pi*k/n) for k in range(1,n)])
    prod_sin_1 = np.prod([np.sin(np.pi*k/n) for k in range(1,n)])
    print(f"    prod sin(m*pi*k/n) = {prod_sin_m:.6f}")
    print(f"    prod sin(pi*k/n) = {prod_sin_1:.6f}")
    print(f"    n/2^(n-1) = {n/2**(n-1):.6f}")
    print(f"    ratio = {prod_sin_m/prod_sin_1:.6f}")
    print(f"    ratio^2 = {(prod_sin_m/prod_sin_1)**2:.6f}")

    # Known: prod_{k=1}^{n-1} sin(pi*k/n) = n / 2^{n-1}
    # For prod_{k=1}^{n-1} sin(m*pi*k/n):
    # Since gcd(m, n) = 1 (n=2m+1 odd, m < n):
    # {mk mod n : k=1,...,n-1} = {1,...,n-1} (permutation)
    # So prod_{k=1}^{n-1} sin(m*pi*k/n) = prod_{k=1}^{n-1} sin(pi*k/n) = n/2^{n-1}
    #
    # THEREFORE: prod Q_k = (n/2^{n-1})^2 / (n/2^{n-1})^2 = 1 !!! QED

    print(f"    => gcd(m,n)={np.gcd(m,n)}: m*k permutes residues mod n")
    print(f"    => prod sin(m*pi*k/n) = prod sin(pi*k/n)")
    print(f"    => prod Q_k = 1  QED")

# ============================================================
print(f"\n{'='*70}")
print("GENERALIZATION: prod(Q_k, k≠0) = 1 for ANY circulant with gcd(|S_generator|,n)=1?")
print("="*70)

# Test for Paley p=7
for p in [7, 11]:
    QR = {a % p for a in range(1, p) if pow(a, (p-1)//2, p) == 1}
    omega = np.exp(2j * np.pi / p)

    prod_Q = 1.0
    for k in range(1, p):
        S_hat = sum(omega**(k*s) for s in QR)
        prod_Q *= abs(S_hat)**2

    print(f"  Paley p={p}, QR={sorted(QR)}: prod(Q_k,k≠0) = {prod_Q:.4f}")

    # For Paley: Q_k = (p+1)/4 for all k≠0
    # So prod = ((p+1)/4)^(p-1)
    q = (p+1)/4
    print(f"    Predicted ((p+1)/4)^(p-1) = {q**(p-1):.4f}")

# For general circulant: prod Q_k = |prod_{k=1}^{n-1} S_hat(k)|^2
# = |Res(Phi_n, f)|^2 where f(x) = sum_{s in S} x^s
# This is the resultant of the n-th cyclotomic polynomial with f.

# ============================================================
print(f"\n{'='*70}")
print("DEEP STRUCTURE: per-eigenspace chi and Morgan-Voyce")
print("="*70)

# For Interval n=2m+1, the per-eigenspace Omega are:
# omega_j^{per} = Omega_j / n
# What are these? Can they be expressed in terms of Morgan-Voyce polynomials?

# B_m(x) = sum_{j=0}^m C(m+j, 2j) x^j is the Morgan-Voyce polynomial
# B_m(1) = F_{2m+1}

# Per-eigenspace generating function?
for n in [5, 7, 9]:
    m = (n-1)//2
    S = set(range(1, m+1))
    A = circulant_tournament(n, S)

    per_eig = []
    for dim in range(n):
        paths = get_regular_paths(A, dim)
        per_eig.append(len(paths) // n)

    print(f"\n  Interval n={n} (m={m}):")
    print(f"    per-eig Omega: {per_eig}")
    print(f"    sum = {sum(per_eig)}")
    print(f"    F_{n} = ", end="")

    # Fibonacci numbers
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    F_n = a
    print(f"{F_n}")

    # Check: sum of per_eig_omega vs various Fibonacci-related quantities
    print(f"    sum/F_n = {sum(per_eig)/F_n:.6f}")

    # Alternating sum
    alt = sum((-1)**j * per_eig[j] for j in range(n))
    print(f"    alt sum = {alt}")

    # Check if per_eig_omega relate to binomial coefficients C(m,j) somehow
    for j in range(n):
        if j < n and per_eig[j] > 0:
            ratio = per_eig[j] / comb(n-1, j) if comb(n-1, j) > 0 else 0
            print(f"      omega_{j}^per = {per_eig[j]}, C({n-1},{j})={comb(n-1,j)}, ratio={ratio:.4f}")

print("\nDONE.")
