#!/usr/bin/env python3
"""
newton_identity_analysis.py -- kind-pasteur-2026-03-13-S60

The eigenvalue moments S_{2k} = sum D_t^{2k} are power sums of {D_1^2,...,D_m^2}.
Newton's identities relate power sums to elementary symmetric polynomials.

Key question: at p=11 (m=5), only 3 moments (S4,S6,S8) determine H.
But 5 power sums are generally needed to determine 5 values.
So either:
(a) The D^2 values satisfy extra constraints reducing the effective dimension, or
(b) H depends on fewer than all e_j's, or
(c) Some combination of these.

At p=13 (m=6), all 5 non-constant moments are needed, suggesting (a) at p=11.

This script:
1. Computes the D^2 multisets at p=7,11,13
2. Checks if D^2 values satisfy polynomial constraints (from S being a half-set of Z_p^*)
3. Tests whether the minimal polynomial of D^2 has lower degree at p=3 mod 4
4. Explores the Gauss sum connection
"""

import cmath
import numpy as np
from collections import defaultdict


def compute_D2_multiset(p, S):
    """Compute the multiset {D_1^2, ..., D_m^2} for a given orientation."""
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D2_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D2_vals.append(lam.imag ** 2)
    return sorted(D2_vals)


def compute_full_eigenvalues(p, S):
    """Compute all eigenvalues lambda_t = S_hat(t) for t=1,...,m."""
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    evals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        evals.append(lam)
    return evals


def gauss_sum_analysis(p, S):
    """Analyze the Gauss sum structure of S."""
    omega = cmath.exp(2j * cmath.pi / p)
    m = (p - 1) // 2

    # S_hat(t) = sum_{s in S} omega^{st}
    S_hat = []
    for t in range(p):
        val = sum(omega ** (s * t) for s in S)
        S_hat.append(val)

    # QR-related Gauss sum
    chi = {}  # Legendre symbol
    for a in range(p):
        if a == 0:
            chi[a] = 0
        elif pow(a, (p-1)//2, p) == 1:
            chi[a] = 1
        else:
            chi[a] = -1

    g = sum(chi[a] * omega ** a for a in range(p))  # quadratic Gauss sum
    g2 = g * g

    return S_hat, g, g2, chi


def analyze_D2_structure(p):
    """Analyze the structure of D^2 multisets at prime p."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"D^2 MULTISET STRUCTURE at p={p}, m={m}")
    print(f"{'='*70}")

    # Gauss sum
    omega = cmath.exp(2j * cmath.pi / p)
    chi = {}
    for a in range(p):
        if a == 0:
            chi[a] = 0
        elif pow(a, (p-1)//2, p) == 1:
            chi[a] = 1
        else:
            chi[a] = -1
    g = sum(chi[a] * omega ** a for a in range(p))
    print(f"  Gauss sum g = {g:.4f}")
    print(f"  g^2 = {g*g:.4f}")
    print(f"  |g|^2 = {abs(g)**2:.4f} (should be p={p})")
    print(f"  g^2/p = {g*g/p:.4f}")
    if p % 4 == 3:
        print(f"  p = 3 mod 4: g^2 = -p = {-p} (pure imaginary g)")
    else:
        print(f"  p = 1 mod 4: g^2 = +p = {p} (real g)")

    # Collect all D^2 multisets
    multisets = defaultdict(list)
    limit = min(N, 128)
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        D2 = tuple(round(d, 6) for d in compute_D2_multiset(p, S))
        multisets[D2].append(bits)

    n_unique = len(multisets)
    print(f"\n  {n_unique} unique D^2 multisets (out of {limit} orientations)")

    # Show all unique multisets
    print(f"\n  Unique D^2 multisets:")
    for D2 in sorted(multisets.keys()):
        cnt = len(multisets[D2])
        s2 = sum(D2)
        s4 = sum(d**2 for d in D2)
        print(f"    {[round(d, 4) for d in D2]}: count={cnt}, "
              f"sum={s2:.4f}, sum_sq={s4:.4f}")

    # Check if D^2 values lie on specific curves
    print(f"\n  --- D^2 value analysis ---")
    all_D2_vals = set()
    for D2 in multisets:
        for d in D2:
            all_D2_vals.add(round(d, 6))
    all_D2_sorted = sorted(all_D2_vals)
    print(f"  {len(all_D2_sorted)} distinct D^2 values across all orientations:")
    for d in all_D2_sorted:
        print(f"    D^2 = {d:.6f}")

    # For the PALEY tournament (if p = 3 mod 4):
    if p % 4 == 3:
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        D2_paley = compute_D2_multiset(p, S_qr)
        print(f"\n  Paley D^2 multiset: {[round(d, 6) for d in D2_paley]}")
        # All D^2 values should be (p+1)/4 for Paley
        target = (p + 1) / 4
        print(f"  Expected uniform value: (p+1)/4 = {target}")
        max_dev = max(abs(d - target) for d in D2_paley)
        print(f"  Max deviation from uniform: {max_dev:.8f}")

    # For the INTERVAL tournament:
    S_int = list(range(1, m + 1))
    D2_int = compute_D2_multiset(p, S_int)
    evals_int = compute_full_eigenvalues(p, S_int)
    print(f"\n  Interval D^2 multiset: {[round(d, 6) for d in D2_int]}")
    print(f"  Interval eigenvalues: {[f'{l:.4f}' for l in evals_int]}")

    # Check: does D_t depend on t in a simple way?
    print(f"\n  --- Eigenvalue-frequency relationship ---")
    for bits in [0, N-1, N//3]:
        if bits >= limit:
            continue
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        evals = compute_full_eigenvalues(p, S)
        print(f"\n    bits={bits}, S={S}:")
        for t in range(m):
            lam = evals[t]
            print(f"      t={t+1}: lambda = {lam.real:.4f} + {lam.imag:.4f}i, "
                  f"|lambda|^2 = {abs(lam)**2:.4f}, D^2 = {lam.imag**2:.4f}")

    # Key test: at p=11, why do 3 moments suffice for 5 D^2 values?
    # Possible: D^2 values satisfy a quadratic constraint
    if p in [11, 13]:
        print(f"\n  --- Polynomial constraints on D^2 values ---")
        for D2 in sorted(multisets.keys()):
            vals = list(D2)
            # Check if values satisfy a degree-d polynomial
            # P(x) = prod (x - D_i^2) should have constrained coefficients
            coeffs = np.polynomial.polynomial.polyfromroots(vals)
            print(f"    Minimal poly coeffs (from roots): "
                  f"{[round(c, 4) for c in coeffs]}")

    # Dimension of the D^2 variety
    print(f"\n  --- Effective dimension of D^2 space ---")
    # Collect all D^2 vectors (unsorted, by frequency order)
    D2_vectors = []
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        D2_unsorted = []
        omega_val = cmath.exp(2j * cmath.pi / p)
        for t in range(1, m + 1):
            lam = sum(omega_val ** (s * t) for s in S)
            D2_unsorted.append(lam.imag ** 2)
        D2_vectors.append(D2_unsorted)

    D2_mat = np.array(D2_vectors)
    # Center the matrix
    D2_centered = D2_mat - D2_mat.mean(axis=0)
    # SVD to find effective rank
    U, s, Vt = np.linalg.svd(D2_centered)
    print(f"  Singular values of centered D^2 matrix: "
          f"{[f'{sv:.4f}' for sv in s]}")
    rank = sum(1 for sv in s if sv > 1e-6)
    print(f"  Effective rank: {rank} (out of m={m})")
    print(f"  Fraction of variance: "
          f"{[f'{sv**2/sum(s**2):.4f}' for sv in s]}")

    return multisets


# ================================================================
# MAIN
# ================================================================
print("=" * 70)
print("NEWTON IDENTITY / D^2 MULTISET ANALYSIS")
print("=" * 70)

for p in [7, 11, 13]:
    analyze_D2_structure(p)

print("\nDONE.")
