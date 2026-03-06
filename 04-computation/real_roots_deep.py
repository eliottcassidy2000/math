#!/usr/bin/env python3
"""
Real Roots Deep Analysis
=========================
Compute full independence polynomial of Omega(T) and analyze roots.

Tests:
1. All roots real and negative? (n=5 exhaustive, n=6-8 random)
2. Newton inequalities: alpha_k^2 >= alpha_{k-1}*alpha_{k+1}*(k+1)*(m-k+1)/((m-k)*k)?
3. Ultra-log-concavity: (alpha_k/C(m,k))^2 >= (alpha_{k-1}/C(m,k-1))*(alpha_{k+1}/C(m,k+1))?
4. Coefficient patterns for Paley tournaments

kind-pasteur-2026-03-05-S16
"""

import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import (tournament_from_bits, all_tournaments,
                             random_tournament, find_odd_cycles,
                             conflict_graph, hamiltonian_path_count)
from itertools import combinations
from math import comb, factorial
import random

def independence_polynomial_coefficients(adj):
    """Compute full independence polynomial coefficients [alpha_0, alpha_1, ..., alpha_d].
    adj is adjacency matrix. Returns list where coeffs[k] = alpha_k."""
    m = len(adj)
    if m == 0:
        return [1]

    # Precompute neighbor bitmasks
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j

    # Count independent sets by size
    max_alpha = m
    coeffs = [0] * (max_alpha + 1)

    for mask in range(1 << m):
        # Check independence
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            k = bin(mask).count('1')
            coeffs[k] += 1

    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()

    return coeffs

def independence_poly_coeffs_fast(cycles):
    """Fast computation using vertex-disjoint cycle collections."""
    m = len(cycles)
    if m == 0:
        return [1]

    vsets = [frozenset(c) for c in cycles]

    # Build adjacency bitmasks
    adj_bits = [0] * m
    for a in range(m):
        for b in range(a + 1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    # For small m, enumerate all subsets
    if m <= 20:
        max_k = 0
        coeffs_dict = {}
        for mask in range(1 << m):
            ok = True
            seen = 0
            temp = mask
            while temp:
                v = (temp & -temp).bit_length() - 1
                if adj_bits[v] & seen:
                    ok = False
                    break
                seen |= 1 << v
                temp &= temp - 1
            if ok:
                k = bin(mask).count('1')
                coeffs_dict[k] = coeffs_dict.get(k, 0) + 1
                max_k = max(max_k, k)

        coeffs = [coeffs_dict.get(k, 0) for k in range(max_k + 1)]
        return coeffs
    else:
        # For larger m, use greedy approach to find max independent set size
        # then enumerate by size
        n_verts = len(set(v for c in cycles for v in c))
        max_k = n_verts // 3

        coeffs = [0] * (max_k + 1)
        coeffs[0] = 1
        coeffs[1] = m

        # Size 2
        for a in range(m):
            for b in range(a+1, m):
                if not (adj_bits[a] & (1 << b)):
                    coeffs[2] = coeffs.get(2, 0) + 1 if isinstance(coeffs, dict) else coeffs[2] + 1

        # Size 3
        if max_k >= 3:
            for a in range(m):
                for b in range(a+1, m):
                    if adj_bits[a] & (1 << b):
                        continue
                    for c_idx in range(b+1, m):
                        if not (adj_bits[a] & (1 << c_idx)) and not (adj_bits[b] & (1 << c_idx)):
                            coeffs[3] += 1

        # Trim
        while len(coeffs) > 1 and coeffs[-1] == 0:
            coeffs.pop()

        return coeffs

def check_real_roots(coeffs):
    """Check if all roots of polynomial with given coefficients are real.
    Returns (all_real, roots)."""
    import numpy as np
    if len(coeffs) <= 1:
        return True, True, []

    # Polynomial is sum alpha_k * x^k
    # numpy.roots expects [a_n, a_{n-1}, ..., a_0]
    poly_coeffs = list(reversed(coeffs))
    roots = np.roots(poly_coeffs)

    # Check if all roots are real (imaginary part < epsilon)
    epsilon = 1e-6
    all_real = all(abs(r.imag) < epsilon for r in roots)
    all_neg = all(r.real < epsilon for r in roots) if all_real else False

    return all_real, all_neg, roots

def check_ulc(coeffs, m):
    """Check ultra-log-concavity: alpha_k/C(m,k) is log-concave."""
    d = len(coeffs) - 1
    if d <= 1:
        return True

    # Compute normalized sequence
    norm = []
    for k in range(d + 1):
        c = comb(m, k)
        if c == 0:
            norm.append(float('inf'))
        else:
            norm.append(coeffs[k] / c)

    # Check log-concavity: norm[k]^2 >= norm[k-1] * norm[k+1]
    for k in range(1, d):
        if norm[k-1] == float('inf') or norm[k] == float('inf') or norm[k+1] == float('inf'):
            continue
        if norm[k] ** 2 < norm[k-1] * norm[k+1] - 1e-10:
            return False

    return True

def check_newton(coeffs):
    """Check Newton inequalities: alpha_k^2 >= alpha_{k-1} * alpha_{k+1} * (k+1)(m-k+1)/((m-k)*k)
    Actually the standard Newton: alpha_k^2 * C(d,k)^{-2} >= alpha_{k-1}*alpha_{k+1} * C(d,k-1)^{-1} * C(d,k+1)^{-1}

    Simplified: for real-rooted poly, alpha_k^2 >= alpha_{k-1}*alpha_{k+1} * (k+1)/(k) * (d-k+1)/(d-k)
    """
    d = len(coeffs) - 1
    if d <= 1:
        return True
    for k in range(1, d):
        if coeffs[k-1] == 0 and coeffs[k+1] == 0:
            continue
        lhs = coeffs[k] ** 2
        rhs = coeffs[k-1] * coeffs[k+1] * (k+1) * (d - k + 1) / (k * (d - k))
        if lhs < rhs - 0.01:
            return False
    return True

# ── Exhaustive n=5 ──
print("=" * 70)
print("REAL ROOTS ANALYSIS — n=5 EXHAUSTIVE")
print("=" * 70)

n = 5
m_tiles = n * (n - 1) // 2
total = 0
real_root_failures = 0
ulc_failures = 0
newton_failures = 0
coeff_histogram = {}

for bits in range(1 << m_tiles):
    T = tournament_from_bits(n, bits)
    cycles = find_odd_cycles(T)
    m = len(cycles)

    if m <= 20:
        cg = conflict_graph(cycles)
        coeffs = independence_polynomial_coefficients(cg)
    else:
        coeffs = independence_poly_coeffs_fast(cycles)

    key = tuple(coeffs)
    coeff_histogram[key] = coeff_histogram.get(key, 0) + 1

    all_real, all_neg, roots = check_real_roots(coeffs)
    if not all_real:
        real_root_failures += 1

    if not check_ulc(coeffs, m):
        ulc_failures += 1

    if not check_newton(coeffs):
        newton_failures += 1

    total += 1

print(f"Total: {total}")
print(f"Real root failures: {real_root_failures}")
print(f"ULC failures: {ulc_failures}")
print(f"Newton failures: {newton_failures}")
print(f"\nDistinct coefficient patterns: {len(coeff_histogram)}")
print("Pattern (alpha_0, alpha_1, ...) : count")
for coeffs, count in sorted(coeff_histogram.items(), key=lambda x: -x[1]):
    h = sum(c * 2**k for k, c in enumerate(coeffs))
    print(f"  {list(coeffs)} -> H={h}, count={count}")

# ── n=6 random ──
print(f"\n{'='*70}")
print("REAL ROOTS ANALYSIS — n=6 (200 random)")
print("='*70")

random.seed(42)
real6 = 0
ulc6 = 0
newton6 = 0
coeff6 = {}

for trial in range(200):
    T = random_tournament(6)
    cycles = find_odd_cycles(T)
    m = len(cycles)
    cg = conflict_graph(cycles)
    coeffs = independence_polynomial_coefficients(cg)

    key = tuple(coeffs)
    coeff6[key] = coeff6.get(key, 0) + 1

    all_real, all_neg, roots = check_real_roots(coeffs)
    if not all_real:
        real6 += 1
    if not check_ulc(coeffs, m):
        ulc6 += 1
    if not check_newton(coeffs):
        newton6 += 1

print(f"Real root failures: {real6}/200")
print(f"ULC failures: {ulc6}/200")
print(f"Newton failures: {newton6}/200")
print(f"Distinct patterns: {len(coeff6)}")
for coeffs, count in sorted(coeff6.items(), key=lambda x: -x[1])[:10]:
    h = sum(c * 2**k for k, c in enumerate(coeffs))
    d = len(coeffs) - 1
    print(f"  {list(coeffs)} -> H={h}, deg={d}, count={count}")

# ── Paley T_7 detailed ──
print(f"\n{'='*70}")
print("PALEY T_7: FULL INDEPENDENCE POLYNOMIAL OF OMEGA")
print("=" * 70)

def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a * a) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

T7 = paley_tournament(7)
cycles7 = find_odd_cycles(T7)
m7 = len(cycles7)
print(f"Omega(T_7): {m7} vertices (cycles)")

# Count by length
len_counts = {}
for c in cycles7:
    L = len(c)
    len_counts[L] = len_counts.get(L, 0) + 1
print(f"By length: {len_counts}")

# Build conflict graph
cg7 = conflict_graph(cycles7)
print(f"Computing independence polynomial coefficients...")

coeffs7 = independence_polynomial_coefficients(cg7)
print(f"I(Omega(T_7), x) coefficients: {coeffs7}")

h7 = sum(c * 2**k for k, c in enumerate(coeffs7))
print(f"I(Omega(T_7), 2) = {h7} (should be 189)")
print(f"H(T_7) = {hamiltonian_path_count(T7)}")

all_real7, all_neg7, roots7 = check_real_roots(coeffs7)
print(f"All roots real? {all_real7}")
print(f"All roots negative? {all_neg7}")
if roots7 is not None and len(roots7) > 0:
    real_roots = sorted([r.real for r in roots7])
    print(f"Roots: {[f'{r:.6f}' for r in real_roots]}")

ulc7 = check_ulc(coeffs7, m7)
print(f"ULC? {ulc7}")

# Normalized coefficients
print(f"\nNormalized alpha_k/C({m7},k):")
for k, a in enumerate(coeffs7):
    c = comb(m7, k)
    print(f"  k={k}: alpha={a}, C({m7},{k})={c}, ratio={a/c:.6f}")

# ── n=7 random sampling ──
print(f"\n{'='*70}")
print("n=7: 50 random tournaments (3-cycle subgraph)")
print("=" * 70)

real7_fail = 0
ulc7_fail = 0

for trial in range(50):
    T = random_tournament(7)
    # Use only 3-cycles for speed
    cycles = [c for c in find_odd_cycles(T) if len(c) == 3]
    m = len(cycles)
    if m == 0:
        continue
    cg = conflict_graph(cycles)
    coeffs = independence_polynomial_coefficients(cg)

    all_real, all_neg, roots = check_real_roots(coeffs)
    if not all_real:
        real7_fail += 1
        print(f"  REAL ROOT FAILURE at trial {trial}! coeffs={coeffs}")

    if not check_ulc(coeffs, m):
        ulc7_fail += 1
        print(f"  ULC FAILURE at trial {trial}! coeffs={coeffs}")

print(f"Real root failures: {real7_fail}/50")
print(f"ULC failures: {ulc7_fail}/50")

print("\nDone.")
