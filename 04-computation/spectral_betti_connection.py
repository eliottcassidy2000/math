#!/usr/bin/env python3
"""
spectral_betti_connection.py - Do eigenvalues of skew-adjacency predict path homology?

The skew-adjacency matrix S of a tournament T has S[i][j] = 1 if i→j, -1 if j→i, 0 on diag.
S is skew-symmetric, so eigenvalues are purely imaginary: ±iλ₁, ±iλ₂, ...

QUESTION: Does the spectral structure correlate with β₁, β₃?

Possible connections:
1. det(S) relates to cycle structure (Pfaffian)
2. Spectral gap might predict "topological phase"
3. Eigenvalue distribution might encode the same "odd-only" structure

Also: the CHARACTERISTIC POLYNOMIAL of S has the form det(xI - S) = x^n + c₂x^{n-2} + c₄x^{n-4} + ...
(only even powers, since S is skew-symmetric). The coefficients c₂ₖ are sums of principal
2k×2k Pfaffians, which count... signed sums of vertex-disjoint cycle covers!

This connects to α_k in the OCF!

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations
import numpy as np
import random

random.seed(42)

def skew_adj(A, n):
    """Skew-adjacency matrix: S[i][j] = A[i][j] - A[j][i]"""
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def spectral_invariants(S, n):
    """Extract eigenvalue-based invariants."""
    eigs = np.linalg.eigvals(S)
    # Eigenvalues should be purely imaginary
    imag_parts = sorted(np.abs(eigs.imag), reverse=True)
    real_parts = np.abs(eigs.real)

    # Pfaffian = sqrt(det(S)) for even n (up to sign)
    det_S = np.linalg.det(S)

    # Characteristic polynomial coefficients (from eigenvalues)
    # det(xI - S) = product of (x - λ_i)

    return {
        'det': det_S,
        'trace_S2': np.trace(S @ S),  # = -2 * sum of edges = -n(n-1) for tournament
        'max_imag': imag_parts[0] if len(imag_parts) > 0 else 0,
        'spectral_gap': imag_parts[0] - imag_parts[1] if len(imag_parts) > 1 else 0,
        'imag_parts': imag_parts[:4],  # top 4
    }

# === n=5: exhaustive ===
print("=" * 60)
print("n=5: SPECTRAL vs β₁ (exhaustive)")
print("=" * 60)

n = 5
m = n*(n-1)//2
spectral_by_beta = {0: [], 1: []}

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0

    S = skew_adj(A, n)
    spec = spectral_invariants(S, n)
    spectral_by_beta[b1].append(spec)

for b1_val in [0, 1]:
    specs = spectral_by_beta[b1_val]
    dets = [s['det'] for s in specs]
    max_imags = [s['max_imag'] for s in specs]
    gaps = [s['spectral_gap'] for s in specs]

    print(f"\nβ₁={b1_val}: {len(specs)} tournaments")
    print(f"  det(S): range [{min(dets):.1f}, {max(dets):.1f}], mean={np.mean(dets):.2f}")
    print(f"  max|λ|: range [{min(max_imags):.2f}, {max(max_imags):.2f}], mean={np.mean(max_imags):.2f}")
    print(f"  spectral gap: range [{min(gaps):.2f}, {max(gaps):.2f}], mean={np.mean(gaps):.2f}")

# Check: does det(S) = 0 iff some topological condition?
# For n odd, det(S) = 0 always (skew-symmetric odd matrix).
# For n even, det(S) = Pf(S)² ≠ 0 generically.
print(f"\nn={n} is odd, so det(S) = 0 for all tournaments (confirmed: {all(abs(s['det']) < 1e-10 for specs in spectral_by_beta.values() for s in specs)})")

# === n=6: spectral vs β₃ ===
print("\n" + "=" * 60)
print("n=6: SPECTRAL vs β₁, β₃ (sampling 500)")
print("=" * 60)

n = 6
spectral_by_type = {'b0': [], 'b1': [], 'b3': []}

for trial in range(500):
    A = random_tournament(n)
    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    S = skew_adj(A, n)
    spec = spectral_invariants(S, n)

    if b3 > 0:
        spectral_by_type['b3'].append(spec)
    elif b1 > 0:
        spectral_by_type['b1'].append(spec)
    else:
        spectral_by_type['b0'].append(spec)

for key in ['b0', 'b1', 'b3']:
    specs = spectral_by_type[key]
    if not specs:
        print(f"\n{key}: 0 tournaments")
        continue
    dets = [abs(s['det']) for s in specs]
    max_imags = [s['max_imag'] for s in specs]
    gaps = [s['spectral_gap'] for s in specs]

    print(f"\n{key}: {len(specs)} tournaments")
    print(f"  |det(S)|: range [{min(dets):.1f}, {max(dets):.1f}], mean={np.mean(dets):.2f}")
    print(f"  max|λ|: range [{min(max_imags):.2f}, {max(max_imags):.2f}], mean={np.mean(max_imags):.2f}")
    print(f"  spectral gap: range [{min(gaps):.2f}, {max(gaps):.2f}], mean={np.mean(gaps):.2f}")

# === Pfaffian connection ===
print("\n" + "=" * 60)
print("PFAFFIAN CONNECTION")
print("=" * 60)

print("""
For even n, det(S) = Pf(S)². The Pfaffian is:
  Pf(S) = sum over perfect matchings M: sgn(M) * product_{(i,j) in M} S[i][j]

Each S[i][j] = ±1 (for i≠j in tournament), so Pf(S) is a signed count of matchings.
This is related to the cycle structure since each matching decomposes into...
wait, matchings in the COMPLETE graph, not the tournament.

Actually: the Pfaffian of the SKEW-ADJACENCY of a tournament counts
"tournament-adapted" perfect matchings. For n=4:
  Pf(S) = S[0][1]*S[2][3] - S[0][2]*S[1][3] + S[0][3]*S[1][2]

Each factor is ±1, so Pf(S) ∈ {-3, -1, 1, 3} for n=4.

KEY INSIGHT: det(S) = Pf(S)² measures the "matching richness" of the tournament.
High |Pf(S)| → many consistent matchings → possible OCF connection?

For n=6: |Pf(S)| can be much larger, and β₃ > 0 might correlate with |Pf(S)|.
""")

# Check Pfaffian distribution by β type at n=6
if spectral_by_type['b3']:
    print("Pfaffian check at n=6 (|det(S)| = Pf(S)²):")
    for key in ['b0', 'b1', 'b3']:
        specs = spectral_by_type[key]
        if specs:
            dets = sorted(set(round(abs(s['det'])) for s in specs))
            print(f"  {key}: |det(S)| values = {dets[:10]}{'...' if len(dets) > 10 else ''}")

# === Characteristic polynomial ===
print("\n" + "=" * 60)
print("CHARACTERISTIC POLYNOMIAL INVARIANTS")
print("=" * 60)

# For skew-symmetric S, char poly = x^n + c₂x^{n-2} + c₄x^{n-4} + ...
# c₂ₖ = sum of Pf²(S_I) where I ranges over 2k-element subsets
# This is related to cycle covers!

# Let me compute char poly coefficients for specific n=6 tournaments
print("\nn=6: Char poly coefficients by β type (sample):")
for key in ['b0', 'b1', 'b3']:
    if spectral_by_type[key]:
        # Take a representative
        spec = spectral_by_type[key][0]
        # Recompute char poly
        # We need the actual tournament, not just the spectral data
        pass

# Actually, let me recompute with stored tournaments
print("\nRecomputing with stored adjacency matrices...")

n = 6
examples = {'b0': None, 'b1': None, 'b3': None}
for trial in range(2000):
    A = random_tournament(n)
    beta = path_betti_numbers(A, n)
    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    key = 'b3' if b3 > 0 else ('b1' if b1 > 0 else 'b0')
    if examples[key] is None:
        S = skew_adj(A, n)
        coeffs = np.polynomial.polynomial.polyfromroots(np.linalg.eigvals(S))
        # char poly of S
        char_poly = np.real(np.round(coeffs, 2))
        examples[key] = {'A': A, 'S': S, 'char_poly': char_poly}
        print(f"  {key}: char poly coeffs (low to high) = {char_poly}")

    if all(v is not None for v in examples.values()):
        break

print("\n" + "=" * 60)
print("INTERPRETATION")
print("=" * 60)
print("""
The skew-adjacency matrix S of a tournament encodes:
1. det(S) = Pf(S)² — "matching complexity"
2. Eigenvalues iλ₁, ..., iλ_n — "spectral shape"
3. Char poly coefficients — signed cycle cover counts

CONJECTURE: β₃ > 0 tournaments have LARGER |Pf(S)| (more consistent matchings).
This would connect path homology (topology) to spectral theory (algebra)
and OCF (combinatorics) through the Pfaffian.

The three pillars would then be FOUR pillars:
1. OCF (I(Ω,2))
2. Cumulant hierarchy (κ_{2k})
3. Path homology (β_{2k-1})
4. Spectral (eigenvalues of skew-adj)
""")
