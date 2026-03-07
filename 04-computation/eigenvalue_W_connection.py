#!/usr/bin/env python3
"""
Connect tournament eigenvalues to W(r) coefficients.

The skew part S of tournament adjacency A = (J-I)/2 + S has purely imaginary eigenvalues.
The Cayley transform Q = (I+S)(I-S)^{-1} has eigenvalues on the unit circle.

Question: do the W(r) coefficients relate to symmetric functions of the eigenvalues?

We know: W(1/2) = H(T), W(-1/2) = (-1)^{n-1} H(T).
Irving-Omar: W_IO(z) = det(I+zA^T)/det(I-zA) = prod(1+z*mu_k)/(1-z*mu_k)

Are the roots of W(r) related to eigenvalues?

kind-pasteur-2026-03-06-S25f
"""

import numpy as np
from itertools import permutations
import random
from math import factorial

def tournament_random(n, seed):
    random.seed(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def paley(p):
    """Paley tournament on F_p."""
    qr = set()
    for x in range(1, p):
        qr.add((x*x) % p)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j-i) % p) in qr:
                A[i][j] = 1
    return A

def W_coefficients(A):
    n = len(A)
    coeffs = [0.0] * n
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        poly = [1.0]
        for si in s:
            new_poly = [0.0] * (len(poly) + 1)
            for j, c in enumerate(poly):
                new_poly[j+1] += c
                new_poly[j] += c * si
            poly = new_poly
        for k in range(min(len(poly), n)):
            coeffs[k] += poly[k]
    return coeffs

def skew_part(A):
    n = len(A)
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            S[i, j] = (A[i][j] - A[j][i]) / 2
    return S

print("=" * 70)
print("EIGENVALUE-W(r) CONNECTION")
print("=" * 70)

for n in [5, 7]:
    print(f"\n--- n={n} ---")

    # Use Paley tournament if prime, else random
    if n in [3, 5, 7, 11]:
        A = paley(n)
        label = "Paley"
    else:
        A = tournament_random(n, 42)
        label = "random"

    S = skew_part(A)
    eigs_S = np.linalg.eigvals(S)
    # Sort by imaginary part
    eigs_S = sorted(eigs_S, key=lambda x: x.imag)

    # Adjacency matrix eigenvalues
    A_np = np.array(A, dtype=float)
    eigs_A = np.linalg.eigvals(A_np)
    eigs_A = sorted(eigs_A, key=lambda x: x.real)

    print(f"  Tournament: {label}")
    print(f"  Eigenvalues of S: {['%.3f*i' % e.imag for e in eigs_S]}")
    print(f"  Eigenvalues of A: {['%.3f+%.3fi' % (e.real, e.imag) for e in eigs_A]}")

    # W coefficients
    wc = W_coefficients(A)
    print(f"  W(r) coefficients: {[f'{c:.1f}' for c in wc]}")

    # Verify W(1/2) = H
    H = sum(1 for p in permutations(range(n))
            if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))
    W_half = sum(wc[k] * (0.5)**k for k in range(n))
    print(f"  H = {H}, W(1/2) = {W_half:.4f}")

    # For Paley: check if eigenvalues have special structure
    if label == "Paley":
        print(f"\n  Paley eigenvalue analysis:")
        # Paley tournament has eigenvalues related to Gauss sums
        # For p prime: eigenvalues of A_Paley are (p-1)/2 (with multiplicity 1)
        # and (-1 +/- sqrt(p))/2 (each with multiplicity (p-1)/2)
        print(f"  Expected: (p-1)/2 = {(n-1)/2}")
        print(f"  Expected: (-1+sqrt(p))/2 = {(-1+np.sqrt(n))/2:.4f}")
        print(f"  Expected: (-1-sqrt(p))/2 = {(-1-np.sqrt(n))/2:.4f}")
        print(f"  Each non-trivial eigenvalue has multiplicity (p-1)/2 = {(n-1)//2}")

# =====================================================================
# Now: relate W(r) roots to eigenvalues
# =====================================================================
print(f"\n{'='*70}")
print("W(r) ROOTS vs EIGENVALUES")
print(f"{'='*70}")

for n in [5, 7]:
    print(f"\n--- n={n} ---")

    for seed in range(5):
        A = tournament_random(n, seed * 13)
        S = skew_part(A)
        eigs_S = np.sort(np.linalg.eigvals(S).imag)

        wc = W_coefficients(A)
        # W(r) is degree n-1 with only even powers
        # So W(r) = w_0 + w_2*r^2 + w_4*r^4 + ...
        # The roots come in +-pairs
        # W(r) as polynomial of r^2: W(u) = w_0 + w_2*u + w_4*u^2 + ...
        # This is a polynomial of degree (n-1)/2 in u = r^2

        # Find roots of W(r) numerically
        poly_coeffs = list(reversed(wc))  # numpy wants highest degree first
        roots = np.roots(poly_coeffs)
        real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-10])
        complex_roots = sorted([r for r in roots if abs(r.imag) > 1e-10], key=lambda r: r.real)

        print(f"  seed={seed}: W roots (real): {[f'{r:.4f}' for r in real_roots]}")
        if complex_roots:
            print(f"          W roots (complex): {[f'{r.real:.4f}+{r.imag:.4f}i' for r in complex_roots[:4]]}")
        print(f"          S eigenvalues (imag): {[f'{e:.4f}' for e in eigs_S if abs(e) > 0.01]}")

# =====================================================================
# Key insight check: at Paley, W(r) should be especially structured
# =====================================================================
print(f"\n{'='*70}")
print("PALEY W(r) STRUCTURE")
print(f"{'='*70}")

for p in [5, 7]:
    A = paley(p)
    wc = W_coefficients(A)
    H = sum(1 for perm in permutations(range(p))
            if all(A[perm[i]][perm[i+1]] == 1 for i in range(p-1)))

    print(f"\n  Paley p={p}:")
    print(f"    H = {H}")
    print(f"    W coefficients: {[f'{c:.4f}' for c in wc]}")

    # Ratios
    if wc[-1] != 0:
        print(f"    Ratios w_k/w_{{n-1}}: {[f'{c/wc[-1]:.6f}' for c in wc]}")

    # Check: are W coefficients related to binomial coefficients?
    # For the "maximally symmetric" tournament, W might be a power of (r + something)
    # Or a Chebyshev polynomial

    # Check if W(r)/n! is a polynomial with nice coefficients
    print(f"    W(r)/n!: {[f'{c/factorial(p):.6f}' for c in wc]}")

    # For Paley T_5: H = 15
    if p == 5:
        print(f"\n    At p=5: W(r) = {wc[0]:.0f} + {wc[2]:.0f}*r^2 + {wc[4]:.0f}*r^4")
        print(f"    W(r)/120 = {wc[0]/120:.6f} + {wc[2]/120:.6f}*r^2 + {wc[4]/120:.6f}*r^4")
        print(f"    = 1 + (1/4)*r^2 + r^4  ??")
        check = [wc[0]/120, wc[2]/120, wc[4]/120]
        print(f"    Actually: [{check[0]:.4f}, {check[1]:.4f}, {check[2]:.4f}]")

print("\nDONE")
