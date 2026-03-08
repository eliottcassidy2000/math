#!/usr/bin/env python3
"""
tournament_markov_chain.py — Tournament as Markov chain; spectrum and mixing.

IDEA: Define a random walk on {0,...,n-1} by: from vertex i, move to a
random out-neighbor of i (each with prob 1/out-degree(i)).

For a REGULAR tournament (out-degree = (n-1)/2 for all), the transition
matrix is P = (2/(n-1)) * A, where A is the adjacency matrix.

The stationary distribution is uniform (for regular tournaments).
The mixing time depends on the spectral gap.

CONNECTION TO F(T,x): The eigenvalues of A determine:
1. det(xI - A) = characteristic polynomial
2. The number of closed walks of length k = tr(A^k)
3. Cycle counts: c_k relates to tr(A^k) via Möbius function

The SIGNED adjacency matrix B = 2A - J + I (where J = all-ones) is
skew-symmetric: B + B^T = 0. Its eigenvalues are purely imaginary.
For Paley tournaments, these eigenvalues involve Gauss sums.

QUESTION: Does the spectrum of B (or A) encode F(T,x)?

F(T,x) involves PATHS, not cycles. The connection between path counting
and spectrum is via the permanent or via inclusion-exclusion.

NOVEL CONNECTION: The TRANSFER MATRIX of a 1D spin chain.
Think of positions 1,...,n as "sites" and the permutation P as a
"spin configuration" where P_i is the "spin" at site i.
The weight x^{A[P_i][P_{i+1}]} is the Boltzmann factor.
F(T,x) = partition function of this 1D model.

But the constraint "each spin value used exactly once" makes this
a HARD-CORE model (exclusion), not a simple transfer matrix model.

Author: opus-2026-03-07-S44
"""
import numpy as np
from itertools import permutations, combinations
import math
import random

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

# ============================================================
# SPECTRUM OF TOURNAMENT ADJACENCY MATRIX
# ============================================================
print("=" * 60)
print("EIGENVALUE SPECTRUM OF TOURNAMENT")
print("=" * 60)

for n in [5, 7]:
    m = n*(n-1)//2
    print(f"\n--- n={n} ---")

    seen = set()
    if n == 5:
        iterator = range(1 << m)
    else:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(100)]

    results = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        H = F[n-1]
        A_np = np.array(A, dtype=float)
        eigenvals = np.linalg.eigvals(A_np)

        # Signed adjacency B = 2A - J + I
        J = np.ones((n, n))
        I = np.eye(n)
        B = 2 * A_np - J + I
        B_eigenvals = np.linalg.eigvals(B)

        # Characteristic polynomial of A
        det_A = np.linalg.det(A_np)

        results.append({
            'H': H, 'F': F,
            'A_eigs': sorted(eigenvals, key=lambda z: -abs(z)),
            'B_eigs': sorted(B_eigenvals, key=lambda z: -abs(z)),
            'det_A': det_A
        })

    # Show examples
    for r in results[:8]:
        A_str = [f"{e.real:.2f}{e.imag:+.2f}i" for e in r['A_eigs']]
        B_str = [f"{e.real:.2f}{e.imag:+.2f}i" for e in r['B_eigs']]
        print(f"  H={r['H']:3d}: det(A)={r['det_A']:6.1f}")
        print(f"    A eigs: {A_str}")
        print(f"    B eigs: {B_str}")

# ============================================================
# NOVEL: det(A) vs H(T)
# ============================================================
print("\n" + "=" * 60)
print("det(A) vs H(T)")
print("=" * 60)

n = 5
m = n*(n-1)//2
seen = set()

det_H_pairs = []
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    A_np = np.array(A, dtype=float)
    det_A = round(np.linalg.det(A_np))
    det_H_pairs.append((H, det_A))

print(f"  n=5: H vs det(A):")
for H, d in sorted(det_H_pairs):
    print(f"    H={H:3d}, det(A)={d:4d}")

# n=7
n = 7
m = n*(n-1)//2
random.seed(42)
seen = set()
det_H_pairs = []

for trial in range(200):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    A_np = np.array(A, dtype=float)
    det_A = round(np.linalg.det(A_np))
    det_H_pairs.append((H, det_A))

# How many distinct det values?
det_vals = sorted(set(d for _, d in det_H_pairs))
print(f"\n  n=7: {len(det_vals)} distinct det(A) values: {det_vals[:20]}...")

# Does det(A) determine H?
det_to_H = {}
for H, d in det_H_pairs:
    if d not in det_to_H:
        det_to_H[d] = set()
    det_to_H[d].add(H)

multi = {d: sorted(hs) for d, hs in det_to_H.items() if len(hs) > 1}
print(f"  det(A) values mapping to multiple H: {len(multi)}")
if multi:
    for d, hs in sorted(multi.items())[:5]:
        print(f"    det(A)={d}: H in {hs}")

# ============================================================
# NOVEL: PFAFFIAN OF B vs H(T)
# ============================================================
print("\n" + "=" * 60)
print("PFAFFIAN OF B = 2A-J+I (SKEW-SYMMETRIC)")
print("=" * 60)

# For even n, Pf(B) is well-defined. Pf(B)^2 = det(B).
# For odd n, det(B) = 0 (skew-symmetric odd matrix).
# Check at n=6.

n = 6
m = n*(n-1)//2
seen = set()

def compute_F_fast(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

pf_H_pairs = []
for bits in range(min(1 << m, 500)):
    A = tournament_from_bits(bits, n)
    F = compute_F_fast(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    A_np = np.array(A, dtype=float)
    J = np.ones((n, n))
    I = np.eye(n)
    B = 2 * A_np - J + I
    det_B = round(np.linalg.det(B))
    # Pfaffian^2 = det_B
    if det_B >= 0:
        pf = round(math.sqrt(det_B))
    else:
        pf = None  # complex Pfaffian

    pf_H_pairs.append((H, det_B, pf))

print(f"  n=6: H vs det(B):")
for H, db, pf in sorted(pf_H_pairs)[:15]:
    pf_str = str(pf) if pf is not None else "complex"
    print(f"    H={H:3d}, det(B)={db:6d}, Pf={pf_str}")

# ============================================================
# NOVEL: ENERGY FUNCTION — SUM OF EIGENVALUES OF W(x)
# ============================================================
print("\n" + "=" * 60)
print("ENERGY = tr(W(x)) AND SPECTRAL CONNECTIONS")
print("=" * 60)
# W(x)[i][j] = x * A[i][j] + (1 - A[i][j]) for i ≠ j
# tr(W) = 0 (diagonal = 0)
# sum of all entries = sum_{i≠j} [x*A + (1-A)] = x*E + (m-E)
# where E = sum A[i][j] = C(n,2) and m = n*(n-1) pairs.
# So sum = x * C(n,2) + C(n,2) = C(n,2)(x+1).
# This is UNIVERSAL (tournament-independent)!

# What about tr(W^2) = sum_j sum_k W[j][k]*W[k][j]?
# W[j][k]*W[k][j] = x for all j≠k (multiplicative constraint!).
# So tr(W^2) = sum_{j≠k} x = n*(n-1)*x. Also universal!

# tr(W^3) = sum_{j,k,l} W[j][k]*W[k][l]*W[l][j]
# This involves 3-cycles and is NOT universal.

n = 5
m = n*(n-1)//2
seen = set()

print(f"\nn=5:")
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    H = F[n-1]
    x = 2  # symbolic would be better but let's use x=2
    W = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                W[i][j] = x * A[i][j] + (1 - A[i][j])

    tr_W2 = np.trace(W @ W)
    tr_W3 = np.trace(W @ W @ W)
    tr_W4 = np.trace(W @ W @ W @ W)

    print(f"  H={H:3d}: tr(W²)={tr_W2:8.1f}, tr(W³)={tr_W3:8.1f}, tr(W⁴)={tr_W4:8.1f}")

print(f"\n  tr(W²) = n*(n-1)*x = {n*(n-1)*2} (universal ✓)")
print(f"  tr(W³) depends on tournament (3-cycle count)")
print(f"  tr(W^k) for k>=3 encodes cycle structure — standard spectral theory")
