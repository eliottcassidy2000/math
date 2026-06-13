#!/usr/bin/env python3
"""
tournament_spectral.py — Spectral decomposition of tournament matrices.

A tournament adjacency matrix A satisfies A + A^T = J - I.
So A = (J-I)/2 + S/2 where S = A - A^T is SKEW-SYMMETRIC.

S[i][j] = A[i][j] - A[j][i] = 2*A[i][j] - 1 for i≠j, 0 for i=i.
So S = 2A - (J-I) is the "signed tournament matrix."

Eigenvalues of S are purely imaginary: ±i*lambda_k.
Pfaffian of S is well-defined at even n.

QUESTIONS:
1. How do eigenvalues of S relate to F(T,x)?
2. Does Pf(S) have a tournament-theoretic meaning?
3. Is there a spectral formula for H(T) or F(T,x)?
4. Connection to the skew-Schur functions?

Also: A = (J-I+S)/2. The eigenvalues of A are:
  (n-1)/2 + pure imaginary (for the J eigenspace)
  -1/2 + pure imaginary (for the J-perp eigenspace)

So the REAL PART of all eigenvalues of A is either (n-1)/2 (once) or -1/2 (n-1 times).
This is a beautiful universal constraint!

Author: opus-2026-03-07-S46
"""
import numpy as np
from itertools import permutations, combinations
from math import factorial

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

def pfaffian(M):
    """Compute Pfaffian of a skew-symmetric matrix (even size)."""
    n = len(M)
    if n % 2 == 1:
        return 0
    if n == 0:
        return 1
    if n == 2:
        return M[0][1]
    # Expand along first row
    total = 0
    for j in range(1, n):
        if M[0][j] == 0:
            continue
        sign = (-1)**(j-1)
        # Remove rows/cols 0 and j
        remaining = [i for i in range(n) if i != 0 and i != j]
        sub = [[M[remaining[a]][remaining[b]] for b in range(len(remaining))] for a in range(len(remaining))]
        total += sign * M[0][j] * pfaffian(sub)
    return total

# ============================================================
# EIGENVALUE STRUCTURE OF TOURNAMENT MATRICES
# ============================================================
print("=" * 60)
print("EIGENVALUE STRUCTURE: A = (J-I+S)/2")
print("=" * 60)

for n in [4, 5, 6, 7]:
    m = n*(n-1)//2
    seen = set()
    real_parts = set()
    imag_data = []

    import random
    random.seed(42)
    num = min(1 << m, 30000)

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

        A = np.array(adj, dtype=float)
        eigs = np.linalg.eigvals(A)

        # Check real parts
        for e in eigs:
            real_parts.add(round(e.real, 6))

        # The imaginary parts (sorted by magnitude)
        imags = sorted([abs(e.imag) for e in eigs], reverse=True)
        imag_data.append((imags, F))

    print(f"\nn={n}: {len(seen)} distinct F-vectors")
    print(f"  Distinct real parts of eigenvalues: {sorted(real_parts)}")
    print(f"  Expected: {(n-1)/2} (once) and {-0.5} ({n-1} times)")

# ============================================================
# PFAFFIAN OF S = 2A - (J-I)
# ============================================================
print("\n" + "=" * 60)
print("PFAFFIAN OF SIGNED TOURNAMENT MATRIX S = 2A - J + I")
print("=" * 60)

for n in [4, 6]:
    m = n*(n-1)//2
    seen = set()
    pf_values = {}

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        # S = 2A - J + I
        S = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    S[i][j] = 0
                else:
                    S[i][j] = 2*adj[i][j] - 1  # +1 if i->j, -1 if j->i

        pf = pfaffian(S)
        H = sum(F)
        if pf not in pf_values:
            pf_values[pf] = []
        pf_values[pf].append(H)

    print(f"\nn={n}: Pfaffian values and corresponding H(T):")
    for pf_val in sorted(pf_values.keys()):
        H_vals = sorted(set(pf_values[pf_val]))
        print(f"  Pf(S) = {pf_val}: H values = {H_vals} ({len(pf_values[pf_val])} tournaments)")

# ============================================================
# SPECTRAL DETERMINANTS AND F(T,x)
# ============================================================
print("\n" + "=" * 60)
print("det(xI - A) vs F(T,x)")
print("=" * 60)

# The characteristic polynomial of A is det(xI - A).
# Can F(T,x) be expressed via the characteristic polynomial?
# F(T,x) has degree n-1, char poly has degree n.

# Alternative: det(I + x*A) or det(I - x*A)?
# det(I + x*A) = prod(1 + x*lambda_i) where lambda_i are eigenvalues.

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")
    seen = set()

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        A = np.array(adj, dtype=float)

        # det(I + x*A) at several x values
        det_vals = {}
        F_vals = {}
        for x in [0.5, 1, 2, 3]:
            det_vals[x] = np.linalg.det(np.eye(n) + x * A)
            F_vals[x] = sum(F[k] * x**k for k in range(n))

        if len(seen) <= 5:
            H = sum(F)
            print(f"  F={F}, H={H}")
            print(f"  det(I+xA) at x=0.5,1,2,3: {[f'{det_vals[x]:.1f}' for x in [0.5,1,2,3]]}")
            print(f"  F(x) at x=0.5,1,2,3: {[f'{F_vals[x]:.1f}' for x in [0.5,1,2,3]]}")

# ============================================================
# KEY QUESTION: Is Im(eigenvalue) determined by t3?
# ============================================================
print("\n" + "=" * 60)
print("IMAGINARY EIGENVALUES vs t3")
print("=" * 60)

for n in [5]:
    m = n*(n-1)//2
    seen = set()
    t3_to_imag = {}

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        A = np.array(adj, dtype=float)
        eigs = np.linalg.eigvals(A)
        imags = tuple(sorted([round(abs(e.imag), 6) for e in eigs], reverse=True))

        t3 = sum(1 for triple in combinations(range(n), 3)
                 if (adj[triple[0]][triple[1]] and adj[triple[1]][triple[2]] and adj[triple[2]][triple[0]]) or
                    (adj[triple[0]][triple[2]] and adj[triple[2]][triple[1]] and adj[triple[1]][triple[0]]))

        if t3 not in t3_to_imag:
            t3_to_imag[t3] = set()
        t3_to_imag[t3].add(imags)

    print(f"\nn={n}:")
    for t3 in sorted(t3_to_imag.keys()):
        imag_sets = t3_to_imag[t3]
        print(f"  t3={t3}: {len(imag_sets)} distinct imaginary spectra")
        for spec in sorted(imag_sets)[:3]:
            print(f"    {spec}")

# ============================================================
# TOURNAMENT AS ELEMENT OF CLIFFORD ALGEBRA?
# ============================================================
print("\n" + "=" * 60)
print("SIGNED ADJACENCY: S^2 EIGENVALUES")
print("=" * 60)

# S = 2A - J + I is skew-symmetric. S^2 is symmetric negative semidefinite.
# Eigenvalues of S^2 are -lambda_k^2 where lambda_k are imag parts of S eigenvalues.
# S^2[i][j] = sum_k S[i][k]*S[k][j]

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

        S = np.array([[2*adj[i][j]-1 if i!=j else 0 for j in range(n)] for i in range(n)], dtype=float)
        S2 = S @ S
        eigs_S2 = sorted(np.linalg.eigvals(S2).real, reverse=True)

        # S^2[i][i] = sum_j S[i][j]^2 = sum_{j≠i} 1 = n-1 (since S[i][j]=±1 for i≠j)
        # So tr(S^2) = n(n-1)... wait, that's -n(n-1) since S^2 is neg semidef
        # Actually S^2[i][i] = -sum_j S[i][j]^2 = -(n-1) ... no.
        # S[i][j] = ±1 for i≠j, S[i][i]=0
        # (S^2)[i][i] = sum_j S[i][j]*S[j][i] = sum_j S[i][j]*(-S[i][j]) = -sum_{j≠i} 1 = -(n-1)
        # So tr(S^2) = -n(n-1), confirmed.

        if len(seen) <= 5:
            H = sum(F)
            print(f"  F={F}, H={H}")
            print(f"  eigs(S^2) = {[f'{e:.2f}' for e in eigs_S2]}")
            print(f"  tr(S^2) = {sum(eigs_S2):.0f} (should be {-n*(n-1)})")
