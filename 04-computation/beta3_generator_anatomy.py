#!/usr/bin/env python3
"""
ANATOMY OF β₃ GENERATORS
=========================
For tournaments with β₃ = 1, find explicit generators of H₃ and
analyze their structure. What makes them non-trivial?

opus-2026-03-15-S72d
"""
import numpy as np
import sys
import random
from collections import Counter

sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')
from path_homology_v2 import (enumerate_allowed_paths, boundary_coeffs,
                                build_full_boundary_matrix, compute_omega_basis,
                                path_betti_numbers)

def tournament_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                wins = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
                if max(wins) < 2:
                    c3 += 1
    return c3

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

def find_h3_generator(A, n):
    """Find a non-trivial H₃ generator if β₃ = 1."""
    allowed = {}
    for p in range(-1, min(n, 6)):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    omega = {}
    for p in range(min(n, 6)):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    # ∂₃: Ω₃ → Ω₂
    dim_omega3 = omega[3].shape[1] if omega[3].ndim == 2 else 0
    if dim_omega3 == 0:
        return None, None, None

    bd3 = build_full_boundary_matrix(allowed[3], allowed[2])
    bd3_omega = bd3 @ omega[3]  # in A₂ coordinates

    # ker(∂₃) on Ω₃
    S3 = np.linalg.svd(bd3_omega, compute_uv=False)
    rank3 = sum(s > 1e-8 for s in S3)
    ker3_dim = dim_omega3 - rank3

    # im(∂₄) on Ω₃
    dim_omega4 = omega[4].shape[1] if 4 in omega and omega[4].ndim == 2 else 0
    if dim_omega4 > 0:
        bd4 = build_full_boundary_matrix(allowed[4], allowed[3])
        bd4_omega = bd4 @ omega[4]
        S4 = np.linalg.svd(bd4_omega, compute_uv=False)
        im4_dim = sum(s > 1e-8 for s in S4)
    else:
        im4_dim = 0

    beta3 = ker3_dim - im4_dim
    if beta3 <= 0:
        return None, None, None

    # Find kernel basis of ∂₃|Ω₃
    U3, S3_full, Vt3 = np.linalg.svd(bd3_omega, full_matrices=True)
    rank3 = sum(s > 1e-8 for s in S3_full)
    ker_basis = Vt3[rank3:].T  # columns in Ω₃ coordinates

    # Convert to A₃ coordinates
    ker_in_A3 = omega[3] @ ker_basis

    # If there's an image from ∂₄, project out
    if im4_dim > 0:
        im4_in_A3 = bd4 @ omega[4]  # in A₃... wait, bd4 maps to A₃
        # Actually im(∂₄) is in A₃. We need im in Ω₃ coordinates.
        # im(∂₄|Ω₄) in A₃ coordinates
        im4_A3 = (bd4 @ omega[4])
        # Project ker_basis out of im(∂₄)
        # We work in A₃ space
        # Find basis of im(∂₄) in A₃
        U_im, S_im, Vt_im = np.linalg.svd(im4_A3, full_matrices=False)
        nz = sum(s > 1e-8 for s in S_im)
        if nz > 0:
            P_im = U_im[:, :nz]  # basis of im(∂₄) in A₃
            # Project out
            for col in range(ker_in_A3.shape[1]):
                v = ker_in_A3[:, col]
                proj = P_im @ (P_im.T @ v)
                ker_in_A3[:, col] = v - proj

    # Find nonzero generator
    for col in range(ker_in_A3.shape[1]):
        v = ker_in_A3[:, col]
        if np.linalg.norm(v) > 1e-8:
            # Normalize
            v = v / np.max(np.abs(v))
            # Express in terms of allowed 3-paths
            gen = {}
            for i, path in enumerate(allowed[3]):
                if abs(v[i]) > 1e-6:
                    gen[path] = round(v[i], 4)
            return gen, allowed[3], beta3

    return None, None, None

# Search for β₃ = 1 tournaments at n=7
print("=" * 70)
print("SEARCHING FOR β₃ = 1 TOURNAMENTS AT n=7")
print("=" * 70)

n = 7
E = n*(n-1)//2
random.seed(42)
found = 0

for trial in range(2000):
    bits = random.randint(0, (1 << E) - 1)
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=4)

    if betti[3] > 0:
        found += 1
        c3 = count_3cycles(A, n)
        sc = score_seq(A, n)
        print(f"\n--- β₃={betti[3]} tournament #{found} ---")
        print(f"  bits={bits}, c₃={c3}, score={sc}")
        print(f"  Full β = {betti}")

        # Find generator
        gen, paths, b3 = find_h3_generator(A, n)
        if gen:
            print(f"  H₃ generator ({len(gen)} terms):")
            # Group by vertex support
            vertex_support = set()
            for path, coeff in sorted(gen.items(), key=lambda x: -abs(x[1])):
                vertex_support.update(path)
                if abs(coeff) > 0.01:
                    print(f"    {coeff:+.3f} * {path}")

            print(f"  Vertex support: {sorted(vertex_support)} ({len(vertex_support)} vertices)")
            print(f"  Number of terms: {len(gen)}")

            # Check: is the generator "localized"?
            supports = Counter(len(set(p)) for p in gen.keys())
            print(f"  Path lengths in generator: {dict(supports)}")

        if found >= 3:
            break

# Check mutual exclusivity with β₁
print("\n" + "=" * 70)
print("MUTUAL EXCLUSIVITY CHECK")
print("=" * 70)

random.seed(123)
b1_and_b3 = 0
b1_only = 0
b3_only = 0
neither = 0

for trial in range(3000):
    bits = random.randint(0, (1 << E) - 1)
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=4)

    b1 = betti[1] if len(betti) > 1 else 0
    b3 = betti[3] if len(betti) > 3 else 0

    if b1 > 0 and b3 > 0:
        b1_and_b3 += 1
        print(f"  *** BOTH β₁={b1} AND β₃={b3}: bits={bits}, β={betti}")
    elif b1 > 0:
        b1_only += 1
    elif b3 > 0:
        b3_only += 1
    else:
        neither += 1

print(f"\nn=7, 3000 samples:")
print(f"  β₁>0 only: {b1_only}")
print(f"  β₃>0 only: {b3_only}")
print(f"  Both >0: {b1_and_b3}")
print(f"  Neither: {neither}")
if b1_and_b3 == 0:
    print(f"  ✓ β₁ and β₃ are MUTUALLY EXCLUSIVE")
