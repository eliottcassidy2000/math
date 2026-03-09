#!/usr/bin/env python3
"""
beta2_z2_in_imd3.py — Is Z₂(Ω) ⊆ im(∂₃|A₃)?

If YES, then combined with im(∂₃|A₃) ∩ Ω₂ = im(∂₃|Ω₃) (verified),
we get Z₂(Ω) ⊆ B₂(Ω), hence β₂ = 0.

Z₂(Ω) ⊆ im(∂₃|A₃) means: for every z ∈ Z₂(Ω), ∃σ ∈ A₃ with ∂₃(σ) = z.

Note: A_* is NOT a chain complex. ∂₃(A₃) is not contained in ker(∂₂).
But Z₂(Ω) IS in ker(∂₂|Ω₂). The question is whether Z₂(Ω) ⊆ ∂₃(A₃).

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


print("=" * 70)
print("Z₂(Ω) ⊆ im(∂₃|A₃)? UNIVERSAL TEST")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    t0 = time.time()

    z2_in_im = 0
    z2_not_in_im = 0
    z2_zero = 0

    for bits in range(total):
        A = build_adj(n, bits)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

        d2 = dim_om(om2)
        if d2 == 0:
            z2_zero += 1
            continue

        # Z₂ in A₂ coordinates
        bd2 = build_full_boundary_matrix(ap2, ap1)
        d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
        rk = sum(s > 1e-8 for s in S)
        z2_dim = d2 - rk

        if z2_dim == 0:
            z2_zero += 1
            continue

        z2_A2 = om2 @ Vt[rk:, :].T  # |A₂| × z2_dim

        # ∂₃: A₃ → A₂ (FULL boundary, not restricted to Ω)
        bd3_full = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

        # Check: is each Z₂ basis vector in im(∂₃|A₃)?
        # im(∂₃|A₃) = column space of bd3_full
        # z ∈ im(∂₃|A₃) iff z is in the column space of bd3_full
        if bd3_full.shape[1] > 0:
            # For each Z₂ basis vector: check if it's in column space of bd3_full
            all_in = True
            for j in range(z2_dim):
                z = z2_A2[:, j]
                # Solve bd3_full @ σ = z
                sigma, resid, _, _ = np.linalg.lstsq(bd3_full, z, rcond=None)
                actual_resid = np.linalg.norm(bd3_full @ sigma - z)
                if actual_resid > 1e-6:
                    all_in = False
                    break

            if all_in:
                z2_in_im += 1
            else:
                z2_not_in_im += 1
                if z2_not_in_im <= 3:
                    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
                    print(f"  COUNTEREXAMPLE at n={n}: bits={bits}, scores={scores}")
                    print(f"    z2_dim={z2_dim}, |A₃|={len(ap3)}, resid={actual_resid:.6f}")
        else:
            # No 3-paths at all: Z₂ can't be in im(∂₃|A₃) unless Z₂=0
            z2_not_in_im += 1

        if (bits+1) % 10000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {bits+1}/{total} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")
    print(f"  Z₂=0: {z2_zero}")
    print(f"  Z₂ ⊆ im(∂₃|A₃): {z2_in_im}")
    print(f"  Z₂ ⊄ im(∂₃|A₃): {z2_not_in_im}")

    if z2_not_in_im == 0:
        print(f"  *** Z₂(Ω) ⊆ im(∂₃|A₃) CONFIRMED for all n={n} tournaments! ***")

print("\nDone.")
