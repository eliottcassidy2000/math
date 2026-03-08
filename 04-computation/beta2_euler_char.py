#!/usr/bin/env python3
"""
beta2_euler_char.py — Euler characteristic approach to β₂=0

The Euler characteristic χ = Σ (-1)^p β_p is also equal to
Σ (-1)^p dim(Ω_p).

If we can compute χ independently and show it equals 1 - β₁,
then β₂ = 0 follows from β₂ = β₁ - χ + 1 + (higher terms).

More precisely:
χ = β₀ - β₁ + β₂ - β₃ + ...
  = 1 - β₁ + β₂ - β₃ + ...

Also: χ = dim(Ω₀) - dim(Ω₁) + dim(Ω₂) - dim(Ω₃) + ...

If we can show dim(Ω_0) - dim(Ω_1) + dim(Ω_2) - dim(Ω_3) + ... = 1 - β₁,
then β₂ - β₃ + β₄ - ... = 0, which combined with β₂ ≥ 0 and
non-negativity constraints might give β₂ = 0.

But this is weak — we'd need β₃ ≥ β₂ which isn't obvious.

BETTER: Use the ALTERNATING SUM FORMULA for ranks.
rk(∂₁) - rk(∂₂) + rk(∂₃) - ... = dim(Ω₀) - χ

This gives exact relationships between ranks.

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

print("="*70)
print("EULER CHARACTERISTIC AND RANK FORMULAS")
print("="*70)

for n in [4, 5]:
    print(f"\n--- n = {n} ---")

    chi_distribution = Counter()
    rank_patterns = Counter()

    for A in all_tournaments(n):
        allowed = {}
        for p in range(n):
            allowed[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

        omega = {}
        for p in range(n):
            if p == 0:
                omega[p] = np.eye(n)
            elif allowed[p]:
                omega[p] = compute_omega_basis(A, n, p, allowed[p],
                                                allowed[p-1] if p >= 1 else [])
            else:
                omega[p] = np.zeros((0, 0))

        dims = []
        ranks = []
        betti = []

        for p in range(n):
            dim_p = omega[p].shape[1] if omega[p].ndim == 2 and omega[p].shape[0] > 0 else (n if p == 0 else 0)
            dims.append(dim_p)

            if p == 0:
                rk = 0  # ∂₀ = 0
            elif dim_p == 0:
                rk = 0
            else:
                bd = build_full_boundary_matrix(allowed[p], allowed[p-1])
                bd_om = bd @ omega[p]
                S = np.linalg.svd(bd_om, compute_uv=False)
                rk = sum(s > 1e-8 for s in S)
            ranks.append(rk)

        # Betti numbers
        for p in range(n):
            ker = dims[p] - ranks[p]
            im_next = ranks[p+1] if p + 1 < n else 0
            betti.append(ker - im_next)

        chi = sum((-1)**p * b for p, b in enumerate(betti))
        chi_alt = sum((-1)**p * d for p, d in enumerate(dims))

        assert abs(chi - chi_alt) < 1e-8, f"χ mismatch: {chi} vs {chi_alt}"

        chi_distribution[int(round(chi))] += 1
        rank_patterns[tuple(ranks)] += 1

    print(f"  χ distribution: {dict(chi_distribution)}")
    print(f"\n  Rank patterns:")
    for pattern, count in sorted(rank_patterns.items()):
        print(f"    rk = {pattern}: {count}")

    # Check: does χ = 1 - β₁ always?
    print(f"\n  Verifying χ = 1 - β₁:")
    verified = True
    for A in all_tournaments(n):
        allowed = {}
        for p in range(3):
            allowed[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

        om1 = compute_omega_basis(A, n, 1, allowed[1], allowed[0])
        d1 = om1.shape[1] if om1.ndim == 2 else 0

        bd1 = build_full_boundary_matrix(allowed[1], allowed[0])
        bd1_om = bd1 @ om1
        S = np.linalg.svd(bd1_om, compute_uv=False)
        rk1 = sum(s > 1e-8 for s in S)
        beta1 = d1 - rk1  # ker(∂₁) = dim(Ω₁) - rk(∂₁), im(∂₂) contributes later

        # Actually β₁ = ker(∂₁) - im(∂₂)
        om2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
        d2 = om2.shape[1] if om2.ndim == 2 else 0
        if d2 > 0:
            bd2 = build_full_boundary_matrix(allowed[2], allowed[1])
            bd2_om = bd2 @ om2
            coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
            S2 = np.linalg.svd(coords, compute_uv=False)
            im2 = sum(s > 1e-8 for s in S2)
        else:
            im2 = 0

        beta1_actual = (d1 - rk1) - im2

        # Compute χ via dims
        dims_full = []
        for p in range(n):
            ap = enumerate_allowed_paths(A, n, p)
            if p == 0:
                dims_full.append(n)
            elif ap:
                op = compute_omega_basis(A, n, p, ap, enumerate_allowed_paths(A, n, p-1))
                dims_full.append(op.shape[1] if op.ndim == 2 else 0)
            else:
                dims_full.append(0)

        chi = sum((-1)**p * d for p, d in enumerate(dims_full))

        if abs(chi - (1 - beta1_actual)) > 0.5:
            print(f"    FAIL: χ={chi}, 1-β₁={1-beta1_actual}")
            verified = False
            break

    if verified:
        print(f"  ✓ χ = 1 - β₁ for all n={n} tournaments")
        print(f"  This means: β₂ - β₃ + β₄ - ... = 0")
        print(f"  Combined with β₂=0 (verified), this gives β₃ - β₄ + ... = 0")

# KEY INSIGHT
print(f"\n{'='*70}")
print("KEY INSIGHT: χ = 1 - β₁ IS NOT ENOUGH FOR β₂ = 0")
print("="*70)
print(f"""
χ = 1 - β₁ + β₂ - β₃ + β₄ - ...

If χ = 1 - β₁, then β₂ - β₃ + β₄ - ... = 0.

This means β₂ = β₃ - β₄ + β₅ - ...

At n=5: β₃ = β₄ = 0, so β₂ = 0. ✓
At n=7: β₃ can be 1, β₄ = 0 for the same tournament.
  So β₂ = β₃ - 0 = β₃ ≥ 0. But β₂ = β₃ only if β₂ ≥ 0.
  Wait: χ = 1 - β₁ gives β₂ = β₃ always (when β₄=β₅=...=0).
  Since β₂ ≥ 0 and β₃ ≥ 0, this is consistent but NOT a proof of β₂=0.

  For β₃=1 (S-phase at n=7): β₂ = β₃ = 1? NO! β₂ = 0 computationally.
  So χ ≠ 1 - β₁ when β₃ > 0? Let me check.
""")

# Verify χ formula at n=7 for some S-phase tournaments
print(f"Checking χ at n=7 for S-phase tournaments...")
n = 7
import random
random.seed(42)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

for _ in range(20):
    mask = random.randint(0, (1 << m) - 1)
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    dims_full = []
    for p in range(n):
        ap = enumerate_allowed_paths(A, n, p)
        if p == 0:
            dims_full.append(n)
        elif ap:
            op = compute_omega_basis(A, n, p, ap, enumerate_allowed_paths(A, n, p-1))
            dims_full.append(op.shape[1] if op.ndim == 2 else 0)
        else:
            dims_full.append(0)

    chi = sum((-1)**p * d for p, d in enumerate(dims_full))

    # Get β₁
    a1 = enumerate_allowed_paths(A, n, 1)
    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))
    bd1 = build_full_boundary_matrix(a1, enumerate_allowed_paths(A, n, 0))
    bd1_om = bd1 @ om1
    S = np.linalg.svd(bd1_om, compute_uv=False)
    rk1 = sum(s > 1e-8 for s in S)

    a2 = enumerate_allowed_paths(A, n, 2)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    if om2.ndim == 2 and om2.shape[1] > 0:
        bd2 = build_full_boundary_matrix(a2, a1)
        bd2_om = bd2 @ om2
        coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
        S2 = np.linalg.svd(coords, compute_uv=False)
        im2 = sum(s > 1e-8 for s in S2)
    else:
        im2 = 0

    beta1 = (om1.shape[1] - rk1) - im2

    if chi != 1 - beta1:
        print(f"  χ={chi}, 1-β₁={1-beta1}, diff={chi-(1-beta1)}")

print(f"\n  If χ ≠ 1-β₁ at n=7, it means higher Betti numbers contribute.")
print(f"  The Euler char approach alone can't prove β₂=0.")

print("\nDONE")
