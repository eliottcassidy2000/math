#!/usr/bin/env python3
"""
beta2_arcflip_exhaustive_n6.py - EXHAUSTIVE verification of beta_2 arc-flip invariance at n=6

Precompute exactness data for all 32768 tournaments on 6 vertices,
then check delta(Z_2) = delta(rk d_3) for all 32768 * 15 = 491520 arc flips.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
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

def compute_z2_and_rk3(A, n):
    """Compute dim(Z_2) and rk(d_3)."""
    paths = {}
    omega = {}
    for p in range(5):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dim_om2 = omega[2].shape[1] if omega[2].ndim == 2 else 0
    dim_om3 = omega[3].shape[1] if omega[3].ndim == 2 else 0

    if dim_om2 == 0:
        return 0, 0, 0

    bd2 = build_full_boundary_matrix(paths[2], paths[1])
    bd2_om = bd2 @ omega[2]
    sv2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk_d2 = int(np.sum(np.abs(sv2) > 1e-8))
    dim_z2 = dim_om2 - rk_d2

    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(paths[3], paths[2])
        bd3_om = bd3 @ omega[3]
        im3_coords, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)
        rk_d3 = np.linalg.matrix_rank(im3_coords, tol=1e-8)
    else:
        rk_d3 = 0

    return dim_z2, rk_d3, dim_z2 - rk_d3


def flip_arc_bits(bits, i, j, n):
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits


n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"EXHAUSTIVE beta_2 ARC-FLIP INVARIANCE at n={n}")
print(f"Total tournaments: {total}")
print(f"Total arc flips to check: {total * n_arcs}")
print("=" * 70)

# Phase 1: Precompute
print(f"\nPhase 1: Precomputing all {total} tournaments...")
t0 = time.time()
data = {}
for bits in range(total):
    if bits % 5000 == 0 and bits > 0:
        dt = time.time() - t0
        eta = dt / bits * (total - bits)
        print(f"  {bits}/{total} ({dt:.0f}s, ETA {eta:.0f}s)")
    A = build_adj(n, bits)
    data[bits] = compute_z2_and_rk3(A, n)

dt = time.time() - t0
print(f"  Precomputed {total} tournaments in {dt:.0f}s")

# Phase 2: Check all flips
print(f"\nPhase 2: Checking all {total * n_arcs} arc flips...")
t0 = time.time()
mismatches = 0
total_flips = 0

for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        dt = time.time() - t0
        print(f"  {bits}/{total} ({dt:.0f}s), mismatches so far: {mismatches}")

    z2, rk3, b2 = data[bits]

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc_bits(bits, i, j, n)
            z2_2, rk3_2, b2_2 = data[bits2]

            if (z2_2 - z2) != (rk3_2 - rk3):
                mismatches += 1
                if mismatches <= 5:
                    print(f"    MISMATCH: bits={bits}, flip ({i},{j}): "
                          f"dZ2={z2_2-z2}, drk3={rk3_2-rk3}")
            total_flips += 1

dt = time.time() - t0
print(f"\n  Total flips checked: {total_flips}")
print(f"  Mismatches: {mismatches}")
print(f"  Time: {dt:.0f}s")

if mismatches == 0:
    print(f"\n  *** CONFIRMED: beta_2 is arc-flip invariant at n={n} (EXHAUSTIVE) ***")

    # Also verify all beta_2 = 0
    b2_nonzero = sum(1 for d in data.values() if d[2] != 0)
    print(f"  beta_2 != 0 count: {b2_nonzero}")

print("\nDone.")
