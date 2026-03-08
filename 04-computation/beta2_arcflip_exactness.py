#!/usr/bin/env python3
"""
beta2_arcflip_exactness.py - Track exactness at Omega_2 under arc flips

Key idea: The transitive tournament has beta_2 = 0 (it's a simplex).
Under each arc flip, dim(Z_2) and rk(d_3) change.
If we can find a formula for these changes, we can prove beta_2 stays 0.

Specifically, for arc flip u->v becoming v->u:
  delta_Z2 = dim(Z_2(T')) - dim(Z_2(T))
  delta_im3 = rk(d_3|Omega_3(T')) - rk(d_3|Omega_3(T))
  beta_2 = 0 iff delta_Z2 = delta_im3 at every step

This script computes these deltas exhaustively at n=5 and analyses their structure.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
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

def compute_exactness_data(A, n):
    """Compute dim(Omega_2), dim(Z_2), rk(d_3), beta_2."""
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
        return {'dim_om2': 0, 'dim_om3': 0, 'dim_z2': 0, 'rk_d3': 0, 'beta2': 0}

    # d_2: Omega_2 -> Omega_1
    bd2 = build_full_boundary_matrix(paths[2], paths[1])
    bd2_om = bd2 @ omega[2]
    sv2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk_d2 = int(np.sum(np.abs(sv2) > 1e-8))
    dim_z2 = dim_om2 - rk_d2

    # d_3: Omega_3 -> Omega_2
    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(paths[3], paths[2])
        bd3_om = bd3 @ omega[3]
        # Project into Omega_2 coordinates
        im3_coords, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)
        rk_d3 = np.linalg.matrix_rank(im3_coords, tol=1e-8)
    else:
        rk_d3 = 0

    beta2 = dim_z2 - rk_d3

    return {
        'dim_om2': dim_om2, 'dim_om3': dim_om3,
        'dim_z2': dim_z2, 'rk_d3': rk_d3, 'beta2': beta2,
    }


def flip_arc(bits, i, j, n):
    """Flip the arc between vertices i and j (i<j).
    Returns new bits value."""
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits  # shouldn't reach here


# ===== Exhaustive n=5 =====
print("=" * 70)
print("ARC FLIP EXACTNESS TRACKING")
print("=" * 70)

for n in [5]:
    print(f"\n--- n={n} exhaustive ---")
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs

    # Precompute exactness data for all tournaments
    t0 = time.time()
    data = {}
    for bits in range(total):
        A = build_adj(n, bits)
        data[bits] = compute_exactness_data(A, n)

    dt = time.time() - t0
    print(f"  Precomputed {total} tournaments in {dt:.1f}s")

    # Verify beta_2 = 0 for all
    b2_nonzero = sum(1 for d in data.values() if d['beta2'] != 0)
    print(f"  beta_2 != 0: {b2_nonzero}")

    # For each tournament and each arc, compute the delta under flip
    delta_stats = Counter()
    delta_z2_vs_im3 = Counter()  # (delta_z2, delta_im3) pairs
    delta_by_score_diff = defaultdict(list)

    for bits in range(total):
        A = build_adj(n, bits)
        d = data[bits]
        scores = [sum(row) for row in A]

        for i in range(n):
            for j in range(i+1, n):
                bits2 = flip_arc(bits, i, j, n)
                d2 = data[bits2]

                delta_z2 = d2['dim_z2'] - d['dim_z2']
                delta_im3 = d2['rk_d3'] - d['rk_d3']
                delta_om2 = d2['dim_om2'] - d['dim_om2']
                delta_om3 = d2['dim_om3'] - d['dim_om3']

                # Score difference: out-degree of source minus out-degree of target
                if A[i][j] == 1:
                    src, tgt = i, j
                else:
                    src, tgt = j, i
                score_diff = scores[src] - scores[tgt]

                delta_z2_vs_im3[(delta_z2, delta_im3)] += 1
                delta_by_score_diff[score_diff].append(
                    (delta_om2, delta_om3, delta_z2, delta_im3))

    print(f"\n  Delta (dim_Z2, rk_d3) distribution under arc flip:")
    for (dz, di), count in sorted(delta_z2_vs_im3.items()):
        match = "OK" if dz == di else "MISMATCH"
        print(f"    (dZ2={dz:+d}, dim3={di:+d}): {count} flips {match}")

    # Verify: delta_Z2 = delta_im3 ALWAYS
    mismatches = sum(c for (dz, di), c in delta_z2_vs_im3.items() if dz != di)
    print(f"\n  Mismatches (dZ2 != dim3): {mismatches}")

    # Analyze by score difference
    print(f"\n  Deltas by score difference (d_src - d_tgt):")
    for sd in sorted(delta_by_score_diff.keys()):
        entries = delta_by_score_diff[sd]
        avg_dom2 = np.mean([e[0] for e in entries])
        avg_dom3 = np.mean([e[1] for e in entries])
        avg_dz2 = np.mean([e[2] for e in entries])
        avg_dim3 = np.mean([e[3] for e in entries])
        print(f"    score_diff={sd:+d}: N={len(entries)}, "
              f"avg(dOm2)={avg_dom2:+.2f}, avg(dOm3)={avg_dom3:+.2f}, "
              f"avg(dZ2)={avg_dz2:+.2f}, avg(dim3)={avg_dim3:+.2f}")

# ===== n=6 sampled =====
print(f"\n{'='*70}")
print(f"n=6 SAMPLED (1000 tournaments, all arcs)")
print("=" * 70)
import random
random.seed(42)
n = 6
n_arcs = n*(n-1)//2

# Sample tournaments
sample_bits = [random.randint(0, (1 << n_arcs) - 1) for _ in range(1000)]

t0 = time.time()
mismatches_n6 = 0
delta_dist_n6 = Counter()

for idx, bits in enumerate(sample_bits):
    if idx % 200 == 0 and idx > 0:
        print(f"  ... {idx}/1000 ({time.time()-t0:.0f}s)")
    A = build_adj(n, bits)
    d = compute_exactness_data(A, n)

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            A2 = build_adj(n, bits2)
            d2 = compute_exactness_data(A2, n)

            delta_z2 = d2['dim_z2'] - d['dim_z2']
            delta_im3 = d2['rk_d3'] - d['rk_d3']

            delta_dist_n6[(delta_z2, delta_im3)] += 1
            if delta_z2 != delta_im3:
                mismatches_n6 += 1

dt = time.time() - t0
print(f"\n  Processed 1000 tournaments x {n_arcs} arcs = {1000*n_arcs} flips in {dt:.0f}s")
print(f"  Mismatches (dZ2 != dim3): {mismatches_n6}")
print(f"\n  Delta distribution:")
for (dz, di), count in sorted(delta_dist_n6.items()):
    match = "OK" if dz == di else "FAIL"
    print(f"    (dZ2={dz:+d}, dim3={di:+d}): {count} {match}")

print("\nDone.")
