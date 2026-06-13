#!/usr/bin/env python3
"""
beta2_dt_deficit_n7.py - DT filling analysis at n=7

At n=6: DT alone fills Z_2 for 97.1%, deficit is always exactly 1.
What happens at n=7? Is the deficit always bounded? Always small?

Sample n=7 tournaments and measure:
  - dim(Z_2), dim(im(d_3|DT)), gap = dim(Z_2) - rk(DT_bd)
  - dim(Omega_3), |DT|
  - Full Omega_3 fills Z_2? (confirming beta_2 = 0)

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
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

def random_tournament(n):
    bits = random.getrandbits(n*(n-1)//2)
    return bits

def find_DT_paths(A, n):
    """Find all DT 4-paths: (a,b,c,d) with a->b->c->d, a->c, b->d."""
    dt = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if not A[a][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c or not A[c][d]: continue
                    if A[b][d]:
                        dt.append((a,b,c,d))
    return dt

def compute_DT_boundary_matrix(dt_paths, allowed_2, n):
    if not dt_paths or not allowed_2:
        return np.zeros((len(allowed_2), 0))
    path_to_idx = {p: i for i, p in enumerate(allowed_2)}
    mat = np.zeros((len(allowed_2), len(dt_paths)))
    for j, (a, b, c, d) in enumerate(dt_paths):
        faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
        signs = [1, -1, 1, -1]
        for face, sign in zip(faces, signs):
            if face in path_to_idx:
                mat[path_to_idx[face], j] += sign
    return mat

n = 7
n_arcs = n*(n-1)//2
N_SAMPLES = 5000

print("=" * 70)
print(f"DT DEFICIT ANALYSIS AT n={n} ({N_SAMPLES} samples)")
print("=" * 70)

gap_counts = {}  # gap -> count
beta2_violations = 0
trivial = 0

t0 = time.time()
for trial in range(N_SAMPLES):
    if trial % 1000 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES} ({dt:.0f}s)")

    bits = random_tournament(n)
    A = build_adj(n, bits)

    allowed_3 = enumerate_allowed_paths(A, n, 3)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    omega2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        trivial += 1
        continue

    # Z_2 = ker(d_2|Omega_2)
    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    bd2_om = bd2 @ omega2
    S_v = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S_v) > 1e-8))
    Z2_dim = dim_O2 - rk2

    if Z2_dim == 0:
        trivial += 1
        continue

    # DT filling
    dt_paths = find_DT_paths(A, n)
    DT_bd = compute_DT_boundary_matrix(dt_paths, allowed_2, n)
    rk_DT = np.linalg.matrix_rank(DT_bd, tol=1e-8)

    # Get Z_2 basis
    U_dt, S_dt, Vt_dt = np.linalg.svd(bd2_om, full_matrices=True)
    Z2_omega = Vt_dt[rk2:].T
    Z2_A2 = omega2 @ Z2_omega

    combined = np.hstack([DT_bd, Z2_A2])
    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)

    gap = rk_combined - rk_DT  # dimensions of Z_2 not in DT image

    gap_counts[gap] = gap_counts.get(gap, 0) + 1

    # Also check full Omega_3 fills Z_2
    if gap > 0:
        omega3 = compute_omega_basis(A, n, 3, allowed_3, allowed_2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

        if dim_O3 > 0:
            bd3 = build_full_boundary_matrix(allowed_3, allowed_2)
            bd3_om = bd3 @ omega3
            full_combined = np.hstack([bd3_om, Z2_A2])
            rk_full = np.linalg.matrix_rank(bd3_om, tol=1e-8)
            rk_full_combined = np.linalg.matrix_rank(full_combined, tol=1e-8)

            if rk_full_combined > rk_full:
                beta2_violations += 1
                scores = tuple(sorted(sum(row) for row in A))
                print(f"  *** BETA_2 > 0 at bits={bits}: scores={scores}, "
                      f"Z2={Z2_dim}, rk(d3|O3)={rk_full}, gap={gap}")

dt = time.time() - t0
print(f"\n  n={n}, {N_SAMPLES} samples ({dt:.0f}s)")
print(f"  Trivial (Z_2=0): {trivial}")
print(f"  beta_2 violations: {beta2_violations}")
print(f"\n  DT gap distribution:")
for gap in sorted(gap_counts.keys()):
    pct = 100*gap_counts[gap]/(N_SAMPLES - trivial)
    print(f"    gap={gap}: {gap_counts[gap]} ({pct:.1f}%)")

# Score analysis for deficit tournaments
print(f"\n  Now analyzing score patterns of high-deficit cases...")

# Rerun a few high-deficit cases for detailed analysis
high_gap_examples = []
random.seed(42)
for trial in range(2000):
    bits = random_tournament(n)
    A = build_adj(n, bits)

    allowed_3 = enumerate_allowed_paths(A, n, 3)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)

    omega2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    if dim_O2 == 0: continue

    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    bd2_om = bd2 @ omega2
    S_v = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S_v) > 1e-8))
    Z2_dim = dim_O2 - rk2
    if Z2_dim == 0: continue

    dt_paths = find_DT_paths(A, n)
    DT_bd = compute_DT_boundary_matrix(dt_paths, allowed_2, n)
    rk_DT = np.linalg.matrix_rank(DT_bd, tol=1e-8)

    U_dt, S_dt, Vt_dt = np.linalg.svd(bd2_om, full_matrices=True)
    Z2_omega = Vt_dt[rk2:].T
    Z2_A2 = omega2 @ Z2_omega

    combined = np.hstack([DT_bd, Z2_A2])
    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
    gap = rk_combined - rk_DT

    if gap >= 2:
        scores = tuple(sorted(sum(row) for row in A))
        omega3 = compute_omega_basis(A, n, 3, allowed_3, allowed_2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
        high_gap_examples.append((bits, scores, Z2_dim, rk_DT, gap, len(dt_paths), dim_O3))

if high_gap_examples:
    print(f"\n  Found {len(high_gap_examples)} cases with gap >= 2:")
    for bits, scores, Z2, rkDT, gap, ndt, dO3 in high_gap_examples[:10]:
        print(f"    bits={bits}: scores={scores}, Z2={Z2}, rk(DT)={rkDT}, "
              f"gap={gap}, |DT|={ndt}, dim(O3)={dO3}")
else:
    print(f"\n  No cases with gap >= 2 found in 2000 samples")

print("\nDone.")
