#!/usr/bin/env python3
"""
beta2_omega3_filling.py - Does the position-cone filling land in Omega_3?

The position cone test showed that for every z in ker(d_2|Omega_2),
there exist weights such that z = d_3(w) for some w in A_3.

BUT we need w in Omega_3, not just A_3!

This script checks:
1. Whether the position-cone image happens to lie in Omega_3
2. Whether ker(d_2|Omega_2) = im(d_3|Omega_3) directly
3. The relationship between im(d_3|A_3) and im(d_3|Omega_3)

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def direct_beta2_check(A, n):
    """Directly verify ker(d_2|Omega_2) = im(d_3|Omega_3).

    Returns True if beta_2 = 0, plus diagnostic data.
    """
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)

    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    if dim_O2 == 0:
        return True, {'ker_d2': 0, 'im_d3': 0}

    # d_2 restricted to Omega_2
    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    S2 = np.linalg.svd(D2_omega, compute_uv=False)
    rank_d2 = sum(s > 1e-8 for s in S2)
    ker_d2 = dim_O2 - rank_d2

    # d_3 restricted to Omega_3 -> image in A_2
    if dim_O3 > 0:
        D3 = build_full_boundary_matrix(paths3, paths2)
        D3_omega = D3 @ omega3  # columns = image vectors in A_2
        S3 = np.linalg.svd(D3_omega, compute_uv=False)
        im_d3 = sum(s > 1e-8 for s in S3)
    else:
        im_d3 = 0

    # Also compute im(d_3|A_3) for comparison
    D3_full = build_full_boundary_matrix(paths3, paths2)
    S3_full = np.linalg.svd(D3_full, compute_uv=False) if len(paths3) > 0 else []
    im_d3_full = sum(s > 1e-8 for s in S3_full) if len(S3_full) > 0 else 0

    beta2 = ker_d2 - im_d3

    return beta2 == 0, {
        'dim_O2': dim_O2, 'dim_O3': dim_O3,
        'dim_A2': len(paths2), 'dim_A3': len(paths3),
        'rank_d2': rank_d2, 'ker_d2': ker_d2,
        'im_d3_omega': im_d3, 'im_d3_full': im_d3_full,
        'beta2': beta2,
    }


# ============================================================
# PART 1: Verify im(d_3|Omega_3) = im(d_3|A_3) for tournaments?
# ============================================================
print("=" * 70)
print("im(d_3|Omega_3) vs im(d_3|A_3)")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m
    all_equal = True
    mismatch = 0

    for bits in range(total):
        A = build_adj(n, bits)
        ok, info = direct_beta2_check(A, n)
        if info['im_d3_omega'] != info['im_d3_full']:
            all_equal = False
            mismatch += 1
            if mismatch <= 5:
                scores = tuple(sorted([sum(row) for row in A]))
                print(f"  MISMATCH n={n} T#{bits} scores={scores}: "
                      f"im(d3|O3)={info['im_d3_omega']}, im(d3|A3)={info['im_d3_full']}")
                print(f"    dim(O3)={info['dim_O3']}, dim(A3)={info['dim_A3']}")

    print(f"\nn={n}: im(d3|Omega_3) == im(d3|A_3) for ALL? {all_equal}")
    if not all_equal:
        print(f"  Mismatches: {mismatch}/{total}")


# ============================================================
# PART 2: Exhaustive n=6 verification of beta_2 = 0
# with im(d_3) comparison
# ============================================================
print(f"\n{'='*70}")
print("n=6: EXHAUSTIVE beta_2 CHECK + im comparison")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

all_beta2_zero = True
im_equal_count = 0
im_differ_count = 0

for bits in range(total):
    A = build_adj(n, bits)
    ok, info = direct_beta2_check(A, n)

    if not ok:
        all_beta2_zero = False
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  beta_2 != 0: T#{bits} scores={scores} beta_2={info['beta2']}")

    if info['im_d3_omega'] == info['im_d3_full']:
        im_equal_count += 1
    else:
        im_differ_count += 1

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s) beta2_zero={all_beta2_zero} im_equal={im_equal_count} im_differ={im_differ_count}")

elapsed = time.time() - t0
print(f"\nn=6: ({elapsed:.0f}s)")
print(f"  beta_2 = 0 for ALL? {all_beta2_zero}")
print(f"  im(d3|O3) == im(d3|A3): {im_equal_count}/{total}")
print(f"  im(d3|O3) != im(d3|A3): {im_differ_count}/{total}")


# ============================================================
# PART 3: Key question - is im(d_3|A_3) cap Omega_2 = im(d_3|Omega_3)?
# i.e., can we "correct" any A_3 preimage to an Omega_3 preimage?
# ============================================================
print(f"\n{'='*70}")
print("CORRECTION: A_3 preimage -> Omega_3 preimage")
print("=" * 70)

# For each z in ker(d_2|Omega_2), we have z = d_3(w_A) for w_A in A_3
# (by the position cone test).
# Question: can we find w_O in Omega_3 with d_3(w_O) = z?
#
# Since d_3(w_A) = z and d_3(w_O) = z, we need w_O = w_A + ker(d_3).
# So: does (w_A + ker(d_3|A_3)) intersect Omega_3?
# i.e., does w_A + ker_d3 contain an Omega_3 element?
#
# Equivalently: is the projection of w_A onto A_3 / Omega_3 in the
# image of ker(d_3) projected to A_3 / Omega_3?
#
# This is automatically true if Omega_3 + ker(d_3|A_3) = A_3
# i.e., if every element of A_3 can be written as omega + kernel.

n = 5
m = n*(n-1)//2
total = 1 << m

for bits in range(total):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
    dim_A3 = len(paths3)

    # ker(d_3|A_3)
    D3 = build_full_boundary_matrix(paths3, paths2)
    if dim_A3 > 0 and len(paths2) > 0:
        U3, S3, V3t = np.linalg.svd(D3, full_matrices=True)
        rank_d3_A3 = sum(s > 1e-8 for s in S3)
        ker_d3_A3_basis = V3t[rank_d3_A3:].T  # columns = ker vectors
        dim_ker_d3 = ker_d3_A3_basis.shape[1] if ker_d3_A3_basis.ndim == 2 else 0
    else:
        dim_ker_d3 = dim_A3

    # Check: Omega_3 + ker(d_3|A_3) = A_3?
    # Stack columns of omega3 and ker_d3_A3_basis
    if dim_O3 > 0 and dim_ker_d3 > 0:
        combined = np.hstack([omega3, ker_d3_A3_basis])
        Sc = np.linalg.svd(combined, compute_uv=False)
        rank_combined = sum(s > 1e-8 for s in Sc)
        spans_all = (rank_combined == dim_A3)
    elif dim_O3 > 0:
        spans_all = (dim_O3 == dim_A3)
    elif dim_ker_d3 > 0:
        spans_all = (dim_ker_d3 == dim_A3)
    else:
        spans_all = (dim_A3 == 0)

    if not spans_all and bits < 20:
        scores = tuple(sorted([sum(row) for row in A]))
        print(f"  T#{bits} scores={scores}: Omega_3+ker(d3) rank={rank_combined}, dim(A3)={dim_A3}")
        print(f"    dim(Omega_3)={dim_O3}, dim(ker_d3)={dim_ker_d3}")

# Overall check
all_span = True
for bits in range(total):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
    dim_A3 = len(paths3)
    D3 = build_full_boundary_matrix(paths3, paths2)
    if dim_A3 > 0 and len(paths2) > 0:
        U3, S3, V3t = np.linalg.svd(D3, full_matrices=True)
        rank_d3 = sum(s > 1e-8 for s in S3)
        ker_d3_basis = V3t[rank_d3:].T
        dim_ker = ker_d3_basis.shape[1] if ker_d3_basis.ndim == 2 else 0
    else:
        dim_ker = dim_A3
        ker_d3_basis = np.eye(dim_A3) if dim_A3 > 0 else np.zeros((0,0))

    if dim_O3 > 0 and dim_ker > 0:
        combined = np.hstack([omega3, ker_d3_basis])
        Sc = np.linalg.svd(combined, compute_uv=False)
        rank_comb = sum(s > 1e-8 for s in Sc)
        if rank_comb < dim_A3:
            all_span = False
            break
    elif dim_O3 + dim_ker < dim_A3:
        all_span = False
        break

print(f"\nn=5: Omega_3 + ker(d_3|A_3) = A_3 for ALL? {all_span}")


# ============================================================
# PART 4: Simplest check — does d_3|A_3 already have ker(d_2|O_2) in image?
# ============================================================
print(f"\n{'='*70}")
print("DOES d_3(A_3) contain ker(d_2|O_2)?")
print("=" * 70)

n = 5
works = 0
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        works += 1
        continue

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    U2, S2, V2t = np.linalg.svd(D2_omega, full_matrices=True)
    rank_d2 = sum(s > 1e-8 for s in S2)
    ker_d2_basis = V2t[rank_d2:]  # rows = ker vectors in O2 coords
    ker_in_A2 = omega2 @ ker_d2_basis.T  # columns in A_2

    # Image of d_3 from A_3
    D3 = build_full_boundary_matrix(paths3, paths2)

    # Check: is each ker_d2 vector in im(D3)?
    # Stack D3 columns and ker vectors, check rank
    if D3.shape[1] > 0:
        combined = np.hstack([D3, ker_in_A2])
        Sc = np.linalg.svd(combined, compute_uv=False)
        rank_comb = sum(s > 1e-8 for s in Sc)
        Sd = np.linalg.svd(D3, compute_uv=False)
        rank_D3 = sum(s > 1e-8 for s in Sd)
        in_image = (rank_comb == rank_D3)
    else:
        in_image = (ker_in_A2.shape[1] == 0)

    if in_image:
        works += 1

print(f"\nn=5: ker(d_2|O_2) in im(d_3|A_3) for ALL? {works == 1024}")
print(f"  {works}/1024 succeed")

# Same for n=4
n = 4
works4 = 0
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        works4 += 1
        continue

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    U2, S2, V2t = np.linalg.svd(D2_omega, full_matrices=True)
    rank_d2 = sum(s > 1e-8 for s in S2)
    ker_d2_basis = V2t[rank_d2:]
    ker_in_A2 = omega2 @ ker_d2_basis.T

    D3 = build_full_boundary_matrix(paths3, paths2)
    if D3.shape[1] > 0:
        combined = np.hstack([D3, ker_in_A2])
        Sc = np.linalg.svd(combined, compute_uv=False)
        rank_comb = sum(s > 1e-8 for s in Sc)
        Sd = np.linalg.svd(D3, compute_uv=False)
        rank_D3 = sum(s > 1e-8 for s in Sd)
        in_image = (rank_comb == rank_D3)
    else:
        in_image = (ker_in_A2.shape[1] == 0)

    if in_image:
        works4 += 1

print(f"n=4: {works4}/64 succeed")


print("\n\nDone.")
