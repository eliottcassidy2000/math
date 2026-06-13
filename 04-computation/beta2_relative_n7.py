#!/usr/bin/env python3
"""
beta2_relative_n7.py - Verify H_2(T, T\v) = 0 at n=7 by sampling

This is the KEY to an inductive proof of beta_2 = 0:
  - Base: beta_2 = 0 for n <= 4 (direct)
  - Induction: H_2(T, T\v) = 0 => H_2(T\v)=0 implies H_2(T)=0

Verified exhaustively at n=5,6 (HYP-213). Now sample n=7.

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
    return random.getrandbits(n*(n-1)//2)

def compute_relative_h2(A, n, v):
    """Compute H_2(T, T\v) numerically.

    Uses the quotient complex: C_p^rel = Omega_p(T) / Omega_p(T \ {v}).
    We identify Omega_p(T \ {v}) with the subspace of Omega_p(T)
    supported on paths not containing v.

    H_2^rel = ker(d_2^rel) / im(d_3^rel).
    """
    # Get allowed paths and Omega bases
    a2 = enumerate_allowed_paths(A, n, 2)
    a1 = enumerate_allowed_paths(A, n, 1)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 and om2.shape[1] > 0 else 0
    if dim_om2 == 0:
        return 0

    om1 = np.eye(n)  # Omega_1 = A_1 for tournaments (all edges present)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 and om3.shape[1] > 0 else 0

    # Boundary maps in Omega coordinates
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2  # |A_1| x dim_om2

    # Indices of v-paths in A_2 and A_3
    v_idx_2 = [i for i in range(len(a2)) if v in a2[i]]
    v_idx_1 = [i for i in range(len(a1)) if v in a1[i]]

    # V_2 = Omega_2 ∩ <no-v paths>: kernel of restriction to v-paths
    R_v2 = om2[v_idx_2, :]  # |v_paths| x dim_om2
    Sv2 = np.linalg.svd(R_v2, compute_uv=False)
    rank_rv2 = int(np.sum(np.abs(Sv2) > 1e-8))
    dim_V2 = dim_om2 - rank_rv2  # dim(Omega_2 ∩ no-v)

    # dim(Omega_2^rel) = rank_rv2
    if rank_rv2 == 0:
        return 0

    # ker(d_2) in Omega_2 coordinates
    U2, S2, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk_d2 = int(np.sum(np.abs(S2) > 1e-8))
    ker2_dim = dim_om2 - rk_d2

    if ker2_dim == 0:
        return 0  # No cycles at all

    ker2_basis = Vt2[rk_d2:].T  # dim_om2 x ker2_dim

    # ker(d_2) ∩ V_2: intersection of ker(d_2) with no-v subspace
    # V_2 basis in Omega_2 coords
    Uv, Sv, Vtv = np.linalg.svd(R_v2, full_matrices=True)
    if rank_rv2 < dim_om2:
        V2_basis = Vtv[rank_rv2:].T  # dim_om2 x dim_V2
    else:
        V2_basis = np.zeros((dim_om2, 0))

    if V2_basis.shape[1] > 0 and ker2_dim > 0:
        combined = np.hstack([ker2_basis, V2_basis])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        intersect_dim = ker2_dim + dim_V2 - rk_combined
    else:
        intersect_dim = 0

    ker2_rel = ker2_dim - intersect_dim

    if ker2_rel == 0:
        return 0

    # im(d_3^rel) in quotient Omega_2 / V_2
    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(a3, a2)
        bd3_om = bd3 @ om3  # |A_2| x dim_om3

        # Express im(d_3) in Omega_2 coordinates
        im3_om2, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)

        # Rank of im3_om2 modulo V_2
        if V2_basis.shape[1] > 0:
            combined_im = np.hstack([im3_om2, V2_basis])
            rk_im_comb = np.linalg.matrix_rank(combined_im, tol=1e-8)
            im3_rel = rk_im_comb - dim_V2
        else:
            im3_rel = np.linalg.matrix_rank(im3_om2, tol=1e-8)
    else:
        im3_rel = 0

    h2_rel = max(0, ker2_rel - im3_rel)
    return h2_rel


# ===== n=7 sampling =====
n = 7
N_SAMPLES = 3000

print("=" * 70)
print(f"RELATIVE HOMOLOGY H_2(T, T\\v) AT n={n} ({N_SAMPLES} samples)")
print("=" * 70)

violations = 0  # tournaments with H_2^rel > 0 for some v
all_zero = 0    # tournaments with H_2^rel = 0 for all v
max_h2_rel = 0

t0 = time.time()
for trial in range(N_SAMPLES):
    if trial % 500 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES} ({dt:.0f}s), violations={violations}")

    bits = random_tournament(n)
    A = build_adj(n, bits)

    all_v_zero = True
    for v_test in range(n):
        h2_rel = compute_relative_h2(A, n, v_test)
        if h2_rel > 0:
            all_v_zero = False
            if h2_rel > max_h2_rel:
                max_h2_rel = h2_rel
            scores = tuple(sorted(sum(row) for row in A))
            print(f"  *** H_2^rel > 0: bits={bits}, v={v_test}, "
                  f"h2_rel={h2_rel}, scores={scores}")

    if all_v_zero:
        all_zero += 1
    else:
        violations += 1

dt = time.time() - t0
print(f"\n  n={n}, {N_SAMPLES} samples ({dt:.0f}s)")
print(f"  All v: H_2^rel=0: {all_zero}/{N_SAMPLES}")
print(f"  Some v: H_2^rel>0: {violations}/{N_SAMPLES}")
print(f"  Max H_2^rel: {max_h2_rel}")

# ===== n=8 quick check =====
n = 8
N_SAMPLES_8 = 500

print(f"\n{'='*70}")
print(f"RELATIVE HOMOLOGY AT n={n} ({N_SAMPLES_8} samples)")
print("=" * 70)

violations_8 = 0
t0 = time.time()
for trial in range(N_SAMPLES_8):
    if trial % 100 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES_8} ({dt:.0f}s)")

    bits = random_tournament(n)
    A = build_adj(n, bits)

    for v_test in range(n):
        h2_rel = compute_relative_h2(A, n, v_test)
        if h2_rel > 0:
            violations_8 += 1
            scores = tuple(sorted(sum(row) for row in A))
            print(f"  *** H_2^rel > 0: bits={bits}, v={v_test}, "
                  f"h2_rel={h2_rel}, scores={scores}")
            break

dt = time.time() - t0
print(f"\n  n={n}, {N_SAMPLES_8} samples ({dt:.0f}s)")
print(f"  Violations: {violations_8}/{N_SAMPLES_8}")

# ===== n=9 very quick check =====
n = 9
N_SAMPLES_9 = 100

print(f"\n{'='*70}")
print(f"RELATIVE HOMOLOGY AT n={n} ({N_SAMPLES_9} samples)")
print("=" * 70)

violations_9 = 0
t0 = time.time()
for trial in range(N_SAMPLES_9):
    if trial % 25 == 0 and trial > 0:
        dt = time.time() - t0
        print(f"  ... {trial}/{N_SAMPLES_9} ({dt:.0f}s)")

    bits = random_tournament(n)
    A = build_adj(n, bits)

    # Test just a few vertices to save time
    for v_test in range(min(4, n)):
        h2_rel = compute_relative_h2(A, n, v_test)
        if h2_rel > 0:
            violations_9 += 1
            scores = tuple(sorted(sum(row) for row in A))
            print(f"  *** H_2^rel > 0: bits={bits}, v={v_test}, "
                  f"h2_rel={h2_rel}, scores={scores}")
            break

dt = time.time() - t0
print(f"\n  n={n}, {N_SAMPLES_9} samples ({dt:.0f}s)")
print(f"  Violations: {violations_9}/{N_SAMPLES_9}")

print("\nDone.")
