#!/usr/bin/env python3
"""
beta2_connecting_map.py - Check injectivity of connecting map delta

From the long exact sequence (assuming H_2(T\\v) = 0 by induction):
  0 -> H_2(T) -> H_2(T, T\\v) ->^delta H_1(T\\v) -> H_1(T) -> ...

H_2(T) = ker(delta). So beta_2(T) = 0 iff delta is injective.
Equivalently: dim(H_2(T, T\\v)) = dim(ker(i_*: H_1(T\\v) -> H_1(T))).

This DOES NOT require H_2(T, T\\v) = 0! Just delta-injectivity.

Verify at n=4,5,6: delta is always injective (confirming beta_2 = 0).

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
    build_full_boundary_matrix, path_betti_numbers
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

def delete_vertex(A, n, v):
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B

def compute_relative_h2_correct(A, n, v):
    """CORRECT relative H_2 computation."""
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
        return 0

    bd2 = build_full_boundary_matrix(paths[2], paths[1])
    bd2_om = bd2 @ omega[2]

    v_arc_idx = [i for i in range(len(paths[1])) if v in paths[1][i]]

    M = bd2_om[v_arc_idx, :]
    sv_M = np.linalg.svd(M, compute_uv=False)
    rk_M = int(np.sum(np.abs(sv_M) > 1e-8))
    dim_preimage = dim_om2 - rk_M

    v_path_idx_2 = [i for i in range(len(paths[2])) if v in paths[2][i]]
    if v_path_idx_2:
        R_v2 = omega[2][v_path_idx_2, :]
        sv_rv2 = np.linalg.svd(R_v2, compute_uv=False)
        rk_rv2 = int(np.sum(np.abs(sv_rv2) > 1e-8))
        dim_V2 = dim_om2 - rk_rv2
    else:
        dim_V2 = dim_om2

    ker2_rel = dim_preimage - dim_V2

    if ker2_rel <= 0:
        return 0

    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(paths[3], paths[2])
        bd3_om = bd3 @ omega[3]
        im3_om2, _, _, _ = np.linalg.lstsq(omega[2], bd3_om, rcond=None)

        if v_path_idx_2 and rk_rv2 < dim_om2:
            Uv, Sv, Vtv = np.linalg.svd(omega[2][v_path_idx_2, :], full_matrices=True)
            V2_basis = Vtv[rk_rv2:].T
        elif not v_path_idx_2:
            V2_basis = np.eye(dim_om2)
        else:
            V2_basis = np.zeros((dim_om2, 0))

        if V2_basis.shape[1] > 0:
            combined = np.hstack([im3_om2, V2_basis])
            rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
            im3_rel = rk_combined - V2_basis.shape[1]
        else:
            im3_rel = np.linalg.matrix_rank(im3_om2, tol=1e-8)
    else:
        im3_rel = 0

    return max(0, ker2_rel - im3_rel)


def compute_h1_kernel(A, n, v):
    """Compute dim(ker(i_*: H_1(T\\v) -> H_1(T))).

    i_* sends a 1-cycle of T\\v to the same cycle viewed in T.
    ker(i_*) = 1-cycles of T\\v that are 1-boundaries in T.
    """
    # H_1(T\\v)
    B = delete_vertex(A, n, v)
    betti_sub = path_betti_numbers(B, n-1, max_dim=2)
    h1_sub = betti_sub[1]

    if h1_sub == 0:
        return 0

    # Get Z_1(T\\v) and B_1(T) to find ker(i_*)
    # Z_1(T\\v): 1-cycles of T\\v
    paths_1_sub = enumerate_allowed_paths(B, n-1, 1)
    paths_0_sub = enumerate_allowed_paths(B, n-1, 0)
    paths_2_sub = enumerate_allowed_paths(B, n-1, 2)

    omega_1_sub = np.eye(n-1) if len(paths_1_sub) > 0 else np.zeros((n-1, 0))
    bd1_sub = build_full_boundary_matrix(paths_1_sub, paths_0_sub)

    U1s, S1s, Vt1s = np.linalg.svd(bd1_sub, full_matrices=True)
    rk1_sub = int(np.sum(np.abs(S1s) > 1e-8))
    z1_sub_dim = len(paths_1_sub) - rk1_sub  # dim Z_1(T\\v)

    if z1_sub_dim == 0:
        return 0

    z1_sub_basis = Vt1s[rk1_sub:].T  # |A_1(T\\v)| x z1_sub_dim

    # B_1(T\\v)
    omega_2_sub = compute_omega_basis(B, n-1, 2, paths_2_sub, paths_1_sub)
    dim_om2_sub = omega_2_sub.shape[1] if omega_2_sub.ndim == 2 else 0

    if dim_om2_sub > 0:
        bd2_sub = build_full_boundary_matrix(paths_2_sub, paths_1_sub)
        bd2_om_sub = bd2_sub @ omega_2_sub
        rk_bd2_sub = np.linalg.matrix_rank(bd2_om_sub, tol=1e-8)
    else:
        rk_bd2_sub = 0

    h1_sub_check = z1_sub_dim - rk_bd2_sub  # should equal betti_sub[1]

    # Now embed Z_1(T\\v) into A_1(T) and check which are in B_1(T)
    # Map: relabel vertices (remove v from numbering)
    # Vertices of T\\v are {0,...,n-1}\\{v}, relabeled as 0,...,n-2
    relabel = []
    for i in range(n):
        if i != v:
            relabel.append(i)
    # relabel[k] = original vertex for T\\v vertex k

    paths_1_T = enumerate_allowed_paths(A, n, 1)
    paths_2_T = enumerate_allowed_paths(A, n, 2)
    omega_2_T = compute_omega_basis(A, n, 2, paths_2_T, paths_1_T)
    dim_om2_T = omega_2_T.shape[1] if omega_2_T.ndim == 2 else 0

    # B_1(T): image of d_2 in A_1(T)
    bd1_T = build_full_boundary_matrix(paths_1_T, enumerate_allowed_paths(A, n, 0))
    # Z_1(T)
    U1T, S1T, Vt1T = np.linalg.svd(bd1_T, full_matrices=True)
    rk1_T = int(np.sum(np.abs(S1T) > 1e-8))

    if dim_om2_T > 0:
        bd2_T = build_full_boundary_matrix(paths_2_T, paths_1_T)
        bd2_om_T = bd2_T @ omega_2_T
        # B_1(T) in A_1(T) coordinates
        B1_T = bd2_om_T  # |A_1(T)| x dim_om2_T, columns span B_1(T)
    else:
        B1_T = np.zeros((len(paths_1_T), 0))

    # Embed Z_1(T\\v) into A_1(T)
    path_to_idx_T = {p: i for i, p in enumerate(paths_1_T)}
    path_to_idx_sub = {p: i for i, p in enumerate(paths_1_sub)}

    embed = np.zeros((len(paths_1_T), len(paths_1_sub)))
    for i, p_sub in enumerate(paths_1_sub):
        p_T = (relabel[p_sub[0]], relabel[p_sub[1]])
        if p_T in path_to_idx_T:
            embed[path_to_idx_T[p_T], i] = 1

    Z1_embedded = embed @ z1_sub_basis  # |A_1(T)| x z1_sub_dim

    # ker(i_*) = Z_1(T\\v) ∩ B_1(T) (in T coordinates), modulo B_1(T\\v)
    # Actually: ker(i_*) = {z in H_1(T\\v) : i_*(z) = 0 in H_1(T)}
    # = {z in Z_1(T\\v)/B_1(T\\v) : embed(z) in B_1(T)}
    # = (embed^{-1}(B_1(T)) ∩ Z_1(T\\v)) / B_1(T\\v)

    # Compute: {z in Z_1(T\\v) : embed(z) in span(B1_T)}
    if B1_T.shape[1] > 0:
        combined = np.hstack([Z1_embedded, B1_T])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        rk_B1T = np.linalg.matrix_rank(B1_T, tol=1e-8)
        # dim(Z_1(T\\v) ∩ B_1(T)) = z1_sub_dim + rk_B1T - rk_combined
        dim_intersect = z1_sub_dim + rk_B1T - rk_combined

        # But we need to account for B_1(T\\v) ⊂ Z_1(T\\v) ∩ B_1(T)
        # ker(i_*) = dim_intersect - rk_bd2_sub
        # Wait, not exactly. Let me think again.

        # Z_1(T\\v) = z1_sub_dim dimensional space
        # B_1(T\\v) = rk_bd2_sub dimensional subspace of Z_1(T\\v)
        # {z in Z_1(T\\v) : embed(z) in B_1(T)} = ?
        # This includes B_1(T\\v) (since embed(B_1(T\\v)) ⊂ B_1(T))
        # ker(i_*) = dim({z in Z_1(T\\v) : embed(z) in B_1(T)}) - dim(B_1(T\\v))
        ker_i = max(0, dim_intersect - rk_bd2_sub)
    else:
        # B_1(T) = 0, so ker(i_*) = {z : embed(z) = 0 in B_1(T)} / B_1(T\\v)
        # = Z_1(T\\v) ∩ ker(embed) / B_1(T\\v)
        # For tournaments embed is injective on Z_1, so this is 0
        ker_i = 0

    return ker_i


# ===== Verify: dim H_2^rel = dim ker(i_*) always? =====
print("=" * 70)
print("CONNECTING MAP INJECTIVITY CHECK")
print("Check: dim H_2(T,T\\v) = dim ker(H_1(T\\v)->H_1(T))")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n={n} exhaustive ---")
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs
    matches = 0
    mismatches = 0

    t0 = time.time()
    for bits in range(total):
        A = build_adj(n, bits)
        for v_test in range(n):
            h2_rel = compute_relative_h2_correct(A, n, v_test)
            ker_i = compute_h1_kernel(A, n, v_test)

            if h2_rel == ker_i:
                matches += 1
            else:
                mismatches += 1
                if mismatches <= 3:
                    scores = tuple(sorted(sum(row) for row in A))
                    print(f"  MISMATCH: bits={bits}, v={v_test}, scores={scores}, "
                          f"H_2^rel={h2_rel}, ker_i={ker_i}")

    dt = time.time() - t0
    print(f"  n={n}: matches={matches}, mismatches={mismatches} ({dt:.0f}s)")
    print(f"  delta injective for all: {'YES' if mismatches == 0 else 'NO'}")

print("\nDone.")
