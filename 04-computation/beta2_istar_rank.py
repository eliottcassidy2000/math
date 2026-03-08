#!/usr/bin/env python3
"""
beta2_istar_rank.py — Verify: β₂=0 ⟺ i_*: H₁(T\v) → H₁(T) has maximal rank

From the LES: β₂(T) = H₂^rel - dim(ker i_*)
We proved H₂^rel = max(0, β₁(T\v) - β₁(T)).

So β₂(T) = max(0, β₁(T\v) - β₁(T)) - dim(ker i_*)
         = max(0, β₁(T\v) - β₁(T)) - (β₁(T\v) - rk(i_*))

Case β₁(T\v) ≥ β₁(T):
  β₂ = β₁(T\v) - β₁(T) - β₁(T\v) + rk(i_*) = rk(i_*) - β₁(T)
  β₂=0 ⟺ rk(i_*) = β₁(T) (surjective)

Case β₁(T\v) < β₁(T):
  β₂ = 0 - β₁(T\v) + rk(i_*) = rk(i_*) - β₁(T\v)
  β₂=0 ⟺ rk(i_*) = β₁(T\v) (injective)

In both cases: β₂=0 ⟺ rk(i_*) = min(β₁(T\v), β₁(T))

This is a CLEAN characterization. Verify it and then prove i_* has max rank.

Author: opus-2026-03-08-S45
"""
import sys
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


def local_to_global_map(n, v):
    return [i for i in range(n) if i != v]


def get_h1_basis(A, n):
    """Return H₁ basis as vectors in Ω₁ coordinates."""
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)

    if not ap1:
        return np.zeros((0, 0)), [], ap1, ap0

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    d1 = om1.shape[1] if om1.ndim == 2 else 0
    if d1 == 0:
        return np.zeros((0, 0)), [], ap1, ap0

    # ∂₁: Ω₁ → Ω₀
    bd1 = build_full_boundary_matrix(ap1, ap0)
    bd1_om = bd1 @ om1
    om0 = np.eye(n)
    coords1 = np.linalg.lstsq(om0, bd1_om, rcond=None)[0]
    U1, S1, Vt1 = np.linalg.svd(coords1, full_matrices=True)
    rk1 = int(sum(s > 1e-8 for s in S1))

    # Z₁ basis (in Ω₁ coords)
    if rk1 < Vt1.shape[0]:
        z1_basis = Vt1[rk1:]  # rows are Z₁ basis vectors
    else:
        return np.zeros((d1, 0)), om1, ap1, ap0

    # B₁ from ∂₂
    if ap2:
        om2 = compute_omega_basis(A, n, 2, ap2, ap1)
        d2 = om2.shape[1] if om2.ndim == 2 else 0
        if d2 > 0:
            bd2 = build_full_boundary_matrix(ap2, ap1)
            bd2_om = bd2 @ om2
            coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
            rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        else:
            rk2 = 0
            coords2 = np.zeros((d1, 0))
    else:
        rk2 = 0
        coords2 = np.zeros((d1, 0))

    # H₁ = Z₁ / B₁
    # Project Z₁ basis modulo B₁
    # H₁ basis vectors in Ω₁ coords: columns of z1_basis.T that are
    # linearly independent mod im(coords2)
    z1_cols = z1_basis.T  # d1 × z1_dim

    beta1 = z1_cols.shape[1] - rk2
    return z1_cols, om1, ap1, ap0  # return Z₁ basis


def compute_istar_rank(A, n, v):
    """Compute rank of i_*: H₁(T\v) → H₁(T)."""
    B = delete_vertex(A, n, v)
    loc2glob = local_to_global_map(n, v)

    # H₁(T) data
    ap0_T = enumerate_allowed_paths(A, n, 0)
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)

    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T) if ap1_T else np.zeros((0,0))
    d1_T = om1_T.shape[1] if om1_T.ndim == 2 and om1_T.shape[0] > 0 else 0

    if d1_T == 0:
        return 0, 0, 0  # rk, beta1_T, beta1_Tv

    bd1_T = build_full_boundary_matrix(ap1_T, ap0_T)
    bd1_om_T = bd1_T @ om1_T
    coords1_T = np.linalg.lstsq(np.eye(n), bd1_om_T, rcond=None)[0]
    rk_d1_T = np.linalg.matrix_rank(coords1_T, tol=1e-8)

    if ap2_T:
        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T)
        d2_T = om2_T.shape[1] if om2_T.ndim == 2 else 0
        if d2_T > 0:
            bd2_T = build_full_boundary_matrix(ap2_T, ap1_T)
            bd2_om_T = bd2_T @ om2_T
            coords2_T = np.linalg.lstsq(om1_T, bd2_om_T, rcond=None)[0]
            rk_d2_T = np.linalg.matrix_rank(coords2_T, tol=1e-8)
        else:
            rk_d2_T = 0
    else:
        rk_d2_T = 0

    beta1_T = (d1_T - rk_d1_T) - rk_d2_T

    # H₁(T\v) data
    ap0_Tv = enumerate_allowed_paths(B, n-1, 0)
    ap1_Tv = enumerate_allowed_paths(B, n-1, 1)
    ap2_Tv = enumerate_allowed_paths(B, n-1, 2)

    om1_Tv = compute_omega_basis(B, n-1, 1, ap1_Tv, ap0_Tv) if ap1_Tv else np.zeros((0,0))
    d1_Tv = om1_Tv.shape[1] if om1_Tv.ndim == 2 and om1_Tv.shape[0] > 0 else 0

    if d1_Tv == 0:
        return 0, beta1_T, 0

    bd1_Tv = build_full_boundary_matrix(ap1_Tv, ap0_Tv)
    bd1_om_Tv = bd1_Tv @ om1_Tv
    coords1_Tv = np.linalg.lstsq(np.eye(n-1), bd1_om_Tv, rcond=None)[0]
    rk_d1_Tv = np.linalg.matrix_rank(coords1_Tv, tol=1e-8)

    if ap2_Tv:
        om2_Tv = compute_omega_basis(B, n-1, 2, ap2_Tv, ap1_Tv)
        d2_Tv = om2_Tv.shape[1] if om2_Tv.ndim == 2 else 0
        if d2_Tv > 0:
            bd2_Tv = build_full_boundary_matrix(ap2_Tv, ap1_Tv)
            bd2_om_Tv = bd2_Tv @ om2_Tv
            coords2_Tv = np.linalg.lstsq(om1_Tv, bd2_om_Tv, rcond=None)[0]
            rk_d2_Tv = np.linalg.matrix_rank(coords2_Tv, tol=1e-8)
        else:
            rk_d2_Tv = 0
    else:
        rk_d2_Tv = 0

    beta1_Tv = (d1_Tv - rk_d1_Tv) - rk_d2_Tv

    if beta1_Tv == 0:
        return 0, beta1_T, 0

    # Now compute rk(i_*: H₁(T\v) → H₁(T))
    # Z₁(T\v) basis in Ω₁(T\v) coords
    U_Tv, S_Tv, Vt_Tv = np.linalg.svd(coords1_Tv, full_matrices=True)
    rk_z = int(sum(s > 1e-8 for s in S_Tv))
    if rk_z < Vt_Tv.shape[0]:
        z1_Tv = Vt_Tv[rk_z:].T  # d1_Tv × z1_dim, cols are Z₁ basis
    else:
        return 0, beta1_T, beta1_Tv

    # Embed Z₁(T\v) into A₁(T) coords
    # First: Z₁(T\v) in A₁(T\v) coords
    z1_Tv_A = om1_Tv @ z1_Tv  # |A₁(T\v)| × z1_dim

    # Inclusion map: A₁(T\v) → A₁(T) with vertex renumbering
    ap1_T_list = [tuple(x) for x in ap1_T]
    ap1_Tv_list = [tuple(x) for x in ap1_Tv]
    T_idx = {p: i for i, p in enumerate(ap1_T_list)}

    incl = np.zeros((len(ap1_T_list), len(ap1_Tv_list)))
    for j, path_local in enumerate(ap1_Tv_list):
        path_global = tuple(loc2glob[k] for k in path_local)
        if path_global in T_idx:
            incl[T_idx[path_global], j] = 1

    # Z₁(T\v) embedded in A₁(T)
    z1_embedded = incl @ z1_Tv_A  # |A₁(T)| × z1_dim

    # Project into Ω₁(T) coords
    z1_in_om1_T = np.linalg.lstsq(om1_T, z1_embedded, rcond=None)[0]  # d1_T × z1_dim

    # Now project into H₁(T) = Z₁(T) / B₁(T)
    # Z₁(T) basis in Ω₁(T) coords
    U_T, S_T, Vt_T = np.linalg.svd(coords1_T, full_matrices=True)
    rk_z_T = int(sum(s > 1e-8 for s in S_T))
    if rk_z_T < Vt_T.shape[0]:
        z1_T = Vt_T[rk_z_T:].T  # d1_T × z1_dim_T
    else:
        return 0, beta1_T, beta1_Tv

    # Project embedded Z₁(T\v) into Z₁(T)
    # z1_in_om1_T should already be in Z₁ if i_* is well-defined
    # Project: z1_T^T @ z1_in_om1_T
    z1_proj = z1_T.T @ z1_in_om1_T  # z1_dim_T × z1_dim_Tv

    # Now mod out B₁(T)
    # B₁(T) in Ω₁(T) coords
    if rk_d2_T > 0:
        b1_T = coords2_T[:, :rk_d2_T]  # columns of B₁
        b1_in_z1 = z1_T.T @ b1_T  # z1_dim_T × rk_d2_T

        # Combined: [b1_in_z1 | z1_proj]
        combined = np.hstack([b1_in_z1, z1_proj])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        rk_b1 = np.linalg.matrix_rank(b1_in_z1, tol=1e-8)
        rk_istar = rk_combined - rk_b1
    else:
        rk_istar = np.linalg.matrix_rank(z1_proj, tol=1e-8)

    return rk_istar, beta1_T, beta1_Tv


# ===== MAIN =====
print("=" * 70)
print("VERIFY: rk(i_*) = min(β₁(T\\v), β₁(T)) FOR ALL (T,v)")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

maxrk_count = 0
nonmaxrk_count = 0
total = 0
details = Counter()

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    for v in range(n):
        rk, b1T, b1Tv = compute_istar_rank(A, n, v)
        expected = min(b1T, b1Tv)
        total += 1

        if rk == expected:
            maxrk_count += 1
            details[('max', b1T, b1Tv, rk)] += 1
        else:
            nonmaxrk_count += 1
            details[('NON-MAX', b1T, b1Tv, rk)] += 1
            if nonmaxrk_count <= 5:
                print(f"  NON-MAX: bits={bits}, v={v}, rk={rk}, "
                      f"β₁(T)={b1T}, β₁(T\\v)={b1Tv}, expected={expected}")

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}, non-max={nonmaxrk_count}")

print(f"\nResults ({total} pairs):")
print(f"  Max rank: {maxrk_count}/{total}")
print(f"  Non-max rank: {nonmaxrk_count}/{total}")

print(f"\nBreakdown:")
for key, count in sorted(details.items()):
    print(f"  {key}: {count}")

if nonmaxrk_count == 0:
    print(f"\n✓ i_* ALWAYS has maximal rank!")
    print(f"✓ This proves β₂(T) = 0 for all n={n} tournaments (by induction)")

print("\nDone.")
