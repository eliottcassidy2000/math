#!/usr/bin/env python3
"""
beta2_delta_failures.py — Identify the 4 δ-injectivity failures at n=5

At n=5, δ: H₂(T,T\\v) → H₁(T\\v) fails injectivity for exactly 4 (T,v) pairs.
Yet for EVERY tournament, SOME v gives injectivity.

Question: What makes these specific pairs fail?

Author: opus-2026-03-08-S49
"""
import sys
import numpy as np
from collections import Counter, defaultdict
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


# Import the checking function from the main script
# (Redefine inline for speed)
def check_delta_detailed(A, n, v):
    """Check δ-injectivity with detailed output."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

    ap0_T = enumerate_allowed_paths(A, n, 0)
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap3_T = enumerate_allowed_paths(A, n, 3)

    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)

    remap = {i: others[i] for i in range(n1)}

    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))

    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

    d2_T = dim_om(om2_T)
    d3_T = dim_om(om3_T)
    d1_sub = dim_om(om1_sub)
    d2_sub = dim_om(om2_sub)

    if d2_T == 0:
        return True, {}

    # Embed Ω₂(T\v) into Ω₂(T)
    ap2_T_list = [tuple(p) for p in ap2_T]
    ap2_sub_list = [tuple(p) for p in ap2_sub] if ap2_sub else []

    if ap2_sub and d2_sub > 0:
        embed = np.zeros((len(ap2_T_list), d2_sub))
        for j in range(d2_sub):
            a2_sub_vec = om2_sub[:, j]
            for k, path_sub in enumerate(ap2_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap2_T_list:
                    idx_T = ap2_T_list.index(path_T)
                    embed[idx_T, j] = a2_sub_vec[k]
        phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
    else:
        phi = np.zeros((d2_T, 0))

    rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)

    # Quotient basis
    if rk_phi > 0:
        U_phi, S_phi, _ = np.linalg.svd(phi, full_matrices=True)
        Q = U_phi[:, rk_phi:]
    else:
        Q = np.eye(d2_T)

    d_rel = Q.shape[1]
    if d_rel == 0:
        return True, {}

    # Boundaries in Ω coords
    coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

    if d3_T > 0:
        coords3_T = np.linalg.lstsq(om2_T, build_full_boundary_matrix(ap3_T, ap2_T) @ om3_T, rcond=None)[0]
    else:
        coords3_T = np.zeros((d2_T, 0))

    # Embed Ω₁(T\v) into Ω₁(T)
    ap1_T_list = [tuple(p) for p in ap1_T]
    d1_T = dim_om(om1_T)

    if d1_sub > 0:
        embed1 = np.zeros((len(ap1_T_list), d1_sub))
        for j in range(d1_sub):
            a1_sub_vec = om1_sub[:, j]
            for k, path_sub in enumerate(ap1_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap1_T_list:
                    idx_T = ap1_T_list.index(path_T)
                    embed1[idx_T, j] = a1_sub_vec[k]
        psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
    else:
        psi = np.zeros((d1_T, 0))

    rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)

    if rk_psi > 0:
        U_psi, _, _ = np.linalg.svd(psi, full_matrices=True)
        R = U_psi[:, rk_psi:]
    else:
        R = np.eye(d1_T)

    # Relative boundary
    coords2_rel_q = R.T @ coords2_T @ Q
    rk_d2_rel = np.linalg.matrix_rank(coords2_rel_q, tol=1e-8)
    z2_rel_dim = d_rel - rk_d2_rel

    if z2_rel_dim == 0:
        return True, {'h2_rel': 0}

    # Relative ∂₃ image
    ap3_T_list = [tuple(p) for p in ap3_T]

    om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
    d3_sub = dim_om(om3_sub)

    if d3_T > 0:
        if d3_sub > 0:
            embed3 = np.zeros((len(ap3_T_list), d3_sub))
            for j in range(d3_sub):
                a3_sub_vec = om3_sub[:, j]
                for k, path_sub in enumerate(ap3_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap3_T_list:
                        idx_T = ap3_T_list.index(path_T)
                        embed3[idx_T, j] = a3_sub_vec[k]
            chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
            rk_chi = np.linalg.matrix_rank(chi, tol=1e-8)
        else:
            chi = np.zeros((d3_T, 0))
            rk_chi = 0

        d3_proj = Q.T @ coords3_T
        if rk_chi > 0:
            U_chi, _, _ = np.linalg.svd(chi, full_matrices=True)
            Q3 = U_chi[:, rk_chi:]
            d3_proj_rel = d3_proj @ Q3
        else:
            d3_proj_rel = d3_proj
    else:
        d3_proj_rel = np.zeros((d_rel, 0))

    rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
    h2_rel = z2_rel_dim - rk_d3_rel

    info = {
        'h2_rel': h2_rel,
        'z2_rel': z2_rel_dim,
        'rk_d3_rel': rk_d3_rel,
        'd_rel': d_rel,
        'd2_T': d2_T, 'd3_T': d3_T,
        'd2_sub': d2_sub, 'd3_sub': d3_sub,
    }

    if h2_rel == 0:
        return True, info

    # Compute δ and check injectivity
    U_rel, S_rel, Vt_rel = np.linalg.svd(coords2_rel_q, full_matrices=True)
    z2_rel_basis = Vt_rel[rk_d2_rel:]

    if d3_proj_rel.shape[1] > 0:
        proj_z2 = z2_rel_basis @ d3_proj_rel
        rk_proj = np.linalg.matrix_rank(proj_z2, tol=1e-8)
        U_p, S_p, _ = np.linalg.svd(proj_z2, full_matrices=True)
        h2_basis_idx = rk_proj
    else:
        h2_basis_idx = 0

    h2_rel_basis = z2_rel_basis[h2_basis_idx:]

    # Compute δ for each H₂^rel basis vector
    delta_images = []
    for k in range(h2_rel_basis.shape[0]):
        h = h2_rel_basis[k]
        bd = coords2_T @ Q @ h
        if d1_sub > 0:
            delta_sub = np.linalg.lstsq(psi, bd, rcond=None)[0]
            delta_images.append(delta_sub)

    if not delta_images or d1_sub == 0:
        return True, info

    # Check in H₁(T\v)
    bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
    rk_d1_sub = np.linalg.matrix_rank(bd1_sub, tol=1e-8)

    if d2_sub > 0:
        coords2_sub = np.linalg.lstsq(om1_sub, build_full_boundary_matrix(ap2_sub, ap1_sub) @ om2_sub, rcond=None)[0]
        rk_d2_sub = np.linalg.matrix_rank(coords2_sub, tol=1e-8)
    else:
        rk_d2_sub = 0

    z1_sub = d1_sub - rk_d1_sub
    beta1_sub = z1_sub - rk_d2_sub
    info['beta1_sub'] = beta1_sub

    U1, S1, Vt1 = np.linalg.svd(bd1_sub @ om1_sub, full_matrices=True)
    z1_basis = Vt1[rk_d1_sub:]

    D = np.array(delta_images)
    D_z1 = D @ z1_basis.T

    if rk_d2_sub > 0:
        B1_z1 = z1_basis @ coords2_sub
        rk_B1 = np.linalg.matrix_rank(B1_z1, tol=1e-8)
        combined = np.hstack([B1_z1, D_z1.T])
        rk_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        delta_rk = rk_combined - rk_B1
    else:
        rk_B1 = 0
        delta_rk = np.linalg.matrix_rank(D_z1, tol=1e-8)

    info['delta_rk'] = delta_rk
    is_injective = (delta_rk == h2_rel)

    return is_injective, info


print("=" * 70)
print("IDENTIFYING δ-INJECTIVITY FAILURES AT n=5")
print("=" * 70)

n = 5
m = n*(n-1)//2

for bits in range(1 << m):
    A = build_adj(n, bits)

    for v in range(n):
        inj, info = check_delta_detailed(A, n, v)
        if not inj:
            scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
            d_out_v = sum(A[v][j] for j in range(n) if j != v)

            print(f"\nFAILURE: T#{bits} scores={scores}, v={v} (d⁺(v)={d_out_v})")
            print(f"  Details: {info}")

            # Show the tournament
            print(f"  Adjacency:")
            for i in range(n):
                row = [A[i][j] for j in range(n)]
                print(f"    {i}: {row}")

            # Which vertices beat v?
            beats_v = [i for i in range(n) if i != v and A[i][v]]
            beaten_by_v = [i for i in range(n) if i != v and A[v][i]]
            print(f"  Vertices beating v: {beats_v}")
            print(f"  Vertices beaten by v: {beaten_by_v}")

            # What is T\v?
            others = [i for i in range(n) if i != v]
            n1 = n - 1
            A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
            scores_sub = tuple(sorted(sum(A_sub[i][j] for j in range(n1) if j!=i) for i in range(n1)))
            c3_sub = sum(1 for i in range(n1) for j in range(i+1,n1) for k in range(j+1,n1)
                         if A_sub[i][j]+A_sub[j][k]+A_sub[k][i] in (0,3))
            print(f"  T\\v scores={scores_sub}, c₃={c3_sub}")

print("\nDone.")
