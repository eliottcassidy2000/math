#!/usr/bin/env python3
"""
beta2_delta_sourcesink.py — Test: δ-injectivity fails ONLY at source/sink vertices

Hypothesis: δ: H₂(T,T\\v) → H₁(T\\v) is injective whenever 1 ≤ d⁺(v) ≤ n-2.
            It can fail when d⁺(v) = 0 (sink) or d⁺(v) = n-1 (source).

This would give the inductive proof of β₂=0:
  - For n≥4, pick any non-source/non-sink vertex v
  - By induction: β₂(T\\v) = 0
  - By δ-injectivity: H₂(T) ↪ ker(δ) = 0
  - QED

Quick test: only check non-source/non-sink vertices.

Author: opus-2026-03-08-S49
"""
import sys, time, random
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


def check_delta_injective_fast(A, n, v):
    """Fast check of δ-injectivity for (T, v)."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    ap0_T = enumerate_allowed_paths(A, n, 0)
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap3_T = enumerate_allowed_paths(A, n, 3)

    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)

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
        return True, 0

    ap2_T_list = [tuple(p) for p in ap2_T]

    # Embed Ω₂(T\v)
    if ap2_sub and d2_sub > 0:
        embed = np.zeros((len(ap2_T_list), d2_sub))
        for j in range(d2_sub):
            for k, path_sub in enumerate(ap2_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap2_T_list:
                    embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
        phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
    else:
        phi = np.zeros((d2_T, 0))

    rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
    if rk_phi > 0:
        U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
        Q = U_phi[:, rk_phi:]
    else:
        Q = np.eye(d2_T)

    d_rel = Q.shape[1]
    if d_rel == 0:
        return True, 0

    coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

    # Embed Ω₁(T\v)
    ap1_T_list = [tuple(p) for p in ap1_T]
    d1_T = dim_om(om1_T)
    if d1_sub > 0:
        embed1 = np.zeros((len(ap1_T_list), d1_sub))
        for j in range(d1_sub):
            for k, path_sub in enumerate(ap1_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap1_T_list:
                    embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
        psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
    else:
        psi = np.zeros((d1_T, 0))

    rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
    if rk_psi > 0:
        U_psi, _, _ = np.linalg.svd(psi, full_matrices=True)
        R = U_psi[:, rk_psi:]
    else:
        R = np.eye(d1_T)

    coords2_rel_q = R.T @ coords2_T @ Q
    rk_d2_rel = np.linalg.matrix_rank(coords2_rel_q, tol=1e-8)
    z2_rel_dim = d_rel - rk_d2_rel

    if z2_rel_dim == 0:
        return True, 0

    # Relative ∂₃
    if d3_T > 0:
        coords3_T = np.linalg.lstsq(om2_T, build_full_boundary_matrix(ap3_T, ap2_T) @ om3_T, rcond=None)[0]
        om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
        d3_sub = dim_om(om3_sub)

        d3_proj = Q.T @ coords3_T
        if d3_sub > 0:
            ap3_T_list = [tuple(p) for p in ap3_T]
            embed3 = np.zeros((len(ap3_T_list), d3_sub))
            for j in range(d3_sub):
                for k, path_sub in enumerate(ap3_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap3_T_list:
                        embed3[ap3_T_list.index(path_T), j] = om3_sub[k, j]
            chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
            rk_chi = np.linalg.matrix_rank(chi, tol=1e-8)
            if rk_chi > 0:
                U_chi, _, _ = np.linalg.svd(chi, full_matrices=True)
                d3_proj_rel = d3_proj @ U_chi[:, rk_chi:]
            else:
                d3_proj_rel = d3_proj
        else:
            d3_proj_rel = d3_proj
    else:
        d3_proj_rel = np.zeros((d_rel, 0))

    rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
    h2_rel = z2_rel_dim - rk_d3_rel

    if h2_rel == 0:
        return True, h2_rel

    # Full δ check
    U_rel, _, Vt_rel = np.linalg.svd(coords2_rel_q, full_matrices=True)
    z2_rel_basis = Vt_rel[rk_d2_rel:]

    if d3_proj_rel.shape[1] > 0:
        proj_z2 = z2_rel_basis @ d3_proj_rel
        rk_proj = np.linalg.matrix_rank(proj_z2, tol=1e-8)
    else:
        rk_proj = 0
    h2_rel_basis = z2_rel_basis[rk_proj:]
    if h2_rel_basis.shape[0] == 0:
        return True, 0

    delta_images = []
    for k in range(h2_rel_basis.shape[0]):
        h = h2_rel_basis[k]
        bd = coords2_T @ Q @ h
        if d1_sub > 0:
            delta_sub = np.linalg.lstsq(psi, bd, rcond=None)[0]
            delta_images.append(delta_sub)

    if not delta_images or d1_sub == 0:
        return True, h2_rel

    bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
    rk_d1_sub = np.linalg.matrix_rank(bd1_sub, tol=1e-8)
    if d2_sub > 0:
        coords2_sub = np.linalg.lstsq(om1_sub, build_full_boundary_matrix(ap2_sub, ap1_sub) @ om2_sub, rcond=None)[0]
        rk_d2_sub = np.linalg.matrix_rank(coords2_sub, tol=1e-8)
    else:
        rk_d2_sub = 0

    _, _, Vt1 = np.linalg.svd(bd1_sub @ om1_sub, full_matrices=True)
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
        delta_rk = np.linalg.matrix_rank(D_z1, tol=1e-8)

    return delta_rk == h2_rel, h2_rel


print("=" * 70)
print("SOURCE/SINK HYPOTHESIS FOR δ-INJECTIVITY")
print("=" * 70)

random.seed(42)

for n in [5, 6, 7]:
    m = n*(n-1)//2
    total = 1 << m

    if n <= 5:
        samples = list(range(total))
    elif n == 6:
        samples = list(range(total))
    else:
        samples = random.sample(range(total), min(500, total))

    fail_interior = 0  # failures at non-source/non-sink
    fail_boundary = 0  # failures at source/sink
    total_interior_checks = 0
    total_boundary_checks = 0
    exist_good_v = 0  # tournaments where SOME non-extremal v gives injectivity

    t0 = time.time()
    for idx, bits in enumerate(samples):
        A = build_adj(n, bits)
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

        has_good_v = False

        for v in range(n):
            dv = scores[v]
            is_extremal = (dv == 0 or dv == n-1)

            inj, h2_rel = check_delta_injective_fast(A, n, v)

            if h2_rel > 0:
                if is_extremal:
                    total_boundary_checks += 1
                    if not inj:
                        fail_boundary += 1
                else:
                    total_interior_checks += 1
                    if not inj:
                        fail_interior += 1
                        sc = tuple(sorted(scores))
                        print(f"  INTERIOR FAIL! T#{bits} scores={sc}, v={v}, d⁺(v)={dv}")
                    else:
                        has_good_v = True
            else:
                if not is_extremal:
                    has_good_v = True

        # Even if no interior v has h2_rel > 0, that's fine (H₂(T,T\v)=0 for all interior v)
        exist_good_v += 1  # we count based on existence

        if (idx + 1) % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {idx+1}/{len(samples)} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {len(samples)} tournaments in {elapsed:.0f}s")
    print(f"  Interior checks (1 ≤ d⁺ ≤ n-2): {total_interior_checks}")
    print(f"  Interior failures: {fail_interior}")
    print(f"  Boundary checks (source/sink): {total_boundary_checks}")
    print(f"  Boundary failures: {fail_boundary}")

    if fail_interior == 0:
        print(f"  ✓ δ is injective for ALL non-source/non-sink vertices at n={n}")
    else:
        print(f"  ✗ {fail_interior} interior failures at n={n}")


print("\nDone.")
