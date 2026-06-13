#!/usr/bin/env python3
"""
beta2_delta_scalar_debug.py — Debug the δ scalar computation

The scalar script showed |λ|=0 for 10 interior vertices, contradicting
the reliable delta-injectivity test. Find and fix the bug.

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
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


def check_delta_reliable(A, n, v):
    """Reliable δ-injectivity check from beta2_delta_injectivity.py"""
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
    d1_T = dim_om(om1_T)

    if d2_T == 0:
        return True, 0

    ap2_T_list = [tuple(p) for p in ap2_T]
    ap1_T_list = [tuple(p) for p in ap1_T]

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
        return True, 0

    # Full δ computation (RELIABLE version from beta2_delta_injectivity.py)
    _, _, Vt_rel = np.linalg.svd(coords2_rel_q, full_matrices=True)
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
print("DEBUG: Cross-checking scalar computation with reliable test")
print("=" * 70)

# Find the 10 interior cases with |λ|≈0 at n=5
n = 5
m = n*(n-1)//2

zero_cases = []
for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue
        if dv != 2:
            continue  # Only d+=2 had zeros

        # Quick check via reliable method
        inj, h2r = check_delta_reliable(A, n, v)
        if h2r > 0 and not inj:
            zero_cases.append((bits, v, dv, scores))
            print(f"  RELIABLE FAILURE: T#{bits} v={v} d+={dv} scores={scores}")

print(f"\nReliable failures for d+=2 interior: {len(zero_cases)}")
if len(zero_cases) == 0:
    print("  ✓ No reliable failures! The scalar computation had a bug.")
    print("  The bug is likely in the H₂(T,T\\v) generator computation.")

# Now identify where the scalar goes wrong
print(f"\n--- Diagnosing the scalar bug ---")
print("The issue: selecting h2_rel_basis = z2_rel_basis[rk_proj:]")
print("This assumes the first rk_proj rows of z2_rel_basis span B₂∩Z₂.")
print("But SVD of proj_z2 doesn't guarantee this ordering!")
print()
print("FIX: Properly compute H₂ = Z₂/(Z₂∩B₂) by projecting Z₂")
print("onto the orthogonal complement of Z₂∩im(∂₃^rel) within Z₂.")

# Let's find one specific zero case and trace through
for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv != 2:
            continue

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
        d1_T = dim_om(om1_T)
        d3_T = dim_om(om3_T)
        d1_sub = dim_om(om1_sub)
        d2_sub = dim_om(om2_sub)

        if d2_T == 0:
            continue

        ap2_T_list = [tuple(p) for p in ap2_T]
        ap1_T_list = [tuple(p) for p in ap1_T]

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
            continue

        coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

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
            continue

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
            continue

        # FOUND ONE — now trace the δ computation
        print(f"\n  TRACE: T#{bits} v={v} d+={dv}")
        print(f"    d_rel={d_rel}, z2_rel={z2_rel_dim}, rk_d3_rel={rk_d3_rel}, h2_rel={h2_rel}")

        # Z₂ basis
        _, _, Vt_rel = np.linalg.svd(coords2_rel_q, full_matrices=True)
        z2_rel_basis = Vt_rel[rk_d2_rel:]

        print(f"    Z₂^rel dim = {z2_rel_basis.shape[0]}")

        if d3_proj_rel.shape[1] > 0:
            proj_z2 = z2_rel_basis @ d3_proj_rel
            print(f"    proj_z2 shape = {proj_z2.shape}")
            print(f"    proj_z2 singular values: {np.linalg.svd(proj_z2, compute_uv=False)}")
            rk_proj = np.linalg.matrix_rank(proj_z2, tol=1e-8)
            print(f"    rk_proj = {rk_proj}")
        else:
            rk_proj = 0

        # The BUGGY selection
        h2_basis_buggy = z2_rel_basis[rk_proj:]
        print(f"    h2_basis_buggy shape = {h2_basis_buggy.shape}")

        # Now the CORRECT way: find the complement of im(∂₃^rel) within Z₂
        if rk_proj > 0:
            # Project out the B₂ directions from Z₂
            U_proj, S_proj, _ = np.linalg.svd(proj_z2.T, full_matrices=True)
            # The first rk_proj rows of z2_rel_basis span directions in B₂
            # But we need to properly orthogonalize
            B2_in_Z2 = proj_z2[:, :rk_proj] if proj_z2.shape[1] >= rk_proj else proj_z2
            # Actually, let's use a different approach: QR decomposition
            # to find which Z₂ directions map into im(∂₃^rel)

            # d3_proj_rel maps Ω₃^rel → rel Ω₂ (in Q-coords).
            # z2_rel_basis @ d3_proj_rel gives the projection of ∂₃ images onto Z₂.
            # A nonzero column means that ∂₃ direction intersects Z₂ (i.e., is in B₂).

            # The image of ∂₃^rel in Z₂ is: column space of proj_z2
            # To get H₂ = Z₂ / (Z₂ ∩ B₂), we need complement of colspace(proj_z2) in Z₂.

            # Since z2_rel_basis rows span Z₂, and proj_z2 = z2_rel_basis @ d3_proj_rel,
            # the column space of proj_z2 gives the B₂ directions within Z₂.

            # CORRECT H₂ generator: any Z₂ vector not in colspace(proj_z2)
            # Use SVD of proj_z2 to find the null space of proj_z2^T (left null space)
            U_z2, S_z2, Vt_z2 = np.linalg.svd(proj_z2, full_matrices=True)
            # U_z2[:, :rk_proj] spans the B₂∩Z₂ directions
            # U_z2[:, rk_proj:] spans the H₂ directions (complement in Z₂)
            h2_basis_correct = U_z2[:, rk_proj:].T @ z2_rel_basis

            print(f"    h2_basis_correct shape = {h2_basis_correct.shape}")

            # Compare
            for k in range(h2_basis_buggy.shape[0]):
                bd_buggy = coords2_T @ Q @ h2_basis_buggy[k]
                print(f"    Buggy δ[{k}] norm in Ω₁(T): {np.linalg.norm(bd_buggy):.6f}")

            for k in range(h2_basis_correct.shape[0]):
                bd_correct = coords2_T @ Q @ h2_basis_correct[k]
                print(f"    Correct δ[{k}] norm in Ω₁(T): {np.linalg.norm(bd_correct):.6f}")

            # Now compute δ correctly
            if d1_sub > 0:
                for k in range(h2_basis_correct.shape[0]):
                    bd = coords2_T @ Q @ h2_basis_correct[k]
                    delta_sub = np.linalg.lstsq(psi, bd, rcond=None)[0]
                    print(f"    Correct δ[{k}] in Ω₁(T\\v): norm={np.linalg.norm(delta_sub):.6f}")
                    print(f"      values: {delta_sub}")
        else:
            print(f"    rk_proj=0, so Z₂=H₂ (no B₂ to quotient)")

        break  # just one case
    else:
        continue
    break

print("\nDone.")
