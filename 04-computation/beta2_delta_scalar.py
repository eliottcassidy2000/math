#!/usr/bin/env python3
"""
beta2_delta_scalar.py — The δ map as a SINGLE SCALAR

Since h₂_rel = 1 and β₁(T\v) = 1 universally (when nonzero),
δ: H₂(T,T\v) → H₁(T\v) is a map R → R, i.e., multiplication by a scalar λ.

δ-injectivity ⟺ λ ≠ 0.

This script computes λ for all (T,v) with interior v and h₂_rel = 1.
We look for:
- Is λ always nonzero? (confirms HYP-249)
- Is |λ| related to any tournament invariant?
- Is λ always ±1 or some other nice value?
- What is λ for source vertices? (should be 0 for the 4 failure cases)

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


def compute_delta_scalar(A, n, v):
    """Compute the scalar value of δ: H₂(T,T\\v) → H₁(T\\v).
    Returns (lambda_val, h2_rel, beta1) or (None, 0, beta1) if h2_rel=0."""
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

    # β₁(T\v)
    bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
    rk_d1 = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)
    z1_sub_dim = d1_sub - rk_d1
    if d2_sub > 0:
        bd2_sub_mat = build_full_boundary_matrix(ap2_sub, ap1_sub)
        bd2_in_om1 = np.linalg.lstsq(om1_sub, bd2_sub_mat @ om2_sub, rcond=None)[0]
        rk_d2_sub = np.linalg.matrix_rank(bd2_in_om1, tol=1e-8)
    else:
        bd2_in_om1 = np.zeros((d1_sub, 0))
        rk_d2_sub = 0
    beta1 = z1_sub_dim - rk_d2_sub

    if d2_T == 0:
        return None, 0, beta1

    ap2_T_list = [tuple(p) for p in ap2_T]
    ap1_T_list = [tuple(p) for p in ap1_T]

    # Relative Ω₂
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
        return None, 0, beta1

    # Relative ∂₂
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
        return None, 0, beta1

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
        return None, 0, beta1

    # Now compute δ as a scalar.
    # H₂(T,T\v) is 1-dimensional; H₁(T\v) is 1-dimensional.
    # Need to find the H₂(T,T\v) generator and compute its δ image in H₁(T\v).

    # Step 1: Find Z₂^rel = ker(∂₂^rel)
    _, _, Vt_rel = np.linalg.svd(coords2_rel_q, full_matrices=True)
    z2_rel_basis = Vt_rel[rk_d2_rel:]  # rows are Z₂^rel basis vectors (in Q-coords)

    # Step 2: Mod out by B₂^rel = im(∂₃^rel)
    if d3_proj_rel.shape[1] > 0 and rk_d3_rel > 0:
        # Project Z₂^rel onto complement of B₂^rel
        proj_z2 = z2_rel_basis @ d3_proj_rel
        rk_proj = np.linalg.matrix_rank(proj_z2, tol=1e-8)
        if rk_proj > 0:
            # Some Z₂ directions are in B₂
            U_z2, _, _ = np.linalg.svd(proj_z2, full_matrices=True)
            h2_generator_in_z2 = z2_rel_basis[rk_proj:]
        else:
            h2_generator_in_z2 = z2_rel_basis
    else:
        h2_generator_in_z2 = z2_rel_basis

    if h2_generator_in_z2.shape[0] == 0:
        return None, 0, beta1

    # h2_generator_in_z2 should be 1D (h2_rel=1)
    h_vec = h2_generator_in_z2[0]  # in relative Ω₂ coords (projected by Q)

    # Step 3: Compute δ(h) = ∂₂(lift(h)) restricted to Ω₁(T\v), modulo B₁(T\v)
    # lift(h) = Q @ h_vec (in Ω₂(T) coords)
    # ∂₂(lift(h)) = coords2_T @ Q @ h_vec (in Ω₁(T) coords)
    bd_image = coords2_T @ Q @ h_vec  # in Ω₁(T) coords

    # Step 4: Project bd_image onto Ω₁(T\v) via psi
    if d1_sub > 0:
        # psi: Ω₁(T) → Ω₁(T\v). We need the component IN Ω₁(T\v).
        # bd_image is in Ω₁(T) coordinates. The T\v part is psi^T @ bd_image... no.
        # psi maps Ω₁(T\v) basis INTO Ω₁(T) basis. To get the T\v component:
        # bd_image_sub = psi^T @ bd_image (but this isn't right either)

        # Actually, the connecting map δ works like this:
        # Take z ∈ Z₂^rel = {x ∈ Ω₂(T) : ∂₂(x) ∈ im(Ω₁(T\v)→Ω₁(T))}
        # Then ∂₂(x) = ι(y) for some y ∈ Ω₁(T\v)
        # δ([z]) = [y] ∈ H₁(T\v)

        # So we need to solve: psi @ y = bd_image for y
        # Wait, psi is in Ω₁(T) coords. ι is the inclusion Ω₁(T\v) → Ω₁(T).
        # In our coordinates: the image of Ω₁(T\v) in Ω₁(T) is column space of psi.
        # So we solve: psi @ y = bd_image

        # But wait, the relative cycle condition says ∂₂(x) ∈ im(Ω₁(T\v)).
        # This means bd_image is NOT necessarily in im(psi).
        # Actually, the relative cycle condition IS that ∂₂^rel = 0,
        # where ∂₂^rel: Ω₂^rel → Ω₁^rel.
        # So R.T @ coords2_T @ Q @ h_vec = 0, which means coords2_T @ Q @ h_vec ∈ im(psi).

        # So bd_image = psi @ y for some y.
        y = np.linalg.lstsq(psi, bd_image, rcond=None)[0]
        residual = np.max(np.abs(psi @ y - bd_image))
        if residual > 1e-6:
            print(f"  WARNING: residual {residual:.6f} in δ solve")
    else:
        return None, h2_rel, beta1

    # Step 5: Project y onto H₁(T\v) = Z₁/B₁
    # Z₁ basis in Ω₁(T\v) coords:
    _, _, Vt1 = np.linalg.svd(bd1_sub @ om1_sub, full_matrices=True)
    z1_basis = Vt1[rk_d1:]  # (z1_dim, d1_sub)

    # B₁ basis in Z₁ coords:
    if rk_d2_sub > 0:
        B1_in_z1 = z1_basis @ bd2_in_om1  # (z1_dim, d2_sub)
        rk_B1_in_z1 = np.linalg.matrix_rank(B1_in_z1, tol=1e-8)
    else:
        B1_in_z1 = np.zeros((z1_sub_dim, 0))
        rk_B1_in_z1 = 0

    # H₁ generator: a direction in Z₁ not in B₁
    if rk_B1_in_z1 > 0:
        U_B1, _, _ = np.linalg.svd(B1_in_z1, full_matrices=True)
        h1_gen = z1_basis[rk_B1_in_z1:]  # (beta1, d1_sub) — rows
    else:
        h1_gen = z1_basis

    # Project y onto Z₁
    y_z1 = z1_basis @ y  # (z1_dim,)

    # Project y_z1 onto H₁ = Z₁/B₁
    y_h1 = h1_gen @ y  # (beta1,) — the δ value in H₁

    # Since h₂_rel = 1 and β₁ = 1, this is a single scalar
    if len(y_h1) == 1:
        lambda_val = y_h1[0]
    else:
        lambda_val = np.linalg.norm(y_h1)

    return lambda_val, h2_rel, beta1


print("=" * 70)
print("δ AS A SCALAR: λ VALUES FOR ALL (T,v)")
print("=" * 70)

# Part 1: n=5 exhaustive
n = 5
m = n*(n-1)//2

print(f"\n--- n={n} exhaustive ---")

interior_lambdas = []
source_lambdas = []
sink_lambdas = []

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        lam, h2r, b1 = compute_delta_scalar(A, n, v)

        if h2r > 0:
            if dv == 0:
                sink_lambdas.append((bits, v, dv, lam, h2r, b1))
            elif dv == n-1:
                source_lambdas.append((bits, v, dv, lam, h2r, b1))
            else:
                interior_lambdas.append((bits, v, dv, lam, h2r, b1))

print(f"\nInterior vertices with h₂_rel > 0: {len(interior_lambdas)}")
if interior_lambdas:
    vals = [abs(x[3]) for x in interior_lambdas]
    print(f"  |λ| range: [{min(vals):.6f}, {max(vals):.6f}]")
    print(f"  |λ| mean: {np.mean(vals):.6f}")
    print(f"  |λ| = 0 count: {sum(1 for v in vals if v < 1e-8)}")

    # Distribution by d+
    by_dplus = {}
    for bits, v, dv, lam, h2r, b1 in interior_lambdas:
        by_dplus.setdefault(dv, []).append(abs(lam))
    for dv in sorted(by_dplus):
        vals_dv = by_dplus[dv]
        print(f"  d⁺={dv}: {len(vals_dv)} cases, |λ| range [{min(vals_dv):.4f}, {max(vals_dv):.4f}]")

    # Show a few specific values
    print(f"\n  Sample interior λ values:")
    for bits, v, dv, lam, h2r, b1 in interior_lambdas[:10]:
        print(f"    T#{bits} v={v} d⁺={dv}: λ={lam:.6f}")

print(f"\nSource vertices with h₂_rel > 0: {len(source_lambdas)}")
for bits, v, dv, lam, h2r, b1 in source_lambdas:
    print(f"  T#{bits} v={v}: λ={lam:.6f} (h₂_rel={h2r}, β₁={b1})")

print(f"\nSink vertices with h₂_rel > 0: {len(sink_lambdas)}")
for bits, v, dv, lam, h2r, b1 in sink_lambdas:
    print(f"  T#{bits} v={v}: λ={lam:.6f} (h₂_rel={h2r}, β₁={b1})")


# Part 2: Look for pattern in |λ|
print(f"\n{'='*70}")
print("PATTERN SEARCH: What determines |λ|?")
print("=" * 70)

if interior_lambdas:
    # Group by score sequence
    by_scores = {}
    for bits, v, dv, lam, h2r, b1 in interior_lambdas:
        A = build_adj(n, bits)
        scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
        by_scores.setdefault(scores, []).append(abs(lam))

    print("\n  By score sequence:")
    for scores, vals in sorted(by_scores.items()):
        print(f"    {scores}: {len(vals)} cases, |λ| values: {sorted(set(round(v, 4) for v in vals))}")


# Part 3: n=6 sample
print(f"\n{'='*70}")
print("n=6 (sampled)")
print("=" * 70)

random.seed(42)
n = 6
m = n*(n-1)//2
total = 1 << m
samples = random.sample(range(total), 2000)

interior_lambdas_6 = []
t0 = time.time()

for idx, bits in enumerate(samples):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        lam, h2r, b1 = compute_delta_scalar(A, n, v)
        if h2r > 0:
            interior_lambdas_6.append((bits, v, dv, lam, h2r, b1))
        break  # one v per tournament

    if (idx + 1) % 500 == 0:
        print(f"  {idx+1}/2000 ({time.time()-t0:.0f}s)")

print(f"\nn=6: {len(interior_lambdas_6)} interior cases with h₂_rel > 0")
if interior_lambdas_6:
    vals = [abs(x[3]) for x in interior_lambdas_6]
    print(f"  |λ| range: [{min(vals):.6f}, {max(vals):.6f}]")
    print(f"  |λ| = 0 count: {sum(1 for v in vals if v < 1e-8)}")
    print(f"  |λ| mean: {np.mean(vals):.6f}")

    # Check if all nonzero
    if all(v > 1e-8 for v in vals):
        print(f"  ✓ ALL |λ| > 0 (δ injective for all interior v)")
    else:
        print(f"  ✗ Some λ = 0!")

    # Unique values
    unique_vals = sorted(set(round(v, 4) for v in vals))
    if len(unique_vals) <= 20:
        print(f"  Unique |λ| values: {unique_vals}")

print("\nDone.")
