#!/usr/bin/env python3
"""
beta2_sigma_bound_algebra.py — Algebraic mechanism for Σ h₂_rel ≤ 3

Key question: WHY is the sum of critical fillers bounded by 3?

Approach:
1. For each T, compute Σ_v β₁(T\v) and Σ h₂_rel
2. The "killing rate" = Σ h₂_rel / Σ β₁(T\v) measures how often
   a nonzero H₁(T\v) is killed by inclusion
3. Study the Z₂ dimension and its relationship to fillers
4. The rank of the collective inclusion map i*: ⊕H₁(T\v) → H₁(T)

Key identity (from LES):
For each v with β₁(T\v) > 0:
  h₂_rel(v) = β₁(T\v) - rk(i*_v) = dim ker(i*_v: H₁(T\v) → H₁(T))

Since β₁(T\v) = 1 universally (HYP-257), h₂_rel(v) ∈ {0,1}, and:
  h₂_rel(v) = 1 ⟺ i*_v = 0 (the H₁ class is killed)

So Σ h₂_rel = #{v : β₁(T\v) = 1 and i*_v = 0}

Author: opus-2026-03-08-S49
"""
import sys, time
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


def analyze_tournament(A, n):
    """Full analysis of β₁ deletions and H₁-killing for tournament A."""
    # Compute T data
    ap = {}
    om = {}
    for p in range(4):
        ap[p] = enumerate_allowed_paths(A, n, p)
    om[0] = np.eye(n)
    om[1] = compute_omega_basis(A, n, 1, ap[1], ap[0])
    om[2] = compute_omega_basis(A, n, 2, ap[2], ap[1]) if ap[2] else np.zeros((0,0))
    om[3] = compute_omega_basis(A, n, 3, ap[3], ap[2]) if ap[3] else np.zeros((0,0))

    # β₁(T)
    bd1 = build_full_boundary_matrix(ap[1], ap[0])
    d1m = np.linalg.lstsq(om[0], bd1 @ om[1], rcond=None)[0]
    rk1 = np.linalg.matrix_rank(d1m, tol=1e-8)
    d1 = dim_om(om[1])

    if dim_om(om[2]) > 0:
        bd2 = build_full_boundary_matrix(ap[2], ap[1])
        d2m = np.linalg.lstsq(om[1], bd2 @ om[2], rcond=None)[0]
        im2 = np.linalg.matrix_rank(d2m, tol=1e-8)
    else:
        im2 = 0
    beta1_T = (d1 - rk1) - im2

    # Z₂ dimension
    d2 = dim_om(om[2])
    if d2 > 0:
        bd2_T = build_full_boundary_matrix(ap[2], ap[1])
        d2_mat = np.linalg.lstsq(om[1], bd2_T @ om[2], rcond=None)[0]
        rk2 = np.linalg.matrix_rank(d2_mat, tol=1e-8)
        z2_dim = d2 - rk2
    else:
        z2_dim = 0

    # For each vertex: β₁(T\v) and h₂_rel
    results = []
    for v in range(n):
        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

        ap0s = enumerate_allowed_paths(A_sub, n1, 0)
        ap1s = enumerate_allowed_paths(A_sub, n1, 1)
        ap2s = enumerate_allowed_paths(A_sub, n1, 2)
        om1s = compute_omega_basis(A_sub, n1, 1, ap1s, ap0s)
        om2s = compute_omega_basis(A_sub, n1, 2, ap2s, ap1s) if ap2s else np.zeros((0,0))

        d1s = dim_om(om1s)
        if d1s == 0:
            results.append({'v': v, 'beta1_sub': 0, 'h2rel': 0, 'score': sum(A[v])})
            continue

        bd1s = build_full_boundary_matrix(ap1s, ap0s)
        rk1s = np.linalg.matrix_rank(bd1s @ om1s, tol=1e-8)
        z1s = d1s - rk1s

        d2s = dim_om(om2s)
        if d2s > 0:
            bd2s = build_full_boundary_matrix(ap2s, ap1s)
            bd2oms = np.linalg.lstsq(om1s, bd2s @ om2s, rcond=None)[0]
            b1s = np.linalg.matrix_rank(bd2oms, tol=1e-8)
        else:
            b1s = 0
        beta1_sub = z1s - b1s

        if beta1_sub == 0:
            results.append({'v': v, 'beta1_sub': 0, 'h2rel': 0, 'score': sum(A[v])})
            continue

        # Check if i*: H₁(T\v) → H₁(T) kills the H₁
        # Find H₁ generator of T\v
        U, S, Vt = np.linalg.svd(bd1s @ om1s, full_matrices=True)
        rk = sum(s > 1e-8 for s in S)
        z1_basis = Vt[rk:, :]

        if b1s > 0:
            b1_in_z1 = z1_basis @ bd2oms
            U_b, S_b, _ = np.linalg.svd(b1_in_z1, full_matrices=True)
            b_rk = sum(s > 1e-8 for s in S_b)
            h1_gen = U_b[:, b_rk:][:, 0]
        else:
            h1_gen = z1_basis[0] if z1_basis.shape[0] > 0 else np.zeros(d1s)

        # Embed into T
        h1_A1_sub = om1s @ (z1_basis.T @ h1_gen)
        ap1_T_list = [tuple(p) for p in ap[1]]
        h1_T = np.zeros(len(ap1_T_list))
        remap = {i: others[i] for i in range(n1)}
        for i, p in enumerate(ap1s):
            path_T = (remap[p[0]], remap[p[1]])
            if path_T in ap1_T_list:
                h1_T[ap1_T_list.index(path_T)] = h1_A1_sub[i]

        # Is h1_T a boundary in T?
        d2_T = dim_om(om[2])
        if d2_T > 0:
            bd2_T = build_full_boundary_matrix(ap[2], ap[1])
            h1_om = np.linalg.lstsq(om[1], h1_T, rcond=None)[0]
            bd2om_T = np.linalg.lstsq(om[1], bd2_T @ om[2], rcond=None)[0]
            preimage = np.linalg.lstsq(bd2om_T, h1_om, rcond=None)[0]
            resid = np.linalg.norm(bd2om_T @ preimage - h1_om)
            h2rel = 1 if resid < 1e-6 else 0
        else:
            h2rel = 0

        score = sum(A[v][j] for j in range(n) if j != v)
        results.append({'v': v, 'beta1_sub': beta1_sub, 'h2rel': h2rel, 'score': score})

    sigma_beta1 = sum(r['beta1_sub'] for r in results)
    sigma_h2rel = sum(r['h2rel'] for r in results)

    return {
        'beta1_T': beta1_T, 'z2_dim': z2_dim,
        'sigma_beta1': sigma_beta1, 'sigma_h2rel': sigma_h2rel,
        'vertices': results,
        'd2': d2,
    }


# ============================================================
# Main analysis
# ============================================================
print("=" * 70)
print("Σ h₂_rel vs Σ β₁(T\\v) ANALYSIS")
print("=" * 70)

for n in [5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    t0 = time.time()

    joint_dist = Counter()  # (sigma_beta1, sigma_h2rel) -> count
    z2_by_sigma = Counter()  # (sigma_h2rel, z2_dim) -> count
    beta1_dist = Counter()  # beta1_T -> count
    d2_by_sigma = Counter()  # (sigma_h2rel, d2) -> count

    for bits in range(total):
        A = build_adj(n, bits)
        info = analyze_tournament(A, n)

        joint_dist[(info['sigma_beta1'], info['sigma_h2rel'])] += 1
        z2_by_sigma[(info['sigma_h2rel'], info['z2_dim'])] += 1
        beta1_dist[info['beta1_T']] += 1
        d2_by_sigma[(info['sigma_h2rel'], info['d2'])] += 1

        if (bits+1) % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {bits+1}/{total} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")

    print(f"\n  β₁(T) distribution: {dict(sorted(beta1_dist.items()))}")

    print(f"\n  Joint (Σβ₁(T\\v), Σh₂_rel):")
    for (sb, sh), cnt in sorted(joint_dist.items()):
        rate = sh/sb if sb > 0 else 0
        print(f"    Σβ₁={sb}, Σh₂_rel={sh}: {cnt} tournaments (rate={rate:.2f})")

    print(f"\n  (Σh₂_rel, dim Z₂):")
    for (sh, z2), cnt in sorted(z2_by_sigma.items()):
        print(f"    Σh₂_rel={sh}, z2_dim={z2}: {cnt}")

    print(f"\n  (Σh₂_rel, dim Ω₂):")
    for (sh, d2), cnt in sorted(d2_by_sigma.items()):
        if cnt >= 10:  # Only show significant entries
            print(f"    Σh₂_rel={sh}, d₂={d2}: {cnt}")

print("\nDone.")
