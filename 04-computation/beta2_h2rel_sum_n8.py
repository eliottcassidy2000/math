#!/usr/bin/env python3
"""
beta2_h2rel_sum_n8.py — Test HYP-262 (Σ h₂_rel ≤ 3) at n=8

Sample random tournaments and compute the sum of h₂_rel over all vertices.
If the bound holds at n=8, it strongly supports universality.

Also: track WHICH vertices contribute h₂_rel > 0 (score distribution).

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


def compute_h2_rel(A, n, v):
    """Compute h₂(T,T\\v)."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap3_T = enumerate_allowed_paths(A, n, 3)
    ap0_T = enumerate_allowed_paths(A, n, 0)

    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)
    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)

    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

    d2_T = dim_om(om2_T); d1_T = dim_om(om1_T); d3_T = dim_om(om3_T)
    d1_sub = dim_om(om1_sub); d2_sub = dim_om(om2_sub)

    if d2_T == 0: return 0

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
    if d_rel == 0: return 0

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
    R = np.linalg.svd(psi, full_matrices=True)[0][:, rk_psi:] if rk_psi > 0 else np.eye(d1_T)

    coords2_rel = R.T @ coords2_T @ Q
    rk_d2_rel = np.linalg.matrix_rank(coords2_rel, tol=1e-8)
    z2_rel = d_rel - rk_d2_rel
    if z2_rel == 0: return 0

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
    return z2_rel - rk_d3_rel


print("=" * 70)
print("HYP-262 TEST: Σ h₂_rel at n=8 (sampled)")
print("=" * 70)

random.seed(42)
n = 8
m = n*(n-1)//2  # 28
total = 1 << m

N_SAMPLES = 200
samples = random.sample(range(total), N_SAMPLES)

sum_dist = Counter()
interior_sum_dist = Counter()
max_sum = 0
score_dist_when_nonzero = Counter()  # d+(v) when h2_rel > 0
vertex_nonzero_counts = Counter()  # how many vertices per tournament have h2_rel > 0

t0 = time.time()
for idx, bits in enumerate(samples):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    total_h2 = 0
    interior_h2 = 0
    nonzero_count = 0
    for v in range(n):
        h2r = compute_h2_rel(A, n, v)
        total_h2 += h2r
        if h2r > 0:
            nonzero_count += 1
            score_dist_when_nonzero[scores[v]] += 1
        if 1 <= scores[v] <= n-2:
            interior_h2 += h2r

    sum_dist[total_h2] += 1
    interior_sum_dist[interior_h2] += 1
    vertex_nonzero_counts[nonzero_count] += 1
    max_sum = max(max_sum, total_h2)

    if (idx + 1) % 20 == 0:
        elapsed = time.time() - t0
        print(f"  n={n}: {idx+1}/{N_SAMPLES} ({elapsed:.0f}s) max_sum={max_sum}")

elapsed = time.time() - t0
print(f"\nn={n}: {N_SAMPLES} tournaments in {elapsed:.0f}s")
print(f"  Σ h₂_rel (all v): {dict(sorted(sum_dist.items()))}")
print(f"  Σ h₂_rel (interior v): {dict(sorted(interior_sum_dist.items()))}")
print(f"  Max Σ h₂_rel = {max_sum}")
print(f"  #vertices with h₂_rel>0: {dict(sorted(vertex_nonzero_counts.items()))}")
print(f"  d⁺(v) when h₂_rel>0: {dict(sorted(score_dist_when_nonzero.items()))}")

if max_sum <= 3:
    print(f"\n  ✓ HYP-262 HOLDS at n=8: Σ h₂_rel ≤ 3")
    print(f"  Interior vertices: at least {n-2} = {n-2} > 3, so HYP-258 follows by pigeonhole")
else:
    print(f"\n  ✗ HYP-262 FAILS at n=8: max Σ = {max_sum}")

print("\nDone.")
