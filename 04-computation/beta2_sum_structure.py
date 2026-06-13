#!/usr/bin/env python3
"""
beta2_sum_structure.py — Algebraic structure of Σ h₂_rel

Goal: understand WHY Σ_v h₂(T,T\v) ≤ 3 for all tournaments.

Key insight to test: the Euler characteristic approach.
χ_rel(T,T\v) = Σ (-1)^p dim(Ω_p^rel) = Σ (-1)^p β_p^rel

So Σ_v h₂_rel = Σ_v β₂^rel relates to Σ_v χ_rel + Σ_v (other β_p^rel terms).

Also test: does Σ h₂_rel relate to β₁(T)?
From HYP-263: β₁(T) > 0 ⟹ Σ = 0. So the sum can only be large when β₁(T) = 0.

Deeper: compute dim(Ω_p^rel) for each v and each p, to see what constrains h₂_rel.

Author: opus-2026-03-08-S49
"""
import sys, time, random
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


def compute_rel_dimensions(A, n, v):
    """Compute dimensions of relative chain complex Ω_p(T)/Ω_p(T\\v) for p=0,1,2,3."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    dims = {}
    for p in range(4):
        ap_T = enumerate_allowed_paths(A, n, p)
        if p > 0:
            apm1_T = enumerate_allowed_paths(A, n, p-1)
            om_T = compute_omega_basis(A, n, p, ap_T, apm1_T) if ap_T else np.zeros((0,0))
        else:
            om_T = np.eye(n) if ap_T else np.zeros((0,0))

        ap_sub = enumerate_allowed_paths(A_sub, n1, p)
        if p > 0:
            apm1_sub = enumerate_allowed_paths(A_sub, n1, p-1)
            om_sub = compute_omega_basis(A_sub, n1, p, ap_sub, apm1_sub) if ap_sub else np.zeros((0,0))
        else:
            om_sub = np.eye(n1) if ap_sub else np.zeros((0,0))

        d_T = dim_om(om_T) if p > 0 else (n if ap_T else 0)
        d_sub = dim_om(om_sub) if p > 0 else (n1 if ap_sub else 0)

        # Relative dimension = d_T - rk(inclusion)
        # For p=0: inclusion maps n-1 vertices into n vertices, rk = n-1
        if p == 0:
            dims[p] = 1  # Always: vertex v contributes 1 to Ω₀^rel
        else:
            if d_sub == 0:
                dims[p] = d_T
            elif d_T == 0:
                dims[p] = 0
            else:
                # Compute inclusion rank
                ap_T_list = [tuple(pp) for pp in ap_T]
                embed = np.zeros((len(ap_T_list), d_sub))
                for j in range(d_sub):
                    for k, path_sub in enumerate(ap_sub):
                        path_T = tuple(remap[x] for x in path_sub)
                        if path_T in ap_T_list:
                            embed[ap_T_list.index(path_T), j] = om_sub[k, j]
                phi = np.linalg.lstsq(om_T, embed, rcond=None)[0]
                rk = np.linalg.matrix_rank(phi, tol=1e-8)
                dims[p] = d_T - rk

    return dims


def compute_h2_rel_fast(A, n, v):
    """Compute h₂(T,T\\v) using the established method."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    ap2_T = enumerate_allowed_paths(A, n, 2)
    if not ap2_T:
        return 0
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap0_T = enumerate_allowed_paths(A, n, 0)
    ap3_T = enumerate_allowed_paths(A, n, 3)

    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))

    d2_T = dim_om(om2_T)
    if d2_T == 0:
        return 0

    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)

    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)

    d2_sub = dim_om(om2_sub)
    d1_T = dim_om(om1_T)
    d1_sub = dim_om(om1_sub)
    d3_T = dim_om(om3_T)

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
        return 0

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
    if z2_rel == 0:
        return 0

    if d3_T > 0:
        ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)
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
print("ALGEBRAIC STRUCTURE OF Σ h₂_rel")
print("=" * 70)

# Part 1: Euler characteristic of relative complex
print("\n--- Part 1: Relative Euler characteristics ---")
n = 5
m = n*(n-1)//2
total = 1 << m

chi_vs_h2_sum = defaultdict(list)
dim_patterns = Counter()

t0 = time.time()
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    h2_sum = 0
    chi_sum = 0
    for v in range(n):
        h2r = compute_h2_rel_fast(A, n, v)
        h2_sum += h2r

        dims = compute_rel_dimensions(A, n, v)
        # χ_rel = dim Ω₀^rel - dim Ω₁^rel + dim Ω₂^rel - dim Ω₃^rel
        chi_rel = dims[0] - dims[1] + dims[2] - dims[3]
        chi_sum += chi_rel

        dim_patterns[(dims[0], dims[1], dims[2], dims[3], h2r)] += 1

    chi_vs_h2_sum[h2_sum].append(chi_sum)

print(f"n={n}: {time.time()-t0:.0f}s")
print(f"\nΣ h₂_rel → Σ χ_rel mapping:")
for h2s in sorted(chi_vs_h2_sum):
    chi_vals = chi_vs_h2_sum[h2s]
    print(f"  Σ h₂_rel = {h2s}: Σ χ_rel in {set(chi_vals)} (count={len(chi_vals)})")

print(f"\nDimension patterns (d0_rel, d1_rel, d2_rel, d3_rel, h2_rel) → count:")
for pat, cnt in sorted(dim_patterns.items(), key=lambda x: -x[1])[:20]:
    print(f"  {pat} → {cnt}")


# Part 2: What is the Ω₂^rel structure when h₂_rel = 1?
# Each path in Ω₂^rel involves v. How many such paths are there?
print(f"\n--- Part 2: Counting v-paths in Ω₂ ---")

v_path_count_dist = defaultdict(Counter)
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    for v in range(n):
        if scores[v] == 0 or scores[v] == n-1:
            continue  # skip source/sink

        # Count allowed 2-paths through v
        ap2 = enumerate_allowed_paths(A, n, 2)
        v_paths = [p for p in ap2 if v in p]

        h2r = compute_h2_rel_fast(A, n, v)
        v_path_count_dist[h2r][len(v_paths)] += 1

print(f"Interior vertices: #(v-paths in A₂) distribution by h₂_rel:")
for h2r in sorted(v_path_count_dist):
    print(f"  h₂_rel = {h2r}: {dict(sorted(v_path_count_dist[h2r].items()))}")


# Part 3: Relationship between h₂_rel and 3-cycles through v
print(f"\n--- Part 3: h₂_rel vs 3-cycles through v ---")

cycles_dist = defaultdict(Counter)
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    for v in range(n):
        if scores[v] == 0 or scores[v] == n-1:
            continue

        # Count 3-cycles through v
        cycle_count = 0
        for i in range(n):
            if i == v:
                continue
            for j in range(n):
                if j == v or j == i:
                    continue
                if A[v][i] and A[i][j] and A[j][v]:
                    cycle_count += 1
        cycle_count //= 2  # each cycle counted twice

        h2r = compute_h2_rel_fast(A, n, v)
        cycles_dist[h2r][cycle_count] += 1

print(f"Interior vertices: #(3-cycles through v) by h₂_rel:")
for h2r in sorted(cycles_dist):
    print(f"  h₂_rel = {h2r}: {dict(sorted(cycles_dist[h2r].items()))}")


# Part 4: The key question — can we express Σ h₂_rel as a function of T?
print(f"\n--- Part 4: Σ h₂_rel as function of tournament invariants ---")

sum_vs_t3 = defaultdict(Counter)
sum_vs_beta1 = defaultdict(Counter)

for bits in range(total):
    A = build_adj(n, bits)

    # Count 3-cycles
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    t3 += 1

    # Compute β₁
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d1 = dim_om(om1)
    d2 = dim_om(om2)
    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    z1 = d1 - rk
    if d2 > 0:
        bd2om = np.linalg.lstsq(om1, build_full_boundary_matrix(ap2, ap1) @ om2, rcond=None)[0]
        b1 = np.linalg.matrix_rank(bd2om, tol=1e-8)
    else:
        b1 = 0
    beta1 = z1 - b1

    h2_sum = sum(compute_h2_rel_fast(A, n, v) for v in range(n))

    sum_vs_t3[h2_sum][t3] += 1
    sum_vs_beta1[h2_sum][beta1] += 1

print(f"Σ h₂_rel vs t₃ (3-cycle count):")
for s in sorted(sum_vs_t3):
    print(f"  Σ={s}: t₃ in {dict(sorted(sum_vs_t3[s].items()))}")

print(f"\nΣ h₂_rel vs β₁(T):")
for s in sorted(sum_vs_beta1):
    print(f"  Σ={s}: β₁ in {dict(sorted(sum_vs_beta1[s].items()))}")

print("\nDone.")
