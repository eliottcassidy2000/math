#!/usr/bin/env python3
"""
beta2_twocone_proof.py — Toward proving 2 cones suffice

Hypothesis: if N⁺(v₁) ∪ N⁺(v₂) ∪ {v₁,v₂} = [n], then the two cones
together cover all of Z₂.

Stronger: just need N⁺(v₁) ∪ N⁺(v₂) to contain all "first vertices" of
paths in Z₂.

Test: for each 2-cone-needed tournament, check the out-neighborhood
coverage of the working pairs vs failing pairs.

Also: analyze the chain homotopy identity ∂₃∘c_v + c_v∘∂₂ = id
and what happens with partial cones.

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from itertools import combinations
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

def compute_cone_image(A, n, v, ap2, ap3_list):
    C = np.zeros((len(ap3_list), len(ap2)))
    for j, path in enumerate(ap2):
        if v in path:
            continue
        if not A[v][path[0]]:
            continue
        cone_path = tuple([v] + list(path))
        if cone_path in ap3_list:
            C[ap3_list.index(cone_path), j] = 1.0
    return C


print("=" * 70)
print("TWO-CONE PROOF INGREDIENTS")
print("=" * 70)

# Test 1: Does out-neighborhood coverage determine cone success?
print("\nTest 1: Out-neighborhood coverage vs cone success")

n = 5
m = n*(n-1)//2
total = 1 << m

coverage_works = Counter()  # (coverage_size, works?)
pair_coverage_works = 0
pair_coverage_fails = 0

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

    d2 = dim_om(om2)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
    U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = d2 - rk

    if z2_dim == 0:
        continue

    z2_om = Vt[rk:, :]
    z2_A2 = om2 @ z2_om.T
    z2_pinv = np.linalg.pinv(z2_A2)

    ap3_list = [tuple(p) for p in ap3]
    bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    projected = {}
    for v in range(n):
        C_v = compute_cone_image(A, n, v, ap2, ap3_list)
        fill_v = bd3 @ C_v @ z2_A2
        projected[v] = z2_pinv @ fill_v

    # For each pair: check coverage and cone success
    for v1, v2 in combinations(range(n), 2):
        coverage = set()
        coverage.update(j for j in range(n) if A[v1][j])
        coverage.update(j for j in range(n) if A[v2][j])
        coverage.add(v1)
        coverage.add(v2)
        cov_size = len(coverage)
        full_coverage = (cov_size == n)

        cols = np.hstack([projected[v1], projected[v2]])
        works = np.linalg.matrix_rank(cols, tol=1e-8) >= z2_dim

        coverage_works[(cov_size, works)] += 1
        if full_coverage:
            if works:
                pair_coverage_works += 1
            else:
                pair_coverage_fails += 1

print(f"\n  Coverage → success table:")
for (cov, w), cnt in sorted(coverage_works.items()):
    print(f"    coverage={cov}, works={w}: {cnt}")
print(f"\n  Full coverage (all n vertices covered) → works: {pair_coverage_works}")
print(f"  Full coverage → fails: {pair_coverage_fails}")

# Test 2: The partial cone identity
print(f"\n{'='*70}")
print("Test 2: PARTIAL CONE IDENTITY")
print("=" * 70)
print("""
For z ∈ Z₂, define z_v = restriction of z to paths cone-able by v.
Then ∂₃(c_v(z_v)) = z_v - c_v(∂₂(z_v)).

Since z is a cycle but z_v is NOT (in general), we get:
∂₃(c_v(z_v)) = z_v - c_v(∂₂(z_v))
r_v := z - ∂₃(c_v(z_v)) = z - z_v + c_v(∂₂(z_v)) = z_rest + c_v(∂₂(z_v))

where z_rest = z - z_v = non-cone-able part.

So the residual r_v involves:
1. z_rest: paths in z that v can't cone
2. c_v(∂₂(z_v)): correction from partial boundary

If we apply a second cone c_w to r_v, we need to cover the residual.
""")

# Test this for a specific tournament
bits = 41  # Known Σ=3 tournament
A = build_adj(n, bits)
scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
print(f"Tournament bits={bits}, scores={scores}")

ap0 = enumerate_allowed_paths(A, n, 0)
ap1 = enumerate_allowed_paths(A, n, 1)
ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)
om1 = compute_omega_basis(A, n, 1, ap1, ap0)
om2 = compute_omega_basis(A, n, 2, ap2, ap1)

bd2 = build_full_boundary_matrix(ap2, ap1)
d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
rk = sum(s > 1e-8 for s in S)
z2_om = Vt[rk:, :]
z2_A2 = om2 @ z2_om.T
z2_dim = z2_A2.shape[1]

ap2_list = [tuple(p) for p in ap2]
ap3_list = [tuple(p) for p in ap3]
bd3 = build_full_boundary_matrix(ap3, ap2)

print(f"  z2_dim = {z2_dim}")

# For each z basis vector and each v: decompose into cone-able and rest
for j in range(min(z2_dim, 2)):
    z = z2_A2[:, j]
    print(f"\n  z_{j}:")

    for v in range(n):
        # Identify cone-able paths
        z_v = np.zeros_like(z)
        z_rest = np.zeros_like(z)
        for i, p in enumerate(ap2):
            if abs(z[i]) < 1e-10:
                continue
            if v not in p and A[v][p[0]] and tuple([v]+list(p)) in ap3_list:
                z_v[i] = z[i]
            else:
                z_rest[i] = z[i]

        # c_v(z_v) in A₃ coords
        cv_zv = np.zeros(len(ap3_list))
        for i, p in enumerate(ap2):
            if abs(z_v[i]) < 1e-10:
                continue
            cone_path = tuple([v] + list(p))
            idx = ap3_list.index(cone_path)
            cv_zv[idx] = z_v[i]

        # ∂₃(c_v(z_v)) in A₂ coords
        d3_cv = bd3 @ cv_zv

        # Residual
        residual = z - d3_cv
        res_norm = np.linalg.norm(residual)

        # Also: c_v(∂₂(z_v)) — partial boundary correction
        d2_zv = bd2 @ z_v  # in A₁ coords
        cv_d2zv = np.zeros(len(ap2_list))
        for i, p in enumerate(ap1):
            if abs(d2_zv[i]) < 1e-10:
                continue
            cone_path = tuple([v] + list(p))
            if cone_path in ap2_list:
                cv_d2zv[ap2_list.index(cone_path)] += d2_zv[i]

        # z_rest + c_v(∂₂(z_v)) should equal residual
        expected_res = z_rest + cv_d2zv
        check = np.linalg.norm(residual - expected_res)

        n_coneable = sum(1 for i in range(len(ap2)) if abs(z_v[i]) > 1e-10)
        n_rest = sum(1 for i in range(len(ap2)) if abs(z_rest[i]) > 1e-10)

        print(f"    v={v}: coneable={n_coneable}, rest={n_rest}, "
              f"||residual||={res_norm:.4f}, identity_check={check:.2e}")

# Test 3: Does coverage of FIRST VERTICES of Z₂ suffice?
print(f"\n{'='*70}")
print("Test 3: FIRST VERTEX COVERAGE")
print("=" * 70)

# For each Z₂ basis vector, which vertices appear as first vertex?
first_vertices_needed = Counter()

for bits in range(total):
    A = build_adj(n, bits)

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

    d2 = dim_om(om2)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
    U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = d2 - rk

    if z2_dim == 0:
        continue

    z2_A2 = om2 @ Vt[rk:, :].T

    # Collect first vertices of ALL nonzero paths in Z₂
    first_verts = set()
    for j in range(z2_dim):
        z = z2_A2[:, j]
        for i, p in enumerate(ap2):
            if abs(z[i]) > 1e-8:
                first_verts.add(p[0])

    first_vertices_needed[len(first_verts)] += 1

print(f"\nn={n}: Number of distinct first vertices in Z₂:")
for k, cnt in sorted(first_vertices_needed.items()):
    print(f"  {k} first vertices: {cnt} tournaments")

print("\nDone.")
