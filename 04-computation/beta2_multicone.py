#!/usr/bin/env python3
"""
beta2_multicone.py — Multi-vertex cone covering of ker(d₂)

For each tournament T, find a set S of vertices such that the cones
c_v for v ∈ S together cover all of ker(d₂).

Concretely: for each z ∈ Z₂ basis, try to decompose z into parts
that can be filled by different cones.

The idea: ∂₃(Σ_v c_v(z_v)) = z if we can decompose z = Σ_v z_v
where each c_v(z_v) ∈ Ω₃ and ∂₃(c_v(z_v)) = z_v.

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


def compute_cone_image(A, n, v, ap2, ap3_list):
    """Compute the image of cone c_v on A₂ paths.
    Returns a matrix C where C[i,j] = coefficient of ap3[i] in c_v(ap2[j]).
    """
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
print("MULTI-VERTEX CONE COVERING OF ker(d₂)")
print("=" * 70)

for n in [5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    if n >= 7:
        break

    t0 = time.time()

    # Statistics
    min_vertices_needed = Counter()
    total_count = 0
    failures = []

    for bits in range(total):
        A = build_adj(n, bits)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

        d2 = dim_om(om2)
        if d2 == 0:
            min_vertices_needed[0] += 1
            total_count += 1
            continue

        # Compute Z₂ = ker(∂₂|Ω₂)
        bd2 = build_full_boundary_matrix(ap2, ap1)
        d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
        rk = sum(s > 1e-8 for s in S)
        z2_dim = d2 - rk

        if z2_dim == 0:
            min_vertices_needed[0] += 1
            total_count += 1
            continue

        # Z₂ basis in A₂ coordinates
        z2_om = Vt[rk:, :]  # z2_dim × d2
        z2_A2 = om2 @ z2_om.T  # |A₂| × z2_dim

        ap3_list = [tuple(p) for p in ap3]
        bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

        # For each vertex, compute the cone filling matrix
        # cone_fill_v: maps Z₂ basis vector z to ∂₃(c_v(z)) in A₂ coords
        cone_fill_matrices = {}
        for v in range(n):
            C_v = compute_cone_image(A, n, v, ap2, ap3_list)
            # c_v maps A₂ → A₃ (subset), then ∂₃ maps back to A₂
            # For z in Z₂ (A₂ coords): c_v(z) = C_v @ z (in A₃), ∂₃(c_v(z)) = bd3 @ C_v @ z
            fill_v = bd3 @ C_v @ z2_A2  # |A₂| × z2_dim
            cone_fill_matrices[v] = fill_v

        # Strategy: try all subsets of vertices, find minimum that spans Z₂
        # For efficiency, represent in Z₂ coordinates: project fill_v into Z₂ basis
        # We want: Σ_v α_v * fill_v_projected = I (on Z₂ basis)
        # Actually: each fill_v gives a map Z₂ → A₂, and we want the SUM to be identity on Z₂.
        # Project to Z₂: z2_A2.T @ fill_v should give z2_dim × z2_dim matrix

        # Better: we need the union of column spaces of projected fill matrices to span Z₂
        # Stack all projected fills and find rank

        # For each v: projected_v = z2_A2.T @ fill_v (z2_dim × z2_dim), but this might not be right
        # Actually z2_A2.T isn't an orthogonal projector in general.

        # Use pseudoinverse: for z in Z₂ (A₂ coords), project fill_v(z) back to Z₂
        # z2_A2 has shape |A₂| × z2_dim, so z2_A2^+ = (z2_A2^T z2_A2)^{-1} z2_A2^T
        z2_pinv = np.linalg.pinv(z2_A2)  # z2_dim × |A₂|

        # For each v: proj_v = z2_pinv @ fill_v = z2_dim × z2_dim
        projected = {}
        for v in range(n):
            projected[v] = z2_pinv @ cone_fill_matrices[v]  # z2_dim × z2_dim

        # Now: we want to find S ⊆ [n] such that Σ_{v∈S} α_v proj_v = I
        # This is asking whether I is in the span of {proj_v : v ∈ [n]} as matrices.

        # Simpler test: does the joint image of all proj_v span Z₂?
        # Stack columns: [proj_0 | proj_1 | ... | proj_{n-1}] and check rank
        all_cols = np.hstack([projected[v] for v in range(n)])  # z2_dim × (n * z2_dim)
        joint_rank = np.linalg.matrix_rank(all_cols, tol=1e-8)

        if joint_rank < z2_dim:
            failures.append(bits)
            min_vertices_needed[-1] += 1
        else:
            # Find minimum subset
            found = False
            for k in range(1, n+1):
                if found:
                    break
                # Try subsets of size k
                from itertools import combinations
                for S in combinations(range(n), k):
                    cols = np.hstack([projected[v] for v in S])
                    if np.linalg.matrix_rank(cols, tol=1e-8) >= z2_dim:
                        min_vertices_needed[k] += 1
                        found = True
                        break

        total_count += 1
        if total_count % 5000 == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {total_count}/{total} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")
    print(f"  Minimum vertices needed distribution: {dict(sorted(min_vertices_needed.items()))}")
    if failures:
        print(f"  FAILURES (joint cone doesn't span Z₂): {len(failures)}")
        for b in failures[:5]:
            scores = [sum(build_adj(n, b)[i][j] for j in range(n) if j!=i) for i in range(n)]
            print(f"    bits={b}, scores={sorted(scores)}")
    else:
        print(f"  ALL tournaments: multi-cone covers Z₂!")

print("\n" + "=" * 70)
print("CONE RESIDUAL ANALYSIS")
print("=" * 70)

# For n=5: study what happens when single cone fails
n = 5
m = n*(n-1)//2
total = 1 << m

print(f"\nn={n}: Analyzing cone residuals for single-cone failures")
residual_patterns = Counter()
count = 0

for bits in range(total):
    A = build_adj(n, bits)

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    if d2 == 0:
        continue

    # Z₂
    bd2 = build_full_boundary_matrix(ap2, ap1)
    d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
    U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = d2 - rk

    if z2_dim == 0:
        continue

    z2_om = Vt[rk:, :]
    z2_A2 = om2 @ z2_om.T

    ap3_list = [tuple(p) for p in ap3]
    bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2), 0))

    # For each vertex: what fraction of Z₂ does its cone cover?
    z2_pinv = np.linalg.pinv(z2_A2)

    cone_ranks = []
    for v in range(n):
        C_v = compute_cone_image(A, n, v, ap2, ap3_list)
        fill_v = bd3 @ C_v @ z2_A2
        proj_v = z2_pinv @ fill_v
        cr = np.linalg.matrix_rank(proj_v, tol=1e-8)
        cone_ranks.append(cr)

    residual_patterns[tuple(sorted(cone_ranks, reverse=True))] += 1
    count += 1

print(f"\n  Total tournaments with Z₂ > 0: {count}")
print(f"  Cone rank patterns (sorted desc):")
for pattern, cnt in sorted(residual_patterns.items(), key=lambda x: -x[1])[:15]:
    print(f"    {pattern}: {cnt}")

print("\nDone.")
