#!/usr/bin/env python3
"""
beta2_cone_allvertex.py — Test cone operator for every vertex

For each tournament T and each vertex v, test if the cone c_v gives
a chain contraction at degree 2: ∂₃∘c_v + c_v∘∂₂ = id on Z₂.

We don't need it to work on all of Ω₂, just on Z₂ = ker(∂₂).
Because β₂ = Z₂/B₂, if we can show c_v maps Z₂ into Ω₃ and
∂₃(c_v(z)) = z for z ∈ Z₂, then Z₂ ⊆ B₂ so β₂ = 0.

Actually: for z ∈ Z₂, ∂₂(z) = 0, so c_v(∂₂(z)) = 0.
So we need ∂₃(c_v(z)) = z for z ∈ Z₂.

This is just: z ∈ im(∂₃), i.e., z ∈ B₂. Which is exactly β₂ = 0.

OK so the cone approach is: for each z ∈ Z₂, find σ ∈ Ω₃ with ∂₃σ = z.
The cone c_v(z) is a candidate for σ.

Test: for each (T, v), does c_v: Z₂ → Ω₃ satisfy ∂₃∘c_v = id|_{Z₂}?

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


def test_cone_on_z2(A, n, v, ap2, ap3, ap1, ap0, om2, om3, om1):
    """Test if cone c_v maps Z₂ into Ω₃ with ∂₃∘c_v = id on Z₂."""
    d2 = dim_om(om2)
    d3 = dim_om(om3)
    d1 = dim_om(om1)

    if d2 == 0:
        return True, 0  # Z₂ = 0, trivially true

    # Compute Z₂ = ker(∂₂|_{Ω₂})
    bd2 = build_full_boundary_matrix(ap2, ap1)
    d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]  # d1 × d2
    U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z2_dim = d2 - rk

    if z2_dim == 0:
        return True, 0  # No cycles, trivially true

    # Z₂ basis in Ω₂ coordinates
    z2_basis = Vt[rk:, :]  # z2_dim × d2

    # For each Z₂ basis vector: compute c_v(z) and check ∂₃(c_v(z)) = z
    ap2_list = [tuple(p) for p in ap2]
    ap3_list = [tuple(p) for p in ap3]
    bd3 = build_full_boundary_matrix(ap3, ap2) if ap3 else np.zeros((len(ap2_list), 0))

    success = True
    for j in range(z2_dim):
        # z in A₂ coordinates
        z = om2 @ z2_basis[j, :]

        # c_v(z) in A₃ coordinates
        cv_z = np.zeros(len(ap3_list))
        for i, path in enumerate(ap2):
            if abs(z[i]) < 1e-10:
                continue
            if v in path:
                continue  # Can't cone if v already in path
            if not A[v][path[0]]:
                continue  # v doesn't point to first vertex
            cone_path = tuple([v] + list(path))
            if cone_path in ap3_list:
                cv_z[ap3_list.index(cone_path)] += z[i]

        # ∂₃(c_v(z)) in A₂ coordinates
        d3_cv = bd3 @ cv_z

        # Check if ∂₃(c_v(z)) = z (in A₂)
        diff = d3_cv - z
        diff_norm = np.linalg.norm(diff)

        if diff_norm > 1e-6:
            success = False

    return success, z2_dim


print("=" * 70)
print("CONE OPERATOR ON Z₂ — TESTING ALL (T, v) PAIRS")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For each tournament: does SOME vertex give a working cone?
some_works = 0
none_works = 0
vertex_works_count = Counter()  # How many vertices work per tournament

t0 = time.time()
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    works_count = 0
    for v in range(n):
        ok, z2_dim = test_cone_on_z2(A, n, v, ap2, ap3, ap1, ap0, om2, om3, om1)
        if ok:
            works_count += 1

    vertex_works_count[works_count] += 1
    if works_count > 0:
        some_works += 1
    else:
        none_works += 1

elapsed = time.time() - t0
print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")
print(f"  SOME vertex gives working cone: {some_works}")
print(f"  NO vertex gives working cone: {none_works}")
print(f"  #working vertices distribution: {dict(sorted(vertex_works_count.items()))}")

# For tournaments where no vertex works: examine them
if none_works > 0:
    print(f"\n--- Examining tournaments with no working cone ---")
    count = 0
    for bits in range(total):
        A = build_adj(n, bits)
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

        any_works = False
        for v in range(n):
            ok, z2_dim = test_cone_on_z2(A, n, v, ap2, ap3, ap1, ap0, om2, om3, om1)
            if ok:
                any_works = True
                break

        if not any_works:
            count += 1
            if count <= 3:
                print(f"  bits={bits}, scores={sorted(scores)}")
                d2 = dim_om(om2); d3 = dim_om(om3)
                d2_mat = np.linalg.lstsq(om1, build_full_boundary_matrix(ap2, ap1) @ om2, rcond=None)[0]
                rk = np.linalg.matrix_rank(d2_mat, tol=1e-8)
                z2 = d2 - rk
                print(f"    dim(Ω₂)={d2}, dim(Z₂)={z2}, dim(Ω₃)={d3}")

                # Check β₂
                if d3 > 0:
                    bd3 = build_full_boundary_matrix(ap3, ap2)
                    d3_mat = np.linalg.lstsq(om2, bd3 @ om3, rcond=None)[0]
                    b2 = np.linalg.matrix_rank(d3_mat, tol=1e-8)
                else:
                    b2 = 0
                print(f"    dim(B₂)={b2}, β₂={z2-b2}")

print("\nDone.")
