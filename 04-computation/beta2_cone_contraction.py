#!/usr/bin/env python3
"""
beta2_cone_contraction.py — Test cone operator as chain contraction

For a tournament T with source vertex v (v → all others):
Define cone operator c_v: A_p → A_{p+1} by c_v(a_0,...,a_p) = (v, a_0,...,a_p).

Does c_v map Ω_p into Ω_{p+1}?
Does ∂_{p+1} ∘ c_v + c_v ∘ ∂_p = id on Ω_p (up to projection)?

If so, this gives a chain contraction making Ω_* acyclic (all β_p = 0 for p ≥ 1).

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


print("=" * 70)
print("CONE OPERATOR AS CHAIN CONTRACTION")
print("=" * 70)

# Test with n=5 tournaments
n = 5
m = n*(n-1)//2

# Find a tournament with a source vertex
# The transitive tournament (bits=0) has vertex 0 as source
bits = 0
A = build_adj(n, bits)
scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
source = scores.index(n-1)  # vertex with d+ = n-1
print(f"\nTransitive T_{n}: source = {source} (d+ = {scores[source]})")
print(f"Scores: {scores}")

# Enumerate paths
ap = {}
for p in range(5):
    ap[p] = enumerate_allowed_paths(A, n, p)
    print(f"  |A_{p}| = {len(ap[p])}")

# Compute Ω bases
om = {}
for p in range(5):
    if p == 0:
        om[p] = np.eye(n)
    elif ap[p]:
        om[p] = compute_omega_basis(A, n, p, ap[p], ap[p-1])
    else:
        om[p] = np.zeros((0,0))
    print(f"  dim(Ω_{p}) = {dim_om(om[p]) if p > 0 else n}")

# For source vertex v: define cone c_v on allowed paths
v = source
print(f"\n--- Cone operator c_{v} ---")

# c_v maps an allowed p-path (a0,...,ap) to (v, a0,...,ap)
# This is a (p+1)-path. It's ALLOWED if v → a0 (true since v is source).

# First: check c_v maps A_p into A_{p+1}
for p in [0, 1, 2, 3]:
    if not ap[p]:
        continue
    cone_paths = []
    cone_in_ap = 0
    for path in ap[p]:
        cone_path = [v] + list(path)
        # Check if this is an allowed (p+1)-path
        is_allowed = True
        for i in range(len(cone_path)-1):
            if not A[cone_path[i]][cone_path[i+1]]:
                is_allowed = False
                break
        if is_allowed:
            cone_paths.append(tuple(cone_path))
            if any(tuple(q) == tuple(cone_path) for q in ap[p+1]):
                cone_in_ap += 1
        else:
            # v → a0 should always work since v is source
            pass

    print(f"  c_{v}(A_{p}): {len(cone_paths)}/{len(ap[p])} paths allowed, {cone_in_ap} in A_{p+1}")

# Now: does c_v map Ω_p into Ω_{p+1}?
# For each Ω_p basis vector, compute c_v(ω) and check if it's in Ω_{p+1}
for p in [1, 2]:
    dp = dim_om(om[p])
    dp1 = dim_om(om[p+1]) if p+1 < 5 else 0
    if dp == 0:
        continue

    ap_list = [tuple(q) for q in ap[p]]
    ap1_list = [tuple(q) for q in ap[p+1]]

    # For each Ω_p basis vector: compute c_v(ω) as a vector in A_{p+1}
    # Then check if it's in Ω_{p+1}
    cone_vectors = []
    for j in range(dp):
        # ω = Σ_i om[p][i,j] * (path_i)
        # c_v(ω) = Σ_i om[p][i,j] * c_v(path_i)
        cv = np.zeros(len(ap1_list))
        for i, path in enumerate(ap[p]):
            if abs(om[p][i,j]) < 1e-10:
                continue
            cone_path = tuple([v] + list(path))
            if cone_path in ap1_list:
                cv[ap1_list.index(cone_path)] += om[p][i,j]

        cone_vectors.append(cv)

    # Project each cone vector into Ω_{p+1}
    if dp1 > 0:
        cone_mat = np.column_stack(cone_vectors)  # |A_{p+1}| × dp
        # Express in Ω_{p+1} coordinates
        coords = np.linalg.lstsq(om[p+1], cone_mat, rcond=None)[0]
        # Residual
        resid = cone_mat - om[p+1] @ coords
        resid_norm = np.linalg.norm(resid)
        rk_cone = np.linalg.matrix_rank(coords, tol=1e-8)
        print(f"\n  c_{v}(Ω_{p}) into Ω_{p+1}: residual = {resid_norm:.6f}, rank = {rk_cone}/{dp}")
    else:
        print(f"\n  c_{v}(Ω_{p}): Ω_{p+1} = 0, can't map there")

# Now: test ∂ ∘ c + c ∘ ∂ = id on Ω_2
# ∂₃(c_v(ω)) + c_v(∂₂(ω)) = ω for ω ∈ Ω₂
print(f"\n--- Chain contraction identity ∂₃∘c_{v} + c_{v}∘∂₂ = id on Ω₂ ---")

d2 = dim_om(om[2])
d1 = dim_om(om[1])
d3 = dim_om(om[3]) if 3 in om else 0

ap2_list = [tuple(q) for q in ap[2]]
ap3_list = [tuple(q) for q in ap[3]]
ap1_list = [tuple(q) for q in ap[1]]

for j in range(d2):
    # ω = j-th Ω₂ basis vector (in A₂ coordinates)
    omega = om[2][:, j]

    # 1. c_v(ω): cone of ω, as vector in A₃
    cv_omega = np.zeros(len(ap3_list))
    for i, path in enumerate(ap[2]):
        if abs(omega[i]) < 1e-10:
            continue
        cone_path = tuple([v] + list(path))
        if cone_path in ap3_list:
            cv_omega[ap3_list.index(cone_path)] += omega[i]

    # 2. ∂₃(c_v(ω)): boundary of the cone, as vector in A₂
    bd3 = build_full_boundary_matrix(ap[3], ap[2])
    d3_cv = bd3 @ cv_omega  # in A₂ coordinates

    # 3. ∂₂(ω): boundary of ω, as vector in A₁
    bd2 = build_full_boundary_matrix(ap[2], ap[1])
    d2_omega = bd2 @ omega  # in A₁ coordinates

    # 4. c_v(∂₂(ω)): cone of boundary, as vector in A₂
    cv_d2 = np.zeros(len(ap2_list))
    for i, path in enumerate(ap[1]):
        if abs(d2_omega[i]) < 1e-10:
            continue
        cone_path = tuple([v] + list(path))
        if cone_path in ap2_list:
            cv_d2[ap2_list.index(cone_path)] += d2_omega[i]

    # 5. ∂₃∘c_v + c_v∘∂₂ should equal ω (on A₂ level)
    lhs = d3_cv + cv_d2
    diff = lhs - omega
    diff_norm = np.linalg.norm(diff)

    if j < 5 or diff_norm > 1e-6:
        print(f"  j={j}: ||∂₃∘c_{v}(ω) + c_{v}∘∂₂(ω) - ω|| = {diff_norm:.6f}")

# Now test with a NON-transitive tournament
print(f"\n{'='*70}")
print("TESTING ON NON-TRANSITIVE TOURNAMENT")
print("=" * 70)

# Find a regular (or near-regular) tournament
for bits in range(1024):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
    if sorted(scores) == [2,2,2,2,2]:  # regular n=5
        print(f"\nRegular T₅: bits={bits}, scores={scores}")

        # This has NO source vertex. Try cone with highest-degree vertex
        v_max = scores.index(max(scores))
        print(f"  Using v={v_max} (d+={scores[v_max]})")

        # Enumerate paths
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap0 = enumerate_allowed_paths(A, n, 0)
        om2_t = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
        om3_t = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))
        om1_t = compute_omega_basis(A, n, 1, ap1, ap0)

        d2_t = dim_om(om2_t)
        ap2_list = [tuple(q) for q in ap2]
        ap3_list = [tuple(q) for q in ap3]
        ap1_list = [tuple(q) for q in ap1]

        # Test c_v on each Ω₂ basis vector
        for j in range(min(d2_t, 5)):
            omega = om2_t[:, j]
            # Cone
            cv_omega = np.zeros(len(ap3_list))
            for i, path in enumerate(ap2):
                if abs(omega[i]) < 1e-10:
                    continue
                cone_path = tuple([v_max] + list(path))
                if cone_path in ap3_list:
                    cv_omega[ap3_list.index(cone_path)] += omega[i]
                elif A[v_max][path[0]]:
                    # Edge exists but cone path not in A₃?
                    pass

            # ∂₃(c_v(ω))
            bd3 = build_full_boundary_matrix(ap3, ap2)
            d3_cv = bd3 @ cv_omega

            # ∂₂(ω)
            bd2 = build_full_boundary_matrix(ap2, ap1)
            d2_omega = bd2 @ omega

            # c_v(∂₂(ω))
            cv_d2 = np.zeros(len(ap2_list))
            for i, path in enumerate(ap1):
                if abs(d2_omega[i]) < 1e-10:
                    continue
                cone_path = tuple([v_max] + list(path))
                if cone_path in ap2_list:
                    cv_d2[ap2_list.index(cone_path)] += d2_omega[i]

            lhs = d3_cv + cv_d2
            diff = lhs - omega
            diff_norm = np.linalg.norm(diff)
            print(f"  j={j}: ||∂₃∘c_{v_max} + c_{v_max}∘∂₂ - id|| = {diff_norm:.6f}")

        break

print("\nDone.")
