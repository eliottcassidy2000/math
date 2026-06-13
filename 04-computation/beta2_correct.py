#!/usr/bin/env python3
"""
β_2 FOR TOURNAMENTS — CORRECT COMPUTATION

Previous computations assumed Ω_2 = span(transitive triples). This is WRONG.
Ω_2 includes linear combinations of non-transitive 2-paths where non-allowed
face terms cancel.

This script computes β_2 using the CORRECT Ω_2 from compute_omega_basis.
"""
import numpy as np
import sys, time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def correct_beta2(A, n):
    """Compute β_2 using correct Ω definitions."""
    a0 = enumerate_allowed_paths(A, n, 0)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    # Ω_2 = ker(proj to non-A_1 faces)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

    if dim_om2 == 0:
        return 0, 0, 0, 0

    # Ω_3 = ker(proj to non-A_2 faces)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # ∂_2: Ω_2 → Ω_1
    # First build ∂_2: A_2 → A_1 (full allowed)
    # Actually, we need ∂_2: Ω_2 → Ω_1 = A_1 (since all 1-paths are ∂-invariant)
    bd2_full = build_full_boundary_matrix(
        [tuple(p) for p in a2],
        [tuple(p) for p in a1]
    )
    # Restrict to Ω_2
    bd2_om = bd2_full @ om2  # len(a1) x dim_om2

    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim2 = dim_om2 - rank2

    if ker_dim2 == 0:
        return 0, dim_om2, dim_om3, 0

    # ∂_3: Ω_3 → Ω_2
    if dim_om3 == 0:
        return ker_dim2, dim_om2, 0, ker_dim2

    bd3_full = build_full_boundary_matrix(
        [tuple(p) for p in a3],
        [tuple(p) for p in a2]
    )
    # Restrict to Ω_3 input
    bd3_om3 = bd3_full @ om3  # len(a2) x dim_om3

    # Now bd3_om3 gives the boundary of Ω_3 elements in A_2 coordinates.
    # We need to express this in Ω_2 coordinates.
    # Since ∂(Ω_3) ⊆ Ω_2, bd3_om3 should be in the column space of om2.
    # Project: solve om2 @ x = bd3_om3 for x
    # x = pinv(om2) @ bd3_om3
    bd3_in_om2, residual, _, _ = np.linalg.lstsq(om2, bd3_om3, rcond=None)

    # Verify projection is exact
    err = np.max(np.abs(om2 @ bd3_in_om2 - bd3_om3))
    if err > 1e-6:
        print(f"  WARNING: projection error = {err:.2e}")

    im_rank = np.linalg.matrix_rank(bd3_in_om2, tol=1e-8)
    beta2 = ker_dim2 - im_rank

    return ker_dim2, dim_om2, dim_om3, beta2

# ===== Exhaustive computation =====
print("=" * 70)
print("CORRECT β_2 COMPUTATION")
print("=" * 70)

for n in [4, 5]:
    t0 = time.time()
    beta2_dist = {}
    count = 0
    for A in all_tournaments_gen(n):
        count += 1
        ker, dom2, dom3, b2 = correct_beta2(A, n)
        key = (ker, dom2, dom3, b2)
        beta2_dist[key] = beta2_dist.get(key, 0) + 1

    t1 = time.time()
    print(f"\nn={n}: {count} tournaments ({t1-t0:.1f}s)")
    print(f"  (ker_dim, dim_Ω_2, dim_Ω_3, β_2): count")
    for key in sorted(beta2_dist.keys()):
        ker, dom2, dom3, b2 = key
        print(f"    ({ker}, {dom2}, {dom3}, β_2={b2}): {beta2_dist[key]}")

    all_zero = all(b2 == 0 for (_, _, _, b2) in beta2_dist.keys())
    print(f"  β_2 = 0 for ALL? {all_zero}")

# ===== n=6 exhaustive =====
print(f"\n--- n=6 exhaustive ---")
n = 6
t0 = time.time()
all_zero = True
count = 0
for A in all_tournaments_gen(n):
    count += 1
    ker, dom2, dom3, b2 = correct_beta2(A, n)
    if b2 != 0:
        all_zero = False
        print(f"  T#{count}: β_2 = {b2} !")
    if count % 5000 == 0:
        print(f"  ... {count}/32768 ({time.time()-t0:.0f}s)", flush=True)

t1 = time.time()
print(f"\nn=6: {count} tournaments ({t1-t0:.1f}s)")
print(f"  β_2 = 0 for ALL? {all_zero}")

# ===== Ω_2 dimensions: TT vs true =====
print(f"\n\n{'='*70}")
print("Ω_2: TRANSITIVE TRIPLES vs TRUE DIMENSION")
print("="*70)

from collections import Counter
n = 5
gap_dist = Counter()
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    tt_count = sum(1 for p in a2 if A[p[0]][p[2]] == 1)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    gap = dim_om2 - tt_count
    gap_dist[gap] += 1

print(f"\nn=5: dim(Ω_2) - |TT| distribution:")
for g in sorted(gap_dist.keys()):
    print(f"  gap={g}: {gap_dist[g]} tournaments")

# What determines the gap?
print(f"\nThe gap comes from non-TT pairs sharing a non-allowed face.")
print(f"Non-allowed face (a,c) appears in (a,b,c) when c→a.")
print(f"If multiple 2-paths share the same non-allowed face,")
print(f"we get cancellation chains.")

# Count: how many non-A_1 face types, and rank
for A in list(all_tournaments_gen(5))[:10]:
    a1 = enumerate_allowed_paths(A, 5, 1)
    a2 = enumerate_allowed_paths(A, 5, 2)
    a1_set = set(tuple(p) for p in a1)

    # Non-A_1 faces
    face_counts = {}
    for p in a2:
        a,b,c = tuple(p)
        faces = [(b,c), (a,c), (a,b)]
        signs = [1, -1, 1]
        for face, sign in zip(faces, signs):
            if face not in a1_set:
                face_counts[face] = face_counts.get(face, 0) + 1

    n_ntt = sum(1 for p in a2 if A[p[0]][p[2]] == 0)
    n_faces = len(face_counts)
    multi = sum(1 for v in face_counts.values() if v > 1)
    dim_om2 = compute_omega_basis(A, 5, 2, a2, a1).shape[1]

    print(f"  |non-TT|={n_ntt}, |non-A_1 faces|={n_faces}, multi={multi}, "
          f"dim_Ω_2={dim_om2}, |TT|={len(a2)-n_ntt}")

print("\nDone.")
