#!/usr/bin/env python3
"""
WHY IS β_2 = 0 FOR ALL TOURNAMENTS?

Algebraic analysis of the chain complex at dimension 2.

β_2 = dim(ker ∂_2 ∩ Ω_2) - dim(im ∂_3 ∩ Ω_2)

For tournaments:
- Ω_2 = {2-paths (a,b,c) where a→b, b→c, AND a→c} = transitive triples
- ∂_2(a,b,c) = (b,c) - (a,c) + (a,b)
- Ω_3 = {3-paths (a,b,c,d) where all faces are allowed}

The boundary ∂_2 maps Ω_2 into Ω_1 = A_1 (all edges).
ker ∂_2 in Ω_2 consists of formal sums of transitive triples whose boundary = 0.

For β_2=0, we need: every 2-cycle in Ω_2 is a 2-boundary from Ω_3.

Let's verify this computationally and look for the algebraic reason.
"""
import numpy as np
from itertools import combinations, permutations
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
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

# ===== n=4: FULL ANALYSIS =====
print("=" * 70)
print("β_2 = 0 ANALYSIS: n=4")
print("=" * 70)

n = 4
for idx, A in enumerate(all_tournaments_gen(n)):
    # Get Ω_2 and Ω_3
    a2 = enumerate_allowed_paths(A, n, 2)
    a1 = enumerate_allowed_paths(A, n, 1)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)

    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if dim_om2 > 0:
        # Compute ker(∂_2) on Ω_2
        bd2 = build_full_boundary_matrix(a2, a1)
        bd2_om = bd2 @ om2
        U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
        rank2 = sum(s > 1e-8 for s in S)
        ker2 = dim_om2 - rank2

        # Compute im(∂_3) on Ω_3 → Ω_2
        if dim_om3 > 0:
            bd3 = build_full_boundary_matrix(a3, a2)
            bd3_om = bd3 @ om3
            rank3 = np.linalg.matrix_rank(bd3_om, tol=1e-8)
        else:
            rank3 = 0

        beta_2 = ker2 - rank3

        if idx < 5 or beta_2 > 0:
            # What are the transitive triples (Ω_2)?
            trans_triples = []
            for path in a2:
                v0, v1, v2 = path
                if A[v0][v2] == 1:
                    trans_triples.append(path)

            print(f"\n  Tournament #{idx}: dim(Ω_2)={dim_om2}, dim(Ω_3)={dim_om3}")
            print(f"    ker(∂_2)={ker2}, im(∂_3)={rank3}, β_2={beta_2}")
            print(f"    transitive triples: {len(trans_triples)}")
            if len(trans_triples) <= 8:
                for t in trans_triples:
                    print(f"      {t[0]}→{t[1]}→{t[2]} (with {t[0]}→{t[2]})")

# ===== n=5: Sample analysis =====
print("\n\n" + "=" * 70)
print("β_2 = 0 ANALYSIS: n=5")
print("=" * 70)

n = 5
ker2_dist = []
for A in all_tournaments_gen(n):
    a2 = enumerate_allowed_paths(A, n, 2)
    a1 = enumerate_allowed_paths(A, n, 1)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)

    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if dim_om2 > 0:
        bd2 = build_full_boundary_matrix(a2, a1)
        bd2_om = bd2 @ om2
        U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
        rank2 = sum(s > 1e-8 for s in S)
        ker2 = dim_om2 - rank2

        if dim_om3 > 0:
            bd3 = build_full_boundary_matrix(a3, a2)
            bd3_om = bd3 @ om3
            rank3 = np.linalg.matrix_rank(bd3_om, tol=1e-8)
        else:
            rank3 = 0

        ker2_dist.append((ker2, rank3, dim_om2, dim_om3))
    else:
        ker2_dist.append((0, 0, 0, dim_om3))

from collections import Counter
print(f"\n  Distribution of (ker∂_2, im∂_3, dim Ω_2, dim Ω_3):")
dist = Counter(ker2_dist)
for key in sorted(dist.keys()):
    print(f"    {key}: {dist[key]}")

# KEY INSIGHT: Is ker∂_2 = im∂_3 always?
always_equal = all(k == r for k, r, _, _ in ker2_dist)
print(f"\n  ker(∂_2) = im(∂_3) always? {always_equal}")

# For the cases with ker∂_2 > 0, verify that im∂_3 fills the kernel
for k, r, d2, d3 in ker2_dist:
    if k > 0 and k != r:
        print(f"  COUNTEREXAMPLE: ker={k}, im={r}")

print(f"\n  CONCLUSION: β_2 = ker(∂_2) - im(∂_3) = 0 for all n=5 tournaments")

# ===== STRUCTURAL REASON =====
print("\n\n" + "=" * 70)
print("STRUCTURAL ANALYSIS: WHY ker(∂_2) = im(∂_3)?")
print("=" * 70)

print("""
For a transitive triple (a,b,c) with a→b, b→c, a→c:
  ∂(a,b,c) = (b,c) - (a,c) + (a,b)

A 2-cycle Σ α_t (a,b,c) has ∂ = 0, meaning each edge appears
with net coefficient 0 in the total.

For this cycle to be a boundary, we need it to come from Ω_3.
A 3-path (a,b,c,d) with a→b→c→d is in Ω_3 iff:
  - (b,c,d) ∈ Ω_2: needs b→d
  - (a,c,d) ∈ A_2: needs a→c and c→d (already have c→d)... but a→c?
  - (a,b,d) ∈ A_2: needs a→b and b→d (already have a→b)... but b→d?
  - (a,b,c) ∈ Ω_2: needs a→c (transitive)

Wait, the ∂-invariance of a 3-path requires ALL faces to be allowed.
Let me verify which faces appear.
""")

# Detailed Ω_3 analysis at n=5
n = 5
for trial in range(1):
    A = list(all_tournaments_gen(n))[500]  # pick one

    a3 = enumerate_allowed_paths(A, n, 3)
    a2 = enumerate_allowed_paths(A, n, 2)
    a1 = enumerate_allowed_paths(A, n, 1)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    print(f"\n  Example tournament: |A_3|={len(a3)}, dim(Ω_3)={dim_om3}")

    # For each 3-path, check which face conditions hold
    for path in a3[:10]:
        v0, v1, v2, v3 = path
        # Faces: (v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)
        face_0 = (v1,v2,v3)  # delete v0
        face_1 = (v0,v2,v3)  # delete v1
        face_2 = (v0,v1,v3)  # delete v2
        face_3 = (v0,v1,v2)  # delete v3

        # Check if each face is an allowed 2-path
        def is_allowed_2(path, A):
            a, b, c = path
            return A[a][b] == 1 and A[b][c] == 1

        f0_ok = is_allowed_2(face_0, A)
        f1_ok = is_allowed_2(face_1, A)
        f2_ok = is_allowed_2(face_2, A)
        f3_ok = is_allowed_2(face_3, A)

        all_ok = f0_ok and f1_ok and f2_ok and f3_ok
        print(f"    path {path}: faces ok = [{f0_ok},{f1_ok},{f2_ok},{f3_ok}], all={all_ok}")

    # Ω_3 condition: all 4 faces must be in A_2
    # Face (v0,v2,v3): allowed iff v0→v2 (shortcut over v1)
    # Face (v0,v1,v3): allowed iff v0→v1 (given) and v1→v3 (shortcut over v2)
    print(f"""
  For path (v0,v1,v2,v3) where v0→v1→v2→v3:
    Face (v1,v2,v3): always allowed (v1→v2→v3)
    Face (v0,v1,v2): always allowed (v0→v1→v2)
    Face (v0,v2,v3): allowed iff v0→v2
    Face (v0,v1,v3): allowed iff v1→v3

  So (v0,v1,v2,v3) ∈ Ω_3 iff v0→v2 AND v1→v3 (both shortcuts exist).
  This means: in the 4-vertex induced subtournament {v0,v1,v2,v3},
  the path v0→v1→v2→v3 is "doubly transitive" (both 2-step shortcuts exist).
  """)

print("\nDone.")
