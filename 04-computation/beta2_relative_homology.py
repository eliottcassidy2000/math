#!/usr/bin/env python3
"""Relative homology H_2(T, T\v) for vertex deletion.

Long exact sequence: H_2(T\v) → H_2(T) → H_2(T, T\v) → H_1(T\v) → H_1(T)

Since H_2(T\v) = 0 by induction (assuming β_2=0 for smaller tournaments),
H_2(T) injects into H_2(T, T\v). So H_2(T) = 0 iff the image of H_2(T)
in H_2(T, T\v) is trivial, which happens iff the connecting map
δ: H_2(T, T\v) → H_1(T\v) is injective, OR H_2(T, T\v) = 0.

Test: is there always a vertex v such that H_2(T, T\v) = 0?
"""
import numpy as np
from itertools import combinations, permutations
import sys
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    path_betti_numbers
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

def delete_vertex(A, n, v):
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B

def compute_relative_h2(A, n, v):
    """Compute H_2(T, T\v) = ker(∂_2^rel) / im(∂_3^rel).

    The relative chain complex is Ω_p(T) / Ω_p(T\v), where
    Ω_p(T\v) embeds into Ω_p(T) as paths not using v.

    Concretely:
    - Ω_p^{rel} = Ω_p(T) / (Ω_p(T) ∩ <paths not using v>)
    - ∂^{rel} is the induced boundary map
    """
    # Compute Ω_p(T) for p=1,2,3
    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(5):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    # For the relative complex, we need to quotient out paths not using v.
    # In practice: Ω_p^{rel} = span of Ω_p elements that have nonzero
    # component on paths using v.

    # More precisely: let V_p = subspace of Ω_p spanned by paths NOT using v.
    # Then Ω_p^{rel} = Ω_p / (Ω_p ∩ V_p).
    # dim(Ω_p^{rel}) = dim(Ω_p) - dim(Ω_p ∩ V_p).

    # Paths using v:
    uses_v = {}
    for p in range(4):
        uses_v[p] = np.array([1 if v in path else 0 for path in a[p]])

    # Ω_p ∩ V_p: subspace of Ω_p consisting of elements with zero
    # component on all v-paths.
    # In Ω_p basis: a vector c is in Ω_p ∩ V_p iff om[p] @ c has zero
    # on all v-path indices.

    rel_dims = {}
    rel_kers = {}
    rel_ims = {}

    for p in [1, 2, 3]:
        dim_om_p = om[p].shape[1] if om[p].ndim == 2 and om[p].shape[1] > 0 else 0
        if dim_om_p == 0:
            rel_dims[p] = 0
            continue

        # v-path indices
        v_idx = [i for i in range(len(a[p])) if v in a[p][i]]
        if not v_idx:
            rel_dims[p] = 0
            continue

        # Restriction matrix: project om[p] onto v-path rows
        R_v = om[p][v_idx, :]  # |v_idx| × dim_om_p

        # dim(Ω_p ∩ V_p) = dim(ker R_v)
        S = np.linalg.svd(R_v, compute_uv=False)
        rank_Rv = sum(s > 1e-8 for s in S)
        dim_inter = dim_om_p - rank_Rv
        rel_dims[p] = dim_om_p - dim_inter  # = rank_Rv

    # Now compute boundaries in the relative complex
    # ∂_2^{rel}: Ω_2^{rel} → Ω_1^{rel}
    # ∂_3^{rel}: Ω_3^{rel} → Ω_2^{rel}

    # For dim_2_rel and dim_3_rel > 0:
    dim_2_rel = rel_dims.get(2, 0)
    dim_3_rel = rel_dims.get(3, 0)

    if dim_2_rel == 0:
        return 0, dim_2_rel, dim_3_rel

    # We need to work in the quotient space.
    # Approach: find a basis for Ω_p^{rel} (complement of Ω_p ∩ V_p in Ω_p)
    # and compute boundary maps in this basis.

    # Actually, simpler approach: compute rank of ∂_p restricted to the
    # "v-essential" part.

    # For ∂_2: Ω_2 → Ω_1
    # In the quotient Ω_2/V_2 → Ω_1/V_1:
    # ker(∂_2^{rel}) = (ker ∂_2 + V_2)/V_2 = ker ∂_2 / (ker ∂_2 ∩ V_2)

    # Hmm, this is the standard diagram chase.
    # ker(∂_2^{rel}) / im(∂_3^{rel}) = H_2(T, T\v)

    # Let me just compute this numerically.
    # Step 1: Find Ω_2, its v-essential part, and the boundary maps.

    dim_om2 = om[2].shape[1] if om[2].ndim == 2 else 0
    dim_om3 = om[3].shape[1] if om[3].ndim == 2 else 0

    if dim_om2 == 0:
        return 0, 0, 0

    # ∂_2 in A coordinates
    bd2 = build_full_boundary_matrix(a[2], a[1])
    bd2_om = bd2 @ om[2]  # |A_1| × dim_om2

    # ker(∂_2) in Ω_2: kernel of bd2_om
    U2, S2, Vt2 = np.linalg.svd(bd2_om)
    r2 = sum(s > 1e-8 for s in S2)
    ker2_dim = dim_om2 - r2
    if ker2_dim == 0:
        return 0, dim_2_rel, dim_3_rel

    ker2_basis = Vt2[r2:].T  # dim_om2 × ker2_dim (in Ω_2 coords)

    # V_2 = Ω_2 ∩ <no-v paths>: kernel of R_v
    v_idx_2 = [i for i in range(len(a[2])) if v in a[2][i]]
    R_v2 = om[2][v_idx_2, :]
    U_rv, S_rv, Vt_rv = np.linalg.svd(R_v2)
    rank_rv = sum(s > 1e-8 for s in S_rv)
    V2_basis = Vt_rv[rank_rv:].T if rank_rv < dim_om2 else np.zeros((dim_om2, 0))

    # ker(∂_2) ∩ V_2: intersection
    if V2_basis.shape[1] > 0:
        # Project ker2 basis onto V2 and find common subspace
        combined = np.hstack([ker2_basis, V2_basis])
        combined_rank = np.linalg.matrix_rank(combined, tol=1e-8)
        intersect_dim = ker2_dim + V2_basis.shape[1] - combined_rank
    else:
        intersect_dim = 0

    # dim(ker ∂_2^{rel}) = ker2_dim - intersect_dim (by quotient formula)
    ker2_rel = ker2_dim - intersect_dim

    # im(∂_3^{rel}): image of ∂_3 in the quotient Ω_2/V_2
    if dim_om3 > 0:
        bd3 = build_full_boundary_matrix(a[3], a[2])
        bd3_om = bd3 @ om[3]  # |A_2| × dim_om3

        # Image of ∂_3 in Ω_2 coordinates: solve om[2] @ x = bd3_om[:, j]
        # for each j. But ∂∂=0 means im(∂_3) ⊂ ker(∂_2) ⊂ Ω_2.
        # So the image lives in Ω_2.

        # Express im(∂_3) in Ω_2 coordinates
        im3_om2, _, _, _ = np.linalg.lstsq(om[2], bd3_om, rcond=None)
        # im3_om2 columns are Ω_2-coord vectors

        # In the quotient Ω_2/V_2, the image is:
        # dim(im ∂_3^{rel}) = rank of im3_om2 modulo V_2

        # Combine im3_om2 with V_2 basis
        if V2_basis.shape[1] > 0:
            combined_im = np.hstack([im3_om2, V2_basis])
            rank_combined_im = np.linalg.matrix_rank(combined_im, tol=1e-8)
            rank_V2 = V2_basis.shape[1]
            im3_rel = rank_combined_im - rank_V2
        else:
            im3_rel = np.linalg.matrix_rank(im3_om2, tol=1e-8)
    else:
        im3_rel = 0

    h2_rel = max(0, ker2_rel - im3_rel)
    return h2_rel, dim_2_rel, dim_3_rel

# ===== Test at n=5 =====
print("=" * 70)
print("RELATIVE HOMOLOGY H_2(T, T\\v)")
print("=" * 70)

for n in [5, 6]:
    print(f"\n--- n={n} ---")
    total = 0
    has_zero_v = 0  # tournaments with some v giving H_2^rel = 0
    all_zero_v = 0  # tournaments with ALL v giving H_2^rel = 0
    max_h2_rel = 0
    h2_rel_distrib = Counter()

    for tidx, A in enumerate(all_tournaments_gen(n)):
        total += 1
        if n == 6 and total % 5000 == 0:
            print(f"  ... {total}", flush=True)

        h2_rels = []
        for v in range(n):
            h2_rel, d2r, d3r = compute_relative_h2(A, n, v)
            h2_rels.append(h2_rel)

        min_h2_rel = min(h2_rels)
        max_h2_rel_v = max(h2_rels)
        max_h2_rel = max(max_h2_rel, max_h2_rel_v)

        if min_h2_rel == 0:
            has_zero_v += 1
        if max_h2_rel_v == 0:
            all_zero_v += 1

        h2_rel_distrib[tuple(sorted(h2_rels))] += 1

    print(f"  Total tournaments: {total}")
    print(f"  ∃v: H_2(T,T\\v)=0: {has_zero_v}/{total}")
    print(f"  ∀v: H_2(T,T\\v)=0: {all_zero_v}/{total}")
    print(f"  Max H_2^rel seen: {max_h2_rel}")

    print(f"\n  H_2^rel patterns (sorted over v):")
    for pattern in sorted(h2_rel_distrib.keys()):
        count = h2_rel_distrib[pattern]
        print(f"    {pattern}: {count}")

# ===== Part 2: Connecting homomorphism =====
print(f"\n\n{'='*70}")
print("CONNECTING HOMOMORPHISM δ: H_2(T,T\\v) → H_1(T\\v)")
print("=" * 70)
print("δ must be injective for H_2(T) = 0 (when H_2(T,T\\v) ≠ 0)")

n = 5
inj_count = 0
not_inj_count = 0
h2rel_zero_count = 0

for A in all_tournaments_gen(n):
    for v in range(n):
        h2_rel, _, _ = compute_relative_h2(A, n, v)
        if h2_rel == 0:
            h2rel_zero_count += 1
            continue

        # H_1(T\v)
        B = delete_vertex(A, n, v)
        betti_sub = path_betti_numbers(B, n-1, max_dim=2)
        h1_sub = betti_sub[1]

        if h1_sub >= h2_rel:
            inj_count += 1  # δ COULD be injective (dimension allows)
        else:
            not_inj_count += 1  # δ CANNOT be injective (target too small)

print(f"  n={n}")
print(f"  H_2^rel = 0: {h2rel_zero_count}")
print(f"  H_2^rel > 0, dim H_1(T\\v) ≥ H_2^rel: {inj_count}")
print(f"  H_2^rel > 0, dim H_1(T\\v) < H_2^rel: {not_inj_count}")

print("\nDone.")
