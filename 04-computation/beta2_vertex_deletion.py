#!/usr/bin/env python3
"""β_2 = 0 via vertex deletion long exact sequence.

GLMY path homology has a Mayer-Vietoris type sequence.
For a vertex v in digraph G, there's a long exact sequence:
... → H_n(G\v) → H_n(G) → H_n(G, G\v) → H_{n-1}(G\v) → ...

where G\v = G with vertex v removed, and H_n(G, G\v) is the relative
homology.

For TOURNAMENTS: removing a vertex v from T_n gives T_{n-1} (still a tournament).
If β_2(T_{n-1}) = 0 by induction, and if the relative term
H_2(T, T\v) = 0, then β_2(T) = 0.

TEST: Is there a vertex v in each tournament such that deletion simplifies β_2?
"""
import numpy as np
from itertools import combinations
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
    """Return adjacency matrix with vertex v removed."""
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B

# ===== Part 1: Verify β_2 of subtournaments =====
print("=" * 70)
print("VERTEX DELETION AND β_2")
print("=" * 70)

for n in [5, 6]:
    print(f"\n--- n={n} ---")
    # For each tournament, delete each vertex and check β_2
    count = 0
    always_zero_sub = 0
    total = 0

    for A in all_tournaments_gen(n):
        count += 1
        if n == 6 and count % 5000 == 0:
            print(f"  ... {count}", flush=True)
        total += 1

        all_sub_b2_zero = True
        for v in range(n):
            B = delete_vertex(A, n, v)
            betti = path_betti_numbers(B, n-1, max_dim=3)
            if betti[2] > 0:
                all_sub_b2_zero = False
                break

        if all_sub_b2_zero:
            always_zero_sub += 1

    print(f"  All {n-1}-vertex subtournaments have β_2=0: {always_zero_sub}/{total}")

# ===== Part 2: Relative homology H_2(T, T\v) =====
print(f"\n\n{'='*70}")
print("RELATIVE HOMOLOGY H_2(T, T\\v)")
print("="*70)
print("H_2(T, T\\v) = Ω_2(T)/Ω_2(T\\v) modulo boundaries")

# For a tournament T on [n] and vertex v:
# Ω_2(T) = paths in T that are in Ω_2
# Ω_2(T\v) = paths in T\v that are in Ω_2(T\v)
# The relative chain space = Ω_2(T) / (Ω_2(T) restricted to T\v)
#   = Ω_2 paths that involve vertex v

# More precisely, the relative chain complex is:
# C_k^{rel} = Ω_k(T) / Ω_k(T\v)
# where Ω_k(T\v) is naturally embedded in Ω_k(T) (paths not using v)

n = 5
print(f"\n--- n={n}, checking relative dimensions ---")

rel_dims = Counter()
for A in all_tournaments_gen(n):
    for v in range(n):
        # Ω_k(T) dimensions
        dims_T = {}
        for d in range(n):
            ap = enumerate_allowed_paths(A, n, d)
            if d <= 1:
                dims_T[d] = len(ap) if d == 1 else n
            else:
                ap_dm1 = enumerate_allowed_paths(A, n, d-1)
                om = compute_omega_basis(A, n, d, ap, ap_dm1)
                dims_T[d] = om.shape[1] if om.ndim == 2 else 0

        # Ω_k(T\v) dimensions
        B = delete_vertex(A, n, v)
        dims_sub = {}
        for d in range(n-1):
            ap = enumerate_allowed_paths(B, n-1, d)
            if d <= 1:
                dims_sub[d] = len(ap) if d == 1 else n-1
            else:
                ap_dm1 = enumerate_allowed_paths(B, n-1, d-1)
                om = compute_omega_basis(B, n-1, d, ap, ap_dm1)
                dims_sub[d] = om.shape[1] if om.ndim == 2 else 0

        # Relative dimensions
        rel = (dims_T[2] - dims_sub.get(2, 0),
               dims_T[3] - dims_sub.get(3, 0))

        # But wait: we need to count Ω_2(T) elements that DON'T use v
        # vs total Ω_2(T) elements
        # The embedding Ω_k(T\v) → Ω_k(T) is not straightforward because
        # the basis might mix paths with and without v.

        # Instead, just count: how many A_k paths use v?
        a2_T = enumerate_allowed_paths(A, n, 2)
        a2_with_v = sum(1 for p in a2_T if v in p)
        a2_without_v = len(a2_T) - a2_with_v

        a3_T = enumerate_allowed_paths(A, n, 3)
        a3_with_v = sum(1 for p in a3_T if v in p)
        a3_without_v = len(a3_T) - a3_with_v

        rel_dims[(a2_with_v, a2_without_v, a3_with_v, a3_without_v)] += 1

    break  # Just first tournament for now

print(f"  (A2_with_v, A2_without, A3_with_v, A3_without): count")
for k in sorted(rel_dims):
    print(f"    {k}: {rel_dims[k]}")

# ===== Part 3: Link of a vertex =====
print(f"\n\n{'='*70}")
print("LINK OF VERTEX v: PATHS STARTING OR ENDING AT v")
print("="*70)

# For the β_2=0 proof, the key observation might be about the "link" of v:
# Paths in Ω_2 that use v contribute to H_2 if they can't be filled.
# The tournament structure ensures every such path CAN be filled because
# v has edges to ALL other vertices.

n = 5
link_data = Counter()
for A in all_tournaments_gen(n):
    for v in range(n):
        out_v = sum(A[v])
        in_v = sum(A[j][v] for j in range(n))

        # 2-paths involving v: (v,*,*), (*,v,*), (*,*,v)
        a2_T = enumerate_allowed_paths(A, n, 2)
        starts_v = sum(1 for p in a2_T if p[0] == v)
        mid_v = sum(1 for p in a2_T if p[1] == v)
        ends_v = sum(1 for p in a2_T if p[2] == v)

        link_data[(out_v, starts_v, mid_v, ends_v)] += 1

print(f"  (out_deg, starts_v, mid_v, ends_v): count")
for k in sorted(link_data):
    print(f"    {k}: {link_data[k]}")

# ===== Part 4: Can we find an acyclic matching on Ω_2? =====
# If there's an acyclic matching (in the discrete Morse theory sense)
# that pairs every Ω_2 element with an Ω_3 element, then H_2 = 0.
print(f"\n\n{'='*70}")
print("DISCRETE MORSE THEORY: ACYCLIC MATCHING Ω_2 ↔ Ω_3")
print("="*70)

# For each tournament at n=5, try to find a matching between
# ker(∂_2) elements and Ω_3 elements via ∂_3.
# If each ker(∂_2) vector is the boundary of a unique Ω_3 element,
# we have an acyclic matching.

n = 5
perfect_match = 0
imperfect = 0
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0: continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        perfect_match += 1  # trivially
        continue

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd3 = build_full_boundary_matrix(a3_list, a2_list)
    if dim_om3 > 0:
        im3 = bd3 @ om3
        im3_in_om2, _, _, _ = np.linalg.lstsq(om2, im3, rcond=None)
        rank3 = np.linalg.matrix_rank(im3_in_om2, tol=1e-8)
    else:
        rank3 = 0

    if rank3 == ker_dim:
        # Check: is rank(∂_3) = dim(Ω_3)? (i.e., ∂_3 injective?)
        if rank3 == dim_om3:
            perfect_match += 1  # ∂_3 is injective, so it's a "perfect" filling
        else:
            imperfect += 1  # ∂_3 has kernel (β_3 related), but still fills all of ker(∂_2)
    else:
        print(f"  UNEXPECTED: rank3={rank3} < ker_dim={ker_dim}")

print(f"  ∂_3 injective (dim(ker ∂_3)=0): {perfect_match}")
print(f"  ∂_3 has kernel but still fills: {imperfect}")

print("\nDone.")
