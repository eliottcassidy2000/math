#!/usr/bin/env python3
"""
beta2_ses_proof.py — Short exact sequence proof attempt for β₂=0

PROOF STRATEGY via vertex deletion:

Given tournament T on n vertices, fix vertex v.
Let T' = T \ v (subtournament on n-1 vertices).

Chain complexes:
  Ω_*(T') ↪ Ω_*(T) → Ω_*(T)/Ω_*(T') =: R_*(v)

This gives a short exact sequence of chain complexes:
  0 → Ω_*(T') → Ω_*(T) → R_*(v) → 0

Long exact sequence in homology:
  ... → H_2(T') → H_2(T) → H_2(R_*(v)) → H_1(T') → H_1(T) → ...

If β_2(T') = 0 (induction hypothesis) then:
  0 → H_2(T) → H_2(R_*(v)) → H_1(T') → H_1(T)

So β_2(T) ≤ dim H_2(R_*(v)).

If H_2(R_*(v)) = 0 for ALL v, then β_2(T) = 0.

WHAT IS R_p(v)?
R_p(v) = Ω_p(T) / Ω_p(T')
= elements of Ω_p(T) that genuinely involve vertex v.

More precisely: R_p(v) is spanned by Ω_p elements that use v in at least one path.

THIS SCRIPT: Compute H_2(R_*(v)) for all tournaments at n=4,5 and verify it's 0.

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def compute_relative_H2(A, n, v):
    """
    Compute H_2 of the relative complex R_*(v) = Ω_*(T)/Ω_*(T\\v).
    
    Method: compute Ω_p(T) and Ω_p(T\\v), then form the quotient complex
    and compute its homology.
    """
    # Full tournament
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    om1 = np.eye(len(a1)) if a1 else np.zeros((0,0))
    om2 = compute_omega_basis(A, n, 2, a2, a1) if a2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, a3, a2) if a3 else np.zeros((0,0))
    
    dim_Om2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0
    dim_Om3 = om3.shape[1] if om3.ndim == 2 and om3.shape[0] > 0 else 0
    
    if dim_Om2 == 0:
        return 0, {'trivial': True}
    
    # Subtournament T' = T\v
    verts_sub = [i for i in range(n) if i != v]
    n_sub = n - 1
    
    # Map: vertex i in T' → position in verts_sub
    # A_sub[remapped_i][remapped_j] = A[verts_sub[i]][verts_sub[j]]
    A_sub = [[A[verts_sub[i]][verts_sub[j]] for j in range(n_sub)] for i in range(n_sub)]
    
    a1_sub = enumerate_allowed_paths(A_sub, n_sub, 1)
    a2_sub = enumerate_allowed_paths(A_sub, n_sub, 2)
    a3_sub = enumerate_allowed_paths(A_sub, n_sub, 3)
    
    om2_sub = compute_omega_basis(A_sub, n_sub, 2, a2_sub, a1_sub) if a2_sub else np.zeros((0,0))
    om3_sub = compute_omega_basis(A_sub, n_sub, 3, a3_sub, a2_sub) if a3_sub else np.zeros((0,0))
    
    dim_Om2_sub = om2_sub.shape[1] if om2_sub.ndim == 2 and om2_sub.shape[0] > 0 else 0
    dim_Om3_sub = om3_sub.shape[1] if om3_sub.ndim == 2 and om3_sub.shape[0] > 0 else 0
    
    # Embed Ω_p(T') into Ω_p(T):
    # A 2-path (a',b',c') in T' maps to (verts_sub[a'], verts_sub[b'], verts_sub[c']) in T.
    # This is an allowed 2-path in T that doesn't use v.
    
    # For the relative complex, we need:
    # R_2 = Ω_2(T) / Ω_2(T'), dim = dim_Om2 - dim_Om2_sub
    # R_3 = Ω_3(T) / Ω_3(T'), dim = dim_Om3 - dim_Om3_sub
    
    dim_R2 = dim_Om2 - dim_Om2_sub
    dim_R3 = dim_Om3 - dim_Om3_sub
    
    if dim_R2 == 0:
        return 0, {'dim_R2': 0, 'dim_R3': dim_R3}
    
    # To compute H_2(R_*), we need the boundary maps on the quotient.
    # This requires constructing the actual quotient bases.
    
    # Step 1: Embed Ω_2(T') into Ω_2(T) via the inclusion map
    # Build the embedding matrix
    
    if dim_Om2_sub > 0:
        # Map T' paths to T paths
        # For each T' 2-path (a',b',c'), find it in T's a2 list
        sub_to_full_2 = {}
        for i, p in enumerate(a2_sub):
            full_path = tuple(verts_sub[x] for x in p)
            if full_path in [tuple(q) for q in a2]:
                j = [tuple(q) for q in a2].index(full_path)
                sub_to_full_2[i] = j
        
        # Embedding: map om2_sub columns to vectors in A₂(T) space
        embed_2 = np.zeros((len(a2), dim_Om2_sub))
        for col in range(dim_Om2_sub):
            for i in range(len(a2_sub)):
                if abs(om2_sub[i, col]) > 1e-10 and i in sub_to_full_2:
                    embed_2[sub_to_full_2[i], col] = om2_sub[i, col]
        
        # Project embedding into Ω_2(T) coordinates
        embed_om2, _, _, _ = np.linalg.lstsq(om2, embed_2, rcond=None)
        # embed_om2: dim_Om2 × dim_Om2_sub
        
        # R_2 basis: complement of im(embed_om2) in Ω_2(T)
        U, S, Vt = np.linalg.svd(embed_om2, full_matrices=True)
        rk_embed = int(np.sum(np.abs(S) > 1e-8))
        R2_basis = U[:, rk_embed:]  # dim_Om2 × dim_R2, in Ω_2(T) coords
    else:
        R2_basis = np.eye(dim_Om2)
    
    actual_dim_R2 = R2_basis.shape[1]
    
    # Step 2: Same for dimension 3
    if dim_Om3_sub > 0 and dim_Om3 > 0:
        sub_to_full_3 = {}
        a3_tuples = [tuple(q) for q in a3]
        for i, p in enumerate(a3_sub):
            full_path = tuple(verts_sub[x] for x in p)
            if full_path in a3_tuples:
                j = a3_tuples.index(full_path)
                sub_to_full_3[i] = j
        
        embed_3 = np.zeros((len(a3), dim_Om3_sub))
        for col in range(dim_Om3_sub):
            for i in range(len(a3_sub)):
                if abs(om3_sub[i, col]) > 1e-10 and i in sub_to_full_3:
                    embed_3[sub_to_full_3[i], col] = om3_sub[i, col]
        
        embed_om3, _, _, _ = np.linalg.lstsq(om3, embed_3, rcond=None)
        U3, S3, Vt3 = np.linalg.svd(embed_om3, full_matrices=True)
        rk_embed3 = int(np.sum(np.abs(S3) > 1e-8))
        R3_basis = U3[:, rk_embed3:]
    elif dim_Om3 > 0:
        R3_basis = np.eye(dim_Om3)
    else:
        R3_basis = np.zeros((0, 0))
    
    actual_dim_R3 = R3_basis.shape[1] if R3_basis.ndim == 2 else 0
    
    # Step 3: Boundary map ∂₂^R: R_2 → R_1
    # In Ω_2(T) coords, ∂₂ is given by bd2 @ om2.
    # Restricted to R_2 basis: ∂₂^R = bd2_om2_restricted_to_R2
    
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    # ∂₂ on R_2: bd2_om @ R2_basis
    d2_R = bd2_om @ R2_basis  # in A_1 space
    
    S_d2R = np.linalg.svd(d2_R, compute_uv=False)
    rk_d2_R = int(np.sum(np.abs(S_d2R) > 1e-8))
    ker_d2_R = actual_dim_R2 - rk_d2_R
    
    # Step 4: Boundary map ∂₃^R: R_3 → R_2
    if actual_dim_R3 > 0 and a3:
        bd3 = build_full_boundary_matrix(a3, a2)
        bd3_om = bd3 @ om3
        # ∂₃ in Ω₂ coords: project bd3_om into Ω₂
        bd3_om2, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
        # ∂₃ restricted to R_3, projected onto R_2:
        d3_R2 = R2_basis.T @ bd3_om2 @ R3_basis  # dim_R2 × dim_R3
        
        S_d3R = np.linalg.svd(d3_R2, compute_uv=False)
        rk_d3_R = int(np.sum(np.abs(S_d3R) > 1e-8))
    else:
        rk_d3_R = 0
    
    # H_2(R_*) = ker(∂₂^R) / im(∂₃^R)
    H2_R = ker_d2_R - rk_d3_R
    
    return H2_R, {
        'dim_R2': actual_dim_R2, 'dim_R3': actual_dim_R3,
        'ker_d2_R': ker_d2_R, 'rk_d3_R': rk_d3_R,
        'rk_d2_R': rk_d2_R,
    }

print("=" * 70)
print("RELATIVE HOMOLOGY H₂(R_*(v)) = H₂(T, T\\v)")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n = {n} ---")
    
    h2_nonzero = 0
    total = 0
    h2_dist = Counter()
    detail_dist = Counter()
    
    for A in all_tournaments(n):
        for v in range(n):
            h2, info = compute_relative_H2(A, n, v)
            total += 1
            h2_dist[h2] += 1
            if h2 > 0:
                h2_nonzero += 1
                if h2_nonzero <= 3:
                    print(f"  NONZERO H₂(R_*({v})): H₂={h2}, info={info}")
            
            if not info.get('trivial', False):
                key = (info['dim_R2'], info['dim_R3'], info['ker_d2_R'], info['rk_d3_R'])
                detail_dist[key] += 1
    
    print(f"  Total (tournament, vertex) pairs: {total}")
    print(f"  H₂(R_*(v)) distribution: {dict(h2_dist)}")
    print(f"  Nonzero: {h2_nonzero}")
    
    if h2_nonzero == 0:
        print(f"  *** H₂(T, T\\v) = 0 for ALL tournaments at n={n}! ***")
    
    print(f"\n  Detail distribution (dim_R2, dim_R3, ker∂₂^R, rk∂₃^R):")
    for key, count in sorted(detail_dist.items()):
        h2 = key[2] - key[3]
        print(f"    {key}: count={count}, H₂={h2}")

# Test at n=6 (sampled)
print(f"\n--- n = 6 (sampled) ---")
import random
random.seed(42)
n = 6
h2_nonzero = 0
total = 0

pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

for _ in range(500):
    mask = random.randint(0, (1 << m) - 1)
    A = [[0]*n for _ in range(n)]
    for idx2, (i,j) in enumerate(pairs):
        if (mask >> idx2) & 1: A[i][j] = 1
        else: A[j][i] = 1
    
    for v in range(n):
        h2, info = compute_relative_H2(A, n, v)
        total += 1
        if h2 > 0:
            h2_nonzero += 1
            print(f"  NONZERO at n=6: v={v}, H₂={h2}, info={info}")

print(f"  Tested: {total}")
print(f"  H₂ nonzero: {h2_nonzero}")
if h2_nonzero == 0:
    print(f"  *** H₂(T, T\\v) = 0 for all sampled n=6 tournaments! ***")

print(f"""
PROOF STRUCTURE:
1. Base case: β₂(T) = 0 for n ≤ 3 (trivial) and n = 4 (exhaustive check: 4 classes)
2. Induction step: Given β₂(T') = 0 for all (n-1)-tournaments:
   - Short exact sequence: 0 → Ω_*(T\\v) → Ω_*(T) → R_*(v) → 0
   - Long exact sequence: 0 → H₂(T) → H₂(R_*(v)) → H₁(T\\v) → ...
   - If H₂(R_*(v)) = 0 for some v, then β₂(T) = 0. ✓
   - REMAINING: Prove H₂(R_*(v)) = 0 algebraically.
   
The relative complex R_*(v) has a concrete description:
   R_2(v): Ω₂ elements using vertex v (quotient by those not using v)
   R_3(v): Ω₃ elements using vertex v
   
H₂(R_*(v)) = 0 means: every relative 2-cycle is a relative 3-boundary.
""")

print("Done.")
