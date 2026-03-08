#!/usr/bin/env python3
"""
beta2_local_to_global.py — Local-to-global approach for β₂=0

IDEA: For each 4-element subset S ⊆ V(T), the restriction T|_S has β₂=0.
Can we glue these local exactness results to get global β₂=0?

The approach: express every z ∈ Z₂ as a sum of "local" 2-cycles,
each supported on 4 vertices. Since local β₂=0, each piece is a boundary.

A 2-cycle z = Σ aᵢ (xᵢ, yᵢ, zᵢ) where each triple uses 3 vertices.
The boundary condition ∂₂(z) = 0 creates relations among the coefficients.

KEY OBSERVATION: Every 2-path (a,b,c) lives in the 3-vertex set {a,b,c}.
Every 3-path (a,b,c,d) lives in the 4-vertex set {a,b,c,d}.
So ∂₃ maps 4-vertex elements to 3-vertex elements.

For β₂=0, we need: for every z ∈ Z₂ (supported on various 3-vertex subsets),
there exists w ∈ Ω₃ (supported on 4-vertex subsets) with ∂₃(w) = z.

The "gluing" question: can we decompose z into local pieces, fill each locally,
and combine?

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
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

print("=" * 70)
print("LOCAL-TO-GLOBAL: DECOMPOSING Z₂ INTO 4-VERTEX PIECES")
print("=" * 70)

n = 5
print(f"\n--- n = {n} ---")

# For each tournament, try to decompose every Z₂ element into 
# boundaries from 4-vertex sub-tournaments
total_nontrivial = 0
local_fills_count = 0
local_fails_count = 0

for idx, A in enumerate(all_tournaments(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a2:
        continue
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_Om2 == 0:
        continue
    
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    dim_Z2 = dim_Om2 - rk2
    
    if dim_Z2 == 0:
        continue
    total_nontrivial += 1
    
    # Get Z₂ in A₂ coordinates
    U2, S2v, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
    Z2_om2 = Vt2[rk2_v:].T
    Z2_a2 = om2 @ Z2_om2
    
    # For each 4-vertex subset, get the 3-paths and their boundaries
    # Build combined boundary matrix from ALL 4-vertex restrictions
    if not a3:
        continue
    
    bd3 = build_full_boundary_matrix(a3, a2)
    
    # Full Ω₃ approach: check if boundaries span Z₂
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_Om3 = om3.shape[1] if om3.ndim == 2 else 0
    if dim_Om3 == 0:
        continue
    
    bd3_om = bd3 @ om3
    bd3_om2, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
    bd3_Z2 = Z2_om2.T @ bd3_om2
    rk_full = np.linalg.matrix_rank(bd3_Z2, tol=1e-8)
    
    if rk_full == dim_Z2:
        local_fills_count += 1
    else:
        local_fails_count += 1

# This is just verifying β₂=0 again. The real question is:
# can we RESTRICT to specific 4-vertex subsets?

print(f"  Total nontrivial: {total_nontrivial}")
print(f"  Full Ω₃ fills Z₂: {local_fills_count}")

# NEW APPROACH: For each 3-path (a,b,c,d), it uses vertices {a,b,c,d}.
# The key question: does the boundary ∂₃(a,b,c,d) in Ω₂ depend ONLY on
# the 4-vertex sub-tournament T|_{a,b,c,d}?
# YES — because the faces are (b,c,d), (a,c,d), (a,b,d), (a,b,c),
# all using subsets of {a,b,c,d}.

# And the Ω₃ membership: (a,b,c,d) ∈ Ω₃ iff all faces are in Ω₂.
# Face (x,y,z) ∈ Ω₂ iff it's TT or part of a cancellation combination.
# But the TT condition for (x,y,z) depends only on x→z, which is determined
# by the 4-vertex sub-tournament.

# SO: Ω₃ membership is determined by the 4-vertex restriction.
# But Z₂ elements can span MULTIPLE 3-vertex subsets.
# The filling requires combining 3-paths from different 4-vertex subsets.

# KEY EXPERIMENT: Can we fill Z₂ using ONLY 3-paths from a SINGLE 4-vertex subset?
print(f"\n--- Single 4-vertex subset filling? ---")

n = 5
single_fills = 0
multi_needed = 0

for idx, A in enumerate(all_tournaments(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a2 or not a3:
        continue
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_Om2 == 0:
        continue
    
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    dim_Z2 = dim_Om2 - rk2
    if dim_Z2 == 0:
        continue
    
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    bd3 = build_full_boundary_matrix(a3, a2)
    bd3_om = bd3 @ om3
    bd3_om2, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
    
    U2, S2v, Vt2 = np.linalg.svd(bd2 @ om2, full_matrices=True)
    rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
    Z2_om2 = Vt2[rk2_v:].T
    
    # For each 4-vertex subset, check if its Ω₃ elements fill Z₂
    filled_by_any_single = False
    for subset in combinations(range(n), 4):
        # Get Ω₃ elements supported on this subset
        subset_set = set(subset)
        om3_indices = [j for j in range(om3.shape[1]) 
                       for i in range(len(a3)) 
                       if abs(om3[i, j]) > 1e-8 and set(a3[i]) <= subset_set]
        # This is wrong — an Ω₃ basis element can span multiple a3 paths
        # Let me just check which a3 paths are in the subset
        a3_in_subset = [j for j, p in enumerate(a3) if set(p) <= subset_set]
        if not a3_in_subset:
            continue
        
        bd3_sub = bd3[:, a3_in_subset]
        bd3_sub_om2, _, _, _ = np.linalg.lstsq(om2, bd3_sub, rcond=None)
        bd3_sub_Z2 = Z2_om2.T @ bd3_sub_om2
        rk = np.linalg.matrix_rank(bd3_sub_Z2, tol=1e-8)
        
        if rk == dim_Z2:
            filled_by_any_single = True
            break
    
    if filled_by_any_single:
        single_fills += 1
    else:
        multi_needed += 1

print(f"  Z₂ filled by single 4-vertex subset: {single_fills}")
print(f"  Need multiple 4-vertex subsets: {multi_needed}")

# How many 4-vertex subsets are needed?
print(f"\n--- Minimum number of 4-vertex subsets needed ---")
n = 5
min_subsets_dist = Counter()

for idx, A in enumerate(all_tournaments(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a2 or not a3:
        continue
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_Om2 == 0:
        continue
    
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    dim_Z2 = dim_Om2 - rk2
    if dim_Z2 == 0:
        min_subsets_dist[0] += 1
        continue
    
    bd3 = build_full_boundary_matrix(a3, a2)
    
    U2, S2v, Vt2 = np.linalg.svd(bd2 @ om2, full_matrices=True)
    rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
    Z2_om2 = Vt2[rk2_v:].T
    
    # Greedy: accumulate a3 paths from different 4-vertex subsets
    all_subsets = list(combinations(range(n), 4))
    accumulated_a3 = []
    for subset in all_subsets:
        subset_set = set(subset)
        new_a3 = [j for j, p in enumerate(a3) if set(p) <= subset_set and j not in accumulated_a3]
        accumulated_a3.extend(new_a3)
        
        if accumulated_a3:
            bd3_acc = bd3[:, accumulated_a3]
            bd3_acc_om2, _, _, _ = np.linalg.lstsq(om2, bd3_acc, rcond=None)
            bd3_acc_Z2 = Z2_om2.T @ bd3_acc_om2
            rk = np.linalg.matrix_rank(bd3_acc_Z2, tol=1e-8)
            
            if rk == dim_Z2:
                # How many subsets did we use?
                n_used = all_subsets.index(subset) + 1
                min_subsets_dist[n_used] += 1
                break

print(f"  Min 4-vertex subsets needed (greedy): {dict(sorted(min_subsets_dist.items()))}")

# THE REAL INSIGHT: Let me count which 4-vertex patterns matter
print(f"\n{'='*70}")
print("THE REAL INSIGHT: BOUNDARY MAP ∂₃ STRUCTURE")
print("=" * 70)
print(f"""
The boundary ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c).

Each face uses 3 of the 4 vertices. The 2-cycle z uses various 3-vertex subsets.
For z to be in im(∂₃), we need 3-paths (a,b,c,d) whose boundaries overlap with z.

The COMPLETENESS of tournaments means:
- For any 3-vertex set {{x,y,z}}, there exist allowed 2-paths
- For any 4-vertex set {{a,b,c,d}}, there exist allowed 3-paths
- The abundance of 3-paths creates enough boundaries to fill any 2-cycle

This is essentially a COVERING argument:
Every 2-path appears as a face of MANY 3-paths (from different 4-vertex sets).
The question is whether these face-appearances create linearly independent 
contributions to Z₂.

At n=4: every 3-path uses ALL 4 vertices, so there's only one 4-vertex subset.
  β₂=0 at n=4 means the single-subset argument works.
At n=5: C(5,4)=5 subsets. We need at most a few of them.
""")

print("Done.")
