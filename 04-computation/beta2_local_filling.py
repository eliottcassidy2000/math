#!/usr/bin/env python3
"""β_2 = 0: LOCAL FILLING APPROACH

Key idea: every 2-cycle z in a tournament can be filled LOCALLY,
meaning that for every 4-vertex subset {a,b,c,d} containing a
nonzero component of z, the restriction of z to that subset is
either already a boundary or can be extended.

The 4-vertex tournament determines the local filling.
There are exactly 4 tournament types on 4 vertices (up to isomorphism):
1. Transitive: 0→1→2→3, 0→2, 0→3, 1→3
2. One 3-cycle: has exactly one 3-cycle
3. Two 3-cycles: "diamond" pattern
4. Not possible at n=4 (tournaments on 4 vertices have 0, 1, or 2 three-cycles)

For each type, compute H_2 and the local filling capacity.

Also: test if every 2-cycle in a tournament at n=5 can be decomposed
as a sum of 2-boundaries coming from 4-vertex subtournaments.
"""
import numpy as np
from itertools import combinations, permutations
import sys
from collections import Counter, defaultdict
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

# ===== Part 1: All n=4 tournament types =====
print("=" * 70)
print("ALL n=4 TOURNAMENT TYPES: CHAIN COMPLEX")
print("=" * 70)

type_data = defaultdict(list)
for A in all_tournaments_gen(4):
    t3 = sum(1 for a, b, c in combinations(range(4), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or
                (A[b][a] and A[a][c] and A[c][b]))
    betti = path_betti_numbers(A, 4, max_dim=3)

    a2 = enumerate_allowed_paths(A, 4, 2)
    a3 = enumerate_allowed_paths(A, 4, 3)
    a1 = enumerate_allowed_paths(A, 4, 1)

    om2 = compute_omega_basis(A, 4, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    om3 = compute_omega_basis(A, 4, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    type_data[t3].append({
        'betti': tuple(betti),
        'dim_om2': dim_om2,
        'dim_om3': dim_om3,
        'n_a2': len(a2),
        'n_a3': len(a3),
    })

for t3 in sorted(type_data):
    examples = type_data[t3]
    r = examples[0]
    print(f"\nt3={t3} ({len(examples)} tournaments):")
    print(f"  |A_2|={r['n_a2']}, dim(Ω_2)={r['dim_om2']}, |A_3|={r['n_a3']}, dim(Ω_3)={r['dim_om3']}")
    print(f"  β = {r['betti']}")

# ===== Part 2: Subtournament boundary decomposition at n=5 =====
print(f"\n\n{'='*70}")
print("SUBTOURNAMENT BOUNDARY DECOMPOSITION AT n=5")
print("="*70)
print("Can every 2-cycle be decomposed into boundaries from 4-vertex subtournaments?")

n = 5
success = 0
failure = 0
total = 0

for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_list = [tuple(p) for p in a1]
    a2_list = [tuple(p) for p in a2]
    a3_list = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0: continue

    bd2 = build_full_boundary_matrix(a2_list, a1_list)
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0: continue

    total += 1

    # ker basis
    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    # For each 4-vertex subset, compute im(∂_3) restricted to paths within that subset
    a2_idx = {p: i for i, p in enumerate(a2_list)}

    # Collect all DT 3-paths from ALL 4-vertex subtournaments
    all_sub_bds = []
    for subset in combinations(range(n), 4):
        s = set(subset)
        # Find DT paths within this subset
        for p in a3_list:
            if set(p) <= s:
                # This is a DT path in the subtournament
                # Its boundary in A_2
                a, b, c, d = p
                bd_vec = np.zeros(len(a2_list))
                faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
                for fi, f in enumerate(faces):
                    if f in a2_idx:
                        bd_vec[a2_idx[f]] += (-1)**fi
                all_sub_bds.append(bd_vec)

    if not all_sub_bds:
        failure += ker_dim
        continue

    sub_bd_matrix = np.column_stack(all_sub_bds)

    # Check: do these boundaries span ker(∂_2|Ω_2)?
    # Project to Ω_2 coordinates
    sub_in_om2, _, _, _ = np.linalg.lstsq(om2, sub_bd_matrix, rcond=None)

    for ci in range(ker_dim):
        z_om2 = Vt[rank2 + ci, :]
        # Check if z is in span of sub_in_om2
        combined = np.column_stack([sub_in_om2, z_om2.reshape(-1,1)])
        rank_without = np.linalg.matrix_rank(sub_in_om2, tol=1e-8)
        rank_with = np.linalg.matrix_rank(combined, tol=1e-8)

        if rank_with == rank_without:
            success += 1
        else:
            failure += 1

print(f"\nTotal tournaments with ker > 0: {total}")
print(f"2-cycles filled by subtournament DT: {success}")
print(f"2-cycles NOT filled: {failure}")

# ===== Part 3: What's special about the boundary of ∂_3 in a tournament? =====
# The key difference from oriented graphs:
# In a tournament, for EVERY pair (a,b), either a→b or b→a.
# This means: for every 2-path (x,y,z) and every vertex w,
# we can form (w,x,y,z) or (x,y,z,w) and check if it's DT.
# The probability of being DT is controlled by the tournament structure.

print(f"\n\n{'='*70}")
print("DT PATH DENSITY AT n=5")
print("="*70)

dt_density = Counter()
for A in all_tournaments_gen(5):
    a3 = enumerate_allowed_paths(A, 5, 3)
    a2 = enumerate_allowed_paths(A, 5, 2)
    a2_set = set(tuple(p) for p in a2)

    dt_count = 0
    for p in a3:
        a, b, c, d = tuple(p)
        if A[a][c] and A[b][d]:
            dt_count += 1

    t3 = sum(1 for a, b, c in combinations(range(5), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or
                (A[b][a] and A[a][c] and A[c][b]))
    dt_density[(t3, dt_count)] += 1

print(f"(t3, #DT): count")
for k in sorted(dt_density):
    print(f"  {k}: {dt_density[k]}")

print("\nDone.")
