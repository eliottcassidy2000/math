#!/usr/bin/env python3
"""Try to prove β_2=0 using Hamiltonian path homotopy.

By Rédei's theorem, every tournament has a Hamiltonian path.
Can we use this to construct a chain homotopy s: Ω_2 → Ω_3
such that ∂_3 ∘ s = id on ker(∂_2)?

Alternative: use the Hamiltonian path to define a FILTRATION
of the chain complex that makes exactness at Ω_2 transparent.

Key idea: the Hamiltonian path v_0 → v_1 → ... → v_{n-1} gives
a total order compatible with many edges. For vertices "close"
in this order, we can extend 2-paths to 3-paths.
"""
import numpy as np
from itertools import permutations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
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

def find_hamiltonian_path(A, n):
    """Find a Hamiltonian path in tournament (always exists by Rédei)."""
    # Simple greedy: insert each vertex into the path
    path = [0]
    for v in range(1, n):
        # Find where to insert v
        inserted = False
        for i in range(len(path)):
            if A[v][path[i]]:
                path.insert(i, v)
                inserted = True
                break
        if not inserted:
            path.append(v)
    return path

def try_filling_map(A, n, hp, z_coeffs, a2_list, a3_list):
    """Try to fill 2-cycle z using Hamiltonian path hp.

    Strategy: for each 2-path (a,b,c) in z with coefficient z_{abc},
    extend to a 3-path using the HP ordering.

    Several strategies to try:
    1. Left-extend: find v before a in HP, use (v,a,b,c)
    2. Right-extend: find v after c in HP, use (a,b,c,v)
    3. Insert: find v between two, use (a,v,b,c) or (a,b,v,c)
    """
    a3_set = set(a3_list)
    hp_pos = {v: i for i, v in enumerate(hp)}

    filled = np.zeros(len(a3_list))
    residual_2 = z_coeffs.copy()

    # For each nonzero term in z, try to extend
    for idx2 in range(len(a2_list)):
        coeff = z_coeffs[idx2]
        if abs(coeff) < 1e-12:
            continue
        a_, b_, c_ = a2_list[idx2]

        # Try right extension: (a,b,c,v) for each v ∉ {a,b,c} with c→v
        best_v = None
        for v in range(n):
            if v in {a_, b_, c_}:
                continue
            if not A[c_][v]:
                continue
            path = (a_, b_, c_, v)
            if path in a3_set:
                # Check if this is DT (all faces in A_2)
                is_dt = A[a_][c_] and A[b_][v]
                if is_dt:
                    best_v = v
                    break
                elif best_v is None:
                    best_v = v

    return None  # placeholder

# ===== Part 1: Analyze Hamiltonian path structure =====
print("=" * 70)
print("HAMILTONIAN PATH AND CYCLE STRUCTURE")
print("=" * 70)

n = 5
for tidx, A in enumerate(all_tournaments_gen(n)):
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    if scores != (2,2,2,2,2):
        continue

    hp = find_hamiltonian_path(A, n)
    print(f"\nTournament #{tidx}: scores={scores}")
    print(f"  Hamiltonian path: {hp}")
    print(f"  HP edges: {[(hp[i], hp[i+1]) for i in range(n-1)]}")

    # What fraction of edges align with HP?
    hp_pos = {v: i for i, v in enumerate(hp)}
    aligned = sum(1 for i in range(n) for j in range(n) if i != j and A[i][j] and hp_pos[i] < hp_pos[j])
    total_edges = n * (n-1) // 2
    print(f"  Aligned edges: {aligned}/{total_edges}")

    # Which 2-paths are "HP-monotone" (all three vertices in HP order)?
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    monotone = sum(1 for (a,b,c) in a2 if hp_pos[a] < hp_pos[b] < hp_pos[c])
    print(f"  HP-monotone 2-paths: {monotone}/{len(a2)}")

    break

# ===== Part 2: Test a specific filling construction =====
print(f"\n\n{'='*70}")
print("FILLING CONSTRUCTION: GREEDY EXTENSION")
print("=" * 70)

# Strategy: for each 2-cycle z, greedily fill it by extending each term
# to a DT 4-path, adjusting for boundary cancellation.

n = 5
total_attempts = 0
total_success = 0

for tidx, A in enumerate(all_tournaments_gen(n)):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)

    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    _, S2, Vt2 = np.linalg.svd(bd2_om)
    rank_d2 = sum(s > 1e-8 for s in S2)
    z2_dim = dim_om2 - rank_d2

    if z2_dim == 0:
        continue

    # Get ∂_3 map in Ω coords
    bd3 = build_full_boundary_matrix(a3, a2)
    bd3_om = bd3 @ om3
    im3_coords, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)

    # Get Z_2 basis
    ker2_basis = Vt2[rank_d2:].T

    # For each basis element of Z_2, check if it can be filled
    for j in range(z2_dim):
        cycle_om = ker2_basis[:, j]
        # Try to solve im3_coords @ w = cycle_om
        w, residual, _, _ = np.linalg.lstsq(im3_coords, cycle_om, rcond=None)
        error = np.linalg.norm(im3_coords @ w - cycle_om)
        total_attempts += 1
        if error < 1e-8:
            total_success += 1

print(f"  n={n}: {total_success}/{total_attempts} 2-cycles filled (should be all)")

# ===== Part 3: Find the SIMPLEST possible filling rule =====
print(f"\n\n{'='*70}")
print("MINIMAL FILLING: DT PATHS ONLY?")
print("=" * 70)

# Question: can every 2-cycle be filled using ONLY DT 4-paths?
# We know this fails at n=6. But can we at least describe the failure mode?

n = 5
dt_only_success = 0
dt_only_total = 0

for tidx, A in enumerate(all_tournaments_gen(n)):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    _, S2, Vt2 = np.linalg.svd(bd2_om)
    rank_d2 = sum(s > 1e-8 for s in S2)
    z2_dim = dim_om2 - rank_d2
    if z2_dim == 0:
        continue

    ker2_basis = Vt2[rank_d2:].T

    # DT paths: all 3-paths where all faces are in A_2
    dt_indices = [i for i, p in enumerate(a3)
                  if A[p[0]][p[2]] and A[p[1]][p[3]]]

    # DT boundary map restricted to ker(∂_2)
    bd3 = build_full_boundary_matrix(a3, a2)
    # Boundary of DT paths in A_2 coords
    dt_bd_a2 = bd3[:, dt_indices]
    # Project into Ω_2 coords
    dt_bd_om2, _, _, _ = np.linalg.lstsq(om2, dt_bd_a2, rcond=None)

    # Check: does im(dt_bd_om2) contain all of Z_2?
    combined = np.hstack([dt_bd_om2, ker2_basis])
    rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
    rank_dt_bd = np.linalg.matrix_rank(dt_bd_om2, tol=1e-8)

    if rank_combined == rank_dt_bd:
        dt_only_success += 1
    dt_only_total += 1

print(f"  n={n}: DT fills all Z_2: {dt_only_success}/{dt_only_total}")

# ===== Part 4: What is the COMBINATORIAL structure of ∂_3? =====
print(f"\n\n{'='*70}")
print("BOUNDARY MAP STRUCTURE AT n=5")
print("=" * 70)

# For each DT 4-path (a,b,c,d), ∂_3(a,b,c,d) = (b,c,d)-(a,c,d)+(a,b,d)-(a,b,c)
# All four faces are TT paths (since a→c and b→d).
# So the DT boundary ONLY involves TT 2-paths.

# For cancellation elements in Ω_3:
# e.g., (a,b,c,d)-(a,b',c,d) where both have the same bad face
# Their boundary involves both TT and non-TT 2-paths.

# Key question: does ∂_3(DT) span all of ker(∂_2|_TT)?
# And does ∂_3(cancellation) span the remaining ker(∂_2)?

# For the regular tournament C_5^{1,2}:
A = [[0]*n for _ in range(n)]
for i in range(n):
    A[i][(i+1)%n] = 1
    A[i][(i+2)%n] = 1

a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

print(f"\nC_5^{{1,2}}: |A_2|={len(a2)}, |A_3|={len(a3)}")

# Classify A_2 paths
tt_idx = [i for i, (a,b,c) in enumerate(a2) if A[a][c]]
non_tt_idx = [i for i, (a,b,c) in enumerate(a2) if not A[a][c]]
print(f"  TT paths: {len(tt_idx)}, non-TT: {len(non_tt_idx)}")

# Classify A_3 paths
dt_idx = [i for i, (a,b,c,d) in enumerate(a3) if A[a][c] and A[b][d]]
non_dt_idx = [i for i in range(len(a3)) if i not in dt_idx]
print(f"  DT paths: {len(dt_idx)}, non-DT: {len(non_dt_idx)}")

# Boundary of DT paths: which faces do they produce?
bd3 = build_full_boundary_matrix(a3, a2)
print(f"\n  DT path boundaries (are they ALL TT 2-paths?):")
dt_all_tt = True
for i in dt_idx:
    path = a3[i]
    faces_idx = [j for j in range(len(a2)) if abs(bd3[j, i]) > 1e-10]
    for fi in faces_idx:
        if fi not in tt_idx:
            dt_all_tt = False
            a,b,c = a2[fi]
            print(f"    DT {path} has non-TT face {a2[fi]} (c→a: {A[c][a]})")
            break

print(f"  All DT boundary faces are TT? {dt_all_tt}")

# ∂_3 restricted to DT in Ω_2 coordinates
om2 = compute_omega_basis(A, n, 2, a2, a1)
dim_om2 = om2.shape[1]

# Compute ∂_3(DT) in Ω_2 coordinates
dt_bd = bd3[:, dt_idx]  # |A_2| × |DT|
dt_bd_om, _, _, _ = np.linalg.lstsq(om2, dt_bd, rcond=None)
rank_dt = np.linalg.matrix_rank(dt_bd_om, tol=1e-8)
print(f"\n  rank(∂_3|_DT) in Ω_2: {rank_dt}")

# Z_2
bd2 = build_full_boundary_matrix(a2, a1)
bd2_om = bd2 @ om2
_, S2, Vt2 = np.linalg.svd(bd2_om)
rank_d2 = sum(s > 1e-8 for s in S2)
z2 = dim_om2 - rank_d2
print(f"  Z_2 = {z2}")
print(f"  DT fills {'all' if rank_dt >= z2 else 'only ' + str(rank_dt) + '/' + str(z2)} of Z_2")

# What about ∂_3 of the FULL Ω_3?
om3 = compute_omega_basis(A, n, 3, a3, a2)
dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
bd3_om = bd3 @ om3
im3_om, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
rank_full = np.linalg.matrix_rank(im3_om, tol=1e-8)
print(f"  rank(∂_3|_Ω3) in Ω_2: {rank_full}")

# ===== Part 5: KEY TEST — is the image of ∂_3(DT) always in the TT subspace of Ω_2? =====
print(f"\n\n{'='*70}")
print("DT → TT BOUNDARY STRUCTURE")
print("=" * 70)

# If ALL DT boundary faces are TT, then im(∂_3|_DT) ⊂ TT ⊂ Ω_2.
# And if im(∂_3|_DT) = ker(∂_2|_TT), then the "TT part" is exact.
# The "cancellation part" would need separate treatment.

n = 5
all_dt_tt_check = True
for A in all_tournaments_gen(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    tt_set = set(i for i, (a,b,c) in enumerate(a2) if A[a][c])

    bd3 = build_full_boundary_matrix(a3, a2)
    dt_indices = [i for i, (a,b,c,d) in enumerate(a3) if A[a][c] and A[b][d]]

    for i in dt_indices:
        faces = [j for j in range(len(a2)) if abs(bd3[j, i]) > 1e-10]
        if not all(f in tt_set for f in faces):
            all_dt_tt_check = False
            break
    if not all_dt_tt_check:
        break

print(f"  n={n}: ∂_3(DT) ⊂ TT? {all_dt_tt_check}")

# Also check at n=4
n = 4
all_dt_tt_4 = True
for A in all_tournaments_gen(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]

    tt_set = set(i for i, (a,b,c) in enumerate(a2) if A[a][c])
    bd3 = build_full_boundary_matrix(a3, a2)
    dt_indices = [i for i, (a,b,c,d) in enumerate(a3) if A[a][c] and A[b][d]]

    for i in dt_indices:
        faces = [j for j in range(len(a2)) if abs(bd3[j, i]) > 1e-10]
        if not all(f in tt_set for f in faces):
            all_dt_tt_4 = False
            break
    if not all_dt_tt_4:
        break

print(f"  n=4: ∂_3(DT) ⊂ TT? {all_dt_tt_4}")

# ===== Part 6: Decomposition — TT subcomplex =====
print(f"\n\n{'='*70}")
print("TT SUBCOMPLEX ANALYSIS")
print("=" * 70)

# If ∂_3(DT) ⊂ TT, then we have a SUBCOMPLEX:
# ... → DT → TT → edges → vertices
# where TT ⊂ Ω_2 and DT ⊂ Ω_3 (or A_3).
# The homology of this subcomplex at degree 2 is ker(∂_2|_TT)/im(∂_3|_DT).
# If this is already 0, then the TT part is exact.

n = 5
tt_exact_count = 0
tt_total = 0

for A in all_tournaments_gen(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    tt_paths = [(a,b,c) for (a,b,c) in a2 if A[a][c]]
    dt_paths = [(a,b,c,d) for (a,b,c,d) in a3 if A[a][c] and A[b][d]]

    if not tt_paths:
        continue
    tt_total += 1

    # ∂_2 on TT paths: (b,c)-(a,c)+(a,b)
    # All faces are edges in A_1
    tt_matrix = np.zeros((len(a1), len(tt_paths)))
    a1_idx = {p: i for i, p in enumerate(a1)}
    for j, (a,b,c) in enumerate(tt_paths):
        bc = (b,c)
        ac = (a,c)
        ab = (a,b)
        if bc in a1_idx: tt_matrix[a1_idx[bc], j] += 1
        if ac in a1_idx: tt_matrix[a1_idx[ac], j] -= 1
        if ab in a1_idx: tt_matrix[a1_idx[ab], j] += 1

    # ker(∂_2|_TT)
    S_tt = np.linalg.svd(tt_matrix, compute_uv=False)
    rank_tt = sum(s > 1e-8 for s in S_tt)
    ker_tt = len(tt_paths) - rank_tt

    # im(∂_3|_DT) in TT coordinates
    if dt_paths:
        a2_idx = {p: i for i, p in enumerate(a2)}
        tt_idx_map = {p: i for i, p in enumerate(tt_paths)}

        dt_bd_tt = np.zeros((len(tt_paths), len(dt_paths)))
        for j, (a,b,c,d) in enumerate(dt_paths):
            # ∂(a,b,c,d) = (b,c,d)-(a,c,d)+(a,b,d)-(a,b,c)
            faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in tt_idx_map:
                    dt_bd_tt[tt_idx_map[face], j] += sign

        rank_dt_bd = np.linalg.matrix_rank(dt_bd_tt, tol=1e-8)
    else:
        rank_dt_bd = 0

    h2_tt = ker_tt - rank_dt_bd
    if h2_tt == 0:
        tt_exact_count += 1

print(f"  n={n}: TT subcomplex exact at degree 2: {tt_exact_count}/{tt_total}")

# ===== Part 7: Flag complex comparison =====
print(f"\n\n{'='*70}")
print("FLAG COMPLEX vs PATH COMPLEX")
print("=" * 70)

# The TT paths form the FLAG COMPLEX (order complex) of the tournament.
# TT = transitive triples = simplices of the flag complex.
# DT = doubly-transitive 4-paths = 3-simplices of the flag complex.
# Is H_2 of the flag complex always 0 for tournaments?

# The flag complex of a tournament = the simplicial complex where
# simplices are transitive subsets (totally ordered by the tournament).
# This is the ORDER COMPLEX of the tournament's partial order
# (which is actually a total pre-order — well, for tournaments it's
# more nuanced since tournaments aren't transitive in general).

# For a TRANSITIVE tournament, the flag complex = full simplex = contractible.
# For general tournaments, the flag complex can have nontrivial homology.

n = 5
flag_h2_count = Counter()

for A in all_tournaments_gen(n):
    # Build flag (simplicial) complex
    # 0-simplices: vertices
    # 1-simplices: edges (all of them)
    # 2-simplices: transitive triples
    # 3-simplices: transitive quadruples

    tt_paths = [(a,b,c) for a in range(n) for b in range(n) for c in range(n)
                if a!=b and b!=c and a!=c and A[a][b] and A[b][c] and A[a][c]]

    # Transitive quadruples: (a,b,c,d) with a→b→c→d, a→c, a→d, b→d
    tq = [(a,b,c,d) for a in range(n) for b in range(n) for c in range(n) for d in range(n)
          if len({a,b,c,d}) == 4 and A[a][b] and A[b][c] and A[c][d]
          and A[a][c] and A[a][d] and A[b][d]]

    # ∂_2 on TT: (b,c)-(a,c)+(a,b) where all are edges
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
    e_idx = {e: i for i, e in enumerate(edges)}

    bd2 = np.zeros((len(edges), len(tt_paths)))
    for j, (a,b,c) in enumerate(tt_paths):
        bd2[e_idx[(b,c)], j] += 1
        bd2[e_idx[(a,c)], j] -= 1
        bd2[e_idx[(a,b)], j] += 1

    rank_bd2 = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_bd2 = len(tt_paths) - rank_bd2

    # ∂_3 on TQ: (b,c,d)-(a,c,d)+(a,b,d)-(a,b,c)
    tt_idx = {p: i for i, p in enumerate(tt_paths)}
    bd3 = np.zeros((len(tt_paths), len(tq)))
    for j, (a,b,c,d) in enumerate(tq):
        if (b,c,d) in tt_idx: bd3[tt_idx[(b,c,d)], j] += 1
        if (a,c,d) in tt_idx: bd3[tt_idx[(a,c,d)], j] -= 1
        if (a,b,d) in tt_idx: bd3[tt_idx[(a,b,d)], j] += 1
        if (a,b,c) in tt_idx: bd3[tt_idx[(a,b,c)], j] -= 1

    rank_bd3 = np.linalg.matrix_rank(bd3, tol=1e-8)
    h2_flag = ker_bd2 - rank_bd3

    flag_h2_count[h2_flag] += 1

print(f"  n={n}: Flag complex H_2 distribution: {dict(sorted(flag_h2_count.items()))}")
print(f"  Flag H_2 = 0 for all? {'YES' if set(flag_h2_count.keys()) == {0} else 'NO'}")

print("\nDone.")
