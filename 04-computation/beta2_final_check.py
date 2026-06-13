#!/usr/bin/env python3
"""
β_2 = 0: FINAL VERIFICATION

Check that im(∂_3|Ω_3) = ker(∂_2|Ω_2) using the CORRECT definition.

Ω_3 = {u ∈ A_3 : ∂u ∈ A_2}
∂_3|Ω_3: Ω_3 → Ω_2 (well-defined since ∂²=0 => ∂(Ω_3) ⊆ Ω_2)

We compute:
1. ker(∂_2|Ω_2) using boundary restricted to transitive triples
2. im(∂_3|Ω_3) using compute_omega_basis for Ω_3, then applying ∂_3

Also compare with:
3. im(projection of ∂(DT paths) to Ω_2) — the "DT sufficiency" test

If (2) = (1) and (3) = (1), both methods give β_2 = 0.
"""
import numpy as np
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

n = 5
print("=" * 70)
print(f"β_2 VERIFICATION at n={n}")
print("=" * 70)

all_pass = True
for t_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    # Ω_2 = transitive triples
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    if len(tt) == 0:
        continue

    # ker(∂_2|Ω_2)
    a1_idx = {tuple(p): i for i, p in enumerate(a1)}
    bd2 = np.zeros((len(a1), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2[a1_idx[(b,c)], j] += 1
        bd2[a1_idx[(a,c)], j] -= 1
        bd2[a1_idx[(a,b)], j] += 1
    rank2 = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(tt) - rank2
    if ker_dim == 0:
        continue

    # Ω_3 (correct definition)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # ∂_3|Ω_3 → A_2, then check it lands in Ω_2
    a2_tuples = [tuple(p) for p in a2]
    a2_idx = {t: i for i, t in enumerate(a2_tuples)}
    a3_tuples = [tuple(p) for p in a3]

    if dim_om3 > 0:
        bd3_full = np.zeros((len(a2_tuples), len(a3_tuples)))
        for j, path in enumerate(a3_tuples):
            v0,v1,v2,v3 = path
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            signs = [1, -1, 1, -1]
            for face, sign in zip(faces, signs):
                if face in a2_idx:
                    bd3_full[a2_idx[face], j] += sign

        # Apply to Ω_3 basis
        bd3_omega = bd3_full @ om3  # columns in A_2

        # Verify: ∂(Ω_3) lands in Ω_2
        tt_set = set(tt)
        non_tt_indices = [i for i, t in enumerate(a2_tuples) if t not in tt_set]
        max_non_tt = np.max(np.abs(bd3_omega[non_tt_indices, :])) if non_tt_indices else 0
        if max_non_tt > 1e-8:
            print(f"  T#{t_idx}: ∂(Ω_3) NOT in Ω_2! max={max_non_tt:.4e}")
            all_pass = False
            continue

        # Project to Ω_2 coordinates
        tt_indices = [a2_idx[t] for t in tt]
        bd3_tt = bd3_omega[tt_indices, :]
        im_rank = np.linalg.matrix_rank(bd3_tt, tol=1e-8)
    else:
        im_rank = 0

    beta2 = ker_dim - im_rank
    if beta2 != 0:
        print(f"  T#{t_idx}: β_2 = {beta2} ≠ 0 !")
        all_pass = False

print(f"\nn={n}: All tournaments β_2 = 0? {all_pass}")

# ===== Now check: does Ω_3 = span(Ω_3 ∩ DT paths) + cancellation chains? =====
print(f"\n{'='*70}")
print("Ω_3 STRUCTURE: INDIVIDUAL DT vs CANCELLATION CHAINS")
print("="*70)

from collections import Counter
om3_structure = Counter()

for t_idx, A in enumerate(all_tournaments_gen(n)):
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a3_tuples = [tuple(p) for p in a3]

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    # DT paths
    dt_indices = [j for j, p in enumerate(a3) if A[p[0]][p[2]]==1 and A[p[1]][p[3]]==1]
    n_dt = len(dt_indices)

    # 4-clique DT paths (a→d too)
    clique_indices = [j for j in dt_indices if A[a3[j][0]][a3[j][3]]==1]
    n_clique = len(clique_indices)

    if dim_om3 > 0:
        # How many Ω_3 basis vectors are pure DT (single DT path)?
        n_pure_dt = 0
        for col in range(dim_om3):
            vec = om3[:, col]
            support = [j for j in range(len(a3_tuples)) if abs(vec[j]) > 1e-8]
            if len(support) == 1 and support[0] in dt_indices:
                n_pure_dt += 1

        om3_structure[(dim_om3, n_dt, n_clique, n_pure_dt)] += 1

print("(dim_Ω_3, #DT, #clique, #pure_DT_basis_vecs): count")
for key in sorted(om3_structure.keys()):
    d, ndt, ncl, npure = key
    print(f"  ({d}, DT={ndt}, clique={ncl}, pure={npure}): {om3_structure[key]}")

# ===== Check: are 4-clique paths always in Ω_3? =====
print(f"\n{'='*70}")
print("ARE 4-CLIQUE DT PATHS ALWAYS INDIVIDUALLY IN Ω_3?")
print("="*70)

always_in = True
for A in all_tournaments_gen(n):
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_set = set(tuple(p) for p in a2)

    for p in a3:
        v0,v1,v2,v3 = p
        if A[v0][v2] and A[v1][v3] and A[v0][v3]:
            # 4-clique DT path — check all faces in A_2
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            for face in faces:
                if face not in a2_set:
                    print(f"  4-clique DT path {tuple(p)} has face {face} NOT in A_2!")
                    always_in = False

print(f"4-clique DT paths always in Ω_3? {always_in}")

# ===== Check: are non-clique DT paths (d→a) ever individually in Ω_3? =====
print(f"\nNon-clique DT paths individually in Ω_3?")
nonclique_in_om3 = 0
nonclique_not_in_om3 = 0
for A in all_tournaments_gen(n):
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_set = set(tuple(p) for p in a2)

    for p in a3:
        v0,v1,v2,v3 = p
        if A[v0][v2] and A[v1][v3] and not A[v0][v3]:
            # Non-clique DT: d→a
            faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
            all_in = all(face in a2_set for face in faces)
            if all_in:
                nonclique_in_om3 += 1
            else:
                nonclique_not_in_om3 += 1

print(f"  In Ω_3: {nonclique_in_om3}")
print(f"  NOT in Ω_3: {nonclique_not_in_om3}")

print("\nDone.")
