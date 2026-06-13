#!/usr/bin/env python3
"""Detailed analysis of DT failures at n=6.

At n=6, 3% of (tournament, 2-cycle) pairs cannot be filled by DT alone.
This script analyzes exactly WHAT the missing fillings look like and
HOW cancellation chains provide them.

Goal: identify the algebraic mechanism that guarantees β_2=0 at n=6.
"""
import numpy as np
from itertools import combinations, permutations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    boundary_coeffs, path_betti_numbers
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

def is_DT(A, a, b, c, d):
    return (A[a][b] == 1 and A[b][c] == 1 and A[c][d] == 1 and
            A[a][c] == 1 and A[b][d] == 1)

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def t3_count(A, n):
    return sum(1 for a, b, c in combinations(range(n), 3)
               if (A[a][b] and A[b][c] and A[c][a]) or
                  (A[b][a] and A[a][c] and A[c][b]))

# ===== Analyze DT failures at n=6 =====
print("=" * 70)
print("DT FAILURE ANALYSIS AT n=6")
print("=" * 70)

n = 6
dt_gap_tournaments = []
total = 0
dt_sufficient = 0
dt_failures = 0

for tidx, A in enumerate(all_tournaments_gen(n)):
    total += 1
    if total % 5000 == 0:
        print(f"  ... {total}", flush=True)

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    # Kernel of ∂_2|Ω_2
    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rank2 = sum(s > 1e-8 for s in S2)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    # DT paths
    dt_idx = [i for i, p in enumerate(a3_t) if is_DT(A, *p)]
    if not dt_idx:
        dt_failures += 1
        dt_gap_tournaments.append((tidx, A, ker_dim, 0, ker_dim))
        continue

    # Image of DT paths
    bd3 = build_full_boundary_matrix(a3_t, a2_t)
    bd3_dt = bd3[:, dt_idx]

    # Does im(∂_3|DT) ⊇ ker(∂_2|Ω_2)?
    _, S2full, Vt2full = np.linalg.svd(bd2_om)
    ker_basis = om2 @ Vt2full[rank2:].T  # columns are kernel generators in A_2 coords

    # Check: is each kernel generator in im(bd3_dt)?
    rank_dt = np.linalg.matrix_rank(bd3_dt, tol=1e-8)
    combined = np.hstack([bd3_dt, ker_basis])
    rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)

    gap = rank_combined - rank_dt
    if gap > 0:
        dt_failures += 1
        dt_gap_tournaments.append((tidx, [row[:] for row in A], ker_dim, rank_dt, gap))
    else:
        dt_sufficient += 1

print(f"\n  Total tournaments: {total}")
print(f"  DT sufficient: {dt_sufficient}")
print(f"  DT failures: {dt_failures}")

# Analyze the failures
print(f"\n--- DT failure details ---")
gap_stats = Counter()
score_stats = Counter()
t3_stats = Counter()

for tidx, A, ker_dim, rank_dt, gap in dt_gap_tournaments:
    ss = score_seq(A, n)
    t3 = t3_count(A, n)
    gap_stats[gap] += 1
    score_stats[ss] += 1
    t3_stats[t3] += 1

print(f"  Gap distribution: {dict(sorted(gap_stats.items()))}")
print(f"  Score sequences: {dict(sorted(score_stats.items()))}")
print(f"  t3 values: {dict(sorted(t3_stats.items()))}")

# ===== Deep dive into a specific failure =====
print(f"\n\n{'='*70}")
print("DEEP DIVE: FIRST DT FAILURE")
print("=" * 70)

if dt_gap_tournaments:
    tidx, A, ker_dim, rank_dt, gap = dt_gap_tournaments[0]

    print(f"Tournament #{tidx}, score={score_seq(A, n)}, t3={t3_count(A, n)}")
    print(f"  ker(∂_2|Ω_2) dim = {ker_dim}, DT image rank = {rank_dt}, gap = {gap}")

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

    print(f"  dim(Ω_2) = {dim_om2}, dim(Ω_3) = {dim_om3}")

    # DT and non-DT in Ω_3
    dt_idx = [i for i, p in enumerate(a3_t) if is_DT(A, *p)]
    non_dt_idx = [i for i, p in enumerate(a3_t) if not is_DT(A, *p)]

    print(f"  |DT| = {len(dt_idx)}, |A_3\\DT| = {len(non_dt_idx)}")

    # Ω_3 elements that are NOT DT
    # Each column of om3 is an Ω_3 basis element in A_3 coordinates
    # Check which basis elements have support on non-DT paths
    pure_dt_count = 0
    mixed_count = 0
    pure_nondt_count = 0

    for j in range(dim_om3):
        col = om3[:, j]
        has_dt = any(abs(col[i]) > 1e-8 for i in dt_idx)
        has_nondt = any(abs(col[i]) > 1e-8 for i in non_dt_idx)
        if has_dt and has_nondt:
            mixed_count += 1
        elif has_dt:
            pure_dt_count += 1
        elif has_nondt:
            pure_nondt_count += 1

    print(f"  Ω_3 basis: {pure_dt_count} pure DT, {mixed_count} mixed, {pure_nondt_count} pure non-DT")

    # Find the unfilled kernel direction
    bd2 = build_full_boundary_matrix(a2_t, a1_t)
    bd2_om = bd2 @ om2
    _, S2full, Vt2full = np.linalg.svd(bd2_om)
    rank2 = sum(s > 1e-8 for s in S2full)
    ker_basis = om2 @ Vt2full[rank2:].T

    bd3 = build_full_boundary_matrix(a3_t, a2_t)
    bd3_dt = bd3[:, dt_idx]

    # Find which kernel generator is NOT in im(DT)
    for k in range(ker_dim):
        z = ker_basis[:, k]
        # Project z onto im(bd3_dt): residual = z - bd3_dt @ (bd3_dt \ z)
        sol, res, _, _ = np.linalg.lstsq(bd3_dt, z, rcond=None)
        residual = z - bd3_dt @ sol
        resid_norm = np.linalg.norm(residual)
        if resid_norm > 1e-6:
            print(f"\n  Kernel generator {k}: NOT in im(DT), residual norm = {resid_norm:.6f}")

            # Show the residual
            z_support = [(a2_t[i], z[i]) for i in range(len(a2_t)) if abs(z[i]) > 1e-8]
            print(f"    Cycle z (support={len(z_support)} paths):")
            for p, c in z_support[:10]:
                print(f"      {c:+.6f} * {p}")
            if len(z_support) > 10:
                print(f"      ... ({len(z_support)} total)")

            res_support = [(a2_t[i], residual[i]) for i in range(len(a2_t))
                           if abs(residual[i]) > 1e-8]
            print(f"    Residual (support={len(res_support)} paths):")
            for p, c in res_support[:10]:
                print(f"      {c:+.6f} * {p}")

            # NOW fill with full Ω_3
            bd3_om = bd3 @ om3
            sol_full, _, _, _ = np.linalg.lstsq(bd3_om, z, rcond=None)
            residual_full = z - bd3_om @ sol_full
            print(f"    Full Ω_3 filling residual: {np.linalg.norm(residual_full):.6e}")

            # What Ω_3 elements are needed?
            w_a3 = om3 @ sol_full
            w_support = [(i, a3_t[i], w_a3[i]) for i in range(len(a3_t))
                         if abs(w_a3[i]) > 1e-8]
            n_dt_used = sum(1 for i, p, _ in w_support if is_DT(A, *p))
            n_nondt_used = len(w_support) - n_dt_used
            print(f"    Filling: {len(w_support)} paths ({n_dt_used} DT + {n_nondt_used} non-DT)")

            # Show non-DT paths in filling
            print(f"    Non-DT paths in filling:")
            for i, p, c in w_support:
                if not is_DT(A, *p):
                    # What faces are bad?
                    bad_faces = []
                    for sign, face in boundary_coeffs(p):
                        if len(set(face)) == len(face):
                            face_t = tuple(face)
                            if face_t not in set(a2_t):
                                bad_faces.append(face_t)
                    print(f"      {c:+.6f} * {p}  bad faces: {bad_faces}")

            break  # Just first unfilled generator

# ===== Part 2: For failure tournaments, analyze the CANCELLATION structure =====
print(f"\n\n{'='*70}")
print("CANCELLATION CHAIN STRUCTURE IN DT FAILURES")
print("=" * 70)

if dt_gap_tournaments:
    # Take a few failure tournaments
    for idx_f, (tidx, A, ker_dim, rank_dt, gap) in enumerate(dt_gap_tournaments[:3]):
        print(f"\n--- Failure #{idx_f}: T#{tidx}, gap={gap} ---")

        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)
        a2_t = [tuple(p) for p in a2]
        a3_t = [tuple(p) for p in a3]
        a2_set = set(a2_t)

        # Find non-DT 3-paths and their bad faces
        non_dt_paths = []
        for p in a3_t:
            if not is_DT(A, *p):
                bad_faces = []
                for sign, face in boundary_coeffs(tuple(p)):
                    face_t = tuple(face)
                    if len(set(face_t)) == len(face_t) and face_t not in a2_set:
                        bad_faces.append((sign, face_t))
                non_dt_paths.append((tuple(p), bad_faces))

        print(f"  Non-DT 3-paths: {len(non_dt_paths)}")

        # Group by bad face
        bf_groups = defaultdict(list)
        for p, bfs in non_dt_paths:
            for sign, bf in bfs:
                bf_groups[bf].append((p, sign))

        print(f"  Bad face groups:")
        for bf in sorted(bf_groups):
            group = bf_groups[bf]
            print(f"    {bf}: {len(group)} paths")
            for p, sign in group[:3]:
                print(f"      {'+' if sign > 0 else '-'} {p}")

        # Key insight: for paths p1, p2 sharing bad face bf,
        # (sign1)*p1 - (sign2)*p2 has bf cancelled.
        # If all OTHER faces of p1 and p2 are allowed, this difference is in Ω_3!
        # Check this:
        print(f"\n  Cancellation pairs (sharing one bad face, all others allowed):")
        cancel_pairs = 0
        for bf, group in bf_groups.items():
            for i in range(len(group)):
                for j in range(i+1, len(group)):
                    p1, s1 = group[i]
                    p2, s2 = group[j]

                    # Difference: s1*p1 - s2*p2 (so bf cancels)
                    # Check other faces
                    diff_faces = set()
                    for sign, face in boundary_coeffs(p1):
                        face_t = tuple(face)
                        if face_t != bf and len(set(face_t)) == len(face_t):
                            diff_faces.add(face_t)
                    for sign, face in boundary_coeffs(p2):
                        face_t = tuple(face)
                        if face_t != bf and len(set(face_t)) == len(face_t):
                            diff_faces.add(face_t)

                    all_allowed = all(f in a2_set for f in diff_faces)
                    if all_allowed:
                        cancel_pairs += 1
                        if cancel_pairs <= 3:
                            print(f"    {p1} - {p2} (bad face {bf})")

        print(f"  Total cancellation pairs: {cancel_pairs}")

# ===== Part 3: Is β_2 related to the simplicial homology of the flag complex? =====
print(f"\n\n{'='*70}")
print("DIMENSION COUNTING ANALYSIS")
print("=" * 70)

n = 6
dim_data = []
for tidx, A in enumerate(all_tournaments_gen(n)):
    if tidx % 5000 == 0:
        print(f"  ... {tidx}", flush=True)

    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a4 = enumerate_allowed_paths(A, n, 4)

    a1_t = [tuple(p) for p in a1]
    a2_t = [tuple(p) for p in a2]
    a3_t = [tuple(p) for p in a3]
    a4_t = [tuple(p) for p in a4]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    om4 = compute_omega_basis(A, n, 4, a4, a3)

    d2 = om2.shape[1] if om2.ndim == 2 else 0
    d3 = om3.shape[1] if om3.ndim == 2 else 0
    d4 = om4.shape[1] if om4.ndim == 2 else 0

    # Boundary ranks
    if d2 > 0:
        bd2 = build_full_boundary_matrix(a2_t, a1_t)
        bd2_om = bd2 @ om2
        S = np.linalg.svd(bd2_om, compute_uv=False)
        r2 = sum(s > 1e-8 for s in S)
    else:
        r2 = 0

    if d3 > 0:
        bd3 = build_full_boundary_matrix(a3_t, a2_t)
        bd3_om = bd3 @ om3
        S = np.linalg.svd(bd3_om, compute_uv=False)
        r3 = sum(s > 1e-8 for s in S)
    else:
        r3 = 0

    if d4 > 0:
        bd4 = build_full_boundary_matrix(a4_t, a3_t)
        bd4_om = bd4 @ om4
        S = np.linalg.svd(bd4_om, compute_uv=False)
        r4 = sum(s > 1e-8 for s in S)
    else:
        r4 = 0

    ker2 = d2 - r2
    ker3 = d3 - r3

    # β_2 = ker2 - r3, should be 0
    beta2 = ker2 - r3

    # β_3 = ker3 - r4
    beta3 = ker3 - r4

    ss = score_seq(A, n)
    t3 = t3_count(A, n)
    dim_data.append((d2, d3, d4, r2, r3, r4, ker2, ker3, beta2, beta3, ss, t3))

# Verify β_2 = 0
beta2_nonzero = sum(1 for x in dim_data if x[8] != 0)
print(f"  β_2 ≠ 0: {beta2_nonzero} / {len(dim_data)}")

# Check the rank formula
# ker2 = d2 - r2 = d2 - (C(n,2) - n + 1 - β_1)
# β_1 = C(n,2) - n + 1 - r2
# So ker2 = d2 - r2

# The key: β_2 = 0 iff r3 = ker2 = d2 - r2
# i.e., rank(∂_3|Ω_3) = dim(Ω_2) - rank(∂_2|Ω_2)

# Check: is r3 always = ker2?
r3_eq_ker2 = sum(1 for x in dim_data if x[6] == x[4])  # ker2 == r3
print(f"  r3 = ker2: {r3_eq_ker2} / {len(dim_data)}")

# Distribution of key dimensions
print(f"\nDimension distributions:")
d2_vals = Counter(x[0] for x in dim_data)
d3_vals = Counter(x[1] for x in dim_data)
ker2_vals = Counter(x[6] for x in dim_data)
r3_vals = Counter(x[4] for x in dim_data)

print(f"  dim(Ω_2): {dict(sorted(d2_vals.items()))}")
print(f"  dim(Ω_3): {dict(sorted(d3_vals.items()))}")
print(f"  ker(∂_2|Ω_2): {dict(sorted(ker2_vals.items()))}")
print(f"  rank(∂_3|Ω_3): {dict(sorted(r3_vals.items()))}")

# Check: is dim(Ω_3) - ker(∂_3) = ker(∂_2)?
# i.e., r3 = ker2, which is β_2 = 0
print(f"\n  dim(Ω_3) - ker(∂_3) = rank(∂_3) = r3")
print(f"  ker(∂_2) = ker2")
print(f"  β_2 = ker2 - r3 = {dict(Counter(x[8] for x in dim_data))}")

# Relationship: d3 = r3 + ker3
# And ker3 = r4 + β_3
# So d3 = r3 + r4 + β_3
# And β_2 = ker2 - r3 = 0 means r3 = ker2 = d2 - r2

# So the exact constraint is:
# dim(Ω_3) ≥ ker(∂_2|Ω_2) = dim(Ω_2) - rank(∂_2|Ω_2)
# AND rank(∂_3|Ω_3) = ker(∂_2|Ω_2) exactly

# Is dim(Ω_3) always ≥ ker(∂_2|Ω_2)?
deficit = [x[1] - x[6] for x in dim_data]  # d3 - ker2
print(f"\n  dim(Ω_3) - ker(∂_2): min={min(deficit)}, max={max(deficit)}")
print(f"  Distribution: {dict(sorted(Counter(deficit).items()))}")

print("\nDone.")
