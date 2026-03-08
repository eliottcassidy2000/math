#!/usr/bin/env python3
"""
β_2 = 0 PROOF VIA CONE CONSTRUCTION

Idea: Given a 2-cycle u ∈ ker(∂_2) ∩ Ω_2, construct v ∈ Ω_3 with ∂_3(v) = u
using a "cone from a source vertex" approach.

KEY INSIGHT: For a SOURCE vertex v (out-degree n-1), every 2-path (a,b,c)
with v→a, v→b gives a DT path (v,a,b,c). Since v beats everyone, this works
for ALL 2-paths.

∂_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)

So: ∂_3(h_v(u)) = u - u_v, where u_v is supported on (v,*,*) paths.
And u_v is still a 2-cycle in Ω_2.

Now iterate: can we kill u_v using DT paths?

For (v,x,y), need (v,x,y,w) DT: requires y→w AND x→w.
Or (z,v,x,y) DT: requires z→v (impossible if v is source).

So the second cone adds w at the END.

TEST: Does this iterated cone construction always work?
"""
import numpy as np
from itertools import combinations
import sys, time
from collections import Counter
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

def cone_reduce(u_coeffs, a2_list, A, n, source_v):
    """Given 2-cycle u = Σ u_coeffs[i] * a2_list[i], compute the residual
    after coning from source_v.

    Returns (3-chain coefficients for DT paths used, residual 2-cycle u_v).
    """
    v = source_v
    a2_idx = {p: i for i, p in enumerate(a2_list)}

    residual = np.copy(u_coeffs)

    # For each term α * (a,b,c) where v beats a and b:
    # Add α * (v,a,b,c) to the 3-chain
    # This subtracts α * (a,b,c) from residual
    # and adds α * (-(v,b,c) + (v,a,c) - (v,a,b)) to residual

    dt_paths = []  # (coefficient, (v,a,b,c))

    for i, p in enumerate(a2_list):
        if abs(u_coeffs[i]) < 1e-12:
            continue
        a, b, c = p
        if a == v or b == v or c == v:
            continue  # skip paths involving v

        if A[v][a] == 1 and A[v][b] == 1:
            alpha = u_coeffs[i]
            dt_paths.append((alpha, (v, a, b, c)))
            residual[i] -= alpha

            # Add contributions from other faces
            if (v, b, c) in a2_idx:
                residual[a2_idx[(v,b,c)]] += alpha  # -(-(v,b,c)) = +(v,b,c)... wait
                # ∂(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)
                # After subtracting α*(a,b,c): residual gets -α*(-(v,b,c) + (v,a,c) - (v,a,b))
                # = +α*(v,b,c) - α*(v,a,c) + α*(v,a,b)
                # Wait no: ∂(h_v(u)) = u - u_v, so u_v = u - ∂(h_v(u))
                # ∂(α*(v,a,b,c)) = α*[(a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)]
                # residual = u - ∂(chain) = u - Σ α*(a,b,c) + Σ α*(v,b,c) - Σ α*(v,a,c) + Σ α*(v,a,b)

    # Let me redo this more carefully
    residual = np.copy(u_coeffs)
    for alpha, path in dt_paths:
        _, a, b, c = path
        # Subtract ∂(alpha * (v,a,b,c)) from residual
        if (a,b,c) in a2_idx:
            residual[a2_idx[(a,b,c)]] -= alpha
        if (v,b,c) in a2_idx:
            residual[a2_idx[(v,b,c)]] += alpha  # -(-alpha) = +alpha
        if (v,a,c) in a2_idx:
            residual[a2_idx[(v,a,c)]] -= alpha  # -(+alpha) = -alpha
        if (v,a,b) in a2_idx:
            residual[a2_idx[(v,a,b)]] += alpha  # -(-alpha) = +alpha

    return dt_paths, residual

def test_cone_reduction(A, n, verbose=False):
    """Test cone reduction on all 2-cycles."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_list = [tuple(p) for p in a2]

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        return None

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        return None

    # Get kernel basis in A_2 coords
    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    # Find source vertex (highest out-degree)
    out_deg = [sum(A[i]) for i in range(n)]
    source = max(range(n), key=lambda v: out_deg[v])

    results = []
    for i in range(ker_dim):
        u = ker_basis[i]
        norm_u = np.linalg.norm(u)

        # First cone reduction from source
        dt1, r1 = cone_reduce(u, a2_list, A, n, source)
        norm_r1 = np.linalg.norm(r1)

        # Check: is residual supported only on (source, *, *) paths?
        non_source_support = sum(abs(r1[j]) > 1e-10
            for j, p in enumerate(a2_list)
            if p[0] != source)

        # Second cone: pick the path that doesn't start with source
        # and try to cone from the END instead
        # (v,x,y) → (v,x,y,w) where y→w and x→w

        results.append({
            'norm_u': norm_u,
            'norm_r1': norm_r1,
            'non_source': non_source_support,
            'source_deg': out_deg[source],
        })

        if verbose:
            print(f"\n  Cycle #{i+1}: ||u|| = {norm_u:.4f}")
            print(f"    After 1st cone from v={source} (deg={out_deg[source]}): ||residual|| = {norm_r1:.4f}")
            print(f"    Non-source support in residual: {non_source_support}")
            if non_source_support > 0:
                # Show what's left
                for j, p in enumerate(a2_list):
                    if abs(r1[j]) > 1e-10 and p[0] != source:
                        print(f"      {r1[j]:+.4f} * {p}")

    return results

# ===== Test at n=5 =====
print("=" * 70)
print("CONE REDUCTION FOR β_2 = 0")
print("=" * 70)

# Test on specific tournaments
print("\n--- Example tournaments at n=5 ---")
count = 0
found = 0
for A in all_tournaments_gen(5):
    count += 1
    result = test_cone_reduction(A, 5)
    if result is None:
        continue
    found += 1
    if found <= 3:
        out_deg = [sum(A[i]) for i in range(5)]
        print(f"\nTournament #{count}, out-degrees = {out_deg}")
        test_cone_reduction(A, 5, verbose=True)

# ===== Statistics =====
print(f"\n\n{'='*70}")
print("CONE REDUCTION STATISTICS (n=5)")
print("="*70)

zero_after_1cone = 0
nonzero_after_1cone = 0
all_source_support = 0
total = 0

for A in all_tournaments_gen(5):
    result = test_cone_reduction(A, 5)
    if result is None:
        continue
    total += 1
    for r in result:
        if r['norm_r1'] < 1e-10:
            zero_after_1cone += 1
        else:
            nonzero_after_1cone += 1
            if r['non_source'] == 0:
                all_source_support += 1

print(f"Total cycles: {zero_after_1cone + nonzero_after_1cone}")
print(f"  Killed by 1st cone: {zero_after_1cone}")
print(f"  Residual nonzero: {nonzero_after_1cone}")
print(f"    Supported only on source-paths: {all_source_support}")
print(f"    Has non-source support: {nonzero_after_1cone - all_source_support}")

# ===== Can residuals be killed by a 2nd cone? =====
print(f"\n\n{'='*70}")
print("ITERATED CONE: 2nd step from END vertex")
print("="*70)

def cone_from_end(r1_coeffs, a2_list, A, n, source):
    """For residual r1 supported on (source,x,y) paths,
    try to cone from the END: (source,x,y) → (source,x,y,w)."""
    a2_idx = {p: i for i, p in enumerate(a2_list)}

    dt_paths = []
    residual = np.copy(r1_coeffs)

    # For (source,x,y), need w with y→w AND x→w
    for i, p in enumerate(a2_list):
        if abs(r1_coeffs[i]) < 1e-12:
            continue
        if p[0] != source:
            continue
        _, x, y = p
        alpha = r1_coeffs[i]

        # Find w beaten by both x and y
        candidates = [w for w in range(n) if w != source and w != x and w != y
                      and A[y][w] == 1 and A[x][w] == 1]
        if candidates:
            w = candidates[0]  # pick first
            dt_paths.append((alpha, (source, x, y, w)))
            # ∂(s,x,y,w) = (x,y,w) - (s,y,w) + (s,x,w) - (s,x,y)
            if (x,y,w) in a2_idx:
                residual[a2_idx[(x,y,w)]] -= alpha
            if (source,y,w) in a2_idx:
                residual[a2_idx[(source,y,w)]] += alpha
            if (source,x,w) in a2_idx:
                residual[a2_idx[(source,x,w)]] -= alpha
            residual[i] += alpha  # -(-(s,x,y)) for the d_3 face

    return dt_paths, residual

killed_2step = 0
survived_2step = 0

for A in all_tournaments_gen(5):
    a1 = enumerate_allowed_paths(A, 5, 1)
    a2 = enumerate_allowed_paths(A, 5, 2)
    a2_list = [tuple(p) for p in a2]

    om2 = compute_omega_basis(A, 5, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    out_deg = [sum(A[i]) for i in range(5)]
    source = max(range(5), key=lambda v: out_deg[v])

    for i in range(ker_dim):
        u = ker_basis[i]
        dt1, r1 = cone_reduce(u, a2_list, A, 5, source)
        if np.linalg.norm(r1) < 1e-10:
            killed_2step += 1
            continue
        dt2, r2 = cone_from_end(r1, a2_list, A, 5, source)
        if np.linalg.norm(r2) < 1e-10:
            killed_2step += 1
        else:
            survived_2step += 1

print(f"After 2-step cone: killed={killed_2step}, survived={survived_2step}")

# ===== What if we try ALL vertices as source, not just max out-degree? =====
print(f"\n\n{'='*70}")
print("TRY ALL VERTICES AS CONE SOURCE")
print("="*70)

total_cycles = 0
killed_best = 0

for A in all_tournaments_gen(5):
    a1 = enumerate_allowed_paths(A, 5, 1)
    a2 = enumerate_allowed_paths(A, 5, 2)
    a2_list = [tuple(p) for p in a2]

    om2 = compute_omega_basis(A, 5, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2_list, [tuple(p) for p in a1])
    bd2_om = bd2 @ om2
    rank2 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_dim = dim_om2 - rank2
    if ker_dim == 0:
        continue

    U, S, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    ker_basis = (om2 @ Vt[rank2:, :].T).T

    for i in range(ker_dim):
        total_cycles += 1
        u = ker_basis[i]

        best_norm = np.linalg.norm(u)
        for v in range(5):
            dt1, r1 = cone_reduce(u, a2_list, A, 5, v)
            norm1 = np.linalg.norm(r1)
            if norm1 < best_norm:
                best_norm = norm1

        if best_norm < 1e-10:
            killed_best += 1

print(f"Total cycles: {total_cycles}")
print(f"Killed by single cone from best vertex: {killed_best}")
print(f"Need multi-step: {total_cycles - killed_best}")

# ===== Check if single-vertex cone EVER suffices =====
# For a tournament with a Hamiltonian path, the source vertex of that path
# beats all others. But not every tournament has such a vertex.
print(f"\n--- Out-degree distribution for source vertex ---")
deg_dist = Counter()
for A in all_tournaments_gen(5):
    out_deg = [sum(A[i]) for i in range(5)]
    max_deg = max(out_deg)
    deg_dist[max_deg] += 1

print(f"  Max out-degree: count")
for d in sorted(deg_dist.keys()):
    print(f"    {d}: {deg_dist[d]}")
# At n=5: max out-degree is at most 4 (source), at least 3 (regular tournament)

print("\nDone.")
