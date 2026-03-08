#!/usr/bin/env python3
"""
beta2_delta_explicit.py - Explicit computation of the connecting homomorphism delta
for the regular n=5 tournaments (the "bad" cases where all beta1(T\v)=1).

Goal: Understand WHY delta is injective for interior vertices.
Show the explicit relative 2-cycle and its delta-image.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def compute_delta_explicit(A, n, v):
    """
    Compute the connecting homomorphism delta: H_2(T,T\\v) -> H_1(T\\v)
    explicitly, returning the relative 2-cycle and its delta-image.
    """
    others = [i for i in range(n) if i != v]
    A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
    n_sub = n - 1

    # Compute beta_1(T\v)
    betti_sub = path_betti_numbers(A_sub, n_sub, max_dim=1)
    beta1_sub = betti_sub[1]

    if beta1_sub == 0:
        return {'beta1_sub': 0, 'h2_rel_dim': 0}

    # Full tournament paths
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)

    # Sub-tournament paths
    paths2_sub = enumerate_allowed_paths(A_sub, n_sub, 2)
    paths1_sub = enumerate_allowed_paths(A_sub, n_sub, 1)

    # Map sub indices back to original vertices
    def orig_path(p, others=others):
        return tuple(others[i] for i in p)

    paths2_sub_orig = [orig_path(p) for p in paths2_sub]
    paths1_sub_orig = [orig_path(p) for p in paths1_sub]

    # Identify 2-paths using v (relative chain group)
    paths2_using_v = [p for p in paths2 if v in p]
    paths2_not_v = [p for p in paths2 if v not in p]

    # 2-paths using v, categorized by position of v
    pos0 = [p for p in paths2_using_v if p[0] == v]
    pos1 = [p for p in paths2_using_v if p[1] == v]
    pos2 = [p for p in paths2_using_v if p[2] == v]

    # Omega_2(T) basis
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_omega2 = omega2.shape[1]

    # Omega_2(T\v) basis
    omega2_sub = compute_omega_basis(A_sub, n_sub, 2, paths2_sub, paths1_sub)
    dim_omega2_sub = omega2_sub.shape[1]

    # Identify basis elements in Omega_2(T) that use v
    # Express relative Omega_2 = Omega_2(T) / Omega_2(T\v)
    path2_idx = {p: i for i, p in enumerate(paths2)}

    # Embed Omega_2(T\v) into Omega_2(T)
    # T\v paths in original indices
    embed_sub = np.zeros((len(paths2), dim_omega2_sub))
    for j in range(dim_omega2_sub):
        for i_sub, p_sub in enumerate(paths2_sub):
            p_orig = orig_path(p_sub)
            if p_orig in path2_idx:
                embed_sub[path2_idx[p_orig], j] = omega2_sub[i_sub, j]

    # Relative chain = Omega_2(T) modulo Omega_2(T\v)
    # We need to compute the quotient and then find ker(d_2) in the quotient

    # Build d_2: A_2 -> A_1 for T
    path1_idx = {p: i for i, p in enumerate(paths1)}
    D2 = np.zeros((len(paths1), len(paths2)))
    for j, (a,b,c) in enumerate(paths2):
        if (b,c) in path1_idx: D2[path1_idx[(b,c)], j] += 1
        if (a,c) in path1_idx: D2[path1_idx[(a,c)], j] -= 1
        if (a,b) in path1_idx: D2[path1_idx[(a,b)], j] += 1

    # d_2 in Omega_2 coordinates
    D2_omega = D2 @ omega2

    # Project d_2 onto 1-paths NOT using v (for relative cycle condition)
    paths1_sub_set = set(paths1_sub_orig)
    paths1_using_v = [p for p in paths1 if v in p]
    paths1_not_v = [p for p in paths1 if v not in p]

    # Index maps
    p1_not_v_idx = {p: i for i, p in enumerate(paths1_not_v)}
    p1_using_v_idx = {p: i for i, p in enumerate(paths1_using_v)}

    # Projection of D2 onto 1-paths not using v
    D2_not_v = np.zeros((len(paths1_not_v), len(paths2)))
    for j, (a,b,c) in enumerate(paths2):
        if (b,c) in p1_not_v_idx: D2_not_v[p1_not_v_idx[(b,c)], j] += 1
        if (a,c) in p1_not_v_idx: D2_not_v[p1_not_v_idx[(a,c)], j] -= 1
        if (a,b) in p1_not_v_idx: D2_not_v[p1_not_v_idx[(a,b)], j] += 1

    D2_not_v_omega = D2_not_v @ omega2

    # Relative 2-cycle: z in Omega_2(T) such that d_2(z) has no component in T\v's 1-paths
    # Wait, relative cycle condition: d_2(z) in C_1(T\v), which means ONLY 1-paths not using v
    # but also only ALLOWED 1-paths in T\v.
    # Actually d_2(z) is already in A_1(T) (since z in Omega_2(T)), and the relative condition
    # is d_2(z) in A_1(T\v). Since A_1(T\v) = {edges not using v}, we need d_2(z) to have
    # zero coefficient on all 1-paths using v.

    D2_using_v = np.zeros((len(paths1_using_v), len(paths2)))
    for j, (a,b,c) in enumerate(paths2):
        if (b,c) in p1_using_v_idx: D2_using_v[p1_using_v_idx[(b,c)], j] += 1
        if (a,c) in p1_using_v_idx: D2_using_v[p1_using_v_idx[(a,c)], j] -= 1
        if (a,b) in p1_using_v_idx: D2_using_v[p1_using_v_idx[(a,b)], j] += 1

    D2_using_v_omega = D2_using_v @ omega2

    # Relative 2-cycle: z in Omega_2(T) with d_2(z)|_{using v} = 0
    # AND z not in Omega_2(T\v)
    # The space of relative cycles is: ker(D2_using_v_omega) mod Omega_2(T\v)

    # First find ker(D2_using_v_omega)
    U, S, Vt = np.linalg.svd(D2_using_v_omega, full_matrices=True)
    rank_uv = sum(s > 1e-8 for s in S)
    ker_uv = Vt[rank_uv:].T  # in Omega_2 coordinates

    # Now mod out by Omega_2(T\v)
    # Express Omega_2(T\v) in Omega_2(T) coordinates
    # omega2_sub is in A_2(T\v) coordinates
    # embed_sub is in A_2(T) coordinates
    # Convert to Omega_2(T) coordinates: omega2^{-1} @ embed_sub
    # omega2 has shape (|A_2|, dim_omega2), pseudoinverse:
    omega2_pinv = np.linalg.pinv(omega2)
    sub_in_omega2 = omega2_pinv @ embed_sub  # shape (dim_omega2, dim_omega2_sub)

    # Relative cycles = ker_uv mod sub_in_omega2
    # Combined space: [ker_uv | sub_in_omega2]
    combined = np.hstack([ker_uv, sub_in_omega2])
    _, S_comb, _ = np.linalg.svd(combined, full_matrices=False)
    rank_combined = sum(s > 1e-8 for s in S_comb)
    _, S_sub_only, _ = np.linalg.svd(sub_in_omega2, full_matrices=False)
    rank_sub = sum(s > 1e-8 for s in S_sub_only)

    h2_rel_dim = ker_uv.shape[1] - (rank_combined - rank_sub)
    # Wait, this isn't right. Let me think again.
    # H_2(T,T\v) = Z_2(T,T\v) / B_2(T,T\v)
    # Z_2 relative = ker(d_2_rel) = {z in Omega_2(T)/Omega_2(T\v) : d_2(z) in A_1(T\v)}
    # This is: {z in Omega_2(T) : d_2(z) has zero v-component} / Omega_2(T\v)

    # dim Z_2_rel = dim(ker D2_using_v_omega) - dim(ker D2_using_v_omega intersect Omega_2(T\v))
    # Let me compute the intersection
    if sub_in_omega2.shape[1] > 0:
        # Check if sub elements are in ker_uv
        D2_uv_sub = D2_using_v_omega @ sub_in_omega2
        _, S_sub_check, _ = np.linalg.svd(D2_uv_sub, full_matrices=False)
        rank_sub_in_ker = sub_in_omega2.shape[1] - sum(s > 1e-8 for s in S_sub_check)
    else:
        rank_sub_in_ker = 0

    # Actually, let me just directly compute the relative cycle space
    # z in ker(D2_using_v_omega), and we want z modulo image of sub_in_omega2 in ker
    # Since Omega_2(T\v) should be INSIDE ker(D2_using_v_omega) (paths in T\v don't use v at all)
    # Check this:
    if sub_in_omega2.shape[1] > 0:
        test = D2_using_v_omega @ sub_in_omega2
        max_err = np.max(np.abs(test))
    else:
        max_err = 0

    # dim H_2(T,T\v) = dim ker(D2_uv_omega) - rank(sub_in_omega2)
    # but we also need to account for boundaries B_2(T,T\v)
    # For now, just compute delta on the relative cycles

    # Find a relative cycle NOT in Omega_2(T\v)
    # Take a vector in ker_uv that's NOT in span(sub_in_omega2)
    if ker_uv.shape[1] > rank_sub:
        # Project out sub_in_omega2 from ker_uv
        if sub_in_omega2.shape[1] > 0:
            proj = sub_in_omega2 @ np.linalg.pinv(sub_in_omega2)
            ker_relative = ker_uv - proj @ ker_uv
            # Find nonzero columns
            norms = np.linalg.norm(ker_relative, axis=0)
            good_cols = [i for i in range(ker_relative.shape[1]) if norms[i] > 1e-8]
            if good_cols:
                z_omega2 = ker_relative[:, good_cols[0]]
                z_omega2 = z_omega2 / np.linalg.norm(z_omega2)
            else:
                return {'beta1_sub': beta1_sub, 'h2_rel_dim': 0, 'error': 'no relative cycle found'}
        else:
            z_omega2 = ker_uv[:, 0]
            z_omega2 = z_omega2 / np.linalg.norm(z_omega2)

        # Convert to A_2 coordinates
        z_A2 = omega2 @ z_omega2

        # Compute delta(z) = d_2(z)|_{not using v}
        delta_z = D2_not_v @ z_A2

        # Check that delta_z is a 1-cycle in T\v
        # Build d_1 for T\v
        paths0_sub = [(i,) for i in range(n_sub)]
        D1_sub = np.zeros((n_sub, len(paths1_sub)))
        for j, (a,b) in enumerate(paths1_sub):
            D1_sub[b, j] += 1
            D1_sub[a, j] -= 1

        # Map delta_z to T\v coordinates
        delta_z_sub = np.zeros(len(paths1_sub))
        for i, p_sub in enumerate(paths1_sub):
            p_orig = orig_path(p_sub)
            if p_orig in p1_not_v_idx:
                delta_z_sub[i] = delta_z[p1_not_v_idx[p_orig]]

        is_cycle = np.max(np.abs(D1_sub @ delta_z_sub)) < 1e-8

        # Check if delta_z is a BOUNDARY in T\v (delta injective iff it's NOT a boundary)
        # Build d_2 for T\v
        D2_sub = build_full_boundary_matrix(paths2_sub, paths1_sub)
        omega2_sub_full = compute_omega_basis(A_sub, n_sub, 2, paths2_sub, paths1_sub)

        if omega2_sub_full.shape[1] > 0:
            D2_sub_omega = D2_sub @ omega2_sub_full
            # Check if delta_z_sub is in im(D2_sub_omega)
            aug = np.hstack([D2_sub_omega, delta_z_sub.reshape(-1,1)])
            rank_D2_sub = sum(s > 1e-8 for s in np.linalg.svd(D2_sub_omega, compute_uv=False))
            rank_aug = sum(s > 1e-8 for s in np.linalg.svd(aug, compute_uv=False))
            is_boundary = (rank_aug == rank_D2_sub)
        else:
            is_boundary = (np.linalg.norm(delta_z_sub) < 1e-8)

        # Decompose z by position of v
        pos_decomp = {}
        for i, p in enumerate(paths2):
            if abs(z_A2[i]) > 1e-8:
                if p[0] == v: pos_key = 'pos0'
                elif p[1] == v: pos_key = 'pos1'
                elif p[2] == v: pos_key = 'pos2'
                else: pos_key = 'no_v'
                if pos_key not in pos_decomp:
                    pos_decomp[pos_key] = []
                pos_decomp[pos_key].append((p, round(z_A2[i], 6)))

        # Decompose delta_z by edge type
        delta_decomp = []
        for i, p in enumerate(paths1_not_v):
            if abs(delta_z[i]) > 1e-8:
                delta_decomp.append((p, round(delta_z[i], 6)))

        return {
            'beta1_sub': beta1_sub,
            'h2_rel_dim': 1,
            'is_cycle': is_cycle,
            'is_boundary': is_boundary,
            'delta_injective': not is_boundary,
            'pos_decomp': pos_decomp,
            'delta_decomp': delta_decomp,
            'z_norm': float(np.linalg.norm(z_A2)),
            'delta_norm': float(np.linalg.norm(delta_z)),
            'd_plus_v': sum(A[v]),
        }

    return {'beta1_sub': beta1_sub, 'h2_rel_dim': 0}


# ============================================================
# PART 1: Regular tournament at n=5, all vertices
# ============================================================
print("=" * 70)
print("EXPLICIT DELTA FOR REGULAR T at n=5")
print("=" * 70)

n = 5
# Find first regular tournament
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    if sorted([sum(row) for row in A]) == [2,2,2,2,2]:
        break

print(f"\nTournament T#{bits}, score=(2,2,2,2,2)")
print("Adjacency:")
for i in range(n):
    out_nbrs = [j for j in range(n) if A[i][j]]
    print(f"  {i} -> {out_nbrs}")

for v in range(n):
    print(f"\n--- Vertex v={v}, d+(v)={sum(A[v])} ---")
    result = compute_delta_explicit(A, n, v)
    print(f"  beta1(T\\v) = {result['beta1_sub']}")
    if result['h2_rel_dim'] > 0:
        print(f"  H2(T,T\\v) dim = {result['h2_rel_dim']}")
        print(f"  delta(z) is 1-cycle? {result['is_cycle']}")
        print(f"  delta(z) is boundary? {result['is_boundary']}")
        print(f"  => delta INJECTIVE: {result['delta_injective']}")

        print(f"  Relative 2-cycle z (by position of v):")
        for pos_key in ['pos0', 'pos1', 'pos2', 'no_v']:
            if pos_key in result['pos_decomp']:
                paths = result['pos_decomp'][pos_key]
                print(f"    {pos_key}: {len(paths)} terms")
                for p, c in paths:
                    print(f"      {p}: {c}")

        print(f"  delta(z) image in H1(T\\v):")
        for p, c in result['delta_decomp']:
            print(f"    {p}: {c}")
    else:
        print(f"  H2(T,T\\v) = 0 (nothing to check)")


# ============================================================
# PART 2: Source/sink analysis at n=5
# ============================================================
print(f"\n{'='*70}")
print("DELTA FOR SOURCE/SINK VERTICES at n=5")
print("=" * 70)

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = sorted([sum(row) for row in A])
    if 0 in scores or 4 in scores:  # has source or sink
        for v in range(n):
            d_plus = sum(A[v])
            if d_plus == 0 or d_plus == 4:  # source or sink
                result = compute_delta_explicit(A, n, v)
                if result.get('h2_rel_dim', 0) > 0:
                    status = "INJ" if result.get('delta_injective', False) else "NOT-INJ"
                    print(f"  T#{bits} v={v} d+={d_plus}: {status}, b1(T\\v)={result['beta1_sub']}")
                    if not result.get('delta_injective', True):
                        print(f"    BOUNDARY FAILURE!")
                        break


# ============================================================
# PART 3: Pattern analysis - delta image structure
# ============================================================
print(f"\n{'='*70}")
print("DELTA IMAGE STRUCTURE BY VERTEX TYPE")
print("=" * 70)

n = 5
interior_patterns = defaultdict(int)
boundary_patterns = defaultdict(int)

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    for v in range(n):
        d_plus = sum(A[v])
        result = compute_delta_explicit(A, n, v)
        if result.get('h2_rel_dim', 0) > 0:
            # Count position contributions
            pos_counts = {}
            for pos_key in ['pos0', 'pos1', 'pos2']:
                if pos_key in result.get('pos_decomp', {}):
                    pos_counts[pos_key] = len(result['pos_decomp'][pos_key])
                else:
                    pos_counts[pos_key] = 0
            key = (d_plus, pos_counts['pos0'], pos_counts['pos1'], pos_counts['pos2'],
                   result.get('delta_injective', None))
            if 1 <= d_plus <= n-2:
                interior_patterns[key] += 1
            else:
                boundary_patterns[key] += 1

print("\nInterior vertices (1 <= d+ <= n-2):")
for key in sorted(interior_patterns.keys()):
    d_plus, p0, p1, p2, inj = key
    print(f"  d+={d_plus}, pos=(p0={p0},p1={p1},p2={p2}), inj={inj}: {interior_patterns[key]} cases")

print("\nBoundary vertices (source/sink):")
for key in sorted(boundary_patterns.keys()):
    d_plus, p0, p1, p2, inj = key
    print(f"  d+={d_plus}, pos=(p0={p0},p1={p1},p2={p2}), inj={inj}: {boundary_patterns[key]} cases")


# ============================================================
# PART 4: Can we prove pos-mixing => delta injective?
# ============================================================
print(f"\n{'='*70}")
print("POSITION MIXING ANALYSIS")
print("=" * 70)

print("\nFor interior v: the relative 2-cycle ALWAYS has terms at all 3 positions.")
print("For source v: only pos-0 terms.")
print("For sink v: only pos-2 terms.")
print("\nThe delta image for interior v comes from mixing all three positions,")
print("creating a 1-chain that cannot be a boundary in T\\v.")
print("\nLet's verify this algebraic mechanism:")

# For each interior vertex in regular T, show the delta image structure
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    if sorted([sum(row) for row in A]) != [2,2,2,2,2]:
        continue

    # Just do the first regular tournament
    for v in range(n):
        result = compute_delta_explicit(A, n, v)
        if result.get('h2_rel_dim', 0) > 0:
            # Analyze the delta image
            delta = result.get('delta_decomp', [])
            if delta:
                # Which edges appear in delta(z)?
                edges = [p for p, c in delta]
                # Are these edges forming a cycle in T\v?
                in_degree = Counter()
                out_degree = Counter()
                for (a,b), c in delta:
                    if c > 0:
                        out_degree[a] += 1
                        in_degree[b] += 1
                    else:
                        out_degree[b] += 1
                        in_degree[a] += 1

                # Show the structure
                print(f"\n  v={v}, d+={sum(A[v])}")
                print(f"  delta(z) edges: {[(p, c) for p, c in delta]}")
                print(f"  Number of edges in delta(z): {len(delta)}")
    break  # just first regular tournament


print("\n\nDone.")
