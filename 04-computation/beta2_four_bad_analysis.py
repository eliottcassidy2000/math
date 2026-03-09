#!/usr/bin/env python3
"""
beta2_four_bad_analysis.py — Investigate why ≤3 bad vertices when β₁(T)=0

A "bad vertex" v has β₁(T\v) = 1. Computationally Σ_v β₁(T\v) ≤ 3 always.
KEY QUESTION: Why can't there be 4 bad vertices?

Analysis:
1. Find tournaments with β₁(T)=0 and exactly 3 bad vertices; examine structure
2. Try to construct tournaments with 4 bad vertices via arc flips
3. Analyze the mechanism preventing 4 bad vertices

Uses path_homology_v2.py functions (imported without running validation).
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import random
import sys
import time

# ============================================================
# Core functions from path_homology_v2.py (copied to avoid validation output)
# ============================================================

def enumerate_allowed_paths(A, n, p):
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))
    if p == 0:
        return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None:
        max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti.append(0)
            continue
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0
        beta_p = ker_dim - im_dim
        betti.append(max(0, beta_p))
    return betti

# ============================================================
# Helper functions
# ============================================================

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def delete_vertex(A, n, v):
    """Return adjacency matrix of T \ v."""
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    ri, ci = 0, 0
    for i in range(n):
        if i == v:
            continue
        ci = 0
        for j in range(n):
            if j == v:
                continue
            B[ri][ci] = A[i][j]
            ci += 1
        ri += 1
    return B, new_n

def vertex_map_after_delete(n, v):
    """Map from old vertex indices to new indices after deleting v."""
    m = {}
    new_idx = 0
    for i in range(n):
        if i != v:
            m[i] = new_idx
            new_idx += 1
    return m

def score_sequence(A, n, vertices=None):
    """Out-degree sequence of vertices (sorted)."""
    if vertices is None:
        vertices = range(n)
    scores = []
    for v in vertices:
        s = sum(A[v][j] for j in range(n))
        scores.append(s)
    scores.sort()
    return tuple(scores)

def subtournament(A, n, vertices):
    """Extract subtournament on given vertices."""
    vlist = sorted(vertices)
    k = len(vlist)
    B = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, k

def is_transitive(A, n):
    """Check if tournament is transitive (acyclic)."""
    # Transitive iff score sequence is (0,1,...,n-1)
    scores = sorted(sum(A[i][j] for j in range(n)) for i in range(n))
    return scores == list(range(n))

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def find_h1_generator(A, n):
    """Find a generator of H_1(T) when β_1=1.
    Returns a list of (coeff, edge) pairs representing the 1-cycle in Ω_1.
    """
    allowed_0 = enumerate_allowed_paths(A, n, 0)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)

    omega_1 = compute_omega_basis(A, n, 1, allowed_1, allowed_0)
    dim_omega_1 = omega_1.shape[1] if omega_1.ndim == 2 else 0
    if dim_omega_1 == 0:
        return None

    # ker(∂_1) on Ω_1
    bd_1 = build_full_boundary_matrix(allowed_1, allowed_0)
    bd_1_omega = bd_1 @ omega_1
    U, S, Vt = np.linalg.svd(bd_1_omega, full_matrices=True)
    rank_1 = sum(s > 1e-8 for s in S)
    ker_1_basis = Vt[rank_1:].T  # in Ω_1 coordinates

    # im(∂_2) in Ω_1
    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0
    if dim_omega_2 > 0:
        bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
        bd_2_omega = bd_2 @ omega_2
        # Express image in A_1 coordinates, project to ker(∂_1)
        # For simplicity, express ker_1 in A_1 coordinates
        ker_1_A1 = omega_1 @ ker_1_basis  # columns in A_1 coordinates

        # im(∂_2) in A_1 coordinates
        im_A1 = bd_2_omega  # columns in A_1

        # Find quotient: H_1 = ker / im
        # Stack im columns and find complement in ker
        if im_A1.shape[1] > 0:
            # Project ker basis vectors to quotient
            S_im = np.linalg.svd(im_A1, compute_uv=False)
            im_rank = sum(s > 1e-8 for s in S_im)
        else:
            im_rank = 0

        h1_dim = ker_1_A1.shape[1] - im_rank
        if h1_dim <= 0:
            return None
    else:
        ker_1_A1 = omega_1 @ ker_1_basis
        if ker_1_A1.shape[1] == 0:
            return None

    # Take first ker vector as representative
    z = ker_1_A1[:, 0]

    # Express as edges
    edges = []
    for i, path in enumerate(allowed_1):
        if abs(z[i]) > 1e-8:
            edges.append((round(z[i], 6), path))
    return edges

def flip_arc(A, n, i, j):
    """Return new tournament with arc (i,j) flipped."""
    B = [row[:] for row in A]
    B[i][j] = 1 - B[i][j]
    B[j][i] = 1 - B[j][i]
    return B

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def tournament_key(A, n):
    """Hashable key for a tournament."""
    return tuple(tuple(row) for row in A)

# ============================================================
# ANALYSIS 1: Find and examine 3-bad-vertex tournaments
# ============================================================

def analyze_bad_vertices(A, n):
    """For tournament A on n vertices, find bad vertices and their β₁(T\v)."""
    bad = []
    good = []
    for v in range(n):
        B, m = delete_vertex(A, n, v)
        betti = path_betti_numbers(B, m, max_dim=1)
        b1 = betti[1] if len(betti) > 1 else 0
        if b1 > 0:
            bad.append((v, b1))
        else:
            good.append(v)
    return bad, good

def run_analysis_1(n_val, tournaments_iter, label="", max_count=None):
    """Analysis 1: examine 3-bad-vertex tournaments."""
    print(f"\n{'='*70}")
    print(f"ANALYSIS 1: Tournaments with β₁(T)=0 and exactly 3 bad vertices (n={n_val}, {label})")
    print(f"{'='*70}")

    sigma_dist = Counter()
    three_bad_examples = []
    total_b1_zero = 0
    total_checked = 0
    max_bad = 0

    for A in tournaments_iter:
        total_checked += 1
        if max_count and total_checked > max_count:
            break

        betti = path_betti_numbers(A, n_val, max_dim=1)
        if betti[1] != 0:
            continue
        total_b1_zero += 1

        bad, good = analyze_bad_vertices(A, n_val)
        sigma = len(bad)
        sigma_dist[sigma] += 1
        if sigma > max_bad:
            max_bad = sigma

        if sigma == 3:
            three_bad_examples.append((A, bad, good))

    print(f"\nTotal checked: {total_checked}")
    print(f"Total with β₁(T)=0: {total_b1_zero}")
    print(f"Max bad vertices seen: {max_bad}")
    print(f"\nΣ distribution (among β₁(T)=0 tournaments):")
    for s in sorted(sigma_dist.keys()):
        print(f"  Σ = {s}: {sigma_dist[s]} tournaments ({100*sigma_dist[s]/total_b1_zero:.1f}%)")

    print(f"\n--- Detailed analysis of 3-bad-vertex cases ---")
    print(f"Found {len(three_bad_examples)} tournaments with exactly 3 bad vertices")

    transitive_count = 0
    sub_score_dist = Counter()

    for idx, (A, bad, good) in enumerate(three_bad_examples):
        bad_verts = [v for v, b1 in bad]

        # Subtournament on bad vertices
        sub_A, sub_n = subtournament(A, n_val, bad_verts)
        trans = is_transitive(sub_A, sub_n)
        if trans:
            transitive_count += 1

        sub_scores = score_sequence(sub_A, sub_n)
        sub_score_dist[sub_scores] += 1

        # Print first few examples in detail
        if idx < 5:
            print(f"\n  Example {idx+1}: bad vertices = {bad_verts}")
            print(f"    Full scores: {score_sequence(A, n_val)}")
            print(f"    Bad vertex scores in T: {[sum(A[v][j] for j in range(n_val)) for v in bad_verts]}")
            print(f"    Subtournament on bad verts transitive? {trans}")
            print(f"    Sub score seq: {sub_scores}")

            # Find H_1 generators for each T\v
            for v, b1 in bad:
                B, m = delete_vertex(A, n_val, v)
                gen = find_h1_generator(B, m)
                vmap = vertex_map_after_delete(n_val, v)
                inv_map = {new: old for old, new in vmap.items()}
                if gen:
                    # Translate back to original vertex labels
                    orig_edges = [(c, (inv_map[p[0]], inv_map[p[1]])) for c, p in gen]
                    print(f"    z_{v} (H₁ gen of T\\{v}): {orig_edges}")

            # Check edge overlap between generators
            all_gen_edges = []
            for v, b1 in bad:
                B, m = delete_vertex(A, n_val, v)
                gen = find_h1_generator(B, m)
                if gen:
                    vmap = vertex_map_after_delete(n_val, v)
                    inv_map = {new: old for old, new in vmap.items()}
                    edges_set = set()
                    for c, p in gen:
                        if abs(c) > 1e-8:
                            edges_set.add((inv_map[p[0]], inv_map[p[1]]))
                    all_gen_edges.append((v, edges_set))

            if len(all_gen_edges) >= 2:
                print(f"    Edge overlaps between generators:")
                for i in range(len(all_gen_edges)):
                    for j in range(i+1, len(all_gen_edges)):
                        v1, e1 = all_gen_edges[i]
                        v2, e2 = all_gen_edges[j]
                        overlap = e1 & e2
                        print(f"      z_{v1} ∩ z_{v2}: {len(overlap)} edges: {overlap}")

    print(f"\n  Summary:")
    print(f"    Transitive subtournament on bad verts: {transitive_count}/{len(three_bad_examples)} ({100*transitive_count/max(1,len(three_bad_examples)):.1f}%)")
    print(f"    Sub score distributions: {dict(sub_score_dist)}")

    return three_bad_examples, sigma_dist

# ============================================================
# ANALYSIS 2: Try to construct 4-bad-vertex tournaments
# ============================================================

def run_analysis_2(n_val, three_bad_examples, label=""):
    """Try arc flips to create a 4th bad vertex."""
    print(f"\n{'='*70}")
    print(f"ANALYSIS 2: Can arc flips create a 4th bad vertex? (n={n_val})")
    print(f"{'='*70}")

    attempts = 0
    successes = 0
    failure_reasons = Counter()

    for ex_idx, (A, bad, good) in enumerate(three_bad_examples[:20]):  # limit for speed
        bad_verts = set(v for v, _ in bad)

        # Try flipping each arc
        for i in range(n_val):
            for j in range(i+1, n_val):
                B = flip_arc(A, n_val, i, j)
                attempts += 1

                # Check β₁(B) = 0
                betti_B = path_betti_numbers(B, n_val, max_dim=1)
                if betti_B[1] != 0:
                    failure_reasons["flip breaks β₁(T)=0"] += 1
                    continue

                # Count bad vertices in B
                bad_B, good_B = analyze_bad_vertices(B, n_val)
                sigma_B = len(bad_B)

                if sigma_B >= 4:
                    successes += 1
                    print(f"\n  *** FOUND 4+ bad vertices! Example {ex_idx}, flip ({i},{j}) ***")
                    print(f"    Bad vertices: {[v for v, _ in bad_B]}")
                    print(f"    Original bad: {[v for v, _ in bad]}")
                elif sigma_B == 3:
                    new_bad = set(v for v, _ in bad_B)
                    if new_bad != bad_verts:
                        failure_reasons["3 bad but different set"] += 1
                    else:
                        failure_reasons["same 3 bad vertices"] += 1
                elif sigma_B < 3:
                    failure_reasons[f"only {sigma_B} bad"] += 1

    print(f"\n  Total arc-flip attempts: {attempts}")
    print(f"  Successes (4+ bad): {successes}")
    print(f"  Failure breakdown:")
    for reason, count in failure_reasons.most_common():
        print(f"    {reason}: {count} ({100*count/max(1,attempts):.1f}%)")

# ============================================================
# ANALYSIS 3: Mechanism analysis
# ============================================================

def run_analysis_3(n_val, three_bad_examples, label=""):
    """Analyze the mechanism preventing 4 bad vertices."""
    print(f"\n{'='*70}")
    print(f"ANALYSIS 3: Mechanism analysis (n={n_val})")
    print(f"{'='*70}")

    for ex_idx, (A, bad, good) in enumerate(three_bad_examples[:10]):
        bad_verts = [v for v, _ in bad]

        if ex_idx < 3:
            print(f"\n  --- Example {ex_idx+1} ---")
            print(f"  Bad vertices: {bad_verts}")
            print(f"  Good vertices: {good}")

        # Check 1: β₁ of T minus 2 bad vertices
        if ex_idx < 3:
            print(f"\n  β₁ after deleting pairs of bad vertices:")
        for pair in combinations(bad_verts, 2):
            remaining = [v for v in range(n_val) if v not in pair]
            sub_A, sub_n = subtournament(A, n_val, remaining)
            betti_sub = path_betti_numbers(sub_A, sub_n, max_dim=1)
            if ex_idx < 3:
                print(f"    T \\ {{{pair[0]},{pair[1]}}}: β = {betti_sub}")

        # Check 2: Common neighbor patterns
        if ex_idx < 3:
            print(f"\n  Neighborhood patterns of bad vertices:")
            for v in bad_verts:
                out_nbrs = [j for j in range(n_val) if j != v and A[v][j] == 1]
                in_nbrs = [j for j in range(n_val) if j != v and A[j][v] == 1]
                print(f"    v={v}: out={out_nbrs}, in={in_nbrs}")

            # Common out-neighbors among bad vertices
            for pair in combinations(bad_verts, 2):
                v1, v2 = pair
                common_out = [j for j in range(n_val) if j not in pair and A[v1][j] == 1 and A[v2][j] == 1]
                common_in = [j for j in range(n_val) if j not in pair and A[j][v1] == 1 and A[j][v2] == 1]
                print(f"    Common out-nbrs of {v1},{v2}: {common_out}")
                print(f"    Common in-nbrs of {v1},{v2}: {common_in}")

        # Check 3: Is there a filling witness in Ω₂?
        if ex_idx < 3:
            print(f"\n  Ω₂ filling analysis:")
            allowed_1 = enumerate_allowed_paths(A, n_val, 1)
            allowed_2 = enumerate_allowed_paths(A, n_val, 2)
            omega_2 = compute_omega_basis(A, n_val, 2, allowed_2, allowed_1)
            dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0
            print(f"    dim(Ω₂) = {dim_omega_2}")
            print(f"    |A₂| = {len(allowed_2)} (transitive triples)")

            # For each bad vertex, find the z_v cycle
            generators = {}
            for v, b1 in bad:
                B, m = delete_vertex(A, n_val, v)
                gen = find_h1_generator(B, m)
                if gen:
                    vmap = vertex_map_after_delete(n_val, v)
                    inv_map = {new: old for old, new in vmap.items()}
                    # Build vector in A_1(T) coordinates (lift the cycle)
                    edge_idx = {path: i for i, path in enumerate(allowed_1)}
                    z_vec = np.zeros(len(allowed_1))
                    for c, p in gen:
                        orig_edge = (inv_map[p[0]], inv_map[p[1]])
                        if orig_edge in edge_idx:
                            z_vec[edge_idx[orig_edge]] = c
                    generators[v] = z_vec

            # Check if each z_v is in im(∂₂) in T
            if dim_omega_2 > 0:
                bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
                bd_2_omega = bd_2 @ omega_2  # image of ∂₂ in A₁

                for v in generators:
                    z = generators[v]
                    # Check if z is in column space of bd_2_omega
                    # Augment and check rank
                    aug = np.column_stack([bd_2_omega, z.reshape(-1,1)])
                    rank_orig = np.linalg.matrix_rank(bd_2_omega, tol=1e-8)
                    rank_aug = np.linalg.matrix_rank(aug, tol=1e-8)
                    in_image = (rank_aug == rank_orig)
                    print(f"    z_{v} ∈ im(∂₂) of T? {in_image}")

                    if in_image and dim_omega_2 > 0:
                        # Find the filling: solve bd_2_omega @ x = z
                        x, residuals, rank, sv = np.linalg.lstsq(bd_2_omega, z, rcond=None)
                        # Express in Ω₂ coordinates
                        filling = omega_2 @ x  # in A₂ coordinates
                        nonzero_terms = [(round(filling[i], 4), allowed_2[i]) for i in range(len(allowed_2)) if abs(filling[i]) > 1e-6]
                        print(f"      Filling of z_{v}: {nonzero_terms[:10]}")

                # Check if there's a SINGLE Ω₂ element that fills all three
                if len(generators) == 3:
                    z_sum = sum(generators.values())
                    aug_sum = np.column_stack([bd_2_omega, z_sum.reshape(-1,1)])
                    rank_sum = np.linalg.matrix_rank(aug_sum, tol=1e-8)
                    sum_in_image = (rank_sum == rank_orig)
                    print(f"\n    z_sum (sum of all 3 generators) ∈ im(∂₂)? {sum_in_image}")

    # Global statistics for all 3-bad examples
    print(f"\n  --- Global statistics across all 3-bad examples ---")

    # Check subtournament transitivity for ALL examples
    trans_count = 0
    sub_3cycle_count = 0
    for A, bad, good in three_bad_examples:
        bad_verts = [v for v, _ in bad]
        sub_A, sub_n = subtournament(A, n_val, bad_verts)
        if is_transitive(sub_A, sub_n):
            trans_count += 1
        else:
            sub_3cycle_count += 1

    print(f"    Transitive sub on bad verts: {trans_count}/{len(three_bad_examples)}")
    print(f"    3-cycle sub on bad verts: {sub_3cycle_count}/{len(three_bad_examples)}")

    # Check: for each good vertex, what happens if we look at T\v?
    print(f"\n  --- Good vertex β₁ analysis ---")
    good_beta1_dist = Counter()
    for A, bad, good_verts in three_bad_examples[:20]:
        for v in good_verts:
            B, m = delete_vertex(A, n_val, v)
            betti = path_betti_numbers(B, m, max_dim=1)
            good_beta1_dist[betti[1]] += 1
    print(f"    β₁(T\\v) for good v: {dict(good_beta1_dist)}")

# ============================================================
# ANALYSIS 4: Exhaustive search for 4-bad-vertex tournaments
# ============================================================

def exhaustive_search_4bad(n_val, tournaments_iter, max_count=None):
    """Exhaustively search for any tournament with 4+ bad vertices and β₁=0."""
    print(f"\n{'='*70}")
    print(f"EXHAUSTIVE SEARCH for 4+ bad vertices with β₁(T)=0 (n={n_val})")
    print(f"{'='*70}")

    total = 0
    b1_zero = 0
    max_bad = 0
    max_bad_example = None
    sigma_dist = Counter()

    for A in tournaments_iter:
        total += 1
        if max_count and total > max_count:
            break
        if total % 1000 == 0:
            print(f"  ... checked {total} tournaments, {b1_zero} with β₁=0, max_bad={max_bad}", flush=True)

        betti = path_betti_numbers(A, n_val, max_dim=1)
        if betti[1] != 0:
            continue
        b1_zero += 1

        bad, good = analyze_bad_vertices(A, n_val)
        sigma = len(bad)
        sigma_dist[sigma] += 1

        if sigma > max_bad:
            max_bad = sigma
            max_bad_example = (A, bad, good)
            if sigma >= 4:
                print(f"\n  *** FOUND {sigma} bad vertices! ***")
                print(f"    A = {A}")
                print(f"    Bad: {bad}")

    print(f"\n  Total checked: {total}")
    print(f"  β₁=0 count: {b1_zero}")
    print(f"  Maximum bad vertices: {max_bad}")
    print(f"  Σ distribution:")
    for s in sorted(sigma_dist.keys()):
        print(f"    Σ={s}: {sigma_dist[s]}")

    return max_bad, sigma_dist

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("="*70)
    print("β₂=0 PROOF INVESTIGATION: Why can't there be 4 bad vertices?")
    print("="*70)
    print(f"Date: 2026-03-08")
    print()

    t0 = time.time()

    # ---- n=5 (exhaustive: 2^10 = 1024 tournaments) ----
    print("\n" + "#"*70)
    print(f"# n = 5 (exhaustive)")
    print("#"*70)

    three_bad_5, sigma_5 = run_analysis_1(5, all_tournaments(5), "exhaustive")
    if three_bad_5:
        run_analysis_2(5, three_bad_5, "exhaustive")
        run_analysis_3(5, three_bad_5, "exhaustive")
    exhaustive_search_4bad(5, all_tournaments(5))

    print(f"\nn=5 done. Elapsed: {time.time()-t0:.1f}s")

    # ---- n=6 (exhaustive: 2^15 = 32768 tournaments) ----
    t1 = time.time()
    print("\n" + "#"*70)
    print(f"# n = 6 (exhaustive)")
    print("#"*70)

    three_bad_6, sigma_6 = run_analysis_1(6, all_tournaments(6), "exhaustive")
    if three_bad_6:
        run_analysis_2(6, three_bad_6, "exhaustive")
        run_analysis_3(6, three_bad_6, "exhaustive")
    exhaustive_search_4bad(6, all_tournaments(6))

    print(f"\nn=6 done. Elapsed: {time.time()-t1:.1f}s")

    # ---- n=7 (sampled) ----
    t2 = time.time()
    print("\n" + "#"*70)
    print(f"# n = 7 (sampled: 5000 random tournaments)")
    print("#"*70)

    random.seed(42)
    sample_7 = [random_tournament(7) for _ in range(5000)]

    three_bad_7, sigma_7 = run_analysis_1(7, iter(sample_7), "sampled 5000")
    if three_bad_7:
        run_analysis_2(7, three_bad_7[:10], "sampled")
        run_analysis_3(7, three_bad_7, "sampled")

    # Also check max Σ in sample
    print(f"\n  Checking max Σ in n=7 sample...")
    max_bad_7 = 0
    sigma_7_full = Counter()
    for A in sample_7:
        betti = path_betti_numbers(A, 7, max_dim=1)
        if betti[1] != 0:
            continue
        bad, good = analyze_bad_vertices(A, 7)
        sigma = len(bad)
        sigma_7_full[sigma] += 1
        if sigma > max_bad_7:
            max_bad_7 = sigma
    print(f"  Max bad vertices in n=7 sample: {max_bad_7}")
    print(f"  Full Σ distribution: {dict(sorted(sigma_7_full.items()))}")

    print(f"\nn=7 done. Elapsed: {time.time()-t2:.1f}s")

    # ---- SUMMARY ----
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("""
Key findings to report:
1. Is the subtournament on 3 bad vertices ALWAYS transitive?
2. How do the H₁ generators relate to each other?
3. What prevents a 4th bad vertex from appearing?
4. What is the mechanism?
""")

    print(f"\nTotal elapsed: {time.time()-t0:.1f}s")
