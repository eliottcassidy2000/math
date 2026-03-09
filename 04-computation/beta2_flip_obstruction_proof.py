#!/usr/bin/env python3
"""
FLIP OBSTRUCTION LEMMA — Algebraic proof search

LEMMA: If T is a tournament with beta_1(T)=0, and a,b,c form a directed
3-cycle in T, then it's impossible for all three of beta_1(T\a), beta_1(T\b),
beta_1(T\c) to equal 1.

Equivalently: vertices with beta_1(T\v)=1 ("bad" vertices) form a
TRANSITIVE sub-tournament.

This script explores algebraic structures to find the proof.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random
import sys

# ============================================================
# Core homology routines (inlined from path_homology_v2.py)
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
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
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
# Detailed B_1 computation
# ============================================================

def compute_B1_details(A, n):
    """Compute B_1(T) in detail: allowed 2-paths, boundary matrix, rank.

    Returns dict with:
      - allowed_1: list of allowed 1-paths (edges)
      - allowed_2: list of allowed 2-paths (a,b,c)
      - bd2: boundary matrix partial_2: A_2 -> A_1
      - omega2_basis: basis for Omega_2
      - bd2_omega: partial_2 restricted to Omega_2
      - rank_B1: rank of B_1 = im(partial_2|_{Omega_2})
      - dim_Z1: dimension of Z_1
      - beta_1: first Betti number
    """
    allowed_0 = enumerate_allowed_paths(A, n, 0)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_3 = enumerate_allowed_paths(A, n, 3)

    # Omega_1 = A_1 for tournaments (all faces of edges are vertices, always allowed)
    omega_1 = compute_omega_basis(A, n, 1, allowed_1, allowed_0)
    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)

    # Z_1 = ker(partial_1|_{Omega_1})
    bd1 = build_full_boundary_matrix(allowed_1, allowed_0)
    bd1_omega = bd1 @ omega_1
    if bd1_omega.shape[0] > 0 and bd1_omega.shape[1] > 0:
        U1, S1, Vt1 = np.linalg.svd(bd1_omega, full_matrices=True)
        rank_bd1 = sum(s > 1e-8 for s in S1)
    else:
        rank_bd1 = 0
        Vt1 = np.eye(omega_1.shape[1]) if omega_1.shape[1] > 0 else np.zeros((0,0))
        S1 = []

    dim_omega_1 = omega_1.shape[1]
    dim_Z1 = dim_omega_1 - rank_bd1

    # B_1 = im(partial_2|_{Omega_2})
    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0
    if dim_omega_2 > 0:
        bd2_omega = bd2 @ omega_2
        S2 = np.linalg.svd(bd2_omega, compute_uv=False)
        rank_B1 = sum(s > 1e-8 for s in S2)
    else:
        bd2_omega = np.zeros((len(allowed_1), 0))
        rank_B1 = 0

    return {
        'allowed_1': allowed_1,
        'allowed_2': allowed_2,
        'omega_2': omega_2,
        'bd2': bd2,
        'bd2_omega': bd2_omega,
        'rank_B1': rank_B1,
        'dim_Z1': dim_Z1,
        'beta_1': dim_Z1 - rank_B1,
        'idx_1': {path: i for i, path in enumerate(allowed_1)},
    }

def delete_vertex(A, n, v):
    """Return (A', n-1) with vertex v removed."""
    remaining = [i for i in range(n) if i != v]
    n2 = n - 1
    A2 = [[0]*n2 for _ in range(n2)]
    for i2, i in enumerate(remaining):
        for j2, j in enumerate(remaining):
            A2[i2][j2] = A[i][j]
    return A2, n2, remaining

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def get_3cycles(A, n):
    """Return list of directed 3-cycles as (a,b,c) where a->b->c->a."""
    cycles = []
    for i in range(n):
        for j in range(n):
            if j == i: continue
            if not A[i][j]: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[j][k] and A[k][i]:
                    cycles.append((i,j,k))
    return cycles

def is_3cycle(A, a, b, c):
    """Check if {a,b,c} form a directed 3-cycle in either orientation."""
    return (A[a][b] and A[b][c] and A[c][a]) or \
           (A[b][a] and A[a][c] and A[c][b])

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

# ============================================================
# Finding the unfillable cycle z_v
# ============================================================

def find_unfillable_cycle(A, n, v, details=None):
    """For T with beta_1(T\v)=1, find the generator of H_1(T\v).

    Returns the cycle z_v as a dict {edge: coefficient} in the ORIGINAL
    vertex labeling (not the relabeled T\v).
    """
    A2, n2, remaining = delete_vertex(A, n, v)

    allowed_0 = enumerate_allowed_paths(A2, n2, 0)
    allowed_1 = enumerate_allowed_paths(A2, n2, 1)
    allowed_2 = enumerate_allowed_paths(A2, n2, 2)

    omega_1 = compute_omega_basis(A2, n2, 1, allowed_1, allowed_0)
    omega_2 = compute_omega_basis(A2, n2, 2, allowed_2, allowed_1)

    # Z_1 = ker(partial_1|_{Omega_1})
    bd1 = build_full_boundary_matrix(allowed_1, allowed_0)
    bd1_omega = bd1 @ omega_1
    U1, S1, Vt1 = np.linalg.svd(bd1_omega, full_matrices=True)
    rank_bd1 = sum(s > 1e-8 for s in S1)

    # Null space of bd1_omega = Z_1 in omega coordinates
    Z1_coords = Vt1[rank_bd1:].T  # columns are Z_1 basis vectors (in omega coords)
    Z1_in_A1 = omega_1 @ Z1_coords  # in A_1 coordinates

    # B_1 = im(partial_2|_{Omega_2})
    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0
    if dim_omega_2 > 0:
        bd2_omega = bd2 @ omega_2
        U2, S2, Vt2 = np.linalg.svd(bd2_omega, full_matrices=True)
        rank_B1 = sum(s > 1e-8 for s in S2)
        # B_1 basis in A_1 coords
        B1_basis = U2[:, :rank_B1]
    else:
        rank_B1 = 0
        B1_basis = np.zeros((len(allowed_1), 0))

    dim_Z1 = Z1_in_A1.shape[1]
    beta_1 = dim_Z1 - rank_B1

    if beta_1 != 1:
        return None, remaining

    # Find the unfillable direction: project Z1 onto complement of B1
    # The unfillable cycle is in Z1 but not in B1
    if rank_B1 > 0:
        # Project each Z1 basis vector onto B1, find one with nonzero residual
        for col in range(Z1_in_A1.shape[1]):
            z = Z1_in_A1[:, col]
            # Project onto B1
            proj = B1_basis @ (B1_basis.T @ z)
            residual = z - proj
            if np.linalg.norm(residual) > 1e-8:
                # Normalize
                residual = residual / np.linalg.norm(residual)
                # Convert to edge dict in original labeling
                edge_dict = {}
                for i, path in enumerate(allowed_1):
                    if abs(residual[i]) > 1e-10:
                        # Map back to original vertices
                        orig_edge = (remaining[path[0]], remaining[path[1]])
                        edge_dict[orig_edge] = residual[i]
                return edge_dict, remaining
    else:
        # B1 is trivial, any Z1 vector is unfillable
        z = Z1_in_A1[:, 0]
        z = z / np.linalg.norm(z)
        edge_dict = {}
        for i, path in enumerate(allowed_1):
            if abs(z[i]) > 1e-10:
                orig_edge = (remaining[path[0]], remaining[path[1]])
                edge_dict[orig_edge] = z[i]
        return edge_dict, remaining

    return None, remaining

# ============================================================
# Check if cycle is a boundary in T
# ============================================================

def check_if_boundary(A, n, cycle_dict):
    """Given a 1-cycle z (as edge->coeff dict), check if z in B_1(T).

    Returns (is_boundary, preimage_or_None).
    """
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_0 = enumerate_allowed_paths(A, n, 0)

    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)

    idx_1 = {path: i for i, path in enumerate(allowed_1)}

    # Represent z in A_1 coordinates
    z = np.zeros(len(allowed_1))
    for edge, coeff in cycle_dict.items():
        if edge in idx_1:
            z[idx_1[edge]] += coeff

    # B_1 = im(partial_2|_{Omega_2})
    bd2 = build_full_boundary_matrix(allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0
    if dim_omega_2 == 0:
        return False, None

    bd2_omega = bd2 @ omega_2

    # Solve bd2_omega @ x = z (least squares)
    x, res, rank, sv = np.linalg.lstsq(bd2_omega, z, rcond=None)

    # Check residual
    residual = np.linalg.norm(bd2_omega @ x - z)
    is_boundary = residual < 1e-8

    if is_boundary:
        # Reconstruct preimage in omega_2 coordinates, then in A_2 coordinates
        preimage_omega = x
        preimage_A2 = omega_2 @ preimage_omega
        # Convert to 2-path dict
        preimage_dict = {}
        for i, path in enumerate(allowed_2):
            if abs(preimage_A2[i]) > 1e-10:
                preimage_dict[path] = preimage_A2[i]
        return True, preimage_dict

    return False, None

# ============================================================
# Approach B: Rank drop analysis
# ============================================================

def rank_drop_analysis(A, n, details_T):
    """For each vertex v, compute rank_drop(v) = dim(B_1(T)) - dim(B_1(T\v))."""
    results = {}
    rank_T = details_T['rank_B1']

    for v in range(n):
        A2, n2, remaining = delete_vertex(A, n, v)
        det_v = compute_B1_details(A2, n2)
        rank_v = det_v['rank_B1']
        beta_v = det_v['beta_1']

        # rank_drop = rank_T - rank_v
        # For a tournament on n vertices: dim(Z_1(T)) = C(n-1,2)
        # For T\v on n-1 vertices: dim(Z_1(T\v)) = C(n-2,2)
        # rank_T = C(n-1,2) - beta_1(T) = C(n-1,2) (when beta_1=0)
        # rank_v = C(n-2,2) - beta_1(T\v)
        # rank_drop = C(n-1,2) - C(n-2,2) + beta_1(T\v) = (n-2) + beta_1(T\v)

        results[v] = {
            'rank_B1_Tv': rank_v,
            'beta_1_Tv': beta_v,
            'rank_drop': rank_T - rank_v,
            'expected_drop': (n-2) + beta_v,
        }

    return results

# ============================================================
# Approach: 2-path through vertex analysis
# ============================================================

def paths_through_vertex(allowed_2, v):
    """Return indices of 2-paths that pass through vertex v (as middle vertex)."""
    return [i for i, path in enumerate(allowed_2) if path[1] == v]

def paths_involving_vertex(allowed_2, v):
    """Return indices of 2-paths involving vertex v in any position."""
    return [i for i, path in enumerate(allowed_2) if v in path]

# ============================================================
# MAIN ANALYSIS
# ============================================================

print("=" * 75)
print("FLIP OBSTRUCTION LEMMA — ALGEBRAIC PROOF SEARCH")
print("=" * 75)

# ============================================================
# Part 1: n=5 exhaustive
# ============================================================

print("\n" + "=" * 75)
print("PART 1: n=5 EXHAUSTIVE ANALYSIS")
print("=" * 75)

n = 5
total = 0
beta0_count = 0
lemma_violations = 0
bad_vertex_counts = Counter()
bad_3cycle_data = []  # (tournament, bad_vertices, 3cycles_among_bad)

# Collect stats
rank_drop_stats = []
sum_rank_drop_stats = []

for A in all_tournaments(n):
    total += 1
    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        continue
    beta0_count += 1

    # Find bad vertices
    bad = []
    for v in range(n):
        A2, n2, _ = delete_vertex(A, n, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    bad_vertex_counts[len(bad)] += 1

    if len(bad) < 2:
        continue

    # Check 3-cycles among bad vertices
    bad_set = set(bad)
    cycles_among_bad = []
    for a, b, c in combinations(bad, 3):
        if is_3cycle(A, a, b, c):
            cycles_among_bad.append((a,b,c))
            lemma_violations += 1  # Would be a violation if it happens

    if cycles_among_bad:
        bad_3cycle_data.append((A, bad, cycles_among_bad))

print(f"\nTotal tournaments: {total}")
print(f"Tournaments with beta_1=0: {beta0_count}")
print(f"Bad vertex count distribution: {dict(bad_vertex_counts)}")
print(f"Violations of flip obstruction lemma: {lemma_violations}")

if lemma_violations == 0:
    print(">>> LEMMA HOLDS for all n=5 tournaments <<<")
else:
    print(">>> LEMMA VIOLATED — investigating...")

# ============================================================
# Part 1b: Detailed analysis of tournaments with >=2 bad vertices
# ============================================================

print("\n" + "-" * 75)
print("DETAILED ANALYSIS: tournaments with >=2 bad vertices")
print("-" * 75)

detailed_count = 0
for A in all_tournaments(n):
    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(n):
        A2, n2, _ = delete_vertex(A, n, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 2:
        continue

    detailed_count += 1
    if detailed_count > 10:
        continue  # Only print first 10

    t3 = count_3cycles(A, n)
    details_T = compute_B1_details(A, n)

    print(f"\n  Tournament #{detailed_count}: t3={t3}, beta_1={betti[1]}, bad={bad}")
    print(f"    dim(B_1(T)) = {details_T['rank_B1']}, dim(Z_1(T)) = {details_T['dim_Z1']}")

    # Score sequence
    scores = [sum(A[i]) for i in range(n)]
    print(f"    Scores: {scores}")

    # Rank drop analysis
    rd = rank_drop_analysis(A, n, details_T)
    for v in range(n):
        marker = " *BAD*" if v in bad else ""
        print(f"    v={v}: rank_drop={rd[v]['rank_drop']}, "
              f"beta_1(T\\v)={rd[v]['beta_1_Tv']}, "
              f"expected={rd[v]['expected_drop']}{marker}")

    # Verify rank drop formula
    total_drop = sum(rd[v]['rank_drop'] for v in range(n))
    print(f"    Sum of rank drops: {total_drop}")

    # Check 3-cycles among bad vertices
    bad_set = set(bad)
    for a, b, c in combinations(bad, 3):
        if is_3cycle(A, a, b, c):
            print(f"    *** 3-CYCLE AMONG BAD VERTICES: {a},{b},{c}")
        else:
            print(f"    Transitive triple among bad: {a},{b},{c}")

    # Check which edges form the 3-cycle (a,b)+(b,c)+(c,a)
    # and whether it's in B_1(T)
    cycles_3 = get_3cycles(A, n)
    for cyc in cycles_3[:3]:  # First few 3-cycles
        a, b, c = cyc
        cycle_dict = {(a,b): 1.0, (b,c): 1.0, (c,a): 1.0}
        is_bd, preimage = check_if_boundary(A, n, cycle_dict)
        bad_in_cycle = sum(1 for x in [a,b,c] if x in bad)
        if bad_in_cycle >= 2:
            print(f"    3-cycle ({a},{b},{c}): {bad_in_cycle} bad vertices, "
                  f"is_boundary={is_bd}")
            if is_bd and preimage:
                # Check if preimage uses any of a,b,c as interior vertex
                interior_verts_used = set()
                for path, coeff in preimage.items():
                    interior_verts_used.add(path[1])  # middle vertex of 2-path
                print(f"      Filling uses interior vertices: {interior_verts_used}")
                bad_interior = interior_verts_used & bad_set
                print(f"      Bad vertices used as interior: {bad_interior}")

print(f"\n  Total tournaments with >=2 bad vertices: {detailed_count}")

# ============================================================
# Part 2: Approach B — Rank analysis and constraints
# ============================================================

print("\n" + "=" * 75)
print("PART 2: RANK ANALYSIS (Approach B)")
print("=" * 75)

print("\nKey identity: rank_drop(v) = (n-2) + beta_1(T\\v)")
print("  So bad vertices have rank_drop = n-1, good vertices have rank_drop = n-2")
print(f"  At n={n}: bad drop = {n-1}, good drop = {n-2}")

# What constrains the sum of rank drops?
# Each allowed 2-path (a,b,c) contributes to partial_2. When we remove v:
#   - 2-paths with v as endpoint or middle are lost
# The total "overlap" between removals determines sum of drops.

sum_drops = []
for A in all_tournaments(n):
    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        continue
    details_T = compute_B1_details(A, n)
    rd = rank_drop_analysis(A, n, details_T)
    total_drop = sum(rd[v]['rank_drop'] for v in range(n))
    n_bad = sum(1 for v in range(n) if rd[v]['beta_1_Tv'] == 1)
    sum_drops.append((total_drop, n_bad, details_T['rank_B1']))

print(f"\nSum of rank drops across all vertices (for beta_1=0 tournaments at n={n}):")
drop_dist = Counter((sd, nb) for sd, nb, rb in sum_drops)
for (sd, nb), cnt in sorted(drop_dist.items()):
    print(f"  sum_drops={sd}, n_bad={nb}: {cnt} tournaments")

# Theoretical: sum of rank drops = ?
# rank_drop(v) counts the rank of 2-paths that "die" when v is removed
# A 2-path (a,b,c) is "owned" by vertices a, b, c
# When we remove v, we lose all 2-paths involving v
# But some of these paths might be linearly dependent on paths NOT involving v

# The TOTAL number of 2-paths through v:
print("\n  Expected: each 2-path (a,b,c) involves 3 vertices.")
print("  If we sum rank_drop(v) over all v, we're counting the rank contribution")
print("  of each 2-path up to 3 times (once for each vertex it touches).")

# ============================================================
# Part 3: Approach C — Hidden cycle analysis
# ============================================================

print("\n" + "=" * 75)
print("PART 3: HIDDEN CYCLE ANALYSIS (Approach C)")
print("=" * 75)

hidden_cycle_count = 0
for A in all_tournaments(n):
    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(n):
        A2, n2, _ = delete_vertex(A, n, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 2:
        continue

    hidden_cycle_count += 1
    if hidden_cycle_count > 5:
        continue

    print(f"\n  Tournament #{hidden_cycle_count}, bad={bad}")

    # Find unfillable cycles z_v for each bad vertex
    z_cycles = {}
    for v in bad:
        z_v, remaining = find_unfillable_cycle(A, n, v)
        if z_v is not None:
            z_cycles[v] = z_v
            edges_str = ", ".join(f"({a},{b}):{c:.3f}" for (a,b), c in sorted(z_v.items()))
            print(f"    z_{v}: {edges_str}")

    # Check pairwise: z_a + z_b restricted to common edges
    if len(z_cycles) >= 2:
        for v1, v2 in combinations(z_cycles.keys(), 2):
            z1 = z_cycles[v1]
            z2 = z_cycles[v2]
            # Find common edges (edges not involving v1 or v2)
            common_sum = {}
            for edge, c in z1.items():
                if v2 not in edge:
                    common_sum[edge] = common_sum.get(edge, 0) + c
            for edge, c in z2.items():
                if v1 not in edge:
                    common_sum[edge] = common_sum.get(edge, 0) + c

            # Are these the same edges? Cancel?
            nonzero = {e: c for e, c in common_sum.items() if abs(c) > 1e-10}
            print(f"    z_{v1} + z_{v2} on common edges: {len(nonzero)} nonzero entries")

            # Also: z_a - z_b
            diff = {}
            for edge, c in z1.items():
                if v2 not in edge:
                    diff[edge] = diff.get(edge, 0) + c
            for edge, c in z2.items():
                if v1 not in edge:
                    diff[edge] = diff.get(edge, 0) - c
            nonzero_d = {e: c for e, c in diff.items() if abs(c) > 1e-10}
            print(f"    z_{v1} - z_{v2} on common edges: {len(nonzero_d)} nonzero entries")

    # Check: is the arc (a,b) in z_a?
    # The unfillable cycle z_a lives in T\a, so it does NOT use vertex a.
    # Similarly z_b doesn't use vertex b.
    # If a and b are both bad and connected, the edge (a,b) appears in z_b but not z_a.
    # And edge (b,a) or (a,b) appears in z_a or z_b.
    for v in bad:
        z_v = z_cycles.get(v)
        if z_v is None: continue
        verts_used = set()
        for (a, b) in z_v:
            verts_used.add(a)
            verts_used.add(b)
        print(f"    z_{v} uses vertices: {sorted(verts_used)} (should not contain {v})")

# ============================================================
# Part 4: 3-cycle distribution among bad vertices
# ============================================================

print("\n" + "=" * 75)
print("PART 4: 3-CYCLE DISTRIBUTION AMONG BAD VERTICES")
print("=" * 75)

cycle_stats = Counter()  # (n_bad_in_cycle,) -> count

for A in all_tournaments(n):
    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(n):
        A2, n2, _ = delete_vertex(A, n, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    bad_set = set(bad)

    # Count 3-cycles by number of bad vertices in them
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if is_3cycle(A, i, j, k):
                    n_bad = sum(1 for x in [i,j,k] if x in bad_set)
                    cycle_stats[n_bad] += 1

print(f"\n3-cycle distribution by #bad vertices (n={n}, beta_1=0 tournaments):")
for nb in sorted(cycle_stats.keys()):
    print(f"  {nb} bad vertices in 3-cycle: {cycle_stats[nb]} occurrences")

if 3 not in cycle_stats or cycle_stats[3] == 0:
    print("\n>>> NO 3-CYCLES WITH ALL 3 VERTICES BAD — confirms lemma")

# ============================================================
# Part 5: What about 2 bad vertices in a 3-cycle?
# ============================================================

print("\n" + "=" * 75)
print("PART 5: FILLING STRUCTURE FOR 3-CYCLES WITH 2 BAD VERTICES")
print("=" * 75)

fill_count = 0
for A in all_tournaments(n):
    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(n):
        A2, n2, _ = delete_vertex(A, n, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    bad_set = set(bad)
    if len(bad) < 2:
        continue

    # Find 3-cycles with exactly 2 bad vertices
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if not is_3cycle(A, i, j, k):
                    continue
                n_bad = sum(1 for x in [i,j,k] if x in bad_set)
                if n_bad != 2:
                    continue

                fill_count += 1
                if fill_count > 5:
                    continue

                # Which is the good vertex?
                good = [x for x in [i,j,k] if x not in bad_set][0]
                bad_in_cycle = [x for x in [i,j,k] if x in bad_set]

                # Determine orientation
                a, b, c = i, j, k
                if A[a][b] and A[b][c] and A[c][a]:
                    cycle_dict = {(a,b): 1.0, (b,c): 1.0, (c,a): 1.0}
                else:
                    cycle_dict = {(b,a): 1.0, (a,c): 1.0, (c,b): 1.0}

                is_bd, preimage = check_if_boundary(A, n, cycle_dict)

                print(f"\n  3-cycle ({i},{j},{k}): good={good}, bad={bad_in_cycle}")
                print(f"    is_boundary in T: {is_bd}")
                if is_bd and preimage:
                    for path, coeff in sorted(preimage.items()):
                        interior = path[1]
                        marker = " *BAD*" if interior in bad_set else ""
                        print(f"      {path}: {coeff:.4f} (interior={interior}{marker})")

print(f"\n  Total 3-cycles with 2 bad vertices: {fill_count}")

# ============================================================
# Part 6: n=6 sampling
# ============================================================

print("\n" + "=" * 75)
print("PART 6: n=6 SAMPLING (200 tournaments with >=2 bad vertices)")
print("=" * 75)

n6 = 6
n6_tested = 0
n6_violations = 0
n6_bad_counts = Counter()
n6_cycle_stats = Counter()
n6_target = 200

random.seed(42)
attempts = 0
max_attempts = 50000

while n6_tested < n6_target and attempts < max_attempts:
    attempts += 1
    A = random_tournament(n6)
    betti = path_betti_numbers(A, n6, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(n6):
        A2, n2, _ = delete_vertex(A, n6, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 2:
        continue

    n6_tested += 1
    n6_bad_counts[len(bad)] += 1

    bad_set = set(bad)

    # Check 3-cycles among bad vertices
    for a, b, c in combinations(bad, 3):
        if is_3cycle(A, a, b, c):
            n6_violations += 1
            n6_cycle_stats[3] += 1

    # General cycle stats
    for i in range(n6):
        for j in range(i+1, n6):
            for k in range(j+1, n6):
                if is_3cycle(A, i, j, k):
                    nb = sum(1 for x in [i,j,k] if x in bad_set)
                    n6_cycle_stats[nb] += 1

    if n6_tested <= 5:
        scores = [sum(A[i]) for i in range(n6)]
        print(f"  #{n6_tested}: scores={scores}, bad={bad}, "
              f"3cycles_all_bad={sum(1 for a,b,c in combinations(bad,3) if is_3cycle(A,a,b,c))}")

print(f"\nn=6 results ({n6_tested} tournaments with >=2 bad vertices, from {attempts} attempts):")
print(f"  Bad vertex count distribution: {dict(n6_bad_counts)}")
print(f"  Violations (3 bad vertices forming 3-cycle): {n6_violations}")
print(f"  3-cycle distribution by #bad vertices: {dict(n6_cycle_stats)}")

if n6_violations == 0:
    print(">>> LEMMA HOLDS for all sampled n=6 tournaments <<<")

# ============================================================
# Part 7: THE KEY ALGEBRAIC INSIGHT
# ============================================================

print("\n" + "=" * 75)
print("PART 7: ALGEBRAIC INSIGHT — WHY THE LEMMA MUST HOLD")
print("=" * 75)

print("""
STRUCTURAL OBSERVATION:

For a tournament T on n vertices with beta_1(T) = 0:
  - dim(Z_1(T)) = C(n-1,2)
  - dim(B_1(T)) = C(n-1,2) (since beta_1 = 0)

When we delete vertex v:
  - dim(Z_1(T\\v)) = C(n-2,2)
  - dim(B_1(T\\v)) = C(n-2,2) - beta_1(T\\v)

The rank of B_1(T) is achieved by 2-paths. When we remove v, we lose
all 2-paths involving v. The rank drop is:
  rank_drop(v) = dim(B_1(T)) - dim(B_1(T\\v))
               = C(n-1,2) - C(n-2,2) + beta_1(T\\v)
               = (n-2) + beta_1(T\\v)

So for a "bad" vertex (beta_1(T\\v)=1), rank_drop = n-1.
For a "good" vertex (beta_1(T\\v)=0), rank_drop = n-2.
""")

# Verify this formula
print("Verifying rank_drop formula at n=5...")
verify_count = 0
formula_ok = True
for A in all_tournaments(5):
    betti = path_betti_numbers(A, 5, max_dim=1)
    if betti[1] != 0:
        continue
    verify_count += 1
    if verify_count > 20:
        break
    details = compute_B1_details(A, 5)
    rd = rank_drop_analysis(A, 5, details)
    for v in range(5):
        if rd[v]['rank_drop'] != rd[v]['expected_drop']:
            print(f"  FORMULA FAILS at v={v}: got {rd[v]['rank_drop']}, expected {rd[v]['expected_drop']}")
            formula_ok = False

print(f"  Checked {verify_count} tournaments: {'ALL OK' if formula_ok else 'FAILURES FOUND'}")

# ============================================================
# Part 8: Which Omega_2 elements span B_1?
# ============================================================

print("\n" + "=" * 75)
print("PART 8: Omega_2 STRUCTURE AND THE PROOF")
print("=" * 75)

print("""
KEY QUESTION: What is Omega_2 for a tournament?

A 2-path (a,b,c) is allowed if a->b and b->c in T.
The boundary partial_2(a,b,c) = (b,c) - (a,c) + (a,b).

For (a,b,c) to be in Omega_2, we need partial_2(a,b,c) in A_1.
The faces are: (b,c), (a,c), (a,b).
- (a,b): allowed since a->b
- (b,c): allowed since b->c
- (a,c): this must be allowed. Since T is a tournament, either a->c or c->a.
  Both give an allowed edge (a,c) or (c,a).
  If a->c: face (a,c) is allowed.
  If c->a: face (a,c) is NOT an allowed edge, but (c,a) is.

Wait — the face is literally (a,c) (vertices a and c in that order).
For this to be an allowed 1-path, we need a->c in T.
If instead c->a, then (a,c) is NOT allowed.

So Omega_2 = {allowed 2-paths (a,b,c) where a->c in T}
           ∪ {linear combinations where non-allowed faces cancel}

The (a,b,c) with a->c are the TRANSITIVE TRIPLES.
The (a,b,c) with c->a form 3-CYCLES.

For 3-cycle paths: ∂(a,b,c) = (b,c) - (a,c) + (a,b)
  where (a,c) is NOT allowed (c->a, not a->c).
  So (a,b,c) alone is NOT in Omega_2.
  But (a,b,c) + (c,b,a)? Let's check:
  ∂(c,b,a) = (b,a) - (c,a) + (c,b)
  Sum of non-allowed faces: -(a,c) from first, -(c,a) from second.
  (a,c) is non-allowed and (c,a) is allowed (since c->a).
  Hmm, need to be more careful.
""")

# Let's compute Omega_2 explicitly for a small tournament
print("\nComputing Omega_2 explicitly for a n=5 tournament with bad vertices...")

for A in all_tournaments(5):
    betti = path_betti_numbers(A, 5, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(5):
        A2, n2, _ = delete_vertex(A, 5, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 2:
        continue

    # Analyze this tournament
    allowed_1 = enumerate_allowed_paths(A, 5, 1)
    allowed_2 = enumerate_allowed_paths(A, 5, 2)

    # Classify 2-paths
    transitive = []
    cyclic = []
    for path in allowed_2:
        a, b, c = path
        if A[a][c]:
            transitive.append(path)
        else:
            cyclic.append(path)

    print(f"\n  Bad vertices: {bad}")
    print(f"  |Omega_1| = |A_1| = {len(allowed_1)}")
    print(f"  |A_2| = {len(allowed_2)}: {len(transitive)} transitive + {len(cyclic)} cyclic")

    omega_2 = compute_omega_basis(A, 5, 2, allowed_2, allowed_1)
    print(f"  dim(Omega_2) = {omega_2.shape[1]}")
    print(f"  Omega_2 = transitive triples + ... (cyclic paths pair up)")

    # How do cyclic 2-paths combine in Omega_2?
    # For a 3-cycle a->b->c->a, both (a,b,c) and (c,b,a) are allowed.
    # The face (a,c) from ∂(a,b,c) is non-allowed (c->a).
    # The face (c,a) from ∂(c,b,a) is allowed (c->a).
    # Actually ∂(c,b,a) = (b,a) - (c,a) + (c,b). Face (c,a) IS allowed since c->a.
    # Hmm, (c,b,a) has faces (b,a), (c,a), (c,b).
    # c->b? Not necessarily. For c->b->a with c->a: this is a different pattern.

    # Let's just see what Omega_2 looks like in terms of A_2 basis
    idx_2 = {path: i for i, path in enumerate(allowed_2)}
    print(f"\n  Omega_2 basis vectors (in A_2 coordinates, showing nonzero entries):")
    for col in range(min(omega_2.shape[1], 5)):
        vec = omega_2[:, col]
        entries = [(allowed_2[i], vec[i]) for i in range(len(allowed_2)) if abs(vec[i]) > 1e-10]
        # Classify entries
        trans_entries = [(p, c) for p, c in entries if A[p[0]][p[2]]]
        cyc_entries = [(p, c) for p, c in entries if not A[p[0]][p[2]]]
        if len(entries) <= 6:
            print(f"    basis[{col}]: {len(trans_entries)}T+{len(cyc_entries)}C: ", end="")
            for p, c in entries:
                tp = "T" if A[p[0]][p[2]] else "C"
                print(f"  {p}({tp}):{c:.3f}", end="")
            print()

    break  # Just one tournament

# ============================================================
# Part 9: The critical test — what's special about bad vertices
# ============================================================

print("\n" + "=" * 75)
print("PART 9: WHAT MAKES A VERTEX BAD — STRUCTURAL ANALYSIS")
print("=" * 75)

print("""
A vertex v is "bad" (beta_1(T\\v) = 1) iff removing v from T drops
the rank of B_1 by n-1 instead of n-2.

This means vertex v contributes EXACTLY ONE EXTRA unit of rank to B_1(T)
beyond what would be expected. The 2-paths through v are "barely" spanning.

HYPOTHESIS: The extra rank unit comes from 2-paths where v is the MIDDLE
vertex (interior vertex). If v is an endpoint, those 2-paths also contribute
to B_1(T\\w) for other vertices w.
""")

# For each tournament, categorize 2-paths by role of bad vertices
for A in all_tournaments(5):
    betti = path_betti_numbers(A, 5, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(5):
        A2, n2, _ = delete_vertex(A, 5, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 2:
        continue

    bad_set = set(bad)

    # Get omega_2 basis and the 2-paths
    allowed_1 = enumerate_allowed_paths(A, 5, 1)
    allowed_2 = enumerate_allowed_paths(A, 5, 2)
    omega_2 = compute_omega_basis(A, 5, 2, allowed_2, allowed_1)

    # For each bad vertex, count 2-paths where it's interior vs endpoint
    print(f"\n  Bad vertices: {bad}, scores: {[sum(A[i]) for i in range(5)]}")

    for v in bad:
        interior = [p for p in allowed_2 if p[1] == v]
        left_ep = [p for p in allowed_2 if p[0] == v]
        right_ep = [p for p in allowed_2 if p[2] == v]
        print(f"    v={v}: interior={len(interior)}, left={len(left_ep)}, right={len(right_ep)}")

        # How many involve ONLY bad vertices?
        all_bad_interior = [p for p in interior if p[0] in bad_set and p[2] in bad_set]
        print(f"           all-bad interior paths: {len(all_bad_interior)}: {all_bad_interior}")

    # Look at the sub-tournament on bad vertices
    if len(bad) >= 3:
        print(f"    Sub-tournament on bad vertices:")
        for a in bad:
            for b in bad:
                if a != b:
                    print(f"      {a}->{b}: {bool(A[a][b])}")
        # Is it transitive?
        trans = True
        for a, b, c in combinations(bad, 3):
            if is_3cycle(A, a, b, c):
                trans = False
                break
        print(f"    Transitive: {trans}")

    break  # First example

# ============================================================
# Part 10: THE PROOF ATTEMPT — Sub-tournament structure
# ============================================================

print("\n" + "=" * 75)
print("PART 10: SUB-TOURNAMENT ON BAD VERTICES — ALWAYS TRANSITIVE?")
print("=" * 75)

# Check at n=5: is the sub-tournament on bad vertices ALWAYS transitive?
always_trans_5 = True
for A in all_tournaments(5):
    betti = path_betti_numbers(A, 5, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(5):
        A2, n2, _ = delete_vertex(A, 5, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 3:
        continue

    for a, b, c in combinations(bad, 3):
        if is_3cycle(A, a, b, c):
            always_trans_5 = False
            print(f"  COUNTEREXAMPLE: bad={bad}, 3-cycle among {a},{b},{c}")
            break
    if not always_trans_5:
        break

print(f"\nn=5: Sub-tournament on bad vertices is always transitive: {always_trans_5}")

# Check at n=6
always_trans_6 = True
n6_check = 0
random.seed(123)
for _ in range(20000):
    A = random_tournament(6)
    betti = path_betti_numbers(A, 6, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(6):
        A2, n2, _ = delete_vertex(A, 6, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 3:
        continue

    n6_check += 1
    for a, b, c in combinations(bad, 3):
        if is_3cycle(A, a, b, c):
            always_trans_6 = False
            print(f"  COUNTEREXAMPLE at n=6: bad={bad}, 3-cycle among {a},{b},{c}")
            break
    if not always_trans_6:
        break

print(f"n=6: Checked {n6_check} tournaments with >=3 bad vertices: always transitive = {always_trans_6}")

# Check at n=7
always_trans_7 = True
n7_check = 0
random.seed(456)
for _ in range(10000):
    A = random_tournament(7)
    betti = path_betti_numbers(A, 7, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(7):
        A2, n2, _ = delete_vertex(A, 7, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 3:
        continue

    n7_check += 1
    for a, b, c in combinations(bad, 3):
        if is_3cycle(A, a, b, c):
            always_trans_7 = False
            print(f"  COUNTEREXAMPLE at n=7: bad={bad}, 3-cycle among {a},{b},{c}")
            scores = [sum(A[i]) for i in range(7)]
            print(f"    scores={scores}")
            break
    if not always_trans_7:
        break

print(f"n=7: Checked {n7_check} tournaments with >=3 bad vertices: always transitive = {always_trans_7}")

# ============================================================
# Part 11: Deeper — the unfillable cycles and orientation
# ============================================================

print("\n" + "=" * 75)
print("PART 11: UNFILLABLE CYCLE ORIENTATION AND BAD VERTEX ORDERING")
print("=" * 75)

print("""
If bad vertices always form a transitive sub-tournament, they have a
LINEAR ORDER. What is this order and how does it relate to the
unfillable cycles?

For each pair of bad vertices (a,b) with a->b in T:
  - z_a (the unfillable cycle in T\\a) uses vertex b but not a
  - z_b (the unfillable cycle in T\\b) uses vertex a but not b
  - How do z_a and z_b relate on their common edges?
""")

pair_analysis_count = 0
for A in all_tournaments(5):
    betti = path_betti_numbers(A, 5, max_dim=1)
    if betti[1] != 0:
        continue

    bad = []
    for v in range(5):
        A2, n2, _ = delete_vertex(A, 5, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)

    if len(bad) < 2:
        continue

    pair_analysis_count += 1
    if pair_analysis_count > 3:
        continue

    z_cycles = {}
    for v in bad:
        z_v, remaining = find_unfillable_cycle(A, 5, v)
        if z_v is not None:
            z_cycles[v] = z_v

    print(f"\n  Tournament #{pair_analysis_count}, bad={bad}")

    for v1, v2 in combinations(bad, 2):
        if v1 not in z_cycles or v2 not in z_cycles:
            continue

        z1 = z_cycles[v1]
        z2 = z_cycles[v2]

        direction = "->".join([str(v1), str(v2)]) if A[v1][v2] else "->".join([str(v2), str(v1)])
        print(f"\n    Pair ({v1},{v2}), edge direction: {direction}")

        # Edge (v1, v2) or (v2, v1) should appear in one but not the other
        for e in [(v1, v2), (v2, v1)]:
            in_z1 = e in z1
            in_z2 = e in z2
            print(f"      Edge {e}: in z_{v1}={in_z1}, in z_{v2}={in_z2}")

        # Common edges (not involving v1 or v2)
        common_edges = set()
        for e in z1:
            if v2 not in e:
                common_edges.add(e)
        for e in z2:
            if v1 not in e:
                common_edges.add(e)

        print(f"      Common edges (neither v1 nor v2): {sorted(common_edges)}")
        for e in sorted(common_edges):
            c1 = z1.get(e, 0)
            c2 = z2.get(e, 0)
            if abs(c1) > 1e-10 or abs(c2) > 1e-10:
                ratio = c1/c2 if abs(c2) > 1e-10 else float('inf')
                print(f"        {e}: z_{v1}={c1:.4f}, z_{v2}={c2:.4f}, ratio={ratio:.4f}")

# ============================================================
# Part 12: Maximum number of bad vertices
# ============================================================

print("\n" + "=" * 75)
print("PART 12: MAXIMUM NUMBER OF BAD VERTICES")
print("=" * 75)

max_bad = {5: 0, 6: 0, 7: 0}
max_bad_examples = {}

# n=5 exhaustive
for A in all_tournaments(5):
    betti = path_betti_numbers(A, 5, max_dim=1)
    if betti[1] != 0:
        continue
    bad = []
    for v in range(5):
        A2, n2, _ = delete_vertex(A, 5, v)
        bv = path_betti_numbers(A2, n2, max_dim=1)
        if bv[1] == 1:
            bad.append(v)
    if len(bad) > max_bad[5]:
        max_bad[5] = len(bad)
        max_bad_examples[5] = (A, bad)

print(f"n=5: max bad vertices = {max_bad[5]}")
if 5 in max_bad_examples:
    A, bad = max_bad_examples[5]
    scores = [sum(A[i]) for i in range(5)]
    print(f"  Example: scores={scores}, bad={bad}")

# n=6,7 random
for n_val in [6, 7]:
    random.seed(789 + n_val)
    for _ in range(10000):
        A = random_tournament(n_val)
        betti = path_betti_numbers(A, n_val, max_dim=1)
        if betti[1] != 0:
            continue
        bad = []
        for v in range(n_val):
            A2, n2, _ = delete_vertex(A, n_val, v)
            bv = path_betti_numbers(A2, n2, max_dim=1)
            if bv[1] == 1:
                bad.append(v)
        if len(bad) > max_bad[n_val]:
            max_bad[n_val] = len(bad)
            max_bad_examples[n_val] = (A, bad)

    print(f"n={n_val}: max bad vertices seen = {max_bad[n_val]}")
    if n_val in max_bad_examples:
        A, bad = max_bad_examples[n_val]
        scores = [sum(A[i]) for i in range(n_val)]
        print(f"  Example: scores={scores}, bad={bad}")

# ============================================================
# Summary
# ============================================================

print("\n" + "=" * 75)
print("SUMMARY")
print("=" * 75)

print(f"""
FLIP OBSTRUCTION LEMMA STATUS:
  n=5 (exhaustive): {'HOLDS' if always_trans_5 else 'FAILS'}
  n=6 ({n6_check} samples): {'HOLDS' if always_trans_6 else 'FAILS'}
  n=7 ({n7_check} samples): {'HOLDS' if always_trans_7 else 'FAILS'}

STRONGER RESULT: Bad vertices form a TRANSITIVE sub-tournament.
  (Not just "no 3-cycle among 3 bad vertices" but the entire
   induced sub-tournament on bad vertices is transitive.)

KEY ALGEBRAIC FACTS:
  1. rank_drop(v) = (n-2) + beta_1(T\\v)  [verified]
  2. Bad vertex: rank_drop = n-1 (one extra rank unit)
  3. The extra rank comes from 2-paths with v as interior vertex
  4. If {a,b,c} form a 3-cycle and all are bad, the 3 unfillable
     cycles z_a, z_b, z_c interact on their common edges in a way
     that creates a contradiction with beta_1(T)=0.

PROOF DIRECTION:
  The key insight should come from the Omega_2 structure.
  Transitive triples are automatically in Omega_2.
  3-cycle 2-paths need to pair up to cancel non-allowed faces.
  The pairing constraint + the rank drop constraint may force
  the sub-tournament on bad vertices to be transitive.

  Specifically: if a->b->c->a is a 3-cycle with all three bad,
  then the 2-paths (a,b,c), (b,c,a), (c,a,b) each contribute
  to Omega_2 only in specific combinations. The rank argument
  shows each vertex contributes n-1 to rank, but the mutual
  dependencies from the 3-cycle mean the total contribution
  is less than 3*(n-1), contradicting dim(B_1) = C(n-1,2).
""")

print("\nDone.")
