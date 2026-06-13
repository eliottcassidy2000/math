#!/usr/bin/env python3
"""
beta1_kernel_analysis.py — Kernel and rank-drop analysis for β₁ heredity

Goal: Prove that if β₁(T)=0 for a tournament T on n≥5 vertices,
then there exists v with β₁(T\v)=0.

Analyzes:
1. Kernel K_v = ker(res_v) in Z¹(T) for each vertex v
2. Rank-drop structure and sum constraints
3. U_v space dimensions and intersection with B¹(T\v)
4. Type A/B boundary decomposition and 3-cycle connections

Uses path_homology_v2.py for core homology computation.
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import sys
import random

# Import core functions from path_homology_v2
# We need: enumerate_allowed_paths, boundary_coeffs, build_full_boundary_matrix,
#           compute_omega_basis, path_betti_numbers

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


def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def delete_vertex(A, n, v):
    """Return adjacency matrix of T\v (tournament with v deleted)."""
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    ri = 0
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
    return B


def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3


def compute_Z1_B1(A, n):
    """Compute Z¹ (cocycles) and B¹ (coboundaries) for path homology.

    Returns:
        Z1_basis: columns are basis vectors of Z¹ = ker(∂₁|_{Ω₁})
        B1_basis: columns are basis vectors of B¹ = im(∂₂|_{Ω₂})
        allowed_1: list of allowed 1-paths (edges)
        dim_Z1, dim_B1, beta_1
    """
    allowed_0 = enumerate_allowed_paths(A, n, 0)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_3 = enumerate_allowed_paths(A, n, 3)

    # Omega_1 basis
    omega_1 = compute_omega_basis(A, n, 1, allowed_1, allowed_0)
    dim_omega_1 = omega_1.shape[1] if omega_1.ndim == 2 else 0

    # ∂₁: Ω₁ → A₀
    bd_1 = build_full_boundary_matrix(allowed_1, allowed_0)
    bd_1_omega = bd_1 @ omega_1

    # Z¹ = ker(∂₁|_{Ω₁})
    if bd_1_omega.shape[0] > 0 and bd_1_omega.shape[1] > 0:
        U, S, Vt = np.linalg.svd(bd_1_omega, full_matrices=True)
        rank_1 = sum(s > 1e-8 for s in S)
        # Null space of bd_1_omega in Ω₁ coordinates
        Z1_omega_coords = Vt[rank_1:].T  # columns
        # Convert to A₁ coordinates
        Z1_basis = omega_1 @ Z1_omega_coords
    else:
        Z1_basis = omega_1

    dim_Z1 = Z1_basis.shape[1] if Z1_basis.ndim == 2 else 0

    # Omega_2 basis
    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0

    # B¹ = im(∂₂|_{Ω₂})
    if dim_omega_2 > 0:
        bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
        bd_2_omega = bd_2 @ omega_2  # columns in A₁ coordinates
        # Rank = dim(B¹)
        S2 = np.linalg.svd(bd_2_omega, compute_uv=False)
        dim_B1 = sum(s > 1e-8 for s in S2)
        # Get a basis for B¹ (column space of bd_2_omega)
        U2, S2full, _ = np.linalg.svd(bd_2_omega, full_matrices=False)
        B1_basis = U2[:, :dim_B1]
    else:
        dim_B1 = 0
        B1_basis = np.zeros((len(allowed_1), 0))

    beta_1 = dim_Z1 - dim_B1

    return Z1_basis, B1_basis, allowed_1, dim_Z1, dim_B1, beta_1


def restriction_map(allowed_1_T, allowed_1_Tv, n, v):
    """Build the restriction map res_v: R^{|A₁(T)|} → R^{|A₁(T\v)|}.

    Maps edge (a,b) in T to:
    - 0 if a=v or b=v (edge involves v)
    - the corresponding edge in T\v (with relabeled vertices)
    """
    # Build vertex relabeling: vertices of T\v are 0..n-2,
    # corresponding to T vertices {0,..,n-1}\{v}
    old_to_new = {}
    new_idx = 0
    for i in range(n):
        if i != v:
            old_to_new[i] = new_idx
            new_idx += 1

    idx_Tv = {path: i for i, path in enumerate(allowed_1_Tv)}

    R = np.zeros((len(allowed_1_Tv), len(allowed_1_T)))
    for j, (a, b) in enumerate(allowed_1_T):
        if a == v or b == v:
            continue  # edge involves v, maps to 0
        new_a = old_to_new[a]
        new_b = old_to_new[b]
        new_edge = (new_a, new_b)
        if new_edge in idx_Tv:
            R[idx_Tv[new_edge], j] = 1.0

    return R


def compute_Uv(A, n, v, allowed_1, allowed_2):
    """Compute U_v = span{∂₂(p) : p ∈ Ω₂(T), v ∈ supp(p)}.

    Returns U_v as column vectors in A₁(T) coordinates, and also
    separates Type A (v is middle vertex) and Type B (v is endpoint).
    """
    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0

    if dim_omega_2 == 0:
        return np.zeros((len(allowed_1), 0)), np.zeros((len(allowed_1), 0)), np.zeros((len(allowed_1), 0))

    bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)

    # Get all Ω₂ basis vectors and their boundaries
    # First, identify which allowed 2-paths involve v
    idx_2 = {path: i for i, path in enumerate(allowed_2)}

    # For each Ω₂ basis vector, check if it has support on v-paths
    # A 2-path (a,b,c) involves v if v ∈ {a,b,c}
    # Type A: v = b (middle), Type B: v = a or v = c (endpoint)

    # We need to work in Ω₂ coordinates
    # Each column of omega_2 is a basis vector of Ω₂ in A₂ coordinates

    # For each allowed 2-path, determine if it involves v and how
    v_paths_all = []  # indices into allowed_2
    v_paths_A = []    # v is middle
    v_paths_B = []    # v is endpoint

    for i, (a, b, c) in enumerate(allowed_2):
        if v == b:
            v_paths_A.append(i)
            v_paths_all.append(i)
        elif v == a or v == c:
            v_paths_B.append(i)
            v_paths_all.append(i)

    # Project Ω₂ onto v-path components
    # U_v = im(∂₂ restricted to v-path part of Ω₂)
    # This is tricky because Ω₂ elements can mix v-paths and non-v-paths.

    # Actually, U_v = span{∂₂(p) : p ∈ Ω₂, v ∈ supp(p)}
    # But Ω₂ elements are linear combinations of allowed 2-paths.
    # We want boundaries of Ω₂ elements that have support on v-containing paths.

    # Simpler approach: compute all ∂₂(omega_2 basis vectors) and separate
    # the v-involved part.

    # Actually the most faithful interpretation:
    # For each basis element of Ω₂ that has ANY nonzero coefficient on a v-path,
    # include its boundary in U_v.
    # But that's too broad. Let me re-read the user's intent.

    # "U_v = span{∂₂(p) : p ∈ Ω₂(T), v ∈ supp(p)}"
    # This means: for each element p of Ω₂ whose support involves v,
    # take its boundary. U_v is the span of all such boundaries.

    # But every Ω₂ element is a linear combination. "v ∈ supp(p)" means
    # the element has a nonzero coefficient on some path containing v.

    # Actually, the simplest interpretation:
    # For each ALLOWED 2-path (a,b,c) that is also in Ω₂ (as a single path,
    # which it is if Ω₂ = A₂, which happens for tournaments):
    # If v ∈ {a,b,c}, include ∂₂(a,b,c) in U_v.

    # For tournaments, Ω₂ = A₂ (all transitive triples), so we can just
    # take boundaries of individual allowed 2-paths involving v.

    all_v_boundaries = []
    typeA_boundaries = []
    typeB_boundaries = []

    idx_1 = {path: i for i, path in enumerate(allowed_1)}

    for i, (a, b, c) in enumerate(allowed_2):
        if v not in (a, b, c):
            continue
        # Boundary: (b,c) - (a,c) + (a,b)
        bvec = np.zeros(len(allowed_1))
        for sign, face in boundary_coeffs((a, b, c)):
            if face in idx_1:
                bvec[idx_1[face]] += sign

        all_v_boundaries.append(bvec)
        if v == b:
            typeA_boundaries.append(bvec)
        else:
            typeB_boundaries.append(bvec)

    def to_matrix(vecs):
        if not vecs:
            return np.zeros((len(allowed_1), 0))
        return np.column_stack(vecs)

    return to_matrix(all_v_boundaries), to_matrix(typeA_boundaries), to_matrix(typeB_boundaries)


def matrix_rank(M, tol=1e-8):
    if M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    S = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(S > tol))


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


# ========================================================================
print("=" * 70)
print("BETA_1 KERNEL ANALYSIS — Heredity Proof Investigation")
print("=" * 70)

# ========================================================================
# SECTION 1: Kernel K_v analysis at n=5
# ========================================================================
print("\n" + "=" * 70)
print("SECTION 1: Kernel K_v = ker(res_v|_{Z¹(T)}) at n=5")
print("=" * 70)

n = 5
beta0_count = 0
beta1_count = 0

kernel_data = []  # Store for analysis

for A in all_tournaments(n):
    betti = path_betti_numbers(A, n, max_dim=1)
    if betti[1] != 0:
        beta1_count += 1
        continue

    beta0_count += 1

    Z1, B1, allowed_1, dim_Z1, dim_B1, beta_1 = compute_Z1_B1(A, n)

    # For each vertex v, compute kernel K_v
    kernels = {}
    beta1_Tv = {}

    for v in range(n):
        Av = delete_vertex(A, n, v)
        betti_v = path_betti_numbers(Av, n-1, max_dim=1)
        beta1_Tv[v] = betti_v[1]

        # Restriction map
        allowed_1_Tv = enumerate_allowed_paths(Av, n-1, 1)
        R = restriction_map(allowed_1, allowed_1_Tv, n, v)

        # res_v restricted to Z¹(T): R @ Z1
        R_Z1 = R @ Z1  # columns: images of Z¹ basis vectors in A₁(T\v)

        # Kernel of R_Z1
        if R_Z1.shape[1] > 0:
            U, S, Vt = np.linalg.svd(R_Z1, full_matrices=True)
            rank_R = sum(s > 1e-8 for s in S)
            ker_coords = Vt[rank_R:].T  # in Z¹ coordinates
            # Convert to A₁(T) coordinates
            if ker_coords.shape[1] > 0:
                ker_A1 = Z1 @ ker_coords
            else:
                ker_A1 = np.zeros((len(allowed_1), 0))
        else:
            ker_A1 = Z1  # everything is in kernel

        kernels[v] = ker_A1

    # Analysis
    bad_verts = [v for v in range(n) if beta1_Tv[v] > 0]
    good_verts = [v for v in range(n) if beta1_Tv[v] == 0]

    entry = {
        'A': A,
        'dim_Z1': dim_Z1,
        'kernels': kernels,
        'beta1_Tv': beta1_Tv,
        'bad_verts': bad_verts,
        'good_verts': good_verts,
        't3': count_3cycles(A, n),
    }
    kernel_data.append(entry)

print(f"\nn=5: {beta0_count} tournaments with β₁=0, {beta1_count} with β₁=1")

# 1a: Are kernels always 1-dimensional?
all_1d = True
for entry in kernel_data:
    for v in range(n):
        dim_k = entry['kernels'][v].shape[1]
        if dim_k != 1:
            all_1d = False
            print(f"  SURPRISE: K_{v} has dim {dim_k} (t3={entry['t3']})")

print(f"\n1a) All kernels K_v are 1-dimensional: {all_1d}")

# 1b: Do any K_a = K_b (same 1-d subspace)?
coincidence_count = 0
total_pairs = 0
for entry in kernel_data:
    for a in range(n):
        for b in range(a+1, n):
            total_pairs += 1
            ka = entry['kernels'][a]
            kb = entry['kernels'][b]
            if ka.shape[1] == 1 and kb.shape[1] == 1:
                # Check if they span the same line
                combined = np.column_stack([ka, kb])
                r = matrix_rank(combined)
                if r == 1:
                    coincidence_count += 1

print(f"\n1b) Kernel coincidences K_a = K_b: {coincidence_count} / {total_pairs} pairs")

# 1c: Linear relations among K_v
print(f"\n1c) Linear relations among the n={n} kernel vectors in Z¹(T) (dim={n-1}):")
# n points in (n-1)-dim space: they span at most (n-1)-dim.
# So n=5 points in 4-dim space can be in general position (rank 4 from 5 vectors)
rank_dist = Counter()
for entry in kernel_data:
    if all(entry['kernels'][v].shape[1] == 1 for v in range(n)):
        K_matrix = np.column_stack([entry['kernels'][v][:, 0] for v in range(n)])
        r = matrix_rank(K_matrix)
        rank_dist[r] += 1

for r in sorted(rank_dist.keys()):
    print(f"  rank({'{K_0,...,K_4}'}) = {r}: {rank_dist[r]} tournaments")

# 1d: For bad vertices, are their kernels collinear?
print(f"\n1d) Kernel geometry of bad vertices:")
bad_rank_dist = Counter()
n_bad_dist = Counter()
for entry in kernel_data:
    bv = entry['bad_verts']
    n_bad_dist[len(bv)] += 1
    if len(bv) >= 2 and all(entry['kernels'][v].shape[1] == 1 for v in bv):
        K_bad = np.column_stack([entry['kernels'][v][:, 0] for v in bv])
        r = matrix_rank(K_bad)
        bad_rank_dist[(len(bv), r)] += 1

print(f"  Number of bad vertices distribution: {dict(n_bad_dist)}")
for (nb, r), cnt in sorted(bad_rank_dist.items()):
    print(f"  {nb} bad vertices, rank of their kernels = {r}: {cnt} tournaments")

# Check if bad vertices always form specific structures
print(f"\n  Bad vertex structure details:")
bad_structure_count = Counter()
for entry in kernel_data:
    bv = entry['bad_verts']
    if len(bv) == 0:
        bad_structure_count['all_good'] += 1
    else:
        # Check if bad vertices form a transitive triple
        A = entry['A']
        if len(bv) == 3:
            a, b, c = bv
            # Check 3-cycle vs transitive
            edges_among = sum(A[i][j] for i in bv for j in bv if i != j)
            is_cycle = (edges_among == 3)  # 3-cycle has 3 directed edges
            bad_structure_count[f'3bad_cycle={is_cycle}'] += 1
        else:
            bad_structure_count[f'{len(bv)}bad'] += 1

for key, cnt in sorted(bad_structure_count.items()):
    print(f"    {key}: {cnt}")


# ========================================================================
# SECTION 2: Rank-drop analysis and "all bad" obstruction
# ========================================================================
print("\n\n" + "=" * 70)
print("SECTION 2: Rank-drop analysis and all-bad obstruction")
print("=" * 70)

for n in [5, 6]:
    print(f"\n--- n={n} ---")

    if n == 5:
        tournaments = list(all_tournaments(n))
    else:
        # Sample for n=6
        tournaments = [random_tournament(n) for _ in range(200)]

    beta0_cases = []

    for A in tournaments:
        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue

        Z1, B1, allowed_1, dim_Z1, dim_B1, beta_1 = compute_Z1_B1(A, n)

        allowed_2 = enumerate_allowed_paths(A, n, 2)

        vertex_data = {}
        for v in range(n):
            Av = delete_vertex(A, n, v)
            betti_v = path_betti_numbers(Av, n-1, max_dim=1)

            # Compute U_v
            Uv_all, Uv_A, Uv_B = compute_Uv(A, n, v, allowed_1, allowed_2)
            dim_Uv = matrix_rank(Uv_all)
            dim_Uv_A = matrix_rank(Uv_A)
            dim_Uv_B = matrix_rank(Uv_B)

            # B¹(T\v) — compute
            Z1v, B1v, allowed_1v, dim_Z1v, dim_B1v, beta1v = compute_Z1_B1(Av, n-1)

            # res_v(Z¹(T)) and its relationship to Z¹(T\v), B¹(T\v)
            allowed_1_Tv = enumerate_allowed_paths(Av, n-1, 1)
            R = restriction_map(allowed_1, allowed_1_Tv, n, v)
            R_Z1 = R @ Z1
            rank_res = matrix_rank(R_Z1)

            # U_v ∩ B¹(T\v): project U_v through restriction, intersect with B¹(T\v)
            # Actually U_v lives in A₁(T), and B¹(T\v) lives in A₁(T\v).
            # The rank drop is: dim(B₁(T)) - dim(B₁(T\v)) when decomposing
            # But more precisely:
            # B₁(T) = B₁(T\v)_lifted + U_v
            # where B₁(T\v)_lifted is B₁(T\v) embedded in A₁(T) via the inclusion.

            # The image of U_v under res_v
            R_Uv = R @ Uv_all
            dim_R_Uv = matrix_rank(R_Uv)

            # Intersection of im(res_v(U_v)) with B¹(T\v)
            if dim_R_Uv > 0 and dim_B1v > 0:
                combined = np.column_stack([R_Uv, B1v])
                r_combined = matrix_rank(combined)
                dim_intersect = dim_R_Uv + dim_B1v - r_combined
            else:
                dim_intersect = 0

            rank_drop = rank_res - (dim_Z1 - 1)  # Should be 0 normally

            vertex_data[v] = {
                'beta1': betti_v[1],
                'dim_Uv': dim_Uv,
                'dim_Uv_A': dim_Uv_A,
                'dim_Uv_B': dim_Uv_B,
                'dim_Z1v': dim_Z1v,
                'dim_B1v': dim_B1v,
                'rank_res': rank_res,
                'dim_R_Uv': dim_R_Uv,
                'dim_intersect': dim_intersect,
            }

        bad_verts = [v for v in range(n) if vertex_data[v]['beta1'] > 0]
        good_verts = [v for v in range(n) if vertex_data[v]['beta1'] == 0]

        beta0_cases.append({
            'A': A,
            'vertex_data': vertex_data,
            'bad_verts': bad_verts,
            'good_verts': good_verts,
            'dim_Z1': dim_Z1,
            'dim_B1': dim_B1,
            'dim_Omega2': len(allowed_2),
            't3': count_3cycles(A, n),
        })

    print(f"  Found {len(beta0_cases)} tournaments with β₁=0")

    # 2a: Distribution of number of bad vertices
    bad_dist = Counter()
    for case in beta0_cases:
        bad_dist[len(case['bad_verts'])] += 1
    print(f"\n  2a) Bad vertex count distribution: {dict(sorted(bad_dist.items()))}")

    # 2b: Sum of rank drops = n(n-2) + Σβ₁(T\v) ?
    print(f"\n  2b) Rank-drop sum analysis:")
    for case in beta0_cases[:5]:  # First 5 cases
        sum_rank_res = sum(case['vertex_data'][v]['rank_res'] for v in range(n))
        sum_beta1 = sum(case['vertex_data'][v]['beta1'] for v in range(n))
        print(f"    t3={case['t3']}: Σ rank_res = {sum_rank_res}, "
              f"dim_Z1 = {case['dim_Z1']}, Σβ₁(T\\v) = {sum_beta1}, "
              f"|bad| = {len(case['bad_verts'])}")

    # 2c: dim(U_v) statistics
    print(f"\n  2c) dim(U_v) statistics:")
    Uv_by_status = {'good': [], 'bad': []}
    for case in beta0_cases:
        for v in range(n):
            status = 'bad' if v in case['bad_verts'] else 'good'
            Uv_by_status[status].append(case['vertex_data'][v]['dim_Uv'])

    for status in ['good', 'bad']:
        vals = Uv_by_status[status]
        if vals:
            c = Counter(vals)
            print(f"    {status} vertices: dim(U_v) distribution = {dict(sorted(c.items()))}")

    # 2d: Type A vs Type B
    print(f"\n  2d) Type A (v middle) vs Type B (v endpoint) dimensions:")
    for status in ['good', 'bad']:
        typeA_vals = []
        typeB_vals = []
        for case in beta0_cases:
            for v in range(n):
                s = 'bad' if v in case['bad_verts'] else 'good'
                if s == status:
                    typeA_vals.append(case['vertex_data'][v]['dim_Uv_A'])
                    typeB_vals.append(case['vertex_data'][v]['dim_Uv_B'])
        if typeA_vals:
            print(f"    {status}: Type A = {dict(sorted(Counter(typeA_vals).items()))}, "
                  f"Type B = {dict(sorted(Counter(typeB_vals).items()))}")

    # 2e: res_v(U_v) ∩ B¹(T\v) dimension
    print(f"\n  2e) dim(res_v(U_v) ∩ B¹(T\\v)) by vertex status:")
    for status in ['good', 'bad']:
        vals = []
        for case in beta0_cases:
            for v in range(n):
                s = 'bad' if v in case['bad_verts'] else 'good'
                if s == status:
                    vals.append(case['vertex_data'][v]['dim_intersect'])
        if vals:
            print(f"    {status}: {dict(sorted(Counter(vals).items()))}")

    # 2f: "All bad" test — theoretical bound
    print(f"\n  2f) If all vertices were bad:")
    total_Omega2 = [case['dim_Omega2'] for case in beta0_cases]
    print(f"    dim(Ω₂) range: {min(total_Omega2)} to {max(total_Omega2)}")
    print(f"    Each 2-path involves 3 vertices, so Σ_v |v-paths in Ω₂| = 3*dim(Ω₂)")
    print(f"    For all-bad: need Σ rank_drop = n*(n-1) = {n*(n-1)}")
    print(f"    Available 'rank budget' from 3*dim(Ω₂): {3*min(total_Omega2)} to {3*max(total_Omega2)}")

    # 2g: Actual Σ dim(U_v) and Σ rank_drop
    print(f"\n  2g) Actual sums:")
    for case in beta0_cases[:10]:
        sum_dim_Uv = sum(case['vertex_data'][v]['dim_Uv'] for v in range(n))
        sum_beta1 = sum(case['vertex_data'][v]['beta1'] for v in range(n))
        sum_rank_res = sum(case['vertex_data'][v]['rank_res'] for v in range(n))
        sum_dim_intersect = sum(case['vertex_data'][v]['dim_intersect'] for v in range(n))

        print(f"    t3={case['t3']}, |bad|={len(case['bad_verts'])}: "
              f"Σdim(U_v)={sum_dim_Uv}, Σrank_res={sum_rank_res}, "
              f"Σβ₁(T\\v)={sum_beta1}, Σintersect={sum_dim_intersect}")


# ========================================================================
# SECTION 3: Direct attempt — per-vertex dimension accounting
# ========================================================================
print("\n\n" + "=" * 70)
print("SECTION 3: Per-vertex dimension accounting")
print("=" * 70)

n = 5
print(f"\n--- n={n} detailed dimension table ---")
print(f"{'v':>3} {'β₁':>3} {'dimZ1v':>6} {'dimB1v':>6} {'dimUv':>5} "
      f"{'rkRes':>5} {'dimRUv':>6} {'intct':>5} {'UvA':>4} {'UvB':>4}")

for case in kernel_data[:10]:  # Use kernel_data from section 1
    A = case['A']
    t3 = case['t3']

    Z1, B1, allowed_1, dim_Z1, dim_B1, _ = compute_Z1_B1(A, n)
    allowed_2 = enumerate_allowed_paths(A, n, 2)

    print(f"\n  Tournament t3={t3}, bad={case['bad_verts']}, good={case['good_verts']}")

    for v in range(n):
        Av = delete_vertex(A, n, v)
        Z1v, B1v, allowed_1v, dim_Z1v, dim_B1v, beta1v = compute_Z1_B1(Av, n-1)

        Uv_all, Uv_A, Uv_B = compute_Uv(A, n, v, allowed_1, allowed_2)
        dim_Uv = matrix_rank(Uv_all)
        dim_Uv_A = matrix_rank(Uv_A)
        dim_Uv_B = matrix_rank(Uv_B)

        allowed_1_Tv = enumerate_allowed_paths(Av, n-1, 1)
        R = restriction_map(allowed_1, allowed_1_Tv, n, v)
        R_Z1 = R @ Z1
        rank_res = matrix_rank(R_Z1)

        R_Uv = R @ Uv_all
        dim_R_Uv = matrix_rank(R_Uv)

        if dim_R_Uv > 0 and dim_B1v > 0:
            combined = np.column_stack([R_Uv, B1v])
            r_combined = matrix_rank(combined)
            dim_intersect = dim_R_Uv + dim_B1v - r_combined
        else:
            dim_intersect = 0

        marker = " BAD" if case['beta1_Tv'][v] > 0 else ""
        print(f"  {v:>3} {case['beta1_Tv'][v]:>3} {dim_Z1v:>6} {dim_B1v:>6} "
              f"{dim_Uv:>5} {rank_res:>5} {dim_R_Uv:>6} {dim_intersect:>5} "
              f"{dim_Uv_A:>4} {dim_Uv_B:>4}{marker}")


# ========================================================================
# SECTION 4: Vertex score and 3-cycle connection
# ========================================================================
print("\n\n" + "=" * 70)
print("SECTION 4: Vertex properties vs good/bad status")
print("=" * 70)

n = 5
print(f"\n--- n={n}: Vertex properties ---")

score_good = Counter()
score_bad = Counter()
t3_through_good = []
t3_through_bad = []

for case in kernel_data:
    A = case['A']
    for v in range(n):
        # Out-degree of v
        out_deg = sum(A[v])

        # Number of 3-cycles through v
        t3v = 0
        for i in range(n):
            if i == v:
                continue
            for j in range(i+1, n):
                if j == v:
                    continue
                verts = [v, i, j]
                edges = sum(A[a][b] for a in verts for b in verts if a != b)
                if edges == 3:  # exactly 3 directed edges = 3-cycle
                    t3v += 1

        is_bad = case['beta1_Tv'][v] > 0
        if is_bad:
            score_bad[out_deg] += 1
            t3_through_bad.append(t3v)
        else:
            score_good[out_deg] += 1
            t3_through_good.append(t3v)

print(f"\n  Out-degree distribution:")
print(f"    Good vertices: {dict(sorted(score_good.items()))}")
print(f"    Bad vertices:  {dict(sorted(score_bad.items()))}")

print(f"\n  3-cycles through vertex:")
print(f"    Good: min={min(t3_through_good)}, max={max(t3_through_good)}, "
      f"mean={np.mean(t3_through_good):.2f}")
if t3_through_bad:
    print(f"    Bad:  min={min(t3_through_bad)}, max={max(t3_through_bad)}, "
          f"mean={np.mean(t3_through_bad):.2f}")

# Check: is good/bad related to vertex being in a "dominant" position?
print(f"\n  Source/sink analysis:")
source_good = 0
source_bad = 0
sink_good = 0
sink_bad = 0
for case in kernel_data:
    A = case['A']
    for v in range(n):
        out_deg = sum(A[v])
        in_deg = n - 1 - out_deg
        is_bad = case['beta1_Tv'][v] > 0
        if out_deg == n-1:  # source
            if is_bad: source_bad += 1
            else: source_good += 1
        if in_deg == n-1:  # sink
            if is_bad: sink_bad += 1
            else: sink_good += 1

print(f"    Sources (out-deg={n-1}): good={source_good}, bad={source_bad}")
print(f"    Sinks (in-deg={n-1}): good={sink_good}, bad={sink_bad}")


# ========================================================================
# SECTION 5: Key relationship — when is res_v surjective onto Z¹(T\v)?
# ========================================================================
print("\n\n" + "=" * 70)
print("SECTION 5: res_v surjectivity analysis")
print("=" * 70)

for n in [5, 6]:
    print(f"\n--- n={n} ---")

    if n == 5:
        tournaments = list(all_tournaments(n))
    else:
        tournaments = [random_tournament(n) for _ in range(200)]

    surj_stats = {'good_surj': 0, 'good_not_surj': 0,
                  'bad_surj': 0, 'bad_not_surj': 0}
    corank_dist = {'good': Counter(), 'bad': Counter()}

    for A in tournaments:
        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue

        Z1, B1, allowed_1, dim_Z1, dim_B1, _ = compute_Z1_B1(A, n)

        for v in range(n):
            Av = delete_vertex(A, n, v)
            betti_v = path_betti_numbers(Av, n-1, max_dim=1)
            is_bad = betti_v[1] > 0

            Z1v, B1v, allowed_1v, dim_Z1v, dim_B1v, beta1v = compute_Z1_B1(Av, n-1)

            allowed_1_Tv = enumerate_allowed_paths(Av, n-1, 1)
            R = restriction_map(allowed_1, allowed_1_Tv, n, v)

            # Image of Z¹(T) under res_v
            R_Z1 = R @ Z1
            rank_image = matrix_rank(R_Z1)

            # Is the image = Z¹(T\v)?
            if dim_Z1v > 0:
                combined = np.column_stack([R_Z1, Z1v]) if R_Z1.shape[1] > 0 else Z1v
                rank_combined = matrix_rank(combined)
                corank = dim_Z1v - (rank_image + dim_Z1v - rank_combined)
                # Actually: corank = dim(Z¹(T\v)) - dim(res_v(Z¹(T)) ∩ Z¹(T\v))
                # Let's be more careful
                # res_v(Z¹(T)) is a subspace of A₁(T\v)
                # We want: does it contain Z¹(T\v)?

                # Project Z¹(T\v) onto the image
                if R_Z1.shape[1] > 0:
                    # Check if Z¹(T\v) ⊆ im(R_Z1)
                    combined = np.column_stack([R_Z1, Z1v])
                    rank_combined = matrix_rank(combined)
                    rank_R_Z1 = matrix_rank(R_Z1)
                    surjective = (rank_combined == rank_R_Z1)
                    corank = rank_combined - rank_R_Z1  # dim of Z¹(T\v) not in image
                else:
                    surjective = (dim_Z1v == 0)
                    corank = dim_Z1v
            else:
                surjective = True
                corank = 0

            status = 'bad' if is_bad else 'good'
            if surjective:
                surj_stats[f'{status}_surj'] += 1
            else:
                surj_stats[f'{status}_not_surj'] += 1
            corank_dist[status][corank] += 1

    print(f"  Surjectivity of res_v: Z¹(T) → Z¹(T\\v):")
    for key, val in sorted(surj_stats.items()):
        print(f"    {key}: {val}")

    print(f"\n  Corank distribution (dim Z¹(T\\v) not hit by res_v(Z¹(T))):")
    for status in ['good', 'bad']:
        print(f"    {status}: {dict(sorted(corank_dist[status].items()))}")


# ========================================================================
# SECTION 6: The critical formula — B₁(T) = B₁(T\v) + U_v decomposition
# ========================================================================
print("\n\n" + "=" * 70)
print("SECTION 6: B₁ decomposition B₁(T) = B₁(T\\v)_lifted + U_v")
print("=" * 70)

n = 5
print(f"\n--- n={n}: B₁ decomposition ---")

for case_idx, case in enumerate(kernel_data[:5]):
    A = case['A']
    t3 = case['t3']

    Z1, B1, allowed_1, dim_Z1, dim_B1, _ = compute_Z1_B1(A, n)
    allowed_2 = enumerate_allowed_paths(A, n, 2)

    print(f"\n  Tournament {case_idx} (t3={t3}, bad={case['bad_verts']}):")
    print(f"    dim B₁(T) = {dim_B1}, dim Z₁(T) = {dim_Z1}")

    for v in range(n):
        Av = delete_vertex(A, n, v)
        Z1v, B1v, allowed_1v, dim_Z1v, dim_B1v, beta1v = compute_Z1_B1(Av, n-1)

        # U_v: boundaries from paths through v
        Uv_all, _, _ = compute_Uv(A, n, v, allowed_1, allowed_2)
        dim_Uv = matrix_rank(Uv_all)

        # Lift B₁(T\v) into A₁(T): for each edge (a',b') in T\v,
        # map to (a,b) in T where a,b are the original labels
        allowed_1_Tv = enumerate_allowed_paths(Av, n-1, 1)

        # Inverse of restriction map: embed T\v edges into T
        new_to_old = []
        for i in range(n):
            if i != v:
                new_to_old.append(i)

        idx_T = {path: i for i, path in enumerate(allowed_1)}

        if dim_B1v > 0:
            # B₁(T\v) lives in A₁(T\v). Lift to A₁(T).
            E = np.zeros((len(allowed_1), len(allowed_1_Tv)))
            for j, (a_new, b_new) in enumerate(allowed_1_Tv):
                a_old = new_to_old[a_new]
                b_old = new_to_old[b_new]
                if (a_old, b_old) in idx_T:
                    E[idx_T[(a_old, b_old)], j] = 1.0

            B1v_lifted = E @ B1v
            dim_B1v_lifted = matrix_rank(B1v_lifted)
        else:
            B1v_lifted = np.zeros((len(allowed_1), 0))
            dim_B1v_lifted = 0

        # Combined space
        parts = []
        if dim_B1v_lifted > 0:
            parts.append(B1v_lifted)
        if dim_Uv > 0:
            parts.append(Uv_all)

        if parts:
            combined = np.column_stack(parts)
            dim_combined = matrix_rank(combined)
        else:
            dim_combined = 0

        # Check: does B₁(T\v)_lifted + U_v = B₁(T)?
        if dim_combined > 0 and dim_B1 > 0:
            full = np.column_stack([combined, B1])
            dim_full = matrix_rank(full)
            equals_B1 = (dim_full == dim_B1) and (dim_combined == dim_B1)
        else:
            equals_B1 = (dim_B1 == 0)

        # Intersection dimension
        dim_intersect = dim_B1v_lifted + dim_Uv - dim_combined

        marker = " BAD" if case['beta1_Tv'][v] > 0 else ""
        print(f"    v={v}: B1v={dim_B1v}, B1v_lift={dim_B1v_lifted}, "
              f"dim_Uv={dim_Uv}, sum={dim_combined}, "
              f"B1v_lift+Uv=B1? {equals_B1}, "
              f"intersect={dim_intersect}{marker}")


print("\n\n" + "=" * 70)
print("DONE — Analysis complete")
print("=" * 70)
