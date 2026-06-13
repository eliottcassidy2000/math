#!/usr/bin/env python3
"""
beta2_transitive_mechanism.py — Investigate WHY bad vertices form transitive triples

When β₁(T)=0 and deleting vertex v gives β₁(T\v)=1, v is "bad".
This script investigates:
1. What is z_v (the H₁ generator in T\v) concretely?
2. What is the filling chain c_v in T with ∂₂(c_v) = z_v?
3. Why are bad vertices always transitive?
4. Why can't there be 4 bad vertices?
5. Dimension/rank analysis of restriction maps

Uses path_homology_v2.py functions (copied to avoid module-level execution).

opus-2026-03-08
"""

import numpy as np
from itertools import permutations, combinations
from math import comb
from collections import defaultdict, Counter
import random
import sys

# ============================================================
# CORE PATH HOMOLOGY FUNCTIONS (from path_homology_v2.py)
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

def path_homology_detailed(A, n, p_target=1):
    """Return detailed info: betti, Z_p basis (cycles), B_p basis (boundaries),
    and a representative generator if β_p > 0.

    Returns dict with keys: betti, allowed_p, omega_basis, cycle_basis, boundary_basis,
    generator_coords, generator_paths
    """
    max_dim = p_target + 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    p = p_target
    dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0

    result = {
        'betti': 0,
        'allowed_p': allowed[p],
        'allowed_p_plus1': allowed[p+1],
        'omega_basis': omega[p],
        'dim_omega': dim_omega_p,
        'cycle_basis': None,
        'boundary_basis': None,
        'generators': [],
    }

    if dim_omega_p == 0:
        return result

    # ∂_p restricted to Ω_p
    bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
    bd_p_omega = bd_p @ omega[p]

    # Find kernel of ∂_p on Ω_p (= Z_p in Ω_p coords)
    if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
        U, S_vals, Vt = np.linalg.svd(bd_p_omega, full_matrices=True)
        rank_p = sum(s > 1e-8 for s in S_vals)
        ker_basis_omega = Vt[rank_p:].T  # columns = kernel vectors in Ω_p coords
    else:
        rank_p = 0
        ker_basis_omega = np.eye(dim_omega_p)

    ker_dim = ker_basis_omega.shape[1] if ker_basis_omega.ndim == 2 else 0

    # Convert kernel basis to A_p coords
    if ker_dim > 0:
        cycle_basis_Ap = omega[p] @ ker_basis_omega  # columns in A_p coords
    else:
        cycle_basis_Ap = np.zeros((len(allowed[p]), 0))

    result['cycle_basis'] = cycle_basis_Ap

    # Image of ∂_{p+1} in A_p
    dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
    if dim_omega_p1 > 0:
        bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
        bd_p1_omega = bd_p1 @ omega[p+1]
        S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
        im_dim = sum(s > 1e-8 for s in S_p1)
        result['boundary_basis'] = bd_p1_omega
        result['im_dim'] = im_dim
    else:
        im_dim = 0
        result['boundary_basis'] = np.zeros((len(allowed[p]), 0))
        result['im_dim'] = 0

    beta = ker_dim - im_dim
    result['betti'] = max(0, beta)

    # Extract generators: project cycle basis onto complement of boundary image
    if beta > 0 and ker_dim > 0:
        # We have cycle_basis_Ap (columns are cycles in A_p coords)
        # and boundary image bd_p1_omega (columns are boundaries in A_p coords)
        # Generator = any cycle not in boundary image

        # Stack boundaries and cycles, find which cycles are NOT in span of boundaries
        if im_dim > 0:
            B = result['boundary_basis']
            # For each cycle basis vector, check if it's in span(B)
            # Use QR of B to get orthogonal basis of im
            Q, R = np.linalg.qr(B, mode='reduced')
            # Project each cycle onto complement of im(B)
            projections = []
            for i in range(ker_dim):
                c = cycle_basis_Ap[:, i]
                c_proj = c - Q @ (Q.T @ c)  # component orthogonal to boundaries
                if np.linalg.norm(c_proj) > 1e-8:
                    projections.append(c)
            result['generators'] = projections[:beta]
        else:
            result['generators'] = [cycle_basis_Ap[:, i] for i in range(min(beta, ker_dim))]

    return result


# ============================================================
# TOURNAMENT UTILITIES
# ============================================================

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def delete_vertex(A, n, v):
    """Return adjacency matrix with vertex v deleted, and the vertex relabeling map."""
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    old_to_new = {}
    j = 0
    for i in range(n):
        if i == v:
            continue
        old_to_new[i] = j
        j += 1
    for i in range(n):
        if i == v:
            continue
        for k in range(n):
            if k == v:
                continue
            B[old_to_new[i]][old_to_new[k]] = A[i][k]
    return B, new_n, old_to_new

def tournament_to_string(A, n):
    """Compact string representation."""
    bits = []
    for i in range(n):
        for j in range(i+1, n):
            bits.append('1' if A[i][j] else '0')
    return ''.join(bits)

def is_transitive_set(A, vertices):
    """Check if the sub-tournament on given vertices is transitive."""
    vlist = sorted(vertices)
    k = len(vlist)
    if k <= 1:
        return True
    # Check all triples
    for i in range(k):
        for j in range(i+1, k):
            for l in range(j+1, k):
                a, b, c = vlist[i], vlist[j], vlist[l]
                # Count 3-cycles
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[b][a] and A[c][b] and A[a][c]):
                    return False
    return True

def cycle_representation(gen_vec, allowed_paths):
    """Convert a generator vector to a readable sum of edges/paths."""
    terms = []
    for i, coeff in enumerate(gen_vec):
        if abs(coeff) > 1e-8:
            terms.append((coeff, allowed_paths[i]))
    return terms

def format_cycle(terms):
    """Pretty-print a 1-cycle as sum of directed edges."""
    parts = []
    for coeff, path in sorted(terms, key=lambda x: x[1]):
        c = round(coeff, 4)
        if abs(c - round(c)) < 1e-6:
            c = int(round(c))
        parts.append(f"{c:+}*({path[0]}→{path[1]})")
    return " ".join(parts)


# ============================================================
# QUESTION 1: What is z_v concretely?
# ============================================================

def analyze_generator(A, n, v, verbose=True):
    """For tournament T on n vertices with β₁(T)=0, β₁(T\\v)=1,
    compute the H₁ generator z_v of T\\v."""
    B, m, old_to_new = delete_vertex(A, n, v)
    new_to_old = {new: old for old, new in old_to_new.items()}

    info = path_homology_detailed(B, m, p_target=1)

    if info['betti'] != 1:
        return None

    if not info['generators']:
        return None

    gen = info['generators'][0]
    allowed_1 = info['allowed_p']

    # Normalize: make the largest absolute coefficient = 1
    max_abs = max(abs(c) for c in gen)
    if max_abs > 1e-10:
        gen = gen / max_abs

    # Round to integers if close
    gen_int = np.round(gen).astype(int)
    if np.allclose(gen, gen_int, atol=1e-6):
        gen = gen_int

    terms = cycle_representation(gen, allowed_1)

    # Translate back to original vertex labels
    terms_orig = []
    for coeff, path in terms:
        orig_path = tuple(new_to_old[x] for x in path)
        terms_orig.append((coeff, orig_path))

    # Which vertices does z_v use?
    verts_used = set()
    for _, path in terms_orig:
        for x in path:
            verts_used.add(x)

    if verbose:
        print(f"  z_{v} (generator of H₁(T\\{v})):")
        print(f"    = {format_cycle(terms_orig)}")
        print(f"    vertices used: {sorted(verts_used)} ({len(verts_used)} of {n-1})")
        print(f"    #nonzero edges: {len(terms_orig)}")

    return {
        'gen_vec': gen,
        'terms': terms_orig,
        'vertices_used': verts_used,
        'allowed_1': allowed_1,
        'old_to_new': old_to_new,
        'new_to_old': new_to_old,
        'info': info,
    }


# ============================================================
# QUESTION 2: Filling chain c_v in T
# ============================================================

def find_filling_chain(A, n, v, z_data):
    """Since β₁(T)=0, z_v (viewed as a 1-cycle in T) is a boundary.
    Find c_v ∈ Ω₂(T) with ∂₂(c_v) = z_v."""

    # Get the full chain complex of T
    allowed_1_T = enumerate_allowed_paths(A, n, 1)
    allowed_2_T = enumerate_allowed_paths(A, n, 2)
    allowed_0_T = enumerate_allowed_paths(A, n, 0)

    omega_2_T = compute_omega_basis(A, n, 2, allowed_2_T, allowed_1_T)

    # z_v as a vector in A_1(T)
    idx_1_T = {path: i for i, path in enumerate(allowed_1_T)}
    z_in_T = np.zeros(len(allowed_1_T))

    for coeff, path in z_data['terms']:
        if path in idx_1_T:
            z_in_T[idx_1_T[path]] += coeff

    # ∂₂: Ω₂(T) → A₁(T)
    bd_2_T = build_full_boundary_matrix(allowed_2_T, allowed_1_T)
    dim_omega_2 = omega_2_T.shape[1] if omega_2_T.ndim == 2 else 0

    if dim_omega_2 == 0:
        return None

    bd_2_omega = bd_2_T @ omega_2_T  # maps Ω₂ coords → A₁ coords

    # Solve bd_2_omega @ x = z_in_T  (least squares)
    x, residuals, rank, sv = np.linalg.lstsq(bd_2_omega, z_in_T, rcond=None)

    # Check solution quality
    resid = np.linalg.norm(bd_2_omega @ x - z_in_T)
    if resid > 1e-6:
        print(f"  WARNING: filling chain residual = {resid:.6e}")
        return None

    # Convert to A₂ coords
    c_in_A2 = omega_2_T @ x

    # Which 2-paths (allowed) have nonzero coefficients?
    terms = []
    for i, coeff in enumerate(c_in_A2):
        if abs(coeff) > 1e-8:
            terms.append((round(coeff, 6), allowed_2_T[i]))

    # Does c_v use paths through v?
    uses_v = any(v in path for _, path in terms)
    verts_used = set()
    for _, path in terms:
        for x in path:
            verts_used.add(x)

    return {
        'terms': terms,
        'uses_v': uses_v,
        'verts_used': verts_used,
        'residual': resid,
    }


# ============================================================
# QUESTION 3: Relate z_a and z_b for bad pair {a,b}
# ============================================================

def compare_generators(A, n, bad_verts, z_data_map):
    """For pairs of bad vertices, compare their generators."""
    for a, b in combinations(bad_verts, 2):
        print(f"\n  Comparing z_{a} and z_{b}:")

        # Direction between a and b
        if A[a][b]:
            print(f"    {a}→{b} (a beats b)")
        else:
            print(f"    {b}→{a} (b beats a)")

        za = z_data_map[a]
        zb = z_data_map[b]

        # Edges in z_a that don't involve b (i.e., edges of T\{a,b})
        za_common = [(c, p) for c, p in za['terms'] if b not in p]
        zb_common = [(c, p) for c, p in zb['terms'] if a not in p]

        # Build edge-coefficient dicts
        za_dict = {p: c for c, p in za_common}
        zb_dict = {p: c for c, p in zb_common}

        all_edges = set(za_dict.keys()) | set(zb_dict.keys())

        # Compute z_a - z_b restricted to common edges
        diff = {}
        for e in all_edges:
            d = za_dict.get(e, 0) - zb_dict.get(e, 0)
            if abs(d) > 1e-8:
                diff[e] = d

        print(f"    z_{a} on T\\{{{a},{b}}}: {len(za_common)} edges")
        print(f"    z_{b} on T\\{{{a},{b}}}: {len(zb_common)} edges")
        print(f"    z_{a}-z_{b} restricted: {len(diff)} nonzero edges")

        if len(diff) == 0:
            print(f"    *** z_{a} = z_{b} on common edges! ***")
        else:
            # Check if this difference is a boundary in T\{a,b}
            C, m_ab, o2n = delete_vertex(A, n, a)
            # Now delete b from C
            b_new = o2n[b]
            D, m_ab2, o2n2 = delete_vertex(C, m_ab, b_new)

            # Map original vertices to T\{a,b} labels
            full_map = {}
            for orig, mid in o2n.items():
                if mid in o2n2:
                    full_map[orig] = o2n2[mid]

            # Check if diff is a boundary in T\{a,b}
            allowed_1_ab = enumerate_allowed_paths(D, m_ab2, 1)
            allowed_2_ab = enumerate_allowed_paths(D, m_ab2, 2)
            omega_2_ab = compute_omega_basis(D, m_ab2, 2, allowed_2_ab, allowed_1_ab)

            idx_1_ab = {path: i for i, path in enumerate(allowed_1_ab)}
            diff_vec = np.zeros(len(allowed_1_ab))

            for edge, coeff in diff.items():
                mapped = tuple(full_map.get(x, -1) for x in edge)
                if -1 not in mapped and mapped in idx_1_ab:
                    diff_vec[idx_1_ab[mapped]] += coeff

            if np.linalg.norm(diff_vec) < 1e-8:
                print(f"    diff maps to zero in T\\{{{a},{b}}} (edges lost)")
            else:
                # Check if it's a boundary
                bd_2_ab = build_full_boundary_matrix(allowed_2_ab, allowed_1_ab)
                dim_om2 = omega_2_ab.shape[1] if omega_2_ab.ndim == 2 else 0
                if dim_om2 > 0:
                    bd_2_om = bd_2_ab @ omega_2_ab
                    x, _, _, _ = np.linalg.lstsq(bd_2_om, diff_vec, rcond=None)
                    resid = np.linalg.norm(bd_2_om @ x - diff_vec)
                    if resid < 1e-6:
                        print(f"    z_{a}-z_{b} IS a boundary in T\\{{{a},{b}}} (residual {resid:.2e})")
                    else:
                        print(f"    z_{a}-z_{b} is NOT a boundary in T\\{{{a},{b}}} (residual {resid:.2e})")
                else:
                    print(f"    Ω₂(T\\{{{a},{b}}}) = 0, cannot be boundary")

            # Also check β₁ of T\{a,b}
            betti_ab = path_betti_numbers(D, m_ab2, max_dim=1)
            print(f"    β₁(T\\{{{a},{b}}}) = {betti_ab[1]}")


# ============================================================
# QUESTION 4: Why not 4 bad vertices?
# ============================================================

def check_four_bad(A, n, bad_verts, z_data_map):
    """If there were 4 bad vertices, look for linear dependence."""
    if len(bad_verts) < 4:
        return

    print(f"\n  4+ bad vertices found: {bad_verts}")

    # Restrict all z_v to common edges (edges not involving any bad vertex)
    common_edges = set()
    for v in bad_verts:
        zd = z_data_map[v]
        edges_v = set(p for _, p in zd['terms'] if all(x not in bad_verts or x == p[0] or x == p[1] for x in p))
        # Actually: edges not involving any OTHER bad vertex
        edges_v_common = set(p for _, p in zd['terms']
                           if not any(x in bad_verts and x != v for x in p)
                           and v not in p)
        if not common_edges:
            common_edges = edges_v_common

    # Get all edges used by any z_v that don't involve any bad vertex
    all_common = set()
    for v in bad_verts:
        for _, p in z_data_map[v]['terms']:
            if not any(x in bad_verts for x in p):
                all_common.add(p)

    edge_list = sorted(all_common)
    edge_idx = {e: i for i, e in enumerate(edge_list)}

    print(f"  Common edges (not involving any bad vertex): {len(edge_list)}")

    # Build matrix: rows = bad vertices, cols = common edges
    M = np.zeros((len(bad_verts), len(edge_list)))
    for i, v in enumerate(bad_verts):
        for coeff, p in z_data_map[v]['terms']:
            if p in edge_idx:
                M[i, edge_idx[p]] += coeff

    rank = np.linalg.matrix_rank(M, tol=1e-8)
    print(f"  Rank of restriction matrix ({len(bad_verts)} x {len(edge_list)}): {rank}")
    if rank < len(bad_verts):
        print(f"  *** LINEAR DEPENDENCE among restricted generators! ***")
    else:
        print(f"  Generators are linearly independent on common edges")


# ============================================================
# QUESTION 5: Dimension analysis
# ============================================================

def dimension_analysis(A, n, v):
    """Analyze dimensions of Z₁, B₁ for T and T\\v."""
    # T
    allowed_1_T = enumerate_allowed_paths(A, n, 1)
    allowed_2_T = enumerate_allowed_paths(A, n, 2)
    allowed_0_T = enumerate_allowed_paths(A, n, 0)

    omega_1_T = compute_omega_basis(A, n, 1, allowed_1_T, allowed_0_T)
    omega_2_T = compute_omega_basis(A, n, 2, allowed_2_T, allowed_1_T)

    dim_omega1_T = omega_1_T.shape[1] if omega_1_T.ndim == 2 else 0
    dim_omega2_T = omega_2_T.shape[1] if omega_2_T.ndim == 2 else 0

    # Z₁(T) = ker(∂₁ on Ω₁)
    bd_1_T = build_full_boundary_matrix(allowed_1_T, allowed_0_T)
    bd_1_om = bd_1_T @ omega_1_T if dim_omega1_T > 0 else np.zeros((0,0))
    if bd_1_om.size > 0:
        sv = np.linalg.svd(bd_1_om, compute_uv=False)
        rank_d1 = sum(s > 1e-8 for s in sv)
    else:
        rank_d1 = 0
    dim_Z1_T = dim_omega1_T - rank_d1

    # B₁(T) = im(∂₂ on Ω₂)
    bd_2_T = build_full_boundary_matrix(allowed_2_T, allowed_1_T)
    bd_2_om = bd_2_T @ omega_2_T if dim_omega2_T > 0 else np.zeros((0,0))
    if bd_2_om.size > 0:
        sv2 = np.linalg.svd(bd_2_om, compute_uv=False)
        dim_B1_T = sum(s > 1e-8 for s in sv2)
    else:
        dim_B1_T = 0

    # T\v
    B, m, o2n = delete_vertex(A, n, v)
    allowed_1_Tv = enumerate_allowed_paths(B, m, 1)
    allowed_2_Tv = enumerate_allowed_paths(B, m, 2)
    allowed_0_Tv = enumerate_allowed_paths(B, m, 0)

    omega_1_Tv = compute_omega_basis(B, m, 1, allowed_1_Tv, allowed_0_Tv)
    omega_2_Tv = compute_omega_basis(B, m, 2, allowed_2_Tv, allowed_1_Tv)

    dim_omega1_Tv = omega_1_Tv.shape[1] if omega_1_Tv.ndim == 2 else 0
    dim_omega2_Tv = omega_2_Tv.shape[1] if omega_2_Tv.ndim == 2 else 0

    bd_1_Tv = build_full_boundary_matrix(allowed_1_Tv, allowed_0_Tv)
    bd_1_om_v = bd_1_Tv @ omega_1_Tv if dim_omega1_Tv > 0 else np.zeros((0,0))
    if bd_1_om_v.size > 0:
        sv = np.linalg.svd(bd_1_om_v, compute_uv=False)
        rank_d1_v = sum(s > 1e-8 for s in sv)
    else:
        rank_d1_v = 0
    dim_Z1_Tv = dim_omega1_Tv - rank_d1_v

    bd_2_Tv = build_full_boundary_matrix(allowed_2_Tv, allowed_1_Tv)
    bd_2_om_v = bd_2_Tv @ omega_2_Tv if dim_omega2_Tv > 0 else np.zeros((0,0))
    if bd_2_om_v.size > 0:
        sv2 = np.linalg.svd(bd_2_om_v, compute_uv=False)
        dim_B1_Tv = sum(s > 1e-8 for s in sv2)
    else:
        dim_B1_Tv = 0

    return {
        'T': {'|A₁|': len(allowed_1_T), 'dim_Ω₁': dim_omega1_T,
              'dim_Z₁': dim_Z1_T, 'dim_B₁': dim_B1_T, 'β₁': dim_Z1_T - dim_B1_T,
              '|Ω₂|': dim_omega2_T},
        'T\\v': {'|A₁|': len(allowed_1_Tv), 'dim_Ω₁': dim_omega1_Tv,
                'dim_Z₁': dim_Z1_Tv, 'dim_B₁': dim_B1_Tv, 'β₁': dim_Z1_Tv - dim_B1_Tv,
                '|Ω₂|': dim_omega2_Tv},
    }


# ============================================================
# MAIN ANALYSIS
# ============================================================

def analyze_tournament(A, n, label="", verbose=True):
    """Full analysis of a single tournament."""
    betti_T = path_betti_numbers(A, n, max_dim=1)

    if betti_T[1] != 0:
        return None  # Only interested in β₁(T) = 0

    # Find bad vertices
    bad = []
    for v in range(n):
        B, m, _ = delete_vertex(A, n, v)
        betti_v = path_betti_numbers(B, m, max_dim=1)
        if betti_v[1] == 1:
            bad.append(v)

    num_bad = len(bad)
    if num_bad == 0:
        return {'num_bad': 0}

    if verbose:
        print(f"\n{'='*70}")
        print(f"Tournament {label}: n={n}, β₁(T)=0, bad vertices={bad}")
        print(f"{'='*70}")

        # Check transitivity
        trans = is_transitive_set(A, bad)
        print(f"Bad vertices transitive? {trans}")

        # Orientation among bad vertices
        for a, b in combinations(bad, 2):
            arrow = f"{a}→{b}" if A[a][b] else f"{b}→{a}"
            print(f"  {arrow}")

    # Q1: Compute generators z_v
    z_data_map = {}
    for v in bad:
        if verbose:
            print()
        zd = analyze_generator(A, n, v, verbose=verbose)
        if zd is not None:
            z_data_map[v] = zd

    # Q2: Find filling chains
    if verbose:
        print(f"\n--- Filling chains in T ---")
    filling_map = {}
    for v in bad:
        if v in z_data_map:
            fc = find_filling_chain(A, n, v, z_data_map[v])
            filling_map[v] = fc
            if verbose and fc is not None:
                print(f"  c_{v}: {len(fc['terms'])} terms, uses v={v}? {fc['uses_v']}, "
                      f"vertices: {sorted(fc['verts_used'])}")

    # Q3: Compare generators
    if verbose and len(z_data_map) >= 2:
        print(f"\n--- Generator comparisons ---")
        compare_generators(A, n, bad, z_data_map)

    # Q4: Check for 4+ bad vertices
    if len(bad) >= 4:
        check_four_bad(A, n, bad, z_data_map)

    # Q5: Dimension analysis
    if verbose:
        print(f"\n--- Dimension analysis ---")
        for v in bad[:3]:  # Limit output
            dims = dimension_analysis(A, n, v)
            print(f"  v={v}:")
            for space, d in dims.items():
                print(f"    {space}: {d}")

    return {
        'num_bad': num_bad,
        'bad': bad,
        'transitive': is_transitive_set(A, bad),
        'z_data': z_data_map,
        'filling': filling_map,
    }


def run_exhaustive(n, max_report=10):
    """Run exhaustive analysis for all tournaments of size n."""
    print(f"\n{'#'*70}")
    print(f"# EXHAUSTIVE ANALYSIS: n={n}")
    print(f"{'#'*70}")

    stats = Counter()
    bad_count_dist = Counter()
    transitive_count = 0
    non_transitive_count = 0
    total_beta0 = 0
    reported = 0

    gen_vertex_counts = []  # How many vertices does z_v use?
    filling_uses_v_count = 0
    filling_total = 0

    for idx, A in enumerate(all_tournaments(n)):
        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue
        total_beta0 += 1

        # Find bad vertices
        bad = []
        for v in range(n):
            B, m, _ = delete_vertex(A, n, v)
            bv = path_betti_numbers(B, m, max_dim=1)
            if bv[1] == 1:
                bad.append(v)

        num_bad = len(bad)
        bad_count_dist[num_bad] += 1

        if num_bad > 0:
            trans = is_transitive_set(A, bad)
            if trans:
                transitive_count += 1
            else:
                non_transitive_count += 1
                print(f"\n*** NON-TRANSITIVE BAD SET FOUND! Tournament #{idx}, bad={bad} ***")

        # Detailed analysis for tournaments with 3 bad vertices
        verbose = (num_bad == 3 and reported < max_report)
        if verbose:
            result = analyze_tournament(A, n, label=f"#{idx}", verbose=True)
            reported += 1

            # Collect stats
            if result and result.get('z_data'):
                for v, zd in result['z_data'].items():
                    gen_vertex_counts.append(len(zd['vertices_used']))
            if result and result.get('filling'):
                for v, fc in result['filling'].items():
                    if fc is not None:
                        filling_total += 1
                        if fc['uses_v']:
                            filling_uses_v_count += 1
        elif num_bad == 3:
            # Still collect stats silently
            result = analyze_tournament(A, n, label=f"#{idx}", verbose=False)
            if result and result.get('z_data'):
                for v, zd in result['z_data'].items():
                    gen_vertex_counts.append(len(zd['vertices_used']))
            if result and result.get('filling'):
                for v, fc in result['filling'].items():
                    if fc is not None:
                        filling_total += 1
                        if fc['uses_v']:
                            filling_uses_v_count += 1

    print(f"\n{'='*70}")
    print(f"SUMMARY for n={n}")
    print(f"{'='*70}")
    print(f"Total tournaments with β₁=0: {total_beta0}")
    print(f"Bad vertex count distribution: {dict(sorted(bad_count_dist.items()))}")
    print(f"Transitive bad sets: {transitive_count}")
    print(f"Non-transitive bad sets: {non_transitive_count}")

    if gen_vertex_counts:
        print(f"\nGenerator z_v vertex counts (for bad v in β₁=0 tournaments):")
        vc = Counter(gen_vertex_counts)
        for k in sorted(vc.keys()):
            print(f"  Uses {k} vertices: {vc[k]} times")

    if filling_total > 0:
        print(f"\nFilling chain c_v uses vertex v: {filling_uses_v_count}/{filling_total} "
              f"({100*filling_uses_v_count/filling_total:.1f}%)")

    return bad_count_dist


def run_sampled(n, num_samples=200, max_report=5):
    """Sample random tournaments and analyze those with β₁=0 and 3 bad vertices."""
    print(f"\n{'#'*70}")
    print(f"# SAMPLED ANALYSIS: n={n}, {num_samples} targets")
    print(f"{'#'*70}")

    stats = Counter()
    bad_count_dist = Counter()
    transitive_count = 0
    non_transitive_count = 0
    total_beta0 = 0
    reported = 0
    attempts = 0
    max_attempts = num_samples * 50

    gen_vertex_counts = []
    filling_uses_v_count = 0
    filling_total = 0

    while total_beta0 < num_samples and attempts < max_attempts:
        attempts += 1
        # Random tournament
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue
        total_beta0 += 1

        bad = []
        for v in range(n):
            B, m, _ = delete_vertex(A, n, v)
            bv = path_betti_numbers(B, m, max_dim=1)
            if bv[1] == 1:
                bad.append(v)

        num_bad = len(bad)
        bad_count_dist[num_bad] += 1

        if num_bad > 0:
            trans = is_transitive_set(A, bad)
            if trans:
                transitive_count += 1
            else:
                non_transitive_count += 1
                print(f"\n*** NON-TRANSITIVE BAD SET! bad={bad} ***")

        verbose = (num_bad == 3 and reported < max_report)
        if num_bad >= 1:
            result = analyze_tournament(A, n, label=f"sample-{total_beta0}", verbose=verbose)
            if verbose:
                reported += 1
            if result and result.get('z_data'):
                for v, zd in result['z_data'].items():
                    gen_vertex_counts.append(len(zd['vertices_used']))
            if result and result.get('filling'):
                for v, fc in result['filling'].items():
                    if fc is not None:
                        filling_total += 1
                        if fc['uses_v']:
                            filling_uses_v_count += 1

        if total_beta0 % 50 == 0:
            print(f"  ... {total_beta0}/{num_samples} β₁=0 tournaments found ({attempts} attempts)")

    print(f"\n{'='*70}")
    print(f"SUMMARY for n={n} (sampled)")
    print(f"{'='*70}")
    print(f"Total tournaments with β₁=0: {total_beta0} (from {attempts} random)")
    print(f"Bad vertex count distribution: {dict(sorted(bad_count_dist.items()))}")
    print(f"Transitive bad sets: {transitive_count}")
    print(f"Non-transitive bad sets: {non_transitive_count}")

    if gen_vertex_counts:
        print(f"\nGenerator z_v vertex counts:")
        vc = Counter(gen_vertex_counts)
        for k in sorted(vc.keys()):
            print(f"  Uses {k} vertices: {vc[k]} times")

    if filling_total > 0:
        print(f"\nFilling chain c_v uses vertex v: {filling_uses_v_count}/{filling_total} "
              f"({100*filling_uses_v_count/filling_total:.1f}%)")

    return bad_count_dist


# ============================================================
# FOCUSED ANALYSIS: The algebraic mechanism
# ============================================================

def algebraic_mechanism_analysis(n, num_examples=5):
    """Deep dive into the algebraic structure for a few examples."""
    print(f"\n{'#'*70}")
    print(f"# ALGEBRAIC MECHANISM DEEP DIVE: n={n}")
    print(f"{'#'*70}")

    count = 0
    for A in all_tournaments(n) if n <= 6 else []:
        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue

        bad = []
        for v in range(n):
            B, m, _ = delete_vertex(A, n, v)
            bv = path_betti_numbers(B, m, max_dim=1)
            if bv[1] == 1:
                bad.append(v)

        if len(bad) != 3:
            continue

        count += 1
        if count > num_examples:
            break

        print(f"\n{'='*60}")
        print(f"Example {count}: bad = {bad}")
        print(f"{'='*60}")

        # Show full adjacency restricted to bad vertices
        print(f"Adjacency among bad vertices:")
        for a in bad:
            for b in bad:
                if a != b:
                    print(f"  {a}→{b}" if A[a][b] else "", end="")
            print()

        # Detailed z_v analysis
        z_data_map = {}
        for v in bad:
            zd = analyze_generator(A, n, v, verbose=True)
            if zd:
                z_data_map[v] = zd

        if len(z_data_map) < 3:
            continue

        # KEY ANALYSIS: Look at the structure of each z_v
        print(f"\n--- Structural comparison ---")

        # For each z_v, decompose into: edges involving other bad vertices vs not
        for v in bad:
            zd = z_data_map[v]
            other_bad = [b for b in bad if b != v]

            edges_no_bad = [(c, p) for c, p in zd['terms'] if not any(x in other_bad for x in p)]
            edges_with_bad = [(c, p) for c, p in zd['terms'] if any(x in other_bad for x in p)]

            print(f"\n  z_{v} decomposition:")
            print(f"    Edges not involving other bad: {len(edges_no_bad)}")
            for c, p in sorted(edges_with_bad, key=lambda x: x[1]):
                print(f"    Edge involving bad: {c:+.0f}*({p[0]}→{p[1]})")

        # CRITICAL: Check z_a + z_b + z_c on the common edges
        print(f"\n--- Sum z_a + z_b + z_c on common edges ---")
        a, b, c = bad

        # Common edges = edges not involving any bad vertex
        all_common_edges = set()
        for v in bad:
            for _, p in z_data_map[v]['terms']:
                if not any(x in bad for x in p):
                    all_common_edges.add(p)

        sum_on_common = {}
        for e in all_common_edges:
            s = 0
            for v in bad:
                for coeff, p in z_data_map[v]['terms']:
                    if p == e:
                        s += coeff
            if abs(s) > 1e-8:
                sum_on_common[e] = round(s, 4)

        print(f"  Common edges: {len(all_common_edges)}")
        print(f"  Nonzero in sum: {len(sum_on_common)}")
        if sum_on_common:
            for e, s in sorted(sum_on_common.items()):
                print(f"    {s:+}*({e[0]}→{e[1]})")
        else:
            print(f"  *** SUM IS ZERO on all common edges! ***")

        # Pairwise differences on common edges
        print(f"\n--- Pairwise analysis ---")
        compare_generators(A, n, bad, z_data_map)

        # Filling chains
        print(f"\n--- Filling chains ---")
        for v in bad:
            fc = find_filling_chain(A, n, v, z_data_map[v])
            if fc:
                # Count terms through each bad vertex
                for b2 in bad:
                    if b2 == v:
                        continue
                    terms_thru = [(c, p) for c, p in fc['terms'] if b2 in p]
                    print(f"  c_{v}: {len(terms_thru)} terms through bad vertex {b2}")
                terms_thru_v = [(c, p) for c, p in fc['terms'] if v in p]
                print(f"  c_{v}: {len(terms_thru_v)} terms through v={v} itself")


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    random.seed(42)
    np.set_printoptions(precision=4, suppress=True, linewidth=120)

    print("β₁ TRANSITIVE MECHANISM INVESTIGATION")
    print("=" * 70)
    print("Goal: Understand why bad vertices form transitive triples")
    print("      and prove Σ_v β₁(T\\v) ≤ 3 when β₁(T)=0")
    print()

    # n=5: exhaustive
    dist_5 = run_exhaustive(5, max_report=5)

    # n=5: algebraic deep dive
    algebraic_mechanism_analysis(5, num_examples=3)

    # n=6: exhaustive (feasible: 2^15 = 32768 tournaments)
    dist_6 = run_exhaustive(6, max_report=3)

    # n=6: algebraic deep dive
    algebraic_mechanism_analysis(6, num_examples=3)

    # n=7: sampled
    dist_7 = run_sampled(7, num_samples=100, max_report=3)

    print(f"\n\n{'#'*70}")
    print(f"# FINAL SUMMARY")
    print(f"{'#'*70}")
    print(f"\nn=5 bad count distribution: {dict(sorted(dist_5.items()))}")
    print(f"n=6 bad count distribution: {dict(sorted(dist_6.items()))}")
    print(f"n=7 bad count distribution: {dict(sorted(dist_7.items()))}")

    print("\nDone.")
