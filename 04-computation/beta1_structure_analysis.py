#!/usr/bin/env python3
"""
beta1_structure_analysis.py — Deep structural analysis of beta_1 under vertex deletion

Goal: Prove that for any tournament T with beta_1(T)=0, sum_v beta_1(T\v) <= 3.

Analysis:
1. Explicit boundary matrix d2 and its image B1
2. Vertex deletion: rank drops and the "v-boundary space"
3. Triple intersections of vertex boundary spaces
4. Double-deletion B1(T\{a,b}) dimensions
5. Matroid perspective

opus-2026-03-08
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
import random
import sys

# ============================================================
# CORE: GLMY path homology computation (from path_homology_v2.py)
# ============================================================

def enumerate_allowed_paths(A, n, p):
    """All sequences of p+1 distinct vertices following directed edges."""
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

def compute_beta1_full(A, n):
    """Compute beta_1 and return all intermediate data."""
    allowed_0 = enumerate_allowed_paths(A, n, 0)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_3 = enumerate_allowed_paths(A, n, 3)

    omega_1 = compute_omega_basis(A, n, 1, allowed_1, allowed_0)
    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)

    # d1: Omega_1 -> A_0
    bd_1 = build_full_boundary_matrix(allowed_1, allowed_0)
    bd_1_omega = bd_1 @ omega_1
    S1 = np.linalg.svd(bd_1_omega, compute_uv=False)
    rank_d1 = sum(s > 1e-8 for s in S1)
    dim_omega1 = omega_1.shape[1]
    ker_d1 = dim_omega1 - rank_d1  # = dim(Z_1)

    # d2: Omega_2 -> A_1
    bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
    dim_omega2 = omega_2.shape[1] if omega_2.ndim == 2 else 0
    if dim_omega2 > 0:
        bd_2_omega = bd_2 @ omega_2
        S2 = np.linalg.svd(bd_2_omega, compute_uv=False)
        rank_d2 = sum(s > 1e-8 for s in S2)
    else:
        bd_2_omega = np.zeros((len(allowed_1), 0))
        rank_d2 = 0

    beta_1 = ker_d1 - rank_d2

    return {
        'beta_1': max(0, beta_1),
        'allowed_1': allowed_1,
        'allowed_2': allowed_2,
        'omega_1': omega_1,
        'omega_2': omega_2,
        'bd_2': bd_2,
        'bd_2_omega': bd_2_omega,
        'dim_omega1': dim_omega1,
        'dim_omega2': dim_omega2,
        'dim_Z1': ker_d1,
        'dim_B1': rank_d2,
        'rank_d1': rank_d1,
    }

def compute_beta1_only(A, n):
    """Just compute beta_1, no extras."""
    info = compute_beta1_full(A, n)
    return info['beta_1']

def subtournament(A, n, remove_vertices):
    """Return adjacency matrix with vertices removed."""
    remove_set = set(remove_vertices)
    keep = [v for v in range(n) if v not in remove_set]
    n2 = len(keep)
    B = [[0]*n2 for _ in range(n2)]
    for i, vi in enumerate(keep):
        for j, vj in enumerate(keep):
            B[i][j] = A[vi][vj]
    return B, n2

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

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def count_transitive_triples(A, n):
    """Count transitive triples: ordered (a,b,c) with a->b->c and a->c."""
    tt = 0
    for a in range(n):
        for b in range(n):
            if a == b or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    tt += 1
    return tt

def rank_of_matrix(M, tol=1e-8):
    if M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    S = np.linalg.svd(M, compute_uv=False)
    return int(sum(s > tol for s in S))

# ============================================================
# PART 1: Understand B1 structure
# ============================================================

def analyze_B1_structure(A, n, label=""):
    """Detailed analysis of the boundary space B1 = im(d2)."""
    info = compute_beta1_full(A, n)
    t3 = count_3cycles(A, n)
    tt = count_transitive_triples(A, n)

    print(f"\n  {label}")
    print(f"    n={n}, t3={t3}, TT(ordered)={tt}")
    print(f"    |A_1| = |edges| = {len(info['allowed_1'])}")
    print(f"    |A_2| = |allowed 2-paths| = {len(info['allowed_2'])}")
    print(f"    dim(Omega_1) = {info['dim_omega1']}")
    print(f"    dim(Omega_2) = {info['dim_omega2']}")
    print(f"    dim(Z_1) = {info['dim_Z1']}")
    print(f"    dim(B_1) = {info['dim_B1']}")
    print(f"    beta_1 = {info['beta_1']}")
    print(f"    C(n-1,2) = {(n-1)*(n-2)//2}")

    # For tournaments: Omega_1 = A_1 = all n(n-1) directed edges
    # Z_1 = flow-balanced edge combos
    # B_1 = im(d2|_{Omega_2})

    # Check: for tournaments, is Omega_2 = A_2?
    if info['dim_omega2'] == len(info['allowed_2']):
        print(f"    Omega_2 = A_2 (all 2-paths are allowed)")
    else:
        print(f"    Omega_2 != A_2 ({info['dim_omega2']} vs {len(info['allowed_2'])})")

    return info

# ============================================================
# PART 2: Vertex deletion analysis
# ============================================================

def analyze_deletion(A, n, v, info_T=None):
    """Analyze what happens when vertex v is deleted.

    Returns dict with:
    - beta1_Tv: beta_1(T\v)
    - rank_drop: dim(B1(T)) - dim(B1(T\v) restricted to edges of T\v)
    - dim_Vv: dimension of v-boundary space
    """
    if info_T is None:
        info_T = compute_beta1_full(A, n)

    # Compute T\v
    Av, nv = subtournament(A, n, [v])
    info_Tv = compute_beta1_full(Av, nv)

    # Identify which 2-paths in Omega_2(T) involve vertex v
    # A 2-path (a,b,c) involves v if v in {a,b,c}
    involved_indices = []
    not_involved_indices = []
    for j, path in enumerate(info_T['allowed_2']):
        if v in path:
            involved_indices.append(j)
        else:
            not_involved_indices.append(j)

    # The boundary matrix bd_2 maps A_2 -> A_1
    # V_v = image of {d2(p) : p in Omega_2(T), v in p}
    # We need to work in Omega_2 coordinates

    # Get the image of v-paths through d2
    bd_2_omega = info_T['bd_2_omega']  # shape (|A_1|, dim_Omega_2)

    # But we need to identify which Omega_2 basis vectors correspond to v-paths
    # omega_2 is a basis of Omega_2 as columns in A_2 coordinates
    # The columns of omega_2 are linear combinations of A_2 paths
    # We can't simply split by "contains v" at the Omega_2 level

    # Instead, let's work directly in A_2:
    # Build the boundary image from v-paths in A_2
    bd_2_full = info_T['bd_2']  # maps A_2 -> A_1

    # Columns of bd_2_full corresponding to v-paths
    if involved_indices:
        V_v_matrix = bd_2_full[:, involved_indices]
        dim_Vv = rank_of_matrix(V_v_matrix)
    else:
        V_v_matrix = np.zeros((bd_2_full.shape[0], 0))
        dim_Vv = 0

    # Columns not involving v
    if not_involved_indices:
        no_v_matrix = bd_2_full[:, not_involved_indices]
        dim_no_v = rank_of_matrix(no_v_matrix)
    else:
        no_v_matrix = np.zeros((bd_2_full.shape[0], 0))
        dim_no_v = 0

    # But wait - we should restrict to Omega_2, not all of A_2
    # For tournaments, Omega_2 = A_2 typically. Let's check and handle both cases.

    # Actually for the intersection analysis we need more care.
    # Let's use a different approach: work with the FULL d2 matrix
    # and project onto the edges of T\v.

    # Edges of T\v: those not involving v
    edge_idx_Tv = {path: i for i, path in enumerate(info_T['allowed_1']) if v not in path}
    edge_indices_Tv = [i for i, path in enumerate(info_T['allowed_1']) if v not in path]

    # B1(T) restricted to edges of T\v
    full_image = bd_2_omega  # image of d2 in A_1(T) coordinates
    if full_image.shape[1] > 0:
        # Project to T\v edges
        B1_T_restricted = full_image[edge_indices_Tv, :]
        rank_B1_restricted = rank_of_matrix(B1_T_restricted)
    else:
        rank_B1_restricted = 0

    return {
        'beta1_Tv': info_Tv['beta_1'],
        'dim_B1_Tv': info_Tv['dim_B1'],
        'dim_Z1_Tv': info_Tv['dim_Z1'],
        'dim_Vv': dim_Vv,
        'dim_no_v': dim_no_v,
        'n_vpaths': len(involved_indices),
        'n_novpaths': len(not_involved_indices),
        'rank_B1_restricted': rank_B1_restricted,
        'info_Tv': info_Tv,
    }

# ============================================================
# PART 3: Triple intersection analysis
# ============================================================

def analyze_triple_intersection(A, n, vertices, info_T=None):
    """Analyze V_a, V_b, V_c boundary spaces and their intersections."""
    if info_T is None:
        info_T = compute_beta1_full(A, n)

    bd_2_full = info_T['bd_2']
    allowed_2 = info_T['allowed_2']

    # For each vertex, get columns of bd_2 from paths involving that vertex
    V_matrices = {}
    for v in vertices:
        cols = [j for j, path in enumerate(allowed_2) if v in path]
        if cols:
            V_matrices[v] = bd_2_full[:, cols]
        else:
            V_matrices[v] = np.zeros((bd_2_full.shape[0], 0))

    # Dimensions
    dims = {v: rank_of_matrix(V_matrices[v]) for v in vertices}

    # Pairwise: V_a + V_b (column-concatenate, take rank)
    pair_dims = {}
    for a, b in combinations(vertices, 2):
        combined = np.hstack([V_matrices[a], V_matrices[b]])
        pair_dims[(a,b)] = rank_of_matrix(combined)

    # Triple: V_a + V_b + V_c
    if len(vertices) >= 3:
        triple_combined = np.hstack([V_matrices[v] for v in vertices[:3]])
        triple_dim = rank_of_matrix(triple_combined)
    else:
        triple_dim = None

    # Intersection dimensions via inclusion-exclusion:
    # dim(V_a cap V_b) = dim(V_a) + dim(V_b) - dim(V_a + V_b)
    intersect_dims = {}
    for a, b in combinations(vertices, 2):
        intersect_dims[(a,b)] = dims[a] + dims[b] - pair_dims[(a,b)]

    return {
        'dims': dims,
        'pair_union_dims': pair_dims,
        'intersect_dims': intersect_dims,
        'triple_union_dim': triple_dim,
    }

# ============================================================
# PART 4: Double deletion analysis
# ============================================================

def analyze_double_deletion(A, n, v1, v2):
    """Compute B1(T\{v1,v2})."""
    Avv, nvv = subtournament(A, n, [v1, v2])
    info = compute_beta1_full(Avv, nvv)
    return {
        'beta1': info['beta_1'],
        'dim_B1': info['dim_B1'],
        'dim_Z1': info['dim_Z1'],
    }

# ============================================================
# MAIN ANALYSIS
# ============================================================

def run_analysis():
    print("=" * 70)
    print("BETA_1 STRUCTURE ANALYSIS — Deletion Bound Investigation")
    print("=" * 70)

    # ============================
    # PART 1: B1 structure at n=4,5
    # ============================
    print("\n" + "=" * 70)
    print("PART 1: B1 STRUCTURE")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")
        print(f"  Total edges = n(n-1) = {n*(n-1)}")
        print(f"  C(n-1,2) = {(n-1)*(n-2)//2}  [expected dim(Z1) for tournament]")

        beta_dist = Counter()
        b1_data = []

        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            t3 = count_3cycles(A, n)
            beta_dist[info['beta_1']] += 1
            b1_data.append({
                't3': t3,
                'beta1': info['beta_1'],
                'dim_omega2': info['dim_omega2'],
                'dim_B1': info['dim_B1'],
                'dim_Z1': info['dim_Z1'],
                'n_2paths': len(info['allowed_2']),
            })

        print(f"\n  beta_1 distribution: {dict(beta_dist)}")

        # Group by t3 and show B1 dimensions
        by_t3 = defaultdict(list)
        for d in b1_data:
            by_t3[d['t3']].append(d)

        print(f"\n  {'t3':>4} | {'#T':>4} | {'b1=0':>5} | {'b1=1':>5} | {'dim_Z1':>7} | {'dim_B1 range':>14} | {'|A_2| range':>14} | {'dim_Om2 range':>14}")
        print(f"  {'-'*4}-+-{'-'*4}-+-{'-'*5}-+-{'-'*5}-+-{'-'*7}-+-{'-'*14}-+-{'-'*14}-+-{'-'*14}")
        for t3 in sorted(by_t3.keys()):
            items = by_t3[t3]
            n_b0 = sum(1 for d in items if d['beta1'] == 0)
            n_b1 = sum(1 for d in items if d['beta1'] == 1)
            z1s = [d['dim_Z1'] for d in items]
            b1s = [d['dim_B1'] for d in items]
            a2s = [d['n_2paths'] for d in items]
            om2s = [d['dim_omega2'] for d in items]
            print(f"  {t3:>4} | {len(items):>4} | {n_b0:>5} | {n_b1:>5} | {z1s[0]:>7} | [{min(b1s):>5},{max(b1s):>5}] | [{min(a2s):>5},{max(a2s):>5}] | [{min(om2s):>5},{max(om2s):>5}]")

    # ============================
    # PART 2: Deletion analysis at n=5
    # ============================
    print("\n" + "=" * 70)
    print("PART 2: VERTEX DELETION ANALYSIS (n=5 exhaustive)")
    print("=" * 70)

    n = 5
    bad_vertex_counts = Counter()  # how many "bad" vertices (beta1(T\v) > 0) per tournament
    sum_beta1_dist = Counter()  # distribution of sum_v beta1(T\v)
    max_sum_beta1 = 0
    worst_cases = []

    all_deletion_data = []

    for idx, A in enumerate(all_tournaments(n)):
        info_T = compute_beta1_full(A, n)
        if info_T['beta_1'] != 0:
            continue  # only consider T with beta_1(T) = 0

        deletion_betas = []
        for v in range(n):
            Av, nv = subtournament(A, n, [v])
            b1v = compute_beta1_only(Av, nv)
            deletion_betas.append(b1v)

        n_bad = sum(1 for b in deletion_betas if b > 0)
        sum_b = sum(deletion_betas)
        bad_vertex_counts[n_bad] += 1
        sum_beta1_dist[sum_b] += 1

        if sum_b > max_sum_beta1:
            max_sum_beta1 = sum_b

        if sum_b >= 2:
            worst_cases.append((A, deletion_betas, idx))

        all_deletion_data.append({
            'A': A,
            'idx': idx,
            'deletion_betas': deletion_betas,
            'n_bad': n_bad,
            'sum_b': sum_b,
            't3': count_3cycles(A, n),
        })

    print(f"\n  Tournaments with beta_1(T)=0: {len(all_deletion_data)}")
    print(f"\n  # bad vertices distribution: {dict(sorted(bad_vertex_counts.items()))}")
    print(f"  sum beta_1(T\\v) distribution: {dict(sorted(sum_beta1_dist.items()))}")
    print(f"  MAX sum beta_1(T\\v) = {max_sum_beta1}")

    # Detail on worst cases
    if worst_cases:
        print(f"\n  Worst cases (sum >= 2):")
        for A, betas, idx in worst_cases[:10]:
            t3 = count_3cycles(A, n)
            bad_verts = [v for v in range(n) if betas[v] > 0]
            print(f"    idx={idx}: t3={t3}, betas={betas}, bad_vertices={bad_verts}")

    # ============================
    # DETAILED ANALYSIS of worst cases at n=5
    # ============================
    print("\n" + "=" * 70)
    print("PART 2b: DETAILED DELETION ANALYSIS (worst cases at n=5)")
    print("=" * 70)

    for A, betas, idx in worst_cases[:5]:
        t3 = count_3cycles(A, n)
        print(f"\n  === Tournament idx={idx}, t3={t3}, deletion betas={betas} ===")
        info_T = compute_beta1_full(A, n)
        print(f"    dim(B1(T)) = {info_T['dim_B1']}, dim(Z1(T)) = {info_T['dim_Z1']}")

        bad_verts = [v for v in range(n) if betas[v] > 0]

        for v in range(n):
            del_info = analyze_deletion(A, n, v, info_T)
            marker = " ***BAD***" if betas[v] > 0 else ""
            print(f"    v={v}: beta1(T\\v)={del_info['beta1_Tv']}, "
                  f"dim_B1(T\\v)={del_info['dim_B1_Tv']}, dim_Z1(T\\v)={del_info['dim_Z1_Tv']}, "
                  f"#v-paths={del_info['n_vpaths']}, dim(V_v)={del_info['dim_Vv']}, "
                  f"dim(no-v)={del_info['dim_no_v']}{marker}")

        # Triple intersection analysis (if >= 3 bad vertices)
        if len(bad_verts) >= 2:
            trip = analyze_triple_intersection(A, n, bad_verts, info_T)
            print(f"\n    Boundary space dimensions: {trip['dims']}")
            print(f"    Pairwise union dims: {trip['pair_union_dims']}")
            print(f"    Pairwise intersection dims: {trip['intersect_dims']}")
            if trip['triple_union_dim'] is not None:
                print(f"    Triple union dim: {trip['triple_union_dim']}")

        # Double deletion
        if len(bad_verts) >= 2:
            print(f"\n    Double deletions:")
            for a, b in combinations(bad_verts, 2):
                dd = analyze_double_deletion(A, n, a, b)
                print(f"      T\\{{{a},{b}}}: beta1={dd['beta1']}, dim_B1={dd['dim_B1']}, dim_Z1={dd['dim_Z1']}")

    # ============================
    # PART 3: n=6 analysis
    # ============================
    print("\n" + "=" * 70)
    print("PART 3: n=6 EXHAUSTIVE (beta_1=0 tournaments)")
    print("=" * 70)

    n = 6
    bad_vertex_counts_6 = Counter()
    sum_beta1_dist_6 = Counter()
    max_sum_6 = 0
    worst_6 = []
    count_b0 = 0
    count_total = 0

    for A in all_tournaments(n):
        count_total += 1
        info_T = compute_beta1_full(A, n)
        if info_T['beta_1'] != 0:
            continue
        count_b0 += 1

        deletion_betas = []
        for v in range(n):
            Av, nv = subtournament(A, n, [v])
            b1v = compute_beta1_only(Av, nv)
            deletion_betas.append(b1v)

        n_bad = sum(1 for b in deletion_betas if b > 0)
        sum_b = sum(deletion_betas)
        bad_vertex_counts_6[n_bad] += 1
        sum_beta1_dist_6[sum_b] += 1

        if sum_b > max_sum_6:
            max_sum_6 = sum_b

        if sum_b >= 3:
            worst_6.append((A, deletion_betas, count_total))

        if count_b0 % 1000 == 0:
            print(f"    ... processed {count_b0} tournaments with beta_1=0 (of {count_total} total)", file=sys.stderr)

    print(f"\n  Total tournaments: {count_total}")
    print(f"  With beta_1(T)=0: {count_b0}")
    print(f"\n  # bad vertices distribution: {dict(sorted(bad_vertex_counts_6.items()))}")
    print(f"  sum beta_1(T\\v) distribution: {dict(sorted(sum_beta1_dist_6.items()))}")
    print(f"  MAX sum beta_1(T\\v) = {max_sum_6}")

    if worst_6:
        print(f"\n  Worst cases (sum >= 3) — showing up to 10:")
        for A, betas, idx in worst_6[:10]:
            t3 = count_3cycles(A, n)
            bad_verts = [v for v in range(n) if betas[v] > 0]
            print(f"    idx={idx}: t3={t3}, betas={betas}, bad_vertices={bad_verts}")

    # Detailed analysis of extreme n=6 cases
    if worst_6:
        print(f"\n  --- Detailed analysis of worst n=6 cases ---")
        for A, betas, idx in worst_6[:3]:
            t3 = count_3cycles(A, n)
            print(f"\n  === Tournament idx={idx}, t3={t3}, deletion betas={betas} ===")
            info_T = compute_beta1_full(A, n)
            print(f"    dim(B1(T)) = {info_T['dim_B1']}, dim(Z1(T)) = {info_T['dim_Z1']}")

            bad_verts = [v for v in range(n) if betas[v] > 0]

            for v in range(n):
                del_info = analyze_deletion(A, n, v, info_T)
                marker = " ***BAD***" if betas[v] > 0 else ""
                print(f"    v={v}: beta1(T\\v)={del_info['beta1_Tv']}, "
                      f"dim_B1(T\\v)={del_info['dim_B1_Tv']}, dim_Z1(T\\v)={del_info['dim_Z1_Tv']}, "
                      f"#v-paths={del_info['n_vpaths']}, dim(V_v)={del_info['dim_Vv']}{marker}")

            if len(bad_verts) >= 2:
                trip = analyze_triple_intersection(A, n, bad_verts[:4], info_T)
                print(f"\n    Boundary space dims: {trip['dims']}")
                print(f"    Pairwise union dims: {trip['pair_union_dims']}")
                print(f"    Pairwise intersect dims: {trip['intersect_dims']}")
                if trip['triple_union_dim'] is not None:
                    print(f"    Triple union dim: {trip['triple_union_dim']}")

            if len(bad_verts) >= 2:
                print(f"\n    Double deletions:")
                for a, b in combinations(bad_verts[:4], 2):
                    dd = analyze_double_deletion(A, n, a, b)
                    print(f"      T\\{{{a},{b}}}: beta1={dd['beta1']}, dim_B1={dd['dim_B1']}, dim_Z1={dd['dim_Z1']}")

    # ============================
    # PART 4: KEY RANK INEQUALITY ANALYSIS
    # ============================
    print("\n" + "=" * 70)
    print("PART 4: RANK INEQUALITY ANALYSIS")
    print("=" * 70)

    print("\n  For a tournament T on n vertices with beta_1(T)=0:")
    print("  dim(Z_1) = dim(B_1) = rank(d2|_{Omega_2})")
    print("  For T\\v: dim(Z_1(T\\v)) = C(n-2,2), dim(B_1(T\\v)) = C(n-2,2) - beta_1(T\\v)")
    print()

    # Verify the rank-drop formula
    print("  Verifying rank-drop formula at n=5:")
    n5 = 5
    for A, betas, idx in worst_cases[:3]:
        info_T = compute_beta1_full(A, n5)
        print(f"\n    Tournament idx={idx}, betas={betas}")
        print(f"    dim(B1(T)) = {info_T['dim_B1']}")

        for v in range(5):
            Av, nv = subtournament(A, 5, [v])
            info_Tv = compute_beta1_full(Av, nv)
            rank_drop = info_T['dim_B1'] - info_Tv['dim_B1']
            expected_drop = info_Tv['dim_Z1'] - info_Tv['dim_B1']  # = beta_1(T\v)... no
            # Actually:
            # dim(B1(T)) = dim(Z1(T)) = C(n-1,2) = C(4,2) = 6 (since beta_1(T)=0)
            # dim(Z1(T\v)) = C(n-2,2) = C(3,2) = 3
            # dim(B1(T\v)) = dim(Z1(T\v)) - beta_1(T\v) = 3 - beta_1(T\v)
            # rank_drop = dim(B1(T)) - dim(B1(T\v)) = 6 - (3 - beta_1(T\v)) = 3 + beta_1(T\v)
            print(f"    v={v}: dim_B1(T\\v)={info_Tv['dim_B1']}, dim_Z1(T\\v)={info_Tv['dim_Z1']}, "
                  f"rank_drop={rank_drop}, "
                  f"expected {info_T['dim_B1'] - info_Tv['dim_Z1'] + info_Tv['beta_1']}")

    # ============================
    # PART 5: Z1 dimension analysis
    # ============================
    print("\n" + "=" * 70)
    print("PART 5: Z_1 DIMENSION FOR TOURNAMENTS")
    print("=" * 70)

    print("\n  Key question: Is dim(Z_1) always = C(n-1,2) for tournaments?")
    for n in [3, 4, 5]:
        print(f"\n  n={n}: C(n-1,2) = {(n-1)*(n-2)//2}")
        z1_values = set()
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            z1_values.add(info['dim_Z1'])
        print(f"    dim(Z_1) values observed: {sorted(z1_values)}")

    # n=6 sample
    print(f"\n  n=6: C(5,2) = 10")
    z1_vals_6 = set()
    cnt = 0
    for A in all_tournaments(6):
        info = compute_beta1_full(A, 6)
        z1_vals_6.add(info['dim_Z1'])
        cnt += 1
        if cnt >= 200:
            break
    print(f"    dim(Z_1) values observed (first 200): {sorted(z1_vals_6)}")

    # ============================
    # PART 6: The Omega_2 = A_2 question
    # ============================
    print("\n" + "=" * 70)
    print("PART 6: IS Omega_2 = A_2 FOR TOURNAMENTS?")
    print("=" * 70)

    for n in [3, 4, 5]:
        all_equal = True
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            if info['dim_omega2'] != len(info['allowed_2']):
                all_equal = False
                break
        print(f"  n={n}: Omega_2 = A_2 for all tournaments? {all_equal}")

    # ============================
    # PART 7: Sum constraint derivation attempt
    # ============================
    print("\n" + "=" * 70)
    print("PART 7: SUM CONSTRAINT ANALYSIS")
    print("=" * 70)

    print("""
  KEY INSIGHT ATTEMPT:

  For tournament T on n vertices with beta_1(T) = 0:
  - dim(Z_1(T)) = dim(B_1(T)) = some value D
  - For T\\v: dim(Z_1(T\\v)) = D' (should be constant across all T\\v if all are tournaments)
  - beta_1(T\\v) = dim(Z_1(T\\v)) - dim(B_1(T\\v))

  The 2-paths of T decompose: each 2-path (a,b,c) involves exactly 3 vertices.
  The boundaries d2(a,b,c) span B_1(T).

  When we delete v, we lose:
  1. All 2-paths involving v → their boundary contributions are lost from B_1
  2. All edges involving v → these are removed from the ambient space
  3. The cycle space Z_1 restricted to T\\v has lower dimension

  The "price" of deleting v:
  rank_drop(v) = dim(B_1(T)) - dim(B_1(T\\v) ∩ B_1(T))

  But B_1(T\\v) is NOT a subspace of B_1(T) in general, because T\\v may have
  2-paths that are NOT 2-paths of T (new "allowed" paths could emerge from removing v).

  Wait — removing v can only REDUCE the set of 2-paths (a,b,c) since we lose all
  paths through v, and paths not through v remain unchanged. So A_2(T\\v) ⊂ A_2(T)
  (after relabeling). Similarly Omega_2(T\\v) ⊂ Omega_2(T).

  Therefore: B_1(T\\v) ⊂ B_1(T) (as subspaces of the edges of T\\v).

  And: beta_1(T\\v) = dim(Z_1(T\\v)) - dim(B_1(T\\v))

  Since B_1(T\\v) ⊂ B_1(T)|_{edges of T\\v}, and B_1(T)|_{edges of T\\v} ⊂ Z_1(T\\v)
  (because cycles of T restricted to T\\v edges are still cycles),
  we get: dim(B_1(T\\v)) ≤ dim(B_1(T)|_{edges of T\\v}) ≤ dim(Z_1(T\\v))
  """)

    # Verify the containment B_1(T\v) subset B_1(T)
    print("  Verifying B_1(T\\v) ⊂ B_1(T) at n=5:")
    for A, betas, idx in worst_cases[:3]:
        info_T = compute_beta1_full(A, 5)
        for v in range(5):
            Av, nv = subtournament(A, 5, [v])
            info_Tv = compute_beta1_full(Av, nv)

            # Get B_1(T\v) generators in T's edge coordinates
            # T\v edges map to T's edges (removing v-edges)
            keep = [u for u in range(5) if u != v]
            # Map: edge (i,j) in T\v (0-indexed) corresponds to (keep[i], keep[j]) in T
            edge_map = {}
            for i, path in enumerate(info_Tv['allowed_1']):
                # path = (a,b) in T\v's labeling, corresponds to (keep[a], keep[b]) in T
                orig_edge = (keep[path[0]], keep[path[1]])
                for j, te in enumerate(info_T['allowed_1']):
                    if te == orig_edge:
                        edge_map[i] = j
                        break

            # Embed B_1(T\v) into T's edge space
            bd_Tv = info_Tv['bd_2_omega']
            if bd_Tv.shape[1] == 0:
                continue

            embedded = np.zeros((len(info_T['allowed_1']), bd_Tv.shape[1]))
            for i in range(bd_Tv.shape[0]):
                if i in edge_map:
                    embedded[edge_map[i], :] = bd_Tv[i, :]

            # Check if columns of embedded are in B_1(T) = column space of info_T['bd_2_omega']
            B1_T = info_T['bd_2_omega']
            if B1_T.shape[1] == 0:
                continue

            combined = np.hstack([B1_T, embedded])
            rank_combined = rank_of_matrix(combined)
            rank_B1T = rank_of_matrix(B1_T)

            contained = (rank_combined == rank_B1T)
            if not contained:
                print(f"    VIOLATION! idx={idx}, v={v}: B_1(T\\v) NOT in B_1(T)! "
                      f"rank_B1T={rank_B1T}, rank_combined={rank_combined}")

    print("    (If no violations printed, containment holds for all checked cases)")

    # ============================
    # PART 8: Counting argument
    # ============================
    print("\n" + "=" * 70)
    print("PART 8: COUNTING ARGUMENT — 2-PATH OVERLAP")
    print("=" * 70)

    print("""
  Each 2-path (a,b,c) in a tournament involves 3 vertices.
  Its boundary d2(a,b,c) = (b,c) - (a,c)_or_(c,a) + (a,b) involves 3 edges.

  When we delete vertex v, we lose all 2-paths containing v.
  A 2-path (a,b,c) contains v iff v ∈ {a,b,c}.

  So the 2-paths containing v have v at position 0, 1, or 2.
  - Position 0: (v,b,c) — v is start
  - Position 1: (a,v,c) — v is middle
  - Position 2: (a,b,v) — v is end

  For distinct vertices v, w, a 2-path can contain BOTH v and w.
  """)

    for n in [5]:
        print(f"  n={n}: 2-path overlap statistics")
        for A in [worst_cases[0][0]] if worst_cases else []:
            info = compute_beta1_full(A, n)
            paths_2 = info['allowed_2']

            # Count paths per vertex
            paths_per_v = defaultdict(int)
            for p in paths_2:
                for v in p:
                    paths_per_v[v] += 1

            # Count paths per vertex pair
            paths_per_pair = Counter()
            for p in paths_2:
                for a, b in combinations(p, 2):
                    paths_per_pair[(a,b)] += 1
                    paths_per_pair[(b,a)] += 1

            # Count paths per vertex triple (= paths containing all 3 vertices of a triple)
            paths_per_triple = Counter()
            for p in paths_2:
                vs = tuple(sorted(p))
                paths_per_triple[vs] += 1

            print(f"    Total 2-paths: {len(paths_2)}")
            print(f"    Paths per vertex: {dict(paths_per_v)}")
            print(f"    Max shared paths (pair): {max(paths_per_pair.values()) if paths_per_pair else 0}")
            print(f"    Paths per triple: min={min(paths_per_triple.values())}, max={max(paths_per_triple.values())}")

    # ============================
    # PART 9: The deficiency formula
    # ============================
    print("\n" + "=" * 70)
    print("PART 9: DEFICIENCY FORMULA — beta_1(T\\v) in terms of rank drops")
    print("=" * 70)

    print("""
  For tournament T with beta_1(T) = 0 on n vertices:

  dim(Z_1(T)) = dim(B_1(T)) =: D
  dim(Z_1(T\\v)) =: D'  (should be constant = C(n-2,2) for all v)

  beta_1(T\\v) = D' - dim(B_1(T\\v))

  Since B_1(T\\v) ⊂ B_1(T) (restricted to T\\v edges):
  dim(B_1(T\\v)) = D' - beta_1(T\\v)

  So beta_1(T\\v) measures the "deficiency": how many independent cycles in T\\v
  are NOT boundaries of 2-paths in T\\v.

  Key: some boundaries of 2-paths in T INVOLVING v, when restricted to T\\v edges,
  give cycles in Z_1(T\\v) that are NOT in B_1(T\\v). These are the "phantom boundaries."
  """)

    # Compute the deficiency more carefully
    for n_val in [5]:
        print(f"\n  n={n_val}: Detailed deficiency analysis")
        for A, betas, idx in worst_cases[:3]:
            print(f"\n    Tournament idx={idx}, betas={betas}")
            info_T = compute_beta1_full(A, n_val)
            D = info_T['dim_B1']

            for v in range(n_val):
                Av, nv = subtournament(A, n_val, [v])
                info_Tv = compute_beta1_full(Av, nv)

                # The v-boundary space restricted to T\v edges
                # 2-paths through v in T: their boundaries restricted to T\v edges
                keep = [u for u in range(n_val) if u != v]
                edge_to_idx_T = {path: i for i, path in enumerate(info_T['allowed_1'])}
                edge_to_idx_Tv = {path: i for i, path in enumerate(info_Tv['allowed_1'])}

                # For each 2-path through v, compute its boundary and extract T\v components
                phantom_vecs = []
                for path in info_T['allowed_2']:
                    if v not in path:
                        continue
                    bd = boundary_coeffs(path)
                    # Extract components on T\v edges
                    vec = np.zeros(len(info_Tv['allowed_1']))
                    has_tv_component = False
                    for sign, face in bd:
                        if v not in face:
                            # This face is an edge in T\v
                            # Map to T\v labeling
                            mapped_face = tuple(keep.index(u) for u in face)
                            if mapped_face in edge_to_idx_Tv:
                                vec[edge_to_idx_Tv[mapped_face]] += sign
                                has_tv_component = True
                    if has_tv_component and np.any(np.abs(vec) > 1e-10):
                        phantom_vecs.append(vec)

                if phantom_vecs:
                    phantom_mat = np.column_stack(phantom_vecs)
                    phantom_rank = rank_of_matrix(phantom_mat)

                    # How many of these are already in B_1(T\v)?
                    B1_Tv = info_Tv['bd_2_omega']
                    if B1_Tv.shape[1] > 0:
                        combined = np.hstack([B1_Tv, phantom_mat])
                        rank_combined = rank_of_matrix(combined)
                        new_rank = rank_combined - rank_of_matrix(B1_Tv)
                    else:
                        new_rank = phantom_rank

                    # How many phantoms are cycles in T\v?
                    bd_1_Tv = build_full_boundary_matrix(info_Tv['allowed_1'], info_Tv['allowed_0'] if hasattr(info_Tv, 'allowed_0') else [(u,) for u in range(nv)])
                    # Actually compute d1 applied to phantom vecs
                    allowed_0_Tv = [(u,) for u in range(nv)]
                    bd_1_Tv = build_full_boundary_matrix(info_Tv['allowed_1'], allowed_0_Tv)
                    phantom_in_Z1 = 0
                    for vec in phantom_vecs:
                        if np.max(np.abs(bd_1_Tv @ vec)) < 1e-8:
                            phantom_in_Z1 += 1

                    print(f"    v={v}: #v-2paths={len(phantom_vecs)}, phantom_rank={phantom_rank}, "
                          f"new (not in B1(T\\v))={new_rank}, in_Z1(T\\v)={phantom_in_Z1}, "
                          f"beta1(T\\v)={betas[v]}")
                else:
                    print(f"    v={v}: no phantom boundaries, beta1(T\\v)={betas[v]}")

    # ============================
    # PART 10: Global constraint
    # ============================
    print("\n" + "=" * 70)
    print("PART 10: GLOBAL CONSTRAINT — TOTAL RANK BUDGET")
    print("=" * 70)

    print("""
  TOTAL RANK BUDGET ARGUMENT:

  Let T have n vertices and beta_1(T)=0.
  Then dim(B_1(T)) = dim(Z_1(T)).

  The 2-paths of T partition into:
  - Paths through exactly one vertex v (for each v)
  - Actually each 2-path goes through exactly 3 vertices.

  Consider the total: sum_v dim(V_v) where V_v = boundary image from v-paths.
  Each 2-path contributes to exactly 3 vertex spaces.
  So sum_v dim(V_v) <= 3 * dim(B_1(T)) (with equality if V_v's are independent).

  But dim(B_1(T)) = C(n-1,2), so sum_v dim(V_v) <= 3*C(n-1,2).

  Now: for each v, dim(V_v) >= rank_drop(v) >= n-2  (?)
  And beta_1(T\\v) relates to how much of Z_1(T\\v) is NOT covered by B_1(T\\v).
  """)

    # Compute sum_v dim(V_v) vs 3*dim(B_1(T))
    for n_val in [5]:
        print(f"\n  n={n_val}:")
        for A, betas, idx in worst_cases[:5]:
            info_T = compute_beta1_full(A, n_val)
            total_Vv = 0
            for v in range(n_val):
                del_info = analyze_deletion(A, n_val, v, info_T)
                total_Vv += del_info['dim_Vv']

            print(f"    idx={idx}: sum dim(V_v)={total_Vv}, "
                  f"3*dim(B1)={3*info_T['dim_B1']}, "
                  f"dim(B1)={info_T['dim_B1']}, "
                  f"sum_beta={sum(betas)}, betas={betas}")

    # ============================
    # FINAL SUMMARY
    # ============================
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)

    print(f"""
  n=5: MAX sum beta_1(T\\v) = {max_sum_beta1} (over all T with beta_1(T)=0)
  n=6: MAX sum beta_1(T\\v) = {max_sum_6} (over all T with beta_1(T)=0)

  The bound sum beta_1(T\\v) <= 3 {'HOLDS' if max(max_sum_beta1, max_sum_6) <= 3 else 'FAILS'} at n=5,6.
  """)

    # Characterize the extremal cases
    print("  Extremal cases at n=6 (sum = max):")
    for A, betas, idx in worst_6[:5]:
        if sum(betas) == max_sum_6:
            t3 = count_3cycles(A, n)
            print(f"    idx={idx}: t3={t3}, betas={betas}")

if __name__ == '__main__':
    run_analysis()
