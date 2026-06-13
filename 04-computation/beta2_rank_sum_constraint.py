#!/usr/bin/env python3
"""
Bad-vertex rank-drop constraints in tournament path homology.

Setup: Tournament T on n vertices with beta_1(T) = 0.
"Bad vertex" v means beta_1(T\v) = 1.
rank_drop(v) = dim(B_1(T)) - dim(B_1(T\v)) = (n-2) + beta_1(T\v)

Investigations:
1. Rank of v-contribution matrix (2-paths through v) vs rank_drop
2. Type A (v=middle) vs Type B+C (v=endpoint) decomposition
3. Bounds on #bad via Omega_2 structure
4. Joint distribution of (c_3, #bad)
5. Can we prove #bad <= n-2 via rank arguments?
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import sys
import time

# ---- Core homology routines (extracted from path_homology_v2.py) ----

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

def compute_beta1(A, n):
    """Compute beta_1 for digraph with adjacency matrix A."""
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_0 = enumerate_allowed_paths(A, n, 0)
    allowed_2 = enumerate_allowed_paths(A, n, 2)

    omega_1 = compute_omega_basis(A, n, 1, allowed_1, allowed_0)
    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)

    dim_omega_1 = omega_1.shape[1] if omega_1.ndim == 2 else 0
    if dim_omega_1 == 0:
        return 0

    bd_1 = build_full_boundary_matrix(allowed_1, allowed_0)
    bd_1_omega = bd_1 @ omega_1
    if bd_1_omega.shape[0] > 0 and bd_1_omega.shape[1] > 0:
        S1 = np.linalg.svd(bd_1_omega, compute_uv=False)
        rank_1 = sum(s > 1e-8 for s in S1)
    else:
        rank_1 = 0
    ker_dim = dim_omega_1 - rank_1

    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0
    if dim_omega_2 > 0:
        bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
        bd_2_omega = bd_2 @ omega_2
        S2 = np.linalg.svd(bd_2_omega, compute_uv=False)
        im_dim = sum(s > 1e-8 for s in S2)
    else:
        im_dim = 0

    return max(0, ker_dim - im_dim)

def compute_B1_rank(A, n):
    """Compute dim(B_1(T)) = rank of boundary map restricted to Omega_2."""
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_0 = enumerate_allowed_paths(A, n, 0)

    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0

    if dim_omega_2 == 0:
        return 0

    bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
    bd_2_omega = bd_2 @ omega_2

    if bd_2_omega.shape[0] > 0 and bd_2_omega.shape[1] > 0:
        S = np.linalg.svd(bd_2_omega, compute_uv=False)
        return sum(s > 1e-8 for s in S)
    return 0

def delete_vertex(A, n, v):
    """Return (A', n-1) with vertex v removed."""
    new_n = n - 1
    A_new = [[0]*new_n for _ in range(new_n)]
    def remap(i):
        return i if i < v else i - 1
    for i in range(n):
        if i == v: continue
        for j in range(n):
            if j == v: continue
            A_new[remap(i)][remap(j)] = A[i][j]
    return A_new, new_n

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
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def classify_2paths_through_v(allowed_2, v):
    """Classify 2-paths (a,b,c) by role of v:
    Type A: v = b (middle)
    Type B: v = a (first)
    Type C: v = c (last)
    """
    type_A = []  # v is middle
    type_B = []  # v is first
    type_C = []  # v is last
    for idx, path in enumerate(allowed_2):
        a, b, c = path
        if b == v:
            type_A.append(idx)
        elif a == v:
            type_B.append(idx)
        elif c == v:
            type_C.append(idx)
    return type_A, type_B, type_C

def matrix_rank(M, tol=1e-8):
    if M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    S = np.linalg.svd(M, compute_uv=False)
    return int(sum(s > tol for s in S))


# ===== MAIN INVESTIGATION =====

def investigate_tournament(A, n, verbose=False):
    """Full rank-drop analysis for one tournament."""
    beta1_T = compute_beta1(A, n)
    if beta1_T != 0:
        return None  # skip tournaments with beta_1 != 0

    # Compute B_1(T) info
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_0 = enumerate_allowed_paths(A, n, 0)

    omega_2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    dim_omega_2 = omega_2.shape[1] if omega_2.ndim == 2 else 0

    bd_2_full = build_full_boundary_matrix(allowed_2, allowed_1)

    if dim_omega_2 > 0:
        bd_2_omega = bd_2_full @ omega_2
        B1_rank_T = matrix_rank(bd_2_omega)
    else:
        bd_2_omega = np.zeros((len(allowed_1), 0))
        B1_rank_T = 0

    # For each vertex, compute rank_drop and type decomposition
    results = []
    for v in range(n):
        A_del, n_del = delete_vertex(A, n, v)
        beta1_del = compute_beta1(A_del, n_del)
        rank_drop = (n - 2) + beta1_del  # by formula

        # Also compute B1_rank directly for T\v
        B1_rank_del = compute_B1_rank(A_del, n_del)
        rank_drop_direct = B1_rank_T - B1_rank_del

        # Classify 2-paths through v
        type_A, type_B, type_C = classify_2paths_through_v(allowed_2, v)

        # Build the boundary sub-matrix for paths through v (in Omega_2)
        # First, identify which Omega_2 basis vectors involve v-paths
        through_v_indices = set(type_A + type_B + type_C)

        # The contribution of v-paths to the boundary
        # In the full A_2 basis, select columns for v-paths
        if through_v_indices and dim_omega_2 > 0:
            # Project omega_2 onto the v-path subspace
            # omega_2 has shape (|A_2|, dim_omega_2)
            # We want the boundary image of the omega_2 subspace
            # restricted to paths through v... but omega_2 mixes paths.

            # Alternative: compute rank contribution directly
            # rank_drop = rank(bd_2_omega) - rank(bd_2_omega restricted to non-v-paths)

            # Identify non-v-path indices in A_2
            not_through_v = [i for i in range(len(allowed_2)) if i not in through_v_indices]

            # Project omega_2 onto non-v-paths: zero out v-path rows
            omega_2_no_v = omega_2.copy()
            for idx in through_v_indices:
                omega_2_no_v[idx, :] = 0

            bd_2_no_v = bd_2_full @ omega_2_no_v
            rank_no_v = matrix_rank(bd_2_no_v)

            rank_from_v = B1_rank_T - rank_no_v
        else:
            rank_from_v = 0

        # Type A only contribution (middle vertex paths)
        if type_A and dim_omega_2 > 0:
            omega_2_no_A = omega_2.copy()
            for idx in type_A:
                omega_2_no_A[idx, :] = 0
            bd_2_no_A = bd_2_full @ omega_2_no_A
            rank_no_A = matrix_rank(bd_2_no_A)
            rank_from_A = B1_rank_T - rank_no_A
        else:
            rank_from_A = 0

        # Type BC only contribution
        type_BC = type_B + type_C
        if type_BC and dim_omega_2 > 0:
            omega_2_no_BC = omega_2.copy()
            for idx in type_BC:
                omega_2_no_BC[idx, :] = 0
            bd_2_no_BC = bd_2_full @ omega_2_no_BC
            rank_no_BC = matrix_rank(bd_2_no_BC)
            rank_from_BC = B1_rank_T - rank_no_BC
        else:
            rank_from_BC = 0

        # Score of v
        score_out = sum(A[v][j] for j in range(n) if j != v)

        results.append({
            'v': v,
            'beta1_del': beta1_del,
            'rank_drop_formula': rank_drop,
            'rank_drop_direct': rank_drop_direct,
            'is_bad': beta1_del > 0,
            'n_type_A': len(type_A),
            'n_type_B': len(type_B),
            'n_type_C': len(type_C),
            'rank_from_v': rank_from_v,
            'rank_from_A': rank_from_A,
            'rank_from_BC': rank_from_BC,
            'score_out': score_out,
        })

    t3 = count_3cycles(A, n)
    n_bad = sum(r['is_bad'] for r in results)
    sum_rank_drop = sum(r['rank_drop_formula'] for r in results)

    return {
        'n': n,
        't3': t3,
        'beta1_T': beta1_T,
        'B1_rank_T': B1_rank_T,
        'dim_omega_2': dim_omega_2,
        'n_2paths': len(allowed_2),
        'n_bad': n_bad,
        'sum_rank_drop': sum_rank_drop,
        'expected_sum': n*(n-2) + n_bad,
        'vertices': results,
    }


def print_detailed(info, show_vertices=True):
    """Print detailed results for one tournament."""
    print(f"  t3={info['t3']}, B1_rank={info['B1_rank_T']}, "
          f"dim(Omega_2)={info['dim_omega_2']}, #2paths={info['n_2paths']}, "
          f"#bad={info['n_bad']}, sum_rank_drop={info['sum_rank_drop']} "
          f"(expected n(n-2)+#bad = {info['expected_sum']})")
    if show_vertices:
        for r in info['vertices']:
            bad_mark = " *BAD*" if r['is_bad'] else ""
            print(f"    v={r['v']}: score={r['score_out']}, "
                  f"rank_drop={r['rank_drop_formula']} (direct={r['rank_drop_direct']}), "
                  f"types A/B/C={r['n_type_A']}/{r['n_type_B']}/{r['n_type_C']}, "
                  f"rank_from_v={r['rank_from_v']}, "
                  f"rank_from_A={r['rank_from_A']}, rank_from_BC={r['rank_from_BC']}"
                  f"{bad_mark}")


def run_exhaustive(n, verbose_limit=5):
    """Run exhaustive analysis for all tournaments of size n."""
    print(f"\n{'='*70}")
    print(f"EXHAUSTIVE ANALYSIS n={n}")
    print(f"{'='*70}")

    all_infos = []
    count = 0
    skip_count = 0
    t0 = time.time()

    for A in all_tournaments(n):
        info = investigate_tournament(A, n)
        count += 1
        if info is None:
            skip_count += 1
            continue
        all_infos.append(info)

        if len(all_infos) <= verbose_limit:
            print(f"\nTournament #{count} (beta_1=0):")
            print_detailed(info, show_vertices=(len(all_infos) <= 3))

    elapsed = time.time() - t0
    print(f"\nProcessed {count} tournaments in {elapsed:.1f}s")
    print(f"  beta_1=0: {len(all_infos)}, beta_1>0: {skip_count}")

    if not all_infos:
        return

    # ---- INVESTIGATION 1: Rank drop consistency ----
    print(f"\n--- INVESTIGATION 1: Rank drop formula consistency ---")
    consistent = all(
        all(r['rank_drop_formula'] == r['rank_drop_direct'] for r in info['vertices'])
        for info in all_infos
    )
    print(f"  rank_drop = (n-2) + beta_1(T\\v) matches direct computation: {consistent}")

    # ---- INVESTIGATION 2: Type A vs B+C contribution ----
    print(f"\n--- INVESTIGATION 2: Type A vs B+C rank contributions ---")
    for info in all_infos[:5]:
        for r in info['vertices']:
            print(f"  t3={info['t3']}, v={r['v']}: "
                  f"rank_from_v={r['rank_from_v']}, "
                  f"rank_from_A={r['rank_from_A']}, "
                  f"rank_from_BC={r['rank_from_BC']}, "
                  f"rank_drop={r['rank_drop_formula']}, "
                  f"bad={r['is_bad']}")

    # ---- INVESTIGATION 3: Sum of rank drops ----
    print(f"\n--- INVESTIGATION 3: Sum of rank drops ---")
    sum_drops = [info['sum_rank_drop'] for info in all_infos]
    n_bads = [info['n_bad'] for info in all_infos]
    print(f"  sum_rank_drop range: [{min(sum_drops)}, {max(sum_drops)}]")
    print(f"  n(n-2) = {n*(n-2)}")
    print(f"  #bad range: [{min(n_bads)}, {max(n_bads)}]")
    print(f"  #bad distribution: {Counter(n_bads)}")

    # ---- INVESTIGATION 4: Joint distribution (c3, #bad) ----
    print(f"\n--- INVESTIGATION 4: Joint distribution (c3, #bad) ---")
    joint = defaultdict(int)
    for info in all_infos:
        joint[(info['t3'], info['n_bad'])] += 1

    print(f"  {'t3':>4} {'#bad':>5} {'count':>6}")
    for (t3, nb) in sorted(joint.keys()):
        print(f"  {t3:>4} {nb:>5} {joint[(t3,nb)]:>6}")

    # ---- INVESTIGATION 5: Is #bad <= n-2? ----
    print(f"\n--- INVESTIGATION 5: #bad <= n-2 bound ---")
    max_bad = max(n_bads) if n_bads else 0
    print(f"  max(#bad) = {max_bad}, n-2 = {n-2}")
    print(f"  #bad <= n-2 holds: {max_bad <= n-2}")

    # Additional: bad vertex scores
    print(f"\n--- INVESTIGATION 5b: Bad vertex properties ---")
    bad_scores = []
    good_scores = []
    for info in all_infos:
        for r in info['vertices']:
            if r['is_bad']:
                bad_scores.append(r['score_out'])
            else:
                good_scores.append(r['score_out'])
    if bad_scores:
        print(f"  Bad vertex out-degree distribution: {Counter(bad_scores)}")
    if good_scores:
        print(f"  Good vertex out-degree distribution: {Counter(good_scores)}")

    # ---- INVESTIGATION 6: Relationship between dim(Omega_2), #2paths, and #bad ----
    print(f"\n--- INVESTIGATION 6: Omega_2 structure vs #bad ---")
    for info in all_infos[:10]:
        # How many 2-paths are in Omega_2 vs total allowed 2-paths?
        print(f"  t3={info['t3']}, #2paths={info['n_2paths']}, "
              f"dim(Omega_2)={info['dim_omega_2']}, "
              f"B1_rank={info['B1_rank_T']}, #bad={info['n_bad']}")

    # ---- INVESTIGATION 7: Bad vertex overlap with scores ----
    print(f"\n--- INVESTIGATION 7: Which vertices are bad? ---")
    for info in all_infos:
        if info['n_bad'] > 0:
            bad_verts = [r['v'] for r in info['vertices'] if r['is_bad']]
            bad_scores_list = [r['score_out'] for r in info['vertices'] if r['is_bad']]
            all_scores = [r['score_out'] for r in info['vertices']]
            print(f"  t3={info['t3']}, scores={all_scores}, "
                  f"bad_verts={bad_verts} (scores {bad_scores_list})")
            if len([i for i in all_infos if i['n_bad'] > 0]) > 20:
                break  # don't flood

    # ---- INVESTIGATION 8: Double-counting argument ----
    print(f"\n--- INVESTIGATION 8: Double-counting of 2-paths ---")
    print("  Each 2-path (a,b,c) touches 3 vertices: a (first), b (middle), c (last)")
    print("  Sum over v of #(2-paths through v) = 3 * #(2-paths)")
    for info in all_infos[:5]:
        total_through = sum(
            r['n_type_A'] + r['n_type_B'] + r['n_type_C']
            for r in info['vertices']
        )
        print(f"  t3={info['t3']}: sum_through_v = {total_through}, "
              f"3*#2paths = {3*info['n_2paths']}, match: {total_through == 3*info['n_2paths']}")

    # ---- INVESTIGATION 9: rank_from_v vs rank_drop ----
    print(f"\n--- INVESTIGATION 9: rank_from_v (zeroing approach) vs rank_drop ---")
    print("  rank_from_v = rank(bd @ omega) - rank(bd @ omega with v-paths zeroed)")
    print("  rank_drop = B1_rank(T) - B1_rank(T\\v)")
    print("  These measure different things: zeroing vs deletion")
    diffs = []
    for info in all_infos:
        for r in info['vertices']:
            diff = r['rank_from_v'] - r['rank_drop_formula']
            diffs.append(diff)
            if abs(diff) > 0 and len(diffs) <= 20:
                print(f"  t3={info['t3']}, v={r['v']}: "
                      f"rank_from_v={r['rank_from_v']}, rank_drop={r['rank_drop_formula']}, "
                      f"DIFF={diff}")
    print(f"  Total vertices: {len(diffs)}")
    print(f"  Diff=0: {sum(1 for d in diffs if d == 0)}")
    print(f"  Diff distribution: {Counter(diffs)}")

    return all_infos


# ===== MAIN =====
print("="*70)
print("BAD VERTEX RANK-DROP CONSTRAINTS IN TOURNAMENT PATH HOMOLOGY")
print("="*70)

# n=5 exhaustive
infos_5 = run_exhaustive(5, verbose_limit=3)

# n=6 exhaustive
infos_6 = run_exhaustive(6, verbose_limit=3)

# n=7: sample (exhaustive is 2^21 = 2M tournaments, too slow)
print(f"\n{'='*70}")
print(f"SAMPLED ANALYSIS n=7")
print(f"{'='*70}")

import random
random.seed(42)
n = 7
edges_7 = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges_7)

sample_size = 500
all_infos_7 = []
skip_7 = 0
t0 = time.time()

for trial in range(sample_size):
    A = [[0]*n for _ in range(n)]
    for (i,j) in edges_7:
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1
    info = investigate_tournament(A, n)
    if info is None:
        skip_7 += 1
        continue
    all_infos_7.append(info)

elapsed = time.time() - t0
print(f"\nSampled {sample_size} tournaments in {elapsed:.1f}s")
print(f"  beta_1=0: {len(all_infos_7)}, beta_1>0: {skip_7}")

if all_infos_7:
    n_bads_7 = [info['n_bad'] for info in all_infos_7]
    print(f"\n--- n=7 #bad distribution ---")
    print(f"  #bad range: [{min(n_bads_7)}, {max(n_bads_7)}]")
    print(f"  #bad <= n-2 = 5 holds: {max(n_bads_7) <= 5}")
    print(f"  Distribution: {Counter(n_bads_7)}")

    # Joint (t3, #bad)
    joint_7 = defaultdict(int)
    for info in all_infos_7:
        joint_7[(info['t3'], info['n_bad'])] += 1
    print(f"\n--- n=7 Joint distribution (c3, #bad) ---")
    t3_vals = sorted(set(t3 for t3, _ in joint_7.keys()))
    nb_vals = sorted(set(nb for _, nb in joint_7.keys()))
    # Print as table
    header = f"  {'t3':>4} | " + " ".join(f"b={nb}" for nb in nb_vals)
    print(header)
    print("  " + "-"*len(header))
    for t3 in t3_vals:
        row = f"  {t3:>4} | " + " ".join(f"{joint_7.get((t3,nb),0):>4}" for nb in nb_vals)
        print(row)

    # Consistency checks
    consistent = all(
        all(r['rank_drop_formula'] == r['rank_drop_direct'] for r in info['vertices'])
        for info in all_infos_7
    )
    print(f"\n  Rank drop formula consistent: {consistent}")

    # Bad vertex scores at n=7
    bad_scores_7 = []
    for info in all_infos_7:
        for r in info['vertices']:
            if r['is_bad']:
                bad_scores_7.append(r['score_out'])
    if bad_scores_7:
        print(f"  Bad vertex out-degree dist: {Counter(bad_scores_7)}")


# ===== SYNTHESIS =====
print(f"\n{'='*70}")
print("SYNTHESIS AND CONJECTURES")
print("="*70)

for n_val, infos in [(5, infos_5), (6, infos_6), (7, all_infos_7)]:
    if not infos:
        continue
    n_bads = [info['n_bad'] for info in infos]
    max_bad = max(n_bads)
    print(f"\nn={n_val}: max(#bad) = {max_bad}, n-2 = {n_val-2}")
    print(f"  #bad distribution: {Counter(n_bads)}")

    # Check: is #bad always even?
    parities = Counter(nb % 2 for nb in n_bads)
    print(f"  #bad parity: {parities}")

    # Check: correlation of #bad with t3
    t3s = [info['t3'] for info in infos]
    if len(set(t3s)) > 1 and len(set(n_bads)) > 1:
        mean_t3 = np.mean(t3s)
        mean_nb = np.mean(n_bads)
        cov = np.mean([(t - mean_t3)*(b - mean_nb) for t, b in zip(t3s, n_bads)])
        std_t3 = np.std(t3s)
        std_nb = np.std(n_bads)
        if std_t3 > 0 and std_nb > 0:
            corr = cov / (std_t3 * std_nb)
            print(f"  corr(t3, #bad) = {corr:.4f}")

print("\n\nDone.")
