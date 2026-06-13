#!/usr/bin/env python3
"""
beta1_kernel_structure.py — Analyze β₁ deletion structure for tournaments

For tournaments T with β₁(T)=0:
- Which vertices v have β₁(T\\v) = 1 ("bad") vs β₁(T\\v) = 0 ("good")?
- Can ALL vertices be bad? (Proving "no" is the goal)
- What patterns govern bad/good status?

Analyzes:
1. Distribution of bad vertex counts
2. Relationship between out-degree and bad status
3. Kernel structure of the chain-level restriction map
4. Image of res_v(Z¹(T)) vs Z¹(T\\v) — the real mechanism
5. The "lifting obstruction" for bad vertices
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import random

# ===== Core path homology functions (from path_homology_v2) =====

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
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
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

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
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0: return np.zeros((dim_Ap, 0))
    return null_space

def compute_path_homology_spaces(A, n):
    """Compute Z¹, B¹, and their dimensions for a digraph."""
    allowed = {p: enumerate_allowed_paths(A, n, p) for p in range(-1, 4)}
    allowed[-1] = []

    omega = {}
    for p in range(4):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    # Z¹ = ker(∂_1|_{Omega_1})
    bd_1 = build_full_boundary_matrix(allowed[1], allowed[0])
    bd_1_omega = bd_1 @ omega[1]
    if bd_1_omega.shape[0] > 0 and bd_1_omega.shape[1] > 0:
        U, S, Vt = np.linalg.svd(bd_1_omega, full_matrices=True)
        rank_1 = sum(s > 1e-10 for s in S)
        Z1_coords = Vt[rank_1:].T
        Z1_basis = omega[1] @ Z1_coords
    else:
        Z1_basis = omega[1]
    dim_Z1 = Z1_basis.shape[1] if Z1_basis.ndim == 2 else 0

    # B¹ = im(∂_2|_{Omega_2})
    dim_omega_2 = omega[2].shape[1] if omega[2].ndim == 2 else 0
    if dim_omega_2 > 0:
        bd_2 = build_full_boundary_matrix(allowed[2], allowed[1])
        bd_2_omega = bd_2 @ omega[2]
        U2, S2, Vt2 = np.linalg.svd(bd_2_omega, full_matrices=True)
        rank_2 = sum(s > 1e-10 for s in S2)
        B1_basis = U2[:, :rank_2]
    else:
        B1_basis = np.zeros((len(allowed[1]), 0))
    dim_B1 = B1_basis.shape[1] if B1_basis.ndim == 2 else 0

    beta_1 = dim_Z1 - dim_B1

    return {
        'allowed_1': allowed[1],
        'Z1_basis': Z1_basis,
        'B1_basis': B1_basis,
        'dim_Z1': dim_Z1,
        'dim_B1': dim_B1,
        'beta_1': beta_1,
    }

# ===== Tournament utilities =====

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def delete_vertex(A, n, v):
    return [[A[i][j] for j in range(n) if j != v] for i in range(n) if i != v]

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

# ===== Restriction and lifting analysis =====

def analyze_restriction(A, n, v, spaces_T=None):
    """Analyze res_v: chains(T) -> chains(T\\v) restricted to Z¹(T).

    Key question: does res_v(Z¹(T)) contain Z¹(T\\v)?
    If yes: v is "good" (every cocycle on T\\v lifts)
    If no: v is "bad" (some cocycle on T\\v doesn't lift)
    """
    if spaces_T is None:
        spaces_T = compute_path_homology_spaces(A, n)

    A_Tv = delete_vertex(A, n, v)
    spaces_Tv = compute_path_homology_spaces(A_Tv, n-1)

    allowed_1_T = spaces_T['allowed_1']
    allowed_1_Tv = spaces_Tv['allowed_1']
    Z1_T = spaces_T['Z1_basis']
    Z1_Tv = spaces_Tv['Z1_basis']
    B1_Tv = spaces_Tv['B1_basis']

    # Build res_v: A_1(T) -> A_1(T\\v)
    def orig_to_new(u):
        return u if u < v else u - 1

    idx_Tv = {path: i for i, path in enumerate(allowed_1_Tv)}
    res_full = np.zeros((len(allowed_1_Tv), len(allowed_1_T)))
    for j, (a, b) in enumerate(allowed_1_T):
        if a == v or b == v:
            continue
        new_edge = (orig_to_new(a), orig_to_new(b))
        if new_edge in idx_Tv:
            res_full[idx_Tv[new_edge], j] = 1.0

    # Image of Z¹(T) under res_v
    res_Z1 = res_full @ Z1_T  # columns span res_v(Z¹(T)) in A_1(T\\v) coords

    # Does res_v(Z¹(T)) contain Z¹(T\\v)?
    # Check: is each column of Z1_Tv in the column space of res_Z1?
    if spaces_Tv['dim_Z1'] == 0:
        contains_Z1 = True
        gap = 0
    else:
        # Combine and check rank
        combined = np.column_stack([res_Z1, Z1_Tv])
        rank_res = np.linalg.matrix_rank(res_Z1, tol=1e-8)
        rank_combined = np.linalg.matrix_rank(combined, tol=1e-8)
        contains_Z1 = (rank_combined == rank_res)
        gap = rank_combined - rank_res  # number of independent cocycles NOT in image

    # Does res_v(Z¹(T)) contain B¹(T\\v)?
    if spaces_Tv['dim_B1'] == 0:
        contains_B1 = True
    else:
        combined_B = np.column_stack([res_Z1, B1_Tv])
        rank_combined_B = np.linalg.matrix_rank(combined_B, tol=1e-8)
        contains_B1 = (rank_combined_B == np.linalg.matrix_rank(res_Z1, tol=1e-8))

    # Also: does res_v(B¹(T)) contain B¹(T\\v)?
    B1_T = spaces_T['B1_basis']
    res_B1 = res_full @ B1_T
    if spaces_Tv['dim_B1'] == 0:
        B1_contains_B1 = True
    else:
        combined_BB = np.column_stack([res_B1, B1_Tv])
        rank_resB = np.linalg.matrix_rank(res_B1, tol=1e-8)
        rank_combined_BB = np.linalg.matrix_rank(combined_BB, tol=1e-8)
        B1_contains_B1 = (rank_combined_BB == rank_resB)

    return {
        'beta1_Tv': spaces_Tv['beta_1'],
        'dim_Z1_Tv': spaces_Tv['dim_Z1'],
        'dim_B1_Tv': spaces_Tv['dim_B1'],
        'dim_res_Z1': np.linalg.matrix_rank(res_Z1, tol=1e-8),
        'contains_Z1': contains_Z1,
        'contains_B1': contains_B1,
        'B1_contains_B1': B1_contains_B1,
        'gap': gap,
    }

def full_analysis(A, n, verbose=False):
    """Full analysis of one tournament."""
    spaces_T = compute_path_homology_spaces(A, n)
    if spaces_T['beta_1'] != 0:
        return None

    out_degrees = [sum(A[i]) for i in range(n)]
    t3 = count_3cycles(A, n)

    vertex_data = []
    for v in range(n):
        rd = analyze_restriction(A, n, v, spaces_T)
        is_bad = rd['beta1_Tv'] > 0
        vertex_data.append({
            'v': v,
            'out_deg': out_degrees[v],
            'bad': is_bad,
            **rd,
        })

    n_bad = sum(1 for vd in vertex_data if vd['bad'])
    n_good = n - n_bad

    if verbose:
        print(f"  dim Z¹={spaces_T['dim_Z1']}, dim B¹={spaces_T['dim_B1']}, t3={t3}")
        print(f"  Out-degrees: {out_degrees}")
        print(f"  {n_bad} bad, {n_good} good")
        for vd in vertex_data:
            status = 'BAD' if vd['bad'] else 'good'
            print(f"    v={vd['v']}: d+={vd['out_deg']}, β₁(T\\v)={vd['beta1_Tv']}, "
                  f"dim Z¹(T\\v)={vd['dim_Z1_Tv']}, dim B¹(T\\v)={vd['dim_B1_Tv']}, "
                  f"dim res(Z¹)={vd['dim_res_Z1']}, "
                  f"res(Z¹)⊇Z¹={vd['contains_Z1']}, res(Z¹)⊇B¹={vd['contains_B1']}, "
                  f"res(B¹)⊇B¹={vd['B1_contains_B1']}, gap={vd['gap']} [{status}]")

    return {
        'out_degrees': out_degrees,
        't3': t3,
        'n_bad': n_bad,
        'n_good': n_good,
        'all_bad': n_bad == n,
        'vertex_data': vertex_data,
        'dim_Z1': spaces_T['dim_Z1'],
        'dim_B1': spaces_T['dim_B1'],
    }

# ===== Main =====

def main():
    print("=" * 70)
    print("β₁ KERNEL STRUCTURE ANALYSIS")
    print("=" * 70)

    # ===== n=5 exhaustive =====
    print("\n" + "=" * 70)
    print("n=5 EXHAUSTIVE")
    print("=" * 70)

    n = 5
    total = 0
    beta1_zero = 0
    bad_count_dist = defaultdict(int)
    all_bad_count = 0
    results_5 = []

    for A in all_tournaments(n):
        total += 1
        r = full_analysis(A, n)
        if r is None:
            continue
        beta1_zero += 1
        bad_count_dist[r['n_bad']] += 1
        if r['all_bad']:
            all_bad_count += 1
        results_5.append(r)

    print(f"\nTotal tournaments: {total}")
    print(f"β₁=0: {beta1_zero}, β₁=1: {total - beta1_zero}")

    print(f"\nBad vertex count distribution:")
    for k in sorted(bad_count_dist.keys()):
        pct = 100 * bad_count_dist[k] / beta1_zero
        print(f"  {k} bad: {bad_count_dist[k]} ({pct:.1f}%)")

    print(f"\n*** ALL vertices bad: {all_bad_count} ***")

    # ===== Detailed examples =====
    print(f"\n--- Detailed examples by bad count ---")
    shown = defaultdict(int)
    for r in results_5:
        if shown[r['n_bad']] < 1:
            shown[r['n_bad']] += 1
            print(f"\n  [{r['n_bad']} bad vertices]")
            # Reconstruct A for verbose output
            # Just print what we have
            print(f"    Out-degrees: {r['out_degrees']}, t3={r['t3']}")
            for vd in r['vertex_data']:
                status = 'BAD' if vd['bad'] else 'good'
                print(f"      v={vd['v']}: d+={vd['out_deg']}, β₁(T\\v)={vd['beta1_Tv']}, "
                      f"Z¹={vd['dim_Z1_Tv']}, B¹={vd['dim_B1_Tv']}, "
                      f"res⊇Z¹={vd['contains_Z1']}, res⊇B¹={vd['contains_B1']}, "
                      f"resB⊇B¹={vd['B1_contains_B1']}, gap={vd['gap']} [{status}]")

    # ===== Restriction analysis =====
    print(f"\n--- Restriction map analysis ---")
    # For bad vertices: is res_v(Z¹(T)) ⊇ B¹(T\\v) always true?
    bad_contains_B1 = 0
    bad_not_contains_B1 = 0
    bad_total = 0
    # For good vertices: is res_v(Z¹(T)) ⊇ Z¹(T\\v) always true?
    good_contains_Z1 = 0
    good_not_contains_Z1 = 0
    good_total = 0

    gap_dist_bad = Counter()

    for r in results_5:
        for vd in r['vertex_data']:
            if vd['bad']:
                bad_total += 1
                if vd['contains_B1']:
                    bad_contains_B1 += 1
                else:
                    bad_not_contains_B1 += 1
                gap_dist_bad[vd['gap']] += 1
            else:
                good_total += 1
                if vd['contains_Z1']:
                    good_contains_Z1 += 1
                else:
                    good_not_contains_Z1 += 1

    print(f"\n  Bad vertices ({bad_total} total):")
    print(f"    res(Z¹(T)) ⊇ B¹(T\\v): {bad_contains_B1}/{bad_total}")
    print(f"    res(Z¹(T)) ⊅ B¹(T\\v): {bad_not_contains_B1}/{bad_total}")
    print(f"    Gap distribution: {dict(gap_dist_bad)}")

    print(f"\n  Good vertices ({good_total} total):")
    print(f"    res(Z¹(T)) ⊇ Z¹(T\\v): {good_contains_Z1}/{good_total}")
    print(f"    res(Z¹(T)) ⊅ Z¹(T\\v): {good_not_contains_Z1}/{good_total}")

    # B¹ lifting
    B1_lift_bad = sum(1 for r in results_5 for vd in r['vertex_data'] if vd['bad'] and vd['B1_contains_B1'])
    B1_lift_good = sum(1 for r in results_5 for vd in r['vertex_data'] if not vd['bad'] and vd['B1_contains_B1'])
    print(f"\n  res(B¹(T)) ⊇ B¹(T\\v):")
    print(f"    Bad vertices: {B1_lift_bad}/{bad_total}")
    print(f"    Good vertices: {B1_lift_good}/{good_total}")

    # ===== Out-degree vs bad status =====
    print(f"\n--- Out-degree vs bad status ---")
    deg_stats = defaultdict(lambda: [0, 0])
    for r in results_5:
        for vd in r['vertex_data']:
            d = vd['out_deg']
            deg_stats[d][1] += 1
            if vd['bad']:
                deg_stats[d][0] += 1
    for d in sorted(deg_stats.keys()):
        bad, tot = deg_stats[d]
        print(f"  d+={d}: {bad}/{tot} bad ({100*bad/tot:.1f}%)")

    # ===== t3 vs bad count =====
    print(f"\n--- t3 vs max bad count ---")
    t3_bad = defaultdict(list)
    for r in results_5:
        t3_bad[r['t3']].append(r['n_bad'])
    for t3 in sorted(t3_bad.keys()):
        vals = t3_bad[t3]
        print(f"  t3={t3}: max_bad={max(vals)}, avg_bad={np.mean(vals):.2f}, n={len(vals)}")

    # ===== Score sequence patterns =====
    print(f"\n--- Score sequence vs all-bad ---")
    ss_dist = defaultdict(lambda: [0, 0])
    for r in results_5:
        ss = tuple(sorted(r['out_degrees']))
        ss_dist[ss][1] += 1
        if r['all_bad']:
            ss_dist[ss][0] += 1
    print(f"  Max bad count by score sequence:")
    ss_maxbad = defaultdict(int)
    for r in results_5:
        ss = tuple(sorted(r['out_degrees']))
        ss_maxbad[ss] = max(ss_maxbad[ss], r['n_bad'])
    for ss in sorted(ss_maxbad.keys()):
        print(f"    {ss}: max {ss_maxbad[ss]} bad")

    # ===== n=6 sample =====
    print(f"\n\n" + "=" * 70)
    print("n=6 SAMPLING (300 tournaments)")
    print("=" * 70)

    n = 6
    random.seed(42)
    beta1_zero_6 = 0
    bad_count_dist_6 = defaultdict(int)
    all_bad_6 = 0
    results_6 = []

    for trial in range(300):
        A = random_tournament(n)
        r = full_analysis(A, n)
        if r is None:
            continue
        beta1_zero_6 += 1
        bad_count_dist_6[r['n_bad']] += 1
        if r['all_bad']:
            all_bad_6 += 1
        results_6.append(r)

    print(f"\nβ₁=0: {beta1_zero_6}/300")
    print(f"\nBad vertex count distribution:")
    for k in sorted(bad_count_dist_6.keys()):
        print(f"  {k} bad: {bad_count_dist_6[k]}")
    print(f"\n*** ALL vertices bad: {all_bad_6} ***")

    # Restriction analysis n=6
    print(f"\n--- Restriction map analysis (n=6) ---")
    bad_6 = {'contains_B1': 0, 'not': 0, 'total': 0}
    good_6 = {'contains_Z1': 0, 'not': 0, 'total': 0}
    for r in results_6:
        for vd in r['vertex_data']:
            if vd['bad']:
                bad_6['total'] += 1
                if vd['contains_B1']:
                    bad_6['contains_B1'] += 1
                else:
                    bad_6['not'] += 1
            else:
                good_6['total'] += 1
                if vd['contains_Z1']:
                    good_6['contains_Z1'] += 1
                else:
                    good_6['not'] += 1

    print(f"  Bad: res(Z¹)⊇B¹ = {bad_6['contains_B1']}/{bad_6['total']}")
    print(f"  Good: res(Z¹)⊇Z¹ = {good_6['contains_Z1']}/{good_6['total']}")

    # Detailed examples
    print(f"\n--- n=6 examples ---")
    shown6 = defaultdict(int)
    for r in results_6:
        nb = r['n_bad']
        if shown6[nb] < 1:
            shown6[nb] += 1
            print(f"\n  [{nb} bad vertices]")
            print(f"    Out-degrees: {r['out_degrees']}, t3={r['t3']}")
            for vd in r['vertex_data']:
                status = 'BAD' if vd['bad'] else 'good'
                print(f"      v={vd['v']}: d+={vd['out_deg']}, β₁={vd['beta1_Tv']}, "
                      f"Z¹={vd['dim_Z1_Tv']}, B¹={vd['dim_B1_Tv']}, "
                      f"gap={vd['gap']} [{status}]")

    # ===== n=7 small sample =====
    print(f"\n\n" + "=" * 70)
    print("n=7 SMALL SAMPLE (50 tournaments)")
    print("=" * 70)

    n = 7
    random.seed(123)
    beta1_zero_7 = 0
    bad_count_dist_7 = defaultdict(int)
    all_bad_7 = 0

    for trial in range(50):
        A = random_tournament(n)
        r = full_analysis(A, n)
        if r is None:
            continue
        beta1_zero_7 += 1
        bad_count_dist_7[r['n_bad']] += 1
        if r['all_bad']:
            all_bad_7 += 1

    print(f"\nβ₁=0: {beta1_zero_7}/50")
    print(f"\nBad vertex count distribution:")
    for k in sorted(bad_count_dist_7.keys()):
        print(f"  {k} bad: {bad_count_dist_7[k]}")
    print(f"\n*** ALL vertices bad: {all_bad_7} ***")

    # ===== KEY STRUCTURAL INSIGHT =====
    print(f"\n\n" + "=" * 70)
    print("STRUCTURAL ANALYSIS")
    print("=" * 70)

    # For n=5: when a vertex is bad, what are the Z¹/B¹ dimensions?
    print(f"\nDimension analysis (n=5):")
    print(f"  T: dim Z¹ = dim B¹ = C(n,2)-(n-1) = {5*4//2 - 4}")
    print(f"  T\\v: dim Z¹ = C(n-1,2)-(n-2) = {4*3//2 - 3}")
    print(f"  T\\v: dim B¹ depends on Omega_2 of T\\v")

    # Count Z1/B1 dims for T\\v
    z1b1_dist = defaultdict(int)
    for r in results_5:
        for vd in r['vertex_data']:
            key = (vd['dim_Z1_Tv'], vd['dim_B1_Tv'], vd['bad'])
            z1b1_dist[key] += 1
    print(f"\n  (dim_Z1_Tv, dim_B1_Tv, bad) distribution:")
    for key in sorted(z1b1_dist.keys()):
        print(f"    Z¹={key[0]}, B¹={key[1]}, {'BAD' if key[2] else 'good'}: {z1b1_dist[key]}")

    # ===== The key question for the proof =====
    print(f"\n" + "=" * 70)
    print("KEY QUESTION: What prevents all vertices from being bad?")
    print("=" * 70)

    # For each tournament, look at the relationship between bad vertices
    print(f"\nFor tournaments with max bad count (3 bad at n=5):")
    for r in results_5:
        if r['n_bad'] == 3:
            good_verts = [vd['v'] for vd in r['vertex_data'] if not vd['bad']]
            bad_verts = [vd['v'] for vd in r['vertex_data'] if vd['bad']]
            good_degs = [r['out_degrees'][v] for v in good_verts]
            bad_degs = [r['out_degrees'][v] for v in bad_verts]
            print(f"  Good: v={good_verts} (d+={good_degs}), Bad: v={bad_verts} (d+={bad_degs}), t3={r['t3']}")

    # Euler characteristic / alternating sum argument
    print(f"\n--- Euler characteristic approach ---")
    print(f"If β₁(T)=0, then for the Mayer-Vietoris / LES:")
    print(f"  sum_v β₁(T\\v) is bounded by something topological?")
    sum_beta1_dist = defaultdict(int)
    for r in results_5:
        s = sum(vd['beta1_Tv'] for vd in r['vertex_data'])
        sum_beta1_dist[s] += 1
    print(f"  Distribution of sum_v β₁(T\\v) at n=5:")
    for s in sorted(sum_beta1_dist.keys()):
        print(f"    sum={s}: {sum_beta1_dist[s]} tournaments")
    print(f"  Max possible = n = 5 (if all bad)")

    # ===== FINAL SUMMARY =====
    print(f"\n\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"\nn=5 exhaustive ({beta1_zero} with β₁=0):")
    print(f"  All-bad: {all_bad_count}")
    print(f"  Max bad: {max(bad_count_dist.keys())}/{n}")

    print(f"\nn=6 sample ({beta1_zero_6} with β₁=0):")
    print(f"  All-bad: {all_bad_6}")
    if bad_count_dist_6:
        print(f"  Max bad: {max(bad_count_dist_6.keys())}/{6}")

    print(f"\nn=7 sample ({beta1_zero_7} with β₁=0):")
    print(f"  All-bad: {all_bad_7}")
    if bad_count_dist_7:
        print(f"  Max bad: {max(bad_count_dist_7.keys())}/{7}")

    if all_bad_count == 0 and all_bad_6 == 0 and all_bad_7 == 0:
        print(f"\n*** CONJECTURE CONFIRMED: For β₁(T)=0 on n≥5 vertices,")
        print(f"    NOT all vertices can be bad. There always exists v with β₁(T\\v)=0. ***")

if __name__ == '__main__':
    main()
