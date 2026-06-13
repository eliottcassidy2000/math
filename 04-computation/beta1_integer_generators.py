#!/usr/bin/env python3
"""
beta1_integer_generators.py — Extract INTEGER generators z_v and analyze structure

The SVD-based extraction in the main script gives floating-point generators.
This script uses a different approach: construct Z_1/B_1 explicitly using
integer linear algebra (row reduction over Q).

Key observations from first run:
- z_v ALWAYS uses ALL n-1 vertices of T\v
- Filling chain c_v ALWAYS uses vertex v (100%)
- Bad vertices ALWAYS transitive (0 exceptions through n=7)
- max bad count = 3 (never 4+)
- At n=5: dim_Z₁=C(n-1,2)=6, dim_B₁=5 when β₁=1 (i.e., B₁ drops by 1)
- At n=6: dim_Z₁=C(n-1,2)=10, dim_B₁=9 when β₁=1

NEW APPROACH: Use the "transitive triple" observation as a CONSTRAINT.
For tournaments, Ω₁ = A₁ (all edges are ∂-invariant) so Z₁ = ker(∂₁).
And β₁ = dim(Z₁) - dim(B₁).

The key formula: dim(Z₁) = C(n-1,2) for any tournament on n vertices.
(This is because ∂₁: A₁ → A₀ has rank n-1 for a tournament, and |A₁|=C(n,2).)

So β₁(T) = C(n-1,2) - dim(B₁(T)).

When β₁(T)=0: dim(B₁) = C(n-1,2)
When β₁(T\v)=1: dim(B₁(T\v)) = C(n-2,2) - 1

The question: why does deleting v lose exactly 1 dimension from B₁?

opus-2026-03-08
"""

import numpy as np
from fractions import Fraction
from itertools import combinations
from collections import Counter
import random

# ============================================================
# PATH HOMOLOGY (copied)
# ============================================================

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
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

def build_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1),0), max(len(allowed_p),0)))
    idx = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx: M[idx[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}; na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count; na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed: P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    ns = Vt[rank:].T
    return ns if ns.shape[1] > 0 else np.zeros((dim_Ap, 0))

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None: max_dim = n - 1
    allowed = {p: [] if p < 0 else enumerate_allowed_paths(A, n, p) for p in range(-1, max_dim+2)}
    omega = {p: compute_omega_basis(A, n, p, allowed[p], allowed[p-1]) for p in range(max_dim+2)}
    betti = []
    for p in range(max_dim + 1):
        dim_om = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_om == 0: betti.append(0); continue
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_om = bd_p @ omega[p]
        if bd_p_om.size > 0:
            sv = np.linalg.svd(bd_p_om, compute_uv=False)
            rk = sum(s > 1e-8 for s in sv)
        else: rk = 0
        ker = dim_om - rk
        dim_om1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_om1 > 0:
            bd1 = build_boundary_matrix(allowed[p+1], allowed[p])
            bd1_om = bd1 @ omega[p+1]
            sv1 = np.linalg.svd(bd1_om, compute_uv=False)
            im = sum(s > 1e-8 for s in sv1)
        else: im = 0
        betti.append(max(0, ker - im))
    return betti


def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def delete_vertex(A, n, v):
    new_n = n - 1
    B = [[0]*new_n for _ in range(new_n)]
    o2n = {}; j = 0
    for i in range(n):
        if i == v: continue
        o2n[i] = j; j += 1
    for i in range(n):
        if i == v: continue
        for k in range(n):
            if k == v: continue
            B[o2n[i]][o2n[k]] = A[i][k]
    return B, new_n, o2n


# ============================================================
# KEY ANALYSIS: The B₁ dimension drop mechanism
# ============================================================

def analyze_B1_drop(A, n, verbose=True):
    """For a tournament T on n vertices with β₁(T)=0,
    analyze which 2-chains get lost when we delete each vertex.

    Key insight: dim(B₁(T)) = C(n-1,2) = dim(Z₁(T)) when β₁=0.
    When we delete v, dim(B₁(T\\v)) = C(n-2,2)-1 (drops by 1 from max).

    The allowed 2-paths in T\\v are a SUBSET of those in T.
    So Ω₂(T\\v) ⊆ Ω₂(T) (after relabeling).
    But B₁(T\\v) = im(∂₂|_{Ω₂(T\\v)}).

    The 2-chains LOST are those in Ω₂(T) that use paths THROUGH v.
    """

    # Compute Ω₂(T) and its boundary image
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a0 = enumerate_allowed_paths(A, n, 0)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    bd2 = build_boundary_matrix(a2, a1)

    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    bd2_om = bd2 @ om2 if dim_om2 > 0 else np.zeros((len(a1), 0))

    # rank of ∂₂ on Ω₂ = dim(B₁(T))
    if bd2_om.size > 0:
        sv = np.linalg.svd(bd2_om, compute_uv=False)
        rank_B1_T = sum(s > 1e-8 for s in sv)
    else:
        rank_B1_T = 0

    if verbose:
        print(f"  T: |A₁|={len(a1)}, |Ω₂|={dim_om2}, dim(B₁)={rank_B1_T}, "
              f"C(n-1,2)={n*(n-1)//2 - (n-1)}")

    # For each vertex v, analyze what happens when we delete it
    idx_a2 = {p: i for i, p in enumerate(a2)}

    results = {}
    for v in range(n):
        # Which 2-paths go through v?
        through_v = [i for i, p in enumerate(a2) if v in p]
        not_through_v = [i for i, p in enumerate(a2) if v not in p]

        # Which Ω₂ basis vectors have support on paths through v?
        if dim_om2 > 0:
            # Project Ω₂ onto "through v" vs "not through v"
            om2_thru = om2[through_v, :] if through_v else np.zeros((0, dim_om2))
            om2_rest = om2[not_through_v, :] if not_through_v else np.zeros((0, dim_om2))

            # Basis vectors that are ENTIRELY supported on non-v paths
            # = those with om2_thru[:, j] == 0
            pure_non_v = []
            for j in range(dim_om2):
                if np.linalg.norm(om2_thru[:, j]) < 1e-10:
                    pure_non_v.append(j)
        else:
            pure_non_v = []

        # Compute β₁(T\v) and dim(B₁(T\v))
        Bv, mv, o2n = delete_vertex(A, n, v)
        bv = path_betti_numbers(Bv, mv, max_dim=1)
        is_bad = (bv[1] == 1)

        # Get dim(B₁(T\v)) directly
        a1v = enumerate_allowed_paths(Bv, mv, 1)
        a2v = enumerate_allowed_paths(Bv, mv, 2)
        a0v = enumerate_allowed_paths(Bv, mv, 0)
        om2v = compute_omega_basis(Bv, mv, 2, a2v, a1v)
        bd2v = build_boundary_matrix(a2v, a1v)
        dim_om2v = om2v.shape[1] if om2v.ndim == 2 else 0
        bd2v_om = bd2v @ om2v if dim_om2v > 0 else np.zeros((len(a1v), 0))
        if bd2v_om.size > 0:
            svv = np.linalg.svd(bd2v_om, compute_uv=False)
            rank_B1v = sum(s > 1e-8 for s in svv)
        else:
            rank_B1v = 0

        dim_Z1v = (mv*(mv-1))//2 - (mv - 1)  # C(mv-1,2) for tournament

        results[v] = {
            'bad': is_bad,
            'beta1': bv[1],
            'dim_B1': rank_B1v,
            'dim_Z1': dim_Z1v,
            'dim_Omega2': dim_om2v,
            'through_v_paths': len(through_v),
            'not_through_v_paths': len(not_through_v),
            'pure_non_v_omega': len(pure_non_v),
        }

        if verbose:
            status = "BAD" if is_bad else "ok"
            print(f"  v={v} [{status}]: β₁={bv[1]}, dim(B₁)={rank_B1v}/{dim_Z1v}, "
                  f"|Ω₂(T\\v)|={dim_om2v}, "
                  f"2-paths thru v: {len(through_v)}, "
                  f"pure non-v Ω₂: {len(pure_non_v)}")

    return results


def analyze_transitive_triple_structure(A, n, bad, verbose=True):
    """Given three bad vertices forming a transitive triple,
    analyze the structural reason for transitivity.

    Key question: If bad={a,b,c} with a→b→c→(and c→a would make a 3-cycle,
    but instead a→c making it transitive), what is special about this ordering?
    """

    if len(bad) != 3:
        return

    a, b, c = bad
    # Find the transitive ordering
    # In a transitive triple, one vertex beats both others ("source"),
    # one loses to both ("sink"), one beats one and loses to one ("middle")
    scores = {}
    for v in bad:
        s = sum(1 for u in bad if u != v and A[v][u])
        scores[v] = s

    source = [v for v in bad if scores[v] == 2][0]
    sink = [v for v in bad if scores[v] == 0][0]
    mid = [v for v in bad if scores[v] == 1][0]

    if verbose:
        print(f"\n  Transitive ordering: {source} → {mid} → {sink}")
        print(f"    (source={source} beats both, sink={sink} loses to both)")

    # Score of each bad vertex in the FULL tournament
    for v in bad:
        full_score = sum(A[v][u] for u in range(n) if u != v)
        non_bad_score = sum(A[v][u] for u in range(n) if u != v and u not in bad)
        if verbose:
            print(f"    Score of {v}: {full_score} total, {non_bad_score} vs non-bad")

    # KEY: Count 3-cycles involving each bad vertex with non-bad vertices
    for v in bad:
        t3_with_v = 0
        for i in range(n):
            if i == v or i in bad:
                continue
            for j in range(i+1, n):
                if j == v or j in bad:
                    continue
                triple = [v, i, j]
                # Check for 3-cycle
                if (A[v][i] and A[i][j] and A[j][v]) or (A[i][v] and A[j][i] and A[v][j]):
                    t3_with_v += 1
        if verbose:
            print(f"    3-cycles involving {v} (with non-bad): {t3_with_v}")

    # Cross-edges: between bad and non-bad
    non_bad = [u for u in range(n) if u not in bad]
    if verbose:
        print(f"\n  Cross-edge pattern (bad→non-bad):")
        for v in [source, mid, sink]:
            pattern = ""
            for u in non_bad:
                if A[v][u]:
                    pattern += f"+{u} "
                else:
                    pattern += f"-{u} "
            print(f"    {v}: {pattern}")


def focused_dimension_analysis(n, max_examples=20):
    """Focus on the dimension drop mechanism."""
    print(f"\n{'#'*70}")
    print(f"# DIMENSION DROP ANALYSIS: n={n}")
    print(f"{'#'*70}")

    count = 0
    for A in all_tournaments(n):
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
        if count > max_examples:
            break

        print(f"\n--- Example {count}: bad = {bad} ---")
        results = analyze_B1_drop(A, n, verbose=True)
        analyze_transitive_triple_structure(A, n, bad, verbose=True)

    print(f"\n  (Analyzed {count} examples)")


def aggregate_dimension_analysis(n):
    """Collect aggregate statistics on dim(Ω₂) and dim(B₁) patterns."""
    print(f"\n{'#'*70}")
    print(f"# AGGREGATE DIMENSION ANALYSIS: n={n}")
    print(f"{'#'*70}")

    # For bad vs non-bad vertices, what is dim(Ω₂(T\v))?
    om2_bad = []
    om2_nonbad = []
    B1_bad = []
    B1_nonbad = []

    for A in all_tournaments(n):
        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue

        bad = set()
        for v in range(n):
            B, m, _ = delete_vertex(A, n, v)
            bv = path_betti_numbers(B, m, max_dim=1)
            if bv[1] == 1:
                bad.add(v)

        if not bad:
            continue

        for v in range(n):
            Bv, mv, _ = delete_vertex(A, n, v)
            a1v = enumerate_allowed_paths(Bv, mv, 1)
            a2v = enumerate_allowed_paths(Bv, mv, 2)
            om2v = compute_omega_basis(Bv, mv, 2, a2v, a1v)
            dim_om2v = om2v.shape[1] if om2v.ndim == 2 else 0

            bd2v = build_boundary_matrix(a2v, a1v)
            bd2v_om = bd2v @ om2v if dim_om2v > 0 else np.zeros((len(a1v), 0))
            if bd2v_om.size > 0:
                svv = np.linalg.svd(bd2v_om, compute_uv=False)
                rk = sum(s > 1e-8 for s in svv)
            else:
                rk = 0

            if v in bad:
                om2_bad.append(dim_om2v)
                B1_bad.append(rk)
            else:
                om2_nonbad.append(dim_om2v)
                B1_nonbad.append(rk)

    print(f"\nBad vertices (β₁(T\\v)=1):")
    print(f"  dim(Ω₂) distribution: {Counter(om2_bad)}")
    print(f"  dim(B₁) distribution: {Counter(B1_bad)}")

    print(f"\nNon-bad vertices (β₁(T\\v)=0):")
    print(f"  dim(Ω₂) distribution: {Counter(om2_nonbad)}")
    print(f"  dim(B₁) distribution: {Counter(B1_nonbad)}")

    expected_Z1 = (n-2)*(n-3)//2  # C(n-2,2) for T\v
    print(f"\n  expected dim(Z₁) = C({n-2},2) = {expected_Z1}")
    print(f"  Bad: dim(B₁) should be {expected_Z1 - 1} (one less)")
    print(f"  Non-bad: dim(B₁) should be {expected_Z1} (full)")


def score_sequence_analysis(n):
    """Analyze score sequences of bad vertices vs the full tournament."""
    print(f"\n{'#'*70}")
    print(f"# SCORE SEQUENCE ANALYSIS: n={n}")
    print(f"{'#'*70}")

    bad_scores = Counter()
    bad_role = Counter()  # source/mid/sink score among bad triple

    for A in all_tournaments(n):
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

        # Score of each bad vertex
        for v in bad:
            s = sum(A[v][u] for u in range(n) if u != v)
            bad_scores[s] += 1

            # Role in triple
            s_in_triple = sum(A[v][u] for u in bad if u != v)
            if s_in_triple == 2:
                bad_role[('source', s)] += 1
            elif s_in_triple == 0:
                bad_role[('sink', s)] += 1
            else:
                bad_role[('mid', s)] += 1

    print(f"\nScore distribution of bad vertices: {dict(sorted(bad_scores.items()))}")
    print(f"\nRole × full score:")
    for key in sorted(bad_role.keys()):
        print(f"  {key}: {bad_role[key]}")


# ============================================================
# THE CRITICAL TEST: Can a 3-cycle among bad vertices exist?
# ============================================================

def three_cycle_obstruction(n, num_samples=None):
    """Try to construct a tournament where the bad set forms a 3-cycle.
    If impossible, the reason IS the transitivity constraint.

    Approach: take a known tournament with 3 bad vertices {a,b,c} in
    transitive order a→b→c, a→c. What if we flip the a→c edge to c→a?
    Does β₁(T') become nonzero?
    """
    print(f"\n{'#'*70}")
    print(f"# 3-CYCLE OBSTRUCTION TEST: n={n}")
    print(f"{'#'*70}")

    count = 0
    flipped_still_beta0 = 0
    flipped_creates_beta1 = 0
    flipped_bad_count = Counter()

    gen = all_tournaments(n) if n <= 6 else None

    if gen is None:
        # Sample
        tested = 0
        while count < (num_samples or 200) and tested < 50000:
            tested += 1
            A = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5: A[i][j] = 1
                    else: A[j][i] = 1

            betti = path_betti_numbers(A, n, max_dim=1)
            if betti[1] != 0: continue

            bad = []
            for v in range(n):
                B, m, _ = delete_vertex(A, n, v)
                bv = path_betti_numbers(B, m, max_dim=1)
                if bv[1] == 1: bad.append(v)

            if len(bad) != 3: continue
            count += 1
            _do_flip_test(A, n, bad, count <= 5)
    else:
        for A in gen:
            betti = path_betti_numbers(A, n, max_dim=1)
            if betti[1] != 0: continue

            bad = []
            for v in range(n):
                B, m, _ = delete_vertex(A, n, v)
                bv = path_betti_numbers(B, m, max_dim=1)
                if bv[1] == 1: bad.append(v)

            if len(bad) != 3: continue
            count += 1

            result = _do_flip_test(A, n, bad, count <= 10)
            if result == 'still_beta0':
                flipped_still_beta0 += 1
            else:
                flipped_creates_beta1 += 1

    print(f"\n  Total with 3 bad: {count}")
    print(f"  After flipping source→sink to create 3-cycle:")
    print(f"    Still β₁=0: {flipped_still_beta0}")
    print(f"    Now β₁≥1: {flipped_creates_beta1}")


def _do_flip_test(A, n, bad, verbose):
    """Flip the source→sink edge among the bad triple to create a 3-cycle."""
    a, b, c = bad
    scores = {v: sum(A[v][u] for u in bad if u != v) for v in bad}
    source = [v for v in bad if scores[v] == 2][0]
    sink = [v for v in bad if scores[v] == 0][0]
    mid = [v for v in bad if scores[v] == 1][0]

    # Flip source→sink
    A2 = [row[:] for row in A]
    A2[source][sink] = 0
    A2[sink][source] = 1

    # Now bad triple has: source→mid, mid→sink, sink→source = 3-cycle!
    betti2 = path_betti_numbers(A2, n, max_dim=1)

    if verbose:
        print(f"\n  bad={bad}, source={source}, mid={mid}, sink={sink}")
        print(f"  Original: β₁(T)=0")
        print(f"  Flipped {source}→{sink}: β₁(T')={betti2[1]}")

        if betti2[1] == 0:
            # Check bad vertices of T'
            bad2 = []
            for v in range(n):
                B, m, _ = delete_vertex(A2, n, v)
                bv = path_betti_numbers(B, m, max_dim=1)
                if bv[1] == 1: bad2.append(v)
            print(f"  Bad vertices of T': {bad2}")

    return 'still_beta0' if betti2[1] == 0 else 'creates_beta1'


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    random.seed(42)
    np.set_printoptions(precision=4, suppress=True, linewidth=120)

    print("β₁ TRANSITIVE MECHANISM — PART 2: INTEGER GENERATORS & DIMENSION ANALYSIS")
    print("=" * 70)

    # n=5: focused dimension analysis
    focused_dimension_analysis(5, max_examples=5)

    # n=5: aggregate dimensions
    aggregate_dimension_analysis(5)

    # n=5: score analysis
    score_sequence_analysis(5)

    # n=5: flip test
    three_cycle_obstruction(5)

    # n=6: aggregate
    aggregate_dimension_analysis(6)
    score_sequence_analysis(6)
    three_cycle_obstruction(6)

    print("\n\nDone.")
