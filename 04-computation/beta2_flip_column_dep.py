#!/usr/bin/env python3
"""
beta2_flip_column_dep.py — Verify Grok's flip obstruction claim

For tournaments with β₁=0 and exactly 3 bad vertices, flip the transitive
triple's a→c edge to c→a, creating T' with β₁=1.

Claim: rank(∂₂') = rank(∂₂) - 1, i.e., "one minimal column dependence is
added to ∂₂, dropping rank by precisely 1."

We decompose: which columns are removed/added, and how does each step
affect rank?

Author: opus-2026-03-09-S1
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
import sys

# ============================================================
# Core functions (self-contained, adapted from path_homology_v2)
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
    result = []
    for i in range(len(path)):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result


def build_boundary_matrix(paths_p, paths_pm1):
    """∂: columns=p-paths, rows=(p-1)-paths."""
    if not paths_p or not paths_pm1:
        return np.zeros((max(len(paths_pm1), 0), max(len(paths_p), 0)))
    idx = {p: i for i, p in enumerate(paths_pm1)}
    M = np.zeros((len(paths_pm1), len(paths_p)))
    for j, path in enumerate(paths_p):
        for sign, face in boundary_coeffs(path):
            if face in idx:
                M[idx[face], j] += sign
    return M


def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))
    if p == 0:
        return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed:
                    non_allowed[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed:
                P[non_allowed[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    return Vt[rank:].T


def path_betti(A, n, max_dim=None):
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
        dim_op = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_op == 0:
            betti.append(0); continue
        bd_p = build_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_o = bd_p @ omega[p]
        if bd_p_o.shape[0] > 0 and bd_p_o.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_o, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_op - rank_p
        dim_op1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_op1 > 0:
            bd_p1 = build_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_o = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_o, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0
        betti.append(max(0, ker_dim - im_dim))
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
        yield [row[:] for row in A]


def mat_rank(M, tol=1e-8):
    if M.size == 0:
        return 0
    sv = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(sv > tol))


def nullspace(M, tol=1e-8):
    """Return matrix whose columns span ker(M)."""
    if M.size == 0:
        return np.eye(M.shape[1]) if M.shape[1] > 0 else np.zeros((0, 0))
    U, S, Vt = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(S > tol))
    return Vt[rank:].T


# ============================================================
# Identify "bad vertices" and transitive triples
# ============================================================

def find_bad_vertices(A, n):
    """A vertex v is 'bad' if it's NOT on every 3-cycle.
    Actually, for β₁: a vertex is bad if it's a source or sink in some
    local sense. Let's use the standard: v is 'bad' if removing v drops β₁.

    Simpler: find the transitive triple — 3 vertices forming a transitive
    tournament (no 3-cycle among them).
    """
    # Find all transitive triples
    tt = []
    for triple in combinations(range(n), 3):
        i, j, k = triple
        # Check if they form a 3-cycle
        cyc = (A[i][j] and A[j][k] and A[k][i]) or \
              (A[j][i] and A[k][j] and A[i][k])
        if not cyc:
            # Transitive — find the ordering a→b→c, a→c
            verts = list(triple)
            for perm in permutations(verts):
                a, b, c = perm
                if A[a][b] and A[b][c] and A[a][c]:
                    tt.append((a, b, c))
                    break
    return tt


def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3


# ============================================================
# TT (transitive triple) = 3-path in tournament context
# For ∂₂, columns are allowed 2-paths (i.e., paths (x,y,z) with x→y→z)
# ============================================================

def flip_analysis(A, n, a, b, c, verbose=False):
    """
    Flip a→c to c→a in tournament A.
    Analyze ∂₂ column changes and rank changes.

    Returns dict with all analysis results.
    """
    # Build T' (flipped)
    A2 = [row[:] for row in A]
    assert A2[a][c] == 1 and A2[c][a] == 0, f"Expected a→c edge: a={a}, c={c}"
    A2[a][c] = 0
    A2[c][a] = 1

    # Allowed 1-paths (edges) for T and T'
    edges_T = enumerate_allowed_paths(A, n, 1)
    edges_T2 = enumerate_allowed_paths(A2, n, 1)

    # Allowed 2-paths (TTs) for T and T'
    TT_T = enumerate_allowed_paths(A, n, 2)
    TT_T2 = enumerate_allowed_paths(A2, n, 2)

    set_TT = set(TT_T)
    set_TT2 = set(TT_T2)

    removed = sorted(set_TT - set_TT2)
    added = sorted(set_TT2 - set_TT)
    kept = sorted(set_TT & set_TT2)

    # Build ∂₂ for T: rows=edges_T, cols=TT_T
    d2_T = build_boundary_matrix(TT_T, edges_T)

    # Build ∂₂ for T': rows=edges_T2, cols=TT_T2
    d2_T2 = build_boundary_matrix(TT_T2, edges_T2)

    rank_T = mat_rank(d2_T)
    rank_T2 = mat_rank(d2_T2)

    # Now decompose: what happens when we remove R then add A?
    # Use common edge set (all edges except (a,c) which becomes (c,a))
    # Actually edges_T and edges_T2 differ: T has (a,c), T' has (c,a)

    # For the decomposition, work in the FULL edge space (union of both)
    all_edges = sorted(set(edges_T) | set(edges_T2))
    edge_idx = {e: i for i, e in enumerate(all_edges)}
    nrows = len(all_edges)

    # Build columns for each TT in the full edge space
    def tt_column(tt_path, edge_idx_map, nrows):
        col = np.zeros(nrows)
        for sign, face in boundary_coeffs(tt_path):
            if face in edge_idx_map:
                col[edge_idx_map[face]] += sign
        return col

    # Full ∂₂ for T in the combined edge space
    cols_T_full = np.column_stack([tt_column(tt, edge_idx, nrows) for tt in TT_T]) if TT_T else np.zeros((nrows, 0))
    cols_T2_full = np.column_stack([tt_column(tt, edge_idx, nrows) for tt in TT_T2]) if TT_T2 else np.zeros((nrows, 0))

    rank_T_full = mat_rank(cols_T_full)
    rank_T2_full = mat_rank(cols_T2_full)

    # Kept columns only
    kept_cols = np.column_stack([tt_column(tt, edge_idx, nrows) for tt in kept]) if kept else np.zeros((nrows, 0))
    rank_kept = mat_rank(kept_cols)

    # Kept + added
    if added:
        added_cols = np.column_stack([tt_column(tt, edge_idx, nrows) for tt in added])
        kept_plus_added = np.column_stack([kept_cols, added_cols]) if kept else added_cols
    else:
        kept_plus_added = kept_cols
    rank_kept_plus_added = mat_rank(kept_plus_added)

    rank_lost = rank_T_full - rank_kept
    rank_gained = rank_kept_plus_added - rank_kept
    net = rank_gained - rank_lost

    # Nullspace analysis
    ns_T = nullspace(d2_T.T)   # ker(∂₂^T) in edge space — wait, we want ker(∂₂) in TT space
    # ker(∂₂): vectors v s.t. ∂₂ v = 0, where v is over TT columns
    ns_T_tt = nullspace(d2_T)    # columns of this = null vectors in TT-coordinate space
    ns_T2_tt = nullspace(d2_T2)

    dim_ker_T = ns_T_tt.shape[1]
    dim_ker_T2 = ns_T2_tt.shape[1]

    # The "new" null vector: project T' null space onto T's TT basis
    # Find which TT_T2 null vectors are NOT in the span of kept TTs' null space

    result = {
        'n': n, 'a': a, 'b': b, 'c': c,
        'n_TT': len(TT_T), 'n_TT2': len(TT_T2),
        'n_removed': len(removed), 'n_added': len(added), 'n_kept': len(kept),
        'removed': removed, 'added': added,
        'rank_T': rank_T, 'rank_T2': rank_T2,
        'rank_T_full': rank_T_full, 'rank_T2_full': rank_T2_full,
        'rank_kept': rank_kept,
        'rank_kept_plus_added': rank_kept_plus_added,
        'rank_lost': rank_lost, 'rank_gained': rank_gained, 'net': net,
        'dim_ker_T': dim_ker_T, 'dim_ker_T2': dim_ker_T2,
        'n_edges_T': len(edges_T), 'n_edges_T2': len(edges_T2),
    }

    if verbose:
        print(f"  Flip ({a},{b},{c}): a→c becomes c→a")
        print(f"    TTs: {len(TT_T)} → {len(TT_T2)}  (removed {len(removed)}, added {len(added)}, kept {len(kept)})")
        print(f"    Removed TTs: {removed}")
        print(f"    Added TTs:   {added}")
        print(f"    rank(∂₂):  {rank_T} → {rank_T2}  (change: {rank_T2 - rank_T})")
        print(f"    rank(∂₂) [full edge space]: {rank_T_full} → {rank_T2_full}")
        print(f"    rank(kept): {rank_kept}")
        print(f"    rank(kept+added): {rank_kept_plus_added}")
        print(f"    rank_lost={rank_lost}, rank_gained={rank_gained}, net={net}")
        print(f"    dim ker(∂₂): {dim_ker_T} → {dim_ker_T2}")

    return result


# ============================================================
# Ω₂ analysis (the actual boundary map on Ω₂, not raw A₂)
# ============================================================

def omega_rank_analysis(A, n, a, b, c, verbose=False):
    """
    Full Ω-level analysis: compute ∂₂ restricted to Ω₂ for T and T'.
    This is what actually determines β₂.
    """
    A2 = [row[:] for row in A]
    A2[a][c] = 0
    A2[c][a] = 1

    # T
    al1_T = enumerate_allowed_paths(A, n, 1)
    al2_T = enumerate_allowed_paths(A, n, 2)
    al3_T = enumerate_allowed_paths(A, n, 3)
    om2_T = compute_omega_basis(A, n, 2, al2_T, al1_T)
    om3_T = compute_omega_basis(A, n, 3, al3_T, al2_T)
    bd2_T = build_boundary_matrix(al2_T, al1_T)
    bd2_om_T = bd2_T @ om2_T if om2_T.shape[1] > 0 else np.zeros((len(al1_T), 0))
    rank_bd2_om_T = mat_rank(bd2_om_T)

    # T'
    al1_T2 = enumerate_allowed_paths(A2, n, 1)
    al2_T2 = enumerate_allowed_paths(A2, n, 2)
    al3_T2 = enumerate_allowed_paths(A2, n, 3)
    om2_T2 = compute_omega_basis(A2, n, 2, al2_T2, al1_T2)
    om3_T2 = compute_omega_basis(A2, n, 3, al3_T2, al2_T2)
    bd2_T2 = build_boundary_matrix(al2_T2, al1_T2)
    bd2_om_T2 = bd2_T2 @ om2_T2 if om2_T2.shape[1] > 0 else np.zeros((len(al1_T2), 0))
    rank_bd2_om_T2 = mat_rank(bd2_om_T2)

    dim_om2_T = om2_T.shape[1]
    dim_om2_T2 = om2_T2.shape[1]
    dim_om3_T = om3_T.shape[1]
    dim_om3_T2 = om3_T2.shape[1]

    # β₂ = dim(ker ∂₂|Ω₂) - dim(im ∂₃|Ω₃)
    ker2_T = dim_om2_T - rank_bd2_om_T
    ker2_T2 = dim_om2_T2 - rank_bd2_om_T2

    # im ∂₃
    if dim_om3_T > 0:
        bd3_T = build_boundary_matrix(al3_T, al2_T)
        bd3_om_T = bd3_T @ om3_T
        im3_T = mat_rank(bd3_om_T)
    else:
        im3_T = 0
    if dim_om3_T2 > 0:
        bd3_T2 = build_boundary_matrix(al3_T2, al2_T2)
        bd3_om_T2 = bd3_T2 @ om3_T2
        im3_T2 = mat_rank(bd3_om_T2)
    else:
        im3_T2 = 0

    beta2_T = ker2_T - im3_T
    beta2_T2 = ker2_T2 - im3_T2

    result = {
        'dim_om2_T': dim_om2_T, 'dim_om2_T2': dim_om2_T2,
        'rank_bd2_om_T': rank_bd2_om_T, 'rank_bd2_om_T2': rank_bd2_om_T2,
        'ker2_T': ker2_T, 'ker2_T2': ker2_T2,
        'im3_T': im3_T, 'im3_T2': im3_T2,
        'beta2_T': beta2_T, 'beta2_T2': beta2_T2,
        'dim_om3_T': dim_om3_T, 'dim_om3_T2': dim_om3_T2,
    }

    if verbose:
        print(f"    Ω₂: dim {dim_om2_T} → {dim_om2_T2}")
        print(f"    rank(∂₂|Ω₂): {rank_bd2_om_T} → {rank_bd2_om_T2}")
        print(f"    ker(∂₂|Ω₂): {ker2_T} → {ker2_T2}")
        print(f"    Ω₃: dim {dim_om3_T} → {dim_om3_T2}")
        print(f"    im(∂₃|Ω₃): {im3_T} → {im3_T2}")
        print(f"    β₂: {beta2_T} → {beta2_T2}")

    return result


# ============================================================
# Main analysis
# ============================================================

def main():
    print("=" * 70)
    print("FLIP OBSTRUCTION COLUMN DEPENDENCY ANALYSIS")
    print("Grok's claim: flipping a→c in transitive triple drops rank(∂₂) by 1")
    print("=" * 70)

    # ---- n=5: ALL tournaments with β₁=0 and exactly 1 transitive triple ----
    print("\n" + "=" * 70)
    print("n=5: Tournaments with β₁=0, scanning for those with transitive triples")
    print("=" * 70)

    n = 5
    results_n5 = []
    count_total = 0
    count_b1_zero = 0
    count_has_tt = 0

    # Statistics accumulators
    net_changes = defaultdict(int)
    rank_lost_gained = defaultdict(int)

    for A in all_tournaments(n):
        count_total += 1
        betti = path_betti(A, n, max_dim=2)
        if betti[1] != 0:
            continue
        count_b1_zero += 1

        # Find transitive triples
        tts = find_bad_vertices(A, n)
        if not tts:
            continue
        count_has_tt += 1

        # Use the first transitive triple
        a, b, c = tts[0]

        # Verify T' has β₁=1
        A2 = [row[:] for row in A]
        A2[a][c] = 0
        A2[c][a] = 1
        betti2 = path_betti(A2, n, max_dim=2)

        r = flip_analysis(A, n, a, b, c, verbose=(count_has_tt <= 3))
        r_omega = omega_rank_analysis(A, n, a, b, c, verbose=(count_has_tt <= 3))

        r['beta1_T'] = betti[1]
        r['beta1_T2'] = betti2[1]
        r['beta2_T'] = r_omega['beta2_T']
        r['beta2_T2'] = r_omega['beta2_T2']
        r['n_trans_triples'] = len(tts)
        r.update(r_omega)
        results_n5.append(r)

        net_changes[r['net']] += 1
        rank_lost_gained[(r['rank_lost'], r['rank_gained'])] += 1

        if count_has_tt <= 3:
            print(f"\n  Tournament #{count_has_tt}: t3={count_3cycles(A, n)}, "
                  f"β={betti}, β'={betti2}, #TT_triples={len(tts)}")
            print()

    print(f"\n--- n=5 Summary ---")
    print(f"  Total tournaments: {count_total}")
    print(f"  With β₁=0: {count_b1_zero}")
    print(f"  With β₁=0 AND transitive triples: {count_has_tt}")

    print(f"\n  NET rank change distribution (rank(∂₂') - rank(∂₂)):")
    actual_rank_changes = defaultdict(int)
    for r in results_n5:
        actual_rank_changes[r['rank_T2'] - r['rank_T']] += 1
    for k in sorted(actual_rank_changes):
        print(f"    Δrank = {k}: {actual_rank_changes[k]} cases")

    print(f"\n  Decomposed rank change (rank_lost, rank_gained) distribution:")
    for (rl, rg), cnt in sorted(rank_lost_gained.items()):
        print(f"    lost={rl}, gained={rg} (net={rg-rl}): {cnt} cases")

    print(f"\n  Ω₂-level rank changes:")
    omega_rank_changes = defaultdict(int)
    for r in results_n5:
        delta = r['rank_bd2_om_T2'] - r['rank_bd2_om_T']
        omega_rank_changes[delta] += 1
    for k in sorted(omega_rank_changes):
        print(f"    Δrank(∂₂|Ω₂) = {k}: {omega_rank_changes[k]} cases")

    print(f"\n  dim(Ω₂) changes:")
    om2_changes = defaultdict(int)
    for r in results_n5:
        delta = r['dim_om2_T2'] - r['dim_om2_T']
        om2_changes[delta] += 1
    for k in sorted(om2_changes):
        print(f"    Δdim(Ω₂) = {k}: {om2_changes[k]} cases")

    print(f"\n  β₁ changes (T → T'):")
    b1_changes = defaultdict(int)
    for r in results_n5:
        b1_changes[(r['beta1_T'], r['beta1_T2'])] += 1
    for k in sorted(b1_changes):
        print(f"    β₁: {k[0]} → {k[1]}: {b1_changes[k]} cases")

    print(f"\n  β₂ changes (T → T'):")
    b2_changes = defaultdict(int)
    for r in results_n5:
        b2_changes[(r['beta2_T'], r['beta2_T2'])] += 1
    for k in sorted(b2_changes):
        print(f"    β₂: {k[0]} → {k[1]}: {b2_changes[k]} cases")

    # Detailed table for first 20
    print(f"\n  Detailed table (first 20 of {len(results_n5)}):")
    hdr = ("  {:>3} {:>4} {:>4} {:>3} {:>3} "
            "{:>3} {:>3} {:>4} {:>4} {:>4} "
            "{:>3} {:>3} {:>5} {:>5} {:>3} {:>3}").format(
            "#", "TTs", "TTs2", "rm", "add",
            "rk", "rk2", "Drk", "lost", "gain",
            "O2", "O2b", "rk|O", "rk|Ob", "b1b", "b2b")
    print(hdr)
    for i, r in enumerate(results_n5[:20]):
        row = ("  {:3d} {:4d} {:4d} {:3d} {:3d} "
               "{:3d} {:3d} {:4d} {:4d} {:4d} "
               "{:3d} {:3d} {:5d} {:5d} {:3d} {:3d}").format(
               i+1, r['n_TT'], r['n_TT2'], r['n_removed'], r['n_added'],
               r['rank_T'], r['rank_T2'], r['rank_T2']-r['rank_T'],
               r['rank_lost'], r['rank_gained'],
               r['dim_om2_T'], r['dim_om2_T2'],
               r['rank_bd2_om_T'], r['rank_bd2_om_T2'],
               r['beta1_T2'], r['beta2_T2'])
        print(row)

    # ---- Check Grok's specific claim ----
    print(f"\n" + "=" * 70)
    print("GROK'S CLAIM VERIFICATION")
    print("Claim: rank(∂₂') = rank(∂₂) - 1 exactly")
    print("=" * 70)

    claim_holds = sum(1 for r in results_n5 if r['rank_T2'] - r['rank_T'] == -1)
    claim_fails = sum(1 for r in results_n5 if r['rank_T2'] - r['rank_T'] != -1)
    print(f"  n=5: Claim holds for {claim_holds}/{len(results_n5)} cases")
    if claim_fails > 0:
        print(f"  COUNTEREXAMPLES ({claim_fails}):")
        for r in results_n5:
            if r['rank_T2'] - r['rank_T'] != -1:
                print(f"    ({r['a']},{r['b']},{r['c']}): rank {r['rank_T']} → {r['rank_T2']} "
                      f"(Δ={r['rank_T2']-r['rank_T']})")

    # Also check at Ω level
    print(f"\n  At Ω₂ level: rank(∂₂|Ω₂)' = rank(∂₂|Ω₂) - 1 ?")
    om_claim_holds = sum(1 for r in results_n5 if r['rank_bd2_om_T2'] - r['rank_bd2_om_T'] == -1)
    om_claim_fails = len(results_n5) - om_claim_holds
    print(f"  Holds for {om_claim_holds}/{len(results_n5)} cases")

    # ---- Nullspace analysis on first few ----
    print(f"\n" + "=" * 70)
    print("NULLSPACE ANALYSIS (first 5 cases)")
    print("=" * 70)

    case_idx = 0
    for A in all_tournaments(n):
        if case_idx >= 5:
            break
        betti = path_betti(A, n, max_dim=2)
        if betti[1] != 0:
            continue
        tts = find_bad_vertices(A, n)
        if not tts:
            continue

        a, b, c = tts[0]
        A2 = [row[:] for row in A]
        A2[a][c] = 0
        A2[c][a] = 1

        # Get TTs and boundary matrices for T and T'
        edges_T = enumerate_allowed_paths(A, n, 1)
        TT_T = enumerate_allowed_paths(A, n, 2)
        d2_T = build_boundary_matrix(TT_T, edges_T)

        edges_T2 = enumerate_allowed_paths(A2, n, 1)
        TT_T2 = enumerate_allowed_paths(A2, n, 2)
        d2_T2 = build_boundary_matrix(TT_T2, edges_T2)

        ns_T = nullspace(d2_T)
        ns_T2 = nullspace(d2_T2)

        print(f"\n  Case {case_idx+1}: flip ({a},{b},{c})")
        print(f"    dim ker(∂₂) for T:  {ns_T.shape[1]}  (over {len(TT_T)} TTs)")
        print(f"    dim ker(∂₂) for T': {ns_T2.shape[1]}  (over {len(TT_T2)} TTs)")

        # Find the H₁ generator of T'
        # β₁(T') should be 1. The H₁ generator is a 1-cycle in Ω₁.
        al0_T2 = enumerate_allowed_paths(A2, n, 0)
        om1_T2 = compute_omega_basis(A2, n, 1, edges_T2, al0_T2)
        bd1_T2 = build_boundary_matrix(edges_T2, al0_T2)
        bd1_om_T2 = bd1_T2 @ om1_T2 if om1_T2.shape[1] > 0 else np.zeros((len(al0_T2), 0))
        ker_bd1 = nullspace(bd1_om_T2)

        if ker_bd1.shape[1] > 0:
            # Express in edge coordinates
            gamma_om = ker_bd1[:, 0]  # first kernel vector in Ω₁ coords
            gamma_edge = om1_T2 @ gamma_om  # in edge coords

            # Which edges have nonzero coefficient?
            nonzero = [(edges_T2[i], gamma_edge[i]) for i in range(len(edges_T2))
                       if abs(gamma_edge[i]) > 1e-10]
            print(f"    H₁ generator γ of T' (in edge coords):")
            for edge, coeff in nonzero:
                print(f"      {edge}: {coeff:+.4f}")

            # Check: is ∂₂'(any TT) related to γ?
            # For each new TT column, compute its overlap with γ
            if ns_T2.shape[1] > ns_T.shape[1]:
                print(f"    New null vectors in ker(∂₂'):")
                # The new null vectors span... we can look at them
                for j in range(min(3, ns_T2.shape[1])):
                    v = ns_T2[:, j]
                    nonzero_tts = [(TT_T2[i], v[i]) for i in range(len(TT_T2))
                                   if abs(v[i]) > 1e-10]
                    print(f"      Null vector {j}: {len(nonzero_tts)} nonzero entries")
                    for tt, coeff in nonzero_tts[:5]:
                        print(f"        {tt}: {coeff:+.6f}")
                    if len(nonzero_tts) > 5:
                        print(f"        ... ({len(nonzero_tts) - 5} more)")

        case_idx += 1

    # ---- n=6 sample ----
    print(f"\n" + "=" * 70)
    print("n=6: SAMPLE (first 50 with β₁=0 and transitive triples)")
    print("=" * 70)

    n = 6
    results_n6 = []
    count_checked = 0

    import random
    random.seed(42)

    for A in all_tournaments(n):
        if len(results_n6) >= 50:
            break
        betti = path_betti(A, n, max_dim=2)
        if betti[1] != 0:
            continue
        tts = find_bad_vertices(A, n)
        if not tts:
            continue

        a, b, c = tts[0]
        A2 = [row[:] for row in A]
        A2[a][c] = 0
        A2[c][a] = 1
        betti2 = path_betti(A2, n, max_dim=2)

        r = flip_analysis(A, n, a, b, c, verbose=False)
        r_omega = omega_rank_analysis(A, n, a, b, c, verbose=False)
        r['beta1_T'] = betti[1]
        r['beta1_T2'] = betti2[1]
        r.update(r_omega)
        results_n6.append(r)

        count_checked += 1
        if count_checked % 10 == 0:
            print(f"  ... checked {count_checked}", file=sys.stderr)

    print(f"  Collected {len(results_n6)} cases")

    # Summary for n=6
    actual_rank_changes_6 = defaultdict(int)
    for r in results_n6:
        actual_rank_changes_6[r['rank_T2'] - r['rank_T']] += 1
    print(f"\n  Raw ∂₂ rank change distribution:")
    for k in sorted(actual_rank_changes_6):
        print(f"    Δrank = {k}: {actual_rank_changes_6[k]} cases")

    rl_rg_6 = defaultdict(int)
    for r in results_n6:
        rl_rg_6[(r['rank_lost'], r['rank_gained'])] += 1
    print(f"\n  Decomposed (lost, gained):")
    for (rl, rg), cnt in sorted(rl_rg_6.items()):
        print(f"    lost={rl}, gained={rg} (net={rg-rl}): {cnt} cases")

    om_rank_6 = defaultdict(int)
    for r in results_n6:
        om_rank_6[r['rank_bd2_om_T2'] - r['rank_bd2_om_T']] += 1
    print(f"\n  Ω₂-level rank changes:")
    for k in sorted(om_rank_6):
        print(f"    Δrank(∂₂|Ω₂) = {k}: {om_rank_6[k]} cases")

    b1_6 = defaultdict(int)
    for r in results_n6:
        b1_6[(r['beta1_T'], r['beta1_T2'])] += 1
    print(f"\n  β₁ changes:")
    for k in sorted(b1_6):
        print(f"    {k[0]} → {k[1]}: {b1_6[k]} cases")

    b2_6 = defaultdict(int)
    for r in results_n6:
        b2_6[(r['beta2_T'], r['beta2_T2'])] += 1
    print(f"\n  β₂ changes:")
    for k in sorted(b2_6):
        print(f"    {k[0]} → {k[1]}: {b2_6[k]} cases")

    claim_holds_6 = sum(1 for r in results_n6 if r['rank_T2'] - r['rank_T'] == -1)
    print(f"\n  Grok's claim (Δrank=-1): {claim_holds_6}/{len(results_n6)}")

    # ---- Final verdict ----
    print(f"\n" + "=" * 70)
    print("FINAL VERDICT")
    print("=" * 70)

    all_results = results_n5 + results_n6
    claim_ok = sum(1 for r in all_results if r['rank_T2'] - r['rank_T'] == -1)
    print(f"  Total cases tested: {len(all_results)} (n=5: {len(results_n5)}, n=6: {len(results_n6)})")
    print(f"  Grok's claim Δrank(∂₂)=-1: {claim_ok}/{len(all_results)}")

    # More nuanced: what IS the rank change pattern?
    all_deltas = defaultdict(int)
    for r in all_results:
        all_deltas[r['rank_T2'] - r['rank_T']] += 1
    print(f"  Actual Δrank distribution: {dict(sorted(all_deltas.items()))}")

    all_rl_rg = defaultdict(int)
    for r in all_results:
        all_rl_rg[(r['rank_lost'], r['rank_gained'])] += 1
    print(f"  (rank_lost, rank_gained) distribution: {dict(sorted(all_rl_rg.items()))}")


if __name__ == '__main__':
    main()
