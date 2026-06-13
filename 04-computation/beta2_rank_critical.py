#!/usr/bin/env python3
"""
Rank-critical transitive triples in tournament path homology.

A transitive triple (TT) τ = (a,b,c) with a→b, b→c, a→c is an element of Ω₂.
Its boundary ∂₂(τ) = (b,c) - (a,c) + (a,b).

A TT is "rank-critical" if removing it from Ω₂ drops rank(∂₂) by 1,
i.e., its boundary is NOT in the span of boundaries of all other TTs.

Analysis for all n=5 (1024) and n=6 (32768) tournaments.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
import sys
import time

# ─── Core tournament enumeration ───

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

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

# ─── Allowed paths and boundary ───

def enumerate_allowed_paths(A, n, p):
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

def get_transitive_triples(A, n):
    """Return all TTs (a,b,c) with a→b, b→c, a→c."""
    tts = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]:
                continue
            for c in range(n):
                if c == a or c == b:
                    continue
                if A[b][c] and A[a][c]:
                    tts.append((a, b, c))
    return tts

def build_boundary_matrix_for_tts(tts, allowed_1paths):
    """Build ∂₂ restricted to given TTs.

    ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)

    Rows = allowed 1-paths (edges), Columns = TTs.
    """
    if not tts or not allowed_1paths:
        return np.zeros((max(len(allowed_1paths), 0), max(len(tts), 0)))

    idx_1 = {path: i for i, path in enumerate(allowed_1paths)}
    M = np.zeros((len(allowed_1paths), len(tts)))

    for j, (a, b, c) in enumerate(tts):
        # ∂(a,b,c) = (b,c) - (a,c) + (a,b)
        if (b, c) in idx_1:
            M[idx_1[(b, c)], j] += 1
        if (a, c) in idx_1:
            M[idx_1[(a, c)], j] -= 1
        if (a, b) in idx_1:
            M[idx_1[(a, b)], j] += 1

    return M

# ─── β₁ computation ───

def compute_beta1(A, n):
    """Compute β₁ for tournament."""
    allowed_0 = enumerate_allowed_paths(A, n, 0)
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)

    # Ω₁ = A₁ for tournaments (all faces of edges are vertices, always allowed)
    # Actually need to check: ∂(a,b) = (b) - (a), both are in A₀. So Ω₁ = A₁.

    # ∂₁: A₁ → A₀
    idx_0 = {v: i for i, v in enumerate(allowed_0)}
    bd1 = np.zeros((len(allowed_0), len(allowed_1)))
    for j, (a, b) in enumerate(allowed_1):
        bd1[idx_0[(b,)], j] += 1
        bd1[idx_0[(a,)], j] -= 1

    # ker(∂₁) dimension
    if bd1.shape[1] == 0:
        return 0
    S1 = np.linalg.svd(bd1, compute_uv=False)
    rank1 = int(np.sum(S1 > 1e-8))
    ker1_dim = len(allowed_1) - rank1

    # Ω₂ and ∂₂
    # For tournaments, need to find Ω₂
    allowed_1_set = set(allowed_1)

    # Check which 2-paths have all faces allowed
    # ∂(a,b,c) = (b,c) - (a,c) + (a,b)
    # Need (b,c), (a,c), (a,b) all in A₁
    omega2_paths = []
    for path in allowed_2:
        a, b, c = path
        if (b, c) in allowed_1_set and (a, c) in allowed_1_set and (a, b) in allowed_1_set:
            omega2_paths.append(path)

    # ∂₂: Ω₂ → A₁
    bd2 = build_boundary_matrix_for_tts(omega2_paths, allowed_1)

    if bd2.shape[1] == 0:
        im2_dim = 0
    else:
        S2 = np.linalg.svd(bd2, compute_uv=False)
        im2_dim = int(np.sum(S2 > 1e-8))

    beta1 = ker1_dim - im2_dim
    return max(0, beta1)

# ─── "Bad vertices" ───

def get_bad_vertices(A, n):
    """A vertex v is 'bad' if it's NOT the endpoint of a transitive triple as middle vertex
    that provides coverage... Actually, let's define bad vertices more carefully.

    In tournament path homology, a 'bad face' of a 2-path (a,b,c) is a face
    that is NOT an allowed 1-path. The face (a,c) is not allowed iff there's
    no edge a→c. But in a tournament, either a→c or c→a, so (a,c) is allowed
    iff A[a][c]=1.

    A 2-path (a,b,c) is in Ω₂ iff ALL its faces are allowed 1-paths.
    Faces: (b,c), (a,c), (a,b). In a tournament, (b,c) is allowed iff b→c,
    (a,b) is allowed iff a→b. These are given since (a,b,c) is allowed (a→b→c).
    The only question is (a,c): allowed iff a→c.

    So Ω₂ = { allowed 2-paths (a,b,c) with a→c } = transitive triples!

    A 'bad' 2-path is one NOT in Ω₂, i.e., a→b→c but c→a (a 3-cycle).

    Let's define 'bad vertices' differently. A vertex v is bad if v participates
    in at least one 3-cycle. In that case, some 2-paths through v have bad faces.

    Actually, let me count bad vertices as vertices that are part of 3-cycles.
    """
    bad = set()
    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][j] and A[j][k] and A[k][i]:
                    bad.add(i)
                    bad.add(j)
                    bad.add(k)
    return bad

# ─── Rank-critical analysis ───

def analyze_rank_critical(A, n):
    """For a tournament, find all rank-critical TTs.

    Returns dict with:
    - tts: list of TTs
    - rank_full: rank(∂₂) with all TTs
    - rank_critical: list of TTs that are rank-critical
    - beta1: β₁ of the tournament
    """
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    tts = get_transitive_triples(A, n)

    if not tts:
        return {
            'tts': tts,
            'rank_full': 0,
            'rank_critical': [],
            'beta1': compute_beta1(A, n),
            'n_tts': 0,
            'n_rc': 0,
        }

    # Full boundary matrix
    bd = build_boundary_matrix_for_tts(tts, allowed_1)
    S_full = np.linalg.svd(bd, compute_uv=False)
    rank_full = int(np.sum(S_full > 1e-8))

    # Check each TT for rank-criticality
    rank_critical = []
    for i in range(len(tts)):
        # Remove column i
        cols = list(range(len(tts)))
        cols.pop(i)
        if len(cols) == 0:
            rank_without = 0
        else:
            bd_reduced = bd[:, cols]
            S_red = np.linalg.svd(bd_reduced, compute_uv=False)
            rank_without = int(np.sum(S_red > 1e-8))

        if rank_without < rank_full:
            rank_critical.append(tts[i])

    beta1 = compute_beta1(A, n)

    return {
        'tts': tts,
        'rank_full': rank_full,
        'rank_critical': rank_critical,
        'beta1': beta1,
        'n_tts': len(tts),
        'n_rc': len(rank_critical),
    }

def analyze_missing_tts_beta1(A, n):
    """For β₁>0 tournaments: which NON-TTs, if they were TTs, would increase rank(∂₂)?

    These are 2-paths (a,b,c) with a→b→c but c→a (3-cycle paths).
    If we flipped c→a to a→c, they'd become TTs. Which of these would fill homology?
    """
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    tts = get_transitive_triples(A, n)

    # Current boundary matrix
    bd = build_boundary_matrix_for_tts(tts, allowed_1)
    S_full = np.linalg.svd(bd, compute_uv=False)
    rank_full = int(np.sum(S_full > 1e-8))

    # Non-TT allowed 2-paths (3-cycle paths)
    tt_set = set(tts)
    non_tts = [p for p in allowed_2 if p not in tt_set]

    # For each non-TT, check if adding it (as if it were a TT) would increase rank
    filling = []
    for path in non_tts:
        a, b, c = path
        # Its hypothetical boundary: (b,c) - (a,c) + (a,b)
        # But (a,c) is NOT an allowed 1-path (since c→a, not a→c)
        # However, what if we imagine flipping that edge?
        # Actually the question is: what TTs could be created by flipping arcs?
        # If we flip c→a to a→c, then (a,b,c) becomes a TT.
        # Let's see what that does to the boundary matrix.
        # New TT boundary uses allowed edges of the NEW tournament.
        # This is complex. Let's instead ask: which edges, if flipped, would reduce β₁?
        pass

    # Simpler: for each 3-cycle {a,b,c}, try flipping each arc and see if β₁ drops
    cycles_3 = set()
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    cycles_3.add((i, j, k))

    # For each edge in each 3-cycle, try flipping and check β₁
    flip_results = []
    edges_tried = set()
    for cyc in cycles_3:
        for i_idx in range(3):
            for j_idx in range(3):
                if i_idx == j_idx:
                    continue
                a, b = cyc[i_idx], cyc[j_idx]
                if A[a][b] and (a, b) not in edges_tried:
                    edges_tried.add((a, b))
                    # Flip a→b to b→a
                    A2 = [row[:] for row in A]
                    A2[a][b] = 0
                    A2[b][a] = 1
                    new_beta1 = compute_beta1(A2, n)
                    if new_beta1 < compute_beta1(A, n):
                        flip_results.append(((a, b), new_beta1))

    return flip_results

# ═══════════════════════════════════════════════════════════════
# MAIN ANALYSIS
# ═══════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("RANK-CRITICAL TRANSITIVE TRIPLES IN TOURNAMENT PATH HOMOLOGY")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        t0 = time.time()

        # Collect data
        all_data = []
        beta1_dist = Counter()

        count = 0
        for A in all_tournaments(n):
            result = analyze_rank_critical(A, n)
            t3 = count_3cycles(A, n)
            bad_verts = get_bad_vertices(A, n)
            result['t3'] = t3
            result['n_bad'] = len(bad_verts)
            result['bad_verts'] = bad_verts
            all_data.append(result)
            beta1_dist[result['beta1']] += 1
            count += 1
            if count % 1000 == 0:
                print(f"  ... processed {count} tournaments ({time.time()-t0:.1f}s)")

        elapsed = time.time() - t0
        print(f"\nProcessed {count} tournaments in {elapsed:.1f}s")
        print(f"\nβ₁ distribution: {dict(sorted(beta1_dist.items()))}")

        # ─── Section 1: Overall statistics ───
        print(f"\n--- Section 1: Overall rank-critical statistics ---")

        rc_counts = [d['n_rc'] for d in all_data]
        tt_counts = [d['n_tts'] for d in all_data]
        print(f"  #TTs range: [{min(tt_counts)}, {max(tt_counts)}]")
        print(f"  #rank-critical range: [{min(rc_counts)}, {max(rc_counts)}]")
        print(f"  Mean #rank-critical: {np.mean(rc_counts):.2f}")
        print(f"  Mean fraction rank-critical: {np.mean([d['n_rc']/d['n_tts'] if d['n_tts']>0 else 0 for d in all_data]):.4f}")

        # Distribution of #rank-critical
        rc_dist = Counter(rc_counts)
        print(f"\n  Distribution of #rank-critical TTs:")
        for k in sorted(rc_dist.keys()):
            print(f"    {k} rank-critical: {rc_dist[k]} tournaments")

        # ─── Section 2: β₁=0, grouped by #bad vertices ───
        print(f"\n--- Section 2: β₁=0 tournaments, by #bad vertices ---")

        beta0_data = [d for d in all_data if d['beta1'] == 0]
        bad_groups = defaultdict(list)
        for d in beta0_data:
            bad_groups[d['n_bad']].append(d)

        for nb in sorted(bad_groups.keys()):
            group = bad_groups[nb]
            rcs = [d['n_rc'] for d in group]
            tts = [d['n_tts'] for d in group]
            fracs = [d['n_rc']/d['n_tts'] if d['n_tts']>0 else 0 for d in group]
            print(f"\n  #bad_vertices = {nb}: {len(group)} tournaments")
            print(f"    #TTs: mean={np.mean(tts):.1f}, range=[{min(tts)},{max(tts)}]")
            print(f"    #rank-critical: mean={np.mean(rcs):.2f}, range=[{min(rcs)},{max(rcs)}]")
            print(f"    fraction RC: mean={np.mean(fracs):.4f}")
            rc_sub = Counter(rcs)
            print(f"    RC distribution: {dict(sorted(rc_sub.items()))}")

        # ─── Section 2b: β₁=0, specifically 3 bad vertices ───
        if 3 in bad_groups:
            print(f"\n--- Section 2b: β₁=0, exactly 3 bad vertices (one 3-cycle) ---")
            group3 = bad_groups[3]
            # Check: is the rank = n(n-1)/2 - n + 1 - 0 = C(n,2) - n + 1?
            # (formula from memory: rank(∂₂|Ω₂) = C(n,2) - n + 1 - β₁)
            expected_rank = n*(n-1)//2 - n + 1
            for d in group3[:5]:
                print(f"    rank_full={d['rank_full']} (expected {expected_rank}), "
                      f"#TTs={d['n_tts']}, #RC={d['n_rc']}")

        # ─── Section 2c: β₁=0, 0 bad vertices (transitive) ───
        if 0 in bad_groups:
            print(f"\n--- Section 2c: β₁=0, 0 bad vertices (transitive tournament) ---")
            group0 = bad_groups[0]
            for d in group0[:5]:
                print(f"    rank_full={d['rank_full']}, #TTs={d['n_tts']}, #RC={d['n_rc']}")
                if d['rank_critical']:
                    print(f"      RC TTs: {d['rank_critical'][:10]}")

        # ─── Section 3: Correlation #RC vs #bad vertices ───
        print(f"\n--- Section 3: Correlation between #RC and #bad_vertices ---")
        rcs_all = np.array([d['n_rc'] for d in all_data], dtype=float)
        bads_all = np.array([d['n_bad'] for d in all_data], dtype=float)
        t3s_all = np.array([d['t3'] for d in all_data], dtype=float)
        tts_all = np.array([d['n_tts'] for d in all_data], dtype=float)

        if np.std(rcs_all) > 0 and np.std(bads_all) > 0:
            corr_rc_bad = np.corrcoef(rcs_all, bads_all)[0, 1]
            print(f"  corr(#RC, #bad_vertices) = {corr_rc_bad:.4f}")

        if np.std(rcs_all) > 0 and np.std(t3s_all) > 0:
            corr_rc_t3 = np.corrcoef(rcs_all, t3s_all)[0, 1]
            print(f"  corr(#RC, t3) = {corr_rc_t3:.4f}")

        if np.std(rcs_all) > 0 and np.std(tts_all) > 0:
            corr_rc_tts = np.corrcoef(rcs_all, tts_all)[0, 1]
            print(f"  corr(#RC, #TTs) = {corr_rc_tts:.4f}")

        # ─── Section 4: rank_full = ? ───
        print(f"\n--- Section 4: rank(∂₂|Ω₂) statistics ---")
        ranks = [d['rank_full'] for d in all_data]
        rank_dist = Counter(ranks)
        for r in sorted(rank_dist.keys()):
            # Cross with β₁
            b1_for_rank = Counter(d['beta1'] for d in all_data if d['rank_full'] == r)
            print(f"  rank={r}: {rank_dist[r]} tournaments, β₁ dist: {dict(b1_for_rank)}")

        # Verify formula: rank = C(n,2) - n + 1 - β₁
        print(f"\n  Verify: rank = C(n,2) - n + 1 - β₁ = {n*(n-1)//2} - {n} + 1 - β₁")
        violations = 0
        for d in all_data:
            expected = n*(n-1)//2 - n + 1 - d['beta1']
            if d['rank_full'] != expected:
                violations += 1
                if violations <= 3:
                    print(f"    VIOLATION: rank={d['rank_full']}, expected={expected}, β₁={d['beta1']}")
        print(f"  Violations: {violations}/{count}")

        # ─── Section 5: β₁>0 tournaments — which flips reduce β₁? ───
        print(f"\n--- Section 5: β₁>0 tournaments ---")
        beta1_pos = [d for d in all_data if d['beta1'] > 0]
        print(f"  {len(beta1_pos)} tournaments with β₁>0")

        if beta1_pos and n <= 5:
            # Detailed analysis for small n
            for idx, d in enumerate(beta1_pos[:5]):
                # Reconstruct adjacency
                pass  # We'd need to store A; let's do it differently

        # For β₁>0, what are the RC stats?
        if beta1_pos:
            rcs_b1 = [d['n_rc'] for d in beta1_pos]
            print(f"  #RC for β₁>0: mean={np.mean(rcs_b1):.2f}, range=[{min(rcs_b1)},{max(rcs_b1)}]")
            rc_b1_dist = Counter(rcs_b1)
            print(f"  RC distribution: {dict(sorted(rc_b1_dist.items()))}")

        # ─── Section 5b: β₁>0 — which TTs, if added, would fill homology? ───
        if n <= 5:
            print(f"\n--- Section 5b: β₁>0 — arc flips that reduce β₁ ---")
            flip_count = 0
            for A in all_tournaments(n):
                beta1 = compute_beta1(A, n)
                if beta1 == 0:
                    continue
                flips = analyze_missing_tts_beta1(A, n)
                if flips and flip_count < 5:
                    t3 = count_3cycles(A, n)
                    print(f"  Tournament with t3={t3}, β₁={beta1}:")
                    print(f"    Flips reducing β₁: {flips}")
                    flip_count += 1

        # ─── Section 6: rank-critical vs rank-full relationship ───
        print(f"\n--- Section 6: #RC as function of rank and #TTs ---")

        # Group by (rank, #TTs)
        rt_groups = defaultdict(list)
        for d in all_data:
            rt_groups[(d['rank_full'], d['n_tts'])].append(d['n_rc'])

        print(f"  (rank, #TTs) → #RC distribution:")
        for key in sorted(rt_groups.keys()):
            vals = rt_groups[key]
            c = Counter(vals)
            if len(vals) <= 20:
                print(f"    rank={key[0]}, #TTs={key[1]}: n={len(vals)}, #RC dist={dict(sorted(c.items()))}")
            else:
                print(f"    rank={key[0]}, #TTs={key[1]}: n={len(vals)}, #RC mean={np.mean(vals):.2f}, range=[{min(vals)},{max(vals)}]")

        # ─── KEY QUESTION: Is #RC always = rank? ───
        print(f"\n--- KEY: Is #rank-critical always equal to rank(∂₂)? ---")
        matches = sum(1 for d in all_data if d['n_rc'] == d['rank_full'])
        print(f"  #RC == rank: {matches}/{count} ({100*matches/count:.1f}%)")

        mismatches_less = [(d['n_rc'], d['rank_full'], d['n_tts'], d['beta1'])
                           for d in all_data if d['n_rc'] < d['rank_full']]
        mismatches_more = [(d['n_rc'], d['rank_full'], d['n_tts'], d['beta1'])
                           for d in all_data if d['n_rc'] > d['rank_full']]
        if mismatches_less:
            print(f"  #RC < rank: {len(mismatches_less)} cases")
            for m in mismatches_less[:5]:
                print(f"    #RC={m[0]}, rank={m[1]}, #TTs={m[2]}, β₁={m[3]}")
        if mismatches_more:
            print(f"  #RC > rank: {len(mismatches_more)} cases")
            for m in mismatches_more[:5]:
                print(f"    #RC={m[0]}, rank={m[1]}, #TTs={m[2]}, β₁={m[3]}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
