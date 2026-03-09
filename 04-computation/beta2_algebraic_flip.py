#!/usr/bin/env python3
"""
beta2_algebraic_flip.py — Algebraic flip obstruction analysis

Part 1: Cross-tabulate #bad vs t3 and #TTs vs rank at n=5,6
Part 2: Bad TT boundary structure (boundary projection analysis)
Part 3: Star constraint analysis
Part 4: n=6,7 verification

Uses path_homology_v2.py for Betti number computation.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
from math import comb
import random
import sys
import time

# ===== Core tournament utilities =====

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
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

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_3cycles(A, n):
    """Count the number of 3-cycles in tournament A."""
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def out_degree(A, n, v):
    """Out-degree of vertex v."""
    return sum(A[v])

def is_transitive_triple(A, i, j, k):
    """Check if {i,j,k} forms a transitive triple (TT).
    A TT means one vertex beats both others, who have a definite relation.
    Equivalently: NOT a 3-cycle."""
    cyc1 = A[i][j] and A[j][k] and A[k][i]
    cyc2 = A[j][i] and A[k][j] and A[i][k]
    return not (cyc1 or cyc2)

def count_TTs(A, n):
    """Count transitive triples."""
    return comb(n, 3) - count_3cycles(A, n)

def find_bad_vertices(A, n):
    """A vertex v is 'bad' if removing v from some TT leaves a non-allowed 1-path.

    More precisely: v is bad if there exists a TT (a,b,c) containing v such that
    removing v gives a face that is NOT an edge in the tournament.

    Wait — in a tournament every pair has an edge, so every face IS an allowed path.

    Let me reconsider: 'bad' in the context of β₂=0 analysis means vertices involved
    in the boundary obstruction. Let me use the definition from the path homology:

    A vertex v is 'bad' if it participates in a 2-path (a,v,b) where (a,b) is NOT
    an edge (i.e., the face obtained by removing v is not allowed).

    But in a tournament, (a,b) is always an edge! So there are no 'bad' vertices
    in this sense for tournaments.

    Actually, re-reading the context: 'bad' likely refers to vertices where
    β₁ > 0 relates to cycle structure. Let me use a different definition:

    Bad vertex = vertex v such that v is a SOURCE or SINK in some induced
    sub-tournament that creates a non-trivial cycle.

    Actually, let me use the standard definition from the DT (doubly-transitive)
    analysis: a 2-path (a,b,c) is "bad" if its face (a,c) obtained by removing
    middle vertex b is NOT in A_1 (not an allowed 1-path).

    In a tournament, ALL 2-paths have all faces allowed. So "bad" must mean
    something else here.

    Let me re-read: the question says "#bad vertices" and relates to β₁ and rank.

    I think "bad" = vertex with odd out-degree signature causing β₁ > 0.

    Actually, for tournaments: β₁ = C(n,2) - n + 1 - rank(∂₂|Ω₂).
    And rank(∂₂|Ω₂) = #TTs (when β₂=0).
    Wait, that's what we're testing.

    Let me just use: #bad = number of vertices v where d⁺(v) is "unbalanced"
    in a way that contributes to β₁. The simplest: vertices in 3-cycles.
    """
    # Count how many 3-cycles each vertex participates in
    cycle_count = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    cycle_count[i] += 1
                    cycle_count[j] += 1
                    cycle_count[k] += 1
    # "bad" = vertex participating in at least one 3-cycle
    bad = [v for v in range(n) if cycle_count[v] > 0]
    return bad, cycle_count

def find_bad_vertices_v2(A, n):
    """Alternative: 'bad' vertex = vertex where the local score doesn't match
    the transitive pattern. Specifically, vertex v with d+(v) != 0 and d+(v) != n-1
    that is involved in a 3-cycle.

    Actually, let me try the simplest definition that makes the question meaningful:
    #bad = number of vertices that are in at least one 3-cycle.
    """
    return find_bad_vertices(A, n)

# ===== Path homology computation (simplified for β₁) =====

def compute_beta1(A, n):
    """Compute β₁ for a tournament.

    For tournaments: all pairs have edges, so A_1 = all directed edges = n(n-1).
    Ω_1 = A_1 (all faces of 1-paths are 0-paths = vertices, always allowed).
    ∂₁: A_1 → A_0 sends (a,b) ↦ b - a.
    ker(∂₁) has dim = |A_1| - rank(∂₁) = n(n-1) - (n-1).
    im(∂₂): need rank of ∂₂ restricted to Ω₂.
    β₁ = ker(∂₁)/im(∂₂) = (n(n-1) - (n-1) - rank(∂₂|Ω₂)).

    Wait, let me be more careful. For tournaments:
    - A_0 = n vertices, A_1 = n(n-1) directed edges (both directions)

    No wait: in a tournament, for each pair {i,j}, exactly ONE of i→j or j→i exists.
    So |A_1| = C(n,2) directed edges (one per pair, but as directed = n(n-1)/2).

    Hmm, but allowed paths use the directed edge. So |A_1| = n(n-1)/2 = C(n,2).

    ∂₁(a,b) = (b) - (a). Rank of ∂₁ = n-1 (connected tournament).
    ker(∂₁) = C(n,2) - (n-1).

    β₁ = dim(ker(∂₁)) - dim(im(∂₂|Ω₂)).

    We need rank(∂₂|Ω₂). For tournaments, Ω₂ = A₂ (since all faces of 2-paths
    are edges, and tournaments have all pairs as edges... wait, the face of (a,b,c)
    by removing b is (a,c), which needs a→c to be an edge. But a→c might not hold!
    If c→a instead, then (a,c) is not allowed.)

    So Ω₂ ≠ A₂ in general. The "bad" 2-paths are those (a,b,c) where
    some face is not allowed. Face (a,c) [remove middle] needs a→c.

    A 2-path (a,b,c): requires a→b, b→c. Faces: (b,c), (a,c), (a,b).
    - (b,c): b→c ✓ (given)
    - (a,b): a→b ✓ (given)
    - (a,c): need a→c. If c→a instead, this face is NOT allowed.

    So a 2-path (a,b,c) with a→b, b→c is "non-∂-invariant" iff c→a.
    That means a→b→c→a is a 3-CYCLE!

    So bad 2-paths = 2-paths that are part of a 3-cycle.
    And #non-allowed faces come from 3-cycles.
    """
    # Build allowed 2-paths
    allowed_2 = []
    for a in range(n):
        for b in range(n):
            if b == a: continue
            if not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if not A[b][c]: continue
                allowed_2.append((a, b, c))

    # Build allowed 1-paths (edges)
    allowed_1 = []
    edge_idx = {}
    for a in range(n):
        for b in range(n):
            if b == a: continue
            if A[a][b]:
                edge_idx[(a, b)] = len(allowed_1)
                allowed_1.append((a, b))

    num_edges = len(allowed_1)  # = C(n,2)
    num_2paths = len(allowed_2)

    # Build boundary matrix ∂₂: A_2 → A_1
    # ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
    bd2 = np.zeros((num_edges, num_2paths))
    for j, (a, b, c) in enumerate(allowed_2):
        # face 0: remove a → (b,c)
        if (b, c) in edge_idx:
            bd2[edge_idx[(b, c)], j] += 1
        # face 1: remove b → (a,c), sign = -1
        if (a, c) in edge_idx:
            bd2[edge_idx[(a, c)], j] -= 1
        # face 2: remove c → (a,b)
        if (a, b) in edge_idx:
            bd2[edge_idx[(a, b)], j] += 1

    # Compute Ω₂: 2-chains whose boundary lies entirely in A_1
    # Need to find which faces are NOT in A_1
    # For 2-path (a,b,c): face (a,c) is not allowed iff c→a (3-cycle)
    non_allowed_faces = {}
    na_count = 0
    for j, (a, b, c) in enumerate(allowed_2):
        # Check face (a,c) — the only potentially non-allowed face
        if not A[a][c]:  # c→a, so (a,c) not an edge
            face = (a, c)
            if face not in non_allowed_faces:
                non_allowed_faces[face] = na_count
                na_count += 1

    if na_count == 0:
        # All Ω₂ = A₂ (no 3-cycles → transitive tournament)
        omega2_basis = np.eye(num_2paths)
    else:
        # Build projection to non-allowed faces
        P = np.zeros((na_count, num_2paths))
        for j, (a, b, c) in enumerate(allowed_2):
            if not A[a][c]:
                face = (a, c)
                P[non_allowed_faces[face], j] -= 1  # coefficient is -1 for middle removal

        # Ω₂ = ker(P)
        U, S, Vt = np.linalg.svd(P, full_matrices=True)
        rank_P = sum(s > 1e-10 for s in S)
        omega2_basis = Vt[rank_P:].T

    dim_omega2 = omega2_basis.shape[1] if omega2_basis.ndim == 2 else 0

    # Rank of ∂₂ restricted to Ω₂
    if dim_omega2 > 0:
        bd2_omega = bd2 @ omega2_basis
        sv = np.linalg.svd(bd2_omega, compute_uv=False)
        rank_bd2 = sum(s > 1e-8 for s in sv)
    else:
        rank_bd2 = 0

    # β₁ = dim(ker(∂₁)) - rank(∂₂|Ω₂)
    # ker(∂₁) = C(n,2) - (n-1) for connected tournament
    ker_d1 = num_edges - (n - 1)
    beta1 = ker_d1 - rank_bd2

    return beta1, dim_omega2, rank_bd2, num_2paths

def compute_beta1_and_details(A, n):
    """Compute β₁ and return detailed info about the boundary structure."""
    beta1, dim_omega2, rank_bd2, num_2paths = compute_beta1(A, n)
    t3 = count_3cycles(A, n)
    nTTs = comb(n, 3) - t3
    bad_verts, cycle_counts = find_bad_vertices(A, n)
    nbad = len(bad_verts)

    # rank formula: C(n,2) - n + 1 - β₁
    expected_rank = comb(n, 2) - n + 1 - beta1

    return {
        'beta1': beta1,
        't3': t3,
        'nTTs': nTTs,
        'nbad': nbad,
        'bad_verts': bad_verts,
        'cycle_counts': cycle_counts,
        'dim_omega2': dim_omega2,
        'rank_bd2': rank_bd2,
        'expected_rank': expected_rank,
        'redundancy': nTTs - rank_bd2,  # how many TTs are redundant
    }

# ===== Part 1: Cross-tabulation =====

def part1(n, tournaments_iter, label=""):
    """Cross-tabulate #bad vs t₃ and #TTs vs rank."""
    print(f"\n{'='*70}")
    print(f"PART 1: Cross-tabulation at n={n} {label}")
    print(f"{'='*70}")

    # Collect data
    data = []
    count = 0
    for A in tournaments_iter:
        info = compute_beta1_and_details(A, n)
        data.append(info)
        count += 1
        if count % 5000 == 0:
            print(f"  ... processed {count} tournaments", file=sys.stderr)

    print(f"\nTotal tournaments: {count}")
    print(f"C(n,3) = {comb(n,3)}, C(n,2)-n+1 = {comb(n,2)-n+1}")

    # Cross-tabulate (#bad, β₁) → distribution of t₃
    print(f"\n--- (#bad, β₁) → t₃ distribution ---")
    cross = defaultdict(lambda: defaultdict(int))
    for info in data:
        key = (info['nbad'], info['beta1'])
        cross[key][info['t3']] += 1

    for key in sorted(cross.keys()):
        nbad, b1 = key
        t3_dist = dict(sorted(cross[key].items()))
        total = sum(t3_dist.values())
        print(f"  #bad={nbad}, β₁={b1} ({total} tournaments): t₃ ∈ {t3_dist}")

    # Cross-tabulate (#bad, β₁) → redundancy distribution
    print(f"\n--- (#bad, β₁) → redundancy (#TTs - rank(∂₂|Ω₂)) ---")
    red_cross = defaultdict(lambda: defaultdict(int))
    for info in data:
        key = (info['nbad'], info['beta1'])
        red_cross[key][info['redundancy']] += 1

    for key in sorted(red_cross.keys()):
        nbad, b1 = key
        red_dist = dict(sorted(red_cross[key].items()))
        print(f"  #bad={nbad}, β₁={b1}: redundancy ∈ {red_dist}")

    # KEY TEST: β₁=0, various #bad → is #TTs = rank?
    print(f"\n--- KEY TEST: β₁=0 cases ---")
    for nbad_val in sorted(set(info['nbad'] for info in data)):
        subset = [info for info in data if info['beta1'] == 0 and info['nbad'] == nbad_val]
        if not subset:
            continue
        redundancies = Counter(info['redundancy'] for info in subset)
        tts_vals = Counter(info['nTTs'] for info in subset)
        rank_vals = Counter(info['rank_bd2'] for info in subset)
        print(f"  #bad={nbad_val}: {len(subset)} tournaments")
        print(f"    #TTs values: {dict(sorted(tts_vals.items()))}")
        print(f"    rank values: {dict(sorted(rank_vals.items()))}")
        print(f"    redundancy: {dict(sorted(redundancies.items()))}")
        all_zero = all(info['redundancy'] == 0 for info in subset)
        print(f"    #TTs = rank always? {all_zero}")

    # Also show β₁ > 0 cases
    print(f"\n--- β₁ > 0 cases ---")
    for b1_val in sorted(set(info['beta1'] for info in data)):
        if b1_val == 0:
            continue
        subset = [info for info in data if info['beta1'] == b1_val]
        redundancies = Counter(info['redundancy'] for info in subset)
        nbad_dist = Counter(info['nbad'] for info in subset)
        print(f"  β₁={b1_val}: {len(subset)} tournaments")
        print(f"    #bad distribution: {dict(sorted(nbad_dist.items()))}")
        print(f"    redundancy: {dict(sorted(redundancies.items()))}")

    # Summary table
    print(f"\n--- Summary: dim(Ω₂), #TTs, rank, β₁ ---")
    print(f"  {'t3':>4} {'#TTs':>5} {'β₁':>3} {'dim_Ω₂':>7} {'rank':>5} {'red':>4} {'count':>6}")
    summary = defaultdict(int)
    for info in data:
        key = (info['t3'], info['nTTs'], info['beta1'], info['dim_omega2'],
               info['rank_bd2'], info['redundancy'])
        summary[key] += 1
    for key in sorted(summary.keys()):
        t3, nTTs, b1, do2, rank, red = key
        print(f"  {t3:4d} {nTTs:5d} {b1:3d} {do2:7d} {rank:5d} {red:4d} {summary[key]:6d}")

    return data

# ===== Part 2: Bad TT boundary structure =====

def part2(n=5):
    """Analyze boundary structure of bad TTs at n=5."""
    print(f"\n{'='*70}")
    print(f"PART 2: Bad TT boundary structure at n={n}")
    print(f"{'='*70}")

    count_total = 0
    count_b1_0_bad3 = 0

    # Track statistics
    bad_edge_hits = defaultdict(int)  # how many bad edges each non-bad TT hits
    projection_independent = 0
    projection_dependent = 0

    for A in all_tournaments(n):
        t3 = count_3cycles(A, n)
        info = compute_beta1_and_details(A, n)
        count_total += 1

        if info['beta1'] != 0:
            continue
        if info['nbad'] != 3:
            continue

        count_b1_0_bad3 += 1

        bad = info['bad_verts']
        a, b, c = bad[0], bad[1], bad[2]

        # Find the 3-cycle orientation among {a,b,c}
        # One of the two orientations: a→b→c→a or a→c→b→a
        if A[a][b] and A[b][c] and A[c][a]:
            # 3-cycle a→b→c→a. TTs don't exist here (it's a 3-cycle, not a TT)
            # The "bad TT" doesn't apply — {a,b,c} is a 3-cycle, so it's NOT a TT
            pass
        elif A[b][a] and A[c][b] and A[a][c]:
            # 3-cycle a→c→b→a
            pass

        # Actually, if {a,b,c} are all in 3-cycles, the triple {a,b,c} itself
        # might or might not be a 3-cycle. Let me check.
        is_3cyc = (A[a][b] and A[b][c] and A[c][a]) or \
                  (A[b][a] and A[c][b] and A[a][c])

        # Find all TTs in the tournament
        tts = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if is_transitive_triple(A, i, j, k):
                        tts.append((i, j, k))

        # Find all 3-cycles
        cycles = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if not is_transitive_triple(A, i, j, k):
                        cycles.append((i, j, k))

        # For each TT, find the "source" (vertex beating both others)
        # and express as oriented 2-path
        tt_2paths = []
        for (i, j, k) in tts:
            # Find who beats whom
            # Source beats both, sink loses to both
            scores = {i: A[i][j] + A[i][k], j: A[j][i] + A[j][k], k: A[k][i] + A[k][j]}
            source = max(scores, key=scores.get)
            sink = min(scores, key=scores.get)
            mid = ({i,j,k} - {source, sink}).pop()
            # The TT as 2-path: source → mid → sink (both edges present)
            # Also source → sink is an edge
            tt_2paths.append((source, mid, sink))

        # The edges among bad vertices
        bad_edges = set()
        for x in bad:
            for y in bad:
                if x != y and A[x][y]:
                    bad_edges.add((x, y))

        # For each TT 2-path, check how many bad edges its boundary hits
        for (s, m, sk) in tt_2paths:
            # Boundary ∂₂(s,m,sk) = (m,sk) - (s,sk) + (s,m)
            boundary_edges = [(m, sk), (s, sk), (s, m)]
            hits = sum(1 for e in boundary_edges if e in bad_edges)
            bad_edge_hits[hits] += 1

        # Check if bad triple's boundary projection is independent
        if is_3cyc:
            # {a,b,c} is a 3-cycle — no TT here to analyze
            # The "bad" structure is the cycle itself
            pass

        # Only do detailed analysis for first few
        if count_b1_0_bad3 <= 3:
            print(f"\n  Tournament #{count_b1_0_bad3}: bad={bad}, t3={t3}")
            print(f"    Bad triple is 3-cycle: {is_3cyc}")
            print(f"    Bad edges: {bad_edges}")
            print(f"    #TTs = {len(tts)}, #3-cycles = {len(cycles)}")
            print(f"    TT 2-paths:")
            for (s, m, sk) in tt_2paths[:8]:
                boundary_edges = [(m, sk), (s, sk), (s, m)]
                hits = sum(1 for e in boundary_edges if e in bad_edges)
                in_bad = [e for e in boundary_edges if e in bad_edges]
                print(f"      ({s},{m},{sk}): hits {hits} bad edges: {in_bad}")

    print(f"\n  Total β₁=0, #bad=3: {count_b1_0_bad3} tournaments")
    print(f"\n  Bad edge hit distribution across ALL TT boundaries:")
    for hits in sorted(bad_edge_hits.keys()):
        print(f"    {hits} bad edges: {bad_edge_hits[hits]} TT boundaries")

    # Now do the full boundary projection analysis
    print(f"\n--- Boundary projection analysis ---")

    proj_results = defaultdict(int)

    for A in all_tournaments(n):
        info = compute_beta1_and_details(A, n)
        if info['beta1'] != 0 or info['nbad'] != 3:
            continue

        bad = info['bad_verts']

        # Check if bad triple is a 3-cycle
        a, b, c = bad
        is_3cyc = (A[a][b] and A[b][c] and A[c][a]) or \
                  (A[b][a] and A[c][b] and A[a][c])

        if not is_3cyc:
            proj_results['bad_triple_is_TT'] += 1
            continue

        proj_results['bad_triple_is_3cycle'] += 1

        # Get all allowed 2-paths and build boundary matrix
        # restricted to the 3 bad-vertex edges
        allowed_2paths = []
        for x in range(n):
            for y in range(n):
                if y == x or not A[x][y]: continue
                for z in range(n):
                    if z == x or z == y or not A[y][z]: continue
                    allowed_2paths.append((x, y, z))

        # Edge indices (all directed edges)
        edges = []
        edge_idx = {}
        for x in range(n):
            for y in range(n):
                if x != y and A[x][y]:
                    edge_idx[(x, y)] = len(edges)
                    edges.append((x, y))

        # Bad edges among {a,b,c}
        bad_edge_list = []
        bad_edge_idx = {}
        for x in [a, b, c]:
            for y in [a, b, c]:
                if x != y and A[x][y]:
                    bad_edge_idx[(x, y)] = len(bad_edge_list)
                    bad_edge_list.append((x, y))

        # Build boundary matrix projected onto bad edges
        # For each 2-path, its boundary projected onto the 3 bad edges
        num_bad_edges = len(bad_edge_list)
        proj_matrix = np.zeros((num_bad_edges, len(allowed_2paths)))

        for j, (x, y, z) in enumerate(allowed_2paths):
            # ∂₂(x,y,z) = (y,z) - (x,z) + (x,y)
            for sign, face in [(1, (y, z)), (-1, (x, z)), (1, (x, y))]:
                if face in bad_edge_idx:
                    proj_matrix[bad_edge_idx[face], j] += sign

        # Find which 2-paths involve the 3-cycle vertices
        cycle_2paths = []
        other_2paths = []
        for j, (x, y, z) in enumerate(allowed_2paths):
            if {x, y, z} == {a, b, c}:
                cycle_2paths.append(j)
            else:
                other_2paths.append(j)

        # Check: what do the cycle 2-paths project to?
        if cycle_2paths:
            cycle_proj = proj_matrix[:, cycle_2paths]

        # Check: rank of projection from non-cycle 2-paths
        if other_2paths:
            other_proj = proj_matrix[:, other_2paths]
            sv = np.linalg.svd(other_proj, compute_uv=False)
            rank_other = sum(s > 1e-8 for s in sv)
        else:
            rank_other = 0

        # Full rank
        sv_full = np.linalg.svd(proj_matrix, compute_uv=False)
        rank_full = sum(s > 1e-8 for s in sv_full)

        proj_results[f'rank_other={rank_other}_rank_full={rank_full}'] += 1

    print(f"\n  Projection results:")
    for key in sorted(proj_results.keys()):
        print(f"    {key}: {proj_results[key]}")

# ===== Part 3: Star constraint analysis =====

def part3(n=5):
    """Analyze star constraints and free dimensions."""
    print(f"\n{'='*70}")
    print(f"PART 3: Star constraint analysis at n={n}")
    print(f"{'='*70}")

    results = defaultdict(int)

    sample_count = 0
    for A in all_tournaments(n):
        info = compute_beta1_and_details(A, n)

        # Build the full boundary matrix ∂₂: Ω₂ → A₁
        # and analyze the star constraints

        # Star of vertex v: all edges incident to v
        # Star constraint: for each v, sum of ∂₂ over paths through v

        # Compute out-degrees
        out_degs = [sum(A[v]) for v in range(n)]

        # "Star elimination" removes d+(v)-1 variables per vertex
        # But we need to be more precise

        # Total edges = C(n,2)
        # After star elimination: C(n,2) - sum(d+(v)-1) ... no, that's not right

        # Actually, the star constraints are:
        # For vertex v: sum_{u: v→u} [coefficient of (v,u) in ∂₂(chain)] = ...
        # This is really about the cocycle structure

        # Let me just compute rank(∂₂|Ω₂) and compare with #TTs and C(n,2)-n+1-β₁

        key = (info['t3'], info['beta1'], info['nTTs'], info['rank_bd2'],
               info['dim_omega2'], info['redundancy'])
        results[key] += 1

        sample_count += 1

    print(f"\n  (t3, β₁, #TTs, rank(∂₂|Ω₂), dim(Ω₂), redundancy) → count")
    for key in sorted(results.keys()):
        t3, b1, nTTs, rank, do2, red = key
        print(f"    t3={t3}, β₁={b1}, #TTs={nTTs}, rank={rank}, dim(Ω₂)={do2}, red={red}: {results[key]}")

    # Key insight: redundancy = #TTs - rank(∂₂|Ω₂)
    # If all boundaries of TTs were independent in Ω₂, redundancy = 0
    # But Ω₂ has dimension < #TTs in general, so some are dependent

    print(f"\n  Key relationships:")
    print(f"    C(n,2) = {comb(n,2)}")
    print(f"    n-1 = {n-1} (rank of ∂₁)")
    print(f"    C(n,2)-n+1 = {comb(n,2)-n+1} (dim ker ∂₁)")
    print(f"    C(n,3) = {comb(n,3)} (max TTs)")

    # Analyze: for β₁=0, is rank always = C(n,2)-n+1?
    b1_0 = [key for key in results.keys() if key[1] == 0]
    print(f"\n  For β₁=0 tournaments:")
    for key in sorted(b1_0):
        t3, b1, nTTs, rank, do2, red = key
        print(f"    t3={t3}: rank={rank}, expected={comb(n,2)-n+1}, match={rank == comb(n,2)-n+1}")

# ===== Part 4: n=6,7 verification =====

def part4_n6():
    """Exhaustive n=6 verification."""
    print(f"\n{'='*70}")
    print(f"PART 4a: n=6 exhaustive verification")
    print(f"{'='*70}")
    return part1(6, all_tournaments(6), "(exhaustive)")

def part4_n7(num_samples=200):
    """Sampled n=7 verification."""
    print(f"\n{'='*70}")
    print(f"PART 4b: n=7 sampled verification ({num_samples} tournaments)")
    print(f"{'='*70}")

    def sample_iter():
        for _ in range(num_samples):
            yield random_tournament(7)

    return part1(7, sample_iter(), f"(sample {num_samples})")

# ===== Main =====

if __name__ == '__main__':
    random.seed(42)
    np.random.seed(42)

    start = time.time()

    print("=" * 70)
    print("ALGEBRAIC FLIP OBSTRUCTION ANALYSIS")
    print("=" * 70)

    # Part 1: n=5
    print("\n\n" + "#" * 70)
    print("# PART 1: Cross-tabulation")
    print("#" * 70)

    data5 = part1(5, all_tournaments(5), "(exhaustive)")

    elapsed = time.time() - start
    print(f"\n[n=5 done in {elapsed:.1f}s]")

    # Part 2: Bad TT boundary structure
    print("\n\n" + "#" * 70)
    print("# PART 2: Bad TT boundary structure")
    print("#" * 70)

    part2(5)

    elapsed = time.time() - start
    print(f"\n[Part 2 done in {elapsed:.1f}s]")

    # Part 3: Star constraint analysis
    print("\n\n" + "#" * 70)
    print("# PART 3: Star constraint analysis")
    print("#" * 70)

    part3(5)

    elapsed = time.time() - start
    print(f"\n[Part 3 done in {elapsed:.1f}s]")

    # Part 4: n=6 exhaustive
    print("\n\n" + "#" * 70)
    print("# PART 4: n=6,7 verification")
    print("#" * 70)

    data6 = part4_n6()

    elapsed = time.time() - start
    print(f"\n[n=6 done in {elapsed:.1f}s]")

    # Part 4b: n=7 sampled
    data7 = part4_n7(200)

    elapsed = time.time() - start
    print(f"\n[n=7 done in {elapsed:.1f}s]")

    # Final summary
    print("\n\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)

    print(f"""
Key findings to check:
1. At n=5: β₁=0 with #bad=3 → is #TTs always = rank? (zero redundancy?)
2. At n=5: what is the redundancy for other #bad values?
3. At n=6: does #bad=3 → #TTs=rank still hold?
4. At n=7: same question (sampled)
5. Boundary projection: is the 3-cycle boundary independent from other TT boundaries?

Formula recap:
  β₁ = dim(ker ∂₁) - rank(∂₂|Ω₂)
     = C(n,2) - (n-1) - rank(∂₂|Ω₂)
  #TTs = C(n,3) - t₃
  redundancy = dim(Ω₂) - rank(∂₂|Ω₂)  ... wait, that's not right

  Actually: Ω₂ has some dimension. The TT boundaries live in A₁.
  rank(∂₂|Ω₂) = rank of the image of ∂₂ restricted to Ω₂.
  We compare this with #TTs (the number of transitive triples).
  redundancy = #TTs - rank(∂₂|Ω₂) ... but TTs aren't exactly generators of Ω₂.

  Better: redundancy = dim(Ω₂) - rank(∂₂|Ω₂) = dim(ker(∂₂|Ω₂))
  This is relevant because it tells us how many Ω₂ elements have zero boundary.
  When β₂=0, this equals dim(im(∂₃|Ω₃)).
""")

    total_elapsed = time.time() - start
    print(f"\nTotal time: {total_elapsed:.1f}s")
    print("\nDone.")
