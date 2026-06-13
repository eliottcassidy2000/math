#!/usr/bin/env python3
"""
ALGEBRAIC PROOF THAT β_2 = 0 FOR ALL TOURNAMENTS — v2

STRATEGY: Show ker(∂_2|Ω_2) ⊆ im(∂_3|Ω_3) using the tournament structure.

KEY STRUCTURAL FACTS:
1. Ω_2 = {transitive triples}: (a,b,c) with a→b, b→c, AND a→c
2. Ω_3 = {doubly transitive 4-paths}: (a,b,c,d) with a→b→c→d AND a→c AND b→d
3. ∂_2(a,b,c) = (b,c) - (a,c) + (a,b)
4. ∂_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

A 2-cycle z = Σ α_t · t (t ∈ Ω_2) satisfies ∂_2 z = 0.
We need: z = ∂_3(w) for some w ∈ Ω_3.

APPROACH: Examine the structure of ker(∂_2|Ω_2) directly.
Each edge appears in net zero in ∂_2 z = 0.
"""
import numpy as np
from itertools import combinations
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== Understand the kernel structure at n=4 =====
print("=" * 70)
print("KERNEL OF ∂_2 ON Ω_2: DETAILED STRUCTURE")
print("=" * 70)

n = 4
count = 0
for A in all_tournaments_gen(n):
    a2 = enumerate_allowed_paths(A, n, 2)
    a1 = enumerate_allowed_paths(A, n, 1)
    a3 = enumerate_allowed_paths(A, n, 3)

    # Ω_2 = transitive triples
    trans_triples = []
    for path in a2:
        v0, v1, v2 = path
        if A[v0][v2] == 1:  # shortcut exists
            trans_triples.append(path)

    if len(trans_triples) == 0:
        continue

    # Build ∂_2 restricted to transitive triples
    edge_list = [(i,j) for i in range(n) for j in range(n) if A[i][j] == 1]
    edge_idx = {e: i for i, e in enumerate(edge_list)}

    bd2 = np.zeros((len(edge_list), len(trans_triples)))
    for j, (a, b, c) in enumerate(trans_triples):
        bd2[edge_idx[(b,c)], j] += 1   # (b,c)
        bd2[edge_idx[(a,c)], j] -= 1   # -(a,c)
        bd2[edge_idx[(a,b)], j] += 1   # (a,b)

    rank = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(trans_triples) - rank

    # Ω_3 = doubly transitive 4-paths
    dt_paths = []
    for path in a3:
        v0, v1, v2, v3 = path
        if A[v0][v2] == 1 and A[v1][v3] == 1:
            dt_paths.append(path)

    # im(∂_3) in terms of transitive triples
    if len(dt_paths) > 0:
        bd3 = np.zeros((len(trans_triples), len(dt_paths)))
        tt_idx = {tuple(t): i for i, t in enumerate(trans_triples)}
        for j, (a, b, c, d) in enumerate(dt_paths):
            # ∂_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
            for face, sign in [((b,c,d), 1), ((a,c,d), -1), ((a,b,d), 1), ((a,b,c), -1)]:
                if face in tt_idx:
                    bd3[tt_idx[face], j] += sign
        im_rank = np.linalg.matrix_rank(bd3, tol=1e-8)
    else:
        im_rank = 0

    if ker_dim > 0 and count < 10:
        count += 1
        print(f"\nTournament with ker(∂_2)={ker_dim}:")
        print(f"  |Ω_2|={len(trans_triples)}, |Ω_3|={len(dt_paths)}")
        print(f"  rank(∂_2)={rank}, ker=dim(Ω_2)-rank={ker_dim}")
        print(f"  rank(∂_3|Ω_2)={im_rank}")
        print(f"  β_2 = {ker_dim - im_rank}")

        # Find the actual kernel vector
        U, S_vals, Vt = np.linalg.svd(bd2, full_matrices=True)
        null_start = rank
        for ki in range(null_start, len(trans_triples)):
            vec = Vt[ki]
            nonzero = [(trans_triples[j], vec[j]) for j in range(len(vec)) if abs(vec[j]) > 1e-8]
            print(f"  Kernel vector {ki-null_start}: ", end="")
            for (a,b,c), coeff in nonzero:
                print(f"  {coeff:+.3f}·({a},{b},{c})", end="")
            print()

        # Show the transitive triples
        print(f"  Transitive triples:", end="")
        for (a,b,c) in trans_triples:
            print(f" ({a},{b},{c})", end="")
        print()

        # Show the doubly-transitive 4-paths
        print(f"  DT 4-paths:", end="")
        for (a,b,c,d) in dt_paths:
            print(f" ({a},{b},{c},{d})", end="")
        print()

# ===== KEY QUESTION: What creates the kernel? =====
print("\n\n" + "=" * 70)
print("CYCLE STRUCTURE IN Ω_2")
print("=" * 70)

print("""
Consider the TRIPLE GRAPH: vertices = transitive triples, edges = shared edges.
A 2-cycle in Ω_2 is a formal sum Σ α_t · t where each edge appears with net 0.

For each transitive triple t = (a,b,c), the boundary has 3 edges:
  ∂t = +(b,c) - (a,c) + (a,b)

For z = Σ α_t · t to be a cycle:
  For each edge (u,v): Σ_{t containing (u,v)} ± α_t = 0

In a tournament, every pair of vertices has exactly one edge.
So the constraint is: for each ordered pair (u,v) with u→v,
the alternating sum of coefficients of triples using edge (u,v) must be 0.
""")

# n=5: Detailed analysis of ker(∂_2|Ω_2) structure
print("=" * 70)
print("n=5 ANALYSIS: CYCLE TYPES")
print("=" * 70)

n = 5
cycle_types = {}
count = 0
for idx, A in enumerate(all_tournaments_gen(n)):
    a2 = enumerate_allowed_paths(A, n, 2)
    a1 = enumerate_allowed_paths(A, n, 1)

    # Ω_2
    trans_triples = [p for p in a2 if A[p[0]][p[2]] == 1]
    if len(trans_triples) == 0:
        continue

    edge_list = [(i,j) for i in range(n) for j in range(n) if A[i][j] == 1]
    edge_idx = {e: i for i, e in enumerate(edge_list)}

    bd2 = np.zeros((len(edge_list), len(trans_triples)))
    for j, (a, b, c) in enumerate(trans_triples):
        bd2[edge_idx[(b,c)], j] += 1
        bd2[edge_idx[(a,c)], j] -= 1
        bd2[edge_idx[(a,b)], j] += 1

    rank = np.linalg.matrix_rank(bd2, tol=1e-8)
    ker_dim = len(trans_triples) - rank

    if ker_dim > 0:
        count += 1

        # Count 3-cycles in tournament
        t3 = 0
        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    if (A[i][j]+A[j][i]) == 1 and (A[i][k]+A[k][i]) == 1 and (A[j][k]+A[k][j]) == 1:
                        s = A[i][j] + A[j][k] + A[k][i]
                        if s == 0 or s == 3:
                            continue
                        # Check if it's a 3-cycle
                        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                            t3 += 1

        key = (len(trans_triples), ker_dim, t3)
        cycle_types[key] = cycle_types.get(key, 0) + 1

print(f"Tournaments with ker(∂_2) > 0: {count}/{2**10}")
print(f"\nDistribution (|Ω_2|, ker_dim, t3):")
for key in sorted(cycle_types.keys()):
    print(f"  {key}: {cycle_types[key]}")

# ===== THE KEY INSIGHT: Look at 4-vertex subtournaments =====
print("\n\n" + "=" * 70)
print("4-VERTEX SUBTOURNAMENT CONNECTION")
print("=" * 70)

print("""
Conjecture: β_2 = 0 because for every 4-vertex subset {a,b,c,d},
if its transitive triples form a cycle, it's always filled by a
doubly-transitive 4-path.

There are 4 non-isomorphic tournaments on 4 vertices:
1. Transitive: 0→1→2→3 with all shortcuts (TT4)
2. One 3-cycle: 0→1→2→0, 3 wins/loses to all (two types)
3. Regular: 0→1→2→3→0 with shortcuts 0→2 and 1→3

Let's check which have ker(∂_2) > 0 and which have enough Ω_3 to fill it.
""")

n = 4
from collections import Counter
patterns = Counter()
for A in all_tournaments_gen(n):
    # Score sequence
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))

    a2 = enumerate_allowed_paths(A, n, 2)
    trans_triples = [p for p in a2 if A[p[0]][p[2]] == 1]
    a3 = enumerate_allowed_paths(A, n, 3)
    dt_paths = [(v0,v1,v2,v3) for (v0,v1,v2,v3) in a3 if A[v0][v2]==1 and A[v1][v3]==1]

    edge_list = [(i,j) for i in range(n) for j in range(n) if A[i][j] == 1]
    edge_idx = {e: i for i, e in enumerate(edge_list)}

    if len(trans_triples) > 0:
        bd2 = np.zeros((len(edge_list), len(trans_triples)))
        for j, (a, b, c) in enumerate(trans_triples):
            bd2[edge_idx[(b,c)], j] += 1
            bd2[edge_idx[(a,c)], j] -= 1
            bd2[edge_idx[(a,b)], j] += 1
        rank = np.linalg.matrix_rank(bd2, tol=1e-8)
        ker_dim = len(trans_triples) - rank
    else:
        ker_dim = 0

    patterns[(scores, len(trans_triples), ker_dim, len(dt_paths))] += 1

print(f"n=4 patterns (scores, |Ω_2|, ker∂_2, |Ω_3|):")
for key in sorted(patterns.keys()):
    print(f"  {key}: {patterns[key]}")

# ===== REGULAR TOURNAMENT ANALYSIS =====
print("\n\n" + "=" * 70)
print("REGULAR 4-TOURNAMENT: THE KEY CASE")
print("=" * 70)

# Build a specific regular tournament on 4 vertices: 0→1→2→3→0, 0→2, 1→3
n = 4
A = [[0]*n for _ in range(n)]
A[0][1] = A[1][2] = A[2][3] = A[3][0] = 1  # cycle 0→1→2→3→0
A[0][2] = A[1][3] = 1  # shortcuts

print(f"Regular tournament: 0→1→2→3→0, 0→2, 1→3")
print(f"Adjacency: ", end="")
for i in range(n):
    for j in range(n):
        if A[i][j]: print(f"{i}→{j} ", end="")
print()

a2 = enumerate_allowed_paths(A, n, 2)
trans_triples = [(a,b,c) for (a,b,c) in a2 if A[a][c] == 1]
print(f"\nTransitive triples (Ω_2):")
for (a,b,c) in trans_triples:
    print(f"  ({a},{b},{c}): ∂ = +({b},{c}) -({a},{c}) +({a},{b})")

a3 = enumerate_allowed_paths(A, n, 3)
dt_paths = [(a,b,c,d) for (a,b,c,d) in a3 if A[a][c]==1 and A[b][d]==1]
print(f"\nDoubly-transitive 4-paths (Ω_3):")
for (a,b,c,d) in dt_paths:
    bd = f"+({b},{c},{d}) -({a},{c},{d}) +({a},{b},{d}) -({a},{b},{c})"
    print(f"  ({a},{b},{c},{d}): ∂ = {bd}")

# Show boundary ∂_2
edge_list = [(i,j) for i in range(n) for j in range(n) if A[i][j] == 1]
edge_idx = {e: i for i, e in enumerate(edge_list)}

bd2 = np.zeros((len(edge_list), len(trans_triples)))
for j, (a, b, c) in enumerate(trans_triples):
    bd2[edge_idx[(b,c)], j] += 1
    bd2[edge_idx[(a,c)], j] -= 1
    bd2[edge_idx[(a,b)], j] += 1

print(f"\n∂_2 matrix ({len(edge_list)} edges × {len(trans_triples)} triples):")
for i, e in enumerate(edge_list):
    row = [int(bd2[i,j]) for j in range(len(trans_triples))]
    print(f"  {e}: {row}")

rank = np.linalg.matrix_rank(bd2, tol=1e-8)
ker = len(trans_triples) - rank
print(f"\nrank(∂_2) = {rank}, ker = {ker}")

if ker > 0:
    U, S_vals, Vt = np.linalg.svd(bd2, full_matrices=True)
    for ki in range(rank, len(trans_triples)):
        vec = Vt[ki]
        print(f"\nKernel vector: ", end="")
        for j, (a,b,c) in enumerate(trans_triples):
            if abs(vec[j]) > 1e-8:
                print(f" {vec[j]:+.3f}·({a},{b},{c})", end="")
        print()

    # Show that this kernel vector is in im(∂_3)
    if len(dt_paths) > 0:
        tt_idx = {t: i for i, t in enumerate(trans_triples)}
        bd3 = np.zeros((len(trans_triples), len(dt_paths)))
        for j, (a, b, c, d) in enumerate(dt_paths):
            for face, sign in [((b,c,d), 1), ((a,c,d), -1), ((a,b,d), 1), ((a,b,c), -1)]:
                if face in tt_idx:
                    bd3[tt_idx[face], j] += sign
        print(f"\n∂_3 matrix ({len(trans_triples)} triples × {len(dt_paths)} DT-paths):")
        for i, t in enumerate(trans_triples):
            row = [int(bd3[i,j]) for j in range(len(dt_paths))]
            print(f"  {t}: {row}")

        # Check: does im(∂_3) contain the kernel?
        print(f"\nColumn space of ∂_3:")
        for j, (a,b,c,d) in enumerate(dt_paths):
            col = bd3[:, j]
            print(f"  ∂_3({a},{b},{c},{d}) = ", end="")
            for i, t in enumerate(trans_triples):
                if abs(col[i]) > 1e-8:
                    print(f" {int(col[i]):+d}·{t}", end="")
            print()

print("\nDone.")
