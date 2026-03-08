#!/usr/bin/env python3
"""
even_betti_vanishing_proof.py - Why β_{2k} = 0 for tournaments?

CONJECTURE (opus-2026-03-07-S46e):
For any tournament T on n vertices, β_{2k}(T) = 0 for all k >= 1,
where β_p denotes the p-th GLMY path homology Betti number.

EVIDENCE:
- Exhaustive n=3,4,5: 0 failures (all 8 + 64 + 1024 tournaments)
- Sampled n=6: 0 failures in 500 random
- Sampled n=7: 0 failures in 600+ random
- β₂ = 0 ALWAYS, β₃ ∈ {0, 1}, β₁ ∈ {0, 1}

WHY SHOULD THIS BE TRUE?

Key fact: tournaments have ONLY ODD directed cycles.
In GLMY path homology, the boundary maps ∂_p: Ω_p → Ω_{p-1} act on
"allowed p-paths" (sequences of vertices where consecutive pairs form edges).

For a tournament T:
- Every allowed path v₀→v₁→...→v_p uses directed edges
- A p-cycle (boundary = 0) must "close up": v₀→v₁→...→v_p→v₀
- This requires p+1 edges forming a directed (p+1)-cycle
- But tournaments have only ODD directed cycles!

Wait, this isn't quite right — GLMY paths are NOT just simple cycles.
An "allowed p-path" in GLMY is a formal sum of elementary p-paths
(sequences of p+1 vertices v₀,...,v_p where each v_i→v_{i+1}).

The key is the boundary operator:
∂_p(v₀,...,v_p) = Σ_{i=0}^{p} (-1)^i (v₀,...,v̂_i,...,v_p)

where (v₀,...,v̂_i,...,v_p) is the path with v_i deleted, BUT only if
this shorter path is also "allowed" (all consecutive edges exist).

For tournaments, this is:
∂_p(v₀,...,v_p) = Σ_{i=0}^{p} (-1)^i [v_i removable] (v₀,...,v̂_i,...,v_p)

where [v_i removable] means T[v_{i-1}][v_{i+1}] = 1 (or boundary cases).

HYPOTHESIS: The even Betti vanishing is related to the ANTISYMMETRY
of the tournament adjacency: T[i][j] + T[j][i] = 1.

When we try to form a 2-cycle in Ω_2 (a formal sum of 2-paths = directed triangles),
the antisymmetry forces ALL boundary terms to be nonzero (every pair is connected),
creating a "rigid" boundary structure that prevents even-dimensional holes.

Let me investigate this algebraically at small sizes.

Author: opus-2026-03-07-S46e
"""
import numpy as np
from itertools import combinations, permutations
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths,
    build_full_boundary_matrix
)

# === Investigate the boundary maps at n=4 ===
print("=" * 70)
print("BOUNDARY MAP ANALYSIS AT n=4")
print("=" * 70)

# Tournament with t3=2 (the one with β₁=1)
# Score (1,1,2,2): 0→2, 0→3, 1→2, 1→3, 2→0, 3→1 ... wait
# Let me build the unique (up to iso) tournament with t3=2 at n=4
# Score sequence (1,1,2,2): the regular tournament C4
# 0→1, 1→2, 2→3, 3→0, and diagonal: 0→2 or 2→0, 1→3 or 3→1

# C4 tournament: 0→1, 1→2, 2→3, 3→0, 0→2 (transitive across diagonal)
# No wait. At n=4 with t3=2, we need exactly 2 directed 3-cycles.

# Let me just test a few specific tournaments
def build_adj(n, edges):
    """Build tournament from list of edges (i→j)."""
    A = [[0]*n for _ in range(n)]
    for i, j in edges:
        A[i][j] = 1
    # Fill in remaining (j→i)
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 0 and A[j][i] == 0:
                A[j][i] = 1
    return A

# n=4: tournament with t3=2
# Try the "rotational" tournament: 0→1, 1→2, 2→3, 3→0, 0→2, 1→3
A4 = [[0]*4 for _ in range(4)]
A4[0][1] = 1; A4[1][2] = 1; A4[2][3] = 1; A4[3][0] = 1
A4[0][2] = 1; A4[1][3] = 1

n = 4
print(f"\nTournament (n=4): rotational")
for i in range(n):
    out = [j for j in range(n) if A4[i][j]]
    print(f"  {i} → {out}")

# Count 3-cycles
t3 = sum(1 for i,j,k in combinations(range(n), 3)
         if (A4[i][j] and A4[j][k] and A4[k][i]) or
            (A4[i][k] and A4[k][j] and A4[j][i]))
print(f"  t3 = {t3}")

beta = path_betti_numbers(A4, n)
print(f"  β = {[int(b) for b in beta]}")

# Now let's look at the chain groups and boundary maps
print("\nAllowed paths:")
for p in range(n):
    paths = enumerate_allowed_paths(A4, n, p)
    print(f"  Ω_{p}: {len(paths)} paths")
    if len(paths) <= 20:
        for path in paths:
            print(f"    {path}")

# Build and show boundary matrices
print("\nBoundary matrices:")
for p in range(1, min(n, 4)):
    paths_p = enumerate_allowed_paths(A4, n, p)
    paths_pm1 = enumerate_allowed_paths(A4, n, p-1)
    if len(paths_p) == 0 or len(paths_pm1) == 0:
        print(f"  ∂_{p}: empty (0 paths)")
        continue
    mat = build_full_boundary_matrix(A4, n, p)
    print(f"  ∂_{p}: {mat.shape[0]}×{mat.shape[1]}, rank={np.linalg.matrix_rank(mat)}")

# === n=5: Look at a tournament with β₁=1 and one with β₁=0 ===
print("\n" + "=" * 70)
print("n=5: COMPARING β₁=1 vs β₁=0")
print("=" * 70)

# Tournament with t3=5 (the regular C5, β₁=1)
A5_reg = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        A5_reg[i][(i+d)%5] = 1

# Tournament with t3=0 (transitive, β₁=0)
A5_trans = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(i+1, 5):
        A5_trans[i][j] = 1

for name, A, n in [("Regular C5 (t3=5)", A5_reg, 5), ("Transitive (t3=0)", A5_trans, 5)]:
    print(f"\n{name}:")
    beta = path_betti_numbers(A, n)
    print(f"  β = {[int(b) for b in beta]}")
    print(f"  Chain groups:")
    for p in range(n):
        paths = enumerate_allowed_paths(A, n, p)
        print(f"    Ω_{p}: {len(paths)}")
    print(f"  Boundary ranks:")
    for p in range(1, n):
        paths_p = enumerate_allowed_paths(A, n, p)
        if len(paths_p) == 0:
            print(f"    ∂_{p}: 0 (empty)")
            continue
        mat = build_full_boundary_matrix(A, n, p)
        if mat.size > 0:
            print(f"    ∂_{p}: {mat.shape}, rank={np.linalg.matrix_rank(mat)}")
        else:
            print(f"    ∂_{p}: empty matrix")

# === KEY ANALYSIS: Why β₂ = 0 ===
print("\n" + "=" * 70)
print("WHY β₂ = 0 FOR TOURNAMENTS?")
print("=" * 70)

# β₂ = dim(ker ∂₂) - dim(im ∂₃)
# = dim(ker ∂₂) - rank(∂₃)

# For β₂ > 0 we need a 2-cycle that is NOT a 2-boundary.
# A 2-cycle is a formal sum of 2-paths (triangles) Σ c_T * (v₀,v₁,v₂)
# such that ∂₂(Σ) = 0.

# In a tournament, every 2-path (v₀,v₁,v₂) with v₀→v₁→v₂ exists.
# The boundary is ∂₂(v₀,v₁,v₂) = [v₁→v₂](v₁,v₂) - [v₀→v₂](v₀,v₂) + [v₀→v₁](v₀,v₁)
# For a tournament, ALL of these exist! So:
# ∂₂(v₀,v₁,v₂) = (v₁,v₂) - (v₀,v₂) + (v₀,v₁)

# This means ∂₂ is FULL: every 2-path maps to a non-degenerate boundary.
# The boundary matrix has no zero rows.

# For a general digraph, ∂₂(v₀,v₁,v₂) might SKIP terms if edges are missing.
# But tournaments have ALL edges, so ∂₂ always has all 3 terms.

# This "fullness" of ∂₂ should make ker(∂₂) small, potentially forcing β₂ = 0.

print("""
For a tournament T on n vertices:

∂₂(v₀,v₁,v₂) = (v₁,v₂) - (v₀,v₂) + (v₀,v₁)

ALWAYS has all 3 terms (since T is a complete digraph).
This is because for any a,b with a≠b, either a→b or b→a — so every
1-path (a,b) is an allowed path.

In contrast, for a sparse digraph, ∂₂ may drop terms (when edges are missing),
creating "slack" that allows 2-cycles to form more easily.

The completeness of tournaments makes ∂₂ maximally constrained,
potentially forcing ker(∂₂)/im(∂₃) = 0.

ANALOGY: This is like the simplicial homology of a complete simplicial complex
(which is contractible). Tournaments are "complete as digraphs" (every pair
connected), making their path complex highly constrained.
""")

# Let me verify: for the regular C5 tournament (β₁=1, β₂=0)
# What does the 2-chain group look like?
A = A5_reg
n = 5
paths_2 = enumerate_allowed_paths(A, n, 2)
paths_1 = enumerate_allowed_paths(A, n, 1)
paths_3 = enumerate_allowed_paths(A, n, 3)

print(f"C5 tournament: |Ω₁|={len(paths_1)}, |Ω₂|={len(paths_2)}, |Ω₃|={len(paths_3)}")

mat2 = build_full_boundary_matrix(A, n, 2)
mat3 = build_full_boundary_matrix(A, n, 3)

r2 = np.linalg.matrix_rank(mat2)
r3 = np.linalg.matrix_rank(mat3)
ker2 = mat2.shape[1] - r2 if mat2.size > 0 else 0

print(f"∂₂: {mat2.shape}, rank={r2}, ker={ker2}")
if mat3.size > 0:
    print(f"∂₃: {mat3.shape}, rank={r3}")
print(f"β₂ = ker(∂₂) - rank(∂₃) = {ker2} - {r3} = {ker2 - r3}")

# Check β₁ as well
mat1 = build_full_boundary_matrix(A, n, 1)
r1 = np.linalg.matrix_rank(mat1)
ker1 = mat1.shape[1] - r1 if mat1.size > 0 else 0
print(f"∂₁: {mat1.shape}, rank={r1}, ker={ker1}")
print(f"β₁ = ker(∂₁) - rank(∂₂) = {ker1} - {r2} = {ker1 - r2}")

print("\n" + "=" * 70)
print("COUNTING ARGUMENT")
print("=" * 70)

# For a tournament on n vertices:
# |Ω₀| = n (vertices)
# |Ω₁| = n(n-1) (directed edges = all ordered pairs)
# |Ω₂| = ? (directed 2-paths v₀→v₁→v₂, all distinct)
# |Ω₃| = ? (directed 3-paths)

# In a tournament, |Ω_p| = number of (p+1)-tuples (v₀,...,v_p) of distinct
# vertices with v_i→v_{i+1} for all i.
# This is exactly the number of directed paths of length p.

# |Ω₁| = n(n-1) for tournaments (every ordered pair is a directed edge)
# |Ω₂| = n(n-1)(n-1)/2? No...
# |Ω₂| = #{(a,b,c) distinct : a→b, b→c} = Σ_b outdeg(b)*(n-1-outdeg(b))? No...
# Actually |Ω₂| = #{directed paths of length 2} = Σ_b outdeg(b) * (n - 1 - 1)? No.

# For each middle vertex b: need a→b (indeg(b) choices) and b→c (outdeg(b) choices
# minus a→b→a which isn't length 2 if a=c, but we need a≠c).
# Actually (a,b,c) must be 3 DISTINCT vertices.
# For middle vertex b: a can be any of indeg(b) predecessors,
# c can be any of outdeg(b) successors, minus the case a=c (impossible since a→b and b→c).
# Wait, a is a predecessor (a→b) and c is a successor (b→c).
# If a=c then a→b and b→a, impossible in tournament. So all (a,c) pairs valid.
# |Ω₂| = Σ_b indeg(b) * outdeg(b)

for test_n in [3, 4, 5, 6, 7]:
    # For a regular tournament (outdeg = indeg = (n-1)/2):
    if test_n % 2 == 1:
        omega_2_reg = test_n * ((test_n-1)//2)**2
        omega_1 = test_n * (test_n - 1)
        omega_0 = test_n
        print(f"n={test_n} (regular): |Ω₀|={omega_0}, |Ω₁|={omega_1}, |Ω₂|={omega_2_reg}")
    else:
        # Score could be mixed
        omega_1 = test_n * (test_n - 1)
        omega_0 = test_n
        print(f"n={test_n}: |Ω₀|={omega_0}, |Ω₁|={omega_1}")
