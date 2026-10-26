#!/usr/bin/env python3
"""
Quick analysis of even Betti vanishing and boundary map structure.
Author: opus-2026-03-07-S46e
"""
import numpy as np
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, build_full_boundary_matrix
)
from itertools import combinations

# Regular C5 tournament
A = [[0]*5 for _ in range(5)]
for i in range(5):
    for d in [1, 2]:
        A[i][(i+d)%5] = 1
n = 5

print("C5 tournament (regular, t3=5, β₁=1):")
beta = path_betti_numbers(A, n)
print(f"β = {[int(b) for b in beta]}")

print("\nChain groups and boundary maps:")
for p in range(n):
    paths = enumerate_allowed_paths(A, n, p)
    print(f"  |Ω_{p}| = {len(paths)}")

for p in range(1, n):
    paths_p = enumerate_allowed_paths(A, n, p)
    paths_pm1 = enumerate_allowed_paths(A, n, p-1)
    if len(paths_p) == 0:
        print(f"  ∂_{p}: empty")
        continue
    mat = build_full_boundary_matrix(paths_p, paths_pm1)
    rk = np.linalg.matrix_rank(mat)
    ker = mat.shape[1] - rk
    print(f"  ∂_{p}: {mat.shape[0]}×{mat.shape[1]}, rank={rk}, ker(∂_{p})={ker}")

# Now compute β_p = ker(∂_p) - rank(∂_{p+1})
print("\nBetti numbers from chain complex:")
for p in range(n):
    paths_p = enumerate_allowed_paths(A, n, p)
    if p > 0:
        paths_pm1 = enumerate_allowed_paths(A, n, p-1)
        mat_p = build_full_boundary_matrix(paths_p, paths_pm1)
        ker_p = mat_p.shape[1] - np.linalg.matrix_rank(mat_p)
    else:
        ker_p = len(paths_p)  # ker(∂_0) = |Ω_0|

    if p < n-1:
        paths_pp1 = enumerate_allowed_paths(A, n, p+1)
        if len(paths_pp1) > 0:
            mat_pp1 = build_full_boundary_matrix(paths_pp1, paths_p)
            rk_pp1 = np.linalg.matrix_rank(mat_pp1)
        else:
            rk_pp1 = 0
    else:
        rk_pp1 = 0

    print(f"  β_{p} = ker(∂_{p}) - rank(∂_{p+1}) = {ker_p} - {rk_pp1} = {ker_p - rk_pp1}")

# === KEY: Look at WHY ker(∂_2) = rank(∂_3) ===
print("\n" + "=" * 60)
print("KEY: Why β₂ = 0?")
print("=" * 60)

paths_2 = enumerate_allowed_paths(A, n, 2)
paths_1 = enumerate_allowed_paths(A, n, 1)
paths_3 = enumerate_allowed_paths(A, n, 3)

mat_2 = build_full_boundary_matrix(paths_2, paths_1)
mat_3 = build_full_boundary_matrix(paths_3, paths_2)

ker_2 = mat_2.shape[1] - np.linalg.matrix_rank(mat_2)
rk_3 = np.linalg.matrix_rank(mat_3)

print(f"|Ω_2| = {len(paths_2)}, |Ω_3| = {len(paths_3)}")
print(f"ker(∂_2) = {ker_2}, rank(∂_3) = {rk_3}")
print(f"β_2 = {ker_2} - {rk_3} = {ker_2 - rk_3}")

if ker_2 == rk_3:
    print("β₂ = 0: Every 2-cycle IS a 2-boundary (no 'true' 2-holes)")

# Compare with transitive tournament
print("\n" + "=" * 60)
print("TRANSITIVE tournament (t3=0):")
print("=" * 60)
A_t = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(i+1, 5):
        A_t[i][j] = 1

beta = path_betti_numbers(A_t, n)
print(f"β = {[int(b) for b in beta]}")

for p in range(n):
    paths = enumerate_allowed_paths(A_t, n, p)
    print(f"  |Ω_{p}| = {len(paths)}")

# === The crucial counting: |Ω_p| for tournaments ===
print("\n" + "=" * 60)
print("COUNTING |Ω_p| FOR TOURNAMENTS")
print("=" * 60)

# For a tournament T: |Ω_p| = number of directed paths of length p
# = number of (p+1)-tuples of distinct vertices forming a directed path

# For TRANSITIVE tournament: every ordering consistent with the total order
# is a directed path. So |Ω_p| = P(n, p+1) = n!/(n-p-1)!

# For REGULAR tournament: depends on structure.

# KEY OBSERVATION for β_2:
# ∂_2: Ω_2 → Ω_1 has matrix size |Ω_1| × |Ω_2|
# For tournament: |Ω_1| = n(n-1) (ALL ordered pairs are edges)
#
# The boundary ∂_2(a,b,c) = (b,c) - (a,c) + (a,b)
# But WAIT: (a,c) is in Ω_1 iff a→c in T.
# If c→a instead, then (a,c) is NOT in Ω_1; instead (c,a) is.
# So the boundary formula becomes:
# ∂_2(a,b,c) = (b,c) ± (a,c) or (c,a) + (a,b)

# Actually in GLMY, the boundary of (v_0,...,v_p) is:
# Σ_{i=0}^{p} (-1)^i * (v_0,...,v̂_i,...,v_p)
# where the term is INCLUDED only if (v_0,...,v̂_i,...,v_p) is allowed.

# For p=2: ∂_2(a,b,c) where a→b→c:
# i=0: +(b,c) — always allowed since b→c
# i=1: -(a,c) — allowed iff a→c
# i=2: +(a,b) — always allowed since a→b

# So: if a→c (same direction as the skip):
#   ∂_2(a,b,c) = (b,c) - (a,c) + (a,b)  [3 terms]
# If c→a (opposite direction):
#   ∂_2(a,b,c) = (b,c) + (a,b)  [2 terms, middle term DROPS]

print("""
GLMY boundary for 2-paths in a tournament:
  ∂_2(a→b→c) = (b,c) + (a,b)           if c→a  [2 terms]
  ∂_2(a→b→c) = (b,c) - (a,c) + (a,b)   if a→c  [3 terms]

When a→c (the "transitive" case): the 3 vertices form a transitive triple.
When c→a (the "cyclic" case): the 3 vertices form a directed 3-cycle.

For a TRANSITIVE tournament: ALL triples are transitive, so ALL ∂_2 have 3 terms.
For a CYCLE-RICH tournament: many ∂_2 have only 2 terms.

The "missing middle term" when c→a creates MORE FLEXIBILITY in the boundary,
which might explain why cycle-rich tournaments have β₃ > 0.
""")

# Let me verify this for the C5 tournament
print("C5 tournament — boundary terms per 2-path:")
two_terms = 0
three_terms = 0
for path in enumerate_allowed_paths(A, 5, 2):
    a, b, c = path
    if A[a][c]:
        three_terms += 1
    else:
        two_terms += 1

print(f"  2-term (cyclic triple): {two_terms}")
print(f"  3-term (transitive triple): {three_terms}")

print("\nTransitive tournament — boundary terms per 2-path:")
two_terms_t = 0
three_terms_t = 0
for path in enumerate_allowed_paths(A_t, 5, 2):
    a, b, c = path
    if A_t[a][c]:
        three_terms_t += 1
    else:
        two_terms_t += 1

print(f"  2-term (cyclic triple): {two_terms_t}")
print(f"  3-term (transitive triple): {three_terms_t}")

# === n=3: The 3-cycle tournament ===
print("\n" + "=" * 60)
print("n=3: 3-CYCLE TOURNAMENT (β₁ = 1)")
print("=" * 60)

A3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
n3 = 3
for p in range(3):
    paths = enumerate_allowed_paths(A3, n3, p)
    print(f"|Ω_{p}| = {len(paths)}: {paths}")

paths_1 = enumerate_allowed_paths(A3, n3, 1)
paths_0 = enumerate_allowed_paths(A3, n3, 0)
paths_2 = enumerate_allowed_paths(A3, n3, 2)

mat_1 = build_full_boundary_matrix(paths_1, paths_0)
print(f"\n∂_1 matrix ({mat_1.shape}):")
print(mat_1)
print(f"rank = {np.linalg.matrix_rank(mat_1)}")

if len(paths_2) > 0:
    mat_2 = build_full_boundary_matrix(paths_2, paths_1)
    print(f"\n∂_2 matrix ({mat_2.shape}):")
    print(mat_2)
    print(f"rank = {np.linalg.matrix_rank(mat_2)}")

# The 1-cycle is the directed 3-cycle itself: (0,1) + (1,2) + (2,0)
# ∂_1((0,1) + (1,2) + (2,0)) = (1)-(0) + (2)-(1) + (0)-(2) = 0 ✓
# This is a 1-cycle. Is it a 1-boundary? Only if there's a 2-chain with this boundary.
# If |Ω_2| = 0 (no 2-paths in C3), then im(∂_2) = 0, so β₁ = ker(∂_1) = 1. ✓

print("\nThe 3-cycle gives β₁=1 because there are no 2-paths (no transitive triples).")
print("The cycle (0,1)+(1,2)+(2,0) is a 1-cycle that cannot be filled.")
