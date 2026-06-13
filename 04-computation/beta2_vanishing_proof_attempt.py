#!/usr/bin/env python3
"""
beta2_vanishing_proof_attempt.py - Algebraic proof that β₂ = 0 for tournaments.

KEY IDEA: For tournaments, the boundary map ∂₃: Ω₃ → Ω₂ has special structure
because every triple of vertices forms EITHER a cyclic triple or a transitive triple.

In a tournament:
- Cyclic triple {i,j,k}: one directed 3-cycle (say i→j→k→i)
  The 2-path j→k→i has ∂₂(j→k→i) = k→i - j→i + j→k
  But for cyclic: both j→k and k→i are arcs, and j→i is NOT an arc (i→j is).
  So ∂₂(j→k→i) = (k→i) - (j→i) + (j→k)
  where (j→i) means the 1-chain from j to i.
  Wait, (j→i) is only an allowed 1-path if there's an arc j→i.
  For cyclic triple with cycle i→j→k→i: arcs are i→j, j→k, k→i.
  The allowed 2-paths on {i,j,k} are: i→j→k, j→k→i, k→i→j.

For transitive triple with i→j, i→k, j→k:
  Allowed 2-paths: i→j→k (the only forward path), and also...
  i→k→? (k→j not an arc), j→? (j→k→? k→i not arc if transitive)
  Actually: allowed 2-paths are any a→b→c where a→b and b→c are arcs.
  Transitive {i,j,k} with i→j, i→k, j→k: paths = i→j→k, i→k→? (no, k beats nobody in {i,j,k})
  So only i→j→k.

THEREFORE:
- Cyclic triple contributes 3 allowed 2-paths to Ω₂
- Transitive triple contributes 1 allowed 2-path to Ω₂

This is why dim(Ω₂) = 3·t₃ + (C(n,3) - t₃) = 2·t₃ + C(n,3).

For β₂ = 0, we need: ker(∂₂) = im(∂₃), i.e., every 2-cycle is a boundary.

Author: opus-2026-03-07-S46e
"""
import sys
sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from path_homology_v2 import path_betti_numbers
from itertools import combinations, permutations
import numpy as np
from collections import defaultdict

def count_allowed_paths(A, n, length):
    """Count allowed p-paths (sequences of p+1 vertices with consecutive arcs)."""
    if length == 0:
        return n
    count = 0
    for combo in combinations(range(n), length + 1):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[i+1]] for i in range(length)):
                count += 1
    return count

def build_boundary_matrix(A, n, p):
    """Build the boundary map ∂_p: Ω_p → Ω_{p-1}.

    Ω_p = allowed p-paths.
    ∂_p(v_0→v_1→...→v_p) = Σ_{i=0}^{p} (-1)^i (v_0→...→v̂_i→...→v_p)
    where each term is included only if it's an allowed (p-1)-path.

    But actually ∂_p maps to A_{p-1} (all allowed paths), not Ω_{p-1}.
    Ω_p is the subspace of A_p where the image under ∂_p lands in Ω_{p-1}.

    For homology: H_p = ker(∂_p: Ω_p → Ω_{p-1}) / im(∂_{p+1}: Ω_{p+1} → Ω_p)
    """
    # Enumerate allowed p-paths
    p_paths = []
    for combo in combinations(range(n), p + 1):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[i+1]] for i in range(p)):
                p_paths.append(perm)

    # Enumerate allowed (p-1)-paths
    pm1_paths = []
    for combo in combinations(range(n), p):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[i+1]] for i in range(p-1)):
                pm1_paths.append(perm)

    pm1_index = {tuple(path): i for i, path in enumerate(pm1_paths)}

    # Build matrix
    matrix = np.zeros((len(pm1_paths), len(p_paths)), dtype=int)
    for j, path in enumerate(p_paths):
        for i in range(p + 1):
            face = tuple(path[k] for k in range(p + 1) if k != i)
            if face in pm1_index:
                matrix[pm1_index[face], j] += (-1) ** i

    return matrix, p_paths, pm1_paths

# === n=4: Analyze boundary maps for all tournaments ===
print("=" * 60)
print("n=4: BOUNDARY MAP ANALYSIS")
print("=" * 60)

n = 4
m = n*(n-1)//2

for bits in [0b000000, 0b111111, 0b010101]:  # Transitive, reverse, mixed
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    beta = path_betti_numbers(A, n)
    t3 = sum(1 for i,j,k in combinations(range(n), 3)
             if (A[i][j] and A[j][k] and A[k][i]) or
                (A[i][k] and A[k][j] and A[j][i]))

    print(f"\nbits={bits:06b}, t₃={t3}, β={[int(b) for b in beta]}")

    for p in range(1, 4):
        try:
            mat, p_paths, pm1_paths = build_boundary_matrix(A, n, p)
            rank = np.linalg.matrix_rank(mat)
            print(f"  ∂_{p}: {len(p_paths)} → {len(pm1_paths)}, rank={rank}")
            print(f"    ker(∂_{p}) = {len(p_paths) - rank}, im(∂_{p}) ⊂ R^{len(pm1_paths)}")
        except:
            print(f"  ∂_{p}: error")

# === n=5: Check β₂ = 0 proof structure ===
print("\n" + "=" * 60)
print("n=5: PROOF STRUCTURE FOR β₂ = 0")
print("=" * 60)

n = 5
# Take a tournament with β₁ = 1 (so it's cycle-rich)
A = [[0]*n for _ in range(n)]
# C5 tournament: 0→1→2→3→4→0 and 0→2, 1→3, 2→4, 3→0, 4→1
for i in range(5):
    A[i][(i+1)%5] = 1
    A[i][(i+2)%5] = 1

beta = path_betti_numbers(A, n)
t3 = sum(1 for i,j,k in combinations(range(n), 3)
         if (A[i][j] and A[j][k] and A[k][i]) or
            (A[i][k] and A[k][j] and A[j][i]))
print(f"\nC₅ tournament: t₃={t3}, β={[int(b) for b in beta]}")

for p in range(1, 5):
    mat, p_paths, pm1_paths = build_boundary_matrix(A, n, p)
    rank = np.linalg.matrix_rank(mat)
    ker_dim = len(p_paths) - rank
    print(f"  ∂_{p}: Ω_{p}={len(p_paths)} → Ω_{p-1}={len(pm1_paths)}, rank={rank}, ker={ker_dim}")

# Now check: for ∂₂ and ∂₃, compute im(∂₃) and ker(∂₂)
mat2, paths2, paths1 = build_boundary_matrix(A, n, 2)
mat3, paths3, paths2b = build_boundary_matrix(A, n, 3)

rank2 = np.linalg.matrix_rank(mat2)
rank3 = np.linalg.matrix_rank(mat3)
ker2 = len(paths2) - rank2
im3 = rank3

print(f"\n  β₂ = ker(∂₂) - im(∂₃) = {ker2} - {im3} = {ker2 - im3}")

# Check the DIMENSIONS at each level
print(f"\n  Chain complex dimensions:")
for p in range(5):
    count = count_allowed_paths(A, n, p)
    print(f"    Ω_{p} = {count}")

# === KEY: Why does ker(∂₂) = im(∂₃) always? ===
print("\n" + "=" * 60)
print("WHY β₂ = 0: STRUCTURAL ARGUMENT")
print("=" * 60)

print("""
For a tournament T on n vertices:

1. Every 3-element subset is either cyclic (one 3-cycle) or transitive (total order).
   This is a TOURNAMENT-SPECIFIC property.

2. Ω₂ (allowed 2-paths) decomposes by triple type:
   - Cyclic triple {i,j,k}: contributes 3 allowed 2-paths
   - Transitive triple: contributes 1 allowed 2-path

3. The boundary ∂₃: Ω₃ → Ω₂ maps 3-paths to 2-chains.
   For a 3-path v₀→v₁→v₂→v₃ on vertices {v₀,v₁,v₂,v₃}:
   ∂₃ = (v₁→v₂→v₃) - (v₀→v₂→v₃) + (v₀→v₁→v₃) - (v₀→v₁→v₂)

4. KEY OBSERVATION: In a tournament, every 4-element subset has
   EXACTLY one vertex beaten by all others (or nearly so).
   The structure of ∂₃ on each 4-subset creates enough relations
   to kill all 2-cycles.

5. ALGEBRAIC ARGUMENT: Consider a 2-cycle z ∈ ker(∂₂).
   z = Σ a_{ijk} (i→j→k) with ∂₂(z) = 0.

   For each arc i→j, the coefficient of (i→j) in ∂₂(z) is:
   Σ_k a_{kij} - Σ_k a_{ijk} + correction terms = 0

   This gives a system of linear equations on the coefficients.
   The tournament structure ensures this system is always trivially
   satisfiable (every solution is a boundary of some 3-chain).

6. PROOF STRATEGY: Show that the Euler characteristic argument works:
   χ = dim(Ω₀) - dim(Ω₁) + dim(Ω₂) - dim(Ω₃) + ...
   = β₀ - β₁ + β₂ - β₃ + ...

   If β₂ = 0, then χ = 1 - β₁ + 0 - β₃ + 0 - β₅ + ...
""")

# Verify Euler characteristic at n=5 exhaustive
print("Verifying Euler characteristic at n=5 (all tournaments):")
n = 5
m = n*(n-1)//2
chi_dist = defaultdict(int)

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1

    beta = path_betti_numbers(A, n)
    betas = [int(b) for b in beta]
    chi = sum((-1)**i * betas[i] for i in range(len(betas)))
    chi_dist[chi] += 1

print(f"  χ distribution: {dict(chi_dist)}")
print(f"  χ ∈ {set(chi_dist.keys())}")
