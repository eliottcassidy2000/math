#!/usr/bin/env python3
"""
beta2_twin_obstruction.py — The twin vertex obstruction theory

KEY INSIGHT: In oriented graphs, β₂ > 0 requires "parallel structure" —
specifically, pairs of vertices with very similar neighborhoods.
In tournaments, no two vertices can be twins (one must beat the other).

This script:
1. Constructs oriented graphs WITH β₂ > 0
2. Analyzes what structural feature they have that tournaments lack
3. Attempts to formalize the obstruction

The idea: β₂ > 0 means there's a 2-cycle z ∈ Z₂ that can't be filled.
z = Σ aᵢ (xᵢ, yᵢ, zᵢ) with ∂₂(z) = 0.
For z NOT to be a boundary, no 3-chain in Ω₃ can map to z.

In a tournament, between any two vertices, there's an edge in SOME direction.
This means the "missing edges" that create unfillable 2-holes don't exist.

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def compute_betti(A, n, max_dim=4):
    """Compute path Betti numbers."""
    betti = []
    allowed = {}
    omega = {}
    rk = {}
    
    for p in range(max_dim + 1):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif allowed[p]:
            omega[p] = compute_omega_basis(A, n, p, allowed[p],
                                            allowed.get(p-1, []))
        else:
            omega[p] = np.zeros((0, 0))
    
    dims = {}
    for p in range(max_dim + 1):
        if p == 0:
            dims[p] = n
        elif omega[p].ndim == 2 and omega[p].shape[0] > 0:
            dims[p] = omega[p].shape[1]
        else:
            dims[p] = 0
    
    for p in range(max_dim + 1):
        if p == 0 or dims[p] == 0:
            rk[p] = 0
        else:
            bd = build_full_boundary_matrix(allowed[p], allowed.get(p-1, []))
            bd_om = bd @ omega[p]
            S = np.linalg.svd(bd_om, compute_uv=False)
            rk[p] = int(np.sum(np.abs(S) > 1e-8))
    
    betti = []
    for p in range(max_dim + 1):
        ker = dims[p] - rk[p]
        im_next = rk.get(p+1, 0)
        betti.append(ker - im_next)
    
    return betti

# PART 1: Build oriented graphs with β₂ > 0
print("=" * 70)
print("ORIENTED GRAPHS WITH β₂ > 0")
print("=" * 70)

# The simplest known: bidirected complete graph K₃
# A[i][j] = 1 AND A[j][i] = 1 for all i≠j
n = 3
A_bidir = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j:
            A_bidir[i][j] = 1

b = compute_betti(A_bidir, n)
print(f"\nBidirected K₃: β = {b}")
print(f"  β₂ = {b[2]} > 0!")

# K₄ bidirected
n = 4
A_bidir4 = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j:
            A_bidir4[i][j] = 1
b4 = compute_betti(A_bidir4, n)
print(f"Bidirected K₄: β = {b4}")

# What about PARTIAL bidirections? Add one bidirected edge to a tournament
print(f"\n--- Adding one bidirected edge to tournaments (n=4) ---")
n = 4
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

beta2_pos_count = 0
for mask in range(1 << m):
    # Tournament
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    
    # Add bidirected edge 0↔1
    A[0][1] = 1
    A[1][0] = 1
    
    b = compute_betti(A, n)
    if b[2] > 0:
        beta2_pos_count += 1

print(f"  Tournaments + bidirected 0↔1 with β₂>0: {beta2_pos_count}/{1<<m}")

# What about removing one edge from a tournament?
print(f"\n--- Removing one edge from tournaments (n=5) ---")
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

beta2_pos_removal = 0
total_tested = 0
for mask in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (mask >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    
    # Remove edge 0→1 (if it exists)
    if A[0][1] == 1:
        B = [row[:] for row in A]
        B[0][1] = 0
        b = compute_betti(B, n)
        total_tested += 1
        if b[2] > 0:
            beta2_pos_removal += 1
            if beta2_pos_removal <= 3:
                print(f"  β₂>0 after removing 0→1: β={b}")

print(f"  β₂>0 after removing one edge: {beta2_pos_removal}/{total_tested}")

# PART 2: What's the MINIMAL β₂>0 oriented graph on n vertices?
print(f"\n{'='*70}")
print("MINIMAL β₂>0 ORIENTED GRAPHS")
print("=" * 70)

# Try all oriented graphs on 4 vertices (each pair: nothing, i→j, j→i, or both)
n = 4
beta2_examples = []
pairs4 = [(i,j) for i in range(n) for j in range(i+1,n)]

for config in range(3**len(pairs4)):  # 3 options: i→j, j→i, both
    A = [[0]*n for _ in range(n)]
    c = config
    edge_count = 0
    for (i,j) in pairs4:
        mode = c % 3
        c //= 3
        if mode == 0:
            A[i][j] = 1; edge_count += 1
        elif mode == 1:
            A[j][i] = 1; edge_count += 1
        else:  # both
            A[i][j] = 1; A[j][i] = 1; edge_count += 2
    
    b = compute_betti(A, n)
    if b[2] > 0:
        beta2_examples.append((edge_count, A, b))

beta2_examples.sort(key=lambda x: x[0])
print(f"\nTotal oriented graphs on 4 vertices with β₂>0: {len(beta2_examples)}")
if beta2_examples:
    min_edges = beta2_examples[0][0]
    print(f"Minimum edges needed: {min_edges}")
    for e, A, b in beta2_examples[:5]:
        if e == min_edges:
            edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
            bidir = [(i,j) for (i,j) in edges if A[j][i] == 1 and i < j]
            print(f"  edges={len(edges)}, bidirected={bidir}, β={b}")

# PART 3: The key structural difference
print(f"\n{'='*70}")
print("KEY STRUCTURAL ANALYSIS")
print("=" * 70)

# For each vertex pair in a β₂>0 graph, compute "neighborhood similarity"
# In a tournament, N⁺(v) and N⁺(w) differ: one contains the other's edge
print("""
THEOREM SKETCH:
- In a tournament T on n vertices, for any two vertices u,v:
  * Either u→v or v→u (say u→v)
  * Then u ∈ N⁺(v)^c, v ∈ N⁺(u) (or vice versa)
  * So N⁺(u) ≠ N⁺(v) (differ on {u,v})
  * This means no two vertices are "twins" (identical out-neighborhoods)

- For β₂ > 0 in a digraph G, we need a 2-cycle z that can't be filled.
  z = Σ aᵢ (xᵢ, yᵢ, zᵢ) ∈ Z₂ means ∂₂(z) = 0.
  
  ∂₂(x,y,z) = (y,z) - (x,z) + (x,y)
  
  So z being a cycle means the coefficients form a "balanced" pattern.
  
- For z to be a boundary, we need (a,b,c,d) ∈ Ω₃ with ∂₃ mapping to z.
  In a tournament, Ω₃ = DT + cancellation elements.
  DT paths exist whenever there are "transitive 4-tuples" — very common.
  
- The twin obstruction: if u,v are twins in G, then
  any 3-path through u can be "shifted" to go through v,
  but the SIGN of the boundary term may not match, creating
  an unfillable 2-cycle.
""")

# PART 4: Test the "edge completion" theory
# Adding MISSING edges should kill β₂
print(f"\n--- Edge completion killing β₂ ---")
n = 4
kills = 0
total_nontrivial = 0

for config in range(3**len(pairs4)):
    A = [[0]*n for _ in range(n)]
    c = config
    for (i,j) in pairs4:
        mode = c % 3
        c //= 3
        if mode == 0: A[i][j] = 1
        elif mode == 1: A[j][i] = 1
        else: A[i][j] = 1; A[j][i] = 1
    
    b = compute_betti(A, n)
    if b[2] == 0:
        continue
    
    total_nontrivial += 1
    
    # Try adding all missing edges (making it a tournament on the non-bidirected pairs)
    # For each missing pair, try adding an edge
    missing = [(i,j) for i in range(n) for j in range(i+1,n) 
               if A[i][j] == 0 and A[j][i] == 0]
    
    if not missing:
        # No missing edges — this is a "complete" oriented graph (possibly with bidirections)
        # Convert bidirected edges to single direction
        B = [row[:] for row in A]
        for i in range(n):
            for j in range(i+1,n):
                if B[i][j] == 1 and B[j][i] == 1:
                    B[j][i] = 0  # keep i→j, remove j→i
        b2 = compute_betti(B, n)
        if b2[2] == 0:
            kills += 1

print(f"  β₂>0 oriented graphs: {total_nontrivial}")
print(f"  Killed by resolving bidirections: {kills}/{total_nontrivial}")

# More systematic: for EVERY β₂>0 graph, does resolving ALL bidirected edges 
# into tournaments always kill β₂?
print(f"\n--- Resolving ALL bidirected edges kills β₂? ---")
from itertools import product

always_killed = 0
sometimes_killed = 0
never_killed = 0

for config in range(3**len(pairs4)):
    A = [[0]*n for _ in range(n)]
    c = config
    for (i,j) in pairs4:
        mode = c % 3
        c //= 3
        if mode == 0: A[i][j] = 1
        elif mode == 1: A[j][i] = 1
        else: A[i][j] = 1; A[j][i] = 1
    
    b = compute_betti(A, n)
    if b[2] == 0:
        continue
    
    # Find bidirected pairs
    bidir = [(i,j) for i in range(n) for j in range(i+1,n)
             if A[i][j] == 1 and A[j][i] == 1]
    
    if not bidir:
        # No bidirected edges but β₂>0 — this is an incomplete graph
        # Try adding edges (making it complete)
        missing = [(i,j) for i in range(n) for j in range(i+1,n)
                   if A[i][j] == 0 and A[j][i] == 0]
        if missing:
            all_killed = True
            for dirs in product([0,1], repeat=len(missing)):
                B = [row[:] for row in A]
                for k, (i,j) in enumerate(missing):
                    if dirs[k] == 0: B[i][j] = 1
                    else: B[j][i] = 1
                for i,j in bidir:
                    B[j][i] = 0  # resolve bidirections
                b2 = compute_betti(B, n)
                if b2[2] > 0:
                    all_killed = False
                    break
            if all_killed:
                always_killed += 1
            else:
                sometimes_killed += 1
        continue
    
    # Try all ways to resolve bidirected edges
    all_killed = True
    some_killed = False
    for dirs in product([0,1], repeat=len(bidir)):
        B = [row[:] for row in A]
        for k, (i,j) in enumerate(bidir):
            if dirs[k] == 0:
                B[j][i] = 0  # keep i→j
            else:
                B[i][j] = 0  # keep j→i
        # Also complete any missing edges
        missing = [(i,j) for i in range(n) for j in range(i+1,n)
                   if B[i][j] == 0 and B[j][i] == 0]
        if missing:
            for dirs2 in product([0,1], repeat=len(missing)):
                C = [row[:] for row in B]
                for k2, (i,j) in enumerate(missing):
                    if dirs2[k2] == 0: C[i][j] = 1
                    else: C[j][i] = 1
                b2 = compute_betti(C, n)
                if b2[2] == 0:
                    some_killed = True
                else:
                    all_killed = False
        else:
            b2 = compute_betti(B, n)
            if b2[2] == 0:
                some_killed = True
            else:
                all_killed = False
    
    if all_killed:
        always_killed += 1
    elif some_killed:
        sometimes_killed += 1
    else:
        never_killed += 1

print(f"  Always killed by completion to tournament: {always_killed}")
print(f"  Sometimes killed: {sometimes_killed}")
print(f"  NEVER killed: {never_killed}")
if never_killed > 0:
    print(f"  WARNING: some β₂>0 survives tournament completion!")

print("\nDone.")
