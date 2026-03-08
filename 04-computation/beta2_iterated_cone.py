#!/usr/bin/env python3
"""
beta2_iterated_cone.py — Iterated cone approach to β₂=0

Key insight: For a cycle z ∈ Z₂(Ω(T)), decompose z based on the 
first vertex of each path:

z = Σ_a z_a where z_a = sum of terms (a, b, c) in z

For each a, try coning from a vertex v_a with v_a → a.
The error term from the incomplete cone can then be handled by
subsequent cones.

Actually, the SIMPLEST approach might work:
Pick the SINK vertex (in-degree n-1). It beats nobody, so the cone
from the sink goes BACKWARDS: σ * v (append, not prepend).

For a tournament with a Hamiltonian path, the iterated cone along
the path should give a chain homotopy.

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

print("=" * 70)
print("ITERATED CONE: Hamiltonian path approach")
print("=" * 70)

# For a tournament T, use a topological sort (Hamiltonian path).
# Every tournament has a Hamiltonian path: v₀ → v₁ → ... → v_{n-1}.
# 
# Define cone operators:
# C_k: Ω_p(T_k) → Ω_{p+1}(T_{k+1})
# where T_k = T restricted to {v₀,...,v_k}
# C_k(σ) = (v_k, σ) or σ prepended with v_k
#
# The filtration:
# T_0 ⊂ T_1 ⊂ ... ⊂ T_{n-1} = T
# where T_k = tournament on {v₀,...,v_k}
#
# Adding v_k: since v₀→v₁→...→v_k is a Hamiltonian path,
# v_k is beaten by ALL previous vertices (in the HP).
# Wait, not necessarily — the HP is v₀→v₁→...→v_{n-1}, so v_{k-1}→v_k.
# But v₀ might not beat v_k (v_k could beat v₀).
#
# Actually, in the Hamiltonian PATH ordering, v_{i} → v_{i+1} for all i.
# But v_i and v_j for j > i+1 can go either way.

# Better idea: use the TOPOLOGICAL ORDER from a Hamiltonian path.
# In this order, vertex v_0 beats v_1, v_1 beats v_2, etc.
# The SOURCE vertex v_0 has the property that v_0 → v_1.
# But v_0 might not beat all later vertices.

# For the transitive tournament: the topological order is 0,1,2,3,4.
# v_0=0 beats everyone. Cone from 0 works for all.

# For general tournaments: the HP order gives v_0 → v_1 → ... → v_{n-1}.
# Cone from v_0 works for paths starting with v_0's out-neighbors.
# The key is that v_0 → v_1, so we can always cone paths starting with v_1.

# PROOF IDEA:
# 1. Order vertices by Hamiltonian path: v_0 → v_1 → ... → v_{n-1}.
# 2. Filter: T_k = T|_{v_0,...,v_k}.
# 3. At step k, T_k has one more vertex v_k than T_{k-1}.
# 4. The LES for (T_k, T_{k-1}) gives:
#    H_2(T_{k-1}) → H_2(T_k) → H_2(T_k, T_{k-1}) → H_1(T_{k-1}) → H_1(T_k)
# 5. Base case: T_2 is a tournament on 3 vertices, H_2(T_2) = 0 trivially.
# 6. Inductive step: if H_2(T_{k-1}) = 0, is H_2(T_k) = 0?

# This IS the vertex deletion LES approach! And we showed:
# H_2(T_k) = 0 iff H_2(T_k, T_{k-1}) = dim(ker i_*: H_1(T_{k-1}) → H_1(T_k))
# which is what we verified computationally at n=5.

# But the induction needs H_2(T') = 0 for ALL (n-1)-vertex sub-tournaments,
# not just T_{k-1}. So this is circular.

# UNLESS: we can show H_2(T_k, T_{k-1}) = dim(ker i_*) ALGEBRAICALLY,
# without assuming H_2(T_k) = 0.

# From the LES:
# H_2(T_k) = H_2(T_k, T_{k-1}) - dim(ker i_*: H_1(T_{k-1}) → H_1(T_k))
# This holds ASSUMING H_2(T_{k-1}) = 0 (induction).
# So H_2(T_k) = 0 iff H_2(T_k, T_{k-1}) = dim(ker i_*).
# And H_2(T_k, T_{k-1}) ≥ dim(ker i_*) always (from the LES).
# So H_2(T_k) = 0 iff H_2(T_k, T_{k-1}) ≤ dim(ker i_*).
# i.e., H_2(T_k) = 0 iff no "extra" relative 2-cycles exist.

# The question: does the STRUCTURE of the quotient complex R for
# vertex deletion in a tournament guarantee no extra cycles?

# Let me examine R specifically for the HP edge case.

# Pick the regular tournament
A_reg = [[0, 1, 1, 0, 0],
         [0, 0, 1, 1, 0],
         [0, 0, 0, 1, 1],
         [1, 0, 0, 0, 1],
         [1, 1, 0, 0, 0]]

# HP: 0→1→2→3→4→0 (cycle). Actually need to find a Hamiltonian PATH.
# 0→1→2→3→4 works if 3→4: yes (A[3][4]=1). And 2→3: yes. 1→2: yes. 0→1: yes.
# So HP = 0,1,2,3,4.

# Delete vertex 4 from T_4 = T to get T_3 = T|_{0,1,2,3}
print("\nRegular T on {0,1,2,3,4}, delete v=4:")
for i in range(n):
    print(f"  {i} → {[j for j in range(n) if A_reg[i][j]]}")

v = 4
B = [[A_reg[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
print(f"\nT\\{v} on {{0,1,2,3}}:")
for i in range(n-1):
    print(f"  {i} → {[j for j in range(n-1) if B[i][j]]}")

# Now study the quotient R for this specific deletion
# and the connecting map structure

# R₁: edges incident to v=4
# v=4 beats {3} (A[4][3]? No: A[4] = [1,1,0,0,0], so 4 beats 0,1)
# Actually: A_reg[4] = [1, 1, 0, 0, 0], so 4→0, 4→1
# And 2→4 (A[2][4]=1), 3→4 (A[3][4]=1)
# So v=4: out-nbrs = {0,1}, in-nbrs = {2,3}

S_out = [j for j in range(n) if j != v and A_reg[v][j]]
S_in = [j for j in range(n) if j != v and A_reg[j][v]]
print(f"\nv={v}: out-nbrs = {S_out}, in-nbrs = {S_in}")
print(f"R₁ edges: out = {[(v,j) for j in S_out]}, in = {[(j,v) for j in S_in]}")
print(f"dim R₁ = {n-1}")

# H₁(T\v) 
from beta2_h2rel_fixed import compute_betti_1

b1_Tv = compute_betti_1(B, n-1)
b1_T = compute_betti_1(A_reg, n)
print(f"\nβ₁(T) = {b1_T}, β₁(T\\v) = {b1_Tv}")
print(f"dim(ker i_*) ≥ max(0, {b1_Tv} - {b1_T}) = {max(0, b1_Tv - b1_T)}")

# Compute H₂^rel directly
from beta2_h2rel_fixed import compute_h2_rel_fixed

h2r = compute_h2_rel_fixed(A_reg, n, v)
print(f"H₂(T, T\\v) = {h2r}")

# Now try ALL vertices for the regular tournament
print(f"\nAll vertex deletions for regular T:")
for v in range(n):
    B = [[A_reg[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    b1_T = compute_betti_1(A_reg, n)
    b1_Tv = compute_betti_1(B, n-1)
    h2r = compute_h2_rel_fixed(A_reg, n, v)
    ker_i = max(0, b1_Tv - b1_T)  # lower bound
    print(f"  v={v}: β₁(T)={b1_T}, β₁(T\\v)={b1_Tv}, H₂^rel={h2r}, ker(i*)≥{ker_i}")

# For the PROOF: we need to show that H₂(T,T\v) ≤ dim(ker i_*)
# for at least one vertex v.
# 
# Combined with H₂(T,T\v) ≥ dim(ker i_*) (from LES), this gives equality.
# And then β₂(T) = H₂^rel - dim(ker i_*) = 0.
#
# KEY OBSERVATION: try v = source vertex of T.
# The source beats at least ⌈(n-1)/2⌉ vertices.
# For the source: T\v might have HIGHER connectivity,
# potentially making H₂^rel small.

# Actually, the right approach might be the FLAG COMPLEX comparison.
# In GLMY theory, the chain complex of the COMPLETE directed graph K_n
# has all homology trivial (except H₀). If we can show the Ω-complex
# of a tournament is "close enough" to K_n, we might get β₂=0.

# The complete tournament would be: for each pair, BOTH directions present.
# That's not a tournament — it's the complete digraph.
# The complete digraph on n vertices has Ω_p = free on all p-paths,
# and the complex is the normalized chain complex of the n-simplex ≅ pt.
# So H_p(K_n) = 0 for p ≥ 1.

# A tournament is a SUBGRAPH of K_n. The Ω complex of T is a 
# SUBCOMPLEX of the Ω complex of K_n (since A_p(T) ⊆ A_p(K_n)).
# But subcomplexes don't inherit acyclicity.

# HOWEVER: if we can show that Ω(T) is a DEFORMATION RETRACT of Ω(K_n),
# or that the inclusion Ω(T) → Ω(K_n) is a quasi-isomorphism at level 2,
# that would suffice.

# For a tournament: every undirected edge of K_n is present in T 
# (in one direction). So T contains "half" of K_n.
# The "other half" (reversed edges) form the complement tournament T^op.
# K_n = T ∪ T^op (disjoint union of edge sets).

print(f"\n{'='*70}")
print("K_n COMPARISON")
print("=" * 70)

# Compute Ω dimensions for K_5
K5 = [[1 if i != j else 0 for j in range(n)] for i in range(n)]

for p in range(n):
    ap = enumerate_allowed_paths(K5, n, p)
    if p == 0:
        print(f"  |A_{p}(K₅)| = {len(ap)}, dim Ω_{p} = {n}")
    elif ap:
        om = compute_omega_basis(K5, n, p, ap, enumerate_allowed_paths(K5, n, p-1))
        d = om.shape[1] if om.ndim == 2 else 0
        print(f"  |A_{p}(K₅)| = {len(ap)}, dim Ω_{p} = {d}")
    else:
        print(f"  |A_{p}(K₅)| = 0, dim Ω_{p} = 0")

print("\nDone.")
