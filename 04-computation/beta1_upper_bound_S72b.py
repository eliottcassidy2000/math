"""
Prove β₁ ≤ 1 for all tournaments.

Strategy: β₁ = dim(ker ∂₁) - rank(∂₂|Ω₂)
        = (C(n,2) - n + 1) - rank(∂₂|Ω₂)

So β₁ ≤ 1 iff rank(∂₂|Ω₂) ≥ C(n,2) - n.

We know Ω₂ ⊆ span{transitive triples} and
∂₂(a,b,c) = (b,c) - (a,c) + (a,b).

The question: does im(∂₂|Ω₂) always have codimension ≤ 1 in ker(∂₁)?

Equivalently: is the "harmonic cycle space" at most 1-dimensional?
A harmonic 1-cycle f satisfies:
  (1) ∂₁f = 0  (flow conservation)
  (2) f(a,b) + f(b,c) = f(a,c) for every transitive triple

opus-2026-03-15-S72b
"""
import numpy as np
from itertools import combinations
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import path_betti_numbers

def tournament_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_beta1_details(A, n):
    """Return β₁ and harmonic cycle space details."""
    # Edge list
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edge_idx[(i,j)] = len(edges)
                edges.append((i,j))
    ne = len(edges)  # = C(n,2) for tournament
    
    # Build ∂₁: edges → vertices
    d1 = np.zeros((n, ne), dtype=float)
    for k, (i,j) in enumerate(edges):
        d1[j][k] += 1   # ∂₁(i,j) = j - i
        d1[i][k] -= 1
    
    # Transitive triples: (a,b,c) with a→b, b→c, a→c
    triples = []
    for a in range(n):
        for b in range(n):
            if not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    triples.append((a,b,c))
    
    # Build ∂₂: triples → edges
    # ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
    nt = len(triples)
    d2 = np.zeros((ne, nt), dtype=float)
    for k, (a,b,c) in enumerate(triples):
        d2[edge_idx[(a,b)]][k] += 1
        d2[edge_idx[(a,c)]][k] -= 1
        d2[edge_idx[(b,c)]][k] += 1
    
    # Actually need Ω₂ (regular 2-paths), not just all triples
    # For tournaments, transitive triples ARE Ω₂ since all faces are edges
    # Wait — need to check the regularity condition more carefully
    # Ω₂ = {c ∈ span(A₂) : all faces of any monomial in c are in A₁}
    # For a single triple (a,b,c): faces are (b,c), (a,c), (a,b)
    # All must be edges in T. Since a→b, a→c, b→c, all three are edges. ✓
    # But Ω₂ also requires that the "junk" (non-regular) part is quotiented out.
    # Junk = elements where ∂ maps to non-allowed paths
    # For degree 2: junk triple = (a,b,c) where face (a,c) is NOT an allowed 1-path
    # i.e., neither a→c nor c→a. But in tournament, exactly one holds.
    # If a→c: (a,b,c) is a transitive triple (already counted)
    # If c→a: (a,b,c) is a 3-cycle path, face (a,c) is c→a (wrong direction)
    #   Wait: face₁(a,b,c) = (a,c) in the sense of the path (a,c).
    #   An allowed 1-path (a,c) requires a→c. If c→a, then (a,c) is NOT allowed.
    #   So (a,b,c) with c→a is a "junk" triple.
    #
    # Ω₂ = span{transitive triples} ∩ ker(junk projection)
    # But actually Ω₂ is defined as the set of 2-chains whose boundary 
    # faces are all allowed. Triples with a→b→c, c→a would have 
    # face (a,c) not allowed (since c→a, not a→c). 
    # These are NOT in Ω₂.
    #
    # So Ω₂ = span{transitive triples (a,b,c) : a→b→c AND a→c}
    # = span{transitive triples in the tournament sense}
    # Which is exactly what we computed! ✓
    
    # But wait, Ω₂ might be a SUBSPACE of span{triples}, not the full span.
    # The regularity condition might impose additional constraints.
    # Let me check with the actual path_homology code.
    
    # Use the actual Ω₂ computation from path_homology_v2
    from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis
    
    allowed_1 = enumerate_allowed_paths(A, n, 1)
    allowed_2 = enumerate_allowed_paths(A, n, 2)
    omega_2_basis = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
    
    # Now compute rank of ∂₂ restricted to Ω₂
    # Build boundary matrix for allowed_2 → allowed_1
    from path_homology_v2 import build_full_boundary_matrix
    bd_2 = build_full_boundary_matrix(allowed_2, allowed_1)
    
    # Restrict to Ω₂
    if omega_2_basis.ndim == 2 and omega_2_basis.shape[1] > 0:
        bd2_omega = bd_2 @ omega_2_basis
        rank_d2 = np.linalg.matrix_rank(bd2_omega, tol=1e-8)
    else:
        rank_d2 = 0
    
    # ker(∂₁)
    rank_d1 = np.linalg.matrix_rank(d1, tol=1e-8)
    ker_d1 = ne - rank_d1  # Should be C(n,2) - n + 1
    
    beta1 = ker_d1 - rank_d2
    
    dim_omega2 = omega_2_basis.shape[1] if omega_2_basis.ndim == 2 else 0
    
    return beta1, ker_d1, rank_d2, dim_omega2, len(triples)

print("="*70)
print("β₁ ≤ 1 INVESTIGATION")
print("="*70)

# Exhaustive check at n=5,6
for n in [5, 6]:
    N = 1 << (n*(n-1)//2)
    max_beta1 = 0
    beta1_dist = {}
    for bits in range(N):
        A = tournament_from_bits(n, bits)
        betti = path_betti_numbers(A, n, max_dim=2)
        b1 = betti[1] if len(betti) > 1 else 0
        max_beta1 = max(max_beta1, b1)
        beta1_dist[b1] = beta1_dist.get(b1, 0) + 1
    print(f"\nn={n}: {N} tournaments")
    print(f"  β₁ distribution: {beta1_dist}")
    print(f"  max β₁ = {max_beta1}")

# Detailed analysis at n=7 (sample)
print(f"\nn=7: sampling 5000 random tournaments...")
import random
random.seed(42)
n = 7
max_beta1 = 0
beta1_dist = {}
for trial in range(5000):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=4)
    b1 = betti[1] if len(betti) > 1 else 0
    max_beta1 = max(max_beta1, b1)
    beta1_dist[b1] = beta1_dist.get(b1, 0) + 1
print(f"  β₁ distribution: {beta1_dist}")
print(f"  max β₁ = {max_beta1}")

# n=8 sample
print(f"\nn=8: sampling 2000 random tournaments...")
n = 8
max_beta1 = 0
beta1_dist = {}
for trial in range(2000):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    b1 = betti[1] if len(betti) > 1 else 0
    max_beta1 = max(max_beta1, b1)
    beta1_dist[b1] = beta1_dist.get(b1, 0) + 1
print(f"  β₁ distribution: {beta1_dist}")
print(f"  max β₁ = {max_beta1}")

# n=9 sample (fewer due to cost)
print(f"\nn=9: sampling 500 random tournaments...")
n = 9
max_beta1 = 0
beta1_dist = {}
for trial in range(500):
    bits = random.randint(0, (1 << (n*(n-1)//2)) - 1)
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    b1 = betti[1] if len(betti) > 1 else 0
    max_beta1 = max(max_beta1, b1)
    beta1_dist[b1] = beta1_dist.get(b1, 0) + 1
print(f"  β₁ distribution: {beta1_dist}")
print(f"  max β₁ = {max_beta1}")

# Now investigate the MECHANISM for β₁ ≤ 1
# Harmonic analysis at n=5
print("\n" + "="*70)
print("HARMONIC CYCLE ANALYSIS (n=5)")
print("="*70)

n = 5
# Pick a β₁=1 tournament
for bits in range(1 << (n*(n-1)//2)):
    A = tournament_from_bits(n, bits)
    betti = path_betti_numbers(A, n, max_dim=2)
    if betti[1] == 1:
        beta1, ker_d1, rank_d2, dim_omega2, n_triples = compute_beta1_details(A, n)
        print(f"\n  bits={bits}: β₁={beta1}")
        print(f"  ker(∂₁)={ker_d1}, rank(∂₂|Ω₂)={rank_d2}, dim(Ω₂)={dim_omega2}")
        print(f"  Transitive triples: {n_triples}")
        print(f"  Codimension of im(∂₂) in ker(∂₁): {ker_d1 - rank_d2}")
        
        # Find the harmonic cycle (nullspace of ∂₂^T restricted to ker ∂₁)
        score = tuple(sorted(A.sum(axis=1).astype(int)))
        print(f"  Score sequence: {score}")
        break

# Check: what is the THEORETICAL lower bound on rank(∂₂)?
# Count transitive triples per vertex
print("\n" + "="*70)
print("TRANSITIVE TRIPLE STRUCTURE")
print("="*70)

for n in [5, 6]:
    N = 1 << (n*(n-1)//2)
    min_triples = float('inf')
    max_triples = 0
    min_omega2 = float('inf')
    max_omega2 = 0
    
    for bits in range(N):
        A = tournament_from_bits(n, bits)
        # Count transitive triples
        count = 0
        for a in range(n):
            for b in range(n):
                if not A[a][b]: continue
                for c in range(n):
                    if c == a or c == b: continue
                    if A[b][c] and A[a][c]:
                        count += 1
        min_triples = min(min_triples, count)
        max_triples = max(max_triples, count)
    
    print(f"\nn={n}: Transitive triple count range: [{min_triples}, {max_triples}]")
    # C(n,3) * 3 (3 orderings per unordered triple that could be transitive)
    # 2*C(n,3) = n(n-1)(n-2)/3 (each unordered triple contributes 2 orderings)
    print(f"  C(n,3) = {n*(n-1)*(n-2)//6}, max possible = {n*(n-1)*(n-2)//3}")

# Key insight: even the tournament with the FEWEST transitive triples
# should still have rank(∂₂) ≥ C(n,2) - n. Let's check.

print("\n" + "="*70)
print("RANK OF ∂₂ vs β₁ STATUS")
print("="*70)

n = 5
N = 1 << (n*(n-1)//2)
for bits in range(min(N, 100)):
    A = tournament_from_bits(n, bits)
    beta1, ker_d1, rank_d2, dim_omega2, n_triples = compute_beta1_details(A, n)
    betti = path_betti_numbers(A, n, max_dim=2)
    b1 = betti[1]
    if bits < 10 or b1 == 1:
        print(f"  bits={bits}: β₁={b1}, rank(∂₂)={rank_d2}, ker(∂₁)={ker_d1}, "
              f"dim(Ω₂)={dim_omega2}, #triples={n_triples}")

print("\n\nDONE")
