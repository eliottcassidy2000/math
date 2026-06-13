#!/usr/bin/env python3
"""
beta2_relative_proof.py — Prove β₂=0 via relative homology

The long exact sequence for the pair (T, T\v):
... → H₂(T\v) → H₂(T) → H₂(T, T\v) → H₁(T\v) → H₁(T) → ...

If β₂(T\v) = 0 (induction) and H₂(T, T\v) = 0, then β₂(T) = 0.

Base: n=3, β₂(T) = 0 for all 3-vertex tournaments (verified).
Step: Given β₂(T') = 0 for all (n-1)-tournaments, show β₂(T) = 0.

The long exact sequence:
H₂(T\v) → H₂(T) → H₂(T,T\v) → H₁(T\v) → H₁(T) → ...

With β₂(T\v) = 0: H₂(T) ↪ H₂(T,T\v).
So β₂(T) = 0 iff β₂(T) ≤ dim(H₂(T,T\v)) = 0.
Actually: β₂(T) = 0 iff the map H₂(T) → H₂(T,T\v) is injective AND H₂(T,T\v) = 0.

More precisely: from 0 → H₂(T) → H₂(T,T\v), we get H₂(T) ↪ H₂(T,T\v).
So β₂(T) ≤ dim(H₂(T,T\v)). If H₂(T,T\v) = 0, then β₂(T) = 0.

So the KEY is: H₂(T,T\v) = 0 for all v.

What is H₂(T,T\v)?
The relative chain complex:
  Ω₃(T)/Ω₃(T\v) → Ω₂(T)/Ω₂(T\v) → Ω₁(T)/Ω₁(T\v) → ...

The relative Ω_p consists of Ω_p elements that USE vertex v.

H₂(T,T\v) = ker(∂₂^rel) / im(∂₃^rel) where these are maps on relative chains.

Let me compute this explicitly and understand WHY it vanishes.

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def compute_relative_homology(A, n, v, p_target=2):
    """Compute H_p(T, T\\v) = relative homology at dimension p."""
    # Get Ω bases for T
    allowed = {}
    omega = {}
    for p in range(p_target + 2):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif allowed[p]:
            omega[p] = compute_omega_basis(A, n, p, allowed[p],
                                            allowed.get(p-1, []))
        else:
            omega[p] = np.zeros((0, 0))
    
    # For T\v: vertex v is removed
    # T\v has n-1 vertices: {0,...,n-1}\{v}
    # Remap vertices: keep same labeling but filter out v
    vertices_sub = [i for i in range(n) if i != v]
    A_sub = [[A[i][j] for j in vertices_sub] for i in vertices_sub]
    n_sub = n - 1
    
    # Get Ω bases for T\v  
    allowed_sub = {}
    omega_sub = {}
    for p in range(p_target + 2):
        allowed_sub[p] = enumerate_allowed_paths(A_sub, n_sub, p)
        if p == 0:
            omega_sub[p] = np.eye(n_sub)
        elif allowed_sub[p]:
            omega_sub[p] = compute_omega_basis(A_sub, n_sub, p, allowed_sub[p],
                                                allowed_sub.get(p-1, []))
        else:
            omega_sub[p] = np.zeros((0, 0))
    
    # Relative chains: Ω_p(T) / Ω_p(T\v)
    # Elements using vertex v
    # We need to identify Ω_p(T\v) inside Ω_p(T)
    
    # Ω_p(T\v) paths (in T\v labeling) map to T paths not using v
    # A p-path in T\v: (a₀,...,a_p) where all aᵢ ∈ vertices_sub
    # maps to the same tuple in T (since we kept labels).
    # But wait — compute_omega_basis uses sequential labeling 0..n_sub-1.
    # We need to map back.
    
    # Actually for the relative homology, the key dimensions are:
    # dim_rel_p = dim(Ω_p(T)) - dim(Ω_p(T\v))
    
    dims_T = {}
    dims_sub = {}
    for p in range(p_target + 2):
        if p == 0:
            dims_T[p] = n
            dims_sub[p] = n_sub
        elif omega[p].ndim == 2 and omega[p].shape[0] > 0:
            dims_T[p] = omega[p].shape[1]
        else:
            dims_T[p] = 0
        if p > 0:
            if omega_sub[p].ndim == 2 and omega_sub[p].shape[0] > 0:
                dims_sub[p] = omega_sub[p].shape[1]
            else:
                dims_sub[p] = 0
    
    # Relative dimensions
    rel_dims = {p: dims_T[p] - dims_sub[p] for p in range(p_target + 2)}
    
    # For the actual homology computation, we need the relative boundary maps
    # This is complex. Let me use rank-nullity directly.
    
    # From the long exact sequence:
    # ... → H₂(T\v) → H₂(T) → H₂(T,T\v) → H₁(T\v) → H₁(T) → ...
    
    # Compute β₂(T) and β₂(T\v)
    # β₂(T)
    bd2_T = build_full_boundary_matrix(allowed[2], allowed[1])
    bd2_T_om = bd2_T @ omega[2] if dims_T[2] > 0 else np.zeros((0,0))
    if bd2_T_om.size > 0:
        S = np.linalg.svd(bd2_T_om, compute_uv=False)
        rk2_T = int(np.sum(np.abs(S) > 1e-8))
    else:
        rk2_T = 0
    
    bd3_T = build_full_boundary_matrix(allowed[3], allowed[2]) if allowed.get(3) else np.zeros((0, len(allowed[2]) if allowed[2] else 0))
    if dims_T[3] > 0 and bd3_T.size > 0:
        bd3_T_om = bd3_T @ omega[3]
        S3 = np.linalg.svd(bd3_T_om, compute_uv=False)
        rk3_T = int(np.sum(np.abs(S3) > 1e-8))
    else:
        rk3_T = 0
    
    ker2_T = dims_T[2] - rk2_T
    beta2_T = ker2_T - rk3_T
    
    return {
        'dims_T': dims_T,
        'dims_sub': dims_sub,
        'rel_dims': rel_dims,
        'beta2_T': beta2_T,
        'rk2_T': rk2_T,
        'rk3_T': rk3_T,
    }

print("=" * 70)
print("RELATIVE HOMOLOGY H₂(T, T\\v)")
print("=" * 70)

# Compute relative dimensions and structure
n = 5
print(f"\n--- n = {n}: relative dimension analysis ---")

rel_dim_dist = Counter()
for idx, A in enumerate(all_tournaments(n)):
    for v in range(n):
        data = compute_relative_homology(A, n, v)
        key = tuple(data['rel_dims'].get(p, 0) for p in range(5))
        rel_dim_dist[key] += 1

print(f"  Relative dimension patterns (Ω₀^rel,...,Ω₄^rel):")
for pattern, count in sorted(rel_dim_dist.items()):
    print(f"    {pattern}: {count}")

# THE KEY COMPUTATION: for n=4, verify β₂(T)=0 directly
# AND compute the relative chain complex explicitly
print(f"\n{'='*70}")
print("n=4: EXPLICIT β₂=0 PROOF")
print("=" * 70)

n = 4
print(f"\nFor n=4, every tournament has 4 vertices.")
print(f"The 3-paths (a,b,c,d) use ALL 4 vertices.")
print(f"So the chain complex is very constrained.\n")

# Classify n=4 tournaments
for idx, A in enumerate(all_tournaments(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if idx >= 8:  # Just show a few
        break
    
    t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if A[i][j] and A[j][k] and A[k][i] or A[i][k] and A[k][j] and A[j][i])
    
    om2 = compute_omega_basis(A, n, 2, a2, a1) if a2 else np.zeros((0,0))
    dim_Om2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0
    
    om3 = compute_omega_basis(A, n, 3, a3, a2) if a3 else np.zeros((0,0))
    dim_Om3 = om3.shape[1] if om3.ndim == 2 and om3.shape[0] > 0 else 0
    
    if dim_Om2 > 0:
        bd2 = build_full_boundary_matrix(a2, a1)
        bd2_om = bd2 @ om2
        S2 = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S2) > 1e-8))
    else:
        rk2 = 0
    
    if dim_Om3 > 0:
        bd3 = build_full_boundary_matrix(a3, a2)
        bd3_om = bd3 @ om3
        S3 = np.linalg.svd(bd3_om, compute_uv=False)
        rk3 = int(np.sum(np.abs(S3) > 1e-8))
    else:
        rk3 = 0
    
    dim_Z2 = dim_Om2 - rk2
    beta2 = dim_Z2 - rk3
    
    # DT paths
    dt_count = sum(1 for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1)
    
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    
    print(f"  T{idx}: scores={scores}, t3={t3}, |A₂|={len(a2)}, |A₃|={len(a3)}, "
          f"dim(Ω₂)={dim_Om2}, dim(Ω₃)={dim_Om3}, rk∂₂={rk2}, rk∂₃={rk3}, "
          f"Z₂={dim_Z2}, β₂={beta2}, DT={dt_count}")

# n=4 case analysis for proof
print(f"\n--- n=4 case analysis ---")
n = 4
case_analysis = Counter()

for A in all_tournaments(n):
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    tt = sum(1 for p in a2 if A[p[0]][p[2]] == 1)
    nt = len(a2) - tt
    
    om2 = compute_omega_basis(A, n, 2, a2, a1 if 'a1' in dir() else enumerate_allowed_paths(A, n, 1)) if a2 else np.zeros((0,0))
    dim_Om2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0
    
    dt = sum(1 for p in a3 if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1)
    
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    case_analysis[(scores, len(a2), tt, nt, len(a3), dt, dim_Om2)] += 1

print(f"  (scores, |A₂|, TT, NT, |A₃|, DT, dim(Ω₂)): count")
for case, count in sorted(case_analysis.items()):
    print(f"    {case}: {count}")

print(f"""
AT n=4:
- Transitive (0,1,2,3): 6 paths, all TT, 1 DT → dim(Ω₂)=6, dim(Z₂)=0 (trivially)
  Actually: |A₂|=6 means all triples are allowed and TT
- 1-cycle (0,1,2,3) with one 3-cycle: 5 paths, 4 TT + 1 NT
  Wait: let me check more carefully.
  
The n=4 case should be provable by exhaustion:
- All 4 isomorphism classes
- For each: verify β₂=0 directly
- This IS the base case for induction
""")

print("Done.")
