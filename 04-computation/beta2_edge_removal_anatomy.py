#!/usr/bin/env python3
"""
beta2_edge_removal_anatomy.py — When does removing an edge create β₂>0?

The key to proving β₂=0 for tournaments is understanding WHY completeness
prevents β₂>0. The cleanest approach: show that for every non-tournament
digraph with β₂>0, adding the "missing" edge kills the 2-cycle.

More precisely: if G is a tournament minus one edge (u→v removed), and
β₂(G)>0, then the unfillable 2-cycle z MUST use the "gap" at (u,v).
Adding u→v back creates a 3-path through (u,v) whose boundary fills z.

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
        yield mask, A

def compute_beta2_with_cycle(A, n):
    """Compute β₂ and return the 2-cycle basis if β₂>0."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a2:
        return 0, None, {}
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_Om2 == 0:
        return 0, None, {}
    
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    U2, S2, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    dim_Z2 = dim_Om2 - rk2
    
    if dim_Z2 == 0:
        return 0, None, {'dim_Om2': dim_Om2}
    
    Z2_basis = Vt2[rk2:].T  # in Ω₂ coords
    
    if not a3:
        im3 = 0
    else:
        om3 = compute_omega_basis(A, n, 3, a3, a2)
        dim_Om3 = om3.shape[1] if om3.ndim == 2 else 0
        if dim_Om3 == 0:
            im3 = 0
        else:
            bd3 = build_full_boundary_matrix(a3, a2)
            bd3_om = bd3 @ om3
            bd3_om2, _, _, _ = np.linalg.lstsq(om2, bd3_om, rcond=None)
            bd3_Z2 = Z2_basis.T @ bd3_om2
            im3 = np.linalg.matrix_rank(bd3_Z2, tol=1e-8)
    
    beta2 = dim_Z2 - im3
    
    if beta2 > 0:
        # Get the unfillable 2-cycles (kernel of bd3_Z2^T if it exists)
        # Actually: the quotient Z₂/im(∂₃)
        # Express in A₂ coordinates
        cycle_A2 = om2 @ Z2_basis  # columns are Z₂ basis in A₂ coords
        
        return beta2, {'Z2_A2': cycle_A2, 'a2': a2, 'om2': om2, 'Z2_basis': Z2_basis}, \
               {'dim_Om2': dim_Om2, 'dim_Z2': dim_Z2, 'im3': im3}
    
    return 0, None, {'dim_Om2': dim_Om2, 'dim_Z2': dim_Z2, 'im3': im3}

print("=" * 70)
print("EDGE REMOVAL ANATOMY: WHEN DOES β₂ APPEAR?")
print("=" * 70)

# For n=5, find ALL (tournament, edge) pairs where removal creates β₂>0
n = 5
print(f"\n--- n = {n}: exhaustive edge removal ---")

creation_events = []
for mask, A in all_tournaments(n):
    for u in range(n):
        for v in range(n):
            if u == v or A[u][v] == 0:
                continue
            # Remove edge u→v
            B = [row[:] for row in A]
            B[u][v] = 0
            
            beta2, cycle_data, info = compute_beta2_with_cycle(B, n)
            if beta2 > 0:
                creation_events.append((mask, u, v, A, B, beta2, cycle_data, info))

print(f"  Total creation events: {len(creation_events)}")

# Analyze the structure
if creation_events:
    # What does the unfillable 2-cycle look like?
    print(f"\n  Analyzing unfillable 2-cycles:")
    
    for idx, (mask, u, v, A, B, beta2, cd, info) in enumerate(creation_events[:10]):
        scores_A = tuple(sorted([sum(A[i]) for i in range(n)]))
        scores_B = tuple(sorted([sum(B[i]) for i in range(n)]))
        
        print(f"\n  Event {idx}: remove {u}→{v}, β₂={beta2}")
        print(f"    T scores: {scores_A}, T-e scores: {scores_B}")
        
        # The 2-cycle in A₂ coordinates
        Z2_A2 = cd['Z2_A2']
        a2 = cd['a2']
        
        for cyc_idx in range(Z2_A2.shape[1]):
            z = Z2_A2[:, cyc_idx]
            active = [(a2[i], z[i]) for i in range(len(a2)) if abs(z[i]) > 1e-8]
            
            # Does the 2-cycle involve vertex u or v?
            uses_u = any(u in p for p, _ in active)
            uses_v = any(v in p for p, _ in active)
            uses_uv_pair = any((u in p and v in p) for p, _ in active)
            
            print(f"    Cycle {cyc_idx}: {len(active)} terms, uses_u={uses_u}, uses_v={uses_v}, uses_(u,v)_pair={uses_uv_pair}")
            for p, c in sorted(active, key=lambda x: -abs(x[1]))[:6]:
                print(f"      {p}: {c:.4f}")

    # KEY QUESTION: Does the 2-cycle always involve the missing edge endpoints?
    print(f"\n  Does the unfillable 2-cycle always involve the removed edge?")
    all_involve_both = True
    all_involve_pair = True
    for mask, u, v, A, B, beta2, cd, info in creation_events:
        Z2_A2 = cd['Z2_A2']
        a2 = cd['a2']
        for cyc_idx in range(Z2_A2.shape[1]):
            z = Z2_A2[:, cyc_idx]
            active_paths = [a2[i] for i in range(len(a2)) if abs(z[i]) > 1e-8]
            uses_u = any(u in p for p in active_paths)
            uses_v = any(v in p for p in active_paths)
            uses_pair = any((u in p and v in p) for p in active_paths)
            if not (uses_u and uses_v):
                all_involve_both = False
            if not uses_pair:
                all_involve_pair = False
    
    print(f"    Always involves both u AND v: {all_involve_both}")
    print(f"    Always has a 2-path containing both u and v: {all_involve_pair}")

# n=4 analysis — more detailed since smaller
print(f"\n{'='*70}")
print(f"n = 4: exhaustive edge removal + 2-cycle anatomy")
print("=" * 70)

n = 4
creation_events_4 = []
for mask, A in all_tournaments(n):
    for u in range(n):
        for v in range(n):
            if u == v or A[u][v] == 0:
                continue
            B = [row[:] for row in A]
            B[u][v] = 0
            beta2, cycle_data, info = compute_beta2_with_cycle(B, n)
            if beta2 > 0:
                creation_events_4.append((mask, u, v, A, B, beta2, cycle_data, info))

print(f"  Total creation events at n=4: {len(creation_events_4)}")

if creation_events_4:
    for idx, (mask, u, v, A, B, beta2, cd, info) in enumerate(creation_events_4[:5]):
        print(f"\n  Event {idx}: remove {u}→{v}")
        print(f"    Tournament edges: ", end="")
        print([(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1])
        print(f"    After removal: ", end="")
        print([(i,j) for i in range(n) for j in range(n) if i!=j and B[i][j]==1])
        
        Z2_A2 = cd['Z2_A2']
        a2 = cd['a2']
        for cyc_idx in range(Z2_A2.shape[1]):
            z = Z2_A2[:, cyc_idx]
            active = [(a2[i], z[i]) for i in range(len(a2)) if abs(z[i]) > 1e-8]
            print(f"    2-cycle: {[(p, round(c, 4)) for p, c in active]}")
            
            # In the TOURNAMENT (before removal), what 3-paths have this as boundary?
            a3_T = enumerate_allowed_paths(A, n, 3)
            bd3_T = build_full_boundary_matrix(a3_T, enumerate_allowed_paths(A, n, 2))
            
            # Check: is z a boundary in the tournament?
            a2_T = enumerate_allowed_paths(A, n, 2)
            # z is in B's A₂ coords; need to check if it's a cycle in T
            # First: is z even in T's Ω₂?
            z_vec = np.zeros(len(a2_T))
            for p, c in active:
                if p in a2_T:
                    z_vec[a2_T.index(p)] = c
            
            # Check boundary in T
            bd2_T = build_full_boundary_matrix(a2_T, enumerate_allowed_paths(A, n, 1))
            boundary = bd2_T @ z_vec
            is_cycle_in_T = np.max(np.abs(boundary)) < 1e-8
            print(f"    Is z a cycle in T? {is_cycle_in_T}")
            
            if is_cycle_in_T and a3_T:
                # Find a 3-chain filling z in T
                om3_T = compute_omega_basis(A, n, 3, a3_T, a2_T)
                if om3_T.ndim == 2 and om3_T.shape[1] > 0:
                    bd3_om = bd3_T @ om3_T
                    # Solve bd3_om @ c = z_vec (approximately)
                    c_fill, res, _, _ = np.linalg.lstsq(bd3_om, z_vec, rcond=None)
                    fill_err = np.max(np.abs(bd3_om @ c_fill - z_vec))
                    print(f"    Filling error in T: {fill_err:.2e}")
                    if fill_err < 1e-8:
                        # Which 3-paths fill it?
                        fill_3paths = []
                        fill_vec = om3_T @ c_fill
                        for i in range(len(a3_T)):
                            if abs(fill_vec[i]) > 1e-8:
                                p = a3_T[i]
                                is_dt = A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1
                                fill_3paths.append((p, fill_vec[i], 'DT' if is_dt else 'non-DT'))
                        print(f"    Filling 3-paths: {fill_3paths}")

# THE KEY INSIGHT
print(f"\n{'='*70}")
print("THE MECHANISM: WHY ADDING AN EDGE KILLS THE 2-CYCLE")
print("=" * 70)
print("""
When we remove edge u→v from tournament T:
1. Some 2-paths through (u,v) are lost from A₂
2. Some 3-paths through (u,v) are lost from A₃  
3. The lost 3-paths may have been filling 2-cycles that used those 2-paths
4. If the 2-cycle survives (still in ker ∂₂) but its filling is lost → β₂>0

When we ADD the edge back:
1. New 3-paths appear that go through the edge u→v
2. These 3-paths have boundaries landing in Z₂
3. The new boundaries fill the previously unfillable 2-cycle

THIS IS EXACTLY the DT mechanism:
- Adding u→v creates new DT paths (a,u,v,d) with a→v and u→d
- These DT paths' boundaries cover the "gap" in Z₂

PROOF STRATEGY: 
- For β₂>0, there exists z ∈ Z₂ not in im(∂₃)
- For any pair (u,v) in the tournament, the DT paths using u→v contribute
  boundaries that "cover" part of Z₂
- The UNION of all these contributions covers ALL of Z₂
- This is because every 2-path (a,b,c) has potential "extensions" to DT paths
  (w,a,b,c) with w→b (from tournament completeness) 
  or (a,b,c,w) with b→w (from tournament completeness)
""")

print("Done.")
