#!/usr/bin/env python3
"""
Investigate the full Omega structure at n=8 for T_661 (the H-maximizer)
and compare with T_657 (the Paley extension).

Key question: Can we compute I(Omega, 2) for the FULL graph (76 vertices)
using inclusion-exclusion or the clique structure?

Actually, since OCF is PROVED, I(Omega, 2) = H(T) = 661.
What we want to understand is the STRUCTURE of Omega.

Instance: opus-2026-03-05-S9
"""
from itertools import combinations, permutations
from collections import Counter

def count_ham_dp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1][v] for v in range(n))

def find_all_odd_cycles(T, n):
    """Find all directed odd cycles of length 3, 5, 7."""
    cycles = []
    for length in range(3, n + 1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(T[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    min_idx = perm.index(min(perm))
                    canon = perm[min_idx:] + perm[:min_idx]
                    cycles.append(canon)
                    break
    return list(set(cycles))

# T_661 (H-maximizer)
T661 = [
    [0,0,0,0,1,1,1,0],
    [1,0,0,0,0,1,0,1],
    [1,1,0,0,0,0,1,1],
    [1,1,1,0,0,0,0,1],
    [0,1,1,1,0,0,0,0],
    [0,0,1,1,1,0,0,0],
    [0,1,0,1,1,1,0,0],
    [1,0,0,0,1,1,1,0],
]

# T_657 (Paley extension)
T657 = [
    [0,0,1,0,0,1,1,0],
    [1,0,0,0,0,1,0,1],
    [0,1,0,0,0,0,1,1],
    [1,1,1,0,0,0,0,0],
    [1,1,1,1,0,0,0,0],
    [0,0,1,1,1,0,0,1],
    [0,1,0,1,1,1,0,0],
    [1,0,0,1,1,0,1,0],
]

for name, T in [("T_661", T661), ("T_657", T657)]:
    n = 8
    print(f"\n{'='*60}")
    print(f"Omega structure for {name}")
    print(f"{'='*60}")
    
    cycles = find_all_odd_cycles(T, n)
    len_dist = Counter(len(c) for c in cycles)
    print(f"Cycle counts: {dict(sorted(len_dist.items()))}")
    print(f"Total vertices in Omega: {len(cycles)}")
    
    # Build Omega adjacency
    m = len(cycles)
    cycle_sets = [set(c) for c in cycles]
    adj = [[0]*m for _ in range(m)]
    edge_count = 0
    for i in range(m):
        for j in range(i+1, m):
            if cycle_sets[i] & cycle_sets[j]:
                adj[i][j] = adj[j][i] = 1
                edge_count += 1
    
    print(f"Edges in Omega: {edge_count}")
    density = 2 * edge_count / (m * (m-1)) if m > 1 else 0
    print(f"Density: {density:.4f}")
    
    # Degree distribution by cycle length
    for length in sorted(len_dist.keys()):
        indices = [i for i in range(m) if len(cycles[i]) == length]
        degs = [sum(adj[i]) for i in indices]
        print(f"  {length}-cycles: degrees min={min(degs)}, max={max(degs)}, mean={sum(degs)/len(degs):.1f}")
    
    # Cross-type edges
    idx_3 = [i for i in range(m) if len(cycles[i]) == 3]
    idx_5 = [i for i in range(m) if len(cycles[i]) == 5]
    idx_7 = [i for i in range(m) if len(cycles[i]) == 7]
    
    e_33 = sum(adj[i][j] for i in idx_3 for j in idx_3 if i < j)
    e_55 = sum(adj[i][j] for i in idx_5 for j in idx_5 if i < j)
    e_77 = sum(adj[i][j] for i in idx_7 for j in idx_7 if i < j)
    e_35 = sum(adj[i][j] for i in idx_3 for j in idx_5)
    e_37 = sum(adj[i][j] for i in idx_3 for j in idx_7)
    e_57 = sum(adj[i][j] for i in idx_5 for j in idx_7)
    
    print(f"  3-3 edges: {e_33}")
    print(f"  5-5 edges: {e_55}")
    print(f"  7-7 edges: {e_77}")
    print(f"  3-5 edges: {e_35}")
    print(f"  3-7 edges: {e_37}")
    print(f"  5-7 edges: {e_57}")
    
    # Independence number (by greedy + small exact for 3-cycle subgraph)
    # For the 3-cycle subgraph only
    m3 = len(idx_3)
    adj3 = [[adj[idx_3[i]][idx_3[j]] for j in range(m3)] for i in range(m3)]
    max_indep_3 = 0
    for mask in range(1 << m3):
        vv = [i for i in range(m3) if mask & (1 << i)]
        ok = True
        for a in range(len(vv)):
            for b in range(a+1, len(vv)):
                if adj3[vv[a]][vv[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            max_indep_3 = max(max_indep_3, len(vv))
    print(f"  Independence number of Omega_3: {max_indep_3}")
    
    # I(Omega_3, 2) 
    i_omega3 = 0
    for mask in range(1 << m3):
        vv = [i for i in range(m3) if mask & (1 << i)]
        ok = True
        for a in range(len(vv)):
            for b in range(a+1, len(vv)):
                if adj3[vv[a]][vv[b]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            i_omega3 += 2**len(vv)
    print(f"  I(Omega_3, 2) = {i_omega3}")
    
    h = count_ham_dp(T, n)
    print(f"  H(T) = {h}")
    print(f"  Deficit (H - I_3) = {h - i_omega3} ({100*(h-i_omega3)/h:.1f}% from 5/7-cycles)")
    
    # Check for C5 in Omega_3
    has_c5 = False
    for combo in combinations(range(m3), 5):
        sub = [[adj3[combo[i]][combo[j]] for j in range(5)] for i in range(5)]
        # Check if it's a C5: each vertex has degree exactly 2
        degs = [sum(sub[i]) for i in range(5)]
        if all(d == 2 for d in degs):
            # Check connected (C5, not pentagon union)
            # Actually degree 2 on 5 vertices with 5 edges = C5
            edges = sum(degs) // 2
            if edges == 5:
                has_c5 = True
                break
    print(f"  C5 in Omega_3: {has_c5}")
