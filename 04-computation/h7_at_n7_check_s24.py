#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Does H=7 occur at n=7?

We found α₁=3 exists at n=7 (c3=3, c5=c7=0).
H = 1 + 2·3 + 4·α₂ = 7 + 4·α₂.
If α₂=0 → H=7. If α₂≥1 → H≥11.

Need: 3 three-cycles, all pairwise sharing vertices (α₂=0), no 5-cycles or 7-cycles.
"""

import random
from itertools import combinations, permutations
from collections import Counter

def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def count_3cycles(T):
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def hamiltonian_path_count(T):
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def find_all_directed_odd_cycles(T):
    n = len(T)
    cycles = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                path = (v0,) + perm
                ok = True
                for i in range(length):
                    if T[path[i]][path[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    min_v = min(path)
                    mi = path.index(min_v)
                    canon = tuple(path[(mi + j) % length] for j in range(length))
                    cycles.add(canon)
                    break
    return list(cycles)

def independence_polynomial_coeffs(adj):
    k = len(adj)
    if k == 0:
        return [1]
    coeffs = [0] * (k + 1)
    for mask in range(2**k):
        verts = [i for i in range(k) if (mask >> i) & 1]
        ok = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            coeffs[len(verts)] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def build_conflict_graph(cycles):
    k = len(cycles)
    adj = [[False]*k for _ in range(k)]
    vsets = [set(c) for c in cycles]
    for i in range(k):
        for j in range(i+1, k):
            if vsets[i] & vsets[j]:
                adj[i][j] = adj[j][i] = True
    return adj

# ─── Method 1: Random sampling at n=7, filter for c3=3 ───
print("="*60)
print("Sampling n=7 for H=7")
print("="*60)

h_dist_c3_3 = Counter()
alpha2_dist = Counter()
n_found_7 = 0
n_alpha1_3 = 0
n_checked = 0
target = 500  # α₁=3 instances to collect

while n_alpha1_3 < target:
    T = random_tournament(7)
    n_checked += 1
    c3 = count_3cycles(T)
    if c3 == 3:
        cycles = find_all_directed_odd_cycles(T)
        alpha1 = len(cycles)
        if alpha1 == 3:
            n_alpha1_3 += 1
            adj = build_conflict_graph(cycles)
            coeffs = independence_polynomial_coeffs(adj)
            H = sum(c * 2**k for k, c in enumerate(coeffs))
            alpha2 = coeffs[2] if len(coeffs) > 2 else 0
            h_dist_c3_3[H] += 1
            alpha2_dist[alpha2] += 1
            if H == 7:
                n_found_7 += 1
                print(f"  H=7 FOUND! α₂={alpha2}, cycles: {cycles}")
                # Print tournament
                scores = tuple(sorted(sum(T[i]) for i in range(7)))
                print(f"    scores={scores}")
                cycle_vsets = [set(c) for c in cycles]
                print(f"    cycle vertices: {cycle_vsets}")
                intersection = cycle_vsets[0] & cycle_vsets[1] & cycle_vsets[2]
                pairwise = [(cycle_vsets[i] & cycle_vsets[j]) for i in range(3) for j in range(i+1, 3)]
                print(f"    common vertex: {intersection}")
                print(f"    pairwise intersections: {pairwise}")

            if n_alpha1_3 % 50 == 0:
                print(f"  ...{n_alpha1_3}/{target} α₁=3 found in {n_checked} trials, "
                      f"H dist so far: {dict(sorted(h_dist_c3_3.items()))}")

print(f"\nResults: checked {n_checked} tournaments, found {n_alpha1_3} with α₁=3")
print(f"H distribution among α₁=3: {dict(sorted(h_dist_c3_3.items()))}")
print(f"α₂ distribution: {dict(sorted(alpha2_dist.items()))}")
print(f"H=7 found: {n_found_7}")

# ─── Method 2: Constructive — build a tournament with 3 sharing-vertex cycles ───
print(f"\n{'='*60}")
print("Constructive approach: 3 three-cycles sharing a vertex")
print("="*60)

def construct_sharing_vertex(n):
    """Try to construct T on n vertices with 3 three-cycles all sharing vertex 0,
    no other odd cycles."""
    # Cycles: 0→1→2→0, 0→3→4→0, 0→5→6→0
    # Need: T[0][1]=T[1][2]=T[2][0]=1, T[0][3]=T[3][4]=T[4][0]=1, T[0][5]=T[5][6]=T[6][0]=1
    # Other arcs: must not create any other 3-cycles or 5-cycles or 7-cycles

    T = [[0]*n for _ in range(n)]

    # Set cycle arcs
    # 0→1→2→0
    T[0][1] = T[1][2] = T[2][0] = 1
    T[1][0] = T[2][1] = T[0][2] = 0
    # 0→3→4→0
    T[0][3] = T[3][4] = T[4][0] = 1
    T[3][0] = T[4][3] = T[0][4] = 0
    # 0→5→6→0
    T[0][5] = T[5][6] = T[6][0] = 1
    T[5][0] = T[6][5] = T[0][6] = 0

    # Remaining arcs: between {1,2} and {3,4}, between {1,2} and {5,6},
    # between {3,4} and {5,6}
    # We need these to be "almost transitive" to avoid creating extra cycles

    # Try all 2^6 = 64 possibilities for the remaining cross-group arcs
    remaining = [(1,3), (1,4), (2,3), (2,4), (1,5), (1,6), (2,5), (2,6),
                 (3,5), (3,6), (4,5), (4,6)]

    found_configs = []
    for bits in range(2**len(remaining)):
        for k, (i, j) in enumerate(remaining):
            if (bits >> k) & 1:
                T[i][j] = 1
                T[j][i] = 0
            else:
                T[j][i] = 1
                T[i][j] = 0

        # Check c3
        c3 = count_3cycles(T)
        if c3 != 3:
            continue

        # Check full cycle structure
        cycles = find_all_directed_odd_cycles(T)
        alpha1 = len(cycles)

        if alpha1 == 3:
            adj = build_conflict_graph(cycles)
            coeffs = independence_polynomial_coeffs(adj)
            H = sum(c * 2**k for k, c in enumerate(coeffs))
            a2 = coeffs[2] if len(coeffs) > 2 else 0
            found_configs.append((bits, H, a2, [c for c in cycles]))

    print(f"Found {len(found_configs)} configurations with α₁=3")
    for bits, H, a2, cyc in found_configs[:5]:
        cycle_vsets = [set(c) for c in cyc]
        common = cycle_vsets[0] & cycle_vsets[1] & cycle_vsets[2]
        pairwise = [(cycle_vsets[i] & cycle_vsets[j]) for i in range(3) for j in range(i+1, 3)]
        print(f"  bits={bits:012b}, H={H}, α₂={a2}, common={common}, pairwise={pairwise}")

    if found_configs:
        # Detailed analysis of the first one
        bits, H, a2, cyc = found_configs[0]
        for k, (i, j) in enumerate(remaining):
            if (bits >> k) & 1:
                T[i][j] = 1
                T[j][i] = 0
            else:
                T[j][i] = 1
                T[i][j] = 0
        print(f"\nDetailed first config: H={H}")
        H_check = hamiltonian_path_count(T)
        print(f"  H verified by DP: {H_check}")
        scores = tuple(sorted(sum(T[i]) for i in range(7)))
        print(f"  scores: {scores}")

construct_sharing_vertex(7)

# ─── Method 3: Exhaustive search for H=7 at n=7 ───
print(f"\n{'='*60}")
print("Exhaustive search for H=7 at n=7 (DP + sampling)")
print("="*60)

n_h7 = 0
n_total = 0
# Can't do all 2^21 with DP at n=7 (would take too long)
# Instead, sample heavily
for trial in range(200000):
    T = random_tournament(7)
    H = hamiltonian_path_count(T)
    if H == 7:
        n_h7 += 1
        scores = tuple(sorted(sum(T[i]) for i in range(7)))
        print(f"  H=7 at n=7! scores={scores}")
        if n_h7 <= 3:
            cycles = find_all_directed_odd_cycles(T)
            print(f"    cycles: {cycles}")

print(f"\n200K random n=7 tournaments: {n_h7} with H=7")

# Also check H distribution near 7
print(f"\nFull H distribution from 50K samples at n=7:")
h_dist_7 = Counter()
for trial in range(50000):
    T = random_tournament(7)
    H = hamiltonian_path_count(T)
    h_dist_7[H] += 1

small_H = {h: c for h, c in h_dist_7.items() if h <= 30}
print(f"H ≤ 30: {dict(sorted(small_H.items()))}")

# Missing odd values
achieved = set(h_dist_7.keys())
max_H = max(achieved)
missing = sorted(set(range(1, 50, 2)) - achieved)
print(f"Missing odd H ≤ 49 (in 50K sample): {missing}")

print("\n\nDONE.")
