#!/usr/bin/env python3
"""
rarity_shape_s20fp.py — Rarity and commonality: the shape of the metagraph
kind-pasteur-2026-03-24-S20fp

pi(C) = tilings(C) / 2^m = H(C) / (|Aut(C)| * 2^m)

The SHAPE of the metagraph is determined by the distribution of pi.
Rare classes (small pi): transitive, high-symmetry tournaments
Common classes (large pi): regular, generic tournaments

Questions:
1. How is pi distributed? (histogram, Gini coefficient, entropy)
2. How does rarity correlate with metagraph position? (degree, centrality)
3. Do rare classes cluster together? Are they central or peripheral?
4. What fraction of the metagraph's "weight" is in the top/bottom 10%?
5. The Lorenz curve: how unequal is the distribution?
"""

import sys
import numpy as np
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  RARITY AND SHAPE OF THE METAGRAPH")
print("  kind-pasteur-2026-03-24-S20fp")
print("=" * 80)

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

for n in [5, 6, 7]:
    t0 = time.time()
    TILES = get_tiles(n)
    m = len(TILES)
    N = n
    VERTS = list(range(n, 0, -1))
    all_perms = list(permutations(range(N)))

    tv = [(VERTS.index(x), VERTS.index(y)) for x, y in TILES]

    def b2a(bits):
        A = [[0]*N for _ in range(N)]
        for k in range(N-1): A[k][k+1] = 1
        for i in range(m):
            xi, yi = tv[i]
            if bits[i] == 0: A[xi][yi] = 1
            else: A[yi][xi] = 1
        return A

    def canonicalize(A):
        best = None
        for p in all_perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(N) for j in range(N))
            if best is None or s < best: best = s
        return best

    def count_hp(A):
        dp = [[0]*N for _ in range(1 << N)]
        for v in range(N): dp[1 << v][v] = 1
        for mask in range(1, 1 << N):
            for v in range(N):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(N):
                    if mask & (1 << u): continue
                    if A[v][u]: dp[mask | (1 << u)][u] += dp[mask][v]
        return sum(dp[(1 << N) - 1])

    def count_aut(A):
        return sum(1 for p in all_perms if all(A[p[i]][p[j]] == A[i][j] for i in range(N) for j in range(N)))

    # Build classes
    canon_map = {}
    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = b2a(bits)
        canon_map[mask] = canonicalize(A)

    classes = sorted(set(canon_map.values()))
    V = len(classes)
    cidx = {cn: i for i, cn in enumerate(classes)}
    class_masks = defaultdict(list)
    for mask, cn in canon_map.items():
        class_masks[cn].append(mask)

    # Properties
    tilings = np.array([len(class_masks[cn]) for cn in classes], dtype=float)
    H_vals = np.zeros(V)
    aut_vals = np.zeros(V)
    for i, cn in enumerate(classes):
        A = b2a([(class_masks[cn][0] >> k) & 1 for k in range(m)])
        H_vals[i] = count_hp(A)
        aut_vals[i] = count_aut(A)

    pi = tilings / tilings.sum()

    # Build adjacency
    W = np.zeros((V, V))
    for mask in range(1 << m):
        i = cidx[canon_map[mask]]
        for wi in range(m):
            j = cidx[canon_map[mask ^ (1 << wi)]]
            W[i, j] += 1
    W_off = W.copy()
    np.fill_diagonal(W_off, 0)
    A_unw = (W_off > 0).astype(float)
    deg = A_unw.sum(axis=1)

    print(f"\n{'#'*60}")
    print(f"  n = {n}, m = {m}, V = {V}, 2^m = {2**m}")
    print(f"{'#'*60}")

    # ================================================================
    # 1. DISTRIBUTION OF PI
    # ================================================================
    print(f"\n  1. DISTRIBUTION OF pi = tilings / 2^m:")
    print(f"    min(pi) = {pi.min():.6f} (rarest class)")
    print(f"    max(pi) = {pi.max():.6f} (most common class)")
    print(f"    ratio max/min = {pi.max()/pi.min():.1f}")
    print(f"    mean(pi) = {pi.mean():.6f} = 1/V = {1/V:.6f}")

    # Entropy
    H_entropy = -np.sum(pi * np.log2(pi))
    max_entropy = np.log2(V)
    print(f"    Entropy = {H_entropy:.4f} bits (max = log2(V) = {max_entropy:.4f})")
    print(f"    Entropy ratio = {H_entropy/max_entropy:.4f} (1.0 = uniform)")

    # Gini coefficient
    sorted_pi = np.sort(pi)
    n_pi = len(sorted_pi)
    gini = 1 - 2 * np.sum(sorted_pi * (n_pi - np.arange(n_pi))) / (n_pi * np.sum(sorted_pi))
    print(f"    Gini coefficient = {gini:.4f} (0 = equal, 1 = maximally unequal)")

    # ================================================================
    # 2. TILINGS DISTRIBUTION BY H
    # ================================================================
    print(f"\n  2. TILINGS BY H VALUE:")
    H_unique = sorted(set(H_vals))
    total_tilings = tilings.sum()

    print(f"    {'H':>4} {'#classes':>8} {'tilings':>10} {'%total':>8} {'cumul%':>8} {'avg_aut':>8}")
    cumul = 0
    for h in H_unique:
        mask_h = H_vals == h
        nc = mask_h.sum()
        t = tilings[mask_h].sum()
        cumul += t
        avg_a = aut_vals[mask_h].mean()
        print(f"    {int(h):4d} {nc:8d} {t:10.0f} {t/total_tilings*100:7.1f}% {cumul/total_tilings*100:7.1f}% {avg_a:8.1f}")

    # ================================================================
    # 3. RARITY vs METAGRAPH POSITION
    # ================================================================
    print(f"\n  3. RARITY vs DEGREE:")
    # Correlation
    corr_pi_deg = np.corrcoef(pi, deg)[0,1]
    print(f"    Correlation(pi, degree) = {corr_pi_deg:.4f}")

    # Rare = bottom 25% of pi, Common = top 25%
    q25 = np.percentile(pi, 25)
    q75 = np.percentile(pi, 75)
    rare = pi <= q25
    common = pi >= q75

    rare_deg = deg[rare].mean() if rare.sum() > 0 else 0
    common_deg = deg[common].mean() if common.sum() > 0 else 0
    print(f"    Rare classes (bottom 25%): avg degree = {rare_deg:.1f}")
    print(f"    Common classes (top 25%): avg degree = {common_deg:.1f}")

    # ================================================================
    # 4. CLUSTERING OF RARE CLASSES
    # ================================================================
    # Do rare classes connect mostly to rare, or to common?
    rare_to_rare = 0
    rare_to_common = 0
    rare_to_mid = 0
    for i in range(V):
        if not rare[i]: continue
        for j in range(V):
            if A_unw[i, j] == 0: continue
            if rare[j]: rare_to_rare += 1
            elif common[j]: rare_to_common += 1
            else: rare_to_mid += 1

    total_rare_edges = rare_to_rare + rare_to_common + rare_to_mid
    if total_rare_edges > 0:
        print(f"\n  4. WHERE DO RARE CLASSES CONNECT?")
        print(f"    Rare -> Rare: {rare_to_rare} ({rare_to_rare/total_rare_edges*100:.1f}%)")
        print(f"    Rare -> Mid:  {rare_to_mid} ({rare_to_mid/total_rare_edges*100:.1f}%)")
        print(f"    Rare -> Common: {rare_to_common} ({rare_to_common/total_rare_edges*100:.1f}%)")

    # ================================================================
    # 5. THE LORENZ CURVE (inequality of tilings)
    # ================================================================
    # Fraction of classes holding fraction of tilings
    sorted_tilings = np.sort(tilings)
    cumulative_classes = np.arange(1, V+1) / V
    cumulative_tilings = np.cumsum(sorted_tilings) / total_tilings

    print(f"\n  5. LORENZ CURVE (inequality):")
    for frac in [0.1, 0.25, 0.5, 0.75, 0.9]:
        idx = int(frac * V) - 1
        if idx >= 0 and idx < V:
            print(f"    Bottom {frac*100:.0f}% of classes hold {cumulative_tilings[idx]*100:.1f}% of tilings")

    # ================================================================
    # 6. THE SHAPE: rare-peripheral vs rare-central
    # ================================================================
    # Betweenness-like: is the transitive class (rarest) central or peripheral?
    from collections import deque
    def bfs_dist(adj, start, V):
        dist = [-1] * V
        dist[start] = 0
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in range(V):
                if adj[u][v] > 0 and dist[v] == -1:
                    dist[v] = dist[u] + 1
                    queue.append(v)
        return dist

    trans_idx = int(np.argmin(H_vals))
    reg_idx = int(np.argmax(H_vals))

    dist_from_trans = bfs_dist(A_unw, trans_idx, V)
    dist_from_reg = bfs_dist(A_unw, reg_idx, V)
    eccentricity_trans = max(d for d in dist_from_trans if d >= 0)
    eccentricity_reg = max(d for d in dist_from_reg if d >= 0)

    # Average distance from transitive vs from regular
    avg_dist_trans = np.mean([d for d in dist_from_trans if d >= 0])
    avg_dist_reg = np.mean([d for d in dist_from_reg if d >= 0])

    print(f"\n  6. CENTRALITY OF RARE vs COMMON:")
    print(f"    Transitive (H={H_vals[trans_idx]:.0f}, pi={pi[trans_idx]:.6f}):")
    print(f"      Eccentricity = {eccentricity_trans}, avg distance = {avg_dist_trans:.2f}")
    print(f"    Regular (H={H_vals[reg_idx]:.0f}, pi={pi[reg_idx]:.6f}):")
    print(f"      Eccentricity = {eccentricity_reg}, avg distance = {avg_dist_reg:.2f}")
    print(f"    Rarest class is {'MORE' if avg_dist_trans > avg_dist_reg else 'LESS'} peripheral than commonest")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
