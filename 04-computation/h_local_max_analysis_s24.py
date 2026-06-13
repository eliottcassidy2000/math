#!/usr/bin/env python3
"""
h_local_max_analysis_s24.py — opus-2026-04-05-S24

Deep analysis of the H=37 local maximum at n=6.

The H-landscape is unimodal at n=4,5 but NOT at n=6.
H=37 is a local maximum (no single arc flip increases H).
The global max is H=45.

Questions:
1. What tournament classes achieve H=37 and are local maxima?
2. What's the "valley depth" — how far down do you need to go to escape?
3. What's the basin size of each attractor?
4. Does this predict similar behavior at n=7?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import time
import random
import sys

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def adj_to_bits(A, n):
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                bits |= (1 << idx)
            idx += 1
    return bits

def count_ham_paths(A, n):
    m = 1 << n
    dp = [[0]*n for _ in range(m)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(m):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            if not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = m - 1
    return sum(dp[full][v] for v in range(n))

def count_3cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c += 1
    return c

def score_sequence(A, n):
    return tuple(sorted(int(A[i].sum()) for i in range(n)))

def flip_arc(A, i, j):
    B = A.copy()
    B[i][j] = 1 - B[i][j]
    B[j][i] = 1 - B[j][i]
    return B

def canonical_form(A, n):
    best = None
    for perm in permutations(range(n)):
        bits = 0
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if A[perm[i]][perm[j]] == 1:
                    bits |= (1 << idx)
                idx += 1
        if best is None or bits < best:
            best = bits
    return best

def is_self_complementary(A, n):
    Ac = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            if i != j:
                Ac[i][j] = 1 - A[i][j]
    return canonical_form(A, n) == canonical_form(Ac, n)

n = 6
m = n * (n - 1) // 2
total = 2 ** m

print("=" * 70)
print(f"DEEP ANALYSIS OF H=37 LOCAL MAXIMUM AT n={n}")
print("=" * 70)

# Step 1: Find ALL local maxima exhaustively at n=6
print("\nStep 1: Exhaustive local maxima search at n=6...")
t0 = time.time()

local_maxima = {}  # H -> list of (bits, A)
H_distribution = Counter()
all_H = {}

for bits in range(total):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    all_H[bits] = H
    H_distribution[H] += 1

    is_local_max = True
    for i in range(n):
        for j in range(i+1, n):
            B = flip_arc(A, i, j)
            new_H = count_ham_paths(B, n)
            if new_H > H:
                is_local_max = False
                break
        if not is_local_max:
            break

    if is_local_max:
        if H not in local_maxima:
            local_maxima[H] = []
        local_maxima[H].append(bits)

t1 = time.time()
print(f"  Time: {t1-t0:.1f}s")

global_max_H = max(local_maxima.keys())
print(f"\n  Global max H = {global_max_H}")
print(f"  Local maxima by H value:")
for H_val in sorted(local_maxima.keys()):
    count = len(local_maxima[H_val])
    pct = count / total * 100
    label = " (GLOBAL MAX)" if H_val == global_max_H else " (LOCAL MAX)"
    print(f"    H = {H_val}: {count} tournaments ({pct:.2f}%){label}")

# Step 2: Characterize the H=37 local maxima
print("\n\nStep 2: Characterizing H=37 local maxima...")
for H_val in sorted(local_maxima.keys()):
    if H_val == global_max_H:
        continue
    print(f"\n  H = {H_val} local maximum ({len(local_maxima[H_val])} tournaments):")

    # Group by canonical form
    canon_groups = defaultdict(list)
    for bits in local_maxima[H_val]:
        A = bits_to_adj(bits, n)
        cf = canonical_form(A, n)
        canon_groups[cf].append(bits)

    print(f"    Iso classes: {len(canon_groups)}")
    for cf, members in sorted(canon_groups.items(), key=lambda x: -len(x[1])):
        A = bits_to_adj(cf, n)
        sc = score_sequence(A, n)
        c3 = count_3cycles(A, n)
        is_sc = is_self_complementary(A, n)
        print(f"    cf={cf}: {len(members)} labeled tournaments, score={sc}, c3={c3}, SC={is_sc}")

        # What's the "valley depth" from this local max?
        # For each local-max tournament, find the minimum H reachable by single flip
        min_neighbor_H = float('inf')
        neighbor_H_dist = Counter()
        for bits in members[:10]:  # sample
            A = bits_to_adj(bits, n)
            for i in range(n):
                for j in range(i+1, n):
                    B = flip_arc(A, i, j)
                    nH = count_ham_paths(B, n)
                    min_neighbor_H = min(min_neighbor_H, nH)
                    neighbor_H_dist[nH] += 1
        print(f"      Neighbor H distribution: {dict(sorted(neighbor_H_dist.items()))}")
        print(f"      Min neighbor H: {min_neighbor_H}, Valley depth: {H_val - min_neighbor_H}")

# Step 3: Basin of attraction analysis
print("\n\nStep 3: Basin of attraction analysis (gradient ascent from all 32768 tournaments)...")
t0 = time.time()
terminal_H_by_start_H = defaultdict(Counter)
terminal_count = Counter()
max_steps_by_terminal = defaultdict(int)

for bits in range(total):
    A = bits_to_adj(bits, n)
    start_H = all_H[bits]
    current = A.copy()
    H = start_H

    for step in range(50):
        best_delta = 0
        best_ij = None
        for i in range(n):
            for j in range(i+1, n):
                B = flip_arc(current, i, j)
                new_H = count_ham_paths(B, n)
                delta = new_H - H
                if delta > best_delta:
                    best_delta = delta
                    best_ij = (i, j)
        if best_ij is None:
            break
        current = flip_arc(current, best_ij[0], best_ij[1])
        H += best_delta

    terminal_H_by_start_H[start_H][H] += 1
    terminal_count[H] += 1
    max_steps_by_terminal[H] = max(max_steps_by_terminal[H], step + 1)

t1 = time.time()
print(f"  Time: {t1-t0:.1f}s")

print(f"\n  Overall basins:")
for H_val in sorted(terminal_count.keys()):
    count = terminal_count[H_val]
    pct = count / total * 100
    print(f"    H={H_val}: {count} ({pct:.1f}%), max {max_steps_by_terminal[H_val]} steps")

print(f"\n  Where do different starting H values end up?")
for start_H in sorted(terminal_H_by_start_H.keys()):
    terminals = terminal_H_by_start_H[start_H]
    total_from_start = sum(terminals.values())
    if total_from_start < 10:
        continue  # skip rare
    parts = []
    for term_H in sorted(terminals.keys()):
        pct = terminals[term_H] / total_from_start * 100
        parts.append(f"H={term_H}:{pct:.0f}%")
    print(f"    Start H={start_H}: {total_from_start} tournaments -> {', '.join(parts)}")

# Step 4: 2-flip escape analysis
print("\n\nStep 4: Can we escape H=37 with 2 flips?")
for bits in local_maxima.get(37, [])[:5]:
    A = bits_to_adj(bits, n)
    # Try all pairs of arc flips
    can_escape = False
    escape_paths = 0
    for i1 in range(n):
        for j1 in range(i1+1, n):
            B1 = flip_arc(A, i1, j1)
            H1 = count_ham_paths(B1, n)
            # Now from B1, try all single flips
            for i2 in range(n):
                for j2 in range(i2+1, n):
                    B2 = flip_arc(B1, i2, j2)
                    H2 = count_ham_paths(B2, n)
                    if H2 > 37:
                        escape_paths += 1
                        can_escape = True
    print(f"  bits={bits}: can escape with 2 flips: {can_escape}, escape paths: {escape_paths}/{m*m}")

# Step 5: What fraction of ALL arc-flip neighbors of H=37 tournaments have H > 37 via 2 flips?
print("\n\nStep 5: Minimum saddle depth (barrier height to escape H=37 basin)")
for bits in local_maxima.get(37, [])[:3]:
    A = bits_to_adj(bits, n)
    best_2step = 0
    best_path = None
    worst_valley = 37
    for i1 in range(n):
        for j1 in range(i1+1, n):
            B1 = flip_arc(A, i1, j1)
            H1 = count_ham_paths(B1, n)
            for i2 in range(n):
                for j2 in range(i2+1, n):
                    B2 = flip_arc(B1, i2, j2)
                    H2 = count_ham_paths(B2, n)
                    if H2 > 37 and H1 > worst_valley:
                        worst_valley = H1
                    if H2 > best_2step:
                        best_2step = H2
                        best_path = (H1, H2, (i1,j1), (i2,j2))
    print(f"  bits={bits}: best H reachable in 2 steps = {best_2step}")
    if best_path:
        print(f"    via H={best_path[0]} -> H={best_path[1]} (flips {best_path[2]}, {best_path[3]})")
    print(f"    Shallowest valley (highest H1 with H2 > 37): {worst_valley}")
    print(f"    Barrier depth: {37 - worst_valley}")
