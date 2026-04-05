#!/usr/bin/env python3
"""
opus-2026-04-05-S25: Tournament dynamics as cellular automaton.

CREATIVE INSIGHT: A tournament on n vertices can be viewed as a state in a
cellular automaton defined on the staircase lattice δ_{n-2}.

1. STATE: Each tile (x,y) with x-y≥2 holds a bit (arc direction). State ∈ {0,1}^m.
2. LOCAL RULE: Single tile flip = one CA step. The "neighborhood" is all tiles
   sharing a vertex with the flipped tile.
3. GLOBAL OBSERVABLE: H(T) = number of Hamiltonian paths.

KEY QUESTION: Can we classify tournaments by their CA-like behavior?
- Wolfram Class I: H → fixed value (transitive: H=1)
- Wolfram Class II: H → periodic orbit (complement pairs)
- Wolfram Class III: chaotic/random H under flips
- Wolfram Class IV: complex/emergent behavior at phase transitions

ALSO: The staircase lattice has natural 1D "scan lines" (rows, columns,
antidiagonals). We can view each scan line as a 1D CA state.

EXPERIMENT: Start from a random tournament. Apply single-tile flips.
Track H at each step. Classify the dynamics.
"""

import random
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import time

def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

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

def count_3cycles(T):
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def flip_arc(T, i, j):
    """Flip arc (i,j) and return new tournament."""
    n = len(T)
    T2 = [row[:] for row in T]
    T2[i][j] = 1 - T2[i][j]
    T2[j][i] = 1 - T2[j][i]
    return T2

def tiles(n):
    """Return all tiles (x,y) with x-y >= 2 (non-base-path arcs)."""
    return [(x, y) for y in range(1, n-1) for x in range(y+2, n+1)]

# ═══════════════════════════════════════════════════════
# EXPERIMENT 1: Random walk dynamics — H trajectory
# ═══════════════════════════════════════════════════════
print("="*70)
print("EXPERIMENT 1: Random tile-flip walk — H trajectory")
print("="*70)

for n in [5, 6, 7]:
    T = random_tournament(n)
    tile_list = [(i, j) for i in range(n) for j in range(i+2, n)]  # skip=2+ arcs
    m = len(tile_list)

    H_trajectory = [hamiltonian_path_count(T)]
    c3_trajectory = [count_3cycles(T)]
    steps = 200

    for step in range(steps):
        # Random tile flip
        i, j = random.choice(tile_list)
        T = flip_arc(T, i, j)
        H_trajectory.append(hamiltonian_path_count(T))
        c3_trajectory.append(count_3cycles(T))

    H_arr = np.array(H_trajectory)
    c3_arr = np.array(c3_trajectory)

    print(f"\nn={n} ({m} tiles, {steps} steps):")
    print(f"  H range: [{H_arr.min()}, {H_arr.max()}], mean={H_arr.mean():.1f}, std={H_arr.std():.1f}")
    print(f"  c3 range: [{c3_arr.min()}, {c3_arr.max()}]")

    # Autocorrelation of H
    H_centered = H_arr - H_arr.mean()
    var = np.var(H_centered)
    if var > 0:
        autocorr_1 = np.mean(H_centered[:-1] * H_centered[1:]) / var
        autocorr_5 = np.mean(H_centered[:-5] * H_centered[5:]) / var
        print(f"  H autocorrelation: lag1={autocorr_1:.3f}, lag5={autocorr_5:.3f}")

    # ΔH distribution
    dH = np.diff(H_arr)
    print(f"  ΔH: mean={dH.mean():.2f}, values={dict(Counter(dH))}")

    # Classification attempt
    if H_arr.std() < 2:
        wolfram_class = "I (fixed/near-fixed)"
    elif autocorr_1 > 0.8:
        wolfram_class = "II (periodic-like)"
    elif autocorr_1 < 0.2:
        wolfram_class = "III (chaotic)"
    else:
        wolfram_class = "IV (complex)"
    print(f"  Classification: {wolfram_class}")

# ═══════════════════════════════════════════════════════
# EXPERIMENT 2: Gradient dynamics — always increase H
# ═══════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("EXPERIMENT 2: Greedy ascent — flip arc that maximizes H")
print("="*70)

for n in [5, 6, 7]:
    results = []
    for trial in range(20):
        T = random_tournament(n)
        H = hamiltonian_path_count(T)
        path = [H]

        for step in range(100):
            # Try all arc flips, pick the one that maximizes H
            best_H = H
            best_flip = None
            for i in range(n):
                for j in range(i+1, n):
                    T2 = flip_arc(T, i, j)
                    H2 = hamiltonian_path_count(T2)
                    if H2 > best_H:
                        best_H = H2
                        best_flip = (i, j)

            if best_flip is None:
                break  # local maximum reached

            T = flip_arc(T, *best_flip)
            H = best_H
            path.append(H)

        results.append((path[-1], len(path)-1, path))

    max_H_vals = [r[0] for r in results]
    steps_vals = [r[1] for r in results]
    print(f"\nn={n}: greedy ascent from 20 random starts")
    print(f"  Max H reached: {Counter(max_H_vals)}")
    print(f"  Steps to max: min={min(steps_vals)}, max={max(steps_vals)}, mean={np.mean(steps_vals):.1f}")

    # Show a few paths
    for r in results[:3]:
        print(f"  Path: {r[2][:15]}{'...' if len(r[2])>15 else ''}")

# ═══════════════════════════════════════════════════════
# EXPERIMENT 3: Scan-line CA — row by row
# ═══════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("EXPERIMENT 3: Staircase scan-line encoding")
print("="*70)

# The staircase has rows (indexed by source vertex) and columns (by target).
# Row r contains tiles (r, y) for y < r-1: the arcs from vertex r to non-adjacent vertices.
# This gives a "1D scan" of the tournament state.

for n in [6, 7]:
    T = random_tournament(n)
    H = hamiltonian_path_count(T)

    # Encode tournament as rows of the staircase
    print(f"\nn={n}, H={H}:")
    print(f"  Staircase encoding (row = source vertex, col = skip distance):")
    for src in range(n):
        row = []
        for tgt in range(n):
            if abs(src - tgt) >= 2:
                row.append(T[src][tgt])
        if row:
            print(f"    vertex {src}: {''.join(map(str, row))}")

    # "Rule" analysis: for each vertex, how does its outgoing arc pattern
    # relate to neighbors' patterns?
    scores = [sum(T[i]) for i in range(n)]
    print(f"  Score sequence: {sorted(scores)}")

    # Local rule test: does flipping arc (i,j) affect neighbors predictably?
    # Neighborhood = arcs sharing vertex i or j
    i, j = 0, 2
    T_orig = [row[:] for row in T]
    H_orig = hamiltonian_path_count(T)
    T_flipped = flip_arc(T_orig, i, j)
    H_flipped = hamiltonian_path_count(T_flipped)
    print(f"  Flip ({i},{j}): H {H_orig} → {H_flipped}, ΔH = {H_flipped - H_orig}")

# ═══════════════════════════════════════════════════════
# EXPERIMENT 4: Ergodic mixing — how fast does random walk mix?
# ═══════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("EXPERIMENT 4: Mixing time of random arc-flip walk")
print("="*70)

for n in [5, 6]:
    # Start from transitive tournament
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            T[i][j] = 1  # all arcs go i→j (transitive)

    H_trajectory = []
    c3_trajectory = []
    steps = 500

    for step in range(steps):
        i = random.randint(0, n-1)
        j = random.randint(0, n-1)
        while j == i:
            j = random.randint(0, n-1)
        if i > j:
            i, j = j, i
        T = flip_arc(T, i, j)
        H_trajectory.append(hamiltonian_path_count(T))
        c3_trajectory.append(count_3cycles(T))

    H_arr = np.array(H_trajectory)
    c3_arr = np.array(c3_trajectory)

    # Detect mixing: when does the running mean stabilize?
    window = 50
    running_mean = np.convolve(H_arr, np.ones(window)/window, mode='valid')

    # Theoretical E[H] = n!/2^{n-1}
    from math import factorial as mfact
    EH = mfact(n) / 2**(n-1)

    print(f"\nn={n} (transitive start, {steps} arc-flip steps):")
    print(f"  E[H] = {EH:.2f}")
    print(f"  Running mean (last 100): {H_arr[-100:].mean():.2f}")
    print(f"  Mixing: first 50 mean={H_arr[:50].mean():.2f}, last 50 mean={H_arr[-50:].mean():.2f}")

    # Mixing time estimate: when does running mean first cross E[H]?
    crossing = None
    for t in range(len(running_mean)):
        if abs(running_mean[t] - EH) < EH * 0.1:
            crossing = t + window
            break
    if crossing:
        print(f"  Approximate mixing time (10% of E[H]): ~{crossing} steps")

# ═══════════════════════════════════════════════════════
# EXPERIMENT 5: "Rule 110" analogy — binary orbit patterns
# ═══════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("EXPERIMENT 5: Tournament evolution rule table")
print("="*70)

# For n=5 (m=6 tiles), the tournament state is a 6-bit string.
# Each tile flip is a single-bit mutation.
# The "rule" is: flip tile → new tournament → new H value.
# This is analogous to a 1D CA where the "rule" maps neighborhood → new state.

n = 5
m = n*(n-1)//2 - (n-1)  # number of non-base-path arcs (tiles)
tile_list = [(i, j) for i in range(n) for j in range(i+2, n)]
print(f"\nn={n}: {m} tiles, {2**m} states (tournaments)")
print(f"  Tiles: {tile_list}")

# For each tournament, compute H and the set of H values reachable in 1 flip
from itertools import product

transition_graph = {}
h_landscape = {}

for bits in range(2**m):
    # Construct tournament from bits
    T = [[0]*n for _ in range(n)]
    # Base path: n-1 → n-2 → ... → 1 → 0 (always present)
    for v in range(n-1):
        T[v+1][v] = 1
    # Tiles
    for k, (i, j) in enumerate(tile_list):
        if (bits >> k) & 1:
            T[j][i] = 1  # backward arc
        else:
            T[i][j] = 1  # forward arc (default)

    H = hamiltonian_path_count(T)
    h_landscape[bits] = H

    # Neighbors (1-flip)
    neighbors = []
    for k in range(m):
        nbr = bits ^ (1 << k)
        neighbors.append(nbr)
    transition_graph[bits] = neighbors

# H landscape statistics
H_vals = list(h_landscape.values())
print(f"\nH landscape ({2**m} states):")
print(f"  H values: {sorted(set(H_vals))}")
print(f"  H distribution: {dict(sorted(Counter(H_vals).items()))}")

# Local maxima (no neighbor has higher H)
local_maxima = []
local_minima = []
for bits, H in h_landscape.items():
    nbr_Hs = [h_landscape[nb] for nb in transition_graph[bits]]
    if all(H >= nh for nh in nbr_Hs):
        local_maxima.append((bits, H))
    if all(H <= nh for nh in nbr_Hs):
        local_minima.append((bits, H))

print(f"\nLocal maxima: {len(local_maxima)}")
for bits, H in sorted(local_maxima, key=lambda x: -x[1])[:5]:
    print(f"  bits={bits:0{m}b}, H={H}")
print(f"\nLocal minima: {len(local_minima)}")
for bits, H in sorted(local_minima, key=lambda x: x[1])[:5]:
    print(f"  bits={bits:0{m}b}, H={H}")

# Basin of attraction analysis
# For each state, follow steepest ascent to find which maximum it reaches
basins = defaultdict(list)
for bits in range(2**m):
    current = bits
    visited = set()
    while current not in visited:
        visited.add(current)
        H = h_landscape[current]
        nbrs = transition_graph[current]
        best = current
        best_H = H
        for nb in nbrs:
            if h_landscape[nb] > best_H:
                best_H = h_landscape[nb]
                best = nb
        if best == current:
            break  # local maximum
        current = best
    basins[current].append(bits)

print(f"\nBasin of attraction sizes (greedy ascent):")
for max_bits, basin in sorted(basins.items(), key=lambda x: -len(x[1])):
    H = h_landscape[max_bits]
    print(f"  max H={H} (bits={max_bits:0{m}b}): basin size = {len(basin)}")

# CA-like property: is the dynamics reversible?
# Each tile flip is its own inverse (involution), so YES — the dynamics is reversible.
# This means the "CA" is a reversible CA (bijective rule).
print(f"\nCA property: dynamics is REVERSIBLE (each flip is self-inverse)")

# Connectivity: is the flip graph connected?
# Yes — the hypercube Q_m is connected.
print(f"Flip graph: hypercube Q_{m}, diameter = {m}, connected = YES")

# Spectral gap of the random walk on the H landscape
# The mixing time is related to the spectral gap of the transition matrix
# For Q_m, the second eigenvalue is 1-2/m, so mixing time ~ m*log(m)
print(f"Mixing time of uniform random walk on Q_{m}: ~{m}*log({m}) = ~{m * int(np.log(m) + 1)} steps")

print("\n\nDONE.")
