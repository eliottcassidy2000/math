#!/usr/bin/env python3
"""
Efficient Hamiltonian path count for n=21 tournament on F_21.

Uses numpy flat array DP: dp[mask * n + v] = count.
Total memory: 2^21 * 21 * 8 bytes = 352 MB (tight but feasible).

If this is too much, we use popcount-layered approach:
only keep two popcount layers at a time.

kind-pasteur-2026-03-06-S25e
"""

import numpy as np
import time
import sys

POW2 = [1, 2, 4]

def mul(a, b):
    return ((a[0] + POW2[a[1]] * b[0]) % 7, (a[1] + b[1]) % 3)

def inv(a):
    j_inv = (-a[1]) % 3
    i_inv = (-POW2[j_inv] * a[0]) % 7
    return (i_inv, j_inv)

elements = [(i, j) for j in range(3) for i in range(7)]
idx_map = {e: k for k, e in enumerate(elements)}
n = 21

nonid = [e for e in elements if e != (0, 0)]
pairs = []
seen = set()
for g in nonid:
    if g not in seen:
        pairs.append((g, inv(g)))
        seen.add(g)
        seen.add(inv(g))

# Non-normal connection set
S_nonnormal = set()
for k, (g, gi) in enumerate(pairs):
    if k < 5:
        S_nonnormal.add(g)
    else:
        S_nonnormal.add(gi)

# Build adjacency
adj = np.zeros((n, n), dtype=np.int8)
for a in elements:
    for b in elements:
        if a == b:
            continue
        diff = mul(inv(a), b)
        if diff in S_nonnormal:
            adj[idx_map[a]][idx_map[b]] = 1

# Precompute: for each vertex v, the bitmask of its in-neighbors
in_mask = np.zeros(n, dtype=np.int32)
for v in range(n):
    for u in range(n):
        if adj[u][v]:
            in_mask[v] |= (1 << u)

print(f"n={n}, F_21 tournament (non-normal S)")
print(f"Computing H(T) via flat-array DP...")
sys.stdout.flush()

# Popcount-layered DP to save memory
# At each popcount p, we only need dp entries for masks with exactly p bits set
# and pass them to popcount p+1.

# We use dictionaries but more efficiently: iterate only over existing entries.
t0 = time.time()

# Current layer: dict {(mask, last): count}
current = {}
for v in range(n):
    current[(1 << v, v)] = 1

for pop in range(2, n + 1):
    next_layer = {}
    for (mask, last), cnt in current.items():
        if cnt == 0:
            continue
        # Extend path: try all vertices not in mask, that last can reach
        can_reach = in_mask[:]  # vertices that have an edge FROM last
        # Actually we want adj[last][nxt] = 1, i.e., nxt in out-neighbors of last
        for nxt in range(n):
            if (mask >> nxt) & 1:
                continue  # already in path
            if not adj[last][nxt]:
                continue  # no edge
            new_mask = mask | (1 << nxt)
            key = (new_mask, nxt)
            if key in next_layer:
                next_layer[key] += cnt
            else:
                next_layer[key] = cnt

    current = next_layer
    elapsed = time.time() - t0

    if pop <= 6 or pop >= n - 3 or pop % 3 == 0:
        print(f"  pop={pop:2d}: {len(current):10d} entries, {elapsed:.1f}s")
    sys.stdout.flush()

# H = sum over all (full_mask, v)
full_mask = (1 << n) - 1
H = sum(current.get((full_mask, v), 0) for v in range(n))
print(f"\nH(T) = {H}")
print(f"H / {n} = {H / n}")
print(f"H mod {n} = {H % n}")
print(f"Total time: {time.time() - t0:.1f}s")

# Also check: is this the same H as for the normal S?
# Normal S = C1 union C3
C1 = {(1,0), (2,0), (4,0)}
C3 = {(a, 1) for a in range(7)}
S_normal = C1 | C3

adj2 = np.zeros((n, n), dtype=np.int8)
for a in elements:
    for b in elements:
        if a == b:
            continue
        diff = mul(inv(a), b)
        if diff in S_normal:
            adj2[idx_map[a]][idx_map[b]] = 1

current2 = {}
for v in range(n):
    current2[(1 << v, v)] = 1

for pop in range(2, n + 1):
    next_layer = {}
    for (mask, last), cnt in current2.items():
        if cnt == 0:
            continue
        for nxt in range(n):
            if (mask >> nxt) & 1:
                continue
            if not adj2[last][nxt]:
                continue
            new_mask = mask | (1 << nxt)
            key = (new_mask, nxt)
            if key in next_layer:
                next_layer[key] += cnt
            else:
                next_layer[key] = cnt
    current2 = next_layer

H2 = sum(current2.get((full_mask, v), 0) for v in range(n))
print(f"\nH(T_normal) = {H2}")
print(f"Same H? {H == H2}")
print(f"\nTotal time: {time.time() - t0:.1f}s")
print("DONE")
