#!/usr/bin/env python3
"""
spurious_maxima_n6.py — opus-2026-03-13-S67j

INVESTIGATING: Spurious local maxima in H landscape at n=6.

At n=6, there are 720 spurious local maxima (local max with H < global max).
This is the first n where the landscape becomes "rough."

KEY QUESTION: What H values do spurious maxima have? What score sequences?
Are they related to the degree-4 Fourier correction?

Also investigate:
- At what H value are the spurious maxima?
- What score sequences do they have?
- How far (in flip distance) are they from the global max?
- Does steepest ascent get STUCK at spurious maxima?
"""

import numpy as np
from itertools import permutations
import time

def score_seq(A):
    return tuple(sorted(A.sum(axis=1).astype(int)))

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)
N = 2**m

print("=" * 70)
print(f"SPURIOUS LOCAL MAXIMA AT n={n}")
print("=" * 70)

# Precompute H for all tournaments
print(f"Computing H for all {N} tournaments...")
t0 = time.time()

H = np.zeros(N, dtype=int)
for bits in range(N):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(edges):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for i_p in range(n-1):
            if A[perm[i_p]][perm[i_p+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    H[bits] = count

print(f"  Done in {time.time()-t0:.1f}s")

global_max = np.max(H)
print(f"  Global max: H={global_max}")
print(f"  H distribution: {sorted(set(H))}")

# Find all local maxima
local_maxima = []
for bits in range(N):
    is_max = True
    for k in range(m):
        if H[bits ^ (1 << k)] > H[bits]:
            is_max = False
            break
    if is_max:
        local_maxima.append(bits)

spurious = [b for b in local_maxima if H[b] < global_max]
genuine = [b for b in local_maxima if H[b] == global_max]

print(f"\n  Local maxima: {len(local_maxima)}")
print(f"  Genuine (H={global_max}): {len(genuine)}")
print(f"  Spurious: {len(spurious)}")

# Analyze spurious maxima
spurious_H_vals = sorted(set(H[b] for b in spurious))
print(f"\n  Spurious max H values: {spurious_H_vals}")
for h_val in spurious_H_vals:
    count = sum(1 for b in spurious if H[b] == h_val)
    print(f"    H={h_val}: {count} spurious local maxima")

# Score sequences of spurious maxima
print(f"\n  Score sequences of spurious maxima:")
score_counts = {}
for b in spurious:
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(edges):
        if (b >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    ss = score_seq(A)
    key = (ss, int(H[b]))
    score_counts[key] = score_counts.get(key, 0) + 1

for (ss, h_val), count in sorted(score_counts.items()):
    print(f"    scores={ss}, H={h_val}: {count} tournaments")

# Score sequences of genuine maxima
print(f"\n  Score sequences of genuine maxima (H={global_max}):")
score_counts_genuine = {}
for b in genuine:
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(edges):
        if (b >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    ss = score_seq(A)
    score_counts_genuine[ss] = score_counts_genuine.get(ss, 0) + 1

for ss, count in sorted(score_counts_genuine.items()):
    print(f"    scores={ss}: {count} tournaments")

# Flip distance from spurious maxima to nearest genuine maximum
print(f"\n  Flip distance from spurious to nearest genuine max:")
# This is Hamming distance in the bit representation

distances = {}
for b in spurious:
    min_dist = m
    for g in genuine:
        d = bin(b ^ g).count('1')
        min_dist = min(min_dist, d)
    h_val = int(H[b])
    if h_val not in distances:
        distances[h_val] = []
    distances[h_val].append(min_dist)

for h_val in sorted(distances.keys()):
    dists = distances[h_val]
    print(f"    H={h_val}: min_dist={min(dists)}, max_dist={max(dists)}, "
          f"mean={np.mean(dists):.1f}")

# Basin sizes: how many tournaments converge to spurious vs genuine maxima?
print(f"\n  Basin analysis (steepest ascent):")
basins_genuine = 0
basins_spurious = 0
stuck_at_h = {}

for bits in range(N):
    current = bits
    while True:
        best_k = -1
        best_delta = 0
        for k in range(m):
            delta = H[current ^ (1 << k)] - H[current]
            if delta > best_delta:
                best_delta = delta
                best_k = k
        if best_k == -1:
            break
        current = current ^ (1 << best_k)

    if H[current] == global_max:
        basins_genuine += 1
    else:
        basins_spurious += 1
        h_val = int(H[current])
        stuck_at_h[h_val] = stuck_at_h.get(h_val, 0) + 1

print(f"    Basins reaching global max: {basins_genuine}/{N} ({100*basins_genuine/N:.1f}%)")
print(f"    Basins stuck at spurious: {basins_spurious}/{N} ({100*basins_spurious/N:.1f}%)")
for h_val, count in sorted(stuck_at_h.items()):
    print(f"      Stuck at H={h_val}: {count}")

# Random restarts: how many random starts needed to find global max?
import random
random.seed(42)
trials = 1000
found_global = 0
steps_to_max = []
for _ in range(trials):
    bits = random.randint(0, N-1)
    current = bits
    steps = 0
    while True:
        best_k = -1
        best_delta = 0
        for k in range(m):
            delta = H[current ^ (1 << k)] - H[current]
            if delta > best_delta:
                best_delta = delta
                best_k = k
        if best_k == -1:
            break
        current = current ^ (1 << best_k)
        steps += 1

    if H[current] == global_max:
        found_global += 1
        steps_to_max.append(steps)

print(f"\n  Random start success rate: {found_global}/{trials} = {100*found_global/trials:.1f}%")
if steps_to_max:
    print(f"    Mean steps when successful: {np.mean(steps_to_max):.1f}")

# What's special about H=43 (or whatever the spurious max is)?
# Is it the "almost regular" score sequence?
print(f"\n  Analysis of near-regular tournaments at n=6:")
print(f"  (n-1)/2 = 2.5, so no regular tournaments exist at n=6)")

# Check: are the spurious maxima related to "almost regular" tournaments?
# Score sequences closest to regular: (2,2,2,3,3,3) and (1,2,3,3,3,3) etc.

print(f"\n  H values by score sequence:")
score_to_H_vals = {}
for bits in range(N):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(edges):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    ss = score_seq(A)
    if ss not in score_to_H_vals:
        score_to_H_vals[ss] = set()
    score_to_H_vals[ss].add(int(H[bits]))

for ss in sorted(score_to_H_vals.keys()):
    h_vals = sorted(score_to_H_vals[ss])
    var_score = sum((s - 2.5)**2 for s in ss) / n
    print(f"    {ss} (Var={var_score:.2f}): H in {h_vals}")

print("\n\nDONE — spurious_maxima_n6.py complete")
