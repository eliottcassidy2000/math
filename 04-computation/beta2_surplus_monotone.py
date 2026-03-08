#!/usr/bin/env python3
"""
beta2_surplus_monotone.py — Test surplus connectivity and monotonicity

KEY RESULT from skeleton_tiling_proof.py at n=5:
  ALL 64 tilings reachable from transitive with surplus≥0 at every step!

This means the "surplus ≥ 0" constraint is CONNECTED in the tile-flip graph.
If we can show this at ALL n, then β₂=0 follows by induction on tile flips
starting from the transitive tournament (surplus = C(n-1,4) >> 0).

This script:
1. Tests connectivity at n=6 (exhaustive)
2. Tests connectivity at n=7 (sampled BFS from transitive)
3. Analyzes the WORST-CASE flip paths (minimum surplus encountered)

Author: opus-2026-03-08-S44
"""
import sys, time, os, random
import numpy as np
from collections import Counter, deque, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj_from_bits(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_surplus(A, n):
    """Compute surplus = dim(Ω₃) - dim(Z₂) for tournament A."""
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    if not a2:
        return 0, 0, 0, 0

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 and om2.shape[0] > 0 else 0

    if dim_om2 == 0:
        return 0, 0, 0, 0

    bd2 = build_full_boundary_matrix(a2, a1)
    om1 = compute_omega_basis(A, n, 1, a1, enumerate_allowed_paths(A, n, 0))
    bd2_om = bd2 @ om2
    coords, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)
    S = np.linalg.svd(coords, compute_uv=False)
    rk2 = int(sum(s > 1e-8 for s in S))

    z2 = dim_om2 - rk2

    if not a3:
        return dim_om2, 0, z2, -z2

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 and om3.shape[0] > 0 else 0

    surplus = dim_om3 - z2
    return dim_om2, dim_om3, z2, surplus

def flip_arc(A, n, u, v):
    """Flip arc u→v to v→u. Returns new adjacency matrix."""
    A2 = [row[:] for row in A]
    A2[u][v] = 0
    A2[v][u] = 1
    return A2

print("=" * 70)
print("SURPLUS CONNECTIVITY TEST")
print("=" * 70)

# ===== N=5: EXHAUSTIVE =====
n = 5
m = n*(n-1)//2
total = 1 << m
print(f"\n--- n={n}: exhaustive ({total} tournaments) ---")

t0 = time.time()
surplus_map = {}
for bits in range(total):
    A = build_adj_from_bits(n, bits)
    _, _, _, surp = compute_surplus(A, n)
    surplus_map[bits] = surp

# BFS from transitive (all upper-triangle edges)
trans_bits = sum(1 << idx for idx in range(m))  # all bits set = i→j for i<j
reachable = {trans_bits}
queue = deque([trans_bits])
while queue:
    bits = queue.popleft()
    # Try flipping each arc
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            bits2 = bits ^ (1 << idx)
            if bits2 not in reachable and surplus_map[bits2] >= 0:
                reachable.add(bits2)
                queue.append(bits2)
            idx += 1

print(f"  Time: {time.time()-t0:.1f}s")
print(f"  Reachable: {len(reachable)}/{total} ({100*len(reachable)/total:.1f}%)")
all_nonneg = all(s >= 0 for s in surplus_map.values())
print(f"  All surplus ≥ 0: {all_nonneg}")
print(f"  Connected (all reachable): {len(reachable) == total}")

# ===== N=6: EXHAUSTIVE =====
n = 6
m = n*(n-1)//2
total = 1 << m
print(f"\n--- n={n}: exhaustive ({total} tournaments) ---")

t0 = time.time()
surplus_map6 = {}
for bits in range(total):
    if bits % 5000 == 0 and bits > 0:
        print(f"  ... {bits}/{total} ({time.time()-t0:.0f}s)")
    A = build_adj_from_bits(n, bits)
    _, _, _, surp = compute_surplus(A, n)
    surplus_map6[bits] = surp

print(f"  Computed all in {time.time()-t0:.1f}s")

# Check all surplus ≥ 0
all_nonneg6 = all(s >= 0 for s in surplus_map6.values())
print(f"  All surplus ≥ 0: {all_nonneg6}")
min_surp6 = min(surplus_map6.values())
print(f"  Min surplus: {min_surp6}")

# BFS from transitive
trans_bits6 = sum(1 << idx for idx in range(m))
t0 = time.time()
reachable6 = {trans_bits6}
queue = deque([trans_bits6])
while queue:
    bits = queue.popleft()
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            bits2 = bits ^ (1 << idx)
            if bits2 not in reachable6 and surplus_map6[bits2] >= 0:
                reachable6.add(bits2)
                queue.append(bits2)
            idx += 1

print(f"  BFS time: {time.time()-t0:.1f}s")
print(f"  Reachable: {len(reachable6)}/{total} ({100*len(reachable6)/total:.1f}%)")
print(f"  Connected (all reachable): {len(reachable6) == total}")

# ===== N=7: SAMPLED =====
n = 7
m = n*(n-1)//2
print(f"\n--- n={n}: sampled from transitive ---")

# Start from transitive, do random walks, track min surplus
trans_A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        trans_A[i][j] = 1

# Random walk from transitive
random.seed(42)
n_walks = 200
walk_length = 50
min_surplus_seen = float('inf')
all_surpluses = []

t0 = time.time()
for walk in range(n_walks):
    A = [row[:] for row in trans_A]
    for step in range(walk_length):
        # Pick random arc to flip
        arcs = [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j] == 1]
        u, v = random.choice(arcs)
        A = flip_arc(A, n, u, v)

    _, _, _, surp = compute_surplus(A, n)
    all_surpluses.append(surp)
    if surp < min_surplus_seen:
        min_surplus_seen = surp
        scores = sorted(sum(A[i]) for i in range(n))
        print(f"  Walk {walk}: new min surplus={surp}, scores={scores}")

    if (walk + 1) % 50 == 0:
        print(f"  ... {walk+1}/{n_walks} ({time.time()-t0:.0f}s), min_surplus={min_surplus_seen}")

print(f"\n  n={n}: {n_walks} random walks of length {walk_length}")
print(f"  Min surplus: {min_surplus_seen}")
print(f"  All surplus ≥ 0: {all(s >= 0 for s in all_surpluses)}")
print(f"  Surplus distribution: {dict(Counter(all_surpluses).most_common(10))}")

# ===== CRITICAL: DIRECT PATH ANALYSIS =====
print(f"\n{'='*70}")
print("DIRECT PATH: FLIPPING FROM TRANSITIVE TO RANDOM TOURNAMENT")
print("="*70)

# For n=7, take a random tournament and construct a greedy path from
# transitive that always picks the flip minimizing surplus drop
n = 7
m = n*(n-1)//2

random.seed(123)
for trial in range(5):
    # Generate random target
    target = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                target[i][j] = 1
            else:
                target[j][i] = 1

    # Find differing arcs from transitive
    current = [row[:] for row in trans_A]
    diffs = [(i,j) for i in range(n) for j in range(n) if i != j and current[i][j] != target[i][j] and current[i][j] == 1]

    path_surpluses = []
    _, _, _, surp = compute_surplus(current, n)
    path_surpluses.append(surp)

    for step_num in range(len(diffs)):
        # Try each remaining diff, pick the one with best surplus
        best_surp = -float('inf')
        best_uv = None
        remaining = [(i,j) for i in range(n) for j in range(n) if i != j and current[i][j] != target[i][j] and current[i][j] == 1]
        if not remaining:
            break

        for u, v in remaining:
            A_try = flip_arc(current, n, u, v)
            _, _, _, s = compute_surplus(A_try, n)
            if s > best_surp:
                best_surp = s
                best_uv = (u, v)

        if best_uv:
            current = flip_arc(current, n, best_uv[0], best_uv[1])
            path_surpluses.append(best_surp)

    min_on_path = min(path_surpluses)
    target_scores = sorted(sum(target[i]) for i in range(n))
    print(f"\n  Trial {trial}: target scores={target_scores}")
    print(f"    Path length: {len(path_surpluses)-1} flips")
    print(f"    Min surplus on path: {min_on_path}")
    print(f"    Surplus sequence: {path_surpluses[:15]}{'...' if len(path_surpluses)>15 else ''}")

# ===== MINIMUM SURPLUS ALONG GREEDY PATHS =====
print(f"\n{'='*70}")
print("GREEDY PATH MIN SURPLUS (n=7, 100 random targets)")
print("="*70)

random.seed(42)
min_path_mins = []
t0 = time.time()

for trial in range(100):
    target = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                target[i][j] = 1
            else:
                target[j][i] = 1

    current = [row[:] for row in trans_A]
    path_surpluses = []
    _, _, _, surp = compute_surplus(current, n)
    path_surpluses.append(surp)

    for _ in range(m):
        remaining = [(i,j) for i in range(n) for j in range(n)
                     if i != j and current[i][j] != target[i][j] and current[i][j] == 1]
        if not remaining:
            break
        best_surp = -float('inf')
        best_uv = None
        for u, v in remaining:
            A_try = flip_arc(current, n, u, v)
            _, _, _, s = compute_surplus(A_try, n)
            if s > best_surp:
                best_surp = s
                best_uv = (u, v)
        if best_uv:
            current = flip_arc(current, n, best_uv[0], best_uv[1])
            path_surpluses.append(best_surp)

    min_path_mins.append(min(path_surpluses))

    if (trial + 1) % 25 == 0:
        print(f"  ... {trial+1}/100 ({time.time()-t0:.0f}s)")

print(f"\n  Min-on-path distribution: {dict(Counter(min_path_mins).most_common(10))}")
print(f"  Global minimum of min-on-path: {min(min_path_mins)}")
print(f"  All paths stay ≥ 0: {all(m >= 0 for m in min_path_mins)}")

print("\nDone.")
