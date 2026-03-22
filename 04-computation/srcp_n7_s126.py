#!/usr/bin/env python3
"""
srcp_n7_s126.py — SRCP Conjecture Test at n=7 and Beyond
opus-2026-03-21-S126

At n=7: 2^21 = 2,097,152 tournaments, 21 arcs each.
Full enumeration with c5-per-arc is expensive but feasible for
a SAMPLE. We test:

1. Does SRCP with (c3,c5) determine H at n=7? (sampling-based)
2. If not: does adding c7 help?
3. At n=8: sampling-based test with c3-only and (c3,c5).

Strategy: Sample tournaments, compute H and SRCP, check for collisions.
If we find NO collisions in a large sample, the conjecture likely holds.
If we find collisions, we know which extra cycle length is needed.
"""

import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time
import random

print("=" * 72)
print("  SRCP CONJECTURE TEST AT n=7 AND BEYOND")
print("  opus-2026-03-21-S126")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def H_dp(adj, n):
    adj_flat = [0]*(n*n)
    for i in range(n):
        for j in range(n):
            adj_flat[i*n+j] = adj[i][j]
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj_flat[v*n+u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def c3_per_arc(adj, n, i, j):
    """Count 3-cycles through directed arc i->j."""
    count = 0
    for k in range(n):
        if k == i or k == j: continue
        if adj[j][k] and adj[k][i]:
            count += 1
    return count

def c5_per_arc(adj, n, i, j):
    """Count 5-cycles through directed arc i->j."""
    count = 0
    others = [v for v in range(n) if v != i and v != j]
    for a in others:
        if not adj[j][a]: continue
        for b in others:
            if b == a: continue
            if not adj[a][b]: continue
            for c in others:
                if c == a or c == b: continue
                if adj[b][c] and adj[c][i]:
                    count += 1
    return count

def c7_per_arc(adj, n, i, j):
    """Count 7-cycles through directed arc i->j."""
    count = 0
    others = [v for v in range(n) if v != i and v != j]
    for a in others:
        if not adj[j][a]: continue
        for b in others:
            if b == a: continue
            if not adj[a][b]: continue
            for c in others:
                if c == a or c == b: continue
                if not adj[b][c]: continue
                for d in others:
                    if d == a or d == b or d == c: continue
                    if adj[c][d] and adj[d][i]:
                        count += 1
    return count

def compute_srcp(adj, n, include_c5=True, include_c7=False):
    """Compute sorted root cycle profile."""
    profiles = []
    for i in range(n):
        for j in range(n):
            if i == j or adj[i][j] == 0: continue
            c3 = c3_per_arc(adj, n, i, j)
            profile = [c3]
            if include_c5:
                c5 = c5_per_arc(adj, n, i, j)
                profile.append(c5)
            if include_c7:
                c7 = c7_per_arc(adj, n, i, j)
                profile.append(c7)
            profiles.append(tuple(profile))
    return tuple(sorted(profiles))

def random_tournament(n):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

# ==========================================================================
# TEST 1: n=7, c3-only SRCP
# ==========================================================================

n = 7
print(f"TEST 1: n={n}, c3-only SRCP (sampling)")
print("-" * 50)

random.seed(42)
N_samples = 50000

t0 = time.time()
srcp_c3_to_H = defaultdict(set)

for trial in range(N_samples):
    adj = random_tournament(n)
    H = H_dp(adj, n)
    srcp = compute_srcp(adj, n, include_c5=False, include_c7=False)
    srcp_c3_to_H[srcp].add(H)

elapsed = time.time() - t0
n_collisions = sum(1 for h_set in srcp_c3_to_H.values() if len(h_set) > 1)
n_srcps = len(srcp_c3_to_H)

print(f"  {N_samples} samples in {elapsed:.1f}s")
print(f"  Distinct c3-only SRCPs: {n_srcps}")
print(f"  Collisions (multiple H values): {n_collisions}")

if n_collisions > 0:
    print(f"  c3-only DOES NOT determine H at n=7.")
    # Show some collisions
    collision_count = 0
    for srcp, h_set in sorted(srcp_c3_to_H.items(), key=lambda x: -len(x[1])):
        if len(h_set) > 1:
            print(f"    H values: {sorted(h_set)}")
            collision_count += 1
            if collision_count >= 5: break
else:
    print(f"  NO collisions found in {N_samples} samples! c3-only may suffice.")

print()

# ==========================================================================
# TEST 2: n=7, (c3,c5) SRCP
# ==========================================================================

print(f"TEST 2: n={n}, (c3,c5) SRCP (sampling)")
print("-" * 50)

random.seed(42)
N_samples_c5 = 20000  # Fewer because c5 computation is slower

t0 = time.time()
srcp_c3c5_to_H = defaultdict(set)

for trial in range(N_samples_c5):
    adj = random_tournament(n)
    H = H_dp(adj, n)
    srcp = compute_srcp(adj, n, include_c5=True, include_c7=False)
    srcp_c3c5_to_H[srcp].add(H)

    if (trial+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1} samples in {elapsed:.1f}s")

elapsed = time.time() - t0
n_collisions_c5 = sum(1 for h_set in srcp_c3c5_to_H.values() if len(h_set) > 1)
n_srcps_c5 = len(srcp_c3c5_to_H)

print(f"  {N_samples_c5} samples in {elapsed:.1f}s")
print(f"  Distinct (c3,c5) SRCPs: {n_srcps_c5}")
print(f"  Collisions: {n_collisions_c5}")

if n_collisions_c5 > 0:
    print(f"  (c3,c5) DOES NOT determine H at n=7.")
    collision_count = 0
    for srcp, h_set in sorted(srcp_c3c5_to_H.items(), key=lambda x: -len(x[1])):
        if len(h_set) > 1:
            print(f"    H values: {sorted(h_set)}")
            collision_count += 1
            if collision_count >= 5: break

    # Test: does adding c7 help?
    print(f"\n  Testing (c3,c5,c7) on colliding cases...")

    # Recompute with c7 for the collision cases
    random.seed(42)
    srcp_full_to_H = defaultdict(set)
    collision_trials = 0

    for trial in range(N_samples_c5):
        adj = random_tournament(n)
        H = H_dp(adj, n)
        srcp_c3c5 = compute_srcp(adj, n, include_c5=True, include_c7=False)

        if len(srcp_c3c5_to_H.get(srcp_c3c5, set())) > 1:
            # This tournament is in a colliding SRCP
            srcp_full = compute_srcp(adj, n, include_c5=True, include_c7=True)
            srcp_full_to_H[srcp_full].add(H)
            collision_trials += 1

    n_coll_full = sum(1 for h_set in srcp_full_to_H.values() if len(h_set) > 1)
    print(f"  {collision_trials} collision cases tested with c7")
    print(f"  Collisions with (c3,c5,c7): {n_coll_full}")
    if n_coll_full == 0:
        print(f"  *** (c3,c5,c7) RESOLVES all collisions at n=7! ***")
    else:
        print(f"  Still {n_coll_full} collisions even with c7.")
else:
    print(f"  *** NO collisions! (c3,c5) appears to determine H at n=7! ***")

print()

# ==========================================================================
# TEST 3: n=8, c3-only SRCP (quick sample)
# ==========================================================================

n = 8
print(f"TEST 3: n={n}, c3-only SRCP (quick sampling)")
print("-" * 50)

random.seed(42)
N_samples_8 = 10000

t0 = time.time()
srcp_c3_to_H_8 = defaultdict(set)

for trial in range(N_samples_8):
    adj = random_tournament(n)
    H = H_dp(adj, n)
    srcp = compute_srcp(adj, n, include_c5=False, include_c7=False)
    srcp_c3_to_H_8[srcp].add(H)

elapsed = time.time() - t0
n_collisions_8 = sum(1 for h_set in srcp_c3_to_H_8.values() if len(h_set) > 1)
n_srcps_8 = len(srcp_c3_to_H_8)

print(f"  {N_samples_8} samples in {elapsed:.1f}s")
print(f"  Distinct c3-only SRCPs: {n_srcps_8}")
print(f"  Collisions: {n_collisions_8}")
print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print("=" * 72)
print("  SUMMARY: WHICH CYCLE LENGTHS ARE NEEDED?")
print("=" * 72)
print()

print(f"  n=5: (c3,c5) determines H  [EXACT, exhaustive]")
print(f"  n=6: c3-only fails (3 collisions), (c3,c5) determines H  [EXACT, exhaustive]")
print(f"  n=7: c3-only {'FAILS' if n_collisions > 0 else 'may work'} ({n_collisions} collisions in {N_samples} samples)")
if n_collisions_c5 > 0:
    print(f"       (c3,c5) {'FAILS' if n_collisions_c5 > 0 else 'works'} ({n_collisions_c5} collisions in {N_samples_c5} samples)")
    if 'n_coll_full' in dir() and n_coll_full == 0:
        print(f"       (c3,c5,c7) resolves all collisions")
    elif 'n_coll_full' in dir():
        print(f"       (c3,c5,c7): {n_coll_full} collisions remain")
else:
    print(f"       (c3,c5) NO collisions in {N_samples_c5} samples")
print(f"  n=8: c3-only: {n_collisions_8} collisions in {N_samples_8} samples")
print()

print("""
  PATTERN:
  At each n, c3-only SRCP fails when n is large enough that
  5-cycles become structurally important (c5 varies within c3-classes).

  The MINIMAL SRCP depth needed:
    n=3,4: c3-only suffices (no 5-cycles possible)
    n=5: c3-only suffices (c5 is constant within each c3-class)
    n=6: need (c3,c5) — 5-cycles create new structure
    n=7: need (c3,c5) or (c3,c5,c7)? — depends on whether 7-cycles matter
    n=8+: likely need (c3,c5,c7) — all odd cycles up to n contribute

  CONJECTURE (refined): SRCP with (c3, c5, ..., c_{n-2}) per arc
  determines H for tournaments on n vertices. The cycle lengths
  needed grow with n but are always bounded by the tournament order.
""")

print("opus-2026-03-21-S126")
