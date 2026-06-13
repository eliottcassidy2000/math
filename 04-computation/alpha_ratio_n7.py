#!/usr/bin/env python3
"""
alpha_ratio_n7.py — opus-2026-03-14-S75

Verify α₁ ≥ α₂ at n=7 (2^21 = 2M tournaments).
Also explore the max ratio α₂/α₁ and the I(-1) distribution.

Key question: Does α₂/α₁ ≤ 1/2 always hold? (Observed at n=6)
If so, this is MUCH stronger than α₁ ≥ α₂.
"""

from itertools import combinations
import sys

def count_ham_paths_dp(adj, n):
    """Count Hamiltonian paths via bitmask DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp:
                continue
            cnt = dp[key]
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    new_key = (mask | (1 << u), u)
                    dp[new_key] = dp.get(new_key, 0) + cnt
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def find_3cycles(adj, n):
    """Find all directed 3-cycles."""
    cycles = set()
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.add(frozenset([i, j, k]))
                elif adj[j][i] and adj[i][k] and adj[k][j]:
                    cycles.add(frozenset([i, j, k]))
    return list(cycles)

def find_5cycles(adj, n):
    """Find all directed 5-cycles (chordless)."""
    cycles = set()
    for combo in combinations(range(n), 5):
        verts = list(combo)
        # Check all possible cyclic orderings (5!/5 = 24, but we check half = 12)
        from itertools import permutations as perms
        found = False
        for perm in perms(verts):
            if found:
                break
            is_cycle = True
            for idx in range(5):
                if not adj[perm[idx]][perm[(idx+1) % 5]]:
                    is_cycle = False
                    break
            if not is_cycle:
                continue
            # Check chordless: no chord should exist
            chordless = True
            for idx in range(5):
                v1 = perm[idx]
                v2 = perm[(idx+2) % 5]
                if adj[v1][v2]:
                    chordless = False
                    break
            if chordless:
                cycles.add(frozenset(combo))
                found = True
    return list(cycles)

def find_all_odd_cycles(adj, n):
    """Find 3-cycles and 5-cycles. At n=7, 7-cycles are also possible but rare
    and irrelevant for conflict graph (full Hamilton cycles)."""
    c3 = find_3cycles(adj, n)
    c5 = find_5cycles(adj, n)
    # 7-cycles at n=7 would be Hamiltonian cycles, which are vertex-disjoint from nothing
    # So they contribute to α₁ but not α₂ (can't be disjoint from any other cycle)
    # Actually they use ALL 7 vertices so they're disjoint from no other cycle.
    # Include them for completeness.
    c7 = []
    if n == 7:
        # Check for directed 7-cycles (Hamiltonian cycles)
        verts = list(range(7))
        from itertools import permutations as perms
        seen = set()
        for perm in perms(verts):
            # Normalize: start from 0
            if perm[0] != 0:
                continue
            is_cycle = True
            for idx in range(7):
                if not adj[perm[idx]][perm[(idx+1) % 7]]:
                    is_cycle = False
                    break
            if is_cycle:
                # Canonical form: direction with perm[1] < perm[6]
                if perm[1] < perm[6]:
                    c7.append(frozenset(verts))
                    break  # All 7-cycles use the same vertex set
    return c3 + c5 + c7

def compute_alpha12(adj, n):
    """Compute α₁ and α₂ of independence polynomial of CG(T)."""
    cycles = find_all_odd_cycles(adj, n)
    alpha1 = len(cycles)
    # Count independent pairs (vertex-disjoint)
    alpha2 = 0
    for i in range(len(cycles)):
        for j in range(i + 1, len(cycles)):
            if len(cycles[i] & cycles[j]) == 0:
                alpha2 += 1
    return alpha1, alpha2

# n=7: 2^21 = 2097152 tournaments
n = 7
edges = [(i, j) for i in range(n) for j in range(i+1, n)]
num_edges = len(edges)
total = 2 ** num_edges

print(f"n={n}: {total} tournaments, {num_edges} edges")
print()

max_ratio = 0
max_ratio_info = None
violations = 0
alpha_dist = {}
im1_dist = {}
h_alpha_pairs = {}

# Process in chunks and report progress
chunk = total // 20
for bits in range(total):
    if bits % chunk == 0 and bits > 0:
        print(f"  Progress: {bits}/{total} ({100*bits//total}%), violations so far: {violations}, max α₂/α₁: {max_ratio:.4f}")

    adj = [[False] * n for _ in range(n)]
    for idx, (i, j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    a1, a2 = compute_alpha12(adj, n)

    if a2 > a1:
        violations += 1
        if violations <= 5:
            h = count_ham_paths_dp(adj, n)
            print(f"  VIOLATION! bits={bits}, α₁={a1}, α₂={a2}, H={h}")

    if a1 > 0:
        ratio = a2 / a1
        if ratio > max_ratio:
            max_ratio = ratio
            max_ratio_info = (bits, a1, a2)

    key = (a1, a2)
    alpha_dist[key] = alpha_dist.get(key, 0) + 1

    im1 = 1 - a1 + a2
    im1_dist[im1] = im1_dist.get(im1, 0) + 1

print()
print(f"RESULTS for n={n}:")
print(f"  Total tournaments: {total}")
print(f"  α₂ > α₁ violations: {violations}")
print(f"  Max ratio α₂/α₁: {max_ratio:.6f}")
if max_ratio_info:
    b, a1, a2 = max_ratio_info
    print(f"    Achieved at bits={b}, α₁={a1}, α₂={a2}")

print()
print("  (α₁, α₂) distribution:")
for key in sorted(alpha_dist.keys()):
    a1, a2 = key
    count = alpha_dist[key]
    ratio = a2/a1 if a1 > 0 else 0
    im1 = 1 - a1 + a2
    print(f"    α₁={a1:3d}, α₂={a2:3d}: {count:8d} tours, ratio={ratio:.4f}, I(-1)={im1}")

print()
print("  I(-1) distribution:")
for im1 in sorted(im1_dist.keys()):
    print(f"    I(-1) = {im1:4d}: {im1_dist[im1]:8d} tournaments")

print()
print("  Max I(-1):", max(im1_dist.keys()))
print("  Min I(-1):", min(im1_dist.keys()))
print(f"  I(-1) ≤ 1 for all? {max(im1_dist.keys()) <= 1}")
