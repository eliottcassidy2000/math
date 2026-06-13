#!/usr/bin/env python3
"""
n6_topology_fast_87.py — opus-2026-03-14-S87

Fast n=6 topology exploration using optimized cycle enumeration.
The key insight: instead of brute-force permutation enumeration,
we find directed cycles via DFS.

Key questions for n=6:
1. What fraction of tournaments have Ω non-complete (α₂ > 0)?
2. What (α₁, α₂) pairs are achievable?
3. Does H still determine α-vector at n=6?
4. The H + 2χ = 3 + 6α₂ identity at n=6.
"""

from itertools import combinations
from collections import defaultdict, Counter
import sys

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj, bits

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            if dp[S][v] == 0: continue
            for w in range(n):
                if S & (1 << w): continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def find_directed_cycles_of_length(adj, n, length):
    """Find all directed cycles of given length using DFS.
    Returns list of tuples (v0, v1, ..., v_{L-1}) normalized so min vertex is first.
    """
    cycles = set()
    for start in range(n):
        # DFS from start, looking for cycles of exact length back to start
        def dfs(path, visited):
            if len(path) == length:
                # Check if there's an arc from last vertex back to start
                if adj[path[-1]][start]:
                    # Normalize: rotate so min vertex first
                    min_idx = path.index(min(path))
                    normalized = tuple(path[min_idx:] + path[:min_idx])
                    cycles.add(normalized)
                return
            last = path[-1]
            for v in range(n):
                if v in visited and v != start:
                    continue
                if v == start and len(path) < length:
                    continue  # Can only return to start at the right length
                if not adj[last][v]:
                    continue
                if v != start:
                    path.append(v)
                    visited.add(v)
                    dfs(path, visited)
                    path.pop()
                    visited.discard(v)
        dfs([start], {start})
    return list(cycles)

def get_all_odd_cycles(adj, n):
    """Get all odd directed cycles of length 3, 5, etc."""
    all_cycles = []
    for length in range(3, n+1, 2):
        cycles = find_directed_cycles_of_length(adj, n, length)
        for c in cycles:
            all_cycles.append((c, frozenset(c)))
    return all_cycles

def compute_alpha_vector(cycles):
    """Compute α-vector (number of independent sets of each size) from conflict graph."""
    nc = len(cycles)
    if nc == 0:
        return {0: 1}

    # Build conflict adjacency as bitsets for speed
    conflict = [0] * nc
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][1] & cycles[j][1]:  # share vertex
                conflict[i] |= (1 << j)
                conflict[j] |= (1 << i)

    # Count independent sets by size
    alpha = defaultdict(int)
    for mask in range(1 << nc):
        # Check independence using bitset
        is_indep = True
        bits_remaining = mask
        while bits_remaining:
            v = (bits_remaining & -bits_remaining).bit_length() - 1
            bits_remaining &= bits_remaining - 1
            if conflict[v] & mask:
                is_indep = False
                break
        if is_indep:
            alpha[bin(mask).count('1')] += 1

    return dict(alpha)

# ══════════════════════════════════════════════════════════════════
# MAIN COMPUTATION
# ══════════════════════════════════════════════════════════════════

print("=" * 70)
print("n=6 TOPOLOGY — FAST ENUMERATION")
print("=" * 70)

n = 6
m = n * (n - 1) // 2  # 15
total = 1 << m  # 32768

print(f"n={n}, C(n,2)={m}, total tournaments = {total}")
print()

# Collect results
H_counter = Counter()
alpha_pair_counter = Counter()  # (α₁, α₂) → count
H_plus_2chi_counter = Counter()
H_to_alpha = defaultdict(list)

count = 0
for adj, bits in all_tournaments(n):
    H = compute_H_dp(adj, n)
    H_counter[H] += 1

    cycles = get_all_odd_cycles(adj, n)
    nc = len(cycles)

    if nc == 0:
        alpha = {0: 1}
        chi = 1
        alpha1 = 0
        alpha2 = 0
    else:
        alpha = compute_alpha_vector(cycles)
        alpha1 = alpha.get(1, 0)
        alpha2 = alpha.get(2, 0)
        chi = sum(alpha.get(k, 0) * ((-1)**k) for k in range(max(alpha.keys())+1))

    alpha_pair_counter[(alpha1, alpha2)] += 1
    H_plus_2chi_counter[H + 2 * chi] += 1
    H_to_alpha[H].append((alpha1, alpha2))

    count += 1
    if count % 2000 == 0:
        print(f"  ... {count}/{total} ({100*count/total:.0f}%)")
        sys.stdout.flush()

print(f"\nDone: {count} tournaments processed.")

# ── Results ──────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("RESULT 1: H DISTRIBUTION AT n=6")
print("=" * 70)
for H in sorted(H_counter.keys()):
    print(f"  H = {H:3d}: {H_counter[H]:5d} tournaments")

print("\n" + "=" * 70)
print("RESULT 2: (α₁, α₂) PAIR DISTRIBUTION")
print("=" * 70)
print(f"  {'α₁':>4} {'α₂':>4} {'count':>6} {'H=1+2α₁+4α₂':>14} {'χ=1-α₁+α₂':>12}")
for pair in sorted(alpha_pair_counter.keys()):
    a1, a2 = pair
    H = 1 + 2*a1 + 4*a2
    chi = 1 - a1 + a2
    print(f"  {a1:>4} {a2:>4} {alpha_pair_counter[pair]:>6} {H:>14} {chi:>12}")

print("\n" + "=" * 70)
print("RESULT 3: H + 2χ VALUES (should be 3 + 6α₂)")
print("=" * 70)
for v in sorted(H_plus_2chi_counter.keys()):
    alpha2 = (v - 3) / 6
    print(f"  H + 2χ = {v:3d} (implies α₂ = {alpha2:.1f}): {H_plus_2chi_counter[v]:5d} tournaments")

print("\n" + "=" * 70)
print("RESULT 4: H VALUES WITH MULTIPLE α-VECTORS")
print("=" * 70)
multi_count = 0
for H in sorted(H_to_alpha.keys()):
    distinct = set(H_to_alpha[H])
    if len(distinct) > 1:
        multi_count += 1
        print(f"  H = {H:3d}: {len(distinct)} distinct (α₁,α₂) pairs: {sorted(distinct)}")
if multi_count == 0:
    print("  NONE — H determines (α₁, α₂) uniquely!")
else:
    print(f"  Total: {multi_count} H values with non-unique α-vector")

print("\n" + "=" * 70)
print("RESULT 5: MISSING α₁ VALUES AND FORBIDDEN H")
print("=" * 70)
alpha1_vals = sorted(set(a1 for (a1, a2) in alpha_pair_counter.keys()))
max_a1 = max(alpha1_vals)
missing_a1 = sorted(set(range(max_a1+1)) - set(alpha1_vals))
print(f"  α₁ range: [0, {max_a1}]")
print(f"  α₁ values present: {alpha1_vals}")
print(f"  Missing α₁: {missing_a1}")
for a1 in missing_a1:
    print(f"    α₁ = {a1} → Would give H = {1+2*a1} (with α₂=0) — FORBIDDEN")

# Check which H values actually appear
all_H = sorted(H_counter.keys())
all_possible = set(range(1, max(all_H)+1, 2))  # all odd values
missing_H = sorted(all_possible - set(all_H))
print(f"\n  All odd H in [1, {max(all_H)}]: {len(all_possible)} values")
print(f"  H values achieved: {len(all_H)}")
print(f"  Forbidden H values at n=6: {missing_H}")

print("\n" + "=" * 70)
print("RESULT 6: THE n=6 TRANSITION — Ω COMPLETENESS")
print("=" * 70)
complete = sum(v for (a1, a2), v in alpha_pair_counter.items() if a2 == 0)
non_complete = total - complete
print(f"  Ω complete (α₂ = 0): {complete} ({100*complete/total:.1f}%)")
print(f"  Ω non-complete (α₂ > 0): {non_complete} ({100*non_complete/total:.1f}%)")
print(f"\n  This is the TOPOLOGICAL PHASE TRANSITION at n=6!")
print(f"  Below n=6: Ω always complete, H+2χ=3 exactly.")
print(f"  At n=6: {100*non_complete/total:.1f}% of tournaments cross the transition.")
