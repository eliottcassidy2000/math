#!/usr/bin/env python3
"""
n7_alpha2_sample_87b.py — opus-2026-03-14-S87b

Sample n=7 tournaments to explore α₂ distribution.
Full enumeration of 2^21 = 2097152 tournaments is too slow for full analysis,
but we can sample efficiently.

Key questions:
1. What α₂ values are achievable at n=7?
2. Is α₂=3 impossible (like at n=6)?
3. What is the maximum α₂?
4. Does the gap structure persist?
"""

import random
from itertools import combinations
from collections import Counter
import sys

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

def find_3cycles(adj, n):
    """Find all directed 3-cycles. Returns list of (tuple, frozenset) pairs."""
    cycles = []
    for a, b, c in combinations(range(n), 3):
        # Check both orientations
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles.append(((a,b,c), frozenset({a,b,c})))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            cycles.append(((a,c,b), frozenset({a,b,c})))
    return cycles

def compute_alpha2(cycles):
    """Count vertex-disjoint pairs of 3-cycles (α₂)."""
    nc = len(cycles)
    if nc < 2:
        return 0
    count = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if not (cycles[i][1] & cycles[j][1]):  # disjoint vertex sets
                count += 1
    return count

def compute_H_dp(adj, n):
    """Hamiltonian path count via Held-Karp DP."""
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

# ══════════════════════════════════════════════════════════════════
# SAMPLING
# ══════════════════════════════════════════════════════════════════

n = 7
N_SAMPLES = 100000

print("=" * 70)
print(f"n=7 TOURNAMENT SAMPLING — {N_SAMPLES:,} random tournaments")
print("=" * 70)

random.seed(42)

alpha2_counter = Counter()
alpha1_counter = Counter()
alpha2_to_H = {}  # α₂ → set of H values seen
max_alpha2 = 0
max_alpha2_tournament = None

for trial in range(N_SAMPLES):
    adj = random_tournament(n)
    cycles = find_3cycles(adj, n)
    nc = len(cycles)  # α₁
    a2 = compute_alpha2(cycles)

    alpha1_counter[nc] += 1
    alpha2_counter[a2] += 1

    if a2 not in alpha2_to_H:
        alpha2_to_H[a2] = set()

    # Only compute H for interesting cases
    if a2 > 3 or (a2 > 0 and len(alpha2_to_H[a2]) < 5):
        H = compute_H_dp(adj, n)
        alpha2_to_H[a2].add(H)

    if a2 > max_alpha2:
        max_alpha2 = a2
        max_alpha2_tournament = [row[:] for row in adj]
        H = compute_H_dp(adj, n)
        print(f"  New max α₂ = {a2} at trial {trial}, H = {H}")

    if (trial + 1) % 20000 == 0:
        print(f"  ... {trial+1}/{N_SAMPLES}")
        sys.stdout.flush()

print(f"\nDone sampling {N_SAMPLES:,} tournaments.")

# ── Results ──────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("RESULT 1: α₂ DISTRIBUTION")
print("=" * 70)
for a2 in sorted(alpha2_counter.keys()):
    pct = 100 * alpha2_counter[a2] / N_SAMPLES
    bar = "#" * max(1, int(pct / 2))
    print(f"  α₂ = {a2:2d}: {alpha2_counter[a2]:6d} ({pct:5.1f}%) {bar}")

print(f"\n  Max α₂ seen: {max_alpha2}")

# Check for gaps
a2_vals = sorted(alpha2_counter.keys())
max_a2 = max(a2_vals)
gaps = [v for v in range(max_a2 + 1) if v not in alpha2_counter]
if gaps:
    print(f"  GAPS in α₂: {gaps}")
else:
    print(f"  No gaps in α₂ range [0, {max_a2}]")

print("\n" + "=" * 70)
print("RESULT 2: α₁ DISTRIBUTION")
print("=" * 70)
for a1 in sorted(alpha1_counter.keys()):
    pct = 100 * alpha1_counter[a1] / N_SAMPLES
    print(f"  α₁ = {a1:2d}: {alpha1_counter[a1]:6d} ({pct:5.1f}%)")

print("\n" + "=" * 70)
print("RESULT 3: H VALUES BY α₂")
print("=" * 70)
for a2 in sorted(alpha2_to_H.keys()):
    H_vals = sorted(alpha2_to_H[a2])
    print(f"  α₂ = {a2}: H ∈ {{{', '.join(map(str, H_vals[:10]))}}}" +
          (f" ... ({len(H_vals)} values)" if len(H_vals) > 10 else ""))

# ── Disjoint 3-cycle pair analysis ──────────────────────────────

print("\n" + "=" * 70)
print("RESULT 4: DISJOINT 3-CYCLE PAIR ANALYSIS")
print("=" * 70)

# How many disjoint 3-subset pairs exist at n=7?
pairs_7 = []
for A in combinations(range(7), 3):
    remaining = [x for x in range(7) if x not in A]
    for B in combinations(remaining, 3):
        pair = frozenset([frozenset(A), frozenset(B)])
        if pair not in pairs_7:
            pairs_7.append(pair)

print(f"  Disjoint 3-subset pairs: {len(pairs_7)}")
print(f"  Each pair leaves 1 vertex uncovered")

# For the max-α₂ tournament, which pairs are both-cyclic?
if max_alpha2_tournament:
    adj = max_alpha2_tournament
    cycles = find_3cycles(adj, n)
    cycle_sets = [c[1] for c in cycles]
    print(f"\n  Max-α₂ tournament: {len(cycles)} 3-cycles")

    both_cyclic = 0
    for pair in pairs_7:
        A, B = list(pair)
        if A in cycle_sets and B in cycle_sets:
            both_cyclic += 1
    print(f"  Both-cyclic pairs: {both_cyclic} (= α₂ = {max_alpha2})")

    # Which vertex is always uncovered?
    uncovered_count = Counter()
    for pair in pairs_7:
        A, B = list(pair)
        if A in cycle_sets and B in cycle_sets:
            uncovered = [v for v in range(7) if v not in A and v not in B]
            uncovered_count[uncovered[0]] += 1
    print(f"  Uncovered vertex distribution: {dict(uncovered_count)}")

# ── Compare with n=6 ──────────────────────────────────────────

print("\n" + "=" * 70)
print("COMPARISON: n=6 vs n=7")
print("=" * 70)

print(f"  n=6: α₂ ∈ {{0,1,2,4}}, gap at 3")
print(f"  n=7: α₂ ∈ {sorted(alpha2_counter.keys())}")
print(f"  n=6: 10 complementary partitions, max α₂ = 4")
print(f"  n=7: 70 disjoint 3-subset pairs, max α₂ = {max_alpha2}")
print(f"  n=6: complete → BIBD")
print(f"  n=7: richer structure (what design?)")
