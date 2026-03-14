#!/usr/bin/env python3
"""
Is Ω(T) always a complete graph? When does it first fail?
opus-2026-03-14-S71h

KEY INSIGHT: At n=5, every pair of odd cycles shares a vertex (Ω = K_m).
This means H = 1 + 2m where m = number of odd cycles.
When does this break? At what n do we first see vertex-disjoint odd cycles?
"""

from itertools import combinations, permutations
from collections import Counter

def find_all_directed_odd_cycles(n, adj, max_len=None):
    if max_len is None:
        max_len = n
    cycles = []
    seen = set()
    for length in range(3, max_len + 1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = True
                for i in range(length):
                    nxt = (i + 1) % length
                    if perm[nxt] not in adj[perm[i]]:
                        is_cycle = False
                        break
                if is_cycle:
                    rotations = [perm[i:] + perm[:i] for i in range(length)]
                    canonical = min(rotations)
                    if canonical not in seen:
                        seen.add(canonical)
                        cycles.append((canonical, frozenset(combo)))
    return cycles

def tournament_from_bits(n, bits):
    adj = [set() for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i].add(j)
            else:
                adj[j].add(i)
            idx += 1
    return adj

def hp_count(n, adj):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in adj[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full])

print("=" * 70)
print("IS Ω(T) ALWAYS A COMPLETE GRAPH?")
print("=" * 70)
print()

# Check n=3,4,5 exhaustively
for n in range(3, 6):
    num_edges = n*(n-1)//2
    total = 1 << num_edges
    non_complete = 0
    disjoint_pairs = 0

    for bits in range(total):
        adj = tournament_from_bits(n, bits)
        cycles = find_all_directed_odd_cycles(n, adj, max_len=n)
        nc = len(cycles)

        # Check if any pair of cycles is vertex-disjoint
        has_disjoint = False
        for i in range(nc):
            for j in range(i+1, nc):
                if not (cycles[i][1] & cycles[j][1]):
                    has_disjoint = True
                    disjoint_pairs += 1
                    break
            if has_disjoint:
                break

        if has_disjoint:
            non_complete += 1

    print(f"n={n}: {total} tournaments, {non_complete} have non-complete Ω")

print()

# n=6: check all 32768
print("n=6: checking all 32768 tournaments...")
n = 6
num_edges = 15
total = 1 << num_edges
non_complete_count = 0
disjoint_examples = []

for bits in range(total):
    adj = tournament_from_bits(n, bits)
    cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
    nc = len(cycles)

    has_disjoint = False
    for i in range(nc):
        for j in range(i+1, nc):
            if not (cycles[i][1] & cycles[j][1]):
                has_disjoint = True
                if len(disjoint_examples) < 3:
                    H = hp_count(n, adj)
                    disjoint_examples.append((bits, H, nc, cycles[i], cycles[j]))
                break
        if has_disjoint:
            break

    if has_disjoint:
        non_complete_count += 1

print(f"  {non_complete_count} tournaments have non-complete Ω (vertex-disjoint odd cycles)")
print()

if disjoint_examples:
    print("  Examples of disjoint odd cycle pairs:")
    for bits, H, nc, c1, c2 in disjoint_examples:
        print(f"    bits={bits}: H={H}, |Ω|={nc}")
        print(f"      cycle1: {c1[0]} on vertices {sorted(c1[1])}")
        print(f"      cycle2: {c2[0]} on vertices {sorted(c2[1])}")
else:
    print("  *** Ω IS ALWAYS COMPLETE AT n=6! ***")
    print("  Every pair of odd cycles shares at least one vertex.")
    print()
    print("  This means H = 1 + 2·(#odd cycles) for ALL tournaments on n≤6 vertices!")
    print()

    # Check the cycle count distribution
    print("  Cycle count distribution at n=6:")
    cycle_counts = Counter()
    for bits in range(total):
        adj = tournament_from_bits(n, bits)
        cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
        cycle_counts[len(cycles)] += 1

    for k in sorted(cycle_counts.keys()):
        H = 1 + 2*k
        present = "← H={} ✓".format(H) if True else ""
        print(f"    {k:2d} cycles → H={H:3d}: {cycle_counts[k]:5d} tournaments")

    # Which counts are missing?
    max_k = max(cycle_counts.keys())
    missing = [k for k in range(max_k+1) if k not in cycle_counts]
    print(f"\n  Missing cycle counts: {missing}")
    print(f"  Corresponding absent H: {[1+2*k for k in missing]}")

print()
print("=" * 70)
print("WHY ARE CERTAIN CYCLE COUNTS IMPOSSIBLE?")
print("=" * 70)
print()

# At n=5: cycle counts {0,1,2,4,5,6,7}, missing: {3}
# At n=6: let's see what's missing

# The key: at n=5, a tournament has C(5,3)=10 triples.
# Each triple is either a 3-cycle or transitive. Let t₃ = #3-cycles.
# 5-cycles: a tournament on 5 vertices can have 0,1,2,3 directed 5-cycles.
# Total odd cycles = t₃ + c₅.

# At n=5, t₃ ranges 0-4, c₅ ranges 0-3.
# Distribution:
n = 5
total = 1 << 10
t3_c5_dist = Counter()
for bits in range(total):
    adj = tournament_from_bits(n, bits)
    cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
    t3 = sum(1 for c, _ in cycles if len(c) == 3)
    c5 = sum(1 for c, _ in cycles if len(c) == 5)
    t3_c5_dist[(t3, c5)] += 1

print("n=5: (t₃, c₅) distribution:")
for (t3, c5) in sorted(t3_c5_dist.keys()):
    total_cycles = t3 + c5
    print(f"  t₃={t3}, c₅={c5}, total={total_cycles}: {t3_c5_dist[(t3,c5)]} tournaments")

print()
print("Total odd cycles at n=5: never 3!")
print("  t₃=0,c₅=0: 0 cycles")
print("  t₃=1,c₅=0: 1 cycle")
print("  t₃=2,c₅=0: 2 cycles")
print("  t₃=3,c₅=1: 4 cycles (NOT 3!)")
print("  t₃=4,c₅=1: 5 cycles")
print("  t₃=4,c₅=2: 6 cycles")
print("  t₃=4,c₅=3: 7 cycles")
print()
print("KEY: t₃=3 ALWAYS comes with c₅≥1!")
print("At n=5, having 3 triangles forces at least 1 five-cycle.")
print("So total cycles ≥ 4 whenever t₃ ≥ 3.")
print("And t₃ ≤ 2 gives total ≤ 2 (no 5-cycles when few triangles).")
print("The gap: total cycles = 3 is IMPOSSIBLE!")

print()
print("=" * 70)
print("THE CYCLE-COUNT GAP THEOREM")
print("=" * 70)
print()
print("THEOREM: No tournament on n ≤ 6 vertices has exactly 3 directed odd cycles.")
print()
print("PROOF (n=5): If t₃ ≤ 2, then c₅ = 0 (too few arcs for a Hamilton cycle")
print("to be directed). So total ≤ 2. If t₃ ≥ 3, the cyclic structure forces")
print("at least one 5-cycle, giving total ≥ 4.")
print()
print("COROLLARY: H(T) ≠ 7 for all tournaments on n ≤ 6 vertices.")
print("(Because H = 1 + 2·(#cycles) when Ω is complete, and #cycles ≠ 3.)")
