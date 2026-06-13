#!/usr/bin/env python3
"""
Cycle Count Gaps — Which total odd-cycle counts are impossible?
opus-2026-03-14-S71h

KEY DISCOVERY: At n=5, Ω is always complete, so H = 1 + 2m where m = #cycles.
H=7 is absent because m=3 is impossible (t₃=3 forces c₅≥1, giving m≥4).

Questions:
1. What cycle counts are impossible at each n?
2. Does the gap at m=3 persist at all n? (If so, H=7 is permanently forbidden.)
3. Is there a gap at m=10 that explains H=21?
4. When Ω is NOT complete (n≥6), what's the actual I(Ω,2) = H relationship?
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict

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
print("CYCLE COUNT GAPS AND IMPOSSIBLE H VALUES")
print("=" * 70)
print()

# n=3: only 3-cycles possible (length 3 on 3 vertices)
for n in range(3, 6):
    num_edges = n*(n-1)//2
    total = 1 << num_edges

    cycle_counts = Counter()
    cycle_detail = Counter()  # (t3, c5, c7, ...) -> count

    for bits in range(total):
        adj = tournament_from_bits(n, bits)
        cycles = find_all_directed_odd_cycles(n, adj, max_len=n)
        m = len(cycles)
        cycle_counts[m] += 1

        # Breakdown by cycle length
        by_len = Counter()
        for c, vs in cycles:
            by_len[len(c)] += 1
        key = tuple(by_len.get(k, 0) for k in range(3, n+1, 2))
        cycle_detail[key] += 1

    max_m = max(cycle_counts.keys())
    missing = [m for m in range(max_m+1) if m not in cycle_counts]

    print(f"n={n}:")
    print(f"  Total tournaments: {total}")
    print(f"  Cycle count range: 0 to {max_m}")
    print(f"  Missing counts: {missing}")
    print(f"  Distribution: {dict(sorted(cycle_counts.items()))}")
    print(f"  Corresponding absent H (if Ω complete): {[1+2*m for m in missing]}")
    print()

    print(f"  Cycle length breakdown (t₃, c₅, ...):")
    for key in sorted(cycle_detail.keys()):
        labels = [f"c{2*i+3}={key[i]}" for i in range(len(key)) if key[i] > 0]
        total_cycles = sum(key)
        print(f"    ({', '.join(labels) if labels else 'none'}), total={total_cycles}: {cycle_detail[key]}")
    print()

# n=6: exhaustive but only 3-cycles and 5-cycles (skip 7-cycles — too slow)
print("n=6: (3-cycles and 5-cycles only — ignoring potential 7-cycles... wait, n=6 has no 7-cycles)")
print("     Actually max odd cycle length at n=6 is 5 (6 vertices, no 7-cycle possible)")
n = 6
num_edges = 15
total = 1 << num_edges

cycle_counts = Counter()
t3_c5_dist = Counter()

for bits in range(total):
    adj = tournament_from_bits(n, bits)
    cycles = find_all_directed_odd_cycles(n, adj, max_len=5)
    m = len(cycles)
    cycle_counts[m] += 1

    t3 = sum(1 for c, _ in cycles if len(c) == 3)
    c5 = sum(1 for c, _ in cycles if len(c) == 5)
    t3_c5_dist[(t3, c5)] += 1

max_m = max(cycle_counts.keys())
missing = [m for m in range(max_m+1) if m not in cycle_counts]

print(f"  Total tournaments: {total}")
print(f"  Cycle count range: 0 to {max_m}")
print(f"  Missing counts: {missing}")
print(f"  Corresponding absent H (if Ω complete): {[1+2*m for m in missing]}")
print(f"  Distribution: {dict(sorted(cycle_counts.items()))}")
print()

print(f"  (t₃, c₅) distribution:")
for (t3, c5) in sorted(t3_c5_dist.keys()):
    total_cycles = t3 + c5
    missing_flag = " ← total=3!" if total_cycles == 3 else ""
    missing_flag += " ← total=10!" if total_cycles == 10 else ""
    print(f"    t₃={t3:2d}, c₅={c5:2d}, total={total_cycles:3d}: {t3_c5_dist[(t3,c5)]:5d}{missing_flag}")

print()
print("=" * 70)
print("WHY IS TOTAL=3 IMPOSSIBLE?")
print("=" * 70)
print()

# At n=5: t₃=3 forces c₅=1.
# At n=6: does t₃=3 still force c₅≥1?

# For total=3, we need:
# (a) t₃=3, c₅=0: 3 triangles, no 5-cycle
# (b) t₃=2, c₅=1: impossible because c₅ at n≤5 needs t₃≥3
#     Actually at n=6, a 5-cycle uses 5 vertices, leaving 1 unused
# (c) t₃=1, c₅=2: impossible similarly
# (d) t₃=0, c₅=3: impossible (no triangles means very transitive)

# The key structural constraint: at n=5, 3 triangles force 1 pentagonic cycle.
# Does this persist at n=6?

print("Checking: does t₃=3 always force c₅≥1 at n=6?")
t3_3 = [(t3, c5, cnt) for (t3, c5), cnt in t3_c5_dist.items() if t3 == 3]
print(f"  Cases with t₃=3 at n=6:")
for t3, c5, cnt in sorted(t3_3):
    print(f"    t₃=3, c₅={c5}: {cnt} tournaments")

if all(c5 >= 1 for _, c5, _ in t3_3):
    print("  YES: t₃=3 always has c₅≥1 at n=6!")
else:
    t3_3_c5_0 = [cnt for _, c5, cnt in t3_3 if c5 == 0]
    print(f"  NO: t₃=3, c₅=0 exists at n=6: {sum(t3_3_c5_0)} tournaments")

print()

# Similarly check: what (t₃, c₅) gives total=10?
print("Cases near total=10 at n=6:")
for (t3, c5), cnt in sorted(t3_c5_dist.items()):
    if 8 <= t3+c5 <= 12:
        print(f"  t₃={t3:2d}, c₅={c5:2d}, total={t3+c5}: {cnt}")

print()
print("=" * 70)
print("THE STRUCTURAL EXPLANATION")
print("=" * 70)
print()
print("At n=5: tournament on 5 vertices has C(5,3)=10 triples.")
print("Each triple is either a 3-cycle (t₃) or transitive (10-t₃).")
print("t₃ ranges from 0 to 4 (max for n=5).")
print()
print("The 5-cycle structure:")
print("  A directed Hamiltonian cycle on 5 vertices uses all 5 vertices.")
print("  The number of directed 5-cycles in a 5-tournament depends on")
print("  the tournament structure, not just t₃ — but t₃≥3 forces c₅≥1.")
print()
print("WHY t₃=3 forces c₅≥1 at n=5:")
print("  3 directed 3-cycles on C(5,3)=10 triples means 7 transitive triples.")
print("  The 3 cyclic triples create a pattern that, combined with the arc")
print("  directions of the transitive triples, must form at least one")
print("  directed Hamiltonian cycle (5-cycle).")
print()
print("This is related to the 'packing' constraint: 3 triangles in a")
print("tournament on 5 vertices create enough cyclic structure to force")
print("a Hamiltonian cycle.")
