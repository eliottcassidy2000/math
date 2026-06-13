#!/usr/bin/env python3
"""
pi_forbidden6_89c.py — Why are H = 7, 21, 35, 39 forbidden at n=6?
opus-2026-03-14-S89c

At n=5: H=7 is forbidden (proved earlier: t₃+t₅=3 structurally forbidden)
At n=6: H ∈ {7, 21, 35, 39} are ALL forbidden

Questions:
1. Is there a mod structure? 7 ≡ 1 mod 2, 21 ≡ 1, 35 ≡ 1, 39 ≡ 1. All odd ✓
2. 7 = 7, 21 = 3×7, 35 = 5×7, 39 = 3×13. Is 7 special?
3. Using OCF: H = 1 + 2t₃ + 4d₂ + 8d₃ + ... where t₃ = # 3-cycles, etc.

The OCF decomposition gives constraints on possible (t₃, t₅, d_pairs, ...) tuples.
Let's find which tuples are realized and which H values they produce.
"""

from itertools import combinations
from collections import Counter

def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(1 << m):
        adj = {v: set() for v in range(n)}
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[i].add(j)
            else:
                adj[j].add(i)
        yield adj

def count_hp(adj, n):
    dp = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for u in adj[v]:
                if mask & (1 << u) == 0:
                    new_mask = mask | (1 << u)
                    dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full].values())

def count_3cycles(adj, n):
    count = 0
    for i in range(n):
        for j in adj[i]:
            for k in adj[j]:
                if i in adj[k]:
                    count += 1
    return count // 3

def count_directed_cycles(adj, n, length):
    """Count directed cycles of given length (canonical: min vertex is start)."""
    count = 0
    for start in range(n):
        stack = [(start, [start], 1 << start)]
        while stack:
            v, path, mask = stack.pop()
            if len(path) == length:
                if start in adj[v]:
                    count += 1
                continue
            for u in adj[v]:
                if u < start:
                    continue
                if mask & (1 << u):
                    continue
                stack.append((u, path + [u], mask | (1 << u)))
    return count

print("=" * 70)
print("PART 1: Complete (H, t₃, t₅) data for n=6")
print("=" * 70)

n = 6
data = []
for T in all_tournaments(n):
    H = count_hp(T, n)
    t3 = count_3cycles(T, n)
    t5 = count_directed_cycles(T, n, 5)
    data.append((H, t3, t5))

# Group by H
by_H = {}
for H, t3, t5 in data:
    if H not in by_H:
        by_H[H] = []
    by_H[H].append((t3, t5))

print(f"\n  {'H':>4} | {'count':>5} | t₃ values | t₅ values")
print(f"  {'-'*4}-+-{'-'*5}-+-{'-'*20}-+-{'-'*20}")
for H in sorted(by_H.keys()):
    entries = by_H[H]
    t3_vals = sorted(set(t for t, _ in entries))
    t5_vals = sorted(set(t for _, t in entries))
    t3_counts = Counter(t for t, _ in entries)
    print(f"  {H:4d} | {len(entries):5d} | {t3_vals} | {t5_vals}")

print()
print("=" * 70)
print("PART 2: Which (t₃, t₅) pairs exist at n=6?")
print("=" * 70)

pair_counts = Counter((t3, t5) for _, t3, t5 in data)
print(f"\n  {len(pair_counts)} distinct (t₃, t₅) pairs")

# Map from (t3, t5) to possible H values
pair_to_H = {}
for H, t3, t5 in data:
    key = (t3, t5)
    if key not in pair_to_H:
        pair_to_H[key] = set()
    pair_to_H[key].add(H)

print(f"\n  {'t₃':>3} {'t₅':>3} | {'count':>5} | H values")
print(f"  {'-'*3} {'-'*3}-+-{'-'*5}-+-{'-'*30}")
for (t3, t5) in sorted(pair_to_H.keys()):
    H_set = sorted(pair_to_H[(t3, t5)])
    count = pair_counts[(t3, t5)]
    print(f"  {t3:3d} {t5:3d} | {count:5d} | {H_set}")

print()
print("=" * 70)
print("PART 3: WHY is H=7 forbidden?")
print("=" * 70)

# H = 1 + 2t₃ + 4·(disjoint pairs) + 8·(disjoint triples) + ...
# For H=7: 7 = 1 + 2t₃ + 4c₂ + 8c₃ + ...
# Since all terms after 1 are positive:
# 2t₃ + 4c₂ + 8c₃ + ... = 6
# Possible: t₃=3, c₂=0 OR t₃=1, c₂=1 OR t₃=0, c₂=1, hmm...
# Actually: 2t₃ ≤ 6 means t₃ ≤ 3.
# But wait — the OCF formula has ALL odd cycles, not just 3-cycles.
# H = IP(G(T), 2) = 1 + 2(#cycles) + 4(#disjoint pairs) + 8(#disjoint triples)
# where "cycles" means ALL directed odd cycles.

# For n=6: total odd cycles = t₃ + t₅
# (no Hamiltonian odd cycles since n=6 is even)

# So: H = 1 + 2(t₃+t₅) + 4c₂ + 8c₃
# For H=7: 2(t₃+t₅) + 4c₂ + 8c₃ = 6
# Since 4c₂ ≥ 0 and 8c₃ ≥ 0:
# 2(t₃+t₅) ≤ 6, so t₃+t₅ ≤ 3
# Options:
# (a) t₃+t₅=3, c₂=0, c₃=0 → H = 1+6 = 7
# (b) t₃+t₅=1, c₂=1, c₃=0 → H = 1+2+4 = 7
# (c) t₃+t₅=0, c₂=0, c₃=... no, 4c₂+8c₃=6 has no solution with c₂,c₃≥0

# So H=7 requires EITHER:
# (a) exactly 3 odd cycles, all pairwise sharing a vertex
# (b) exactly 1 odd cycle and 1 disjoint pair... wait, that's impossible.
#     1 cycle means c₂=0 (can't have a pair with just 1 cycle).
#     Actually c₂ counts disjoint PAIRS of cycles. With 1 cycle, c₂=0.
#     So (b) is impossible.
# Only option: t₃+t₅=3, c₂=0.

# From the data above: which (t₃, t₅) pairs have t₃+t₅=3?
print("\n  For H=7, need t₃+t₅=3 and NO disjoint pairs among the 3 cycles.")
print(f"\n  (t₃, t₅) pairs with t₃+t₅=3:")
found = False
for (t3, t5) in sorted(pair_to_H.keys()):
    if t3 + t5 == 3:
        print(f"    t₃={t3}, t₅={t5}: exists, H values = {sorted(pair_to_H[(t3,t5)])}")
        found = True
if not found:
    print("    NONE EXIST!")

# Let's also check t₃+t₅=3 directly
total_cycle_counts = Counter(t3 + t5 for _, t3, t5 in data)
print(f"\n  Distribution of t₃+t₅:")
for total in sorted(total_cycle_counts.keys()):
    print(f"    t₃+t₅={total}: {total_cycle_counts[total]} tournaments")

print()
print("=" * 70)
print("PART 4: The H = 1 + 2Σcycles + 4Σpairs + ... decomposition")
print("=" * 70)

# For each tournament, compute the FULL IP decomposition
# IP_0 = 1, IP_1 = total odd cycles, IP_2 = disjoint pairs, etc.

def full_ip_decomposition(adj, n):
    """Compute IP coefficients for the odd-cycle disjointness graph."""
    # Find all odd cycles
    cycles = []
    for length in range(3, n+1, 2):
        for start in range(n):
            stack = [(start, [start], 1 << start)]
            while stack:
                v, path, mask = stack.pop()
                if len(path) == length:
                    if start in adj[v]:
                        cycles.append(frozenset(path))
                    continue
                for u in adj[v]:
                    if u < start:
                        continue
                    if mask & (1 << u):
                        continue
                    stack.append((u, path + [u], mask | (1 << u)))

    nc = len(cycles)
    if nc == 0:
        return [1]

    # Build disjointness graph
    g_adj = {i: set() for i in range(nc)}
    for i in range(nc):
        for j in range(i+1, nc):
            if not cycles[i] & cycles[j]:
                g_adj[i].add(j)
                g_adj[j].add(i)

    # Max packing ≤ floor(n/3) = 2 for n=6
    # So IP has degree ≤ 2 (possibly 3 if 3+3+... fits, but 3+3+3=9>6)
    # Actually 3+3=6 fits for n=6. So max packing = 2.
    # So we only need IP_0, IP_1, IP_2.

    c1 = nc
    c2 = sum(1 for i in range(nc) for j in g_adj[i] if j > i)

    return [1, c1, c2]

# Sample: compute for a few tournaments with each H value
print(f"\n  IP decomposition for representative tournaments:")

seen_H = set()
for T in all_tournaments(n):
    H = count_hp(T, n)
    if H in seen_H:
        continue
    seen_H.add(H)

    ip = full_ip_decomposition(T, n)
    # Verify
    ip_at_2 = sum(c * 2**i for i, c in enumerate(ip))

    t3 = count_3cycles(T, n)
    t5 = count_directed_cycles(T, n, 5)

    print(f"    H={H:3d}: IP = {ip}, IP(2) = {ip_at_2}, t₃={t3}, t₅={t5}, "
          f"t₃+t₅={t3+t5}, {'✓' if ip_at_2 == H else '✗'}")

print()
print("=" * 70)
print("PART 5: The forbidden values and parity constraints")
print("=" * 70)

# Missing: 7, 21, 35, 39
# H = 1 + 2c₁ + 4c₂ where c₁ = t₃+t₅, c₂ = disjoint pairs
# H is always odd ✓ (since 1 + 2c₁ + 4c₂ is odd)
# H ≡ 1 mod 2 always.
# H mod 4: 1 + 2c₁ mod 4 = 1 if c₁ even, 3 if c₁ odd.
# H mod 8: depends on c₁ mod 4 and c₂ mod 2.

# Check H mod 4 for all values
print("\n  H mod 4 analysis:")
for H in sorted(set(h for h, _, _ in data)):
    count = sum(1 for h, _, _ in data if h == H)
    print(f"    H={H:3d} ≡ {H%4} mod 4, ≡ {H%8} mod 8, count={count}")

# Missing: 7 ≡ 3 mod 4, 21 ≡ 1 mod 4, 35 ≡ 3 mod 4, 39 ≡ 3 mod 4
# Some mod 4 = 1, some = 3. No clean mod pattern.

# Try grouping:
# 7 = 1 + 2×3, needs c₁=3, c₂=0
# 21 = 1 + 2×c₁ + 4×c₂ → 20 = 2c₁ + 4c₂ → 10 = c₁ + 2c₂
# So c₁=10, c₂=0 or c₁=8, c₂=1 or c₁=6, c₂=2 etc.
# 35 = 1 + 34 → 17 = c₁ + 2c₂
# 39 = 1 + 38 → 19 = c₁ + 2c₂

print(f"\n  For missing H values, required (c₁, c₂) with c₁+2c₂ = (H-1)/2:")
for H_miss in [7, 21, 35, 39]:
    target = (H_miss - 1) // 2
    print(f"    H={H_miss}: c₁+2c₂ = {target}")
    # What (c₁, c₂) values actually occur?
    actual = set()
    for T in all_tournaments(n):
        ip = full_ip_decomposition(T, n)
        c1 = ip[1] if len(ip) > 1 else 0
        c2 = ip[2] if len(ip) > 2 else 0
        if c1 + 2*c2 == target:
            actual.add((c1, c2))
    if actual:
        print(f"      Realized (c₁,c₂): {sorted(actual)}")
        # But H for these should be H_miss. Check:
        for T in all_tournaments(n):
            ip = full_ip_decomposition(T, n)
            c1 = ip[1] if len(ip) > 1 else 0
            c2 = ip[2] if len(ip) > 2 else 0
            if c1 + 2*c2 == target:
                H_actual = count_hp(T, n)
                if H_actual != H_miss:
                    print(f"      FOUND: c₁={c1}, c₂={c2} but H={H_actual} ≠ {H_miss}")
                    break
    else:
        print(f"      NO tournament has c₁+2c₂ = {target}!")

print()
print("=" * 70)
print("PART 6: Structural impossibility")
print("=" * 70)

# For H=7: need c₁+2c₂ = 3. We know c₁ = t₃+t₅.
# From the data: t₃+t₅=3 DOES occur (at n=6), and c₂ must be 0.
# But let's check: when t₃+t₅=3, is c₂ always > 0?

print("\n  Tournaments with t₃+t₅=3:")
count_target = 0
for T in all_tournaments(n):
    t3 = count_3cycles(T, n)
    t5 = count_directed_cycles(T, n, 5)
    if t3 + t5 == 3:
        count_target += 1
        if count_target <= 5:
            ip = full_ip_decomposition(T, n)
            H = count_hp(T, n)
            c2 = ip[2] if len(ip) > 2 else 0
            print(f"    t₃={t3}, t₅={t5}, c₂={c2}, H={H}")

if count_target > 5:
    print(f"    ... ({count_target} total)")

# Complete check: for all t₃+t₅=3 tournaments, what is c₂?
c2_at_3 = []
for T in all_tournaments(n):
    t3 = count_3cycles(T, n)
    t5 = count_directed_cycles(T, n, 5)
    if t3 + t5 == 3:
        ip = full_ip_decomposition(T, n)
        c2 = ip[2] if len(ip) > 2 else 0
        c2_at_3.append(c2)

if c2_at_3:
    print(f"\n  When t₃+t₅=3, c₂ values: {sorted(set(c2_at_3))}")
    print(f"  c₂=0 ever? {'YES' if 0 in c2_at_3 else 'NO — THIS IS WHY H=7 IS IMPOSSIBLE!'}")

print("\n\nDone!")
