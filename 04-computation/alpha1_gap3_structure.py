#!/usr/bin/env python3
"""
Structural proof that alpha_1=3 is impossible in any tournament.

Key finding at n=5:
- c3 <= 2 => c5 = 0, so alpha_1 <= 2
- c3 = 3 => c5 >= 1, so alpha_1 >= 4
Therefore alpha_1 = 3 is impossible at n=5.

Question: does this extend to all n? At n>=6, alpha_1 = c3 + c5 + c7 + ...
Need: no combination of cycle counts can sum to exactly 3.

Hypothesis: c3 = 3 ALWAYS forces c5 >= 1 in ANY tournament on ANY n vertices.
If true, then alpha_1 = 3 requires c3 <= 2 and c5+c7+... = 1 or 3.
But at n >= 5, c5 depends on the tournament structure...

Let me check at n=6.

kind-pasteur-2026-03-06-S21
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import combinations

def count_cycles_by_length(T, n):
    result = {}
    for k in range(3, n+1, 2):
        total = 0
        for combo in combinations(range(n), k):
            verts = list(combo)
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if T[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and T[verts[vi]][verts[0]]:
                    total += c
        result[k] = total
    return result

# n=6: check all (c3, c5) pairs and whether c3=3 forces c5>=1
print("n=6: All (c3, c5) pairs:")
n = 6
m = n*(n-1)//2
pairs = {}
alpha1_achievable = set()
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles = count_cycles_by_length(T, n)
    c3, c5 = cycles[3], cycles[5]
    pair = (c3, c5)
    pairs[pair] = pairs.get(pair, 0) + 1
    alpha1_achievable.add(sum(cycles.values()))

for pair in sorted(pairs.keys()):
    c3, c5 = pair
    a1 = c3 + c5
    print(f"  (c3={c3}, c5={c5}): {pairs[pair]} tournaments, alpha_1>={a1}")

print(f"\nAchievable alpha_1 at n=6: {sorted(alpha1_achievable)}")
gaps = [x for x in range(max(alpha1_achievable)+1) if x not in alpha1_achievable]
print(f"Gaps: {gaps}")

# Check: does c3=3 force c5>=1 at n=6?
c3_3_has_c5_0 = any((3, 0) == pair for pair in pairs)
print(f"\n(c3=3, c5=0) exists at n=6? {c3_3_has_c5_0}")

# Check: what is min c5 when c3=3?
if any(pair[0] == 3 for pair in pairs):
    min_c5_at_c3_3 = min(pair[1] for pair in pairs if pair[0] == 3)
    print(f"Min c5 when c3=3 at n=6: {min_c5_at_c3_3}")

# More generally: for each c3, what is the minimum c5?
print(f"\nFor each c3, minimum c5 at n=6:")
c3_values = sorted(set(pair[0] for pair in pairs))
for c3 in c3_values:
    min_c5 = min(pair[1] for pair in pairs if pair[0] == c3)
    max_c5 = max(pair[1] for pair in pairs if pair[0] == c3)
    min_a1 = c3 + min_c5
    print(f"  c3={c3}: c5 in [{min_c5}, {max_c5}], min alpha_1 >= {min_a1}")

print("\nDone.")
