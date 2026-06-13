#!/usr/bin/env python3
"""
At n=7, alpha_1=10 has two compositions:
  (t3=5, t5=4, t7=1) and (t3=6, t5=4, t7=0)

Check i_2 for each composition separately.
Also understand: the 7-cycle uses ALL 7 vertices, so it shares
with EVERY other cycle (which must use >= 3 of the 7 vertices).
So the 7-cycle never contributes to a disjoint pair.

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations
from collections import Counter
import random

def find_directed_cycles_dp(adj, n, k):
    result = []
    for verts in combinations(range(n), k):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        full = (1 << k) - 1
        for S in range(1, full + 1):
            for i in range(k):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(k):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, k):
            if (full, j) in dp and adj[v[j]][v[0]]:
                count += dp[(full, j)]
        if count > 0:
            result.append((frozenset(verts), count))
    return result


def get_all_cycles(adj, n):
    all_cycles = []
    for k in range(3, n + 1, 2):
        if k > n:
            break
        for vs, d in find_directed_cycles_dp(adj, n, k):
            for _ in range(d):
                all_cycles.append(vs)
    return all_cycles


def compute_i2(all_cycles):
    i2 = 0
    for a in range(len(all_cycles)):
        for b in range(a+1, len(all_cycles)):
            if not (all_cycles[a] & all_cycles[b]):
                i2 += 1
    return i2


def main():
    n = 7
    random.seed(42)

    comp_541 = []
    comp_640 = []

    for trial in range(100000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        all_cyc = get_all_cycles(adj, n)
        if len(all_cyc) != 10:
            continue

        c3 = [c for c in all_cyc if len(c) == 3]
        c5 = [c for c in all_cyc if len(c) == 5]
        c7 = [c for c in all_cyc if len(c) == 7]

        i2 = compute_i2(all_cyc)

        if len(c3) == 5 and len(c5) == 4 and len(c7) == 1:
            comp_541.append(i2)
        elif len(c3) == 6 and len(c5) == 4 and len(c7) == 0:
            comp_640.append(i2)

    print(f"n=7: alpha_1=10 composition analysis (100k samples)")
    print(f"  (5,4,1): {len(comp_541)} found, i_2 distribution: {dict(Counter(comp_541))}")
    print(f"  (6,4,0): {len(comp_640)} found, i_2 distribution: {dict(Counter(comp_640))}")

    # The 7-cycle shares with everything (uses all 7 vertices),
    # so it doesn't add disjoint pairs.
    # In (5,4,1): 5 three-cycles + 4 five-cycles + 1 seven-cycle.
    # Disjoint pairs can only be (3,3) (since 3+5=8>7, 3+7=10>7, 5+7=12>7).
    # So i_2 = number of disjoint (3,3) pairs among the 5 three-cycles.
    # With 5 three-cycles on 7 vertices: can they have 0 disjoint pairs?
    # 5 pairwise-intersecting 3-subsets of [7].
    # EKR: C(6,2)=15 is max through common element. 5 << 15. Easy.
    # Hilton-Milner: 13 without common. Also easy.

    print(f"\n  At n=7: disjoint (3,3) pairs are the ONLY possible disjoint type")
    print(f"  3+5=8>7, 3+7=10>7, 5+5=10>7, 5+7=12>7 => all non-3,3 pairs share")
    print(f"  So i_2 = # disjoint pairs among 3-cycles only")

    # For (5,4,1): 5 three-cycles. i_2 distribution shows the answer.
    # For (6,4,0): 6 three-cycles. i_2 distribution shows the answer.

    # KEY QUESTION: Can 5 pairwise-intersecting directed 3-cycles exist at n=7?
    # If yes, then (5,4,1) with i_2=0 is possible (at the 3-cycle level).
    # But we also need the total alpha_1=10, meaning exactly 4 five-cycles and 1 seven-cycle.

    # From data: (5,4,1) always has i_2=? Let's check.
    if comp_541:
        min_i2 = min(comp_541)
        print(f"\n  (5,4,1) min i_2: {min_i2}")
    if comp_640:
        min_i2 = min(comp_640)
        print(f"  (6,4,0) min i_2: {min_i2}")


if __name__ == "__main__":
    main()
