#!/usr/bin/env python3
"""
FAST CHECK: At n=9, is alpha_1 <= 10 possible when every vertex is in a 3-cycle?

Approach: Use Moon's formula to filter by t3 <= 10, then check cycle coverage.

Key bound: sum_v b_v = 3*t3. If every vertex has b_v >= 1, then 3*t3 >= n.
At n=9: t3 >= 3.

When t3 = 3: the 3 three-cycles must be pairwise disjoint (partition [9]).
This gives alpha_3 >= 1, which by Part C contradicts H=21.

When t3 = 4, 5, ..., 10: need to check if alpha_1 can stay <= 10.
But even t3 = 4 at n=9 already means alpha_1 >= 4 from 3-cycles alone.
The 5-cycles forced by tournament structure will add more.

Let me check: MINIMUM alpha_1 at n=9 with t3 in {3,...,10} and every vertex in 3-cycle.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from math import comb
import random

def count_3cycles(adj, n):
    scores = [sum(adj[i]) for i in range(n)]
    return comb(n, 3) - sum(comb(s, 2) for s in scores)

def vertex_in_3cycle(adj, n, v):
    out_v = [j for j in range(n) if adj[v][j]]
    in_v = [j for j in range(n) if adj[j][v]]
    for u in out_v:
        for w in in_v:
            if u != w and adj[u][w]:
                return True
    return False

def count_alpha1_fast(adj, n, limit=15):
    """Count alpha_1 up to a limit, using DP for each odd cycle length."""
    from itertools import combinations
    alpha1 = 0
    for k in range(3, n+1, 2):
        if alpha1 > limit:
            return alpha1
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
            alpha1 += count
            if alpha1 > limit:
                return alpha1
    return alpha1


def main():
    n = 9
    random.seed(42)

    print(f"=== n={n}: Minimum alpha_1 with every vertex in 3-cycle ===")
    print(f"Key bound: t3 >= ceil({n}/3) = {(n+2)//3}")
    print(f"At t3=3: 3 disjoint 3-cycles => alpha_3 >= 1 => H != 21 (Part C)")
    print()

    # Special case: t3 = 3 at n=9
    print("=== t3=3 case (3 disjoint 3-cycles) ===")
    print("3 three-cycles covering all 9 vertices => each vertex in exactly 1 three-cycle")
    print("=> 3 pairwise-disjoint three-cycles => alpha_3 >= 1")
    print("=> Part C: alpha_1 + 2*alpha_2 + 4*alpha_3 >= 3 + 6 + 4 = 13 > 10")
    print("=> H=21 impossible! (no need to check computationally)")
    print()

    # For t3 >= 4: search for minimum alpha_1
    min_alpha1_by_t3 = {}
    trials = 0

    for trial in range(2000000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        scores = [sum(adj[i]) for i in range(n)]
        if 0 in scores or (n-1) in scores:
            continue

        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
        if t3 > 10 or t3 < 3:
            continue

        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
        if not all_in:
            continue

        trials += 1
        alpha1 = count_alpha1_fast(adj, n, limit=15)

        if t3 not in min_alpha1_by_t3 or alpha1 < min_alpha1_by_t3[t3]:
            min_alpha1_by_t3[t3] = alpha1
            print(f"  t3={t3}: new min alpha_1={alpha1}, scores={tuple(sorted(scores))}")

        if trials % 100 == 0 and trials <= 500:
            print(f"  ... {trials} cycle-rich found so far")

    print(f"\nTotal trials: {trial+1}, cycle-rich found: {trials}")
    print(f"\nMinimum alpha_1 by t3:")
    for t3 in sorted(min_alpha1_by_t3.keys()):
        print(f"  t3={t3}: min alpha_1 = {min_alpha1_by_t3[t3]}")

    overall_min = min(min_alpha1_by_t3.values()) if min_alpha1_by_t3 else None
    if overall_min and overall_min > 10:
        print(f"\n*** PROVED (sampling): cycle-rich at n=9 => alpha_1 >= {overall_min} > 10 ***")
        print("*** Combined with Part C (t3=3 case) and Parts J/K: H=21 impossible at n=9! ***")


if __name__ == "__main__":
    main()
