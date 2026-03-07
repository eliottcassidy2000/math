#!/usr/bin/env python3
"""
What is the maximum n for which alpha_1 <= 10 is possible
when every vertex is in some 3-cycle?

Moon's formula: t3 = C(n,3) - sum_v C(s_v, 2)
For t3 <= 10: sum_v C(s_v, 2) >= C(n,3) - 10.

By convexity: sum C(s_v, 2) >= n * C(S/n, 2) where S = n(n-1)/2.
S/n = (n-1)/2. So sum >= n * C((n-1)/2, 2) = n*(n-1)/2*(n-3)/2 / ... hmm.

Actually C((n-1)/2, 2) = (n-1)(n-3)/8 (for n odd).
So sum >= n*(n-1)*(n-3)/8.
And C(n,3) = n(n-1)(n-2)/6.

t3 >= C(n,3) - sum C(s_v,2). For regular: t3 = C(n,3) - n*C((n-1)/2, 2).

For t3 to be small we need scores near transitive. But "every vertex in 3-cycle"
means no vertex has score 0 or n-1 (Part K), and more: each vertex must participate
in at least one 3-cycle. Actually Part J is stronger: if v not in any 3-cycle,
then v not in any cycle at all, meaning all arcs from N-(v) go to N+(v).

The constraint "every vertex in some 3-cycle" means: for each v, there exist
u in N+(v), w in N-(v) with u -> w (creating w -> v -> u -> w).

For vertex v with out-degree s_v: the number of 3-cycles through v is
s_v * (n-1-s_v) - (# arcs from N+(v) to N-(v)).
Wait, more precisely: 3-cycles through v = s_v*(n-1-s_v) - e(N+(v), N-(v))
where e(A,B) = # arcs from A to B.

Actually 3-cycles through v = # edges u->w with u in N+(v), w in N-(v), u->w.
Total possible = s_v * (n-1-s_v). But some go from N+(v) to N-(v) (the "wrong" way).
Let f_v = # arcs from N-(v) to N+(v) (the "transitive" direction).
Then 3-cycles through v = s_v*(n-1-s_v) - f_v? No...

The arcs between N+(v) and N-(v):
- N+(v) has s_v vertices, N-(v) has n-1-s_v vertices.
- Total pairs = s_v*(n-1-s_v).
- Let b_v = arcs from N+(v) to N-(v) (the "backward" direction, creating 3-cycles).
- Then arcs from N-(v) to N+(v) = s_v*(n-1-s_v) - b_v.
- 3-cycles through v = b_v (each such arc u->w with u in N+(v), w in N-(v) gives w->v->u->w).

For v to be in at least one 3-cycle: b_v >= 1, meaning at least one arc goes
from N+(v) back to N-(v).

If b_v = 0 for some v, then all arcs go N-(v) -> N+(v), which is Part J's "layered" structure.

KEY QUESTION: For large n, can we have t3 small AND b_v >= 1 for all v?

t3 = sum_v b_v / ... no. Actually by a well-known identity:
t3 = (1/3) * sum_v [s_v*(n-1-s_v) - (n-1-s_v choose 2) + (something)]

Hmm, let me just use Moon's formula directly and check score constraints.

For "every vertex in a 3-cycle" at large n with t3 small:
The key bound comes from the vertex with the most extreme score.
If min score = 1, then vertex v with s_v=1 has b_v = 1*(n-2) - f_v.
For b_v >= 1: f_v <= n-3. Since f_v <= 1*(n-2) = n-2, this is easy.
Actually f_v = # arcs from N-(v) to N+(v). |N+(v)|=1, |N-(v)|=n-2.
Arcs from N-(v) to N+(v): each of the n-2 vertices in N-(v) either beats
or loses to the single vertex in N+(v). So f_v can be 0 to n-2.
b_v = n-2 - f_v. For b_v >= 1: f_v <= n-3. This fails only if f_v = n-2,
meaning the unique out-neighbor of v beats ALL of N-(v).

For the score constraint at the OTHER end:
If max score = n-2, then vertex w with s_w=n-2 has |N+(w)|=n-2, |N-(w)|=1.
b_w = arcs from N+(w) to N-(w). N+(w) has n-2 vertices, N-(w) has 1.
b_w = # of the n-2 out-neighbors that LOSE to the single in-neighbor.
For b_w >= 1: at least one out-neighbor loses to the in-neighbor. Easy.

So the score constraints alone don't prevent "every vertex in 3-cycle" at any n.
The issue is: with t3 <= 10 AND every vertex in a 3-cycle, can we have
many 5-cycles, 7-cycles etc that keep alpha_1 low?

At large n, tournaments with few 3-cycles are nearly transitive, and nearly
transitive tournaments have even fewer long odd cycles. So alpha_1 should be
dominated by t3 for nearly-transitive tournaments.

Let me verify: for nearly-transitive n=10 tournaments with t3 ~ 5-10,
what is alpha_1?

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from collections import Counter
from math import comb
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


def count_3cycles(adj, n):
    """Fast 3-cycle count using Moon's formula."""
    scores = [sum(adj[i]) for i in range(n)]
    return comb(n, 3) - sum(comb(s, 2) for s in scores)


def vertex_in_3cycle(adj, n, v):
    """Check if vertex v is in at least one 3-cycle."""
    out_v = [j for j in range(n) if adj[v][j]]
    in_v = [j for j in range(n) if adj[j][v]]
    for u in out_v:
        for w in in_v:
            if u != w and adj[u][w]:
                return True
    return False


def main():
    print("=== Maximum n where alpha_1 <= 10 with every vertex in a 3-cycle ===")

    # Strategy: generate nearly-transitive tournaments at n=9,10,11,...
    # with t3 small and check alpha_1 and "every vertex in 3-cycle".

    for n in [8, 9, 10, 11, 12]:
        print(f"\nn={n}:")
        random.seed(42 + n)

        # Approach: start from transitive tournament and flip a few arcs
        # Transitive: vertex i beats vertices 0,...,i-1 (so i has out-degree i... wait)
        # Standard: vertex i beats j iff i > j. Score of i = i. Score seq = (0,1,...,n-1).

        found_alpha_le10 = 0
        found_all_in_3cycle = 0
        found_both = 0
        alpha_vals = []

        for trial in range(50000):
            # Start from transitive and flip k random arcs
            adj = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i):
                    adj[i][j] = 1  # i beats j

            # Flip 2-4 random arcs
            num_flips = random.randint(1, 4)
            flipped = set()
            for _ in range(num_flips):
                i = random.randint(0, n-1)
                j = random.randint(0, n-1)
                if i == j:
                    continue
                if (min(i,j), max(i,j)) in flipped:
                    continue
                flipped.add((min(i,j), max(i,j)))
                # Flip arc between i and j
                if adj[i][j]:
                    adj[i][j] = 0
                    adj[j][i] = 1
                else:
                    adj[j][i] = 0
                    adj[i][j] = 1

            t3 = count_3cycles(adj, n)
            if t3 > 10:
                continue

            # Check all vertices in 3-cycle
            all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))

            # Count alpha_1 (full cycle count)
            alpha1 = 0
            for k in range(3, n+1, 2):
                if k > n:
                    break
                for vs, d in find_directed_cycles_dp(adj, n, k):
                    alpha1 += d
                if alpha1 > 10:
                    break

            if alpha1 <= 10:
                found_alpha_le10 += 1
                if all_in:
                    found_both += 1
                    alpha_vals.append(alpha1)
                    scores = tuple(sorted([sum(adj[i]) for i in range(n)]))
                    if trial < 100 or found_both <= 5:
                        print(f"  alpha_1={alpha1}, t3={t3}, scores={scores}, all_in={all_in}")
            if all_in:
                found_all_in_3cycle += 1

        print(f"  Trials: 50000, alpha_1<=10: {found_alpha_le10}, all_in_3cycle: {found_all_in_3cycle}, both: {found_both}")
        if alpha_vals:
            print(f"  alpha_1 values when both: {Counter(alpha_vals)}")

    # Also check: for large n, what's the MINIMUM alpha_1 when every vertex in 3-cycle?
    print("\n=== Minimum alpha_1 with every vertex in 3-cycle (random sampling) ===")
    for n in [9, 10, 11, 12]:
        print(f"\nn={n}:")
        random.seed(1000 + n)
        min_alpha1 = float('inf')
        min_scores = None

        for trial in range(20000):
            adj = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1

            # Quick filter: check scores first
            scores = [sum(adj[i]) for i in range(n)]
            if 0 in scores or (n-1) in scores:
                continue  # has source/sink, skip

            t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
            if t3 > 15:  # relaxed threshold
                continue

            all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
            if not all_in:
                continue

            alpha1 = 0
            for k in range(3, n+1, 2):
                if k > n:
                    break
                for vs, d in find_directed_cycles_dp(adj, n, k):
                    alpha1 += d
                if alpha1 > 15:
                    break

            if alpha1 < min_alpha1:
                min_alpha1 = alpha1
                min_scores = tuple(sorted(scores))
                print(f"  New min alpha_1={alpha1}, t3={t3}, scores={min_scores}")

        print(f"  Overall min alpha_1: {min_alpha1}")


if __name__ == "__main__":
    main()
