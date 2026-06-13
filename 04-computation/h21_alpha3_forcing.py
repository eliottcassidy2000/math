#!/usr/bin/env python3
"""
CRITICAL QUESTION: If m three-cycles cover all n vertices (n >= 9),
must there be 3 pairwise-disjoint ones?

This would prove H=21 impossible for n >= 9 via Part C.

The greedy argument: Take any 3-cycle C_1. Remove its 3 vertices.
On the remaining n-3 >= 6 vertices, the coverage condition means
some 3-cycle C_2 uses only remaining vertices (i.e., is disjoint from C_1).
Remove C_2's vertices. On the remaining n-6 >= 3 vertices, find C_3.

But does the coverage GUARANTEE that C_2 exists?

Counterexample attempt: All m three-cycles could share vertex v.
Then removing v kills all cycles. But n >= 9 with all cycles through v:
at most C(n-1, 2) three-cycles through v, covering {v} union pairs.
To cover all n vertices from pairs through v: need all n-1 other vertices
paired with v. But not all triples {v, a, b} are 3-cycles.
The coverage requires: for each w != v, some 3-cycle contains {v, w, x}.
This is possible.

So the issue is: CAN all 3-cycles share a common vertex?

If all m three-cycles share vertex v, and they cover all n vertices:
each other vertex appears in some 3-cycle {v, a, b}.
But then removing v: NO 3-cycle survives.
The remaining n-1 vertices have no 3-cycles? That depends on the tournament.

Actually, "the m three-cycles" are ALL directed 3-cycles, not just some.
If we say "every vertex is in some 3-cycle" and "all 3-cycles share v",
then removing v might create new 3-cycles? No, removing vertices can't create
new 3-cycles. Actually it can remove edges, potentially breaking existing
non-3-cycles... no. A 3-cycle on {a, b, c} with a, b, c != v exists in T
iff it exists in T-v.

So: if all 3-cycles contain v, then on T-v (n-1 vertices), there are NO 3-cycles.
By Part J: every vertex in T-v is not in any cycle. So T-v is acyclic = transitive.
And every vertex not equal to v is cycle-free in T-v.
But they ARE in cycles in T (they're in 3-cycles with v).

In this case: all directed odd cycles in T contain v. Every cycle
has length 3 (since the sub-tournament T-v is transitive, restricting
longer cycles to those involving v and the transitive rest).

Actually, can there be 5-cycles through v when T-v is transitive?
A 5-cycle through v on {v, a, b, c, d}: it visits v and 4 vertices
from the transitive part. In the transitive part, the ordering is
a linear order, say a < b < c < d. The 5-cycle needs arcs that go
"backwards" (from higher to lower in the ordering), and the only
backward arcs involve v. So the 5-cycle uses 2 arcs from v (one out, one in)
and 3 arcs among {a,b,c,d}. The 3 arcs among transitive vertices must all
go forward. But a directed 5-cycle on 5 vertices needs at least 2 backward
arcs (in the linear ordering). With only v providing backward connections...

Let me just check computationally.

Instance: kind-pasteur-2026-03-07-S33
"""

import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

from itertools import combinations
from collections import Counter
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

def find_3cycles(adj, n):
    """Find all directed 3-cycles as vertex frozensets (with multiplicity for direction)."""
    cycles = []
    for vs in combinations(range(n), 3):
        a, b, c = vs
        # Check a->b->c->a
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles.append(frozenset(vs))
        # Check a->c->b->a
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cycles.append(frozenset(vs))
    return cycles

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


def main():
    # Test 1: Can all 3-cycles share a common vertex?
    print("=== Test: All 3-cycles through one vertex ===")

    for n in [7, 8, 9]:
        print(f"\nn={n}:")

        # Construct: vertex 0 has N+(0) = {1,...,k}, N-(0) = {k+1,...,n-1}
        # Make T-{0} transitive: i beats j iff i > j (for i,j in {1,...,n-1})
        # Then ALL backward arcs go through vertex 0.
        # 3-cycles through 0: need u in N+(0), w in N-(0), u->w (backward).
        # In the transitive ordering of {1,...,n-1}: if i > j, then i beats j.
        # N+(0) = {1,...,k}, N-(0) = {k+1,...,n-1}.
        # Arc from u in {1,...,k} to w in {k+1,...,n-1}: u beats w iff u > w (in transitive order).
        # But u <= k and w >= k+1, so u < w, meaning w beats u.
        # So ALL arcs from N+(0) go TO N+(0)'s lower vertices... wait.
        #
        # In transitive tournament on {1,...,n-1}: vertex i beats j if i > j.
        # N+(0) = {1,...,k}: these are the LOWEST-score vertices.
        # An arc from u in N+(0) to w in N-(0): if u < w, then w beats u (w is higher).
        # So there are NO backward arcs from N+(0) to N-(0) in the transitive tournament.
        # That means b_0 = 0 and vertex 0 is in NO 3-cycle!
        #
        # To get 3-cycles through 0: need some arc u->w with u in N+(0), w in N-(0).
        # In the transitive tournament, ALL cross arcs go w -> u (higher beats lower).
        # Need to flip some arcs.

        # Let me construct it differently: make T-{0} transitive,
        # then add vertex 0 with specific arcs.
        # Transitive order: 1 < 2 < ... < n-1 (so n-1 beats all others).
        # Set N+(0) = {k+1, ..., n-1} (0 beats the high-ranked), N-(0) = {1,...,k}
        # Then 3-cycles through 0: w -> 0 -> u -> w for w in N-(0), u in N+(0), u->w.
        # u->w means higher beats lower, which happens in the transitive T.
        # u in {k+1,...,n-1}, w in {1,...,k}, u > w so u beats w. YES!
        # So b_0 = |N-(0)| * |N+(0)| cross arcs that go u->w. ALL of them go u->w.
        # b_0 = k * (n-1-k). Each gives a 3-cycle w -> 0 -> u -> w.

        # Set k = n//2 (balanced). Then b_0 ~ n^2/4 three-cycles through 0.
        # But we need FEW three-cycles. So make k small.

        # k = 1: N-(0) = {1}, N+(0) = {2,...,n-1}. b_0 = 1*(n-2) = n-2.
        # 3-cycles: {1, 0, j} for j = 2,...,n-1. That's n-2 three-cycles.
        # But these cover vertices {0, 1, 2, ..., n-1} = all n vertices.
        # All through vertex 0 and vertex 1!

        # Other 3-cycles? T-{0} is transitive, so NO other 3-cycles.
        # Total t3 = n-2.

        # Check: does each vertex v != 0,1 participate in a 3-cycle?
        # Vertex j (j >= 2): is in 3-cycle {1, 0, j}. Yes!
        # Vertex 1: is in 3-cycle {1, 0, 2}. Yes!
        # Vertex 0: is in 3-cycle {1, 0, 2}. Yes!
        # So ALL vertices in 3-cycles, all through vertex 0.

        # No source/sink? Score of 0 = n-2 (out to {2,...,n-1}).
        # Score of 1 = 0 in transitive part, but 1 beats 0? No, 1 is in N-(0), so 1 -> 0 is wrong...
        # Wait: N-(0) = {1}, meaning 1 -> 0 (1 beats 0). But in transitive part, 1 is the lowest.
        # Score of 1: beats 0, loses to 2,3,...,n-1. Score = 1.
        # Score of 0: loses to 1, beats 2,...,n-1. Score = n-2.
        # Score of j (j >= 2): beats 1,...,j-1 and loses to j+1,...,n-1 and beats 0? No, 0 beats j.
        # Actually: 0 -> j for j >= 2 (N+(0) = {2,...,n-1}).
        # And j -> i for i < j in transitive part.
        # Score of j: beats {1,...,j-1} (that's j-1 vertices), loses to {j+1,...,n-1} and to 0.
        # Wait, 0 beats j, so j loses to 0. Score of j = j-1.
        # Score sequence: 1 has score 1, 0 has score n-2, j has score j-1.
        # Sorted: (1-1, 2-1, ..., (n-1)-1, 1, n-2) = (0, 1, ..., n-2, 1, n-2)
        # Wait, vertex 1: score 1. Vertex 0: score n-2. Vertex j: score j-1.
        # Sorted scores: vertex 2 has score 1, vertex 1 has score 1, vertex 3 has score 2, ...
        # vertex n-1 has score n-2, vertex 0 has score n-2.
        # No source (min score = 1, from vertices 1 and 2), no sink.
        # Great!

        # So: tournament with ALL 3-cycles through vertex 0, every vertex in 3-cycle,
        # no source/sink, t3 = n-2. Alpha_1 = ?

        adj = [[0]*n for _ in range(n)]
        # Transitive part on {1,...,n-1}: i beats j if i > j
        for i in range(1, n):
            for j in range(1, i):
                adj[i][j] = 1
        # Vertex 0: loses to 1, beats 2,...,n-1
        adj[1][0] = 1  # 1 -> 0
        for j in range(2, n):
            adj[0][j] = 1  # 0 -> j

        t3 = count_3cycles(adj, n)
        scores = tuple(sorted([sum(adj[i]) for i in range(n)]))
        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))

        # Count all cycles
        alpha1 = 0
        comp = {}
        for k in range(3, n+1, 2):
            ck = 0
            for vs, d in find_directed_cycles_dp(adj, n, k):
                ck += d
            if ck > 0:
                comp[k] = ck
            alpha1 += ck

        print(f"  All-through-0: t3={t3}, scores={scores}, all_in={all_in}")
        print(f"  alpha_1={alpha1}, composition={comp}")

        # Check if any 3 pairwise-disjoint 3-cycles exist
        c3_sets = []
        for vs in combinations(range(n), 3):
            a, b, c = vs
            if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
                if frozenset(vs) not in c3_sets:
                    c3_sets.append(frozenset(vs))

        # Find max independent set in 3-cycle intersection graph
        max_indep = 0
        for size in range(min(len(c3_sets), 5), 0, -1):
            found = False
            for subset in combinations(range(len(c3_sets)), size):
                pairwise_disjoint = True
                for a in range(len(subset)):
                    for b in range(a+1, len(subset)):
                        if c3_sets[subset[a]] & c3_sets[subset[b]]:
                            pairwise_disjoint = False
                            break
                    if not pairwise_disjoint:
                        break
                if pairwise_disjoint:
                    max_indep = size
                    found = True
                    break
            if found:
                break

        print(f"  3-cycle vertex sets: {len(c3_sets)}")
        print(f"  Max pairwise-disjoint 3-cycles: {max_indep}")

    # Test 2: Can we have ALL 3-cycles sharing a vertex AND alpha_1 <= 10?
    print(f"\n=== Can all-sharing-vertex + alpha_1 <= 10 occur at n >= 9? ===")
    # From above: the all-through-0 construction has t3 = n-2.
    # At n=9: t3 = 7. alpha_1 = t3 + t5 + ... = ?
    # The transitive sub-tournament might still give 5-cycles through vertex 0.
    # Need alpha_1 <= 10, so t5 + t7 + ... <= 3.
    # This might work if the tournament is "nearly transitive" enough.

    # But the key: even if all 3-cycles share vertex 0,
    # there could be 5-cycles that DON'T go through vertex 0.
    # Wait: T-{0} is transitive, so no cycles at all in T-{0}.
    # Any cycle in T must go through vertex 0. So alpha_1 = total cycles through 0.

    # 5-cycles through 0 on {0, a, b, c, d}: the path a -> b -> c -> d (or similar)
    # must return to 0 via some arc d -> 0.
    # Since N-(0) = {1}: only vertex 1 has arc to 0. So d = 1.
    # 5-cycle: 0 -> a -> b -> c -> 1 -> 0 (with a, b, c in {2,...,n-1}).
    # Need: 0->a, a->b, b->c, c->1, 1->0.
    # 0->a: a in {2,...,n-1}. YES.
    # a->b: a > b in transitive order. So a > b.
    # b->c: b > c in transitive order. So b > c.
    # c->1: c > 1 in transitive order. So c >= 2.
    # 1->0: YES.
    # So: 5-cycle = 0 -> a -> b -> c -> 1 -> 0 with a > b > c >= 2.
    # Number of such 5-cycles: C(n-2, 3) = C(n-2, 3).
    # At n=9: C(7, 3) = 35. So t5 >= 35.
    # alpha_1 >= t3 + t5 >= 7 + 35 = 42 >> 10.

    # Wait, but these are only 5-cycles of ONE type.
    # There could also be: 0 -> a -> 1 -> b -> c -> 0.
    # Need: 0->a, a->1 (a > 1 in transitive, so a >= 2: a beats 1), 1->b... but 1 is lowest
    # in transitive order, so 1 beats nobody except 0. So 1->b requires b = 0, not allowed.
    # So only one type of 5-cycle exists.

    # At n=9: alpha_1 >= 42 >> 10. So the all-through-0 construction has alpha_1 >> 10.
    # To get alpha_1 <= 10 at n=9, we CAN'T have all 3-cycles through one vertex!

    print(f"\nFor all-through-0 construction:")
    for n in [7, 8, 9, 10]:
        t3 = n - 2
        t5 = comb(n-2, 3)
        print(f"  n={n}: t3={t3}, t5>={t5}, alpha_1>={t3+t5}")


if __name__ == "__main__":
    main()
