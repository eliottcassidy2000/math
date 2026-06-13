#!/usr/bin/env python3
"""
KEY QUESTION: Does "every vertex in a 3-cycle" force t3 to be large enough?

At n=9 with scores (1,1,1,3,4,5,7,7,7): t3=2.
But do the vertices with score 1 necessarily participate in those 2 three-cycles?

Vertex with score 1: out-degree 1, in-degree 7.
For it to be in a 3-cycle: its single out-neighbor must beat some in-neighbor of v.
If the out-neighbor beats ALL in-neighbors, v is in no 3-cycle (Part J).

Actually Part J says: if v is in no 3-cycle, then v is in no cycle at all.
So "every vertex in a 3-cycle" is equivalent to "every vertex in some cycle."

For a vertex with score 1: it's in a 3-cycle iff its out-neighbor has an arc
going "back" to some in-neighbor.

With t3=2 at n=9: only 2 three-cycles. Each covers 3 vertices.
They can cover at most 6 of the 9 vertices. So at least 3 vertices
are NOT in any 3-cycle. By Part J, those vertices are not in any cycle.

BUT: they might still be in 5-cycles or 7-cycles!
Wait no: Part J says NOT in 3-cycle => NOT in any cycle. Period.

So: if t3 covers only 6 of 9 vertices, then 3 vertices are cycle-free.
Part J allows removing them: H(T) = H(T'), where T' is the subtournament
on the 6 vertices in cycles. By Part G: H(T') != 21 for |V(T')| <= 8.

So: the CRITICAL QUESTION IS:
Can the vertices of all 3-cycles span ALL n vertices when t3 is small?

If t3 = k three-cycles, they cover at most 3k vertices.
For all n vertices covered: n <= 3k, so k >= ceil(n/3).
At n=9: k >= 3. At n=10: k >= 4. etc.

But we also need alpha_1 <= 10, meaning t3 + t5 + ... <= 10.
So t3 <= 10.

With t3 three-cycles covering all n vertices: n <= 3*t3 <= 30.
So n <= 30 is the absolute max.

But MORE IMPORTANTLY: the three-cycles through a vertex v with score s_v:
the number of 3-cycles through v = b_v (backward arcs from N+(v) to N-(v)).
For v to be in a 3-cycle: b_v >= 1.

Sum of b_v over all v = 3 * t3 (each 3-cycle contributes 1 to each of 3 vertices).
So sum b_v = 3 * t3.

If every vertex has b_v >= 1: sum b_v >= n, so 3*t3 >= n, t3 >= ceil(n/3).

COMBINED WITH alpha_1 <= 10:
t3 >= ceil(n/3) and t3 <= 10 gives n <= 30.

But actually t3 + t5 + ... = alpha_1 <= 10, and t3 >= ceil(n/3).
So ceil(n/3) <= 10, giving n <= 30.

For the induction to work, we need the base case n <= N_0 proved exhaustively,
and then for n > N_0, either some vertex not in 3-cycle (Part J reduces),
OR t3 >= ceil(n/3) is "high enough" that alpha_1 > 10.

But ceil(n/3) only exceeds 10 at n >= 31! We can't verify exhaustively
up to n=30.

DIFFERENT APPROACH: Can we show alpha_1 >= f(n) for some f(n) > 10 when n >= 9?

Actually: the 3-cycles must cover all n vertices, AND they're in a tournament.
The tournament structure forces additional cycles beyond just the 3-cycles.

KEY INSIGHT: If k three-cycles cover all n vertices, and n is large,
then the tournament has MANY additional 5-cycles forced by the structure.

At n=9 with every vertex in a 3-cycle:
t3 >= 3, and the three-cycles cover all 9 vertices.
But how many 5-cycles does this force?

Let me check computationally: at n=9, find tournaments with:
- Every vertex in a 3-cycle
- t3 <= 10
- alpha_1 <= 10

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
    print("=== At n=9: search for alpha_1 <= 10 with every vertex in 3-cycle ===")
    n = 9
    random.seed(42)

    found = 0
    checked = 0
    alpha_vals = []

    for trial in range(500000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        # Fast filter: score sequence
        scores = sorted([sum(adj[i]) for i in range(n)])
        if 0 in scores or (n-1) in scores:
            continue  # source/sink

        t3 = comb(n, 3) - sum(comb(s, 2) for s in scores)
        if t3 > 10:
            continue

        # Check every vertex in 3-cycle
        all_in = all(vertex_in_3cycle(adj, n, v) for v in range(n))
        if not all_in:
            checked += 1
            continue

        checked += 1

        # Count ALL cycles
        alpha1 = 0
        comp = {}
        for k in range(3, n+1, 2):
            if k > n:
                break
            ck = 0
            for vs, d in find_directed_cycles_dp(adj, n, k):
                ck += d
            comp[k] = ck
            alpha1 += ck
            if alpha1 > 15:
                break

        found += 1
        alpha_vals.append(alpha1)

        if alpha1 <= 10:
            print(f"  FOUND alpha_1={alpha1}! scores={tuple(scores)}, comp={comp}")

        if found <= 10 or (found % 100 == 0):
            if found <= 10:
                print(f"  #{found}: alpha_1={alpha1}, t3={t3}, scores={tuple(scores)}, comp={comp}")

    print(f"\nTrials: {trial+1}, with t3<=10 & every-in-3cycle: {found}")
    if alpha_vals:
        print(f"Alpha_1 range: [{min(alpha_vals)}, {max(alpha_vals)}]")
        print(f"Alpha_1 distribution: {dict(Counter(alpha_vals).most_common(20))}")
        if min(alpha_vals) > 10:
            print(f"\n*** PROVED: at n=9, every-vertex-in-3-cycle + t3<=10 => alpha_1 >= {min(alpha_vals)} > 10 ***")
            print("*** Combined with Part J: H=21 impossible at n=9! ***")
    else:
        print("No tournaments found with t3<=10 and every vertex in 3-cycle!")
        print("This means: at n=9, t3<=10 => some vertex NOT in 3-cycle")
        print("=> Part J applies => reduce to n<=8 => H != 21. QED for n=9!")


if __name__ == "__main__":
    main()
