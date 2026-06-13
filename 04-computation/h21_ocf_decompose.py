#!/usr/bin/env python3
"""
Verify: H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

When alpha_1=10, i_2=0, this means all 10 cycles pairwise share a vertex.
But alpha_2 = # independent sets of size 2 = # vertex-disjoint pairs = i_2 = 0.
So H = 1 + 20 + 0 + 0 + ... = 21? Or is alpha_2 != i_2?

Let me verify the relationship.

In the independence polynomial I(G,x) = sum alpha_k x^k:
- alpha_0 = 1
- alpha_1 = |V(Omega)| = number of odd cycles
- alpha_2 = number of independent sets of size 2 = number of non-adjacent pairs
  = number of vertex-disjoint cycle pairs (since adjacency = sharing vertex)

So YES, alpha_2 = i_2 = number of vertex-disjoint pairs.

And H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

So if alpha_1=10, alpha_2=0, alpha_3=0, then H=21.
But we found alpha_1=10, alpha_2=0 at n=8 with H=81!
This means alpha_k > 0 for some k >= 2, contradicting our i_2=0 finding.

BUT WAIT: i_2 is the number of disjoint PAIRS, which IS alpha_2 of the
independence polynomial. If i_2=0, then alpha_2=0 and alpha_3=0
(since alpha_3 >= 1 requires 3 mutually disjoint cycles, which requires
at least 1 disjoint pair, contradicting i_2=0).

So if alpha_1=10, i_2=0, then H = 1 + 20 = 21. EXACTLY.

But we measured H=81 for these tournaments! Something is wrong.
Let me verify by computing I(Omega, 2) directly.

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations
from collections import Counter
import random

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def find_5cycles_dp(adj, n):
    cycles = []
    for verts in combinations(range(n), 5):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        for S in range(1, 32):
            for i in range(5):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(5):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, 5):
            if (31, j) in dp and adj[v[j]][v[0]]:
                count += dp[(31, j)]
        num = count // 5
        for _ in range(num):
            cycles.append(frozenset(verts))
    return cycles

def find_7cycles_dp(adj, n):
    cycles = []
    for verts in combinations(range(n), 7):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        for S in range(1, 128):
            for i in range(7):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(7):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, 7):
            if (127, j) in dp and adj[v[j]][v[0]]:
                count += dp[(127, j)]
        num = count // 7
        for _ in range(num):
            cycles.append(frozenset(verts))
    return cycles


def independence_polynomial(cycles, x):
    """Compute I(Omega, x) where Omega is the conflict graph on cycles."""
    n = len(cycles)
    if n == 0:
        return 1

    # Build adjacency (conflict = share vertex)
    adj_omega = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if cycles[i] & cycles[j]:  # share vertex
                adj_omega[i][j] = True
                adj_omega[j][i] = True

    # Enumerate all independent sets
    total = 0
    for mask in range(1 << n):
        # Check if this is an independent set
        independent = True
        verts_in_set = []
        for i in range(n):
            if mask & (1 << i):
                verts_in_set.append(i)
        for a in range(len(verts_in_set)):
            for b in range(a+1, len(verts_in_set)):
                if adj_omega[verts_in_set[a]][verts_in_set[b]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            k = len(verts_in_set)
            total += x**k

    return total


def held_karp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full])


def main():
    n = 8
    random.seed(42)
    print(f"n={n}: Verify I(Omega,2) for alpha_1=10, i_2=0 tournaments")

    count = 0
    for trial in range(100000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        c3 = find_3cycles(adj, n)
        if len(c3) > 12:
            continue

        c5 = find_5cycles_dp(adj, n)
        partial = c3 + c5

        if len(partial) == 10:
            # Check i_2
            i2 = 0
            for a in range(10):
                for b in range(a+1, 10):
                    if not (partial[a] & partial[b]):
                        i2 += 1
            if i2 == 0:
                # Compute I(Omega, 2) using ONLY 3+5 cycles
                I_partial = independence_polynomial(partial, 2)
                # Now add 7-cycles
                c7 = find_7cycles_dp(adj, n)
                full = c3 + c5 + c7
                I_full = independence_polynomial(full, 2)
                H = held_karp(adj, n)

                print(f"  t3={len(c3)}, t5={len(c5)}, t7={len(c7)}")
                print(f"  I(Omega_35, 2) = {I_partial}")
                print(f"  I(Omega_full, 2) = {I_full}")
                print(f"  H = {H}")
                print(f"  H == I_full? {H == I_full}")
                print()

                count += 1
                if count >= 5:
                    break

    if count == 0:
        print("No alpha_1(3+5)=10, i2=0 found in 100k samples")


if __name__ == "__main__":
    main()
