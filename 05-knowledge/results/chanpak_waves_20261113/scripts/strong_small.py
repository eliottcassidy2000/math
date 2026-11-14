#!/usr/bin/env python3
"""Hamiltonian-path counts of strong tournaments of orders 3..6 (all labelled tournaments)."""
import itertools
def H(n, adj):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for S in range(1<<n):
        for v in range(n):
            c = dp[S][v]
            if c:
                for w in range(n):
                    if not S>>w & 1 and adj[v]>>w & 1: dp[S|1<<w][w] += c
    return sum(dp[(1<<n)-1])
def strong(n, adj):
    # reachability from 0 forward and backward
    def reach(nb):
        seen = 1; st = [0]
        while st:
            v = st.pop()
            for w in range(n):
                if nb(v, w) and not seen>>w & 1: seen |= 1<<w; st.append(w)
        return seen == (1<<n)-1
    return reach(lambda v,w: adj[v]>>w & 1) and reach(lambda v,w: adj[w]>>v & 1)
for n in range(3, 7):
    pairs = list(itertools.combinations(range(n), 2)); vals = set()
    for mask in range(1 << len(pairs)):
        adj = [0]*n
        for k, (i, j) in enumerate(pairs):
            if mask >> k & 1: adj[i] |= 1<<j
            else: adj[j] |= 1<<i
        if strong(n, adj): vals.add(H(n, adj))
    print(f"order {n}: strong H values {sorted(vals)}")
