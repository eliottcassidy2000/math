#!/usr/bin/env python3
"""Exhaustive check on ALL labelled tournaments of order <= 6 (H) and <= 5 (Aut):
H(T) = prod H(strong components) and |Aut(T)| = prod |Aut(strong components)|."""
import itertools
import networkx as nx

def tourns(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        out = [set() for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: out[i].add(j)
            else: out[j].add(i)
        yield out

def H(vs, out):
    vs = list(vs); idx = {v: k for k, v in enumerate(vs)}; m = len(vs)
    dp = [[0] * m for _ in range(1 << m)]
    for k in range(m): dp[1 << k][k] = 1
    for mask in range(1, 1 << m):
        for k in range(m):
            c = dp[mask][k]
            if not c: continue
            for w in out[vs[k]]:
                if w in idx and not mask >> idx[w] & 1:
                    dp[mask | 1 << idx[w]][idx[w]] += c
    return sum(dp[(1 << m) - 1])

def aut(vs, out):
    vs = list(vs); c = 0
    for p in itertools.permutations(vs):
        mp = dict(zip(vs, p))
        if all((mp[w] in out[mp[v]]) == (w in out[v]) for v in vs for w in vs if v != w):
            c += 1
    return c

def comps(n, out):
    G = nx.DiGraph(); G.add_nodes_from(range(n))
    G.add_edges_from((v, w) for v in range(n) for w in out[v])
    return [sorted(c) for c in nx.strongly_connected_components(G)]

badH = badA = tot = 0
for n in range(1, 7):
    for out in tourns(n):
        cs = comps(n, out)
        hp = 1
        for c in cs: hp *= H(c, out)
        if hp != H(range(n), out): badH += 1
        if n <= 5:
            ap = 1
            for c in cs: ap *= aut(c, out)
            if ap != aut(range(n), out): badA += 1
        tot += 1
print("labelled tournaments checked:", tot, " H-multiplicativity failures:", badH, " |Aut|-product failures (n<=5):", badA)
