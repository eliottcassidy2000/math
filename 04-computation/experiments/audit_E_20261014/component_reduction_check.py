#!/usr/bin/env python3
"""Audit E (2026-10-14): S35's component reduction used by THM-4655 1(d) and 4:  delta(T) = max over strong components
(singletons contribute 0), checked on every labelled tournament with n <= 6 by explicit Hamiltonian-path enumeration;
also checks H(T) = product of H(components)."""
from fractions import Fraction as F
from itertools import combinations

def paths(vs, adj):
    vs = list(vs); res = []
    def rec(p, used):
        if len(p) == len(vs):
            res.append(tuple(p)); return
        for v in vs:
            if v not in used and adj[p[-1]][v]:
                p.append(v); used.add(v); rec(p, used); used.discard(v); p.pop()
    for s in vs:
        rec([s], {s})
    return res

def delta(vs, adj):
    P = paths(vs, adj); H = len(P); best = F(0)
    if H == 0: return None, 0
    pos = [{v: k for k, v in enumerate(p)} for p in P]
    for x, y in combinations(vs, 2):
        c = sum(1 for q in pos if q[x] < q[y]); best = max(best, min(F(c, H), 1 - F(c, H)))
    return best, H

bad = 0; tot = 0
for n in range(2, 7):
    E = list(combinations(range(n), 2))
    for code in range(1 << len(E)):
        adj = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(E):
            if code >> k & 1: adj[i][j] = 1
            else: adj[j][i] = 1
        R = [[adj[i][j] or i == j for j in range(n)] for i in range(n)]
        for k in range(n):
            for i in range(n):
                if R[i][k]:
                    for j in range(n):
                        if R[k][j]: R[i][j] = True
        comps = []; seen = set()
        for v in range(n):
            if v in seen: continue
            C = [u for u in range(n) if R[v][u] and R[u][v]]; seen |= set(C); comps.append(C)
        d, H = delta(range(n), adj)
        dm = F(0); Hp = 1
        for C in comps:
            dc, Hc = delta(C, adj) if len(C) > 1 else (F(0), 1)
            dm = max(dm, dc); Hp *= Hc
        tot += 1
        if d != dm or H != Hp: bad += 1
print(f"labelled tournaments n = 2..6 checked: {tot}; violations of delta(T) = max_C delta(C) or H(T) = prod H(C): {bad}")
