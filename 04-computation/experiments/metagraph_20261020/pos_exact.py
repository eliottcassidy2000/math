#!/usr/bin/env python3
"""Exact palindromic counts via ALL involutory anti-automorphisms (isomorphisms T -> T^rev with sigma^2 = 1, found by
networkx), half-path DP. (A) circulant tournaments of composite odd order 9 and 15 (the affine-only count missed non-affine
anti-automorphisms at n = 9); (B) every strong tournament on 7 vertices (gentourng -c): which are points of symmetry;
(C) every strong tournament on 6 and 8 vertices? (8: 6008 strong classes, feasible) -- even-order strong points of symmetry."""
import itertools, subprocess, sys
from math import gcd
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
def ham_count(n, adj):
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in adj[v]:
                if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[full])
def pal_count(n, adjm, adj):
    D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from((u, v) for u in range(n) for v in range(n) if adjm[u][v])
    R = D.reverse(copy=True)
    total = 0; nanti = 0
    for mp in DiGraphMatcher(D, R).isomorphisms_iter():
        sig = [mp[i] for i in range(n)]
        if any(sig[sig[x]] != x for x in range(n)): continue
        nanti += 1
        fixed = [x for x in range(n) if sig[x] == x]
        if n % 2 == 1 and len(fixed) != 1: continue
        if n % 2 == 0 and len(fixed) != 0: continue
        orbits = []; seen = set()
        for x in range(n):
            if sig[x] == x or x in seen: continue
            seen |= {x, sig[x]}; orbits.append((x, sig[x]))
        k = len(orbits); oid = {}
        for i, (a, b) in enumerate(orbits): oid[a] = i; oid[b] = i
        dp = {}
        for a, b in orbits:
            for v in (a, b): dp[(1 << oid[v], v)] = dp.get((1 << oid[v], v), 0) + 1
        full = (1 << k) - 1; cnt = 0
        for mask in range(1, 1 << k):
            for v in range(n):
                c0 = dp.get((mask, v), 0)
                if not c0: continue
                if mask == full:
                    if n % 2 == 1:
                        if adjm[v][fixed[0]]: cnt += c0
                    else:
                        if adjm[v][sig[v]]: cnt += c0   # middle arc v_{n/2} -> sigma(v_{n/2})
                    continue
                for w in adj[v]:
                    if sig[w] != w and not (mask >> oid[w]) & 1:
                        dp[(mask | (1 << oid[w]), w)] = dp.get((mask | (1 << oid[w]), w), 0) + c0
        total += cnt
    return total, nanti
def circulant_sets(n):
    units = [a for a in range(1, n) if gcd(a, n) == 1]
    classes = list(range(1, (n + 1) // 2)); seen = set(); reps = []
    for signs in itertools.product([0, 1], repeat=len(classes)):
        S = frozenset((x if s == 0 else n - x) for x, s in zip(classes, signs))
        canon = min(tuple(sorted(a * x % n for x in S)) for a in units)
        if canon in seen: continue
        seen.add(canon); reps.append(S)
    return reps
print("(A) circulants of composite odd order, exact anti-automorphism enumeration")
for n in (9, 15):
    for S in circulant_sets(n):
        adjm = [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
        adj = [[j for j in range(n) if adjm[i][j]] for i in range(n)]
        H = ham_count(n, adj); pal, nanti = pal_count(n, adjm, adj)
        print(f"  n = {n} S = {sorted(S)}: H = {H}, palindromic = {pal}, involutory anti-automorphisms {nanti}{'  <-- POS' if pal == H else ''}", flush=True)
print("(B)/(C) strong tournaments on 6, 7, 8 vertices: points of symmetry")
for n in (6, 7, 8):
    out = subprocess.run(['gentourng', '-c', '-q', str(n)], capture_output=True, text=True).stdout.split()
    pos = []
    for line in out:
        adjm = [[0] * n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[k] == '1': adjm[i][j] = 1
                else: adjm[j][i] = 1
                k += 1
        adj = [[j for j in range(n) if adjm[i][j]] for i in range(n)]
        # quick filter: self-converse needed
        D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from((u, v) for u in range(n) for v in range(n) if adjm[u][v])
        if not nx.is_isomorphic(D, D.reverse(copy=True)): continue
        H = ham_count(n, adj); pal, nanti = pal_count(n, adjm, adj)
        if pal == H: pos.append((H, nanti))
    print(f"  n = {n}: strong classes {len(out)}, strong points of symmetry: {len(pos)} {pos}", flush=True)
