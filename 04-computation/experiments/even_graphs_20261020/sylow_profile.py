#!/usr/bin/env python3
"""Sylow-p profiles (THM-4723 Theorem P / the 2026-10-20 packet's Sylow action theorem): for p = 3, 5, 7 and n <= NMAX,
the distribution over isomorphism types of the S_n-conjugacy class of a Sylow p-subgroup of Aut (invariant: order and
multiset of cycle types) for tournaments versus even graphs (and odd graphs for comparison)."""
import sys, subprocess
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher, DiGraphMatcher
from collections import Counter
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
def mul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def cycle_type(g):
    n = len(g); seen = [False] * n; ct = []
    for i in range(n):
        if not seen[i]:
            L = 0; j = i
            while not seen[j]: seen[j] = True; j = g[j]; L += 1
            ct.append(L)
    return tuple(sorted(ct, reverse=True))
def closure(gens, n):
    idt = tuple(range(n)); S = {idt}; fr = [idt]
    while fr:
        x = fr.pop()
        for g in gens:
            y = mul(g, x)
            if y not in S: S.add(y); fr.append(y)
    return S
def vp(m, p):
    k = 0
    while m % p == 0: m //= p; k += 1
    return k
def sylow_invariant(A, p, n):
    target = p ** vp(len(A), p)
    if target == 1: return (1, ())
    pel = [g for g in A if all(L in (1, p) or L % p == 0 for L in cycle_type(g)) and all(L % p == 0 or L == 1 for L in cycle_type(g))]
    # p-elements: order a power of p  <=> every cycle length is a power of p
    def is_ppow(L):
        while L % p == 0: L //= p
        return L == 1
    pel = [g for g in A if all(is_ppow(L) for L in cycle_type(g))]
    P = {tuple(range(n))}
    for g in pel:
        if g in P: continue
        Q = closure(list(P) + [g], n) if len(P) > 1 else closure([g], n)
        if len(Q) == p ** vp(len(Q), p) and len(Q) <= target:
            P = Q
            if len(P) == target: break
    assert len(P) == target, (len(P), target)
    return (len(P), tuple(sorted(Counter(cycle_type(g) for g in P).items())))
for n in range(3, NMAX + 1):
    gs = [nx.from_graph6_bytes(s.encode()) for s in subprocess.run(['geng', '-q', str(n)], capture_output=True, text=True).stdout.split()]
    ts = []
    for line in subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split():
        D = nx.DiGraph(); D.add_nodes_from(range(n)); k = 0
        for i in range(n):
            for j in range(i + 1, n):
                D.add_edge(i, j) if line[k] == '1' else D.add_edge(j, i); k += 1
        ts.append(D)
    dist = {p: {'T': Counter(), 'E': Counter(), 'O': Counter()} for p in (3, 5, 7)}
    for G in gs:
        E = [(min(u, v), max(u, v)) for u, v in G.edges()]
        A = {tuple(m[i] for i in range(n)) for m in GraphMatcher(G, G).isomorphisms_iter()}
        even = True
        for g in A:
            s = 0
            for u, v in E:
                if g[u] > g[v]: s ^= 1
            if s: even = False; break
        for p in (3, 5, 7):
            dist[p]['E' if even else 'O'][sylow_invariant(A, p, n)] += 1
    for T in ts:
        A = {tuple(m[i] for i in range(n)) for m in DiGraphMatcher(T, T).isomorphisms_iter()}
        for p in (3, 5, 7): dist[p]['T'][sylow_invariant(A, p, n)] += 1
    print(f"\n=== n = {n} ===")
    for p in (3, 5, 7):
        keys = sorted(set(dist[p]['T']) | set(dist[p]['E']), key=lambda k: (k[0], k[1]))
        ok = all(dist[p]['T'][k] == dist[p]['E'][k] for k in keys)
        print(f" p = {p}: Sylow-p profile tournaments = even graphs? {ok}")
        for k in keys:
            if k[0] > 1:
                print(f"    |P| = {k[0]} types {k[1]}: T {dist[p]['T'][k]}  E {dist[p]['E'][k]}  (odd graphs {dist[p]['O'][k]})")
