#!/usr/bin/env python3
"""tau(H) vs eps(H) for every odd-order subgroup class H of S_n, n <= NMAX, by enumerating the H-invariant labelled
objects (unions of H-orbits on pairs), canonicalising with nauty's labelg, and counting isomorphism types
(tournaments) resp. even isomorphism types (graphs; evenness from the automorphism group via networkx)."""
import sys, subprocess, itertools
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from collections import Counter
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8

def cyc(perm_cycles, n):
    g = list(range(n))
    for c in perm_cycles:
        for i in range(len(c)): g[c[i]] = c[(i + 1) % len(c)]
    return tuple(g)
def mul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def closure(gens, n):
    idt = tuple(range(n)); S = {idt}; fr = [idt]
    while fr:
        x = fr.pop()
        for g in gens:
            y = mul(g, x)
            if y not in S: S.add(y); fr.append(y)
    return S
def odd_partitions(n):
    def parts(n, m):
        if n == 0: yield (); return
        for k in range(min(n, m), 0, -1):
            if k % 2:
                for r in parts(n - k, k): yield (k,) + r
    return [p for p in parts(n, n) if p != tuple([1] * n)]
def perm_of_type(lam, n):
    cycles = []; v = 0
    for L in lam:
        cycles.append(list(range(v, v + L))); v += L
    return cyc(cycles, n)
def subgroup_reps(n):
    reps = {}
    for lam in odd_partitions(n):
        g = perm_of_type(lam, n); reps[f"Z{len(closure([g],n))} type {lam}"] = closure([g], n)
    if n >= 6:
        a = cyc([[0, 1, 2]], n); b = cyc([[3, 4, 5]], n); reps["Z3xZ3 = <(012),(345)>"] = closure([a, b], n)
    if n >= 7:
        a = cyc([[0, 1, 2, 3, 4, 5, 6]], n); b = cyc([[1, 2, 4], [3, 6, 5]], n)  # x->x+1, x->2x on Z_7
        H = closure([a, b], n); assert len(H) == 21; reps["Z7:Z3 (AGL(1,7) odd part)"] = H
    if n >= 9:
        a = cyc([[0, 1, 2]], n); b = cyc([[3, 4, 5], [6, 7, 8]], n); reps["Z3xZ3 = <(012),(345)(678)>"] = closure([a, b], n)
        a = cyc([[0, 1, 2], [3, 4, 5]], n); b = cyc([[0, 1, 2], [6, 7, 8]], n); reps["Z3xZ3 = <(012)(345),(012)(678)>"] = closure([a, b], n)
        a = cyc([[0, 1, 2], [3, 4, 5], [6, 7, 8]], n); b = cyc([[0, 3, 6], [1, 4, 7], [2, 5, 8]], n); reps["Z3xZ3 regular"] = closure([a, b], n)
        a = cyc([[0, 1, 2]], n); b = cyc([[3, 4, 5]], n); c = cyc([[6, 7, 8]], n); reps["Z3^3 = <(012),(345),(678)>"] = closure([a, b, c], n)
        d = cyc([[0, 3, 6], [1, 4, 7], [2, 5, 8]], n); H = closure([a, b, c, d], n); assert len(H) == 81; reps["Z3 wr Z3 (Sylow-3 of S_9, order 81)"] = H
        e = cyc([[0, 1, 2], [3, 5, 4]], n); H = closure([a, b, c, d], n)
        H27 = closure([cyc([[0,1,2],[3,4,5],[6,7,8]], n), d, cyc([[0,1,2],[3,5,4]], n)], n); reps[f"order-{len(H27)} subgroup <(012)(345)(678), d, (012)(354)>"] = H27
    return reps
def pair_orbits(H, n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    seen = set(); orbs = []
    for p in pairs:
        if p in seen: continue
        O = set()
        for g in H:
            q = (min(g[p[0]], g[p[1]]), max(g[p[0]], g[p[1]])); O.add(q)
        seen |= O; orbs.append(sorted(O))
    return orbs
def graph6(n, edges):
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(edges)
    return nx.to_graph6_bytes(G, header=False).decode().strip()
def labelg(lines):
    out = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
    assert len(out) == len(lines), (len(out), len(lines))
    return out
def digraph6_tournament(n, arcs):
    D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from(arcs)
    return nx.to_graph6_bytes(D, header=False).decode().strip() if False else to_d6(n, arcs)
def to_d6(n, arcs):
    # digraph6 encoding: '&' + N(n) + bits of the n*n adjacency matrix row by row
    A = [[0] * n for _ in range(n)]
    for u, v in arcs: A[u][v] = 1
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        s += chr(63 + int(''.join(map(str, bits[k:k + 6])), 2))
    return s

for n in range(3, NMAX + 1):
    # even/odd table of all graphs by canonical graph6
    out = subprocess.run(['geng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    canon = labelg(out)
    even_of = {}
    for s, c in zip(out, canon):
        G = nx.from_graph6_bytes(s.encode())
        edges = [(min(u, v), max(u, v)) for u, v in G.edges()]
        ev = True
        for m in GraphMatcher(G, G).isomorphisms_iter():
            s2 = 0
            for u, v in edges:
                if m[u] > m[v]: s2 ^= 1
            if s2: ev = False; break
        even_of[c] = ev
    print(f"\n=== n = {n}: graphs {len(out)}, even {sum(even_of.values())} ===")
    for name, H in subgroup_reps(n).items():
        orbs = pair_orbits(H, n)
        # tournaments invariant under H: orient each orbit consistently -- possible iff no element swaps a pair (odd order: ok)
        # choose for each orbit a transversal of oriented pairs: orient p=(i,j) as i->j, propagate by H
        tlines = []; glines = []
        orb_arcs = []
        for O in orbs:
            p = O[0]
            arcs = set()
            for g in H:
                arcs.add((g[p[0]], g[p[1]]))
            orb_arcs.append(arcs)
        for choice in itertools.product([0, 1], repeat=len(orbs)):
            arcs = []
            edges = []
            for c, O, A in zip(choice, orbs, orb_arcs):
                if c: arcs += list(A); edges += O
                else: arcs += [(v, u) for u, v in A]
            tlines.append(to_d6(n, arcs)); glines.append(graph6(n, edges))
        tcan = set(labelg(tlines)); gcan = labelg(glines)
        gtypes = set(gcan)
        e_even = sum(1 for c in gtypes if even_of[c])
        flag = "" if len(tcan) == e_even else "  <-- DIFFER"
        print(f"  H = {name} (|H|={len(H)}, pair orbits {len(orbs)}): tau = {len(tcan)}, eps = {e_even}, all graph types {len(gtypes)}{flag}")
