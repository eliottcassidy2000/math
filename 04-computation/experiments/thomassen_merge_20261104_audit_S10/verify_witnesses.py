#!/usr/bin/env python3
"""Direct (Python, no shared code) checks of single objects quoted in THM-4870 / HYP-9462:
  * the 14-edge n = 10 witness  C_10 + {1,6},{2,8},{3,7},{4,9};
  * the n = 13 counterexample quoted in HYP-9462 Results (rows from thomassen_merge_20261104_confirm.out);
  * planarity of the square-sum graphs Q_N, N <= 40 (networkx), for 'Q_N planar iff N <= 24'.
"""
import itertools, math
import networkx as nx

def ham_cycles(vertices, adj):
    """all Hamiltonian cycles (as frozensets of edges) by DFS from the smallest vertex"""
    vs = sorted(vertices)
    if len(vs) < 3:
        return []
    s = vs[0]
    res = set()
    def dfs(path, visited):
        if len(path) == len(vs):
            if s in adj[path[-1]]:
                res.add(frozenset(frozenset((path[i], path[(i + 1) % len(path)])) for i in range(len(path))))
            return
        for w in adj[path[-1]]:
            if w in vertices and w not in visited:
                path.append(w); visited.add(w)
                dfs(path, visited)
                path.pop(); visited.remove(w)
    dfs([s], {s})
    return list(res)

def analyse(n, edges, label):
    adj = {v: set() for v in range(n)}
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    V = set(range(n))
    hc = ham_cycles(V, adj)
    nc = []
    for v in range(n):
        sub = {u: adj[u] - {v} for u in V - {v}}
        if ham_cycles(V - {v}, sub):
            nc.append(v)
    cyc_ok = all(((i + 1) % n) in adj[i] for i in range(n))
    consec = [(v, (v + 1) % n) for v in nc if (v + 1) % n in nc]
    print("%s: n=%d edges=%d contains cycle 0..n-1: %s  #Ham cycles=%d  NC=%s  cycle-consecutive NC pairs=%s"
          % (label, n, len(edges), cyc_ok, len(hc), nc, consec))
    return adj, nc

# 1. the 14-edge witness
n = 10
E = [(i, (i + 1) % n) for i in range(n)] + [(1, 6), (2, 8), (3, 7), (4, 9)]
analyse(n, E, "THM-4870 (D) witness")

# 2. the n = 13 counterexample of HYP-9462 (B2)
rows = "1012 45 a 514 129 850 a2 140 298 500 a08 1420 801".split()
n = 13
R = [int(x, 16) for x in rows]
sym = all(((R[i] >> j) & 1) == ((R[j] >> i) & 1) for i in range(n) for j in range(n))
E = sorted({(min(i, j), max(i, j)) for i in range(n) for j in range(n) if R[i] >> j & 1})
print("n=13 rows symmetric:", sym, " loops:", any(R[i] >> i & 1 for i in range(n)))
adj, nc = analyse(n, E, "HYP-9462 n=13 witness")
for v, cyc in ((4, [0, 1, 2, 3, 10, 9, 8, 7, 6, 5, 11, 12]), (5, [0, 4, 3, 2, 1, 6, 7, 8, 9, 10, 11, 12])):
    ok = sorted(cyc) == sorted(set(range(n)) - {v}) and all(cyc[(i + 1) % len(cyc)] in adj[cyc[i]] for i in range(len(cyc)))
    print("  quoted Hamiltonian cycle of G - %d valid: %s" % (v, ok))

# 3. planarity of Q_N
def Q(N):
    G = nx.Graph()
    G.add_nodes_from(range(1, N + 1))
    for i in range(1, N + 1):
        for j in range(i + 1, N + 1):
            r = math.isqrt(i + j)
            if r * r == i + j:
                G.add_edge(i, j)
    return G
pl = [N for N in range(1, 41) if nx.check_planarity(Q(N))[0]]
print("planar Q_N for N in 1..40:", pl)
