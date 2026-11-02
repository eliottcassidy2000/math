#!/usr/bin/env python3
"""prime_path_scan.py (klein-2026-11-01) -- incremental scan of the three necessary conditions for a Hamiltonian path
through the primes below N in the square-sum graph (consecutive sums perfect squares): no isolated prime, at most two
leaves, and bipartite balance | |A| - |B| | <= 1 with A = {2} u {p = 1 mod 4}, B = {p = 3 mod 4}, i.e. Chebyshev's bias
pi(N;3,4) - pi(N;1,4) in {0, 1, 2}.  N runs over primes (N = the prime following the last vertex).  Reports the first
passing N, the number of passing N below the bound, and the bias statistics; then searches the first few passing N for
a path with a node budget.  usage: python3 prime_path_scan.py [BOUND=120000]"""
import math, sys
sys.setrecursionlimit(20000)
BOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 120000
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
P = primes_below(BOUND); Ps = set(P)
deg = {}; iso = 0; leaves = 0; bias = 0
def bump(v):
    global iso, leaves
    d = deg[v]
    if d == 0: iso -= 1
    elif d == 1: leaves -= 1
    deg[v] = d + 1
    if d + 1 == 1: leaves += 1
passing = []; cycle_cands = []; bias_small = 0; first_bias_zero = None
for i, p in enumerate(P):
    deg[p] = 0; iso += 1
    if p % 4 == 3: bias += 1
    elif p % 4 == 1: bias -= 1
    s = math.isqrt(p) + 1
    while s * s < 2 * p:
        q = s * s - p
        if q in Ps: bump(q); bump(p)
        s += 1
    if i + 1 < len(P):
        N = P[i + 1]
        if bias in (0, 1, 2):
            bias_small += 1
            if bias == 0 and first_bias_zero is None: first_bias_zero = N
        if iso == 0 and leaves <= 2 and bias in (0, 1, 2):
            passing.append((N, i + 1, leaves, bias))
        if iso == 0 and leaves == 0 and bias == 1:
            cycle_cands.append(N)
print('bound %d: %d primes; N with bias in {0,1,2}: %d (first bias 0 at N = %s); N passing all three conditions: %d' % (
    BOUND, len(P), bias_small, first_bias_zero, len(passing)))
print('first passing N (N, #primes, leaves, bias):', passing[:12])
print('last passing N below the bound:', passing[-1] if passing else None, '; Hamiltonian-cycle candidates (bias 1, no leaf, no isolated):', len(cycle_cands), cycle_cands[:5])
which_leaf = sorted(v for v in deg if deg[v] == 1)
print('leaf primes at the bound:', which_leaf)
print('primes still isolated at the bound:', end=' ')
# which primes are still isolated at the bound
print(sorted(v for v in deg if deg[v] == 0)[:10], '... count', sum(1 for v in deg if deg[v] == 0), '; leaves at the bound:', leaves, '; bias at the bound:', bias)
# search the first passing N
def sq_graph(Pn):
    S = set(Pn); adj = {p: set() for p in Pn}; top = max(Pn)
    for p in Pn:
        s = math.isqrt(p + 2)
        while s * s <= p + top:
            q = s * s - p
            if q > p and q in S: adj[p].add(q); adj[q].add(p)
            s += 1
    return adj
def ham_path(adj, budget):
    n = len(adj); best = [None]; nodes = [0]
    verts = sorted(adj, key=lambda v: len(adj[v]))
    starts = [v for v in verts if len(adj[v]) == 1] or verts[:8]
    def dfs(path, used):
        nodes[0] += 1
        if nodes[0] > budget: return 'budget'
        if len(path) == n: best[0] = list(path); return 'found'
        v = path[-1]
        cand = [w for w in adj[v] if w not in used]
        ends = [w for w in cand if len(adj[w] - used) <= 1]
        if len(ends) > 1: return None
        if ends: cand = ends
        cand.sort(key=lambda w: len(adj[w] - used))
        for w in cand:
            path.append(w); used.add(w)
            r = dfs(path, used)
            if r: return r
            path.pop(); used.discard(w)
        return None
    for s in starts:
        r = dfs([s], {s})
        if r == 'found': return best[0], nodes[0]
        if r == 'budget': return None, nodes[0]
    return None, nodes[0]
for N, k, lv, b in passing[:6]:
    if k > 3000:
        print('N = %d (%d primes): search deferred to prime_path_posa.py (rotation-extension)' % (N, k)); continue
    Pn = P[:k]; adj = sq_graph(Pn)
    path, nodes = ham_path(adj, 2_000_000)
    if path:
        assert len(set(path)) == len(Pn) and all(math.isqrt(a + c) ** 2 == a + c for a, c in zip(path, path[1:]))
        print('N = %d (%d primes): HAMILTONIAN PATH FOUND after %d nodes; first 12 terms %s ... last 6 %s' % (N, k, nodes, path[:12], path[-6:]))
        open('prime_path_N%d.txt' % N, 'w').write(' '.join(map(str, path)) + '\n')
    else:
        print('N = %d (%d primes): no path within %d nodes (UNKNOWN)' % (N, k, nodes))
