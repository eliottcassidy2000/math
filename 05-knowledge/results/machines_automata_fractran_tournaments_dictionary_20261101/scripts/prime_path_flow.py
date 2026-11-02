#!/usr/bin/env python3
"""prime_path_flow.py (klein-2026-11-01) -- Hamiltonian path through the primes below N in the square-sum graph by
2-MATCHING + CYCLE MERGING.  The graph is bipartite (THM-4860 (a)): A = {2} u {p = 1 mod 4}, B = {p = 3 mod 4}.  A
Hamiltonian path has n - 1 edges and all degrees <= 2, so step 1 computes a maximum spanning subgraph with degrees <= 2
as a maximum flow (source -> A, capacity min(deg, 2); graph edges, capacity 1; B -> sink, capacity min(deg, 2)); if its
size nu_2 is below n - 1 there is no Hamiltonian path (THM-4860 (e)); otherwise its components are cycles and at most
one path (a 2-factor has one edge removed).  Step 2 merges every cycle into the path (or cycles together) by 2-opt
exchanges: remove (c, c+) from the cycle and (p, p+) from the path, add (c, p) and (c+, p+), or add (p, c+) and (c, p+).
The degree constraints make the restricted residue classes of THM-4860 (d'),(e) automatic, which is what the
rotation-extension solvers could not do.  The path is verified and written to data/prime_path_N<N>.txt.
usage: python3 prime_path_flow.py N [seed=1]"""
import math, sys, random, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow
N = int(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
rng = random.Random(seed); T0 = time.time()
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
P = primes_below(N); Ps = set(P); top = P[-1]
adj = {p: set() for p in P}
for p in P:
    s = math.isqrt(p + 2)
    while s * s <= p + top:
        q = s * s - p
        if q > p and q in Ps: adj[p].add(q); adj[q].add(p)
        s += 1
n = len(P); deg = {p: len(adj[p]) for p in P}
A = [p for p in P if p == 2 or p % 4 == 1]; B = [p for p in P if p % 4 == 3]
leaves = [p for p in P if deg[p] == 1]; assert len(leaves) <= 2
leaf = leaves[0] if leaves else None
print('N = %d: %d primes (A %d, B %d), %d edges, leaves %s; building the flow network' % (N, n, len(A), len(B), sum(deg.values()) // 2, leaves), flush=True)
# flow network: node 0 = source, 1..nA = A, nA+1..nA+nB = B, last = sink
nA, nB = len(A), len(B); SRC = 0; SNK = nA + nB + 1
ia = {p: 1 + i for i, p in enumerate(A)}; ib = {p: 1 + nA + i for i, p in enumerate(B)}
rows = []; cols = []; caps = []
for p in A: rows.append(SRC); cols.append(ia[p]); caps.append(min(2, deg[p]))
edges = [(a, b) for a in A for b in adj[a]]
rng.shuffle(edges)
for a, b in edges: rows.append(ia[a]); cols.append(ib[b]); caps.append(1)
for p in B: rows.append(ib[p]); cols.append(SNK); caps.append(min(2, deg[p]))
M = csr_matrix((np.array(caps, dtype=np.int32), (np.array(rows), np.array(cols))), shape=(SNK + 1, SNK + 1))
res = maximum_flow(M, SRC, SNK, method='dinic')
print('nu_2 = max flow = %d, n - 1 = %d, %.0f s' % (res.flow_value, n - 1, time.time() - T0), flush=True)
if res.flow_value < n - 1:
    print('NO HAMILTONIAN PATH at N = %d: nu_2 < n - 1 by %d' % (N, n - 1 - res.flow_value)); sys.exit(3)
F = res.flow.tocoo() if hasattr(res, 'flow') else res.residual.tocoo()
nxt = {p: [] for p in P}
Ainv = {v: p for p, v in ia.items()}; Binv = {v: p for p, v in ib.items()}
for r, c, f in zip(F.row, F.col, F.data):
    if f > 0 and r in Ainv and c in Binv:
        a, b = Ainv[r], Binv[c]; nxt[a].append(b); nxt[b].append(a)
assert all(len(nxt[p]) <= 2 for p in P) and sum(len(nxt[p]) for p in P) == 2 * res.flow_value
# components
seen = set(); comps = []
for p in P:
    if p in seen: continue
    # walk
    start = p
    if len(nxt[p]) == 2:            # inside a path component? walk back to an endpoint first
        prev, cur = None, p; steps = 0
        while len(nxt[cur]) == 2 and steps <= n:
            nx = [x for x in nxt[cur] if x != prev][0] if prev is not None else nxt[cur][0]
            if nx == p: break
            prev, cur = cur, nx; steps += 1
        if len(nxt[cur]) == 1: start = cur
    order = [start]; seen.add(start); prev = None; cur = start
    while True:
        cands = [x for x in nxt[cur] if x != prev]
        if not cands: break
        nx = cands[0]
        if nx == start: break
        if nx in seen: break
        order.append(nx); seen.add(nx); prev, cur = cur, nx
    comps.append(order)
paths = [c for c in comps if len(c) == 1 or len(nxt[c[0]]) == 1 or len(nxt[c[-1]]) == 1]
cycles = [c for c in comps if c not in paths]
if not paths:                       # a 2-factor: open the largest cycle
    cycles.sort(key=len, reverse=True); paths = [cycles.pop(0)]
assert len(paths) == 1, 'expected one path component, got %d' % len(paths)
path = paths[0]
if leaf is not None and path[-1] == leaf: path.reverse()
print('2-matching: path of %d vertices + %d cycles (lengths %s%s), %.0f s' % (len(path), len(cycles), sorted(len(c) for c in cycles)[:12], '...' if len(cycles) > 12 else '', time.time() - T0), flush=True)
# verification helper
def ok_edge(a, b): return b in adj[a]
assert all(ok_edge(a, b) for a, b in zip(path, path[1:]))
for c in cycles: assert all(ok_edge(c[i], c[(i + 1) % len(c)]) for i in range(len(c)))
# merging
def merge_into_path(C):
    """merge cycle C into the global path by a 2-opt exchange (successor or predecessor variant); True on success."""
    global path
    pos = {v: i for i, v in enumerate(path)}
    L = len(path); m = len(C); cpos = {v: i for i, v in enumerate(C)}
    for c in C:
        i = cpos[c]; cp = C[(i + 1) % m]
        back = [C[(i - k) % m] for k in range(m)]      # c, c-, ..., c+
        fwd = [C[(i + 1 + k) % m] for k in range(m)]   # c+, ..., c
        for p in adj[c]:
            j = pos.get(p)
            if j is None: continue
            pp = path[j + 1] if j + 1 < L else None
            if pp is None or pp in adj[cp]:            # remove (p,p+): p, c, c-, ..., c+, p+
                path = path[:j + 1] + back + path[j + 1:]; return True
            if pp is not None and cp in adj[p] and pp in adj[c]:   # p, c+, ..., c, p+
                path = path[:j + 1] + fwd + path[j + 1:]; return True
            pm = path[j - 1] if j >= 1 else None
            if pm is not None and pm in adj[cp]:       # remove (p-,p): p-, c+, ..., c, p
                path = path[:j] + fwd + path[j:]; return True
    return False
def merge_3opt(C):
    """insert cycle C between two cuts of the path: P1 + C + rev(P2) + P3, with new edges (p,c), (c+,q), (p+,q+)
    (cuts (p,p+) and (q,q+), q after p+), and the mirror P1 + rev(P2) + C + P3 with (p,q), (p+,c), (c+,q+)."""
    global path
    pos = {v: i for i, v in enumerate(path)}; L = len(path); m = len(C)
    for i, c in enumerate(C):
        cp = C[(i + 1) % m]
        seg = [C[(i - k) % m] for k in range(m)]       # from c back round to c+
        for p in adj[c]:
            j = pos.get(p)
            if j is None or j + 1 >= L: continue
            pp = path[j + 1]
            for q in adj[cp]:
                k = pos.get(q)
                if k is None or k <= j + 1: continue
                qp = path[k + 1] if k + 1 < L else None
                if qp is None or qp in adj[pp]:
                    # P1 = path[:j+1], P2 = path[j+1:k+1] reversed, P3 = path[k+1:]
                    path = path[:j + 1] + seg + path[j + 1:k + 1][::-1] + path[k + 1:]; return True
        # mirror: (p,q), (p+,c), (c+,q+) with c adjacent to p+ and c+ adjacent to q+
        for pp in adj[c]:
            j1 = pos.get(pp)
            if j1 is None or j1 == 0: continue
            j = j1 - 1; p = path[j]
            for qp in adj[cp]:
                k1 = pos.get(qp)
                if k1 is None or k1 <= j1: continue
                q = path[k1 - 1]
                if q in adj[p]:
                    # P1 = path[:j+1], rev(P2) = path[j+1:k1][::-1] (from q down to p+), C from c to c+ ... wait need (p+,c): seg2 starts at c
                    seg2 = [C[(i - t) % m] for t in range(m)]   # c, ..., c+
                    path = path[:j + 1] + path[j + 1:k1][::-1] + seg2 + path[k1:]; return True
    return False
def merge_cycles(C, D):
    dpos = {v: i for i, v in enumerate(D)}; m = len(C); k = len(D)
    for i, c in enumerate(C):
        cp = C[(i + 1) % m]
        for d in adj[c]:
            j = dpos.get(d)
            if j is None: continue
            dp = D[(j + 1) % k]
            if dp in adj[cp]: return C[:i + 1] + [D[(j - t) % k] for t in range(k)] + C[i + 1:]
            dm = D[(j - 1) % k]
            if dm in adj[cp]: return C[:i + 1] + [D[(j + t) % k] for t in range(k)] + C[i + 1:]
    return None
def check_path(pth):
    return all(b in adj[a] for a, b in zip(pth, pth[1:])) and len(set(pth)) == len(pth)
cycles.sort(key=len, reverse=True)
rounds = 0
while cycles:
    rounds += 1; rest = []; progress = False
    for C in cycles:
        if merge_into_path(C) or merge_3opt(C): progress = True
        else: rest.append(C)
    cycles = rest
    assert check_path(path)
    if cycles and not progress:
        cycles.sort(key=len, reverse=True); merged = False
        for a in range(len(cycles)):
            for b in range(len(cycles)):
                if a != b:
                    newc = merge_cycles(cycles[a], cycles[b])
                    if newc is not None:
                        cycles = [c for t, c in enumerate(cycles) if t not in (a, b)] + [newc]; merged = True; break
            if merged: break
        if not merged:
            print('STUCK: %d cycles cannot be merged (lengths %s); restart with another seed' % (len(cycles), sorted(len(c) for c in cycles))); sys.exit(4)
    print('  round %d: path %d, cycles left %d, %.0f s' % (rounds, len(path), len(cycles), time.time() - T0), flush=True)
assert len(path) == n and len(set(path)) == n and (leaf is None or path[0] == leaf or path[-1] == leaf)
assert all(math.isqrt(a + b) ** 2 == a + b for a, b in zip(path, path[1:]))
print('HAMILTONIAN PATH FOUND: %d primes, %.0f s; starts %s ends %s' % (n, time.time() - T0, path[:6], path[-6:]))
import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'prime_path_N%d.txt' % N)
open(out, 'w').write(' '.join(map(str, path)) + '\n'); print('written', os.path.relpath(out))
