#!/usr/bin/env python3
"""prime_path_2matching.py (klein-2026-11-01) -- the 2-matching obstruction to prime square-sum Hamiltonian paths.
A Hamiltonian path of G_N has n - 1 edges and every vertex of degree <= 2, so nu_2(G_N) >= n - 1, where nu_2 is the
maximum number of edges of a spanning subgraph with all degrees <= 2.  For a bipartite graph nu_2 is a maximum flow
(source -> A capacity 2, graph edges capacity 1, B -> sink capacity 2), computed exactly here with scipy's Dinic.
When nu_2 < n - 1 the minimum cut gives a Hall certificate: the set S of A-vertices reachable from the source in the
residual network has 2|S| - cap(N(S)) = n - 1 - nu_2 > 0 ... more precisely nu_2 = 2|A| - max_S (2|S| - cap(N(S))),
cap(b) = min(deg b, 2), and the certificate is re-verified from the adjacency lists, independently of the flow.
usage: python3 prime_path_2matching.py N            -- one bound, with the certificate and its residue classes
       python3 prime_path_2matching.py window LO HI  -- every prime bound N in [LO, HI] passing the count/leaf/isolated
                                                        conditions of THM-4860 (d),(d'); reports nu_2 - (n - 1)"""
import math, sys, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow, breadth_first_order
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
def graph(P):
    Ps = set(P); top = P[-1]; adj = {p: set() for p in P}
    for p in P:
        s = math.isqrt(p + 2)
        while s * s <= p + top:
            q = s * s - p
            if q > p and q in Ps: adj[p].add(q); adj[q].add(p)
            s += 1
    return adj
def nu2(P, adj, want_cert=False):
    A = [p for p in P if p == 2 or p % 4 == 1]; B = [p for p in P if p % 4 == 3]
    nA, nB = len(A), len(B); SRC = 0; SNK = nA + nB + 1
    ia = {p: 1 + i for i, p in enumerate(A)}; ib = {p: 1 + nA + i for i, p in enumerate(B)}
    rows = []; cols = []; caps = []
    for p in A: rows.append(SRC); cols.append(ia[p]); caps.append(min(2, len(adj[p])))
    for a in A:
        for b in adj[a]: rows.append(ia[a]); cols.append(ib[b]); caps.append(1)
    for p in B: rows.append(ib[p]); cols.append(SNK); caps.append(min(2, len(adj[p])))
    M = csr_matrix((np.array(caps, dtype=np.int32), (np.array(rows), np.array(cols))), shape=(SNK + 1, SNK + 1))
    res = maximum_flow(M, SRC, SNK, method='dinic')
    val = int(res.flow_value)
    if not want_cert: return val, None
    # residual reachability from the source
    F = res.flow if hasattr(res, 'flow') else res.residual
    R = (M - F)           # residual capacity on forward arcs; backward arcs = positive flow
    Rt = F.T.tocsr(); Rt.data = np.maximum(Rt.data, 0)
    Res = (R + Rt).tocsr(); Res.eliminate_zeros()
    reach = set(breadth_first_order(Res, SRC, directed=True, return_predecessors=False).tolist())
    Ainv = {v: p for p, v in ia.items()}
    S = [Ainv[v] for v in reach if v in Ainv]
    NS = set(); 
    for a in S: NS |= adj[a]
    capNS = sum(min(2, len(adj[b])) for b in NS)
    demand = sum(min(2, len(adj[a])) for a in S)
    return val, (S, NS, demand, capNS)
if sys.argv[1] == 'window':
    LO, HI = int(sys.argv[2]), int(sys.argv[3])
    P = primes_below(HI + 1000); import bisect
    # candidate bounds: recompute the count conditions (cheap) and the leaf/isolated conditions (incremental degrees)
    Ps = set(P); deg = {}; iso = 0; leaves = 0; c = {1: 0, 5: 0, 7: 0, 11: 0}; cands = []
    def bump(v):
        global iso, leaves
        d = deg[v]
        if d == 0: iso -= 1
        elif d == 1: leaves -= 1
        deg[v] = d + 1
        if d + 1 == 1: leaves += 1
    for i, p in enumerate(P[:-1]):
        deg[p] = 0; iso += 1
        if p > 3: c[p % 12] += 1
        s = math.isqrt(p) + 1
        while s * s < 2 * p:
            q = s * s - p
            if q in Ps and q < p: bump(q); bump(p)
            s += 1
        N = P[i + 1]
        if N < LO or N > HI: continue
        Aa = c[1] + c[5] + 1; Bb = c[7] + c[11] + 1
        if iso == 0 and leaves <= 2 and abs(Aa - Bb) <= 1 and c[1] <= c[11] + 2 and c[7] <= c[5] + 2:
            cands.append((N, i + 1))
    print('window [%d, %d]: %d candidate bounds' % (LO, HI, len(cands)), flush=True)
    fails = 0; passes = []; worst = 0; T0 = time.time()
    for k, (N, cnt) in enumerate(cands):
        Pn = P[:cnt]; adj = graph(Pn); val, _ = nu2(Pn, adj)
        gap = val - (cnt - 1)
        if gap < 0: fails += 1; worst = min(worst, gap)
        else: passes.append(N)
        if k % 10 == 0 or gap >= 0:
            print('  N = %d (%d primes): nu_2 - (n-1) = %d   [%d done, %d excluded, %.0f s]' % (N, cnt, gap, k + 1, fails, time.time() - T0), flush=True)
    print('RESULT: %d candidate bounds, %d excluded by nu_2 < n - 1 (most negative gap %d), %d NOT excluded: %s' % (len(cands), fails, worst, len(passes), passes[:20]))
else:
    N = int(sys.argv[1]); P = primes_below(N); adj = graph(P); n = len(P); T0 = time.time()
    val, cert = nu2(P, adj, want_cert=True)
    print('N = %d: n = %d, nu_2 = %d, n - 1 = %d, gap %d (%.0f s)' % (N, n, val, n - 1, val - (n - 1), time.time() - T0))
    S, NS, demand, capNS = cert
    cls = lambda p: 2 if p == 2 else (3 if p == 3 else p % 12)
    hist = lambda L: {k: sum(1 for x in L if cls(x) == k) for k in (1, 5, 7, 11, 2, 3) if any(cls(x) == k for x in L)}
    print('Hall certificate: |S| = %d A-vertices with demand %d; |N(S)| = %d with capacity %d; deficiency %d (= 2|A| - nu_2 = %d)' % (
        len(S), demand, len(NS), capNS, demand - capNS, 2 * sum(1 for p in P if p == 2 or p % 4 == 1) - val))
    print('  classes mod 12 of S:', hist(S), '; of N(S):', hist(NS))
    print('  degrees in S: min %d, median %d, max %d; N(S) vertices of degree 1: %d' % (min(len(adj[a]) for a in S), sorted(len(adj[a]) for a in S)[len(S) // 2], max(len(adj[a]) for a in S), sum(1 for b in NS if len(adj[b]) == 1)))
    print('  S (sorted, first 40):', sorted(S)[:40])
