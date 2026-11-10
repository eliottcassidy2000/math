#!/usr/bin/env python3
# Checks intentionally survive python -O.
import sys
if "--fail" in sys.argv:
    raise RuntimeError("deliberate checker failure")
"""Exact tests: all simple labelled G on <=5 vertices, no graph filters."""
import itertools, json, hashlib, math
from fractions import Fraction
from pathlib import Path

def construct(m, edges):
    shared = {e: k for k, e in enumerate(edges)}
    n = len(edges)
    cycles = []
    tris = []
    forced = {}

    def arc(a, b):
        key = tuple(sorted((a, b)))
        if not (key not in forced or forced[key] == (a, b)):
            raise RuntimeError('live mathematical check failed')
        forced[key] = (a, b)
    for v in range(m):
        es = [k for e, k in shared.items() if v in e]
        ps = list(range(n, n + 5 + len(es) % 2))
        n += len(ps)
        C = es + ps
        cycles.append(C)
        tris.append(ps[:3])
        for a, b in zip(C, C[1:] + C[:1]):
            arc(a, b)
        arc(ps[2], ps[0])
    adj = [[False] * n for _ in range(n)]
    for a, b in itertools.combinations(range(n), 2):
        u, v = forced.get((a, b), (a, b))
        adj[u][v] = True
    return (adj, cycles, tris, forced)

def independence(m, edges):
    return [s for s in range(1 << m) if all((not (s >> a & 1 and s >> b & 1) for a, b in edges))]

def pack_poly(sets):
    co = [0] * (len(sets) + 1)
    for mask in range(1 << len(sets)):
        used = set()
        ok = True
        for i, C in enumerate(sets):
            if mask >> i & 1:
                if used.intersection(C):
                    ok = False
                    break
                used.update(C)
        if ok:
            co[mask.bit_count()] += 1
    return co

def hp(adj):
    n = len(adj)
    dp = [{} for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in dp[mask].items():
            for w in range(n):
                if not mask >> w & 1 and adj[v][w]:
                    d = dp[mask | 1 << w]
                    d[w] = d.get(w, 0) + count
    return sum(dp[-1].values())

def all_cycles(adj):
    n = len(adj)
    cs = []
    for k in range(3, n + 1, 2):
        for vs in itertools.combinations(range(n), k):
            for tail in itertools.permutations(vs[1:]):
                C = (vs[0],) + tail
                if all((adj[a][b] for a, b in zip(C, C[1:] + C[:1]))):
                    cs.append(C)
    return cs

def run():
    count = 0
    hrows = []
    maxn = 0
    for m in range(1, 6):
        pairs = list(itertools.combinations(range(m), 2))
        for mask in range(1 << len(pairs)):
            edges = [e for i, e in enumerate(pairs) if mask >> i & 1]
            adj, C, D, forced = construct(m, edges)
            n = len(adj)
            maxn = max(maxn, n)
            if not n == len(edges) + 5 * m + sum((sum((v in e for e in edges)) % 2 for v in range(m))):
                raise RuntimeError('live mathematical check failed')
            if not n <= len(edges) + 6 * m:
                raise RuntimeError('live mathematical check failed')
            for cyc in C + D:
                if not all((adj[a][b] for a, b in zip(cyc, cyc[1:] + cyc[:1]))):
                    raise RuntimeError('live mathematical check failed')
            for i, j in pairs:
                if not bool(set(C[i]) & set(C[j])) == ((i, j) in edges):
                    raise RuntimeError('live mathematical check failed')
            for i in range(m):
                for j in range(m):
                    if not bool(set(C[i]) & set(D[j])) == (i == j):
                        raise RuntimeError('live mathematical check failed')
            for i, j in pairs:
                if not not set(D[i]) & set(D[j]):
                    raise RuntimeError('live mathematical check failed')
            for i in range(m):
                a = D[i][0] + 3
                local = [(a, b, c) for b, c in itertools.permutations(C[i], 2) if len({a, b, c}) == 3 and adj[a][b] and adj[b][c] and adj[c][a]]
                if not (local and all((set(t) != set(D[i]) for t in local))):
                    raise RuntimeError('live mathematical check failed')
            poly = pack_poly(C + D)
            inds = independence(m, edges)
            formula = [0] * len(poly)
            for S in inds:
                k = S.bit_count()
                for j in range(m - k + 1):
                    formula[k + j] += math.comb(m - k, j)
            if not poly == formula:
                raise RuntimeError('live mathematical check failed')
            for x in (1, 2, 3, 4):
                left = sum((c * x ** k for k, c in enumerate(poly)))
                right = sum((x ** s.bit_count() * (1 + x) ** (m - s.bit_count()) for s in inds))
                if not left == right:
                    raise RuntimeError('live mathematical check failed')
            if m <= 2:
                H = hp(adj)
                bound = sum((2 ** s.bit_count() * 3 ** (m - s.bit_count()) for s in inds))
                if not H >= bound + 2 * m:
                    raise RuntimeError('live mathematical check failed')
                hrows.append(dict(m=m, edges=edges, n=n, H=H, corona_bound=bound, strengthened_bound=bound + 2 * m))
                if m == 1:
                    cs = all_cycles(adj)
                    p = pack_poly(cs)
                    if not H == sum((a * 2 ** k for k, a in enumerate(p))):
                        raise RuntimeError('live mathematical check failed')
                    brute = sum((all((adj[a][b] for a, b in zip(P, P[1:]))) for P in itertools.permutations(range(n))))
                    if not brute == H:
                        raise RuntimeError('live mathematical check failed')
            count += 1
    adj, C, D, forced = construct(1, [])
    free = [e for e in itertools.combinations(range(5), 2) if e not in forced]
    completion_values = []
    for mask in range(1 << len(free)):
        A = [row[:] for row in adj]
        for j, (a, b) in enumerate(free):
            A[a][b] = bool(mask >> j & 1)
            A[b][a] = not A[a][b]
        cs = all_cycles(A)
        poly = pack_poly(cs)
        if not (poly[0] == 1 and poly[1] >= 3):
            raise RuntimeError('live mathematical check failed')
        H = hp(A)
        if not (H == sum((a * 2 ** k for k, a in enumerate(poly))) and H >= 7):
            raise RuntimeError('live mathematical check failed')
        completion_values.append(H)
    hostile = sum((2 ** s.bit_count() for s in independence(3, list(itertools.combinations(range(3), 2)))))
    if not hostile == 7:
        raise RuntimeError('live mathematical check failed')
    return dict(status='PASS', universe='all labelled simple graphs with 1..5 vertices', graphs=count, max_constructed_order=maxn, hamiltonian_controls=hrows, hostile_K3_value=hostile, one_vertex_completions=len(completion_values), one_vertex_completion_H_values=sorted(set(completion_values)))
if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
