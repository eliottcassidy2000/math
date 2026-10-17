#!/usr/bin/env python3
"""Exact maximum acyclic subgraph MAS(P_p) of Paley tournaments, with an independent brute-force cross-check (S40).

Subset DP with vertex 0 first (P_p is vertex-transitive) and, for the largest case, a second reduction: the stabiliser
of 0 (multiplication by squares) is transitive on the residues and on the non-residues, so some optimal order begins
0, 1 or 0, n0 (n0 the least non-residue); two DPs over subsets of the remaining p - 2 vertices.
The brute force enumerates all (p-1)! orders with 0 first (p <= 11) by Heap's algorithm.
Usage: python paley_mas_exact_20261017.py P [--brute]
Session opus-2026-10-17-S40.
"""
import sys, os, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import kronecker_table, class_number, mas_dp


@njit(cache=True)
def mas_dp_two(adj, n, second):
    """Max forward arcs over orders starting 0, second; DP over subsets of the other n - 2 vertices."""
    rest = np.empty(n - 2, dtype=np.int64)
    k = 0
    for v in range(1, n):
        if v != second:
            rest[k] = v
            k += 1
    m = n - 2
    pos = np.full(n, -1, dtype=np.int64)
    for i in range(m):
        pos[rest[i]] = i
    inmask = np.zeros(m, dtype=np.int64)
    head = np.zeros(m, dtype=np.int64)          # arcs from {0, second} into v
    for i in range(m):
        v = rest[i]
        msk = 0
        for j in range(m):
            u = rest[j]
            if (adj[u] >> v) & 1:
                msk |= 1 << j
        inmask[i] = msk
        head[i] = ((adj[0] >> v) & 1) + ((adj[second] >> v) & 1)
    base0 = (adj[0] >> second) & 1
    size = 1 << m
    dp = np.full(size, -1, dtype=np.int16)
    dp[0] = base0
    for S in range(size):
        b = dp[S]
        if b < 0:
            continue
        for i in range(m):
            bit = 1 << i
            if S & bit:
                continue
            x = S & inmask[i]
            c = 0
            while x:
                x &= x - 1
                c += 1
            val = b + c + head[i]
            T = S | bit
            if val > dp[T]:
                dp[T] = val
    return dp[size - 1]


@njit(cache=True)
def brute(adj, n):
    """All orders with 0 first (Heap's algorithm on 1..n-1)."""
    a = np.arange(1, n)
    m = n - 1
    c = np.zeros(m, dtype=np.int64)

    def score(a):
        s = 0
        seen = 1                      # vertex 0 placed first
        for i in range(m):
            v = a[i]
            x = seen & ~(adj[v])       # placed vertices u with u -> v  <=> not (v -> u)
            x &= ~(1 << v)
            cnt = 0
            while x:
                x &= x - 1
                cnt += 1
            s += cnt
            seen |= 1 << v
        return s

    best = score(a)
    i = 0
    while i < m:
        if c[i] < i:
            if i % 2 == 0:
                a[0], a[i] = a[i], a[0]
            else:
                a[c[i]], a[i] = a[i], a[c[i]]
            s = score(a)
            if s > best:
                best = s
            c[i] += 1
            i = 0
        else:
            c[i] = 0
            i += 1
    return best


def tournament(p):
    chi = kronecker_table(-p)
    adj = np.zeros(p, dtype=np.int64)
    for u in range(p):
        m = 0
        for v in range(p):
            if v != u and chi[(v - u) % p] == 1:
                m |= 1 << v
        adj[u] = m
    return chi, adj


def main():
    p = int(sys.argv[1])
    chi, adj = tournament(p)
    arcs = p * (p - 1) // 2
    h = class_number(-p)
    natural = (arcs + (2 * p * h) // (6 if p == 3 else 2)) // 2
    res = dict(p=p, arcs=arcs, h=h, natural_forward=natural)
    if '--brute' in sys.argv:
        res['brute'] = int(brute(adj, p))
        res['dp_vertex0'] = int(mas_dp(adj, p))
        print('P_%d: brute force %d, vertex-0 DP %d' % (p, res['brute'], res['dp_vertex0']), flush=True)
    n0 = next(a for a in range(2, p) if chi[a] == -1)
    v1 = int(mas_dp_two(adj, p, 1)); v2 = int(mas_dp_two(adj, p, n0))
    res.update(MAS=max(v1, v2), MAS_second_residue=v1, MAS_second_nonresidue=v2)
    print('P_%d: arcs %d, h(-p) = %d, natural forward %d, MAS %d (second vertex residue: %d, non-residue %d: %d)'
          % (p, arcs, h, natural, max(v1, v2), v1, n0, v2), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(res, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_mas_exact_p%d_20261017.json' % p), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
