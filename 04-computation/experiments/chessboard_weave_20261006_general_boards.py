#!/usr/bin/env python3
"""General 2m x 2m boards: ring sizes, slider mobility by ring, knight ring incidence.
Checks the closed forms used in the chessboard-weave note (2026-10-06)."""
from collections import Counter
from itertools import product

def ring(i, j, n):
    c = (n - 1) / 2
    return int(max(abs(i - c), abs(j - c)) - 0.5)

def mobility(i, j, n, dirs):
    tot = 0
    for di, dj in dirs:
        a, b = i + di, j + dj
        while 0 <= a < n and 0 <= b < n:
            tot += 1; a += di; b += dj
    return tot

B4 = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
R4 = [(1, 0), (-1, 0), (0, 1), (0, -1)]
KN = [(a * x, b * y) for (a, b) in ((1, 2), (2, 1)) for x in (1, -1) for y in (1, -1)]

for n in (2, 4, 6, 8, 10, 12, 14, 16):
    m = n // 2
    sq = [(i, j) for i in range(n) for j in range(n)]
    rs = Counter(ring(i, j, n) for i, j in sq)
    ok_sizes = all(rs[k] == 4 * (2 * k + 1) for k in range(m))
    ok_bishop = all(mobility(i, j, n, B4) == 2 * n - 3 - 2 * ring(i, j, n) for i, j in sq)
    ok_rook = all(mobility(i, j, n, R4) == 2 * (n - 1) for i, j in sq)
    # knight edges by ring pair
    E = set()
    for i, j in sq:
        for di, dj in KN:
            a, b = i + di, j + dj
            if 0 <= a < n and 0 <= b < n:
                E.add(frozenset(((i, j), (a, b))))
    within = Counter(ring(*tuple(e)[0], n) for e in E if ring(*tuple(e)[0], n) == ring(*tuple(e)[1], n))
    dr = Counter(abs(ring(*tuple(e)[0], n) - ring(*tuple(e)[1], n)) for e in E)
    # is each ring's within-ring knight subgraph a matching?
    deg = Counter()
    for e in E:
        u, v = tuple(e)
        if ring(*u, n) == ring(*v, n):
            deg[u] += 1; deg[v] += 1
    matching = all(d <= 1 for d in deg.values())
    print(f"n={n:2d}: ring sizes 4(2k+1) {ok_sizes}; bishop mobility = 2n-3-2*ring {ok_bishop}; rook const {ok_rook}; "
          f"knight edges {len(E)} (formula 4(n-1)(n-2)={4*(n-1)*(n-2)}); within-ring by ring {dict(sorted(within.items()))}; "
          f"|dring| {dict(sorted(dr.items()))}; within-ring edges form a matching: {matching} covering {len(deg)} squares")
