#!/usr/bin/env python3
"""The square-sum problem: for which n can 1..n be ordered so that every two consecutive numbers sum to a square?
(Hamiltonian paths in the square-sum graph.)  Exhaustive DFS with a degree heuristic, n <= NMAX.
Known (OEIS A090461 and the folklore 'square-sum problem'): paths exist for n = 15, 16, 17, 23 and every n >= 25 checked;
this script reproduces the exceptional set below NMAX.  Also: Hamiltonian CYCLES (closing sum also a square).
Session opus-2026-10-14-S36.
"""
import sys
from math import isqrt

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
sys.setrecursionlimit(10000)


def graph(n):
    sq = {k * k for k in range(2, isqrt(2 * n) + 2)}
    return {i: [j for j in range(1, n + 1) if j != i and (i + j) in sq] for i in range(1, n + 1)}


def ham_path(n, cycle=False):
    G = graph(n)
    if any(len(G[i]) == 0 for i in G):
        return None
    order = sorted(G, key=lambda v: len(G[v]))
    used = [False] * (n + 1)
    path = []

    def dfs(v):
        path.append(v); used[v] = True
        if len(path) == n:
            if not cycle or path[0] in G[v]:
                return True
        else:
            for w in sorted(G[v], key=lambda u: sum(1 for x in G[u] if not used[x])):
                if not used[w] and dfs(w):
                    return True
        path.pop(); used[v] = False
        return False

    starts = order if not cycle else [order[0]]
    for s in starts:
        if dfs(s):
            return list(path)
    return None


if __name__ == '__main__':
    has_path = [n for n in range(1, NMAX + 1) if ham_path(n) is not None]
    print(f'n <= {NMAX} with a square-sum ordering of 1..n: {has_path}')
    print('exceptions (no ordering):', [n for n in range(2, NMAX + 1) if n not in has_path])
    print('example n = 15:', ham_path(15))
    cyc = [n for n in range(2, min(NMAX, 40) + 1) if ham_path(n, cycle=True) is not None]
    print(f'n <= {min(NMAX, 40)} with a square-sum CYCLE: {cyc}')
