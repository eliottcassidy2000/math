#!/usr/bin/env python3
"""S36 referee audit, part H: square-sum orderings of 1..n (Hamiltonian paths/cycles of the square-sum graph).

Independent method: exact bitmask dynamic programming over vertex subsets for Hamiltonian PATHS when n <= 26 (so the
exceptional set below 27 is decided without search heuristics), explicit verified witnesses for 27 <= n <= 40;
Hamiltonian CYCLES: minimum-degree obstruction or exhaustive DFS with the start fixed at vertex 1, n <= 40.
"""
import sys, time
from math import isqrt
import numpy as np
import numba as nb


def sqgraph(n):
    sqs = {k * k for k in range(2, isqrt(2 * n) + 2)}
    return [[j for j in range(1, n + 1) if j != i and i + j in sqs] for i in range(1, n + 1)]


@nb.njit(cache=True)
def ham_path_dp(n, nbr):
    full = (1 << n) - 1
    reach = np.zeros(1 << n, dtype=np.int32)
    for v in range(n):
        reach[1 << v] = 1 << v
    for mask in range(1, full + 1):
        e = reach[mask]
        if e == 0:
            continue
        for v in range(n):
            if (e >> v) & 1:
                cand = nbr[v] & ~mask
                while cand:
                    lowbit = cand & (-cand)
                    u = 0
                    while (lowbit >> u) != 1:
                        u += 1
                    reach[mask | lowbit] |= 1 << u
                    cand ^= lowbit
    return reach[full] != 0


def find_path(n, cycle=False, limit=5_000_000):
    G = sqgraph(n)
    adj = [set(u - 1 for u in x) for x in G]
    used = [False] * n
    path = []
    steps = [0]

    def dfs(v):
        steps[0] += 1
        if steps[0] > limit:
            raise TimeoutError
        path.append(v); used[v] = True
        if len(path) == n:
            if not cycle or path[0] in adj[v]:
                return True
        else:
            for w in sorted((u - 1 for u in G[v]), key=lambda u: sum(1 for x in G[u] if not used[x - 1])):
                if not used[w] and dfs(w):
                    return True
        path.pop(); used[v] = False
        return False

    starts = [0] if cycle else sorted(range(n), key=lambda v: len(G[v]))
    for s0 in starts:
        if dfs(s0):
            return [x + 1 for x in path]
    return None


def check(seq, n, cycle):
    if sorted(seq) != list(range(1, n + 1)):
        return False
    pairs = list(zip(seq, seq[1:])) + ([(seq[-1], seq[0])] if cycle else [])
    return all(isqrt(a + b) ** 2 == a + b for a, b in pairs)


if __name__ == '__main__':
    sys.setrecursionlimit(20000)
    t0 = time.time()
    paths_dp = []
    for n in range(2, 27):
        G = sqgraph(n)
        nbr = np.array([sum(1 << (u - 1) for u in G[i]) for i in range(n)], dtype=np.int64)
        if ham_path_dp(n, nbr):
            paths_dp.append(n)
    print(f'Hamiltonian square-sum PATH exists (exact subset DP), 2 <= n <= 26: {paths_dp}  [{time.time() - t0:.1f} s]')
    wit = {}
    for n in range(27, 41):
        p = find_path(n)
        wit[n] = p is not None and check(p, n, False)
    print('verified explicit witnesses for 27 <= n <= 40:', all(wit.values()))
    cyc = []
    for n in range(2, 41):
        G = sqgraph(n)
        if min(len(x) for x in G) < 2:
            continue
        try:
            p = find_path(n, cycle=True, limit=20_000_000)
        except TimeoutError:
            cyc.append((n, 'timeout'))
            continue
        if p is not None:
            assert check(p, n, True)
            cyc.append(n)
    print('square-sum CYCLE exists for n <= 40:', cyc, f'  [{time.time() - t0:.1f} s]')
    print('n <= 40 with min degree < 2 (no cycle possible):', [n for n in range(2, 41) if min(len(x) for x in sqgraph(n)) < 2])
