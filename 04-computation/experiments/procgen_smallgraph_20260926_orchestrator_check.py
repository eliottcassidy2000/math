#!/usr/bin/env python3
"""Orchestrator audit of lane `smallgraph`, written from the note's statements;
the lane's scripts were not read.

  1. Theorem 1: the leaf table N(x) and "min degree >= 2 iff n >= 31" for the
     square-sum graph Q_n, n <= 300.
  2. Theorem 2: Q_n has no Hamiltonian path for 2 <= n <= 14; Q_15 has exactly
     one (up to reversal), from 8 to 9.
  3. Theorem 5 (zigzag threshold): for targets j^2 + k(4-k) (j >= 1) the first n
     with a Hamiltonian path is 4k-1, the path is unique, ends 2k and 2k+1
     (k = 2..9, exhaustive path counting).
  4. Proposition F4 (Sundaram): complete i x j grids with n fences <=> 2n+1 odd
     composite (n <= 2000).
"""
import math


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def graph(n, targets):
    tset = set(targets)
    adj = {x: [y for y in range(1, n + 1) if y != x and (x + y) in tset] for x in range(1, n + 1)}
    return adj


def count_ham_paths(n, adj, limit=10 ** 7):
    """number of Hamiltonian paths (each path counted once per direction) and one example."""
    order = sorted(range(1, n + 1), key=lambda v: len(adj[v]))
    total = 0
    ends = set()
    example = None
    visited = [False] * (n + 1)
    path = []

    def dfs(v, depth):
        nonlocal total, example
        if depth == n:
            total += 1
            ends.add((path[0], path[-1]))
            if example is None:
                example = list(path)
            return
        for w in adj[v]:
            if not visited[w]:
                visited[w] = True
                path.append(w)
                dfs(w, depth + 1)
                path.pop()
                visited[w] = False

    for s in range(1, n + 1):
        visited[s] = True
        path.append(s)
        dfs(s, 1)
        path.pop()
        visited[s] = False
    return total, ends, example


# 1. leaf table
squares = [j * j for j in range(1, 60)]
Ntab = {1: 7, 2: 13, 3: 5, 4: 11, 5: 10, 6: 9, 7: 8, 8: 16, 9: 15, 10: 14, 11: 13, 12: 12, 16: 19, 17: 18, 18: 30}
for n in range(1, 301):
    sq = set(j * j for j in range(1, int(math.isqrt(2 * n)) + 2))
    for x in range(1, n + 1):
        deg = sum(1 for y in range(1, n + 1) if y != x and (x + y) in sq)
        predicted_leaf = (x in Ntab and x <= n <= Ntab[x])
        assert (deg <= 1) == predicted_leaf, (n, x, deg)
mindeg_ok = []
for n in range(2, 301):
    sq = set(j * j for j in range(1, int(math.isqrt(2 * n)) + 2))
    md = min(sum(1 for y in range(1, n + 1) if y != x and (x + y) in sq) for x in range(1, n + 1))
    mindeg_ok.append((n, md >= 2))
check(all(ok == (n >= 31) for n, ok in mindeg_ok), "Theorem 1: the leaf table is exact for n <= 300, and min degree >= 2 iff n >= 31")

# 2. Theorem 2
sq = [j * j for j in range(1, 10)]
for n in range(2, 15):
    t, e, ex = count_ham_paths(n, graph(n, sq))
    assert t == 0, (n, t)
t, e, ex = count_ham_paths(15, graph(15, sq))
check(t == 2 and e == {(8, 9), (9, 8)}, f"Theorem 2: no square-sum path for 2 <= n <= 14; exactly one at n = 15 (up to reversal), ends 8 and 9: {ex}")

# 3. Theorem 5
for k in range(2, 10):
    c = k * (4 - k)
    targets = [j * j + c for j in range(1, 60) if j * j + c >= 3]
    first = None
    for n in range(2, 4 * k):
        t, e, ex = count_ham_paths(n, graph(n, targets))
        if t > 0:
            first = (n, t, e)
            break
    assert first is not None and first[0] == 4 * k - 1 and first[1] == 2 and first[2] == {(2 * k, 2 * k + 1), (2 * k + 1, 2 * k)}, (k, first)
check(True, "Theorem 5: for targets j^2 + k(4-k), the first Hamiltonian path appears at n = 4k-1, unique, with ends 2k and 2k+1 (k = 2..9)")

# 4. Sundaram
grid_n = set()
for i in range(1, 1001):
    for j in range(i, 1001):
        grid_n.add(2 * i * j + i + j)
for n in range(1, 2001):
    m = 2 * n + 1
    composite = any(m % d == 0 for d in range(3, int(math.isqrt(m)) + 1, 2))
    assert (n in grid_n) == composite, n
check(True, "Proposition F4: a complete grid with n unit fences exists iff 2n+1 is an odd composite (n <= 2000; Sundaram's sieve)")
