#!/usr/bin/env python3
"""Explicit one-defect witnesses (THM-5071 (4)): a tournament with ONE arc deleted (one mirrored 00 pair, no loop,
no reciprocal pair) with H = 7 (order 5) and H = 21 (order 6); a semicomplete digraph with exactly one reciprocal
pair and H = 21 (order 6); and the sink back-arc construction with H = h for every h (h - 1 in binary)."""
import itertools

def ham(n, A):
    N = 1 << n
    dp = [[0] * n for _ in range(N)]
    for v in range(n): dp[1 << v][v] = 1
    for s in range(1, N):
        for v in range(n):
            c = dp[s][v]
            if c:
                for w in range(n):
                    if not s >> w & 1 and A[v][w]: dp[s | 1 << w][w] += c
    return sum(dp[N - 1])

def first(n, target, mode):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: A[i][j] = 1
            else: A[j][i] = 1
        for (i, j) in pairs:
            a, b = A[i][j], A[j][i]
            if mode == "delete": A[i][j] = A[j][i] = 0
            else: A[i][j] = A[j][i] = 1
            if ham(n, A) == target:
                return ["".join(map(str, r)) for r in A], (i, j)
            A[i][j], A[j][i] = a, b

for (n, t, mode) in ((5, 7, "delete"), (6, 21, "delete"), (6, 21, "add")):
    rows, pair = first(n, t, mode)
    print(f"order {n}, H = {t}, one {'missing' if mode == 'delete' else 'reciprocal'} pair {pair}: rows {' '.join(rows)}")

def sink_construction(h):
    """TT_n plus back arcs from the sink n-1 to {x : bit (n-2-x) of h-1 is 1}: H = 1 + sum 2^(n-2-x)."""
    m = h - 1
    n = max(2, m.bit_length() + 1)
    A = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
    for x in range(n - 1):
        if m >> (n - 2 - x) & 1:
            A[n - 1][x] = 1
    return n, A
for h in range(1, 300):
    n, A = sink_construction(h)
    assert ham(n, A) == h, h
print("sink back-arc construction: H = h for every h in 1..299 (checked)")
