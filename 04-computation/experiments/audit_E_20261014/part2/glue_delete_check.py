#!/usr/bin/env python3
"""Audit E part 2 (2026-10-15): HYP-9279 'Deleting suitable vertices from it keeps delta = 592/1521 at n = 12':
delete each vertex of the n = 13 tournament M6b + M6a glued at a cut vertex (delta 592/1521, H = 13689) and compute delta
of every strong 12-vertex result by DFS."""
from fractions import Fraction as F
def d6dec(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63; bits += [(v >> (5 - t)) & 1 for t in range(6)]
    return [[bits[i * n + j] for j in range(n)] for i in range(n)]
def delta(A):
    n = len(A); P = []
    def rec(p, used):
        if len(p) == n: P.append(tuple(p)); return
        for v in range(n):
            if not used >> v & 1 and A[p[-1]][v]:
                p.append(v); rec(p, used | 1 << v); p.pop()
    for s in range(n): rec([s], 1 << s)
    H = len(P); bef = [[0] * n for _ in range(n)]
    for p in P:
        for i in range(n):
            for j in range(i + 1, n): bef[p[i]][p[j]] += 1
    return H, max(min(F(bef[x][y], H), 1 - F(bef[x][y], H)) for x in range(n) for y in range(x + 1, n))
def strong(A):
    n = len(A); R = [[bool(A[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]: R[i][j] = True
    return all(all(r) for r in R)
A = d6dec('&LPJfNx^{^}f}`~o^oE?@{?[@C?`?Q?')
H, d = delta(A); print(f"n=13 tournament: H={H} delta={d}")
for x in range(13):
    B = [[A[i][j] for j in range(13) if j != x] for i in range(13) if i != x]
    if not strong(B): print(f"  delete {x}: not strong"); continue
    H, d = delta(B); print(f"  delete {x}: H={H} delta={d} = {float(d):.6f}{'  <-- 592/1521' if d == F(592, 1521) else ''}")
