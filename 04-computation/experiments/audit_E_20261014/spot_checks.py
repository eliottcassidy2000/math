#!/usr/bin/env python3
"""Audit E (2026-10-14): single-tournament spot checks by plain DFS (gentourng ascii: upper triangle row by row, 1 iff i->j).
(1) HYP-9279 refuted variant (H): the n = 9 tournament 110111111101001111111111111111111101 (details as claimed);
(2) the claimed n = 10 minimiser 111000111101111111111111111111110001111111101 (delta = 62/157)."""
from fractions import Fraction as F

def analyse(s, n):
    A = [[0] * n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1': A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    paths = []
    def rec(p, used):
        if len(p) == n: paths.append(tuple(p)); return
        for v in range(n):
            if not used >> v & 1 and A[p[-1]][v]:
                p.append(v); rec(p, used | 1 << v); p.pop()
    for s0 in range(n): rec([s0], 1 << s0)
    H = len(paths); bef = [[0] * n for _ in range(n)]
    for p in paths:
        for i in range(n):
            for j in range(i + 1, n): bef[p[i]][p[j]] += 1
    R = [[A[i][j] or i == j for j in range(n)] for i in range(n)]
    for kk in range(n):
        for i in range(n):
            for j in range(n):
                if R[i][kk] and R[kk][j]: R[i][j] = True
    bal = lambda x, y: min(F(bef[x][y], H), 1 - F(bef[x][y], H))
    delta = max(bal(x, y) for x in range(n) for y in range(n) if x != y)
    vb = min(max(bal(v, w) for w in range(n) if w != v) for v in range(n))
    cyc = [p for p in paths if p[0] == 0 and A[p[-1]][0]]
    nobal = [c for c in cyc if not any(F(1, 3) <= F(bef[c[i]][c[(i + 1) % n]], H) <= F(2, 3) for i in range(n))]
    print(f"n={n} {s}: strong={all(all(r) for r in R)} H={H} delta={delta} min vertex balance={vb} scores={[sum(r) for r in A]} "
          f"backward arcs={[(i, j) for i in range(n) for j in range(n) if A[i][j] and i > j]} Hamiltonian cycles={len(cyc)} "
          f"without a balanced arc={len(nobal)}")
    for c in nobal:
        print("   cycle", c, "P(tail before head):", [str(F(bef[c[i]][c[(i + 1) % n]], H)) for i in range(n)])

analyse('110111111101001111111111111111111101', 9)
analyse('111000111101111111111111111111110001111111101', 10)
