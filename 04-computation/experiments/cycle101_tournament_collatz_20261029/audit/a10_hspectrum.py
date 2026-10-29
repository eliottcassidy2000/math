#!/usr/bin/env python3
"""Spot checks of results-note section 1: H = 169 absent at n = 7, present at n = 8 (all classes from gentourng);
the skew characteristic polynomial of R7 (session101 controls)."""
import subprocess, sys
from sympy import Matrix, symbols, factor
def ham(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for m in range(1, 1 << n):
        for v in range(n):
            c = dp[m][v]
            if c:
                for u in range(n):
                    if not (m >> u) & 1 and adj[v][u]: dp[m | (1 << u)][u] += c
    return sum(dp[-1])
def parse(line, n):
    # gentourng ascii: upper triangle row by row, '1' means i -> j
    bits = line.strip()
    adj = [[0] * n for _ in range(n)]
    p = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[p] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            p += 1
    return adj
for n in (7, 8):
    out = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True).stdout.split()
    Hs = sorted(set(ham(parse(l, n), n) for l in out))
    print(f"n = {n}: {len(out)} classes; 169 among H values: {169 in Hs}; H range {Hs[0]}..{Hs[-1]}")
x = symbols('x')
A = Matrix(7, 7, lambda i, j: 1 if (j - i) % 7 in (1, 2, 3) else 0)
print("det(xI - (A - A^T)) for R7 =", factor((x * Matrix.eye(7) - (A - A.T)).det()))
