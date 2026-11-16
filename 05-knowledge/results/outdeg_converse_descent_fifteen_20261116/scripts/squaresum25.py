#!/usr/bin/env python3
"""squaresum25.py: the square-sum graph on 1..25 (i ~ j iff i + j is a perfect square) in rows of five; degrees; which rows
contain edges; the row-pair adjacency; the number of Hamiltonian paths on 1..15 (DP) and the endpoints; for 1..25 the number
of Hamiltonian paths by bitmask DP in a compact form (2^25 states, counted with endpoints)."""
import math
def adj(n):
    sq = {k*k for k in range(2, 10)}
    return [[1 if i != j and (i + j) in sq else 0 for j in range(1, n+1)] for i in range(1, n+1)]
n = 25; A = adj(n)
rows = [list(range(5*r, 5*r + 5)) for r in range(5)]
print("degrees 1..25:", [sum(A[i]) for i in range(n)])
print("rows of five with an internal edge:", [r+1 for r, R in enumerate(rows) if any(A[i][j] for i in R for j in R if i < j)])
pairs = [[sum(A[i][j] for i in rows[a] for j in rows[b]) for b in range(5)] for a in range(5)]
print("edges between rows (matrix):", pairs)
print("row pairs with no edge:", [(a+1, b+1) for a in range(5) for b in range(a+1, 5) if pairs[a][b] == 0])
# Hamiltonian paths on 1..15 and endpoints
def hampaths(A, n):
    N = 1 << n; dp = [[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v] = 1
    for m in range(1, N):
        for l in range(n):
            a = dp[m][l]
            if not a: continue
            for v in range(n):
                if not (m >> v) & 1 and A[l][v]: dp[m|(1<<v)][v] += a
    ends = [dp[N-1][l] for l in range(n)]
    return sum(ends) // 2, [l+1 for l in range(n) if ends[l]]
for m in (15, 16, 17, 23):
    c, e = hampaths(adj(m), m); print(f"square-sum graph on 1..{m}: Hamiltonian paths {c}, endpoints {e}")
