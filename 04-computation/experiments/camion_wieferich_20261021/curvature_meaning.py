# Combinatorial meaning of the curvature statistic X (THM-4750 (3)):  X(m-2) = m! [x^(m-2)] T''/T^2 ... i.e. X(n) = n![x^n] T''(x)/T(x)^2,
# claimed to equal the number of labelled tournaments on [m] (m = n+2) that become strongly connected when the arc at the
# fixed pair {0,1} is made two-way.  Brute force for m <= 7.  mac-mini-2026-10-21-fanopsp.
from fractions import Fraction
from math import comb, factorial
from itertools import product
N = 8
T = [Fraction(2 ** comb(n, 2), factorial(n)) for n in range(N + 3)]
inv = [Fraction(0)] * (N + 3); inv[0] = Fraction(1)
for n in range(1, N + 3): inv[n] = -sum(T[k] * inv[n - k] for k in range(1, n + 1))
inv2 = [sum(inv[k] * inv[n - k] for k in range(n + 1)) for n in range(N + 1)]
Tpp = [T[n + 2] * (n + 2) * (n + 1) for n in range(N + 1)]
X = [factorial(n) * sum(Tpp[k] * inv2[n - k] for k in range(n + 1)) for n in range(N + 1)]
def strong(m, adj):
    seen = {0}; st = [0]
    while st:
        u = st.pop()
        for v in range(m):
            if adj[u][v] and v not in seen: seen.add(v); st.append(v)
    if len(seen) < m: return False
    seen = {0}; st = [0]
    while st:
        u = st.pop()
        for v in range(m):
            if adj[v][u] and v not in seen: seen.add(v); st.append(v)
    return len(seen) == m
for m in range(2, 8):
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    cnt = 0
    for bits in product((0, 1), repeat=len(pairs)):
        adj = [[0] * m for _ in range(m)]
        for (i, j), b in zip(pairs, bits):
            if b: adj[i][j] = 1
            else: adj[j][i] = 1
        adj[0][1] = adj[1][0] = 1
        if strong(m, adj): cnt += 1
    print(f"m={m}: X(m-2) = {X[m-2]}   brute force (strong after making pair 01 two-way) = {cnt}   {'OK' if cnt == X[m-2] else 'MISMATCH'}")
