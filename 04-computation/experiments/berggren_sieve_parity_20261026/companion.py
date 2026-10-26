# Companion vector of a Berggren word: M_w (1,0)^T = (al, ga), with det M_w = al*n - ga*m = (-1)^#A.
import sys
from math import gcd
U = ((2, -1), (1, 0)); A = ((2, 1), (1, 0)); D = ((1, 2), (0, 1))
def mul(X, Y): return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def word(m, n):
    w = []
    while (m, n) != (2, 1):
        if m > 3 * n: w.append('D'); m -= 2 * n
        elif m > 2 * n: w.append('A'); m, n = n, m - 2 * n
        else: w.append('U'); m, n = n, 2 * n - m
    return ''.join(reversed(w))
G = {'U': U, 'A': A, 'D': D}
stats = {}
ex = []
for m in range(2, 120):
    for n in range(1, m):
        if (m - n) % 2 == 0 or gcd(m, n) != 1: continue
        w = word(m, n)
        M = ((1, 0), (0, 1))
        for ch in w: M = mul(G[ch], M)          # apply first letter first
        assert (M[0][0] * 2 + M[0][1], M[1][0] * 2 + M[1][1]) == (m, n), (m, n, w, M)
        al, ga = M[0][0], M[1][0]
        det = al * n - ga * m
        assert det == (-1) ** w.count('A')
        key = (al >= 0, ga >= 0, al <= m, ga <= n, 2 * al <= m)
        stats.setdefault(key, [0, 0])[w.count('A') % 2] += 1
        if len(ex) < 12: ex.append((m, n, w, al, ga, det, pow(n, -1, m) if m > 1 else 0))
for e in ex: print(e)
print(stats)
