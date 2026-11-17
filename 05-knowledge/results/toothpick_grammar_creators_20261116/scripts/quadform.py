#!/usr/bin/env python3
"""quadform.py (R-1): the count of a one-vertex extension is a quadratic form on the cube.  For a tournament T with
h = H(T), e_x = paths ending at x, B_x = root response, and M(y, x) = number of ordered pairs (P_1, P_2) of nonempty paths
covering T with P_1 ending at y and P_2 starting at x (y != x), and v a new vertex with out-neighbourhood X:
    H(T + v_X) = h + sum_{x in X} (B_x - e_x) - sum_{y, x in X, y != x} M(y, x).
Exhaustive check on all tournaments of orders 1-7 and all X; also the specialisations X = {} (sink), X = V (source),
X = {t} (sole exit, THM-4862 (a)), and the image sizes |{H(T + v_X)}| per class."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
from collections import Counter
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def profile(T):
    n = len(T); N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T)); full = N - 1
    e = [dp[full][x] for x in range(n)]; s = [dpr[full][x] for x in range(n)]
    M = [[0]*n for _ in range(n)]
    for m in range(1, N):
        c = full & ~m
        if c == 0: continue
        for y in range(n):
            a = dp[m][y]
            if not a or not (m >> y) & 1: continue
            for x in range(n):
                if (c >> x) & 1 and dpr[c][x]: M[y][x] += a * dpr[c][x]
    B = [s[x] + sum(M[y][x] for y in range(n)) for x in range(n)]
    return H[full], e, s, B, M
def extend(T, X):
    n = len(T); U = [r[:] + [0 if (X >> i) & 1 else 1] for i, r in enumerate(T)] + [[1 if (X >> j) & 1 else 0 for j in range(n)] + [0]]
    return U
def classes1(n): return [[[0]]] if n == 1 else ([[[0, 1], [0, 0]]] if n == 2 else classes(n))
bad = tot = 0; sizes = Counter(); maxsize = {}
for n in range(1, 8):
    for T in classes1(n):
        h, e, s, B, M = profile(T); img = set()
        for X in range(1 << n):
            pred = h + sum(B[x] - e[x] for x in range(n) if (X >> x) & 1) - sum(M[y][x] for y in range(n) for x in range(n) if y != x and (X >> y) & 1 and (X >> x) & 1)
            real = ham_all(extend(T, X))[1][-1]; tot += 1; img.add(real)
            if pred != real: bad += 1
        sizes[(n, len(img))] += 1; maxsize[n] = max(maxsize.get(n, 0), len(img))
print(f"(R-1) quadratic-form lemma on all tournaments of orders 1-7, all out-neighbourhoods: {tot} extensions, {bad} mismatches")
print("      largest number of distinct child counts produced by one class, by order of the parent:", dict(sorted(maxsize.items())))
print("      distribution of |image| at order 7:", sorted((k[1], v) for k, v in sizes.items() if k[0] == 7))
