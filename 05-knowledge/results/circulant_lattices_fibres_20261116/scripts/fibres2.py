#!/usr/bin/env python3
"""fibres2.py (H-3, H-5): per order n <= 8, the number of distinct diagonal polynomials P(q) (= distinct (h, P_2, P_3) for
n <= 8) and of distinct bivariate packing polynomials G (packings enumerated directly); the classes whose diagonal is
shared with a class of different G (invisible to every uniform read, per the concurrent packet); the realised packing
vectors at orders 5-7."""
import numpy as np
from collections import Counter, defaultdict
exec(open('nearjoin_check.py').read().split('# (1)')[0])
P = np.load('profiles8.npz')
def odd_cycles(T):
    n = len(T); cyc = []
    def dfs(start, v, mask, length):
        for w in range(start + 1, n):
            if T[v][w] and not (mask >> w) & 1:
                m2 = mask | (1 << w)
                if T[w][start] and (length + 1) % 2 == 1 and length + 1 >= 3: cyc.append(m2)
                dfs(start, w, m2, length + 1)
    for s in range(n): dfs(s, s, 1 << s, 1)
    return cyc
def packing_poly(T):
    n = len(T); cyc = odd_cycles(T); N = Counter()
    def rec(i, used, c):
        if i == len(cyc): N[(c, n - bin(used).count('1'))] += 1; return
        rec(i + 1, used, c)
        if not (used & cyc[i]): rec(i + 1, used | cyc[i], c + 1)
    rec(0, 0, 0); return tuple(sorted(N.items()))
def diag(G, n):
    return tuple(sum(v * (2 * q) ** c * q ** s for (c, s), v in G) for q in (1, 2, 3))
for n in range(3, 9):
    A = P[f'A{n}']; Gs = [packing_poly(A[i].tolist()) for i in range(len(A))]; Ds = [diag(G, n) for G in Gs]
    byD = defaultdict(set)
    for G, D in zip(Gs, Ds): byD[D].add(G)
    hidden = sum(1 for G, D in zip(Gs, Ds) if len(byD[D]) >= 2)
    ex = next(((D, sorted(byD[D])) for D in byD if len(byD[D]) >= 2), None)
    print(f"order {n}: classes {len(A)}, distinct diagonals P(q) {len(set(Ds))}, distinct bivariate G {len(set(Gs))}, classes whose diagonal is shared with a different G: {hidden}; example: {ex}")
    if n in (5, 6, 7):
        vecs = Counter(Gs); print(f"   realised packing vectors: {len(vecs)}; first: {list(vecs.items())[:3]}")
