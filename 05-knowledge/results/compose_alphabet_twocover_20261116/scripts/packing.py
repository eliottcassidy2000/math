#!/usr/bin/env python3
"""packing.py (C-7): the odd-cycle packing polynomial G_T(x, y) = sum over spanning packings by directed odd cycles
(length >= 3) and singletons of x^cycles y^singletons, enumerated directly for every class of orders 3-7, against the
q-slot cover counts P_q computed by DP: h = G(2, 1) (THM-002), P_2 = G(4, 2), P_3 = G(6, 3); and the signed two-cover
identity sum_W (-1)^|W| h(W) h(W^c) = 0."""
import numpy as np
from collections import Counter
exec(open('nearjoin_check.py').read().split('# (1)')[0])
P = np.load('profiles8.npz')
def covers(T):
    n = len(T); N = 1 << n; full = N - 1; dp, Hm = ham_all(T)
    P2 = sum(Hm[W] * Hm[full & ~W] for W in range(N)); alt = sum((-1) ** bin(W).count('1') * Hm[W] * Hm[full & ~W] for W in range(N)); P3 = 0
    for W1 in range(N):
        c = full & ~W1; h1 = Hm[W1]; sub = c
        while True:
            P3 += h1 * Hm[sub] * Hm[c & ~sub]
            if sub == 0: break
            sub = (sub - 1) & c
    return Hm[full], P2, P3, alt
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
        if i == len(cyc):
            N[(c, n - bin(used).count('1'))] += 1; return
        rec(i + 1, used, c)
        if not (used & cyc[i]): rec(i + 1, used | cyc[i], c + 1)
    rec(0, 0, 0); return N
def G(N, x, y): return sum(v * x ** c * y ** s for (c, s), v in N.items())
bad = Counter(); tot = 0; alt_bad = 0
for n in range(3, 8):
    A = P[f'A{n}']
    for i in range(len(A)):
        T = A[i].tolist(); h, p2, p3, alt = covers(T); N = packing_poly(T); tot += 1
        if G(N, 2, 1) != h: bad['q=1'] += 1
        if G(N, 4, 2) != p2: bad['q=2'] += 1
        if G(N, 6, 3) != p3: bad['q=3'] += 1
        if alt != 0: alt_bad += 1
print(f"(C-7) classes of orders 3-7: {tot}; P_q = G(2q, q) failures: {dict(bad) if bad else 0}; signed two-cover sum nonzero: {alt_bad}")
for name, T in (('C_3', [[0,1,0],[0,0,1],[1,0,0]]), ('strong4', [[0,1,1,0],[0,0,1,1],[0,0,0,1],[1,0,0,0]]), ('RT_5', [[1 if (j - i) % 5 in (1, 2) else 0 for j in range(5)] for i in range(5)])):
    print(f"      {name}: packing polynomial {dict(packing_poly(T))}")
