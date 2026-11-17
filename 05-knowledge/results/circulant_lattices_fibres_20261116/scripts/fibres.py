#!/usr/bin/env python3
"""fibres.py (H-3): distinct packing polynomials per order (orders <= 8: equal to distinct (h, P_2, P_3)), and the
realised coefficient vectors (N_{c,s}) at orders 5-7 with h as the linear form 1 + sum 2^c N_{c,s}."""
import numpy as np
from collections import Counter, defaultdict
exec(open('nearjoin_check.py').read().split('# (1)')[0])
P = np.load('profiles8.npz')
def covers(T):
    n = len(T); N = 1 << n; full = N - 1; dp, Hm = ham_all(T)
    P2 = sum(Hm[W] * Hm[full & ~W] for W in range(N)); P3 = 0
    for W1 in range(N):
        c = full & ~W1; h1 = Hm[W1]; sub = c
        while True:
            P3 += h1 * Hm[sub] * Hm[c & ~sub]
            if sub == 0: break
            sub = (sub - 1) & c
    return Hm[full], P2, P3
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
    rec(0, 0, 0); return N
for n in range(3, 9):
    A = P[f'A{n}']; trip = set(covers(A[i].tolist()) for i in range(len(A)))
    print(f"order {n}: classes {len(A)}, distinct (h, P_2, P_3) = distinct packing polynomials: {len(trip)}")
for n in (5, 6, 7):
    A = P[f'A{n}']; vecs = Counter()
    for i in range(len(A)):
        N = packing_poly(A[i].tolist()); vecs[tuple(sorted(N.items()))] += 1
    keys = sorted(set(k for v in vecs for k, _ in v))
    print(f"order {n}: realised packing vectors {len(vecs)} (types (cycles, singletons) {keys}); h = 1 + sum 2^c N_(c,s); examples: {list(vecs.items())[:4]}")
