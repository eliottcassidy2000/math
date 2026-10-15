#!/usr/bin/env python3
"""Audit E (2026-10-14): independent Python re-check of the (Z2), k = 2 witnesses found by strong_dfs_bloc_z.c:
z is a cut vertex with T - z = D_1 => D_2 (two strong components), and z's best balanced partner has balance
min(P, 1-P) = 10/29 (n = 7), 18/53 (n = 8), 34/101 (n = 9) < 6/17, contradicting 'For k = 2, the best-partner balance of z
stayed >= 6/17' in HYP-9279 (Z) as a description of the n <= 9 data."""
from fractions import Fraction as F
def d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63; bits += [(v >> (5 - t)) & 1 for t in range(6)]
    return n, [[bits[i * n + j] for j in range(n)] for i in range(n)]
def reach(A, S, s):
    r = {s}; st = [s]
    while st:
        v = st.pop()
        for w in S:
            if A[v][w] and w not in r: r.add(w); st.append(w)
    return r
for s, z in [('&FZfpk[EO@?', 0), ('&G\\r{^BwSAAOO', 5), ('&H]x~Fw^@{D?OH?_', 6)]:
    n, A = d6(s)
    V = set(range(n)); R = V - {z}
    strongT = all(reach(A, V, v) == V for v in V)
    comps = []
    for v in R:
        C = frozenset(u for u in reach(A, R, v) if v in reach(A, R, u))
        if C not in comps: comps.append(C)
    paths = []
    def rec(p, used):
        if len(p) == n: paths.append(tuple(p)); return
        for v in range(n):
            if not used >> v & 1 and A[p[-1]][v]:
                p.append(v); rec(p, used | 1 << v); p.pop()
    for s0 in range(n): rec([s0], 1 << s0)
    H = len(paths)
    best = max(min(F(sum(1 for p in paths if p.index(z) < p.index(v)), H), 1 - F(sum(1 for p in paths if p.index(z) < p.index(v)), H)) for v in R)
    print(f"{s} n={n} z={z}: T strong {strongT}; components of T - z: {[sorted(C) for C in comps]} (k = {len(comps)}); H = {H}; "
          f"best-partner balance of z = {best} = {float(best):.5f}; below 6/17 = {float(F(6,17)):.5f}: {best < F(6, 17)}")
