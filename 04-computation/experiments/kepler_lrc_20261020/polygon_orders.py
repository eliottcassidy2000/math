#!/usr/bin/env python3
"""Polygonal order N_min(D) = least N with a lonely vertex time j/N of the regular N-gon (observer lonely: all ||d j/N|| >= 1/n).
Reports per k: extremal sets, the law 'the n-gon time 1/n is lonely iff no speed is divisible by n', tight sets (kappa = 1/n)
and their N_min, the fraction settled by polygons of size <= n, by N in {n, 2n, 3n, 4n}, and the growth of max N_min with M."""
import itertools, math, sys
from fractions import Fraction as F
from collections import Counter
sys.path.insert(0, '.')
from kepler_polygon_lrc import kappa, dist
def lonely_vertex_time(D, n, N):
    return any(all(dist(F(d * j, N)) >= F(1, n) for d in D) for j in range(1, N))
def N_min(D, n, Nmax=600):
    for N in range(2, Nmax + 1):
        if lonely_vertex_time(D, n, N): return N
    return None
for k, M in ((2, 30), (3, 20), (4, 15), (5, 12), (6, 11)):
    n = k + 1; rows = []
    for D in itertools.combinations(range(1, M + 1), k):
        if math.gcd(*D) != 1: continue
        Nm = N_min(D, n); kap, t = kappa(D)
        rows.append((D, Nm, kap, t))
    # law: n-gon works iff no speed divisible by n
    law_ok = all(lonely_vertex_time(D, n, n) == (not any(d % n == 0 for d in D)) for D, _, _, _ in rows)
    tight = [(D, Nm) for D, Nm, kap, t in rows if kap == F(1, n)]
    maxN = max(Nm for _, Nm, _, _ in rows); ext = [D for D, Nm, _, _ in rows if Nm == maxN]
    le_n = sum(1 for _, Nm, _, _ in rows if Nm <= n)
    mult = {m: sum(1 for D, _, _, _ in rows if lonely_vertex_time(D, n, m * n)) for m in (1, 2, 3, 4)}
    cum = {m: sum(1 for D, _, _, _ in rows if any(lonely_vertex_time(D, n, mm * n) for mm in range(1, m + 1))) for m in (1, 2, 3, 4)}
    print(f"k={k} (n={n}), speeds <= {M}: {len(rows)} primitive sets; law 'N=n lonely iff n divides no speed': {law_ok}; tight sets {tight}; max N_min = {maxN} at {ext[:5]}; N_min <= n: {le_n}; settled by N = m n: {mult}; cumulative m n, m <= 1..4: {cum}", flush=True)
    # growth of max N_min with M for k = 2
    if k == 2:
        for MM in (10, 20, 30, 40, 60):
            mx = max((N_min(D, n), D) for D in itertools.combinations(range(1, MM + 1), 2) if math.gcd(*D) == 1)
            print(f"   k=2 max N_min for speeds <= {MM}: {mx}", flush=True)
