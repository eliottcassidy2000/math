#!/usr/bin/env python3
"""Is the polygonal order N_min(D) bounded in terms of k? Census of max N_min over all primitive k-sets of speeds <= M, with
the argmax sets and the N_min histogram; k = 2 (M = 200), k = 3 (M = 45), k = 4 (M = 24). N_min computed up to 2000."""
import itertools, math, sys
from fractions import Fraction as F
from collections import Counter
def dist(x):
    x = x - (x.numerator // x.denominator); return min(x, 1 - x)
def lonely_vertex_time(D, n, N):
    thr = F(1, n)
    for j in range(1, N):
        if all(dist(F(d * j, N)) >= thr for d in D): return True
    return False
def N_min(D, n, Nmax=2000):
    for N in range(2, Nmax + 1):
        if lonely_vertex_time(D, n, N): return N
    return None
for k, M in ((2, 200), (3, 45), (4, 24)):
    n = k + 1; hist = Counter(); mx = (0, []); cnt = 0
    for D in itertools.combinations(range(1, M + 1), k):
        if math.gcd(*D) != 1: continue
        cnt += 1; Nm = N_min(D, n); hist[Nm] += 1
        if Nm > mx[0]: mx = (Nm, [D])
        elif Nm == mx[0] and len(mx[1]) < 8: mx[1].append(D)
    print(f"k={k} (n={n}), speeds <= {M}: {cnt} primitive sets; max N_min = {mx[0]} at {mx[1]}; histogram {dict(sorted(hist.items()))}", flush=True)
