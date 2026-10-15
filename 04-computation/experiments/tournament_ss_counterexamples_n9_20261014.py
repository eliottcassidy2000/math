#!/usr/bin/env python3
"""Counterexamples to the stochastic symmetry (SS) at n = 9 (HYP-9280): random sampling (seed 2026) of strong tournaments
n = 8, 9 finds two with a vertex w (in-degree 7) where Pr[N >= 4] < Pr[N <= 3]; for each, E[N] > d/2, the best
in-neighbour has p > 1/2, and M(T) is not transitive.  Session opus-2026-10-14-S36.
"""
import sys, time, random
import numpy as np
from fractions import Fraction
import importlib.util, os
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ss', os.path.join(HERE, 'tournament_inneighbour_symmetry_20261014.py'))
ss = importlib.util.module_from_spec(spec); sys.argv = ['x']; spec.loader.exec_module(ss)
spec2 = importlib.util.spec_from_file_location('mc', os.path.join(HERE, 'tournament_majority_cycle_20261014.py'))
mc = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(mc)
random.seed(2026)
found = 0
for n, trials in ((8, 6000), (9, 1500)):
    for _ in range(trials):
        adj = np.zeros(n, dtype=np.int64)
        p = random.choice([0.5, 0.8, 0.9])
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < p: adj[i] |= 1 << j
                else: adj[j] |= 1 << i
        if not ss.is_strong(n, adj): continue
        F, B = ss.fb_tables(n, adj)
        for w in range(n):
            dist, d = ss.n_distribution(n, adj, w, F, B)
            le = 0; bad = None
            for k in range(d + 1):
                le += dist[k]; ge = int(dist[d - k:d + 1].sum())
                if ge < le: bad = k; break
            if bad is not None:
                H, C = mc.pair_counts(n, adj)
                I = [u for u in range(n) if u != w and adj[u] >> w & 1]
                EN = Fraction(sum(int(dist[k]) * k for k in range(d + 1)), int(dist.sum()))
                maxp = max(Fraction(int(C[u, w]), int(H)) for u in I)
                trans, ninv, fwd3 = mc.majority_info(n, adj)
                print(f'n={n} adj={list(map(int, adj))} w={w} d={d} dist={list(map(int, dist[:d+1]))} k={bad}: E[N]={EN} (d/2={Fraction(d,2)}), max_u p(u<w)={maxp}, M transitive={trans}, #inverted={ninv}, forward 3-cycle={fwd3}')
                found += 1
print('violations listed:', found)
