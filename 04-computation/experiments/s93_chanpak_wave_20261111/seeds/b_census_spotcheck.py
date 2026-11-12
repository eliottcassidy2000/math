#!/usr/bin/env python3
"""S93 lane seeds: independent spot check of the order-10 census (b_census.py, min-score extensions).
Random labelled tournaments of order 10 (uniform arcs, and near-transitive ones: a transitive tournament with k
random arcs reversed, k = 3..12), every O, by the plain subset recursion of seedlib (not the census code).  Every unit
pair found must be recorded in b_census_minv.npz: odd part c of the slope with minv[sign][c] <= v_2(slope)."""
import numpy as np
from seedlib import dp_tables, pair_mats, subset_sums

rng = np.random.default_rng(1010)
mv = np.load("b_census_minv.npz")["minv_10"]
n = 10
fe = np.zeros((1 << n, n), np.int64); gs = np.zeros((1 << n, n), np.int64); hS = np.zeros(1 << n, np.int64)
Q = np.zeros((n, n), np.int64); U = np.zeros(n, np.int64)
Wsub = np.zeros(1 << n, np.int64); Msub = np.zeros(1 << n, np.int64)
found = [0, 0]; missing = 0; ntour = 0
for trial in range(30000):
    out = np.zeros(n, np.int64)
    if trial % 2 == 0:
        bits = rng.integers(0, 2, 45)
        b = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits[b]: out[i] |= 1 << j
                else: out[j] |= 1 << i
                b += 1
    else:
        perm = rng.permutation(n)
        adj = np.zeros((n, n), np.int8)
        for i in range(n):
            for j in range(i + 1, n):
                adj[perm[i], perm[j]] = 1
        for _ in range(int(rng.integers(3, 13))):
            i, j = rng.choice(n, 2, replace=False)
            adj[i, j], adj[j, i] = adj[j, i], adj[i, j]
        for i in range(n):
            for j in range(n):
                if adj[i, j]: out[i] |= 1 << j
    ntour += 1
    dp_tables(n, out, fe, gs, hS); pair_mats(n, fe, gs, hS, Q, U); subset_sums(n, Q, U, Wsub, Msub)
    H = int(hS[-1])
    for sg, target in ((0, H - 1), (1, H + 1)):
        for O in np.flatnonzero(Wsub == target):
            if O == 0: continue
            M = int(Msub[O]); c = M; v = 0
            while c % 2 == 0: c //= 2; v += 1
            found[sg] += 1
            if not (mv[sg][c >> 1] <= v):
                missing += 1
                print("NOT IN CENSUS", sg, M)
print(f"{ntour} random tournaments of order 10: {found[0]} seeds and {found[1]} anti-seeds found by the independent recursion; not recorded in the census: {missing}")
