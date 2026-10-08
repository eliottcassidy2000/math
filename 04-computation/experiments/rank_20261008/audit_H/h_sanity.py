#!/usr/bin/env python3
"""audit H: sanity for h_greedy on composite d = 4: level-5 table entries and level-6 blocks from Tables.blocks_at agree
with the memoized residue recursion hcore.Map.A_rec (itself checked against brute-force exact simulation)."""
import random, numpy as np
from hdp import Tables
from hcore import Map
rnd = random.Random(5)
for d, m, K in [(4, [1, 3, 5, 7], 6), (5, [1, 6, 4, 31, 1], 5)]:
    tb = Tables(d, m); mp = Map(d, m)
    # prime basis vs position basis: for (1,6,4,31,1) they differ, so compare in the position basis via a change of basis
    from hdp import position_coords
    rho, v, den = position_coords(m)
    mp.v = [tuple(x) for x in v]; mp.n = rho
    H = tb.Hs(K); memo = {}; bad = 0
    idx = [rnd.randrange(len(H)) for _ in range(400)]
    Mv = np.array([H[i] for i in idx], dtype=np.int64); Ev = np.array([rnd.randrange(d ** K) for _ in idx], dtype=np.int64)
    A = tb.blocks_at(K, Mv, Ev)
    for t in range(len(idx)):
        ref = np.array(mp.A_rec(K, int(Mv[t]), int(Ev[t]), memo))
        if not np.array_equal(ref, A[t]): bad += 1
    print(f"Z_{d} m = {m}: 400 random level-{K} blocks from blocks_at vs memoized recursion: mismatches {bad}")
