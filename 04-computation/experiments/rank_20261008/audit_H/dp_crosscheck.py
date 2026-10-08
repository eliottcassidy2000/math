#!/usr/bin/env python3
"""audit H: cross-check the vectorized DP (hdp.Tables) against the memoized residue recursion (hcore.Map.A_rec), which
a_recursion_check.py verified against brute-force exact simulation.  Maps whose position basis = prime basis."""
import random, numpy as np
from hdp import Tables
from hcore import Map
rnd = random.Random(1)
for d, m, kmax in [(5, [1, 2, 3, 7, 1], 4), (5, [1, 1, 2, 3, 7], 4), (5, [1, 2, 3, 7, 11], 3), (7, [1, 1, 1, 1, 2, 3, 5], 3),
                   (7, [1, 1, 1, 2, 3, 5, 11], 3), (4, [1, 3, 5, 7], 5)]:
    tb = Tables(d, m); mp = Map(d, m)
    assert [tuple(x) for x in tb.v] == [tuple(a - b for a, b in zip(mp.v[j], mp.v[0])) for j in range(d)], "basis mismatch"
    bad = 0; n = 0
    for k in range(1, kmax + 1):
        T = tb.table(k); H = tb.Hs(k)
        memo = {}
        for _ in range(300):
            iM = rnd.randrange(len(H)); e = rnd.randrange(d ** k)
            ref = np.array(mp.A_rec(k, int(H[iM]), e, memo))
            n += 1
            if not np.array_equal(ref, T[iM, e]): bad += 1
    print(f"Z_{d} m = {m}: {n} random (k, h) compared, mismatches {bad}")
