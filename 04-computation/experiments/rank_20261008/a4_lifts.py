#!/usr/bin/env python3
"""Full-rank maps on Z_5: the 2-step block matrices depend on the multipliers only through m mod 25 and the log-geometry
(universal at full rank), and on the constants through r mod 25.  For each map below and each of the 5^5 lifts of the
constants r_i mod 25 (r_i = -m_i i mod 5 fixed), test whether the A_4 form Q = 5I - J balances every 2-step block exactly
(block_lmi2.exact_check_vec).  Usage: python3 a4_lifts.py"""
import itertools, math
import numpy as np
from block_balance import int_coords, std_r
from block_lmi2 import last_level_unique, exact_check_vec
Q = [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]
MAPS = [[1, 2, 3, 7, 11], [1, 12, 2, 17, 7], [1, 23, 2, 3, 21], [1, 7, 19, 11, 2], [1, 2, 3, 31, 14], [1, 7, 2, 22, 6], [1, 6, 11, 7, 2]]
for m in MAPS:
    d = 5; rho, v = int_coords(m); assert rho == 4 and v[0] == [0, 0, 0, 0] and all(sorted(map(abs, v[j])) == [0, 0, 0, 1] for j in range(1, 5)), (m, v)
    r0 = std_r(m); ok = 0; bad = []
    for s in itertools.product(range(5), repeat=5):
        r = [r0[i] + 5 * s[i] for i in range(5)]
        U, n, fr = last_level_unique(d, m, r, v, rho, 2)
        good, _ = exact_check_vec(Q, U)
        if good: ok += 1
        else: bad.append(tuple(s))
    print(f"m = {m} (mod 25: {[x % 25 for x in m]}): 5I - J certifies k = 2 for {ok} of 3125 lifts r mod 25; failing lift offsets (first 6): {bad[:6]}", flush=True)
