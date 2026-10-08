#!/usr/bin/env python3
"""Audit G, THM-4607: expanding Matthews-Watts maps (Lambda > 0) on Z_3, Z_5: estimate q(e) = P(y, y+e merge at
equal time) by direct integer orbits (y uniform mod d^(T+64)), for |e| from 1 to 4096.  Prediction: q(e) -> 0.
Also records merges at nonzero debt (must be 0: such meetings are null) and the latest merge time."""
import time, math
from mwsim import MW, run_pairs
maps = [
    ('Z3 (1,5,7) r=(0,1,1) rank 2', 3, [1, 5, 7], [0, 1, 1]),
    ('Z5 (3,4,6,7,9) r=(0,1,3,4,4) rank 3', 5, [3, 4, 6, 7, 9], [0, 1, 3, 4, 4]),
    ('Z3 (1,4,16) r=(0,2,1) rank 1, translation-only', 3, [1, 4, 16], [0, 2, 1]),
    ('Z5 (1,1,6,11,56) r=(0,4,3,2,1)? rank 3 translation-only expanding', 5, [1, 1, 6, 11, 56], None),
]
t0 = time.time()
for name, d, m, r in maps:
    if r is None: r = [(-m[i] * i) % d for i in range(d)]
    mw = MW(d, m, r)
    Lam = sum(math.log(x / d) for x in m) / d
    print(f"{name}: r = {r}, Lambda = {Lam:+.4f}", flush=True)
    for e in (1, 2, 3, 4, 8, 16, 32, 64, 256, 1024, 4096, 65536):
        T = 1024; N = 2000
        res = run_pairs(mw, e, T, N, 7000 + e + 13 * d, [64, 256, 1024], window_visits=False)
        q = res['merged_by'][1024] / N; se = (q * (1 - q) / N) ** .5
        mt = sorted(res['merge_times'])
        print(f"   e = {e:6d}: q(<=64) {res['merged_by'][64]/N:.4f}  q(<=1024) {q:.4f} +- {se:.4f}  latest merge {mt[-1] if mt else None}"
              f"  nonzero-debt merges {res['nonzero_debt_merges']}  [{time.time()-t0:.0f}s]", flush=True)
