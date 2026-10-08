#!/usr/bin/env python3
"""Audit G, THM-4608: direct 2-adic orbits of (y, y+e), y uniform mod 2^(T+64), for maps in every class.
Prediction (THM-4608): (a) (1,1): a.s. iff s0 = r1 - r0 | e, geometric tail, else NEVER;
(b) (1,3)/(3,1): a.s. iff s/3^v3(s) | e (s = r0 + r1), tail ~ T^(-1/2), else NEVER;
(c) m0 m1 >= 5: q(e) -> 0 as |e| -> inf; q(e) < 1 proved only for px+1 (and offsets in sZ of the (1,p) class)."""
import time, math
from mwsim import MW, run_pairs

def v3(n):
    n = abs(n); k = 0
    while n % 3 == 0: n //= 3; k += 1
    return k

def predict(m0, m1, r0, r1, e):
    if (m0, m1) == (1, 1):
        s0 = r1 - r0; return 'a.s.' if e % s0 == 0 else 'never'
    if (m0, m1) in ((1, 3), (3, 1)):
        s = r0 + r1; sp = abs(s) // 3 ** v3(s); return 'a.s.' if e % sp == 0 else 'never'
    return 'expanding'

t0 = time.time()
rows = [
    # (m0, m1, r0, r1, e, T, N)
    (1, 1, 0, 1, 1, 256, 2000), (1, 1, 0, 1, 7, 256, 2000), (1, 1, 0, 1, 1000, 256, 2000),
    (1, 1, 0, 3, 3, 256, 2000), (1, 1, 0, 3, 6, 256, 2000), (1, 1, 0, 3, 1, 2048, 500), (1, 1, 0, 3, 2, 2048, 500),
    (1, 1, 2, 7, 5, 256, 2000), (1, 1, 2, 7, 1, 2048, 500), (1, 1, 2, 7, 3, 2048, 500), (1, 1, 4, -1, 5, 256, 2000),
    (1, 1, 4, -1, 2, 2048, 500), (1, 1, -2, 1, 3, 256, 2000), (1, 1, -2, 1, 1, 2048, 500),
    (1, 3, 0, 1, 1, 16384, 2000), (1, 3, 0, 1, -1, 4096, 1000), (1, 3, 0, 1, 2, 4096, 1000), (1, 3, 0, 1, 3, 4096, 1000),
    (1, 3, 0, 5, 5, 4096, 1000), (1, 3, 0, 5, 1, 4096, 500), (1, 3, 0, 5, 2, 4096, 500), (1, 3, 0, 5, 3, 4096, 500),
    (1, 3, 2, 1, 1, 16384, 2000), (1, 3, 0, 3, 1, 4096, 1000), (1, 3, 0, 15, 5, 4096, 1000), (1, 3, 0, 15, 3, 4096, 500),
    (1, 3, 0, 15, 1, 4096, 500), (1, 3, -2, 1, 1, 4096, 1000), (1, 3, 6, 3, 1, 4096, 1000),
    (3, 1, 0, 1, 1, 16384, 2000), (3, 1, 2, 3, 5, 4096, 1000), (3, 1, 2, 3, 1, 4096, 500), (3, 1, 4, 5, 1, 4096, 1000),
    (3, 1, 0, 7, 7, 4096, 1000), (3, 1, 0, 7, 1, 4096, 500),
    (1, 5, 0, 1, 1, 2048, 8000), (1, 5, 0, 1, 3, 2048, 4000), (1, 5, 0, 1, 100, 2048, 4000),
    (1, 5, 0, 5, 1, 2048, 8000), (1, 5, 0, 5, 5, 2048, 4000), (1, 5, 2, 1, 1, 2048, 4000),
    (5, 1, 0, 1, 1, 2048, 2000), (5, 1, 0, 1, 3, 2048, 4000),
    (1, 7, 0, 1, 1, 2048, 4000), (1, 7, 0, 7, 1, 2048, 4000), (7, 1, 0, 1, 1, 2048, 4000),
    (3, 3, 0, 1, 1, 2048, 1000), (3, 3, 0, 3, 1, 2048, 4000), (3, 3, 0, 3, 2, 2048, 4000),
    (3, 5, 0, 1, 1, 2048, 4000), (3, 5, 0, 1, 2, 2048, 4000), (3, 5, 0, 1, 3, 2048, 4000), (3, 5, 2, 3, 1, 2048, 4000),
    (5, 3, 0, 1, 1, 2048, 4000), (3, 7, 0, 1, 1, 2048, 4000), (1, 9, 0, 1, 1, 2048, 4000), (9, 1, 0, 1, 1, 2048, 4000),
    (3, 9, 0, 1, 1, 2048, 4000), (5, 5, 0, 1, 1, 2048, 1000), (1, 15, 0, 1, 1, 2048, 4000),
]
for (m0, m1, r0, r1, e, T, N) in rows:
    mw = MW(2, [m0, m1], [r0, r1])
    if predict(m0, m1, r0, r1, e) == 'expanding': T = min(T, 1024)
    cps = sorted({c for c in (16, 64, 256, 1024, 4096, 16384) if c <= T} | {T})
    res = run_pairs(mw, e, T, N, 1000 + 17 * m0 + 31 * m1 + 7 * r0 + 3 * r1 + 101 * e, cps, window_visits=False)
    q = [res['merged_by'][c] / N for c in cps]
    pred = predict(m0, m1, r0, r1, e)
    tail = ''
    if pred == 'a.s.':
        tail = '  sqrt(T)*P(no merge): ' + ', '.join(f"{c}:{math.sqrt(c) * (1 - res['merged_by'][c] / N):.2f}" for c in cps if c >= 64)
    mt = sorted(res['merge_times'])
    print(f"({m0},{m1}) r=({r0},{r1}) e={e:5d} pred={pred:9s} N={N}: merged-by {dict(zip(cps, [round(x, 4) for x in q]))}"
          f"  latest merge {mt[-1] if mt else None}  nonzero-debt merges {res['nonzero_debt_merges']}{tail}  [{time.time()-t0:.0f}s]", flush=True)
# Negative multipliers (allowed by the Matthews-Watts definition: m_i nonzero integers, |prod m_i| vs d^d).
# THM-4608's setting restricts to positive multipliers; these test the scope of "every 2-adic Matthews-Watts map"
# and of "3x+1 is the only 2-adic map (up to affine conjugacy) with diffusive coalescence".
# (the debt-vector bookkeeping ignores signs, so only the merge counts are meaningful here)
print("--- negative multipliers (outside THM-4608's positive setting) ---", flush=True)
for (m0, m1, r0, r1, e, T, N) in [(1, -3, 0, 1, 1, 16384, 2000), (1, -3, 0, 1, 2, 4096, 1000), (-3, 1, 0, 1, 1, 4096, 1000),
                                  (1, -1, 0, 1, 1, 1024, 2000), (-1, -1, 0, 1, 1, 1024, 2000), (-1, 1, 0, 1, 1, 1024, 2000),
                                  (-1, 3, 0, 1, 1, 4096, 1000)]:
    mw = MW(2, [m0, m1], [r0, r1])
    cps = sorted({c for c in (16, 64, 256, 1024, 4096, 16384) if c <= T} | {T})
    res = run_pairs(mw, e, T, N, 5000 + 17 * abs(m0) + 31 * abs(m1) + 101 * e, cps, window_visits=False)
    q = [res['merged_by'][c] / N for c in cps]
    tail = '  sqrt(T)*P(no merge): ' + ', '.join(f"{c}:{math.sqrt(c) * (1 - res['merged_by'][c] / N):.2f}" for c in cps if c >= 64)
    mt = sorted(res['merge_times'])
    print(f"({m0},{m1}) r=({r0},{r1}) e={e:5d} N={N}: merged-by {dict(zip(cps, [round(x, 4) for x in q]))}  latest merge {mt[-1] if mt else None}{tail}  [{time.time()-t0:.0f}s]", flush=True)
