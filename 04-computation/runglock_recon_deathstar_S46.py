#!/usr/bin/env python3
"""Recon for HYP-7225 (death-star-S46): the rung lock + divisor collapse.

(A) RUNG LOCK: runners i,j failing at (q,p) with v_j = m*v_i, m <= 12
    => witnesses lock exactly: w_j = m*w_i.  Predicted from
    |w_j - m*w_i|*q < (m+1)*q/14 < q.  Also probe m = 13 for sharpness
    (the bound (m+1)/14 hits 1: lock may fail).
(B) DIVISOR COLLAPSE: on towers, deep instants p/q reduce to ray denominators
    dividing the tower scale; count distinct components K per family window.
Failing runner: frac(v*p/q) < 1/14 or > 13/14  <=>  exists w: 14|v*p - w*q| < q.
Witness: w = round(v*p/q).
"""
from math import gcd
from fractions import Fraction as F
import random, itertools
random.seed(46)

def failing_witness(v, q, p):
    """Return witness w if runner fails the safe band at p/q, else None."""
    r = (v * p) % q
    if 14 * r < q:
        return (v * p) // q
    if 14 * r > 13 * q:
        return (v * p) // q + 1
    return None

# ---------- (A) the rung lock ----------
lock_checked, lock_fails, sharp13_breaks = 0, 0, []
for trial in range(200000):
    m = random.randint(2, 13)
    v1 = random.randint(1, 400)
    v2 = m * v1
    q = random.randint(2, 1000)
    p = random.randint(1, q - 1)
    w1 = failing_witness(v1, q, p)
    w2 = failing_witness(v2, q, p)
    if w1 is None or w2 is None:
        continue
    if m <= 12:
        lock_checked += 1
        if w2 != m * w1:
            lock_fails += 1
            if lock_fails <= 3:
                print("LOCK FAIL", m, v1, q, p, w1, w2)
    else:
        if w2 != m * w1 and len(sharp13_breaks) < 3:
            sharp13_breaks.append((v1, q, p, w1, w2))
print(f"(A) rung lock m<=12: {lock_checked} double-failing pairs, fails={lock_fails}",
      "PASS" if lock_fails == 0 else "FAIL")
print(f"(A') m=13 sharpness: lock-breaking examples found: {len(sharp13_breaks)}",
      sharp13_breaks[:2])

# ---------- (B) divisor collapse on 7-towers ----------
# family: 13 runners v_i = v_bot * 7^i is astronomically steep; use the
# DENSE-CORE-like shape instead: a single comparable block of 13 speeds
# v_bot * m_1 * ... * m_i with consecutive integer ratios in [2,7].
fam_stats = []
for fam in range(60):
    # keep towers moderate so the exhaustive window loop is tractable:
    # v_top <= ~2e5, D = floor(sqrt(7 v_top)) <= ~1200, capped at 300.
    v_bot = random.randint(1, 4)
    ratios = [random.choice([2, 2, 2, 3, 3, 7]) for _ in range(12)]
    v = [v_bot]
    for m in ratios:
        v.append(v[-1] * m)
    v_top = v[-1]
    # window: D^2 <= 7*v_top (super-ladder in the top runner's width)
    D = min(int((7 * v_top) ** 0.5), 300)
    comps = set()          # distinct deep instants (components)
    carriers = set()
    n_deep = 0
    denom_ok = True
    for q in range(2, D + 1):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            fails = [i for i in range(13) if failing_witness(v[i], q, p) is not None]
            if len(fails) >= 6:
                n_deep += 1
                comps.add(F(p, q))
                carriers.add(q)
                # divisor collapse: reduced denominator q divides the scale
                # v_bot * prod(ratios up to top failing rung)?
                scale = v[max(fails)]
                if scale % q != 0:
                    denom_ok = False
    fam_stats.append((v_bot, tuple(ratios), D, len(carriers), len(comps), n_deep, denom_ok))

worst_K = max(f[4] for f in fam_stats)
worst_car = max(f[3] for f in fam_stats)
all_denom = all(f[6] for f in fam_stats)
zero_deep = sum(1 for f in fam_stats if f[5] == 0)
print(f"(B) 60 block families: zero-deep {zero_deep}/60; worst #carriers={worst_car}; "
      f"worst #components K={worst_K}; carriers divide top-failing scale: "
      f"{'PASS' if all_denom else 'FAIL'}")
for f in sorted(fam_stats, key=lambda x: -x[4])[:5]:
    print("   worst-K families:", f[:1], "ratios", f[1][:4], "... D", f[2],
          "carriers", f[3], "K", f[4], "deep", f[5], "divOK", f[6])
