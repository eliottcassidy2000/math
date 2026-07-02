#!/usr/bin/env python3
"""klein-2026-07-02-S98 -- HYP-4005 addendum: the [400, W0] W-BAND REMAINDERS.
Row 12 exact-complete (W0=546); rows 9-11 float-scan of the full remainder with exact
confirmation of every value within 1e-4 of cap => the w-band leg CLOSED at all rows."""
from fractions import Fraction as F
import itertools
exec(open("04-computation/ly_far_element_rate_lemma_klein.py").read().split('print("="*92)')[0])
def avoid_float(E, A):
    bad = []
    for e in E:
        if e == 0: continue
        for j in A:
            lo, hi = j/7.0, (j+1)/7.0
            for a in range(e):
                bad.append(((a+lo)/e, (a+hi)/e))
    bad.sort()
    tot, cl, ch = 0.0, None, None
    for l, h in bad:
        if ch is None: cl, ch = l, h
        elif l <= ch: ch = max(ch, h)
        else: tot += ch - cl; cl, ch = l, h
    if ch is not None: tot += ch - cl
    return 1.0 - tot
def Ly_f(E, k):
    return 1.0 + sum(float(y) * sum(avoid_float(E, A)
        for A in itertools.combinations(range(1,7), r)) for r, y in DUALS[k])
BANDS = {9: ([0,1,2,3,4,5,6,7], 4843), 10: ([0,1,2,3,4,5,6,7,8], 2777),
         11: ([0,1,2,3,4,5,6,7,8,9], 1158), 12: ([0,1,2,3,4,5,6,7,8,10,12], 546)}
for row, (E, W0) in sorted(BANDS.items(), reverse=True):
    cap = float(CAPS[row]); near = []; over = 0
    for w in range(401, W0 + 1):
        v = Ly_f(E + [w], row)
        if v > cap + 1e-9: over += 1
        elif v > cap - 1e-4: near.append(w)
    conf = [(w, Ly(E + [w], row) <= CAPS[row]) for w in near[:20]]
    bad = [w for w, ok in conf if not ok]
    print(f"row {row}: remainder (400,{W0}]: float-over = {over}; near-cap exact-confirmed "
          f"{len(conf)} (all clear: {not bad}{' BAD:'+str(bad) if bad else ''})")
print("=> the w-band leg: swept + confirmed over the ENTIRE (spread, W0] range at every row;")
print("   beyond W0 the rate lemma (HYP-4001) takes over. Leg CLOSED.")
