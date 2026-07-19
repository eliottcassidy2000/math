#!/usr/bin/env python3
"""lrc14_maximiser_sojourn_bound_deathstar_S58b.py
Exact verification of the CORRECTED maximiser bound (fixing MISTAKE-173):
  sojourn(d) = meas{u in [0,1): g(u)=(frac(-d2 u),frac(-d3 u),frac(-d4 u)) in B} <= 2/21,
with B the six-polytope bad set. Also checks the two proof ingredients:
  (A) exactly-2-mirror-runs structure (#runs<=2),
  (B) single-run <= 1/21.
B (one ordering region, sorted g): min<=2/7, consecutive g-diffs<=2/7, max>=5/7  (verified == bad_from_g).
Exact via affine cells (partition [0,1) at k/d_i and ordering crossings).
"""
import sys
from fractions import Fraction as F
from math import gcd
N = int(sys.argv[1]) if len(sys.argv) > 1 else 16

def bad_intervals(DD):
    bps = {F(0), F(1)}
    for d in DD:
        for k in range(d + 1):
            bps.add(F(k, d))
    for i in range(3):
        for j in range(i + 1, 3):
            dd = abs(DD[i] - DD[j])
            if dd > 0:
                for k in range(dd + 1):
                    bps.add(F(k, dd))
    bps = sorted(bps)
    ivs = []
    for idx in range(len(bps) - 1):
        lo, hi = bps[idx], bps[idx + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        gmid = [(-DD[i] * mid) % 1 for i in range(3)]
        a = [gmid[i] + DD[i] * mid for i in range(3)]
        order = sorted(range(3), key=lambda i: gmid[i])
        s0, s1, s2 = order
        cons = [("ge", a[s0] - F(2, 7), DD[s0]), ("le", a[s2] - F(5, 7), DD[s2])]
        feas = True
        for (x, y) in [(s1, s0), (s2, s1)]:
            c = a[x] - a[y]; dc = DD[x] - DD[y]
            if dc > 0: cons.append(("ge", c - F(2, 7), dc))
            elif dc < 0: cons.append(("le", c - F(2, 7), dc))
            elif not (c <= F(2, 7)): feas = False
        if not feas:
            continue
        ulo, uhi = lo, hi
        for typ, const, dd in cons:
            b = const / dd
            if typ == "ge": ulo = max(ulo, b)
            else: uhi = min(uhi, b)
        if uhi > ulo:
            ivs.append((ulo, uhi))
    ivs.sort(); merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged

maxsoj = F(0); dsoj = None; maxruns = 0; over2 = 0; maxrun = F(0); drun = None
for d2 in range(1, N + 1):
    for d3 in range(1, N + 1):
        for d4 in range(1, N + 1):
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            ivs = bad_intervals([d2, d3, d4])
            soj = sum(hi - lo for lo, hi in ivs)
            if soj > maxsoj: maxsoj = soj; dsoj = (d2, d3, d4)
            if len(ivs) > maxruns: maxruns = len(ivs)
            if len(ivs) > 2: over2 += 1
            for lo, hi in ivs:
                if hi - lo > maxrun: maxrun = hi - lo; drun = (d2, d3, d4)
print("primitive directions with components 1..%d:" % N)
print("  MAX sojourn  = %s = %.7f at %s  (2/21=%.7f) -> ceiling %s"
      % (maxsoj, float(maxsoj), dsoj, 2 / 21, "HOLDS" if maxsoj <= F(2, 21) else "BROKEN"))
print("  MAX #runs    = %d ; directions with >2 runs = %d" % (maxruns, over2))
print("  MAX one run  = %s = %.7f at %s  (1/21=%.7f) -> %s"
      % (maxrun, float(maxrun), drun, 1 / 21, "run<=1/21 HOLDS" if maxrun <= F(1, 21) else "BROKEN"))
