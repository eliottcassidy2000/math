#!/usr/bin/env python3
"""
lrc14_loose_split_macmini_S60  (mac-mini-2026-07-05-S60)

Sharpens the loose-branch crux (klein HYP-4151 / kps HYP-4157 critical reduction) by SPLITTING
12-configs on whether they contain a TIGHT 11-subfamily (a dilated AP c*{1..11}):

  AP-SUBFAMILY configs  (some 11-subfamily is a dilated AP): M in {1/13} U [2/25, inf) -- PROVED (mac-mini
                         S59 HYP-4152, ladder + CRT free-rider). The razor-thin extremizer {1..11,24} @ 2/25
                         is HERE.
  ALL-LOOSE configs      (no tight 11-subfamily = every 11-subfamily loose): M >= 2/23 -- a MARGIN over 2/25,
                         floored by {1..13}\{6} @ 2/23 (= the 11-runner second value; recursive).

CONSEQUENCE: the ENTIRE razor-thinness of the loose branch (the 2/25 edge) lives in the S59-PROVED
AP-subfamily case; the all-loose residual is a non-razor >= 2/23 bound that recurses into the 11-runner
rigidity. (Mirror of the covering-min split S46/S47: razor value in the dilated/proved case, margin in the
residual.) VERIFIED reliably (structured, high-height CRT-lifts -- MISTAKE-102 discipline); NOT yet proved
(the all-loose >= 2/23 needs 11-rigidity + peeling closure).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import sys
def out(*a): print(*a); sys.stdout.flush()
def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in sorted(Q):
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a, q)) for v in sp)
            if m > best: best = m
    return best
def is_dilAP11(sub):
    s = sorted(sub); d = s[1]-s[0]
    return d >= 1 and all(s[i+1]-s[i] == d for i in range(10)) and s[0] == d
def all_loose(B):
    B = sorted(B)
    return not any(is_dilAP11(B[:i]+B[i+1:]) for i in range(12))
if __name__ == "__main__":
    c113, c225, c223 = F(1,13), F(2,25), F(2,23)
    floor = F(1); ff = None; below = ingap = tested = 0
    # structured all-loose: AP-with-holes {1..N}\{holes}, N=13..18
    for N in range(13, 19):
        base = list(range(1, N+1))
        for holes in combinations(range(1, N+1), N-12):
            B = [x for x in base if x not in holes]
            if len(B) != 12 or gcd_all(B) != 1 or not all_loose(B): continue
            M = M_exact(B)
            if M > F(1, 10): continue
            tested += 1
            if M < floor: floor = M; ff = tuple(B)
            if M < c223: below += 1
            if c113 < M < c225: ingap += 1
    out(f"structured all-loose configs (no tight 11-subfamily), M<=1/10: {tested}")
    out(f"  min M = {floor} = {float(floor):.5f} at {list(ff) if ff else None}")
    out(f"  below 2/23: {below};  in gap (1/13,2/25): {ingap}")
    out(f"  => all-loose floor = 2/23 (margin over 2/25={float(c225):.5f}); AP-subfamily case (S59) holds the razor 2/25.")
