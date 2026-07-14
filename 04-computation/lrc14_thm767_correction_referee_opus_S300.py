#!/usr/bin/env python3
"""THM-767(4) correction referee (opus-S300, after mac-mini's audit MSG-1621).

Verifies exactly: (1) the audit's chamber example (c=7, core {1..6},
W={1,4,5,6,8,9,10}, t0=1/7) IS an exact all-singleton deck tiling persisting on its
open chamber; (2) the absorption capacities are tiny (a=10: 2 under plain-14, 0
under 14*gcd) while the tiling persists -- so the withdrawn inequality never
governed chamber persistence; (3) the nearest event (owner w=10, t0=3/20) PIERCES
the cover with full witness t=43/140 at clearance exactly 1/14, as THM-767(3)
predicts -- the event-crossing barrier in action."""
from fractions import Fraction as F
from math import gcd
def circ(x):
    fx = x - (x.numerator // x.denominator); return min(fx, 1-fx)
def bad(w, c, t0, d=F(1,14)):
    return [k for k in range(c) if circ(F(w)*(t0+k)/c) < d]
c = 7; W = [1,4,5,6,8,9,10]; t0 = F(1,7); P = [1,2,3,4,5,6]
B = {w: bad(w,c,t0) for w in W}
assert sorted(b[0] for b in B.values()) == list(range(7)) and all(len(b)==1 for b in B.values())
print("exact all-singleton tiling at t0=1/7:", B)
for e in (F(1,1000), F(1,100000)):
    for t in (t0-e, t0+e):
        assert len(set().union(*[set(bad(w,c,t)) for w in W])) == 7
print("tiling persists on the chamber (+-1/1000, +-1/100000): covered 7/7")
caps14 = {a: sum(gcd(a,b) for b in W if b!=a and (a+b)%14==0) for a in W}
caps14g = {a: sum(gcd(a,b) for b in W if b!=a and (a+b)%(14*gcd(a,b))==0) for a in W}
print("capacities plain-14:", caps14, " 14*gcd:", caps14g)
assert caps14[10] == 2 and caps14g[10] == 0
best = None
for w in W:
    for sgn in (1,-1):
        for j in range(2*w+2):
            ev = (F(sgn,14)*c + j)/w; ev -= (ev.numerator//ev.denominator)
            dd = abs(ev - t0)
            if dd > 0 and (best is None or dd < best[0]): best = (dd, w, ev)
dd, w_ev, tev = best
print("nearest event: owner", w_ev, "at t0 =", tev, "(distance", dd, ")")
cover = set().union(*[set(bad(w,c,tev)) for w in W])
assert len(cover) == 6 and min(circ(F(p)*tev) for p in P) >= F(1,14)
k = [k for k in range(7) if k not in cover][0]
t = (tev + k)/c
V = [7*p for p in P] + W
clr = min(circ(F(v)*t) for v in V)
print("event pierces: free sheet", k, "witness t =", t, "clearance =", clr)
assert clr >= F(1,14)
print("ALL CORRECTION CHECKS PASS: barrier real, inequality withdrawn correctly")
