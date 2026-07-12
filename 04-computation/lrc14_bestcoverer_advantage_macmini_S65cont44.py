#!/usr/bin/env python3
"""cont.44: quantify 'consec is the best coverer' (cont.41 crux). p0(consec_k) vs iid
occupancy p0_iid (k phases uniform, P(all 7 sectors hit)). If consec p0 > iid, the
three-gap orbit covers BETTER than random -- the Steinhaus advantage, exact per k.
Also: the k=7 phase transition (p0=0 for k<7 -- can't cover 7 sectors)."""
from fractions import Fraction as F
from math import comb
def p0_consec(k):
    E=list(range(k))
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p0=F(0)
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        if len(hit)==7: p0+=b-a
    return p0
def p0_iid(k):
    # P(k iid uniform phases hit all 7 sectors) = occupancy: sum_{j} (-1)^j C(7,j) ((7-j)/7)^k
    return sum((-1)**j*comb(7,j)*F(7-j,7)**k for j in range(8))
print("consec (three-gap orbit) vs iid coverer: P(all 7 sectors hit) = p0")
print(f"{'k':>3s} {'p0(consec)':>12s} {'p0(iid)':>12s} {'advantage':>11s}  (consec better?)")
for k in range(5,13):
    pc=p0_consec(k); pi=p0_iid(k)
    adv=pc-pi
    print(f"{k:>3d} {float(pc):12.5f} {float(pi):12.5f} {float(adv):+11.5f}  {'YES (three-gap wins)' if adv>0 else 'iid better' if adv<0 else 'p0=0 (k<7: impossible)'}")
print()
print("KEY: consec's coverage ADVANTAGE over random quantifies the three-gap/Steinhaus")
print("structure. If consec >= iid for all k, that's the mechanism behind 'AP is best coverer'")
print("(the POS/T1 dominant term, cont.41). The k=7 transition = apex prime 7 = 14/2.")
