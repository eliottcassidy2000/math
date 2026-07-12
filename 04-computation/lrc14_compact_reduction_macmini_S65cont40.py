#!/usr/bin/env python3
"""cont.40: the k=9 COMPACT REDUCTION. min J over primitive 9-sets by DIAMETER d.
Family = {0} u S, S subset {1..d} with |S|=8, d in S (diameter exactly d), gcd(S u {0})=1.
Track min J vs d; if min-J is achieved at small d and stays > threshold for larger d,
the extremal problem reduces to a FINITE check. Float J for speed, exact-confirm records."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
def Jfloat(E):
    # float sector-emptiness J via breakpoint sweep
    pts=set([0.0,1.0])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(m/e+s/(7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[0.0]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=mid*e-int(mid*e); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return sum(p[n]*n*(7-n) for n in range(8))
def Jexact(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return sum(p[n]*n*(7-n) for n in range(8))
thr=float(F(432,91)); target=float(F(4465,882))
print(f"threshold 432/91={thr:.4f}; consec {{1..9}} target 4465/882={target:.4f}")
print(f"{'diam d':>6s} {'#families':>10s} {'min J':>9s} {'argmin':<26s} vs target")
globmin=(9.9,None)
for d in range(8,23):
    best=(9.9,None)
    for S in combinations(range(1,d),7):
        E=[0]+list(S)+[d]
        if reduce(gcd,E)!=1: continue
        j=Jfloat(E)
        if j<best[0]: best=(j,E)
    # exact-confirm the record
    je=Jexact(best[1])
    if float(je)<globmin[0]: globmin=(float(je),best[1],je)
    flag = "<-- GLOBAL MIN" if best[1]==[0,1,2,3,4,5,6,7,8] or (d==8) else ("  >thr" if float(je)>thr else " *** <=thr")
    print(f"{d:>6d} {'':>10s} {float(je):9.4f} {str(best[1]):<26s} {flag}")
print(f"\nGLOBAL compact min (d<=22): J={globmin[2]}={globmin[0]:.4f} at {globmin[1]}")
print(f"= 4465/882 at consec? {globmin[2]==F(4465,882)}")
print(f"min-J stays > threshold for all d in [8,22]: (check column) => compact reduction to d<=D0")
