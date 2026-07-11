#!/usr/bin/env python3
"""THM-709: enumerate the doubling-tight locus. For each e in {8..13}: the family
AP[e->2e] = {1..13} \ {e} u {2e}. Exact M-scan over all p/q, q < Q. Also the e=7 case
(2e = 14: covering flips) and double-doublings. Mechanism check: the ghost-runner
constraint ||2e t|| >= 1/14 forces ||e t|| >= 1/28 (half-clearance ghost)."""
from fractions import Fraction as F
def Mscan(v, Q=500):
    best=F(0); arg=None
    for q in range(2,Q):
        for p in range(1,q):
            f=F(p,q)
            if f.denominator!=q: continue
            m=None
            for x in v:
                fr=x*f-(x*f).__floor__()
                d=min(fr,1-fr)
                if m is None or d<m: m=d
                if m is not None and m<best: break
            if m is not None and m>best: best,arg=m,f
    return best,arg
print("THE DOUBLING FAMILIES AP[e -> 2e], exact M (q < 500):")
for e in range(7,14):
    v=sorted([x for x in range(1,14) if x!=e]+[2*e])
    M,arg=Mscan(v)
    tag="TIGHT (=1/14)" if M==F(1,14) else f"NOT tight ({M} at t={arg})"
    cov = all(any(x%q==0 for x in v) for q in range(2,15))
    print(f"  e={e:2d} (2e={2*e}): {tag}   covering={cov}")
print()
print("double-doublings (two ghosts):")
for e1,e2 in [(11,12),(12,13),(10,12)]:
    v=sorted([x for x in range(1,14) if x not in (e1,e2)]+[2*e1,2*e2])
    M,arg=Mscan(v)
    print(f"  e={e1},{e2}: M={M} {'TIGHT' if M==F(1,14) else f'at t={arg}'}")
