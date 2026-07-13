#!/usr/bin/env python3
"""cont.55: the far-element M-floor for the 12-AP. M({1..12, f}) as f varies -- find min.
Conjecture: min at f=182 (deep well, M=14/183), i.e. adding a far element to {1..12} lowers
M from 1/13 to at most 14/183. Pins the deep well as the crux extremal."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
def M(v, Q=400):
    best=F(0); qs=set()
    for i in range(len(v)):
        for j in range(len(v)):
            if i!=j: qs.add(v[i]+v[j]); qs.add(abs(v[i]-v[j]))
    for q in sorted(set(list(qs)+list(range(2,Q)))):
        if q<2: continue
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            m=None
            for e in v:
                r=(e*p)%q; d=min(r,q-r); dd=F(d,q)
                if m is None or dd<m: m=dd
            if m>best: best=m
    return best
base=list(range(1,13))  # {1..12}, M=1/13
print(f"M({{1..12}}) alone would be 1/13={float(F(1,13)):.5f}; add far element f:")
print(f"{'f':>5s} {'M({1..12,f})':>13s} {'prim?':>6s} {'cov?':>5s}  vs 14/183={float(F(14,183)):.5f}")
def is_cov(v): return all(any(x%q==0 for x in v) for q in range(2,15))
results=[]
for f in [13,14,26,27,40,84,90,120,156,168,169,182,183,196,364,546]:
    v=sorted(base+[f])
    if len(set(v))!=13: continue
    m=M(v); prim=reduce(gcd,v)==1; cov=is_cov(v)
    results.append((m,f,prim,cov))
    star="" if m!=F(14,183) else "  <-- 14/183 (DEEP WELL)"
    print(f"  {f:>5d} {str(m):>13s}={float(m):.5f} {str(prim):>6s} {str(cov):>5s}{star}")
# find the min over a dense f-scan
print("\ndense scan f=13..400 for the MINIMUM M({1..12,f}) (primitive+covering only):")
gmin=(F(1),None)
for f in range(13,400):
    v=sorted(base+[f])
    if len(set(v))!=13: continue
    if reduce(gcd,v)!=1 or not is_cov(v): continue
    m=M(v)
    if m<gmin[0]: gmin=(m,f)
print(f"  MIN M = {gmin[0]}={float(gmin[0]):.6f} at f={gmin[1]}  (14/183={float(F(14,183)):.6f})")
print(f"  => far-element floor = {'14/183 at f=182 (deep well) CONFIRMED' if gmin[0]==F(14,183) else 'DIFFERENT: '+str(gmin[0])}")
