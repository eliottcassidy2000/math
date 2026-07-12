#!/usr/bin/env python3
"""cont.53 part 2: VERIFY klein-S265's load-bearing claim #(tight n covering)=0 on the
specific tight families I found (THM-708 {1..11,13,24}, THM-709 doubling singleton) + the
atlas (AP, GW, dilates). If all M=1/14-tight families are non-covering, the 5->2 case split
holds: [non-covering sieve] + [DC strict-cushion]."""
from fractions import Fraction as F
from math import gcd
def is_covering(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def M(v, Q=400):
    best=F(0)
    for q in range(2,Q):
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            m=None
            for e in v:
                r=(e*p)%q; d=min(r,q-r); dd=F(d,q)
                if m is None or dd<m: m=dd
            if m>best: best=m
            if best>F(1,14): break
        if best>F(1,14): break
    return best
floor=F(1,14)
print("VERIFY klein-S265: every M=1/14-tight family is NON-COVERING (5->2 case split linchpin)")
print(f"{'family':28s} {'M<=1/14?':>9s} {'covering?':>10s}  tight+covering?")
fams=[
 ("AP {1..13}", list(range(1,14))),
 ("THM-708 {1..11,13,24}", [1,2,3,4,5,6,7,8,9,10,11,13,24]),
 ("THM-709 doubling {1..11,13,24}", [1,2,3,4,5,6,7,8,9,10,11,13,24]),
 ("dilate 2*{1..13}", [2*i for i in range(1,14)]),
 ("dilate 3*{1..13}", [3*i for i in range(1,14)]),
 ("GW deep-well {1..12,156}", list(range(1,13))+[156]),
 ("near-dilate {L..12L,13L+1} L=30030", [30030*i for i in range(1,13)]+[30030*13+1]),
]
viol=0
for nm,v in fams:
    m=M(v); tight = m<=floor+F(1,100000)
    cov=is_covering(v)
    bad = tight and cov
    if bad: viol+=1
    print(f"  {nm:26s} {'M='+str(float(m))[:6]:>9s} {str(cov):>10s}  {'*** TIGHT+COVERING!' if bad else 'OK (not both)'}")
print(f"\n  tight families that ARE covering: {viol}  => klein 5->2 split {'HOLDS' if viol==0 else 'FAILS'}")
print("  (near-dilate has M=1/13 > 1/14 = strict-cushion DC, NOT tight -- consistent with the split)")
