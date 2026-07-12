#!/usr/bin/env python3
"""cont.54: the EXACT min M over primitive COVERING families = the LRC(14) crux value.
Test the candidates: near-dilate {L..12L,13L+1} (1/13), covering deep-wells (14/183?),
and adversarial covering families. Covering = mult of each q in 2..14 (INCL 14)."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
def is_cov(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def M_exact(v, Q=600):
    best=F(0)
    qs=set()
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
print(f"1/14={float(F(1,14)):.6f}, 1/13={float(F(1,13)):.6f}, 14/183={float(F(14,183)):.6f}")
print("candidate covering families + exact M:")
def show(nm,v):
    prim = reduce(gcd,v)==1; cov=is_cov(v)
    m=M_exact(v) if (prim and cov and max(v)<3000) else None
    print(f"  {nm:34s} prim={prim} cov={cov}  M={float(m) if m else 'skip':.6}")
    return m
# near-dilate (small L that is primitive+covering)
for L in [840, 2520]:  # L with even structure
    v=sorted([L*i for i in range(1,13)]+[13*L+1])
    show(f"near-dilate L={L}", v)
# covering deep-well analogs (must include mult of 14)
show("{1..12,14} (drop 13, add 14?)", sorted(set(list(range(1,13))+[14])) )
show("{1..11,13,14}", [1,2,3,4,5,6,7,8,9,10,11,13,14])
show("{2..14}", list(range(2,15)))
show("{1..13} AP (non-cov ref)", list(range(1,14)))
# adversarial: minimize M over primitive covering (hill-climb)
random.seed(54)
def rand_cov(Vmax):
    for _ in range(3000):
        v=sorted(set(random.sample(range(1,Vmax),13)))
        if len(v)==13 and reduce(gcd,v)==1 and is_cov(v): return v
    return None
print("\nadversarial MIN-M over primitive covering (Vmax<=60, hill-climb):")
gmin=(F(1),None)
for s in range(10):
    v=rand_cov(60)
    if v is None: continue
    m=M_exact(v)
    for _ in range(120):
        i=random.randrange(13); nd=random.choice([-1,1])
        v2=sorted(set(v[:i]+[v[i]+nd]+v[i+1:]))
        if len(v2)!=13 or min(v2)<1 or reduce(gcd,v2)!=1 or not is_cov(v2): continue
        r=M_exact(v2)
        if r<m: m=r; v=v2
    if m<gmin[0]: gmin=(m,v)
print(f"  adversarial min-M = {gmin[0]} = {float(gmin[0]):.6f} at {gmin[1]}")
print(f"  > 1/14: {gmin[0]>F(1,14)}, margin {float(gmin[0]-F(1,14)):+.6f}; vs 1/13={float(F(1,13)):.6f}")
