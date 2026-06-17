from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import random

def maxmin_exact(S):
    cands=set()
    Sl=list(S)
    for v in S:
        for k in range(v+1):
            cands.add(Fr(2*k+1,2*v))
    # crossings restricted to plausible range
    for ai in range(len(Sl)):
        a=Sl[ai]
        for bi in range(ai+1,len(Sl)):
            b=Sl[bi]
            for ka in range(a+1):
                for kb in range(b+1):
                    t=Fr(ka+kb,a+b)
                    if 0<t<1: cands.add(t)
                    if a!=b:
                        t2=Fr(ka-kb,a-b)
                        if 0<t2<1: cands.add(t2)
    best=Fr(0)
    for t in cands:
        if not (0<t<1): continue
        m=None
        for v in S:
            x=v*t; lo=x.numerator//x.denominator
            d1=x-lo; d2=(lo+1)-x
            dv=d1 if d1<d2 else d2
            if m is None or dv<m: m=dv
            if m<best: break
        if m>best: best=m
    return best

random.seed(123)
target=Fr(1,14)
mn=None; below=0; spectrum=set()
N=4000
for _ in range(N):
    S=sorted(random.sample(range(1,121),13))
    g=reduce(gcd,S); S=[x//g for x in S]
    mm=maxmin_exact(S)
    spectrum.add(mm)
    if mn is None or mm<mn[0]: mn=(mm,tuple(S))
    if mm<target: below+=1
print(f"Random 13-sets ({N}): min max-min = {mn[0]}={float(mn[0]):.6f}; # with max-min<1/14: {below}")
print(f"  min == 1/14? {mn[0]==Fr(1,14)}")
# verify AP and sporadic
print(f"  AP max-min = {maxmin_exact(list(range(1,14)))}")
print(f"  sporadic max-min = {maxmin_exact([1,2,3,4,5,6,7,8,9,10,11,13,24])}")
small=sorted([s for s in spectrum if s<=Fr(1,12)])[:6]
print("  smallest max-min spectrum at/above 1/14:")
for s in small: print(f"    {s} = {float(s):.6f}")
