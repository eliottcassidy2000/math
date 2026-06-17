from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import random

def maxmin_exact(S):
    cands=set([Fr(0)])
    Sl=list(S)
    for v in S:
        for k in range(2*v+1):
            t=Fr(2*k+1,2*v)
            if 0<=t<1: cands.add(t)
    for a in Sl:
        for b in Sl:
            for ka in range(-1,a+1):
                for kb in range(-1,b+1):
                    if a!=b:
                        t=Fr(ka-kb,a-b)
                        if 0<=t<1: cands.add(t)
                    t2=Fr(ka+kb,a+b)
                    if 0<=t2<1: cands.add(t2)
    best=Fr(0)
    for t in cands:
        vals=[]
        for v in S:
            x=v*t; lower=x.numerator//x.denominator
            d1=x-Fr(lower); d2=Fr(lower+1)-x
            vals.append(d1 if d1<d2 else d2)
        m=min(vals)
        if m>best: best=m
    return best

# Random hunt for max-min < 1/14 (a genuine LRC(14) counterexample)
random.seed(123)
target=Fr(1,14)
mn=None; below=0; spectrum=set()
for _ in range(15000):
    S=sorted(random.sample(range(1,151),13))
    g=reduce(gcd,S)
    S=[x//g for x in S]
    mm=maxmin_exact(S)
    spectrum.add(mm)
    if mn is None or mm<mn[0]: mn=(mm,tuple(S))
    if mm<target: below+=1
print(f"Random 13-sets (15000): min max-min = {mn[0]}={float(mn[0]):.6f}; # with max-min<1/14: {below}")
print(f"  min == 1/14? {mn[0]==Fr(1,14)}")
small=sorted([s for s in spectrum if s<=Fr(1,12)])[:8]
print("  smallest max-min values observed (spectrum above/at 1/14):")
for s in small: print(f"    {s} = {float(s):.6f}")
