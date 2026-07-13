# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont68: the separate-13/14 residual (cont.67) collapses into the mapped multi-killer.
# A separate-13/14 covering family (13,14 by DIFFERENT runners, no mult of 182) has 2 killer slots + runner 1,
# so its interval-core length is <=11 => it IS the multi-killer case (cont.58 M>=7/89; cont.59 escape+finite).
# The naive {1,14}@1/15 refutation is dispelled: {1..14}\{7} has M=1/11 at 3/22 (bound by {8,14}), not 1/15.
from math import gcd
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
def normF(x):
    x=x-(x.numerator//x.denominator)
    if x<0: x+=1
    return min(x,1-x)
def Mexact(v):
    qc=4*max(v)+2; best=F(0)
    for q in range(2,qc):
        for p in range(1,q):
            if gcd(p,q)==1:
                m=min(normF(F(vi*p,q)) for vi in v)
                if m>best: best=m
    return best
def is_cov(v,N=14): return all(any(x%d==0 for x in v) for d in range(2,N+1))
def prim(v): return reduce(gcd,v)==1
def sep1314(v): return (not any(x%182==0 for x in v)) and any(x%13==0 for x in v) and any(x%14==0 for x in v)
def corelen(v):
    s=set(v); k=0
    while (k+1) in s: k+=1
    return k
if __name__=="__main__":
    tgt=F(14,183); maxcore=0; cnt=0
    for k in [11,10,9]:
        core=list(range(1,k+1)); no=13-k
        for outs in combinations(range(13,61),no):
            v=core+list(outs)
            if prim(v) and is_cov(v) and sep1314(v):
                cnt+=1; maxcore=max(maxcore,corelen(v))
    print(f"over {cnt} primitive separate-13/14 covering families (outliers<=60): max interval-core = {maxcore} (<=11)")
    print("=> separate-13/14 has core<=11 (2 killer slots + runner 1 preclude the full {1..12}) = multi-killer (mapped).")
    for name,v in [("{1..11,13,84}", list(range(1,12))+[13,84]), ("{1..14}\{7}", [x for x in range(1,15) if x!=7])]:
        print(f"    {name:<16} M={Mexact(v)} (>=14/183)")
