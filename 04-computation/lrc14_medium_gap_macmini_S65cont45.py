#!/usr/bin/env python3
"""cont.45: close the medium-spread gap d=21..26 for the k=9 base. Exhaustive per-diameter
min J over primitive 9-sets {0}uS, S subset {1..d}, d in S. Float J with early-exit, exact
confirm of the per-diameter record. Closes the gap between compact (d<=20 exh, cont.42) and
klein's far-element tail (HYP-6070, wide families J>=5.677). floor=432/91, compact-min=5.199."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import sys
def Jf(E):
    pts=set([0.0,1.0])
    for e in E:
        if e==0: continue
        ie=1.0/e
        for m in range(e):
            base=m*ie
            for s in range(8): pts.add(base+s*ie/7)
    pts=sorted(x for x in pts if 0.0<=x<=1.0); tot=0.0
    for i in range(len(pts)-1):
        a=pts[i]; b=pts[i+1]; mid=(a+b)*0.5; hit=0; seen=0
        for e in E:
            if e==0:
                if not (seen&1): seen|=1; hit+=1
                continue
            sec=int((mid*e-int(mid*e))*7)
            bit=1<<sec
            if not (seen&bit): seen|=bit; hit+=1
        N=7-hit; tot+=(b-a)*N*(7-N)
    return tot
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
compact_min=float(F(1019,196)); floor=float(F(432,91))
print(f"compact-min J({{0..8}}) = {compact_min:.4f}; floor 432/91 = {floor:.4f}; klein tail (wide) >= 5.677")
print(f"medium-gap exhaustive min-J per diameter d=21..26 (float, exact-confirm record):")
for d in range(21,27):
    best=9.9; barg=None; cnt=0
    for S in combinations(range(1,d),7):
        E=[0]+list(S)+[d]
        if reduce(gcd,E)!=1: continue
        cnt+=1
        j=Jf(E)
        if j<best: best=j; barg=E
    je=Jexact(barg)
    tag = "OK >= compact-min" if float(je)>=compact_min-1e-9 else ("OK >= floor" if float(je)>=floor else "*** BELOW FLOOR")
    print(f"  d={d}: {cnt} primitive fams, min J = {float(je):.4f} at {barg}  {tag}", flush=True)
print("\nif all d=21..26 >= compact-min (5.199), the medium gap is CLOSED and k=9 base verified")
print("end-to-end: [compact d<=20 exh] + [medium 21-26 exh] + [tail >26 klein HYP-6070/THM-688].")
