#!/usr/bin/env python3
"""cont.45 part 2: finish d=24,25,26 of the medium-gap check."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
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
        a=pts[i]; b=pts[i+1]; mid=(a+b)*0.5; seen=0; hit=0
        for e in E:
            if e==0:
                if not (seen&1): seen|=1; hit+=1
                continue
            sec=int((mid*e-int(mid*e))*7); bit=1<<sec
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
cm=float(F(1019,196))
for d in [24,25,26]:
    best=9.9; barg=None; cnt=0
    for S in combinations(range(1,d),7):
        E=[0]+list(S)+[d]
        if reduce(gcd,E)!=1: continue
        cnt+=1; j=Jf(E)
        if j<best: best=j; barg=E
    je=Jexact(barg)
    print(f"  d={d}: {cnt} fams, min J = {float(je):.4f} at {barg}  {'OK >= compact-min' if float(je)>=cm-1e-9 else '*** CHECK'}", flush=True)
