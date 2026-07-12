#!/usr/bin/env python3
"""cont.46: FINISH the k=8 base. Exhaustive MAX deg-3-majorant bound per diameter d=7..26
over primitive 8-cores {0}uS, S subset {1..d}, |S|=6, d in S. Bound (UPPER bound on Phi)
must stay <= cap_9 = 1979/4004; hardest at consec {0..7}, decreases with diameter. Float
bound with exact-confirm of the per-diameter MAX record."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
def boundf(E):  # float deg-3 majorant value
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
                if not(seen&1): seen|=1; hit+=1
                continue
            sec=int((mid*e-int(mid*e))*7); bit=1<<sec
            if not(seen&bit): seen|=bit; hit+=1
        N=7-hit
        q3=1 - (2/3)*N + (47/252)*N*(N-1) - (5/252)*N*(N-1)*(N-2)
        tot+=(b-a)*q3
    return tot
def bounde(E):
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
    def q3(N): return 1 - F(2,3)*N + F(47,252)*N*(N-1) - F(5,252)*N*(N-1)*(N-2)
    return sum(p[N]*q3(N) for N in range(8))
cap9=float(F(1979,4004))
print(f"k=8: deg-3 bound must stay <= cap_9 = {cap9:.4f}; hardest consec {{0..7}}=0.4380")
print("exhaustive MAX bound per diameter d=7..26 (float, exact-confirm record):")
gmax=(0.0,None)
for d in range(7,27):
    best=-1.0; barg=None; cnt=0
    for S in combinations(range(1,d),6):
        E=[0]+list(S)+[d]
        if reduce(gcd,E)!=1: continue
        cnt+=1; v=boundf(E)
        if v>best: best=v; barg=E
    ve=bounde(barg)
    if float(ve)>gmax[0]: gmax=(float(ve),barg)
    tag="OK <= cap_9" if float(ve)<=cap9 else "*** EXCEEDS cap_9"
    print(f"  d={d:2d}: {cnt:>7d} fams, MAX bound = {float(ve):.4f} at {barg}  {tag}", flush=True)
print(f"\nGLOBAL max bound (d<=26) = {gmax[0]:.4f} at {gmax[1]}  vs cap_9 = {cap9:.4f}")
print(f"k=8 base HOLDS (max <= cap_9): {gmax[0]<=cap9}, margin {cap9-gmax[0]:+.4f}")
print("=> [compact+medium exhaustive d<=26] + [tail: bound DECREASES with diameter, klein two-scale]")
