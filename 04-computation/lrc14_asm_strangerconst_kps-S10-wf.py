#!/usr/bin/env python3
"""
Probe the claim's STEP4 'stranger sup L_y' constants 0.25704/0.48729/0.56241.
The pure 1-stranger consec_(k-1)+N converged to ~0.213/0.447/0.523 (lower).
Search which stranger family realizes (or upper-bounds) 0.257/0.487/0.562:
  (a) consec_(k-1) + N over ALL N (already ~0.213) -> recheck max over wider N grid
  (b) ANY (k-1)-subset-of-small + 1 far stranger:  base in bounded box, satellite far
  (c) 2-stranger: bounded core + 2 far satellites
  (d) the sup over the WHOLE wide region (max L_y found in the big hunt was
      0.237/0.459/0.535 -- compare).
Report the max L_y each family reaches and whether 0.257/0.487/0.562 is hit/exceeded/over-cap.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd
from functools import reduce

SECTORS=[(F(j,7),F(j+1,7)) for j in range(1,7)]
def fast_profile(E):
    E=sorted(set(E)); bp=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(1,7):
            for end in (F(j,7),F(j+1,7)):
                m=0
                while True:
                    xv=(end+m)/e
                    if xv>=1:break
                    if xv>=0:bp.add(xv)
                    m+=1
    bp=sorted(b for b in bp if 0<=b<1); Sacc=[F(0)]*5
    nz=[e for e in E if e!=0]
    for lo,hi in zip(bp,bp[1:]+[F(1)]):
        if hi<=lo:continue
        L=hi-lo; mid=(lo+hi)/2; fr=[(e*mid)%1 for e in nz]; free=0
        for (a,b) in SECTORS:
            if not any(a<v<b for v in fr): free+=1
        for r in range(5):
            if free>=r: Sacc[r]+=L*comb(free,r)
    return Sacc
DUAL={8:([F(1),F(-1),F(1),F(-9,10),F(3,5)],4),
      9:([F(1),F(-13,18),F(4,9),F(-1,6)],3),
      10:([F(1),F(-13,18),F(4,9),F(-1,6)],3)}
def Ly(E,k):
    y,R=DUAL[k]; S=fast_profile(E); return sum(y[r]*S[r] for r in range(R+1))
CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91)}
CLAIM={8:0.25704,9:0.48729,10:0.56241}
def prim(E): return reduce(gcd,[e for e in E if e>0])==1

for k in (8,9,10):
    best=(F(0),None)
    # (b) base = any (k-2)-subset of {1..k+2} (bounded), plus the 0, plus 1 far stranger N
    bases=[]
    for sub in combinations(range(1,k+3), k-2):
        E0=[0]+list(sub)
        bases.append(E0)
    Ns=[40,80,160,320,640]
    for E0 in bases:
        for N in Ns:
            E=E0+[N]
            if len(set(E))!=k: continue
            if not prim(E): continue
            ly=Ly(E,k)
            if ly>best[0]: best=(ly,E)
    # (c) 2-stranger: consec_(k-2) + two far
    base=list(range(k-1))
    for N1 in [40,80,160]:
        for N2 in [N1+40,N1+80,2*N1,3*N1]:
            E=base[:k-2]+[N1,N2]
            E=sorted(set(E))
            if len(E)!=k or not prim(E): continue
            ly=Ly(E,k)
            if ly>best[0]: best=(ly,E)
    bf=float(best[0])
    print(f"k={k}: max L_y over bounded-core+far-stranger family = {bf:.5f} ({best[0]})  E={best[1]}")
    print(f"      claim sup {CLAIM[k]}  cap={float(CAP[k]):.5f}  "
          f"reached>=claim? {bf>=CLAIM[k]-1e-4}  claim<cap? {CLAIM[k]<float(CAP[k])}  realized<cap? {bf<float(CAP[k])}")
