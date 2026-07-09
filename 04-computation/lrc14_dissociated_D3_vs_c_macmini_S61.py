"""
mac-mini-2026-07-09-S61 -- the CLEAN per-cluster a-priori inequality closing the dissociated branch.

ARC-COUNT existence: good period exists iff #arcs < rho*.Vmax. Two a-priori facts:
  rho* >= D3(E)              [THM-661: the degree-3 covering-moment bound, EXACT via Farey]
  #arcs = c.spread <= c.Vmax [spread <= Vmax]
So a good period EXISTS whenever  c := #arcs/spread < D3(E).  (Then #arcs <= c.Vmax < D3.Vmax <= rho*.Vmax.)
BOTH sides are a-priori/exact -- no equidistribution, no resonance sum. This tests the single clean
inequality  #arcs/spread < D3(E)  over dissociated (L<=k-6) spread-dense clusters.
"""
import numpy as np
from fractions import Fraction as F
from math import floor, gcd
from functools import reduce
import random
random.seed(61123)

TH = F(1,7); M = F(6,7)
def int_pow(A,B,p,lo,hi):
    if B==0: return A**p*(hi-lo)
    return ((A+B*hi)**(p+1)-(A+B*lo)**(p+1))/(B*(p+1))
def D3_exact(E):
    E=sorted(E); k=len(E); ds=set(e for e in E if e)
    for i in range(k):
        for j in range(i+1,k): ds.add(E[j]-E[i])
    bps=set([F(0),F(1)])
    for d in ds:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(bps); m1=m2=m3=F(0)
    for c in range(len(bps)-1):
        x0,x1=bps[c],bps[c+1]; xm=(x0+x1)/2
        lin=[(F(-floor(e*xm)),F(e)) for e in E]
        order=sorted(range(k),key=lambda j:lin[j][0]+lin[j][1]*xm); sp=[lin[j] for j in order]
        gaps=[(sp[i+1][0]-sp[i][0],sp[i+1][1]-sp[i][1]) for i in range(k-1)]
        gaps.append((F(1)+sp[0][0]-sp[k-1][0],sp[0][1]-sp[k-1][1]))
        subs=set([x0,x1])
        for (a,b) in gaps:
            if b!=0:
                xs=(TH-a)/b
                if x0<xs<x1: subs.add(xs)
        subs=sorted(subs)
        for s in range(len(subs)-1):
            u0,u1=subs[s],subs[s+1]; um=(u0+u1)/2; A=B=F(0)
            for (a,b) in gaps:
                if a+b*um>TH: A+=(a-TH); B+=b
            m1+=int_pow(A,B,1,u0,u1); m2+=int_pow(A,B,2,u0,u1); m3+=int_pow(A,B,3,u0,u1)
    den=m2-m3/M
    return float(m1/M + (m1-m2/M)**2/den) if den>0 else float(m1/M)
def arc_count(E,GRID):
    x=(np.arange(GRID)+0.5)/GRID; Ea=np.array(sorted(E),float)
    ph=np.mod(np.outer(x,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    good=(g.max(axis=1)>1/7).astype(int); return int(np.sum((good-np.roll(good,1))==1))
def prim(E):
    E=sorted(E); return len(E)>=2 and reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best

print("CLEAN a-priori inequality  c=#arcs/spread < D3(E)  over dissociated (L<=k-6) clusters:\n")
import sys
for k in (11, 13):
    diss=k-6; worst_ratio=(0,None); nfail=0; tested=0; minmargin=9.9
    for _ in range(4000):
        if tested>=250: break
        s=random.choice([80,120,160])
        mid=sorted(random.sample(range(1,s),k-2)); E=[0]+mid+[s]
        if len(set(E))!=k or not prim(E): continue
        if longest_ap(E)>diss: continue     # dissociated only
        tested+=1
        d3=D3_exact(E); na=arc_count(E,80*s); c=na/s
        margin=d3-c
        if c>=d3: nfail+=1
        if margin<minmargin: minmargin=margin
        if c/d3>worst_ratio[0]: worst_ratio=(c/d3,(tuple(E[:5]),s,round(c,3),round(d3,3)))
    print(f"k={k} (dissociated L<={diss}): {tested} clusters; c>=D3 (FAILURES): {nfail}; "
          f"min margin (D3-c) = {minmargin:.4f}")
    print(f"   worst c/D3 ratio = {worst_ratio[0]:.3f} at {worst_ratio[1]}")
    print(f"   => {'c < D3 ALWAYS: dissociated branch CLOSES by the clean a-priori inequality' if nfail==0 else 'FAILS -- needs the tighter bound'}\n")
