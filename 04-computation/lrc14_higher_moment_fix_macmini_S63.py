"""
mac-mini-2026-07-09-S63 -- FIX route (c) after MISTAKE-128: D3=B_3 is too weak a mu-lower-bound
for 7-structured co-offsets (c/D3=1.40). But B_d -> mu as d grows. Does a HIGHER-degree
covering-moment bound B_d (d=4,5,6) exceed c, restoring the a-priori certificate  c < B_d <= mu?

B_d(E) = max{ sum_i c_i E[W^i] : sum_i c_i w^i <= 1_{w>0} on [0,6/7] }  (the degree-d moment LP,
a rigorous lower bound on mu = P(W>0)). Test on the MISTAKE-128 counterexample and a zoo of
7-structured / dissociated hard clusters: is c < B_d for d<=6?
"""
import numpy as np
from fractions import Fraction as F
from math import floor, gcd
from functools import reduce
from scipy.optimize import linprog

TH = F(1,7); Mf = 6/7
def W_moments(E, maxp=6):
    E=sorted(E); k=len(E); ds=set(e for e in E if e)
    for i in range(k):
        for j in range(i+1,k): ds.add(E[j]-E[i])
    bps=set([F(0),F(1)])
    for d in ds:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(bps); mom=[F(0)]*(maxp+1)
    def ip(A,B,p,lo,hi):
        if B==0: return A**p*(hi-lo)
        return ((A+B*hi)**(p+1)-(A+B*lo)**(p+1))/(B*(p+1))
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
            for p in range(1,maxp+1): mom[p]+=ip(A,B,p,u0,u1)
    return [float(x) for x in mom]
def mu_exact(E):
    # mu = P(W>0) = meas(good) via the same Farey cells
    return None  # use float below
def mu_float(E, GRID=400000):
    x=(np.arange(GRID)+0.5)/GRID; Ea=np.array(sorted(E),float)
    ph=np.mod(np.outer(x,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return float((g.max(axis=1)>1/7).mean())
_w=np.linspace(1e-4,Mf,1200);
def Bd(mom,d):
    Aub=np.column_stack([_w**i for i in range(1,d+1)])
    res=linprog([-mom[i] for i in range(1,d+1)],A_ub=Aub,b_ub=np.ones(len(_w)),
                bounds=[(None,None)]*d,method='highs')
    return -res.fun if res.success else None
def arc_c(E,GRID=600000):
    x=(np.arange(GRID)+0.5)/GRID; Ea=np.array(sorted(E),float)
    ph=np.mod(np.outer(x,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    good=(g.max(axis=1)>1/7).astype(int); return int(np.sum((good-np.roll(good,1))==1))/(max(E)-min(E))

E0=[0,7,14,21,26,29,37,44,51,58,67,75,82]  # MISTAKE-128 counterexample
mom=W_moments(E0); c=arc_c(E0); mu=mu_float(E0)
print(f"MISTAKE-128 set E={E0}")
print(f"  c=#arcs/spread={c:.4f}, mu={mu:.4f}")
for d in (3,4,5,6):
    b=Bd(mom,d); print(f"  B_{d}={b:.4f}   c<B_{d}? {'YES (certificate restored!)' if c<b else 'no'}")
