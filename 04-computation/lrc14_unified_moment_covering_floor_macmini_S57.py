"""mac-mini-S57 (THM-661): the UNIFIED covering-moment density floor.
mu_{1/7}(E)=P(W>0), W=uncovered measure (THM-657). Degree-d one-sided moment bound
B_d = max{sum c_i E[W^i] : sum c_i w^i <= 1_{w>0} on [0,6/7]} is a rigorous diameter-free
lower bound on mu. Degree 4 clears k=8,9,10 (block 0.761/0.645/0.553 >= bars 0.675/0.562/0.452);
degree 2 (=PZ, THM-660) clears k=11,12,13. Exact block moments via Farey-cell integration."""
from fractions import Fraction as F
from math import floor, comb
from itertools import combinations
import numpy as np
from scipy.optimize import linprog
TH=F(1,7)
def block_moments(k,maxp=4):
    E=list(range(k)); ds=set(e for e in E if e)|set(E[j]-E[i] for i in range(k) for j in range(i+1,k))
    bps=set([F(0),F(1)])
    for d in ds:
        for m in range(d+1): bps.add(F(m,d))
    bps=sorted(bps); M=[F(0)]*(maxp+1)
    for c in range(len(bps)-1):
        x0,x1=bps[c],bps[c+1];xm=(x0+x1)/2
        lin=[(F(-floor(e*xm)),F(e)) for e in E]
        sp=[lin[j] for j in sorted(range(k),key=lambda j:lin[j][0]+lin[j][1]*xm)]
        gaps=[(sp[i+1][0]-sp[i][0],sp[i+1][1]-sp[i][1]) for i in range(k-1)]+[(F(1)+sp[0][0]-sp[k-1][0],sp[0][1]-sp[k-1][1])]
        subs=set([x0,x1])
        for (a,b) in gaps:
            if b!=0 and x0<(TH-a)/b<x1: subs.add((TH-a)/b)
        subs=sorted(subs)
        for s in range(len(subs)-1):
            u0,u1=subs[s],subs[s+1];um=(u0+u1)/2;A=F(0);B=F(0)
            for (a,b) in gaps:
                if a+b*um>TH: A+=a-TH;B+=b
            for p in range(maxp+1):
                M[p]+=sum(comb(p,r)*A**(p-r)*B**r*(u1**(r+1)-u0**(r+1))/(r+1) for r in range(p+1))
    return M
mP=F(14249,252252)
def GPm(P):
    bad=[(max(F(0),F(j,p)-F(1,14*p)),min(F(1),F(j,p)+F(1,14*p))) for p in P for j in range(p+1)]
    bad=[(l,h) for l,h in bad if h>l];bad.sort();mm=[]
    for l,h in bad:
        if mm and l<=mm[-1][1]: mm[-1]=(mm[-1][0],max(mm[-1][1],h))
        else: mm.append((l,h))
    return 1-sum(h-l for l,h in mm)
def Bd(M,N,Wmax=F(6,7),ng=2000):
    ts=np.linspace(1e-7,float(Wmax),ng)
    A=[[w**i for i in range(N+1)] for w in ts]+[[1.0]+[0.0]*N]; b=[1.0]*ng+[0.0]
    r=linprog([-float(M[i]) for i in range(N+1)],A_ub=np.array(A),b_ub=np.array(b),bounds=[(None,None)]*(N+1),method='highs')
    c=[F(ci).limit_denominator(10**6) for ci in r.x]
    return sum(c[i]*M[i] for i in range(N+1)), max(sum(float(c[i])*w**i for i in range(N+1)) for w in ts)
print("k : bar     B_2(PZ)  B_4     verdict")
for k in (8,9,10,11,12,13):
    M=block_moments(k,4); bar=mP+1-min(GPm(P) for P in combinations(range(1,14),13-k)) if k<13 else mP
    b2,_=Bd(M,2); b4,f4=Bd(M,4)
    print(f"{k} : {float(bar):.4f}  {float(b2):.4f}  {float(b4):.4f}  "
          f"{'deg2 clears' if b2>=bar else ('deg4 clears' if b4>=bar else 'SHORT')} (feas {f4:.3f})")
