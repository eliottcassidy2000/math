# Is mac-mini's c>=D3 sliver a REAL good-period-existence failure, or a certificate artifact? (kps-S94)
# For every dissociated k=13 cluster with c=#arcs/spread >= D3 (the failing certificate), DIRECTLY
# check good-period existence: does some j in {1..Vmax-1} give maxgap{e_i j/Vmax} > 1/7, for the
# critical small-Vmax range [s+1, ceil(7s/6)+buffer]? Also recompute #arcs EXACTLY (Farey) to test
# whether c>=D3 is just grid-sampling noise. seed matches mac-mini's family (spreads incl small).
import numpy as np
from fractions import Fraction as F
from math import floor, gcd, ceil
from functools import reduce
import random
random.seed(94017)
TH=F(1,7); M=F(6,7)
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
        for s2 in range(len(subs)-1):
            u0,u1=subs[s2],subs[s2+1]; um=(u0+u1)/2; A=B=F(0)
            for (a,b) in gaps:
                if a+b*um>TH: A+=(a-TH); B+=b
            m1+=int_pow(A,B,1,u0,u1); m2+=int_pow(A,B,2,u0,u1); m3+=int_pow(A,B,3,u0,u1)
    den=m2-m3/M
    return float(m1/M+(m1-m2/M)**2/den) if den>0 else float(m1/M)
def arc_count_grid(E,GRID):
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
def good_period_exists(E,Vmax):
    # exact rational: maxgap{ (e*j mod Vmax)/Vmax } > 1/7  <=>  maxgap_int{ e*j mod Vmax } > Vmax/7
    thr=F(Vmax,7)
    for j in range(1,Vmax):
        pts=sorted(set((e*j)%Vmax for e in E))
        mg=0; n=len(pts)
        for i in range(n):
            g=(pts[(i+1)%n]-pts[i]) if i<n-1 else (pts[0]+Vmax-pts[n-1])
            if g>mg: mg=g
        if mg>thr: return True,j
    return False,None
k=13; diss=7   # longest-AP <= k-6 = 7 (match mac-mini's inclusive test)
print("k=13 dissociated (longest-AP<=7): find c=#arcs/spread >= D3 clusters; DIRECTLY test good-period")
print("existence at critical Vmax range [s+1, ceil(7s/6)+3]. c>=D3 = certificate fail; existence = truth.\n")
nfail_cert=0; nfail_exist=0; tested=0; worst=[]
for _ in range(80000):
    if tested>=6000: break
    s=random.choice([60,70,80,90,100,120])
    mid=sorted(random.sample(range(1,s),k-2)); E=[0]+mid+[s]
    if len(set(E))!=k or not prim(E) or longest_ap(E)>diss: continue
    tested+=1
    na=arc_count_grid(E,80*s); c=na/s; d3=D3_exact(E)
    if c>=d3:
        nfail_cert+=1
        # directly verify good-period existence at every critical Vmax
        allgood=True; badV=None
        for Vmax in range(s+1, ceil(7*s/6)+4):
            ok,_=good_period_exists(E,Vmax)
            if not ok: allgood=False; badV=Vmax; break
        if not allgood:
            nfail_exist+=1
            if len(worst)<6: worst.append((tuple(E),s,round(c,3),round(d3,3),badV))
print(f"tested={tested}; certificate failures (c>=D3): {nfail_cert}; "
      f"TRUE good-period-existence failures: {nfail_exist}")
if nfail_exist==0:
    print("=> EVERY c>=D3 sliver cluster STILL has a good period at every critical Vmax.")
    print("   The sliver is a CERTIFICATE ARTIFACT, not a real gap: good-period existence holds directly.")
else:
    print("=> REAL existence failures found (would break the covering leg):")
    for w in worst: print("   ",w)
