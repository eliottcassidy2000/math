from fractions import Fraction as F
from math import gcd
def norm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return min(r,1-r)
def minnorm(fam,t): return min(norm(F(v)*t) for v in fam)
def Pmax_set(P,level):  # all rational t (denom<=2maxP) with min_P ||.|| >= level, with the max value
    pts=[]; best=F(0)
    for q in range(2,2*max(P)+2):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            val=minnorm(P,F(a,q))
            if val>best: best=val
            if val>=level: pts.append((F(a,q),val))
    return best, pts

def cluster_lift(P,K):
    S=sorted(P+K); pmax=max(P)
    muP,_ = Pmax_set(P, F(1))  # M(P) exactly (best)
    # use P-good points at level 1/14 (bigger freedom), sorted by P-value desc
    _,goodP = Pmax_set(P, F(1,14))
    goodP.sort(key=lambda z:-z[1])
    margin = muP - F(1,14)
    if margin<=0: return None,muP
    # for each good t0, search s in [-margin/pmax, margin/pmax] on a fine rational grid; slide the cluster
    for t0,vP in goodP[:20]:
        # s must keep core safe: |p*s|<=vP-1/14 for all p, i.e. |s|<=(vP-1/14)/pmax
        smax=(vP-F(1,14))/pmax
        if smax<=0: continue
        # sample s to slide the cluster across the circle (rate ~ k0); step so k0*step ~ 1/200
        k0=min(K); steps=400
        for i in range(-steps,steps+1):
            s=smax*F(i,steps)
            t=t0+s
            if minnorm(S,t)>=F(1,14):
                return (t0,vP,t,minnorm(S,t)), muP
    return None,muP

battery=[
 ("P={1..10}",list(range(1,11)),[257,258,263]),
 ("P={1..10}",list(range(1,11)),[300,301,302]),          # very tight cluster spread 2
 ("P={1..10}",list(range(1,11)),[500,502,505,509]),      # need |P|=9 -> adjust
 ("P={1..9}",list(range(1,10)),[500,502,505,509]),
 ("P={1..11}",list(range(1,12)),[257,258]),
 ("P={1..8,11,12}",[1,2,3,4,5,6,7,8,11,12],[257,258,263]),
 ("P={1..10}",list(range(1,11)),[157,158,159]),          # slower cluster (near 13*maxP=130)
 ("P={1..10}",list(range(1,11)),[131,132,133]),          # barely fast (>130)
]
print("cluster-lift witness search (proves M(S)>=1/14 if a witness is found):")
for name,P,K in battery:
    if len(P)+len(K)!=13:
        print(f"  {name}+{K}: wrong size skip"); continue
    res,muP=cluster_lift(P,K)
    if res:
        t0,vP,t,mn=res
        print(f"  {name}+{K}: WITNESS t={t} minnorm={mn}={float(mn):.4f}>=1/14  [M(P)={float(muP):.4f}, t0={t0}] PROVED M(S)>=1/14")
    else:
        print(f"  {name}+{K}: no witness in construction (M(P)={float(muP):.4f}) -- outside regime")
