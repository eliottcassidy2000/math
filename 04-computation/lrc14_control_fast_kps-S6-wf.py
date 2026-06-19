"""
Fast control + focused exhaustive k=8 min-mu_{1/7}.
(1) 2/7 control: random search beats consec for mu_{2/7} at k=9,10  (search is potent).
(2) Exhaustive k=8 spread<=12: min mu_{1/7}, argmin, # below thr_8.
"""
from fractions import Fraction as F
from itertools import combinations
from math import floor, gcd
from functools import reduce
import random

def mu_theta(E,theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[floor(E[order[t]]*mid) for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t];o2=order[(t+1)%n];k1=ks[t];k2=ks[(t+1)%n];wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1];c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s)
                if lo<b: subs.append((lo,b))
            else:
                hi=min(b,(theta-c)/s)
                if a<hi: subs.append((a,hi))
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

def prim(E):
    E=sorted(set(E));
    if len(E)<2: return E
    g=reduce(gcd,[e-E[0] for e in E[1:]])
    return [(e-E[0])//g for e in E] if g>1 else E

# (1) 2/7 control
print("=== CONTROL: mu_{2/7}, does random search beat consec? ===")
th27=F(2,7)
random.seed(123)
for k in [9,10]:
    consec=mu_theta(list(range(k)),th27)
    best=consec; bestE=tuple(range(k))
    for _ in range(8000):
        spread=random.randint(k-1,35); pts={0}
        while len(pts)<k: pts.add(random.randint(1,spread))
        E=prim(sorted(pts))
        if len(E)!=k: continue
        m=mu_theta(E,th27)
        if m<best: best=m;bestE=tuple(E)
    print(f"k={k}: consec_2/7={float(consec):.4f} best={float(best):.4f} at {bestE} BEATEN:{best<consec}")

# (2) exhaustive k=8 min mu_{1/7}
print("\n=== Exhaustive k=8 spread<=12 min mu_{1/7} ===")
th=F(1,7); thr8=F(3637,5880)
best=None;bestE=None;below=0;cnt=0
for combo in combinations(range(1,13),7):
    E=(0,)+combo; cnt+=1
    m=mu_theta(list(E),th)
    if best is None or m<best: best=m;bestE=E
    if m<thr8: below+=1
print(f"sets={cnt} min mu={best}={float(best):.4f} at {bestE} consec-argmin:{bestE==tuple(range(8))} below-thr8:{below}")
print("consec mu == min:", best==F(691,735))
