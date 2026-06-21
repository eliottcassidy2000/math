# Test: balanced 2-cluster where each far cluster is a BLOCK of size 2.
# E = B u {f1,f1+1, f2,f2+1}, B bounded, f2/f1 in (1,2.15). Is p0(E)<=cap_k?
# k = |B| + 4. Take B=[0,2,4,6,8] (size5), far blocks at scale -> k=9. cap_9=0.494.
from fractions import Fraction as F
from math import gcd
def sector(fr):
    s=int(fr*7); return 6 if s>=7 else s
def measS7(E):
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for m in range(0,e+1):
            for j in range(7):
                x=F(7*m+j,7*e)
                if 0<=x<=1: bps.add(x)
    bps=sorted(bps); total=F(0)
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2; hit=set(sector(((e*mid)-int(e*mid))+(1 if ((e*mid)-int(e*mid))<0 else 0)) for e in E)
        if len(hit)==7: total+=b-a
    return total

cap9=F(1979,4004)
print("k=9 block-pair far={f1,f1+1,f2,f2+1}, B=[0,2,4,6,8]; cap_9=%.5f"%float(cap9))
worst=(F(0),None)
for f1 in (30,50,80,120):
    for ratio in [(2,1),(3,2),(5,3),(7,4)]:
        p,q=ratio
        f2=f1*p//q
        if f2<=f1+1: continue
        E=sorted(set([0,2,4,6,8,f1,f1+1,f2,f2+1]))
        if len(E)!=9: continue
        g=0
        for x in E[1:]: g=gcd(g,x)
        v=measS7(E)
        prim = (g==1)
        if v>worst[0]: worst=(v,(f1,ratio,prim))
        print(f"  f1={f1} ratio={p}/{q} f2={f2} E={E} prim={prim}: p0={float(v):.5f} <=cap:{v<=cap9}")
print(f"worst p0 = {float(worst[0]):.5f} at {worst[1]}  (cap_9={float(cap9):.5f})")
