from fractions import Fraction
import itertools, sys
from math import gcd
from functools import reduce
def lm(C, comps=False, theta=Fraction(1,14)):
    arcs=[]
    for d in C:
        w=theta/d
        for m in range(0,d+1): arcs.append((Fraction(m,d)-w,Fraction(m,d)+w))
    segs=[]
    for lo,hi in arcs:
        for s in(-1,0,1):
            a2=max(lo+s,Fraction(0)); b2=min(hi+s,Fraction(1))
            if a2<b2: segs.append((a2,b2))
    segs.sort(); cur=Fraction(-1); u=[]
    for a,b in segs:
        if a>cur: u.append([a,b]); cur=b
        elif b>cur: u[-1][1]=b; cur=b
    safe=Fraction(1)-sum(b-a for a,b in u)
    if comps:
        nc=0; prev=Fraction(0)
        for a,b in u:
            if a>prev: nc+=1
            prev=b
        if prev<Fraction(1): nc+=1
        return safe,nc
    return safe
def cf(x): return x.numerator//x.denominator+(1 if x.numerator%x.denominator else 0)
thr2=Fraction(426,35035)
# a=3 k=2: process largest tail r2 via comb on base'=B+r1. Saturating r1 with safety.
bit=8; below=[]; checked=0; maxr1=0
for extra in itertools.combinations([d for d in range(1,14) if d!=bit],2):
    holes=set(extra)|{bit}
    B=tuple(d for d in range(1,14) if d not in holes)
    r1=14
    while r1<=400:
        if r1 in B: r1+=1; continue
        Bp=tuple(sorted(B+(r1,)))
        Mp,cp=lm(Bp,comps=True)
        denom=6*Mp-7*thr2
        if denom<=0:
            print("DENOM<=0", holes, r1, float(Mp)); sys.exit(2)
        Rp=cf((2*cp)/denom)
        if r1>=Rp: break  # residue (r1,Rp) empty for this r1 and all larger
        for r2 in range(r1+1,Rp):
            if r2 in Bp: continue
            C=tuple(sorted(Bp+(r2,)))
            if len(C)!=12: continue
            if reduce(gcd,C)!=1: continue
            checked+=1
            if lm(C)<thr2: below.append((holes,(r1,r2)))
        maxr1=max(maxr1,r1); r1+=1
print(f"a=3 k=2 DONE: maxr1={maxr1} checked={checked} below_thr2={len(below)}")
for h,t in below[:8]: print("  BELOW",h,t)
print("PROVED" if not below else "FAIL")
