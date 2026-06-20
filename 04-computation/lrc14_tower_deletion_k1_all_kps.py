from fractions import Fraction
import itertools
from math import gcd
from functools import reduce
def lm(C, theta=Fraction(1,14), comps=False):
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
        return safe, nc
    return safe
def cf(x): return x.numerator//x.denominator+(1 if x.numerator%x.denominator else 0)
thr2=Fraction(426,35035)
print("k=1 layer (delete 2^a + one hole + one tail), comb cutoff + exact residue, ALL a:")
for a in [0,1,2,3]:
    bit=2**a; worstR=0; below=[]; ncores=0
    for h in [d for d in range(1,14) if d!=bit]:
        B=tuple(d for d in range(1,14) if d not in (bit,h))
        M,c=lm(B,comps=True)
        denom=6*M-7*thr2
        assert denom>0, (a,h,M)
        R=cf((2*c)/denom); worstR=max(worstR,R)
        for r in range(14,R):
            if r in B: continue
            C=tuple(sorted(B+(r,)))
            if len(C)!=12: continue
            if reduce(gcd,C)!=1: continue
            ncores+=1
            L=lm(C)
            if L<thr2: below.append((L,h,r))
    st="PROVED" if not below else f"{len(below)} BELOW"
    print(f"  a={a} (del {bit}): comb cutoff R={worstR}, residue cores={ncores}, below_thr2={len(below)}  => {st}")
    for L,h,r in sorted(below)[:5]: print(f"      BELOW meas={L} hole={h} tail={r}")
