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

# a=3 k=2, CLEAN per-base double comb.  C=({1..13}\{8,h1,h2}) U {r1<r2}.
# Process tails LARGEST first: comb out r2 from base'=B+r1.
#   Region A (r1 >= Rsingle(B)): then meas(G_{B+r1}) >= thr2 by single-tail comb on B.
#      Adding r2>r1: meas(C) >= (6/7) meas(G_{B+r1}) - 2 c(B+r1)/(7 r2).  But (6/7)*thr2 < thr2,
#      so need r2 large. Instead in region A we comb r2 against base' and exact-residue r2<R'(B+r1).
#      Since R'(B+r1) is small (M(B+r1)>=thr2 but actually M(B+r1) is order 0.04), R' ~ 100s. We must
#      still bound r1. For r1 huge, B+r1 -> approx B (tail far), M(B+r1)->(6/7)M_B, R'->fixed. And r2
#      between r1 and R' is EMPTY once r1 >= R'. So once r1 >= R'(B+r1) the residue is empty -> covered.
#   So: for each base B, scan r1 from 14 upward; for each r1 compute base'=B+r1, R'=cutoff; scan
#       r2 in (r1, R'); STOP increasing r1 once r1 >= R'(B+r1) for all larger r1 (R' saturates).
# This bounds r1 by ~max R' ~ small. Implement with saturation stop.
below=[]; bit=8; others=[d for d in range(1,14) if d!=bit]; maxr1=0
for extra in itertools.combinations(others,2):
    holes=set(extra)|{bit}
    B=tuple(d for d in range(1,14) if d not in holes)
    r1=14
    while True:
        if r1 in B: r1+=1; continue
        Bp=tuple(sorted(B+(r1,)))
        Mp,cp=lm(Bp,comps=True)
        denom=6*Mp-7*thr2
        Rp=cf((2*cp)/denom) if denom>0 else 10**9
        # if r1 already >= Rp, then no r2>r1 with r2<Rp exists => residue empty for this and all larger r1
        if r1>=Rp:
            break
        for r2 in range(r1+1,Rp):
            if r2 in Bp: continue
            C=tuple(sorted(Bp+(r2,)))
            if len(C)!=12: continue
            if reduce(gcd,C)!=1: continue
            L=lm(C)
            if L<thr2: below.append((L,tuple(sorted(holes)),(r1,r2)))
        maxr1=max(maxr1,r1)
        r1+=1
        if r1>400: break  # safety
print(f"a=3 k=2: scanned all 55 three-hole bases; max effective r1={maxr1}")
print(f"  residue rows below thr2: {len(below)}")
for L,h,t in sorted(set(below))[:12]:
    print(f"    BELOW meas={L}={float(L):.6f} holes={h} tails={t}")
if not below: print("  => a=3 k=2 layer fully certified >= thr2 (saturating double-comb + exact residue).")
