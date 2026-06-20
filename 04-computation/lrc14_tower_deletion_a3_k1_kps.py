from fractions import Fraction
import itertools
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
thr2=Fraction(426,35035)

# STRATEGY (mirror THM-543, single-tail layer first, a=3):
#  k=1: C = ({1..13}\{8,h}) U {r}, r>=14.  base B={1..13}\{8,h}, M_B, c_B exact (h != 8).
#  comb: meas(C) >= (6/7) M_B - 2 c_B/(7 r).  This >= thr2 iff  r >= R_h := ceil( 2 c_B / (7*((6/7)M_B - thr2)) )
#        = ceil( 2 c_B / (6 M_B - 7 thr2) ).  Valid only when 6 M_B - 7 thr2 > 0 (always, M_B>=0.045>>thr2).
#  Then finite residue: r in [14, R_h) exact-scan.  Show all >= thr2.
print("=== a=3, k=1 layer (delete 8 + one more hole + one tail) ===")
worstR=0
residue_fail=[]
allbelow=[]
for h in [d for d in range(1,14) if d!=8]:
    B=tuple(d for d in range(1,14) if d not in (8,h))
    M,c=lm(B,comps=True)
    denom=6*M-7*thr2  # >0
    Rh=( (2*c) / denom )
    import math
    Rh_ceil=Rh.numerator//Rh.denominator + (1 if Rh.numerator%Rh.denominator else 0)
    worstR=max(worstR,Rh_ceil)
    # residue scan r in [14, Rh_ceil)
    for r in range(14, Rh_ceil):
        if r in B: continue
        C=tuple(sorted(B+(r,)))
        if len(C)!=12: continue
        from math import gcd; from functools import reduce
        if reduce(gcd,C)!=1: continue
        L=lm(C)
        if L<thr2: allbelow.append((L,h,r))
print(f"  worst comb cutoff R over all h: {worstR}")
print(f"  finite-residue rows below thr2: {len(allbelow)}")
for L,h,r in sorted(allbelow):
    print(f"    BELOW: meas={L}={float(L):.6f}  hole={h} tail={r}")
if not allbelow:
    print("    => k=1 layer for a=3 fully certified >= thr2 (comb cutoff + exact residue).")
