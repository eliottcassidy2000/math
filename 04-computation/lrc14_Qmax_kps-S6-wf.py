#!/usr/bin/env python3
"""
lrc14_Qmax_kps-S6-wf.py
Maximize Q(E)=meas{x: sum_t g_t(x)^2 <= 1/7} over E (0 in E,|E|=8), to test the
candidate uniform proof:  N(E) subset Q(E) and  Q(E) <= Q(consec_8) < cap_8.

Q's per-cell condition is QUADRATIC; roots may be irrational, so Q(E) may be
irrational.  We compute a RIGOROUS UPPER BOUND on Q(E) by, in each cell, taking the
sub-interval where A2 x^2 + A1 x + A0 <= 0 (A0 includes -1/7) and OVER-counting via
an outer rational interval (rounding root outward).  For the MAX search we want an
upper bound on Q; we use exact rational roots when disc is a perfect square and a
*rigorous outer* rational enclosure otherwise (round roots so the interval grows).

Also report Q exactly for consec_8 (it IS rational there) and scan spreads.
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
from math import gcd, isqrt
from functools import reduce
ONE7=F(1,7)

def cells(E):
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        yield a,b

def sqrt_lower(fr):
    """rational r <= sqrt(fr), tight (within 1e-30)."""
    if fr<0: return None
    scale=10**30
    v=int(fr*scale*scale)
    r=isqrt(v)
    return F(r, scale*scale//isqrt(scale*scale)) if False else F(r,scale)  # r/scale <= sqrt(fr)? check
def sqrt_bounds(fr):
    """return (lo,hi) rationals with lo<=sqrt(fr)<=hi, hi-lo<=1e-25."""
    if fr<0: return (F(0),F(0))
    scale=10**26
    v=int(fr*scale*scale)   # = fr*scale^2
    r=isqrt(v)              # r <= sqrt(fr)*scale < r+1
    lo=F(r,scale); hi=F(r+1,scale)
    return lo,hi

def Q_bounds(E):
    """Return (Qlo, Qhi): rigorous rational lower & upper bounds for Q(E)."""
    E=sorted(set(E)); n=len(E); lo_list=[]; hi_list=[]
    for a,b in cells(E):
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        S=[]; C=[]
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=F(E[o2]-E[o1]); c=F(floors[t]-floors[(t+1)%n]+wrap)
            S.append(s); C.append(c)
        A2=sum(s*s for s in S); A1=sum(2*s*c for s,c in zip(S,C)); A0=sum(c*c for c in C)-ONE7
        # region in [a,b] where A2 x^2+A1 x+A0<=0
        if A2==0:
            if A1==0:
                if A0<=0: lo_list.append((a,b)); hi_list.append((a,b))
            elif A1>0:
                xr=-A0/A1; hi=min(b,xr)
                if a<hi: lo_list.append((a,hi)); hi_list.append((a,hi))
            else:
                xr=-A0/A1; loo=max(a,xr)
                if loo<b: lo_list.append((loo,b)); hi_list.append((loo,b))
        else:
            disc=A1*A1-4*A2*A0
            if disc<0:
                continue
            slo,shi=sqrt_bounds(disc)
            # roots r1=(-A1-sqrt)/(2A2), r2=(-A1+sqrt)/(2A2), A2>0 => region [r1,r2]
            # inner (lower bound on Q): use sqrt small => [r1_inner=(-A1-slo)/2A2 (larger), r2_inner=(-A1+slo)/2A2 (smaller)]
            r1_in=(-A1+ (- (slo)) )/(2*A2)  # = (-A1 - slo)/(2A2)
            r1_in=(-A1 - slo)/(2*A2); r2_in=(-A1 + slo)/(2*A2)
            r1_out=(-A1 - shi)/(2*A2); r2_out=(-A1 + shi)/(2*A2)
            # inner enclosure (subset): [max(a,r1_in), min(b,r2_in)]
            lo1=max(a,r1_in); hi1=min(b,r2_in)
            if lo1<hi1: lo_list.append((lo1,hi1))
            # outer enclosure (superset): [max(a,r1_out), min(b,r2_out)]
            lo2=max(a,r1_out); hi2=min(b,r2_out)
            if lo2<hi2: hi_list.append((lo2,hi2))
    def mtot(L):
        L.sort(); out=[]; t=F(0)
        for a,b in L:
            if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
            else: out.append((a,b))
        for a,b in out: t+=b-a
        return t
    return mtot(lo_list), mtot(hi_list)

if __name__=="__main__":
    cap8=F(2243,5880)
    print(f"cap_8={cap8}={float(cap8):.6f}")
    qlo,qhi=Q_bounds(list(range(8)))
    print(f"\nQ(consec_8): in [{float(qlo):.6f}, {float(qhi):.6f}]  (rational since equispaced?)")
    print(f"   cap_8 - Q(consec_8)_hi = {float(cap8-qhi):.6f}  (need >0 for headroom)")

    print("\nMax Q (UPPER bound Qhi) over primitive k=8 by spread W=max(E):")
    gmax=F(0); gmaxE=None
    for W in range(7,17):
        best=F(0); bestE=None
        cnt=0
        for body in itertools.combinations(range(1,W+1),7):
            E=(0,)+body
            if E[-1]!=W: continue
            if reduce(gcd,E)!=1: continue
            cnt+=1
            _,qh=Q_bounds(list(E))
            if qh>best: best=qh; bestE=E
        if best>gmax: gmax=best; gmaxE=bestE
        print(f"   W={W:2d}: max Qhi={float(best):.6f} ({cnt} sets)  E={bestE}  {'<=cap8' if best<=cap8 else '*** >cap8 ***'}")
    print(f"\n  GLOBAL max Qhi (W<=16) = {float(gmax):.6f} at {gmaxE}")
    print(f"  cap_8={float(cap8):.6f}; margin = {float(cap8-gmax):.6f}  "
          f"{'Q-RELAXATION VIABLE so far' if gmax<=cap8 else 'Q EXCEEDS cap8 -> relaxation lossy'}")
