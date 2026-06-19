#!/usr/bin/env python3
"""
lrc14_PZ_kps-S6-wf.py
Second-moment / Paley-Zygmund route to bound Q(E)=meas{V<=1/56}, where
   V(x) = sum_t (g_t(x)-1/8)^2 = sum_t g_t(x)^2 - 1/8 >= 0,
and N(E) subset {maxgap<=1/7} subset {sum g^2<=1/7} = {V<=1/56}.

We need an UPPER bound on Q(E)=meas{V<=1/56} (uniform in E) that beats
   cap_8 = 2243/5880 = 0.381463.
Equivalently a LOWER bound on meas{V>1/56} >= 1-cap_8 = thr_8 = 3637/5880=0.618537.

Tools (all give meas{V<=t} upper bounds from moments of V):
  - One-sided Chebyshev (Cantelli):  meas{V<=t} <= Var(V)/(Var(V)+(E[V]-t)^2)  for t<E[V].
  - Paley-Zygmund:  meas{V>θE[V]} >= (1-θ)^2 E[V]^2/E[V^2].
We compute E[V]=int V dx, E[V^2]=int V^2 dx EXACTLY per E (piecewise polynomial:
V is degree-2 in x per cell, V^2 degree-4), then evaluate the bounds.

Goal: see if Cantelli/PZ with the EXACT moments gives meas{V<=1/56} <= cap_8 for
consec_8 and (more importantly) UNIFORMLY.  We compute the moments for many E and
the resulting bound; we look for the WORST (largest) Cantelli bound over E.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(99)
from math import gcd
from functools import reduce
ONE8=F(1,8); T56=F(1,56)

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

def cell_lin(E,a,b):
    """Return list of (s_t,c_t) with gap_t(x)=s_t x + c_t on cell (a,b)."""
    E=sorted(set(E)); n=len(E); mid=(a+b)/2
    order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
    floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
    out=[]
    for t in range(n):
        o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
        s=F(E[o2]-E[o1]); c=F(floors[t]-floors[(t+1)%n]+wrap)
        out.append((s,c))
    return out

def moments_V(E):
    """E[V], E[V^2] exact, where V=sum_t(g_t-1/8)^2.
       Per cell: h_t(x)=g_t-1/8 = s_t x + (c_t-1/8). V=sum h_t^2 = quadratic p(x).
       Integrate p over cell for E[V]; integrate p^2 (quartic) for E[V^2]."""
    EV=F(0); EV2=F(0)
    for a,b in cells(E):
        lins=cell_lin(E,a,b)
        # p(x)=sum (s x + d)^2 with d=c-1/8 -> P2 x^2 + P1 x + P0
        P2=sum(s*s for s,c in lins)
        P1=sum(2*s*(c-ONE8) for s,c in lins)
        P0=sum((c-ONE8)**2 for s,c in lins)
        # int p over [a,b]
        def intpoly(coefs):  # coefs high->low? use explicit
            pass
        # E[V] contribution: int (P2 x^2+P1 x+P0)
        EV += P2*(b**3-a**3)/3 + P1*(b**2-a**2)/2 + P0*(b-a)
        # p^2 = P2^2 x^4 + 2P2P1 x^3 + (P1^2+2P2P0) x^2 + 2P1P0 x + P0^2
        Q4=P2*P2; Q3=2*P2*P1; Q2=P1*P1+2*P2*P0; Q1=2*P1*P0; Q0=P0*P0
        EV2 += Q4*(b**5-a**5)/5 + Q3*(b**4-a**4)/4 + Q2*(b**3-a**3)/3 + Q1*(b**2-a**2)/2 + Q0*(b-a)
    return EV,EV2

def cantelli_bound(EV,EV2,t):
    """meas{V<=t} <= Var/(Var+(EV-t)^2) for t<EV (V>=0).  Returns None if t>=EV."""
    Var=EV2-EV*EV
    if t>=EV: return None
    return Var/(Var+(EV-t)**2)

if __name__=="__main__":
    cap8=F(2243,5880); thr8=F(3637,5880)
    print(f"cap_8={float(cap8):.6f}, thr_8={float(thr8):.6f}, t=1/56={float(T56):.6f}, 1/8={float(ONE8):.6f}")
    # consec_8
    EV,EV2=moments_V(list(range(8)))
    Var=EV2-EV*EV
    cb=cantelli_bound(EV,EV2,T56)
    print(f"\nconsec_8: E[V]={float(EV):.6f} ({EV}), E[V^2]={float(EV2):.6f}, Var={float(Var):.6f}")
    print(f"   Cantelli meas{{V<=1/56}} <= {float(cb):.6f}  vs cap_8={float(cap8):.6f}  "
          f"{'OK (<cap8)' if cb<=cap8 else 'TOO WEAK (>cap8)'}")
    # broad scan: worst (largest) Cantelli bound over primitive E
    print("\nScanning primitive k=8 for WORST (largest) Cantelli bound and min E[V]:")
    worst=F(0); worstE=None; minEV=F(10); minEVE=None
    cnt=0
    for W in range(7,15):
        for body in itertools.combinations(range(1,W+1),7):
            E=(0,)+body
            if E[-1]!=W: continue
            if reduce(gcd,E)!=1: continue
            cnt+=1
            EV,EV2=moments_V(list(E))
            if EV<minEV: minEV=EV; minEVE=E
            cb=cantelli_bound(EV,EV2,T56)
            if cb is None:
                # t>=E[V]: Cantelli inapplicable (E[V] very small) -> bound is trivial(1). FLAG.
                worst=F(1); worstE=E
                continue
            if cb>worst: worst=cb; worstE=E
    print(f"   scanned {cnt} primitive sets (W<=14)")
    print(f"   min E[V] = {float(minEV):.6f} at {minEVE}  (if E[V] dips below t=1/56, Cantelli fails)")
    print(f"   WORST Cantelli bound = {float(worst):.6f} at {worstE}")
    print(f"   cap_8={float(cap8):.6f}  => {'UNIFORM CANTELLI WORKS (worst<=cap8)' if worst<=cap8 else 'CANTELLI INSUFFICIENT for some E'}")
