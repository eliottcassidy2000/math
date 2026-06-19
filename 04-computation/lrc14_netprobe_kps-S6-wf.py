#!/usr/bin/env python3
"""
lrc14_netprobe_kps-S6-wf.py  -- exploratory probe for the direct net bound.
For the binding consec_8 (and a few hard E), compute N(E)={x: all 8 gaps<=1/7}
EXACTLY as a union of intervals, then for each difference d=e_i-e_j see how N(E)
sits inside {x : ||d x|| in [lo,hi]}.  Goal: find a clean single/few-difference
covering with total measure <= cap_8.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
ONE7=F(1,7)

def net_intervals(E):
    """N(E) = {x in [0,1): every cyclic gap of {frac(e_i x)} <= 1/7}, EXACT."""
    E=sorted(set(E)); n=len(E)
    # breakpoints where the cyclic ORDER of frac(e_i x) changes: x=m/(e_i-e_j)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
    good=[]
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        # On (a,b) order fixed; gap_t(x)=frac(e_{o(t+1)}x)-frac(e_{o(t)}x).
        # gap_t(x) = (e_{o(t+1)}-e_{o(t)}) x - (floor_{t+1}-floor_t), with wrap +1 on last.
        # Condition ALL gaps<=1/7 -> intersection of linear ineqs -> subinterval of (a,b).
        lo=a; hi=b; feasible=True
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(floors[t]-floors[(t+1)%n]+wrap)  # gap = s*x + c
            # want s*x + c <= 1/7
            if s==0:
                if not (c<=ONE7): feasible=False; break
            elif s>0:
                hi=min(hi,(ONE7-c)/s)
            else:
                lo=max(lo,(ONE7-c)/s)
            if lo>=hi: feasible=False; break
        if feasible and lo<hi:
            good.append((lo,hi))
    # merge
    good.sort(); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out

def meas(iv): return sum((b-a for a,b in iv),F(0))

def norm_interval_of_d(N, d):
    """For each x in N, compute ||d x|| and report the min/max range observed
    by sampling cell midpoints (||dx|| is piecewise linear in x)."""
    vals=[]
    for a,b in N:
        # subdivide at points where d x crosses integer/half: x=m/d, x=(m+1/2)/d
        pts=set([a,b])
        # crossings of frac(dx)=0 and =1/2 (where ||dx|| kinks)
        # ||dx|| kinks at dx in Z and dx in Z+1/2
        import math
        m0=math.floor(a*d*2); m1=math.ceil(b*d*2)
        for m in range(m0,m1+1):
            xc=F(m,2*d)
            if a<xc<b: pts.add(xc)
        pts=sorted(pts)
        for u,v in zip(pts,pts[1:]):
            mid=(u+v)/2
            dx=(d*mid)%1
            nx=min(dx,1-dx)
            vals.append(nx)
    return (min(vals),max(vals)) if vals else (None,None)

if __name__=="__main__":
    cap8=F(2243,5880)
    print(f"cap_8 = {cap8} = {float(cap8):.6f}")
    E=list(range(8))
    N=net_intervals(E)
    print(f"\nE=consec_8={E}; meas(N)={meas(N)}={float(meas(N)):.6f}  (cap_8={float(cap8):.4f})")
    print(f"N intervals ({len(N)}):")
    for a,b in N:
        print(f"   [{float(a):.5f},{float(b):.5f}] len={float(b-a):.5f}")
    # For each difference d, what is the range of ||dx|| over N?
    diffs=sorted({E[j]-E[i] for i in range(8) for j in range(i+1,8)})
    print("\nFor each difference d: range of ||dx|| over N(E), and meas{x:||dx|| in that range}:")
    for d in diffs:
        lo,hi=norm_interval_of_d(N,d)
        # meas{x in [0,1): ||dx|| in [lo,hi]} = (number of full periods)*length-per-period
        # ||dx|| in [lo,hi] within one period [0,1/d): length 2*(hi-lo); d periods -> 2(hi-lo)... but careful at edges
        per=2*(hi-lo)  # fraction of [0,1) where ||dx|| in [lo,hi] (since d full periods each contributing 2(hi-lo)/... )
        # Actually meas over [0,1) of {||dx|| in [lo,hi]} = 2*(hi-lo) regardless of d (each period length 1/d, contributes 2(hi-lo)/d... times d = 2(hi-lo))
        print(f"   d={d:2d}: ||dx|| in [{float(lo):.5f},{float(hi):.5f}]  "
              f"=> meas{{||dx|| in range}} = {float(2*(hi-lo)):.5f}  "
              f"{'*** <= cap8 ***' if 2*(hi-lo)<=cap8 else ''}")
