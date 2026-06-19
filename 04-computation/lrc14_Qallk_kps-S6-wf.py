#!/usr/bin/env python3
"""
lrc14_Qallk_kps-S6-wf.py
Check the Sum-of-squares RELAXATION at ALL binding k=8..12:
  N(E) subset {x: sum_t g_t^2 <= 1/7}=:Q(E)  (since maxgap>=sum g^2).
For each k, compute Q(consec_k) (EXACT, rational) and compare to cap_k.
If Q(consec_k) < cap_k AND consec maximizes Q (tested broadly), the relaxation
gives a clean sufficient route at that k.  Report margins.

Also: is Q(consec_k) itself below cap_k?  (necessary for the relaxation to help.)
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
from math import gcd
from functools import reduce
ONE7=F(1,7)

# caps from upstream (EXACT)
CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}

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
    E=sorted(set(E)); n=len(E); mid=(a+b)/2
    order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
    floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
    out=[]
    for t in range(n):
        o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
        s=F(E[o2]-E[o1]); c=F(floors[t]-floors[(t+1)%n]+wrap)
        out.append((s,c))
    return out

def Q_bounds(E):
    """rigorous (Qlo,Qhi) for meas{sum g^2 <= 1/7}."""
    from math import isqrt
    def sqrt_bounds(fr):
        if fr<0: return (F(0),F(0))
        scale=10**26; v=int(fr*scale*scale); r=isqrt(v)
        return F(r,scale),F(r+1,scale)
    E=sorted(set(E)); lo_list=[]; hi_list=[]
    for a,b in cells(E):
        lins=cell_lin(E,a,b)
        A2=sum(s*s for s,c in lins); A1=sum(2*s*c for s,c in lins); A0=sum(c*c for s,c in lins)-ONE7
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
            if disc<0: continue
            slo,shi=sqrt_bounds(disc)
            r1_in=(-A1-slo)/(2*A2); r2_in=(-A1+slo)/(2*A2)
            r1_out=(-A1-shi)/(2*A2); r2_out=(-A1+shi)/(2*A2)
            lo1=max(a,r1_in); hi1=min(b,r2_in)
            if lo1<hi1: lo_list.append((lo1,hi1))
            lo2=max(a,r1_out); hi2=min(b,r2_out)
            if lo2<hi2: hi_list.append((lo2,hi2))
    def mtot(L):
        L.sort(); out=[]; t=F(0)
        for x,y in L:
            if out and x<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],y))
            else: out.append((x,y))
        for x,y in out: t+=y-x
        return t
    return mtot(lo_list),mtot(hi_list)

if __name__=="__main__":
    print("SUM-OF-SQUARES RELAXATION at k=8..12:  N(E) subset {sum g^2<=1/7}=Q(E).")
    print(f"{'k':>3} {'Q(consec)_lo':>13} {'Q(consec)_hi':>13} {'cap_k':>10} {'margin(cap-Qhi)':>16} {'relax helps?':>13}")
    for k in range(8,13):
        ql,qh=Q_bounds(list(range(k)))
        cap=CAP[k]; margin=cap-qh
        helps = qh<=cap
        print(f"{k:>3} {float(ql):>13.6f} {float(qh):>13.6f} {float(cap):>10.6f} {float(margin):>16.6f} "
              f"{'YES' if helps else 'NO':>13}")
    print("\nNow: is consec the MAX of Q over primitive E?  (broad exact scan per k, small spread)")
    for k in range(8,13):
        cap=CAP[k]; _,qhc=Q_bounds(list(range(k)))
        best=qhc; bestE=tuple(range(k)); over=False
        Wcap = {8:13,9:12,10:12,11:13,12:14}[k]
        cnt=0
        for body in itertools.combinations(range(1,Wcap+1),k-1):
            E=(0,)+body
            if reduce(gcd,E)!=1: continue
            cnt+=1
            _,qh=Q_bounds(list(E))
            if qh>best+F(1,10**20):
                best=qh; bestE=E; over=True
        print(f"  k={k:2d}: scanned {cnt} primitive (spread<={Wcap}); max Qhi={float(best):.6f} at {bestE}; "
              f"Q(consec)hi={float(qhc):.6f}; {'consec is max' if not over else '*** consec NOT max ***'}; "
              f"max Qhi {'<=cap' if best<=cap else '> cap'}")
