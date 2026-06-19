#!/usr/bin/env python3
"""
lrc14_netscale_kps-S6-wf.py
Establish & verify EXACT invariances of meas(N(E)) that underpin a rigorous
'consecutive maximizes net' argument.

Invariance 1 (SCALING):  meas(N(c*E)) = meas(N(E)) for any positive integer c
    with gcd(c, .)=... actually for ANY positive integer c (since x->cx mod 1 is
    a measure-preserving c-to-1 cover and the gap structure of {frac(c e_i x)}
    over x in [0,1) equals that of {frac(e_i y)} over y in [0,1)).
    => WLOG gcd(E)=1, AND constant gap vectors all reduce to consecutive.

Invariance 2 (TRANSLATION): N(E)=N(E - min E); WLOG 0 in E (already assumed).

Invariance 3 (REFLECTION): meas(N(E)) = meas(N(maxE - E))  (x -> -x = 1-x sym).

We VERIFY these exactly on many E, and we verify the crucial reduction:
   meas(N(E)) is invariant under E -> c*E, so the family of relevant E is the
   PRIMITIVE sets (gcd=1).  Combined with M3 (consecutive maximizes among small-gap
   vectors), the remaining task is a SPREAD BOUND: show large-spread primitive E
   have net < consec.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(7)
ONE7=F(1,7)
from math import gcd
from functools import reduce

def net_intervals(E):
    E=sorted(set(E)); n=len(E)
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
        lo=a; hi=b; feasible=True
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(floors[t]-floors[(t+1)%n]+wrap)
            if s==0:
                if not (c<=ONE7): feasible=False; break
            elif s>0: hi=min(hi,(ONE7-c)/s)
            else: lo=max(lo,(ONE7-c)/s)
            if lo>=hi: feasible=False; break
        if feasible and lo<hi: good.append((lo,hi))
    good.sort(); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(iv): return sum((b-a for a,b in iv),F(0))
def netmeas(E): return meas(net_intervals(E))

if __name__=="__main__":
    print("[INV1] SCALING: meas(N(c*E)) == meas(N(E))?")
    bad=0; tot=0
    for trial in range(40):
        k=random.choice([8,9])
        E=sorted(set([0]+random.sample(range(1,12),k-1)))
        if len(E)<k: continue
        base=netmeas(E)
        for c in (2,3,5):
            Ec=[c*e for e in E]
            v=netmeas(Ec); tot+=1
            if v!=base:
                bad+=1
                print(f"   VIOLATION E={E} c={c}: {float(base):.6f} vs {float(v):.6f}")
    print(f"   scaling checks: {tot}, violations: {bad}  => {'INVARIANT (PROVED-pattern)' if bad==0 else 'NOT invariant'}")

    print("\n[INV3] REFLECTION: meas(N(E)) == meas(N(maxE - E))?")
    bad=0; tot=0
    for trial in range(60):
        k=random.choice([8,9,10])
        E=sorted(set([0]+random.sample(range(1,14),k-1)))
        if len(E)<k: continue
        Er=sorted(max(E)-e for e in E)
        tot+=1
        if netmeas(E)!=netmeas(Er):
            bad+=1; print(f"   VIOLATION E={E}")
    print(f"   reflection checks: {tot}, violations: {bad}  => {'INVARIANT' if bad==0 else 'NOT'}")

    print("\n[REDUCTION] So WLOG gcd(E)=1. The 'constant gap' family all collapses to consec.")
    print("   Remaining: bound net for PRIMITIVE (gcd=1) E of large spread.")
    # quick: among PRIMITIVE k=8 sets with max(E)<=14, what is the max net?
    print("\n[PRIMITIVE k=8, max(E)<=14] scanning for max net (exact)...")
    best=-1; bestE=None
    cnt=0
    for body in itertools.combinations(range(1,15),7):
        E=(0,)+body
        if reduce(gcd,E)!=1: continue
        cnt+=1
        nm=netmeas(list(E))
        if nm>best: best=nm; bestE=E
    print(f"   scanned {cnt} primitive sets; MAX net = {best} = {float(best):.6f} at E={bestE}")
    print(f"   consec_8 net = {F(44,735)} = {float(F(44,735)):.6f}")
    print(f"   cap_8 = {F(2243,5880)} = {float(F(2243,5880)):.6f}")
