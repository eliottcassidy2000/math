#!/usr/bin/env python3
"""
lrc14_netnonzero_kps-S6-wf.py
Characterize WHEN net(E)>0 and bound the max over ALL primitive k=8 E (any spread)
via a structural argument.

KEY IDEA for a rigorous SPREAD BOUND:
  net(E) = meas{x: 8 integer points frac(e_i x) form a 1/7-net of circle}.
  By SCALING invariance, WLOG gcd(E)=1.
  Claim to test: net(E) is determined by E only through its image in a BOUNDED
  quotient.  Specifically the net intervals have endpoints at m/(e_i-e_j); the
  RELEVANT differences for the binding constraint are small.

We test a concrete REDUCTION:
  (R) net(E) <= max over (E mod something) ...
Instead, the practical rigorous plan: bound net(E) by an EXPLICIT formula that we
can maximize in closed form.  Test the candidate:
     net(E) <= sum_{i} 2*(1/7 - 1/d_i)_+    ???
where d_i are consecutive gaps.  Probe relationships between net and gap vector.

Mainly: get the GLOBAL max of net over ALL primitive k=8 (large spread sweep,
including 'spread-out' adversaries like arithmetic-progression-with-hole, Sidon,
B_2 sets) to confirm 44/735 is the true global max with margin, and find the
structural reason net stays small.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(31337)
from math import gcd
from functools import reduce
ONE7=F(1,7)

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
    consec=F(44,735)
    print(f"consec_8 net=44/735={float(consec):.6f}; cap_8=0.381463")
    gmax=-1; gmaxE=None
    # 1) large random primitive sweep
    print("\n[1] random primitive k=8 sweep over many spreads (exact):")
    for spread in [8,9,10,12,15,20,30,50,80,120]:
        best=-1; bestE=None
        for _ in range(4000):
            E=sorted(set([0]+random.sample(range(1,spread+1),7)))
            if len(E)<8: continue
            if reduce(gcd,E)!=1: continue
            nm=netmeas(E)
            if nm>best: best=nm; bestE=E
        if best>gmax: gmax=best; gmaxE=bestE
        print(f"   spread<={spread:3d}: max net={float(best):.6f} ({best})  E={bestE}")
    # 2) structured adversaries
    print("\n[2] structured adversaries (AP-with-hole, Sidon, geometric-ish):")
    adv=[]
    # AP with one hole
    for hole in range(1,8):
        E=[i for i in range(9) if i!=hole]  # 8 of 0..8
        adv.append(("AP-hole%d"%hole,E))
    # Sidon-like / B2 sets
    adv.append(("Sidon",[0,1,3,7,12,20,30,44]))
    adv.append(("Mian-Chowla",[0,1,3,7,12,22,34,46]))
    adv.append(("powers2",[0,1,2,4,8,16,32,64]))
    adv.append(("two-clusters",[0,1,2,3,20,21,22,23]))
    adv.append(("three-clusters",[0,1,10,11,20,21,30,31]))
    for name,E in adv:
        if reduce(gcd,E)!=1:
            g=reduce(gcd,E); E=[e//g for e in E]
        nm=netmeas(E)
        if nm>gmax: gmax=nm; gmaxE=E
        print(f"   {name:16s} E={E}: net={float(nm):.6f} ({nm})")
    print(f"\n[GLOBAL] max net found = {float(gmax):.6f} ({gmax}) at E={gmaxE}")
    print(f"   consec=44/735={float(consec):.6f}; cap_8=0.381463")
    print(f"   {'all <= consec' if gmax<=consec else 'EXCEEDS consec'}; margin to cap_8 = {float(F(2243,5880)-gmax):.6f}")
