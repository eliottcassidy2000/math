#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C (mac-mini) -- CORRECTED finding: the CDF difference is SINGLE-CROSSING (+then-).

HONEST REVERSAL: full stochastic dominance N_consec <=_st N_E is FALSE (it fails
already for the nearest neighbor E=[0,1,2,3,4,5,6,8]: F_consec-F_E =
[+233/1960, -3/140, -43/5880, -11/840, -1/280, -1/392, 0,0]).  consec WINS at t=0
(measS7) but LOSES the middle CDF.  This matches HYP-2635 (FOSD eliminated).

What IS true (the real structure): the CDF difference
   Phi_E(t) = F_consec(t) - F_E(t) = P(N_consec<=t) - P(N_E<=t)
is SINGLE-CROSSING in t: it starts POSITIVE at t=0 (=measS7 margin>0) and then is
<=0 for all t>=1 (one sign change, + then -).  Equivalently:
   consec has MORE mass on {N=0} and LESS on every up-set {N>=1,...} cumulatively
   beyond the first cell -- a SINGLE-CUT reshaping toward N=0.

This is the EXACT probabilistic content of "measS7 is a CUT/occupancy extremality"
(HYP-2722c): consec pushes the empty-count law toward the COVERED end (N=0) by a
single CDF crossing, NOT by uniform st-dominance.

CONSEQUENCE (the lever):  Phi_E single-crossing (+ then -) with Phi_E(7)=0 is
EXACTLY the statement that  E[psi(N_consec)] <= E[psi(N_E)]  for every function psi
that is INCREASING ... no.  Precisely: if Phi_E(t)>=0 for t<c and <=0 for t>=c
(single down-crossing), then for any psi with a single sign change in its forward
difference matching... -> the clean consequence is:
   sum_t Phi_E(t)*( -Delta psi(t) ) ... ; we TEST which test-functions psi are
ordered.  In particular psi=1[N=0] gives measS7 (the crux).  We test:
  (T0) Phi_E(0)>=0  (crux)         -- expect HOLD
  (T1) single down-crossing of Phi_E (sign pattern +...+,-...-, one change) -- the HYP
  (T2) the IMPLIED order: consec MAXIMIZES E[psi(N)] for every psi that is
       NON-INCREASING and CONVEX-toward-0... we just test the concrete family
       psi_a(N)=1[N=0] is the only sharply-tested one; also test E[a^N] decreasing a.
We report sign-pattern of Phi_E for ALL competitors and the count of multi-crossings
(which would refute the clean single-cut picture).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce

def occupancy_pi(E):
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi
def N_law(E):
    pi=occupancy_pi(E); nl=[F(0)]*8
    for h in range(8): nl[7-h]+=pi[h]
    return nl
def cdf(nl):
    F_=[F(0)]*8; acc=F(0)
    for t in range(8): acc+=nl[t]; F_[t]=acc
    return F_
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def sign_pattern(seq):
    nz=[(1 if x>0 else -1) for x in seq if x!=0]
    sc=sum(1 for a,b in zip(nz,nz[1:]) if a!=b)
    return sc,nz

def run(k,W):
    C=consec(k); FC=cdf(N_law(C))
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E) and list(E)!=C]
    t0_fail=0; multi=[]; not_plusfirst=[]; onecross=0
    for E in bank:
        FE=cdf(N_law(list(E)))
        Phi=[FC[t]-FE[t] for t in range(7)]  # t=7 is 0 trivially, drop
        if Phi[0]<0: t0_fail+=1
        sc,nz=sign_pattern(Phi)
        if sc<=1:
            onecross+=1
            if nz and nz[0]<0: not_plusfirst.append((list(E),nz))
        else:
            multi.append((list(E),nz,[str(x) for x in Phi]))
    return dict(k=k,W=W,n=len(bank),t0_fail=t0_fail,onecross=onecross,
                multi=multi,not_plusfirst=not_plusfirst)

if __name__=="__main__":
    print("="*78)
    print("THREAD C CORRECTED: CDF difference Phi_E(t)=F_consec(t)-F_E(t) single-crossing?")
    print("  HYP-C3: Phi_E starts >=0 (t=0, the crux) and has at most ONE sign change (+then-).")
    print("="*78)
    for k,W in [(8,13),(9,13),(10,12),(8,15),(9,15),(10,15),(8,17)]:
        r=run(k,W)
        print(f"\n--- k={k}, span<={W}, competitors={r['n']} ---")
        print(f"  Phi(0)>=0 (crux measS7) holds on all but: {r['t0_fail']}")
        print(f"  single-crossing (<=1 sign change): {r['onecross']}  / {r['n']}")
        print(f"  MULTI-crossing (refutes HYP-C3): {len(r['multi'])}")
        if r['multi']:
            for E,nz,phi in r['multi'][:5]:
                print(f"      E={E}  signs={nz}  Phi={phi}")
        print(f"  single-cross but starts NEGATIVE: {len(r['not_plusfirst'])}")
    print("\nDONE.")
