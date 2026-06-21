#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C (mac-mini) -- the SINGLE-CROSSING / cut structure of the N-law.

Finding so far:
 * measS7(E) = P(N_E=0) = pgf_E(0); consec maximizes it (VERIFIED k=8,9,10).
 * pgf order has NO universal crossover z* (z* shrinks with span) -> the
   "consec maxes E[z^N] on a fixed [0,z*]" certificate is REFUTED.
 * BUT ~half the competitors are NEVER caught (consec wins pgf on all (0,1]).

NEW STRUCTURAL HYP (THREAD C deliverable). Let
   delta_E(n) = P_consec(N=n) - P_E(N=n),  n=0..7   (sum_n delta=0).
HYP-C2 (SINGLE-CROSSING / sign-regular cut):  for EVERY competitor E (binding rows),
the sequence delta_E(n) changes sign AT MOST ONCE as n goes 0->7, and when it does,
it goes +...+ then -...- ... NO: anchored at delta_E(0)=p0_consec-p0_E >= 0, the
pattern is  (>=0 block at small n)  then  (<=0 block at large n)  with exactly one
sign change.  Equivalently the consec N-law and the competitor N-law are ordered in
the SINGLE-CROSSING (cut) order: consec puts MORE mass on small N (the covered end)
and LESS on large N (the uncovered end), crossing once.

WHY THIS IS THE RIGHT OBJECT:
 * single-crossing of pmfs at the LEFT (n=0) end  ==>  P(N=0) larger for consec,
   AND  P(N<=t) larger for consec for all t below the crossing  (a partial
   stochastic dominance of the LEFT tail).
 * it is the EXACT probabilistic content of "measS7 is a cut/occupancy extremality"
   (HYP-2722c): the empty-count law is reshaped by a single cut toward 0.
 * a single sign change in delta is the hypothesis of the VARIATION-DIMINISHING /
   sign-regular (TP2 / monotone-likelihood-ratio) order.

TEST (exact Fractions): for all binding competitors, count sign changes of
delta_E(n) over n=0..7, and verify the +...+ / -...- (one-crossing, left-positive)
shape.  Report any multi-crossing competitors (these would refute C2 and be the
true obstruction).  Also test the implied LEFT-TAIL dominance:
   F_consec(t) = P(N<=t) >= F_E(t)   for all t <= crossing point.
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

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def sign_changes(seq):
    """count sign changes ignoring zeros; also return the sign pattern (list of +/-)."""
    nz=[(1 if x>0 else -1) for x in seq if x!=0]
    sc=sum(1 for a,b in zip(nz,nz[1:]) if a!=b)
    return sc, nz

def run(k,W):
    C=consec(k); nlC=N_law(C)
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E) and list(E)!=C]
    multi=[]            # >1 sign change
    not_leftpos=[]      # crossing but not +...- (i.e. starts negative)
    left_tail_fail=[]   # F_consec(t)<F_E(t) for some t <= crossing
    onecross=0; zerocross=0
    p0fail=0
    for E in bank:
        nl=N_law(list(E))
        delta=[nlC[n]-nl[n] for n in range(8)]
        if delta[0]<0: p0fail+=1
        sc,nz=sign_changes(delta)
        if sc==0:
            zerocross+=1
        elif sc==1:
            onecross+=1
            if not (nz[0]>0):     # must start positive (left = covered end heavier for consec)
                not_leftpos.append((list(E),nz))
            # left-tail dominance up to crossing index
            # find crossing index c = first n where cumulative sign flips
            # F_consec(t)-F_E(t) = sum_{n<=t} delta(n); must be >=0 up to crossing
            cum=F(0); ok=True; cumlist=[]
            for n in range(8):
                cum+=delta[n]; cumlist.append(cum)
            # single crossing in delta, left-positive => cum is unimodal, rises then falls to 0;
            # left-tail dominance F_consec>=F_E for ALL t holds iff cum>=0 for all t.
            if any(c<0 for c in cumlist):
                left_tail_fail.append((list(E),[str(c) for c in cumlist]))
        else:
            multi.append((list(E),nz,[str(d) for d in delta]))
    return dict(k=k,W=W,n=len(bank),zerocross=zerocross,onecross=onecross,
                multi=multi,not_leftpos=not_leftpos,left_tail_fail=left_tail_fail,p0fail=p0fail)

if __name__=="__main__":
    print("="*78)
    print("THREAD C HYP-C2: single-crossing (cut) structure of delta_E(n)=P_consec(N=n)-P_E(N=n)")
    print("="*78)
    for k,W in [(8,13),(9,13),(10,12),(8,15),(9,15),(10,15),(8,17)]:
        res=run(k,W)
        print(f"\n--- k={k}, span<={W}, competitors={res['n']} ---")
        print(f"  P(N=0) fails: {res['p0fail']}")
        print(f"  0 sign changes (consec dominates pmf-wise trivially? no, sum=0): {res['zerocross']}")
        print(f"  EXACTLY 1 sign change (single crossing): {res['onecross']}")
        print(f"  >1 sign change (MULTI-crossing, refutes C2): {len(res['multi'])}")
        if res['multi']:
            for E,nz,d in res['multi'][:4]:
                print(f"      E={E}  signs={nz}  delta={d}")
        print(f"  1-crossing but NOT left-positive (+...-): {len(res['not_leftpos'])}")
        print(f"  LEFT-TAIL dominance F_consec(t)>=F_E(t) FAILS: {len(res['left_tail_fail'])}")
        if res['left_tail_fail']:
            for E,cl in res['left_tail_fail'][:4]:
                print(f"      E={E}  cum(F_consec-F_E)={cl}")
    print("\nDONE.")
