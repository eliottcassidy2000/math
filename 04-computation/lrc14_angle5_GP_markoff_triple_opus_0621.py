#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 5 -- MARKOFF TRIPLE on the G_P side (opus 2026-06-21).

The binding configurations: k=11 -> P={1,13}; k=10 -> P={1,12,13}; k=12 -> P={1}.
The closed form D_P = G_P/(P p q) for a 2-element P={p,q} uses three legs
||p||, ||q||, ||pq||. We test the Markoff/continued-fraction structure of the
binding pair directly.

For P={p,q}, meas(G_P) = meas{x: ||px||>=1/14 and ||qx||>=1/14}. The structure is
governed by the THREE-DISTANCE / Stern-Brocot mediant of p/q. The binding (minimal
meas(G_P)) pair should be the one whose p/q has the SIMPLEST continued fraction
(largest 'Markoff-like' uniformity) -- i.e. the AP-like / Fibonacci-extremal pair.

FALSIFIABLE HYPOTHESES:
  H_A: among all pairs {p,q}, 1<=p<q<=13, meas(G_{p,q}) is MINIMIZED exactly at the
       pair whose ratio q/p has the LONGEST/simplest CF (the 'most rational/AP' pair).
  H_B: meas(G_{p,q}) = closed form in (||p/q||-type legs); test the Markoff form
       m = (p^2+q^2+ (q-p)^2) or the mediant denominator p+q controls it.
  H_C: the binding triple {1,12,13} satisfies 1+12=13 (a degenerate Markoff/Frobenius
       'sum' relation -- a Stern-Brocot mediant: 13/1 mediant... ). Test whether ALL
       binding P at every k are 'mediant-closed' (each element = sum of two others or
       extreme of an AP).
"""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def measGP(P):
    bad=[]
    for p in P:
        p=int(p)
        for kk in range(0,p+1):
            lo=max((F(kk)-F(1,14))/p,F(0)); hi=min((F(kk)+F(1,14))/p,F(1))
            if lo<hi: bad.append((lo,hi))
    bad.sort(); tot=F(0); cur=cb=None
    for lo,hi in bad:
        if cur is None: cur,cb=lo,hi
        elif lo<=cb: cb=max(cb,hi)
        else: tot+=cb-cur; cur,cb=lo,hi
    if cur is not None: tot+=cb-cur
    return F(1)-tot

def cf(p,q):
    # continued fraction of p/q
    a=[]
    while q:
        a.append(p//q); p,q=q,p%q
    return a

if __name__=="__main__":
    print("#"*78); print("# ANGLE 5: MARKOFF TRIPLE / CF STRUCTURE OF THE BINDING G_P"); print("#"*78)

    # H_A: minimal meas(G_{p,q}) over all pairs, and its CF
    print("\nH_A: meas(G_{p,q}) over pairs 1<=p<q<=13, sorted ascending (binding first):")
    rows=[]
    for p,q in combinations(range(1,14),2):
        m=measGP([p,q]); rows.append((m,p,q))
    rows.sort()
    for m,p,q in rows[:8]:
        g=gcd(p,q); c=cf(q//g,p//g) if p else None
        print(f"   P={{{p},{q}}}: meas={m}={float(m):.5f}  q/p CF={cf(q,p)}  q-p={q-p}  p+q={p+q}  pq={p*q}")
    print(f"   ...")
    mmin,pm,qm=rows[0]
    print(f"   MINIMAL pair = {{{pm},{qm}}}, meas={mmin}={float(mmin):.5f}")
    print(f"   Does {{1,13}} (the k=11 binder) achieve the global pair-min? {(pm,qm)==(1,13)}")

    # The binding for k=11 was {1,13}. Check: is {1,13} the pair-min, or does the
    # 5-element structure differ? Recompute binders and look at their additive str.
    print("\nH_C: ADDITIVE/MEDIANT structure of binding P at each k:")
    binders={8:(1,5,7,8,9),9:(1,11,12,13),10:(1,12,13),11:(1,13),12:(1,)}
    for k,P in binders.items():
        P=sorted(P)
        # check mediant/AP/sum relations
        rels=[]
        for i in range(len(P)):
            for j in range(len(P)):
                for l in range(len(P)):
                    if i<j and P[i]+P[j]==P[l]: rels.append(f"{P[i]}+{P[j]}={P[l]}")
        # AP detection
        diffs=set(P[i+1]-P[i] for i in range(len(P)-1))
        print(f"   k={k} P={P}: sum-relations={rels}  consec-diffs={sorted(diffs)}  contains 1:{1 in P}  contains 13:{13 in P}")

    # H_B: is meas(G_{1,q}) a clean function of q? (the 1-and-large structure)
    print("\nH_B: meas(G_{1,q}) for q=2..13 (the '1 + large' binding family):")
    for q in range(2,14):
        m=measGP([1,q])
        # candidate closed forms
        print(f"   q={q}: meas={m}={float(m):.5f}  6/7-? {m-F(6,7)}  num*den={m.numerator},{m.denominator}")
