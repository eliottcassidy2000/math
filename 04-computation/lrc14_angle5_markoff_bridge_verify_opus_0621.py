#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 5 -- BRIDGE & ROBUSTNESS of the Markoff/CF G_P extremality. opus 2026-06-21.

Established (exact): meas(G_{p,q}) = 5/7 + overlap; overlap({1,q})=1/(7q); binding
2-element pair = {1,13} (simplest CF, max denominator). 

Now (1) does the SAME 'simplest-CF / mediant-extremal' principle pick the binding P
at sizes 3,4,5 (the actual thr_k binders)? and (2) is there a GENERAL Markoff-form
closed form meas(G_P) = 1 - |P|/7 + (resonance corrections from pairwise forms)?
And (3) the genuine TEST: does the Markoff-extremal structure of the binder give a
LOWER BOUND certificate for consec (the wall), or only a re-description?
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
    a=[]
    while q: a.append(p//q); p,q=q,p%q
    return a

def cf_complexity(P):
    """sum of partial quotients of all consecutive ratios -- 'CF simplicity' proxy.
    SIMPLER CF (fewer, larger partial quotients near the rational end) <-> binding."""
    P=sorted(P)
    # use ratio of largest to others; proxy: total number of CF terms across pairs
    terms=0
    for p,q in combinations(P,2):
        terms+=len(cf(q,p))
    return terms

if __name__=="__main__":
    print("#"*78); print("# ANGLE 5: MARKOFF/CF EXTREMALITY -- BRIDGE & ROBUSTNESS"); print("#"*78)

    # (1) for each size, the binding P and its CF/additive signature
    print("\n(1) Binding P at each size; CF & additive structure:")
    for sz in range(1,6):
        rows=[]
        for P in combinations(range(1,14),sz):
            rows.append((measGP(list(P)),P))
        rows.sort()
        m,P=rows[0]
        P=sorted(P)
        cont1 = 1 in P; contmax = 13 in P
        # additive (mediant) relations
        rels=[f"{P[i]}+{P[j]}={P[l]}" for i in range(len(P)) for j in range(i+1,len(P)) for l in range(len(P)) if P[i]+P[j]==P[l]]
        cfterms=cf_complexity(P) if sz>=2 else 0
        print(f"   |P|={sz}: binder={P} meas={float(m):.5f}  contains{{1,13}}={cont1 and contmax}  "
              f"mediant-rels={rels}  CFterms={cfterms}")
        # is binder among the FEWEST CFterms?
        if sz>=2:
            allP=[sorted(Q) for Q in combinations(range(1,14),sz)]
            mincf=min(cf_complexity(Q) for Q in allP)
            print(f"        min possible CFterms at this size={mincf}; binder achieves min: {cfterms==mincf}")

    # (2) general additive closed form attempt:
    # meas(G_P) = 1 - |P|/7 + sum_{pairs} overlap(p,q) - (triple overlaps) + ...
    # For {1,q1,q2,...} test if meas = 1 - |P|/7 + sum_{i<j} overlap(pi,pj) (inclusion-excl, pairs only)
    print("\n(2) Additive (pairwise inclusion-exclusion) approximation to meas(G_P):")
    def overlap(p,q):
        return measGP([p,q])-F(5,7)
    for P in [(1,12,13),(1,11,12,13),(1,5,7,8,9),(1,2,3),(1,13)]:
        P=list(P)
        approx=F(1)-F(len(P),7)+sum(overlap(p,q) for p,q in combinations(P,2))
        act=measGP(P)
        print(f"   P={P}: pairwise-IE approx={float(approx):.5f} act={float(act):.5f} exact-match={approx==act} diff={float(act-approx):.6f}")

    # (3) the BRIDGE: thr_k=1-min meas(G_P). consec must satisfy mu_{1/7}(consec)>=thr_k.
    # The Markoff reading says: the binder is the simplest-CF/max-denominator config.
    # Honest scope: this DESCRIBES the binder, giving an exact thr_k value. It does NOT
    # by itself prove consec clears it -- that is the SEPARATE measS7(consec)>=thr_k check.
    print("\n(3) HONEST BRIDGE STATUS:")
    print("   The Markoff/CF reading gives an EXACT, structurally-explained binder")
    print("   (simplest CF, max denominator) and the closed form meas(G_{1,q})=5/7+1/(7q).")
    print("   It explains WHY the binder is {1,..,13}-extremal but is a re-description of")
    print("   thr_k, NOT yet a proof that consec (measS7-max) clears thr_k.")
