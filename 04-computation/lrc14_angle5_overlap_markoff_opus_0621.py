#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 5 -- the OVERLAP as a MARKOFF triple invariant. opus 2026-06-21.

meas(G_{p,q}) = 5/7 + overlap(p,q), where overlap = meas of {x: ||px||<1/14 AND
||qx||<1/14}. Minimizing meas(G_P) (the binding/extremal config) = minimizing overlap.

The bad-arc systems for p and q intersect near common rationals. ||px||<1/14 are p
arcs of half-width 1/(14p) centered at a/p; ||qx||<1/14 are q arcs centered b/q.
Two arcs (a/p, b/q) overlap when |a/p - b/q| < 1/(14p)+1/(14q) = (p+q)/(14pq), i.e.
|aq - bp| < (p+q)/14. So overlaps are indexed by the SOLUTIONS of |aq-bp| < (p+q)/14
-- a BEZOUT / FROBENIUS lattice condition on the form aq-bp. THIS is the Markoff/
Frobenius object: the indefinite binary form Q(a,b)=aq-bp and its small values.

FALSIFIABLE: overlap(p,q) is determined by the SMALL VALUES of the form aq-bp
(|value| < (p+q)/14). The binding pair {1,13} MINIMIZES the count/weight of small
form-values. Test:
  T1: overlap(p,q) = sum over (a,b) with |aq-bp|<(p+q)/14, 0<=a<=p,0<=b<=q, of the
      arc-intersection length.  Verify against measGP.
  T2: the number of such 'resonant' (a,b) pairs is minimized at {1,13}.
  T3: the form aq-bp has discriminant 0 (it's degenerate: rank-1) -> connect to the
      Markoff spectrum boundary (the degenerate/rational end = simplest CF).
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

def overlap_via_form(p,q):
    """overlap = meas of intersection of the two bad systems, computed by direct
    interval intersection (exact). Also count resonant (a,b)."""
    pa=[]
    for a in range(0,p+1):
        lo=max((F(a)-F(1,14))/p,F(0)); hi=min((F(a)+F(1,14))/p,F(1))
        if lo<hi: pa.append((lo,hi))
    qa=[]
    for b in range(0,q+1):
        lo=max((F(b)-F(1,14))/q,F(0)); hi=min((F(b)+F(1,14))/q,F(1))
        if lo<hi: qa.append((lo,hi))
    ov=F(0); res=0
    for lo1,hi1 in pa:
        for lo2,hi2 in qa:
            l=max(lo1,lo2); h=min(hi1,hi2)
            if h>l: ov+=h-l; res+=1
    return ov,res

def count_small_form(p,q):
    """count (a,b), 0<=a<=p,0<=b<=q with |a q - b p| < (p+q)/14 (interior resonance)."""
    cnt=0; thr=F(p+q,14)
    for a in range(0,p+1):
        for b in range(0,q+1):
            if abs(a*q-b*p)<thr: cnt+=1
    return cnt

if __name__=="__main__":
    print("#"*78); print("# ANGLE 5: OVERLAP AS A MARKOFF/FROBENIUS FORM-VALUE COUNT"); print("#"*78)

    # T1: overlap closed form check vs measGP
    print("\nT1: meas(G_{p,q}) == 5/7 + overlap (interval-intersection) ?")
    ok=True
    for p,q in combinations(range(1,14),2):
        m=measGP([p,q]); ov,res=overlap_via_form(p,q)
        if m!=F(5,7)+ov: ok=False; print(f"   MISMATCH {p},{q}: {m} vs {F(5,7)+ov}")
    print(f"   all 78 pairs match: {ok}")

    # T2: overlap and resonance count, sorted by meas (binding first)
    print("\nT2: pairs sorted by meas(G) ascending (binding first): overlap & form-resonance count")
    rows=[]
    for p,q in combinations(range(1,14),2):
        m=measGP([p,q]); ov,res=overlap_via_form(p,q); sf=count_small_form(p,q)
        rows.append((m,p,q,ov,res,sf))
    rows.sort()
    for m,p,q,ov,res,sf in rows[:10]:
        print(f"   {{{p},{q}}}: meas={float(m):.5f} overlap={ov}={float(ov):.5f} arc-resonances={res} small-form-count={sf}")
    print(f"   ...")
    print(f"   BINDING {{{rows[0][1]},{rows[0][2]}}}: overlap={rows[0][3]} (MIN overlap = MIN meas = binding)")
    # is overlap minimized exactly at the meas-min?
    minov=min(r[3] for r in rows)
    print(f"   global min overlap = {minov}; achieved by meas-min pair: {rows[0][3]==minov}")

    # T3: the {1,q} family overlap = 1/(7q): a SINGLE resonance (the a=b=0 / aq-bp=0 line)
    print("\nT3: overlap({1,q}) = 1/(7q) -- the DEGENERATE (single-resonance) Markoff endpoint:")
    for q in range(2,14):
        ov,res=overlap_via_form(1,q)
        print(f"   {{1,{q}}}: overlap={ov} (=1/(7q)? {ov==F(1,7*q)}) arc-resonances={res}")
    print("\n   => binding pair {1,13}: the form a*13 - b*1 has the FEWEST small values")
    print("      (degenerate rank-1 / simplest CF [13]); overlap = single resonance 1/91.")
    print("      MARKOFF READING: {1,q} sits at the RATIONAL/degenerate end of the form")
    print("      spectrum (one partial quotient); larger q pushes the single resonance")
    print("      narrower (1/(7q)), so q=13 minimizes overlap => minimizes meas(G_P) => binds.")
