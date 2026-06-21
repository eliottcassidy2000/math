#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 7 follow-up: the RAW magnitude vector is the WRONG majorization variable
(test S20 part C: measS7 not Schur-monotone in |e|).  Huffer-Shepp's variable is
ARC LENGTHS, not speeds.  Find the RIGHT variable.

Two candidate "Huffer-Shepp arc-length" analogues for the LRC survival problem:

  CANDIDATE 1 -- the SURVIVAL-WIDTH VECTOR itself W=(W_1..W_6).
    HYP-2757: Z/7* permutes the phases a->ca; consec is the dilation-symmetric
    config => its W vector is the most BALANCED (orbit-averaged).  If sum_a W_a is
    being maximized subject to a fixed "budget", and the per-phase widths trade
    off, a Schur-CONVEX constraint on W with balanced=max would explain consec.
    TEST: is consec's W vector majorization-MINIMAL (most equal) among the
    stratum, AND does measS7 (=sum) correlate with balance?  [sum is invariant
    under majorization, so this alone can't drive sum-max -- check.]

  CANDIDATE 2 -- the COVERED-ARC LENGTHS at the binding instant (the real
    Huffer-Shepp transfer): at each phase a, the 7 clocks paint 7 arcs of the
    cell as x sweeps; W_a = the sub-measure where the UNION covers all sectors.
    The Huffer-Shepp theorem is about the lengths of the COVERING arcs.  In LRC
    the covering "arc" contributed by clock e in cell a has half-width
    proportional to 1/(7|e|) (drift speed 7e => time to cross a sector).
    So the HS arc-length analogue is the vector L = (1/(7|e|))_e, i.e. the
    RECIPROCAL magnitudes.  consec = {0..7} => L ~ (inf,1,1/2,...,1/7): VERY
    spread out (the fixed clock e=0 has infinite arc).  Huffer-Shepp:
    spread-out arc lengths => MORE coverage.  THIS MATCHES consec-max!
    TEST: is measS7 Schur-CONVEX in the reciprocal-magnitude (arc-length) vector?
    HYP-HS2: more spread in L=(1/|e|) => more measS7; consec is L-maximal-spread.
"""
import itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

def occupancy_pi(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi
def measS7(E): return occupancy_pi(E)[7]
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def majorizes(x,y):
    xs=sorted(x,reverse=True); ys=sorted(y,reverse=True)
    if len(xs)!=len(ys) or sum(xs)!=sum(ys): return False
    px=py=F(0)
    for i in range(len(xs)-1):
        px+=xs[i]; py+=ys[i]
        if px<py: return False
    return True

def W_a(E,a):
    import math
    E=sorted(set(E)); lo=F(a,7)-F(1,14); hi=F(a,7)+F(1,14); bps={lo,hi}
    for e in E:
        if e==0: continue
        d=7*abs(e); jlo=lo*d; jhi=hi*d
        j0=math.floor(jlo); j1=math.ceil(jhi)
        for j in range(j0-1,j1+2):
            x=F(j,d)
            if lo<=x<=hi: bps.add(x)
    bps=sorted(bps); tot=F(0)
    for l,h in zip(bps,bps[1:]):
        if h<=l: continue
        xm=(l+h)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=h-l
    return tot

if __name__=="__main__":
    print("="*78)
    print("CANDIDATE 2: Schur-CONVEXITY in the RECIPROCAL-magnitude (arc-length) vector")
    print("  L = (1/|e|)_{e in E, e!=0}.  HS predicts: more SPREAD in L => more measS7.")
    print("  consec={0..7}: L=(1,1/2,...,1/7) + the e=0 fixed clock (infinite arc).")
    print("="*78)
    # To keep finite, drop the e=0 clock (always present in full-residue) and
    # majorize the reciprocal vector of the NONZERO speeds.
    W=14
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),7)]
    bank=[E for E in bank if primitive(E) and set(e%7 for e in E)=={0,1,2,3,4,5,6}]
    print(f"  full-residue shapes span<= {W}: {len(bank)}")
    cache={}
    def m_(E):
        if E not in cache: cache[E]=measS7(list(E))
        return cache[E]
    # group by sum of reciprocals (the L-budget) -- but reciprocals rarely share sum.
    # Instead: test ALL pairs; whenever L(Ea) majorizes L(Eb) (Ea more spread, same
    # sum of reciprocals), check sign.  Few will share sum -> also do an APPROX
    # version: count, over ALL comparable-by-weak-majorization pairs, the sign.
    convL=0; concL=0; compL=0
    Ls={E:sorted((F(1,abs(e)) for e in E if e!=0),reverse=True) for E in bank}
    for Ea,Eb in itertools.combinations(bank,2):
        la,lb=Ls[Ea],Ls[Eb]
        if sum(la)!=sum(lb): continue
        if majorizes(la,lb):
            compL+=1; d=m_(Ea)-m_(Eb)
            if d>F(1,10**15): convL+=1
            if d<-F(1,10**15): concL+=1
        elif majorizes(lb,la):
            compL+=1; d=m_(Eb)-m_(Ea)
            if d>F(1,10**15): convL+=1
            if d<-F(1,10**15): concL+=1
    print(f"  same-reciprocal-sum majorization-comparable pairs = {compL}")
    print(f"  Schur-CONVEX  evidence (more L-spread beat balanced) = {convL}")
    print(f"  Schur-CONCAVE evidence (balanced beat more L-spread) = {concL}")

    print()
    print("="*78)
    print("CANDIDATE 1: is consec's W=(W_1..W_6) vector majorization-MINIMAL (balanced)?")
    print("  Z/7* permutes phases (HYP-2757) so consec's W is a single orbit-balanced")
    print("  vector.  Test whether consec's W is majorized by EVERY same-sum competitor")
    print("  (i.e. consec is the most-equal W => Schur-concave symmetric fn of W max).")
    print("="*78)
    C=consec(8); WC=sorted((W_a(C,a) for a in range(1,7)),reverse=True); sC=sum(WC)
    print(f"  consec W (sorted) = {[str(w) for w in WC]}  sum={sC}={float(sC):.6f}")
    # find competitors with the SAME sum (= same measS7) and compare W-balance
    same_sum=[E for E in bank if m_(E)==sC and list(E)!=C]
    print(f"  competitors with identical measS7 = {len(same_sum)}")
    # KEY realization: among ALL shapes, does HIGHER measS7 come with MORE-BALANCED W?
    # measure W-balance by variance (lower = more balanced). Correlate with measS7.
    import statistics
    rows=[]
    for E in bank:
        wv=[float(W_a(E,a)) for a in range(1,7)]
        rows.append((float(m_(E)), statistics.pvariance(wv), E))
    rows.sort(reverse=True)
    print("  TOP 8 by measS7 (measS7, W-variance, shape):")
    for r in rows[:8]:
        print(f"     measS7={r[0]:.5f}  Wvar={r[1]:.3e}  {r[2]}")
    print("  BOTTOM 4 by measS7:")
    for r in rows[-4:]:
        print(f"     measS7={r[0]:.5f}  Wvar={r[1]:.3e}  {r[2]}")
    # correlation
    ms=[r[0] for r in rows]; vs=[r[1] for r in rows]
    try:
        import math as _m
        n=len(ms); mm=sum(ms)/n; mv=sum(vs)/n
        cov=sum((a-mm)*(b-mv) for a,b in zip(ms,vs))/n
        sd_m=_m.sqrt(sum((a-mm)**2 for a in ms)/n); sd_v=_m.sqrt(sum((b-mv)**2 for b in vs)/n)
        corr=cov/(sd_m*sd_v) if sd_m*sd_v else 0
        print(f"  Pearson corr(measS7, W-variance) = {corr:.4f}")
        print("  (HYP CAND-1: NEGATIVE corr => higher measS7 <-> more balanced W)")
    except Exception as e:
        print("  corr failed",e)
