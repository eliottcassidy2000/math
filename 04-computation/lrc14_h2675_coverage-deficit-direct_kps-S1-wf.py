#!/usr/bin/env python3
"""
LRC(14) HYP-2675 -- "coverage-deficit-direct" angle (kps-S1-wf).

GOAL: prove  span(E) > B  ==>  p0(E) <= cap_k  for primitive E, 0 in E,
|E|=k in {8..12}, with an EXPLICIT B that glues to the DONE span<=14 finite check.

p0(E) = meas(S7(E)) = meas{ x : every inner sector [j/7,(j+1)/7), j=1..6, is hit
by some frac(e*x), e in E }.  (Sector 0 always hit by e=0.)

THE COVERAGE-DEFICIT IDENTITY (exact, elementary):
  Let A_j(E) = { x in [0,1) : NO e in E has frac(e*x) in sector j }   (= sector j empty).
  Then  S7(E) = [0,1) \ union_{j=1}^6 A_j(E),  so
      p0(E) = 1 - meas( union_j A_j )  <=  1 - max_j meas(A_j(E)).            (DEFICIT-1)
  More generally, by inclusion-exclusion / Bonferroni,
      p0(E) <= 1 - meas(A_a) - meas(A_b) + meas(A_a cap A_b)                  (DEFICIT-2)
  for any pair a,b.  These are the L=7 specialisations of THM-534/THM-535.

STRATEGY for the WIDE regime:
  meas(A_j(E)) = meas{ x : for all e, frac(e*x) not in sector j }
               = meas{ x : for all e, {e*x} in [0,1) \ [j/7,(j+1)/7) }.
  Each runner e forbids a 1/7-measure periodic strip.  For a SINGLE runner e,
  meas{x: {e x} notin sector j} = 6/7 exactly.  We want a LOWER bound on the
  intersection over all e in E that stays bounded BELOW by a positive constant
  as span grows -- the "uncovered sector survives".

  KEY:  p0 small  <=>  some sector survives empty on a large set.

This script:
  (1) Exact engine for p0, p1, and per-sector empty measures meas(A_j).
  (2) Verifies DEFICIT-1 / DEFICIT-2 are valid exact upper bounds.
  (3) Measures, over wide primitive E (k=8..12), the quantity
        D(E) := max_j meas(A_j(E))   (the single-sector deficit)
      and checks  1 - D(E) <= cap_k  (sufficient for p0<=cap via DEFICIT-1),
      i.e.  D(E) >= 1 - cap_k.   Target thresholds 1-cap_k:
        k=8: 3637/5880=0.61854   k=9: 2025/4004=0.50575   k=10: 36/91=0.39560
        k=11: 25/91=0.27473      k=12: 1/7=0.14286
  (4) Finds the cross-over span B above which DEFICIT-1 alone certifies, per k,
      and where it does NOT, escalates to DEFICIT-2 (pair) and the full moment
      bound L_y (THM-534).  Reports the minimal B that glues to span<=14.

Marked PROVED only for the exact identities; the threshold facts are
VERIFIED-numerically (exact rationals over sampled/scanned E) unless an
all-E argument is given.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
import random, sys

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7)}
ONE_MINUS_CAP = {k: 1-CAP[k] for k in CAP}

def primitive(E):
    nz=[abs(x) for x in E if x]
    return bool(nz) and reduce(gcd,nz)==1

def breakpoints(E):
    """Exact wall breakpoints: all j/(7e) for e in E nonzero, plus 0,1."""
    E=[e for e in E if e]
    L=reduce(lambda a,b:a*b//gcd(a,b), E, 1)
    d=7*L
    bps=set([0,d])
    for e in E:
        step=d//(7*e)
        x=0
        for _ in range(7*e+1):
            bps.add(x); x+=step
    return d, sorted(bps)

def full_analysis(E):
    """Return (p0, p1, [meas(A_j) for j=1..6]) all as exact Fractions.
    A_j = sector j empty.  Uses midpoint sampling of each elementary interval."""
    E=sorted(set(E))
    nz=[e for e in E if e]
    d,bps=breakpoints(E)
    den2=2*(d//7)   # = 2L ; midpoint num*e//den2 %7 gives sector of frac(e*mid)
    p0=F(0); p1=F(0)
    empty=[F(0)]*7   # empty[j] = measure where sector j is empty (j=0..6)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        midnum=lo+hi
        hitmask=0
        for e in nz:
            s=(e*midnum//den2)%7
            hitmask |= (1<<s)
        length=F(hi-lo,d)
        # inner-sector missed count
        innermask = hitmask & 0b1111110
        missed = 6 - bin(innermask).count('1')
        if missed==0: p0+=length
        elif missed==1: p1+=length
        for j in range(7):
            if not (hitmask & (1<<j)):
                empty[j]+=length
    return p0,p1,empty

def deficit1(E):
    """1 - max_j meas(A_j); valid UPPER bound on p0 (DEFICIT-1)."""
    p0,p1,empty=full_analysis(E)
    D=max(empty[1:7])  # only inner sectors j=1..6
    return p0, 1-D, D

def pair_empty(E,a,b):
    """meas(A_a cap A_b) exact."""
    E=sorted(set(E)); nz=[e for e in E if e]
    d,bps=breakpoints(E); den2=2*(d//7)
    tot=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        midnum=lo+hi
        hitmask=0
        for e in nz:
            hitmask |= (1<<((e*midnum//den2)%7))
        if not (hitmask&(1<<a)) and not (hitmask&(1<<b)):
            tot+=F(hi-lo,d)
    return tot

def deficit2(E):
    """best pair Bonferroni UB: min over (a,b) of 1-mA-mB+mAB."""
    p0,p1,empty=full_analysis(E)
    best=F(10)
    arg=None
    for a in range(1,7):
        for b in range(a+1,7):
            mab=pair_empty(E,a,b)
            ub=1-empty[a]-empty[b]+mab
            if ub<best:
                best=ub; arg=(a,b)
    return best,arg

# THM-534 moment dual L_y (PROVED upper bound on p0)
def moments(E):
    """factorial moments S_r = sum_{|A|=r} J(A,E), r=0..4, A subset {1..6}.
    S_r = E[C(N,r)] where N=#missed inner sectors."""
    E=sorted(set(E)); nz=[e for e in E if e]
    d,bps=breakpoints(E); den2=2*(d//7)
    # we need distribution of N over elementary intervals
    from math import comb
    S=[F(0)]*5
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        midnum=lo+hi
        hitmask=0
        for e in nz:
            hitmask|=(1<<((e*midnum//den2)%7))
        inner=hitmask&0b1111110
        Nt=6-bin(inner).count('1')
        length=F(hi-lo,d)
        for r in range(5):
            S[r]+=comb(Nt,r)*length
    return S

def Ly(E,k):
    S=moments(E)
    if k in (8,):
        return 1 - S[1] + S[2] - F(9,10)*S[3] + F(3,5)*S[4]
    if k in (9,10):
        return 1 - F(13,18)*S[1] + F(4,9)*S[2] - F(1,6)*S[3]
    # k=11,12
    return 1 - F(1,2)*S[1] + F(1,6)*S[2]

def consec(k): return tuple(range(k))

if __name__=="__main__":
    try: sys.stdout.reconfigure(encoding='utf-8')
    except: pass
    print("="*78)
    print("HYP-2675 coverage-deficit-direct: exact engine + identity verification")
    print("="*78)
    print("Targets 1-cap_k (need single-sector deficit D(E) >= this for DEFICIT-1):")
    for k in sorted(CAP):
        print(f"  k={k}: cap={CAP[k]}={float(CAP[k]):.5f}  1-cap={ONE_MINUS_CAP[k]}={float(ONE_MINUS_CAP[k]):.5f}")

    print("\n--- Sanity: consec_k (the finite-check extremal) ---")
    for k in sorted(CAP):
        E=consec(k)
        p0,ub1,D=deficit1(E)
        ub2,arg=deficit2(E)
        ly=Ly(E,k)
        print(f"k={k} consec: p0={float(p0):.5f}  D(maxemptysector)={float(D):.5f}  "
              f"DEFICIT1_ub={float(ub1):.5f}  DEFICIT2_ub={float(ub2):.5f}({arg})  "
              f"Ly={float(ly):.5f}  cap={float(CAP[k]):.5f}")
        assert p0<=ub1+F(1,10**9), "DEFICIT-1 violated!"
        assert p0<=ub2+F(1,10**9), "DEFICIT-2 violated!"
        assert p0<=ly+F(1,10**9), "Ly violated!"
    print("  [PROVED identities hold: p0 <= DEFICIT1 <= ... and p0<=DEFICIT2, p0<=Ly all OK on consec]")

    # Wide scan: for each k, scan/sample wide primitive E, record D(E), p0, which bound certifies.
    print("\n--- WIDE regime: does single-sector deficit D(E) clear 1-cap_k? ---")
    print("    (if D(E) >= 1-cap_k for all wide E, DEFICIT-1 ALONE proves p0<=cap)")
    rng=random.Random(20260620)
    for k in sorted(CAP):
        tgt=ONE_MINUS_CAP[k]
        worst_p0=(F(0),None)
        worst_margin=(F(10),None)     # min over E of (D(E)-(1-cap))  -- want >=0
        d1_fail=0; d1_total=0
        d2_fail=0; ly_fail=0
        # structured wide families + random wide
        fams=[]
        # clusters
        for M in (15,20,30,50,100,200):
            for csz in ([3,3,2] if k==8 else [3,3,3] if k==9 else None or
                        ([3,3,2,2] if k==10 else [3,3,3,2] if k==11 else [3,3,3,3])):
                pass
        # simpler: random wide
        for _ in range(4000):
            mode=rng.choice(['clusters','spread','dyadic','random'])
            E=[0]
            if mode=='clusters':
                nc=rng.choice([2,3,4]); sc=rng.choice([15,20,30,50,100,200])
                for c in range(nc):
                    for j in range(rng.choice([2,3])):
                        E.append(c*sc+j)
            elif mode=='spread':
                step=rng.choice([2,3,4,5]); E=[i*step for i in range(k)]
            elif mode=='dyadic':
                E=[0,1,2,4,8,16,32,64,128,256,512,1024][:k]
            else:
                while len(set(E))<k: E.append(rng.randint(1, rng.choice([30,60,150,400])))
            E=sorted(set(E))
            while len(E)<k: E.append(E[-1]+rng.randint(1,30))
            E=sorted(set(E))[:k]
            if len(E)!=k or not primitive(E): continue
            span=E[-1]
            if span<=14: continue   # WIDE only
            p0,ub1,D=deficit1(E)
            d1_total+=1
            if p0>worst_p0[0]: worst_p0=(p0,tuple(E))
            marg=D-tgt
            if marg<worst_margin[0]: worst_margin=(marg,tuple(E),float(D),float(p0))
            if ub1>CAP[k]:   # DEFICIT-1 fails to certify
                d1_fail+=1
                ub2,arg=deficit2(E)
                if ub2>CAP[k]:
                    d2_fail+=1
                    ly=Ly(E,k)
                    if ly>CAP[k]:
                        ly_fail+=1
        print(f"k={k}: wide sampled={d1_total}")
        print(f"   worst p0 = {float(worst_p0[0]):.5f}  cap={float(CAP[k]):.5f}  margin={float(CAP[k]-worst_p0[0]):.5f}")
        print(f"   min single-sector deficit margin D-(1-cap) = {float(worst_margin[0]):.5f}  "
              f"(D={worst_margin[2]:.5f}, p0={worst_margin[3]:.5f}) at {worst_margin[1]}")
        print(f"   DEFICIT-1 fails to certify on {d1_fail}/{d1_total} wide; of those DEFICIT-2 fails {d2_fail}; of those Ly fails {ly_fail}")
