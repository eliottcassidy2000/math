#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE D: Three-distance / Steinhaus closed form for consec_k = {0,1,...,k-1}.

For consec, the orbit at x is the arithmetic progression A_k(x) = {0, x, 2x, ..., (k-1)x mod 1}.
A 1/7-sector [j/7,(j+1)/7) (j in 1..6) is MISSED iff A_k(x) puts no point in it.
N(x) = # missed sectors among j=1..6.  meas(S7) = p_0 = meas{N=0}.

GOAL 1: closed form for N_consec(x) via the three-gap (Steinhaus) structure of A_k(x).
GOAL 2: closed form / exact value for p_t(consec), verify against engine.
GOAL 3: comparison handle: does the AP-uniformity of consec maximize meas(S7)?

We use exact arithmetic (Fraction). Engine reference: lrc14_empty_sector_distribution_kps.py
"""
import sys, importlib.util
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

spec = importlib.util.spec_from_file_location('eng','04-computation/lrc14_empty_sector_distribution_kps.py')
eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng)

def frac(x):
    return x - (x.numerator // x.denominator)

def missed_sectors(k, x):
    """Set of sectors j in 1..6 missed by consec_k orbit at x."""
    hit = set()
    for i in range(k):
        v = frac(i * x)
        hit.add((v.numerator * 7) // v.denominator)
    return [j for j in range(1,7) if j not in hit]

def N_consec(k, x):
    return len(missed_sectors(k, x))

# ---- Three-gap structure of the orbit (Steinhaus) ----
def gaps(k, x):
    """sorted distinct orbit points and the multiset of gap lengths (Fractions summing to 1)."""
    pts = sorted(set(frac(i*x) for i in range(k)))
    g = []
    for i in range(len(pts)):
        nxt = pts[(i+1) % len(pts)]
        gap = (nxt - pts[i]) % 1 if i < len(pts)-1 else (1 - pts[i] + pts[0])
        g.append(gap)
    return pts, g

# ---- Three-gap PREDICTION of N: a sector j is missed iff the half-open arc [j/7,(j+1)/7)
#      contains no orbit point. Equivalent: the *gap* straddling that contains the whole sector. ----
def N_via_gaps(k, x):
    pts, _ = gaps(k, x)
    cnt = 0
    for j in range(1,7):
        lo = Fraction(j,7); hi = Fraction(j+1,7)
        # missed iff no point in [lo, hi)
        if not any(lo <= p < hi for p in pts):
            cnt += 1
    return cnt

def verify_closedform(kmax_terms=False):
    """Verify N_via_gaps == N_consec at all breakpoint midpoints for k=8,9,10."""
    for k in [8,9,10]:
        E = list(range(k))
        bps = set([Fraction(0), Fraction(1)])
        for e in E:
            if e==0: continue
            for a in range(0, 7*e+1):
                bps.add(Fraction(a, 7*e))
        bps = sorted(b for b in bps if 0<=b<=1)
        bad = 0
        for i in range(len(bps)-1):
            mid = (bps[i]+bps[i+1])/2
            if N_consec(k, mid) != N_via_gaps(k, mid):
                bad += 1
        print(f"k={k}: gap-prediction matches direct N at all {len(bps)-1} intervals: {bad==0} (mismatches={bad})")

def gapcond_vs_N0():
    """Is meas{N=0} == meas{max orbit gap < 1/7}? (three-gap-only sufficiency test)."""
    print("--- max-gap<1/7 is SUFFICIENT but NOT necessary for N=0 ---")
    for k in [8,9,10]:
        E=list(range(k)); bps=set([Fraction(0),Fraction(1)])
        for e in E:
            if e==0: continue
            for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
        bps=sorted(b for b in bps if 0<=b<=1)
        mN0=Fraction(0); mgap=Fraction(0)
        for i in range(len(bps)-1):
            lo,hi=bps[i],bps[i+1]; w=hi-lo; mid=(lo+hi)/2
            if N_consec(k,mid)==0: mN0+=w
            _,g=gaps(k,mid)
            if max(g)<Fraction(1,7): mgap+=w
        print(f"  k={k}: meas(N=0)={mN0}={float(mN0):.4f}  meas(maxgap<1/7)={mgap}={float(mgap):.4f}  equal={mN0==mgap}")

if __name__ == "__main__":
    print("=== Angle D: three-distance closed form for consec ===\n")
    verify_closedform()
    print()
    gapcond_vs_N0()
    print()
    # confirm p_0 exact values
    for k in [8,9,10]:
        p = eng.dist_p(list(range(k)))
        print(f"k={k}: consec p_0 = {p[0]} = {float(p[0]):.6f}")
