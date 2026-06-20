#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEPENDENT adversarial verification of HYP-2661-generalized.
Own exact-rational lonely-measure engine (no import of kps scripts).
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd

THETA = Fraction(1,14)
THR1 = Fraction(7,858)      # drop-6 mouth
THR2 = Fraction(426,35035)  # AP one-hole second value

def lonely_measure(C, theta=THETA):
    """meas{ t in [0,1): ||d t|| > theta for all d in C }, exact rational.

    For speed d, the BAD set (||d t|| <= theta) is union over integers m of
    [m/d - theta/d, m/d + theta/d]. ||x|| = distance to nearest integer.
    We build the union of all bad arcs on the circle [0,1), then complement.
    """
    bad = []
    for d in C:
        if d == 0:
            continue
        w = theta / d
        # m ranges so that m/d in [0,1): m = 0..d (m=d gives 1.0 == 0.0)
        for m in range(0, d+1):
            lo = Fraction(m,d) - w
            hi = Fraction(m,d) + w
            # wrap into [0,1): add translates -1,0,+1 clipped to [0,1)
            for shift in (-1,0,1):
                a = lo+shift; b = hi+shift
                a2 = max(a, Fraction(0)); b2 = min(b, Fraction(1))
                if a2 < b2:
                    bad.append((a2,b2))
    bad.sort()
    cur = Fraction(-1)
    union = []
    for a,b in bad:
        if a > cur:
            union.append([a,b]); cur = b
        elif b > cur:
            union[-1][1] = b; cur = b
    covered = sum(b-a for a,b in union)
    return Fraction(1) - covered

def carry_profile(C):
    d = {}
    for e in C:
        if e == 0: continue
        m=e; a=0
        while m%2==0: m//=2; a+=1
        d[m] = d.get(m,0)+2**a
    return dict(sorted(d.items()))

def has_tower(C):
    s=set(C)
    return {1,2,4,8} <= s

def primitive(C):
    return reduce(gcd, C) == 1

# ---- sanity: reproduce the four known exact values ----
def ap_tail_core(holes, tails):
    return tuple(sorted([d for d in range(1,14) if d not in holes]+list(tails)))

if __name__ == "__main__":
    print("== KNOWN-VALUE SANITY (independent engine) ==")
    checks = {
      "drop-6":   (ap_tail_core({6},[]),   Fraction(7,858)),
      "{5:+2} exc":(ap_tail_core({6,10},[20]), Fraction(3859,420420)),
      "drop-12":  (ap_tail_core({12},[]),  Fraction(426,35035)),
      "two-tail min":(ap_tail_core({4,6,10},[20,46]), Fraction(50189,3223220)),
    }
    ok=True
    for name,(C,exp) in checks.items():
        L = lonely_measure(C)
        match = (L==exp)
        ok = ok and match
        print(f"  {name:13s}: meas={L} expected={exp}  MATCH={match}  ({float(L):.7f})")
    print(f"  THR1=7/858={float(THR1):.7f}  THR2=426/35035={float(THR2):.7f}")
    print(f"  ALL KNOWN VALUES MATCH: {ok}")
    print(f"  THR1<THR2? {THR1<THR2}  (mouth below AP-2nd)")

# ---- DECISIVE EXHAUSTIVE CENSUS ----
def census(cap):
    """All primitive 12-subsets of [1,cap]; find sub-THR2 cores and tower status."""
    below = []
    tower_broken_min = None
    total = 0; primtot = 0
    for C in itertools.combinations(range(1,cap+1), 12):
        total += 1
        if not primitive(C):
            continue
        primtot += 1
        L = lonely_measure(C)
        if L < THR2:
            below.append((L, C, has_tower(C)))
        if not has_tower(C):
            if tower_broken_min is None or L < tower_broken_min[0]:
                tower_broken_min = (L, C)
    return total, primtot, below, tower_broken_min

if __name__ == "__main__" and len(sys.argv)>1 and sys.argv[1]=="census":
    print("== EXHAUSTIVE CENSUS: all primitive 12-subsets of [1,cap] ==")
    overall_counterexample = None
    for cap in range(13,21):
        total, primtot, below, tbmin = census(cap)
        below.sort()
        print(f"\n cap={cap}: {total} 12-subsets, {primtot} primitive.")
        print(f"   sub-THR2 cores: {len(below)}")
        for L,C,tw in below:
            tag = "TOWER" if tw else "*** TOWER-BROKEN (COUNTEREXAMPLE) ***"
            print(f"     meas={L}={float(L):.6f}  C={C}  {tag}")
            if not tw:
                overall_counterexample = (cap,L,C)
        if tbmin:
            L,C = tbmin
            print(f"   min over TOWER-BROKEN primitive cores: meas={L}={float(L):.6f}")
            print(f"     at C={C}  >=THR2? {L>=THR2}  margin={float(L-THR2):.6f}")
    print("\n== VERDICT ==")
    if overall_counterexample:
        print(f"  COUNTEREXAMPLE FOUND: {overall_counterexample}")
    else:
        print("  NO tower-broken core dips below THR2 anywhere in [1,20].")
        print("  Every sub-THR2 primitive core contains the full tower {1,2,4,8}.")
