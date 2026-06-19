#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL re-verification of the single-far-element Weyl-limit reduction.
Independent reimplementation (exact Fractions).

Checks:
  (A) meas_S7 reimplemented from scratch, cross-checked.
  (B) Weyl limit formula LIM(B0)=meas(B0)+(1/7)*P(exactly-1-missed) verified exactly
      against meas(B0 u {w}) for large w, k=8 (AP base) AND a WIDE non-AP base.
  (C) max_{B0} LIM(B0) over ALL primitive (k-1)-bases with 0, bounded by a window.
      Compare to M_k. Is M_k really the bounded global max? Search for a bounded
      set (AP or not) that beats consec.
  (D) Hunt: does ANY single-far set B0 u {w} (w large) exceed M_k? counterexample search.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def missed(E, x):
    hit = set()
    for e in E:
        v = e*x
        v = v - (v.numerator//v.denominator)
        hit.add((v.numerator*7)//v.denominator)
    return [j for j in range(1,7) if j not in hit]

def dist(E):
    """exact measure of {x in [0,1): exactly t of sectors 1..6 missed}, t=0..6."""
    E = sorted(set(E))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e + 1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        t = len(missed(E, mid))
        p[t] += (hi-lo)
    return p

def meas(E):
    return dist(E)[0]

def primitive(E):
    nz = [e for e in E if e != 0]
    return reduce(gcd, nz) == 1 if nz else False

def LIM(B0):
    p = dist(B0)
    return p[0] + p[1]/7

if __name__ == "__main__":
    print("=== (A) baseline meas (exact) ===")
    for k in (8,9,10):
        C = list(range(k))
        m = meas(C)
        print(f"  k={k} consec meas_S7 = {m} = {float(m):.6f}")
    print()

    print("=== (B) Weyl limit exact-vs-large-w, AP and WIDE non-AP bases (k=8) ===")
    test_bases = [
        ("AP-7",      [0,1,2,3,4,5,6]),
        ("wide-nonAP",[0,1,3,7,12,20,33]),   # genuinely spread, primitive
        ("nonAP-2",   [0,2,3,5,6,9,13]),
    ]
    for name,B0 in test_bases:
        lim = LIM(B0)
        print(f"  {name} B0={B0}  LIM={float(lim):.6f} (={lim})")
        for w in [503, 1009, 2003, 5003]:
            E = B0 + [w]
            if not primitive(E): continue
            m = meas(E)
            print(f"     w={w}: meas={float(m):.6f}  diff_from_LIM={float(m-lim):+.6f}  |diff|*w={float(abs(m-lim))*w:.3f}")
    print()

    print("=== (C) Is consec the BOUNDED global max? search small windows ===")
    # exhaustive over primitive k-subsets of {0..W} containing 0, for modest W.
    for k, W in [(8, 12), (9, 12), (10, 13)]:
        Mk = meas(list(range(k)))
        best = Mk; bestset = list(range(k)); beaters = []
        rest = list(range(1, W+1))
        cnt = 0
        for combo in itertools.combinations(rest, k-1):
            S = (0,)+combo
            if reduce(gcd, S[1:]) != 1:  # need primitive
                continue
            cnt += 1
            m = meas(list(S))
            if m > best:
                best = m; bestset = list(S)
            if m > Mk:
                beaters.append((m, list(S)))
        print(f"  k={k} window {{0..{W}}}: tested {cnt} primitive sets")
        print(f"    consec M_{k} = {float(Mk):.6f}   best-in-window = {float(best):.6f} at {bestset}")
        if beaters:
            beaters.sort(reverse=True)
            print(f"    *** {len(beaters)} bounded sets BEAT consec! top: ", end="")
            print([(float(m), s) for m,s in beaters[:3]])
        else:
            print(f"    no bounded set in window beats consec.")
    print()

    print("=== (D) Hunt: any single-far set exceeds M_k? ===")
    for k, W in [(8, 9), (9, 10)]:
        Mk = meas(list(range(k)))
        rest = list(range(1, W+1))
        worst = Fraction(0); worstE = None
        viol = []
        for combo in itertools.combinations(rest, k-2):  # base B0 has k-1 elements incl 0
            B0 = [0]+list(combo)
            if reduce(gcd, B0[1:]) != 1 and reduce(gcd,B0[1:])!=0:
                pass  # base primitivity not required since far w fixes it
            for w in [50, 101, 211, 503, 1009]:
                E = B0 + [w]
                if not primitive(E): continue
                m = meas(E)
                if m > worst:
                    worst = m; worstE = list(E)
                if m > Mk:
                    viol.append((float(m), list(E)))
        print(f"  k={k}: M_{k}={float(Mk):.6f}  worst single-far meas={float(worst):.6f} at {worstE}")
        if viol:
            print(f"    *** {len(viol)} VIOLATIONS of M_k! {viol[:3]}")
        else:
            print(f"    no single-far set reached M_k. (reduction holds on this grid)")
