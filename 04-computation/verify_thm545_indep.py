#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEPENDENT adversarial verification of THM-545 (even-graph parity-completion).
Own exact-rational lonely-measure engine; cross-check vs known values; hunt counterexamples.
"""
import itertools
from fractions import Fraction
from functools import reduce
from math import gcd

THETA = Fraction(1, 14)

def lonely_measure(C, theta=THETA):
    """
    meas{ t in [0,1) : ||d*t|| > theta  for all d in C }.
    ||x|| = distance to nearest integer.
    Danger set for speed d: { t : ||d t|| <= theta } = union over m of [m/d - theta/d, m/d + theta/d].
    We compute the danger union on [0,1) by collecting interval endpoints as exact fractions
    and merging. We INCLUDE endpoints conservatively but since we use strict '>' the boundary
    has measure zero, so closed vs open does not change the measure.
    Returns 1 - |danger union|.
    """
    segs = []
    for d in C:
        if d == 0:
            continue
        w = theta / d
        # centers m/d for m = 0 .. d ; arcs [m/d - w, m/d + w]; wrap into [0,1)
        for m in range(0, d + 1):
            lo = Fraction(m, d) - w
            hi = Fraction(m, d) + w
            for shift in (-1, 0, 1):
                a = lo + shift
                b = hi + shift
                a2 = a if a > 0 else Fraction(0)
                b2 = b if b < 1 else Fraction(1)
                if a2 < b2:
                    segs.append((a2, b2))
    segs.sort()
    covered = Fraction(0)
    cur_lo = None
    cur_hi = None
    for a, b in segs:
        if cur_hi is None:
            cur_lo, cur_hi = a, b
        elif a <= cur_hi:
            if b > cur_hi:
                cur_hi = b
        else:
            covered += cur_hi - cur_lo
            cur_lo, cur_hi = a, b
    if cur_hi is not None:
        covered += cur_hi - cur_lo
    return Fraction(1) - covered

def carry_profile(C):
    d = {}
    for e in C:
        if e == 0:
            continue
        m = e; a = 0
        while m % 2 == 0:
            m //= 2; a += 1
        d[m] = d.get(m, 0) + 2**a
    return dict(sorted(d.items()))

def ap_tail_core(holes, tails):
    return tuple(sorted([x for x in range(1, 14) if x not in holes] + list(tails)))

THR1 = Fraction(7, 858)
THR2 = Fraction(426, 35035)

def main():
    print("=== Cross-check known exact values (independent engine) ===")
    checks = {
        "drop-6":      (ap_tail_core({6}, []),        Fraction(7, 858)),
        "{5:+2} exc":  (ap_tail_core({6, 10}, [20]),  Fraction(3859, 420420)),
        "drop-12":     (ap_tail_core({12}, []),       Fraction(426, 35035)),
        "two-tail min":(ap_tail_core({4,6,10},[20,46]),Fraction(50189, 3223220)),
    }
    all_ok = True
    for name, (C, expected) in checks.items():
        L = lonely_measure(C)
        ok = (L == expected)
        all_ok &= ok
        print(f"  {name:13s}: meas={L}  expected={expected}  MATCH={ok}  carry={carry_profile(C)}")
    print(f"  ALL KNOWN VALUES MATCH: {all_ok}")
    print(f"  thr1=7/858={float(THR1):.8f}  thr2=426/35035={float(THR2):.8f}")

    # 20->40 composition claim
    C2040 = tuple(sorted((set(ap_tail_core({6,10},[20])) - {20}) | {40}))
    print(f"\n  20->40 composition: {C2040}")
    print(f"    len={len(C2040)} gcd={reduce(gcd,C2040)} meas={lonely_measure(C2040)}={float(lonely_measure(C2040)):.8f} (expect 2669/194040={float(Fraction(2669,194040)):.8f})")

    # single bit-shift table
    print("\n=== Single bit-shifts e->2e from drop-6 ===")
    drop6 = set(ap_tail_core({6}, []))
    res = []
    for m in [1,3,5,7,9,11,13]:
        for a in range(0,6):
            e = (2**a)*m; e2 = (2**(a+1))*m
            if e in drop6 and e2 not in drop6:
                C = tuple(sorted((drop6 - {e}) | {e2}))
                if len(C) == 12 and reduce(gcd, C) == 1:
                    res.append((lonely_measure(C), e, e2))
    res.sort()
    for L, e, e2 in res:
        tag = "GAUGE (<thr2)" if L < THR2 else "SPENDING (>=thr2)"
        print(f"    {e:2d}->{e2:2d}: meas={L}={float(L):.8f}  {tag}")

    return all_ok

if __name__ == "__main__":
    main()
