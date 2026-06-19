#!/usr/bin/env python3
"""
Probe the B_7 <= mu_1/7 lower-bound claim, which the main verification found
VIOLATED for consec_9 (B7=0.88563 > mu=0.84014).

The prompt's logic: "empty length-1/7 arc => maxgap >= 1/7 => contributes to mu
a.e." But mu_1/7(E) = meas{maxgap > 1/7} (STRICT). An empty FIXED arc of width
exactly 1/7 only gives maxgap >= 1/7, NOT > 1/7. The set where a fixed width-1/7
arc is empty but maxgap == 1/7 exactly can have POSITIVE measure (e.g. when two
points sit exactly on the arc endpoints over an interval of x). That is the bug.

We confirm by:
 (1) computing B7 and mu for consec_9 EXACTLY and locating the x-interval where
     a fixed arc is empty but maxgap == 1/7 (not > 1/7);
 (2) computing the corrected minorant B7_strict = meas{some fixed arc is empty
     AND maxgap > 1/7}, and checking B7_strict <= mu (should hold).
"""
from fractions import Fraction as F

def maxgap_at(E, x):
    pts = sorted(set((e*x) % 1 for e in E))
    if len(pts) == 1:
        return F(1)
    gaps = []
    for i in range(len(pts)):
        nxt = pts[(i+1) % len(pts)]
        gaps.append(nxt - pts[i] if i+1 < len(pts) else (1-pts[i])+pts[0])
    return max(gaps)

def some_fixed_arc_empty(E, x):
    arcs = [(F(2*i+1,14), F(2*i+3,14)) for i in range(7)]  # arc6 wraps past 1
    pts = [(e*x) % 1 for e in E]
    for (lo, hi) in arcs:
        # membership lo <= p < hi, with hi possibly >1 (arc 6 = [13/14,15/14))
        has = False
        for p in pts:
            if hi <= 1:
                if lo <= p < hi: has = True; break
            else:
                if p >= lo or p < (hi-1): has = True; break
        if not has:
            return True
    return False

def common_breakpoints(E):
    E = sorted(set(E)); n = len(E)
    bp = set([F(0), F(1)])
    # order-cell breakpoints
    for i in range(n):
        for j in range(i+1, n):
            d = E[j]-E[i]
            for m in range(0, d+1): bp.add(F(m, d))
    # fixed-arc-endpoint breakpoints: e*x = m + a, a in {(2i+1)/14}
    for e in E:
        if e == 0: continue
        for m in range(0, e+1):
            for i in range(7):
                bp.add((F(m)+F(2*i+1,14))/e)
                bp.add((F(m)+F(2*i+3,14))/e)
    return sorted(b for b in bp if 0 <= b <= 1)

def measure_pred(E, pred):
    bp = common_breakpoints(E)
    total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2
        if pred(E, mid):
            total += (b-a)
    return total

def main():
    for E in [list(range(9)), list(range(8)), list(range(10))]:
        b7      = measure_pred(E, some_fixed_arc_empty)
        b7strict= measure_pred(E, lambda E,x: some_fixed_arc_empty(E,x) and maxgap_at(E,x) > F(1,7))
        mu      = measure_pred(E, lambda E,x: maxgap_at(E,x) > F(1,7))
        muge    = measure_pred(E, lambda E,x: maxgap_at(E,x) >= F(1,7))
        # region where a fixed arc is empty but maxgap == 1/7 exactly (the leak):
        leak    = measure_pred(E, lambda E,x: some_fixed_arc_empty(E,x) and maxgap_at(E,x) == F(1,7))
        print(f"E=consec_{len(E)}:")
        print(f"   B7 (some fixed arc empty)            = {b7} = {float(b7):.6f}")
        print(f"   B7_strict (empty AND maxgap>1/7)     = {b7strict} = {float(b7strict):.6f}")
        print(f"   mu_1/7 (maxgap > 1/7)                = {mu} = {float(mu):.6f}")
        print(f"   meas(maxgap >= 1/7)                  = {muge} = {float(muge):.6f}")
        print(f"   LEAK (arc empty & maxgap == 1/7)     = {leak} = {float(leak):.6f}")
        print(f"   B7 <= mu ?  {b7 <= mu}    B7_strict <= mu ?  {b7strict <= mu}")
        print(f"   B7 <= meas(maxgap>=1/7) ? {b7 <= muge}")
        print()

if __name__ == "__main__":
    main()
