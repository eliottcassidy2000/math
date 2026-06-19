#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL re-verify of the AP-to-cap MARGIN on the uniform binding row |P|=5.
Independent reimplementation of meas_GP (lonely set) and meas_S (sector cover).

cap_k(n) = min_{|P|=n-1-k} meas(G_P), G_P={x: ||p x||>=1/n for all p in P}.
For |P|=5: k = (n-1)-5 = n-6.  margin = cap - meas_S(consec_k, s=n/2).

Claimed |P|=5 margins: n=8 +0.073, 10 +0.204, 12 +0.192, 14 +0.054(tightest),
16 +0.128, 18 +0.100, 20 +0.136.  Only n=14 < 0.10.

I implement meas_GP with my OWN breakpoint derivation:
||p x|| >= 1/n  <=>  frac(p x) in [1/n, 1-1/n].  The set in x is a union of arcs.
Breakpoints in x where frac(p x) = 1/n or 1-1/n: p x = m + 1/n or m + (1-1/n),
i.e. x = (m + 1/n)/p and x = (m + 1 - 1/n)/p for m=0..p-1.
"""
import sys, itertools
from fractions import Fraction
from math import comb

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')


def fp(v):
    return v - (v.numerator // v.denominator)


def meas_GP(P, n):
    """measure of {x in [0,1): for all p in P, frac(p x) in [1/n, 1-1/n]}."""
    P = sorted(set(p for p in P if p > 0))
    if not P:
        return Fraction(1)
    low = Fraction(1, n)
    high = Fraction(n - 1, n)
    bps = {Fraction(0), Fraction(1)}
    for p in P:
        for m in range(0, p + 1):
            for off in (low, high):
                x = Fraction(m, 1) + off
                x = x / p
                if Fraction(0) <= x <= Fraction(1):
                    bps.add(x)
    bps = sorted(bps)
    tot = Fraction(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        x = (lo + hi) / 2
        ok = True
        for p in P:
            v = fp(p * x)
            if v < low or v > high:
                ok = False
                break
        if ok:
            tot += (hi - lo)
    return tot


def meas_S(E, s):
    """measure of x where ALL sectors 0..s-1 hit by {frac(e x): e in E}."""
    E = sorted(set(E))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, s * e + 1):
            bps.add(Fraction(a, s * e))
    bps = sorted(b for b in bps if Fraction(0) <= b <= Fraction(1))
    p0 = Fraction(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        x = (lo + hi) / 2
        hit = set()
        for e in E:
            v = fp(e * x)
            hit.add((v.numerator * s) // v.denominator)
        if all(j in hit for j in range(s)):
            p0 += (hi - lo)
    return p0


def min_cap(n, psz):
    """min meas_GP over P subset of speeds {1..n-1}, |P|=psz. Exhaustive over real speed set."""
    speeds = list(range(1, n))
    best = None
    bestP = None
    nc = comb(len(speeds), psz)
    restricted = False
    if nc > 60000:
        # restrict to small windows but be generous
        speeds = list(range(1, min(psz + 8, n)))
        restricted = True
    for P in itertools.combinations(speeds, psz):
        m = meas_GP(P, n)
        if best is None or m < best:
            best = m
            bestP = P
    return best, bestP, restricted


def main():
    print("VALIDATION: reproduce cap_8(n=14)=2243/5880 over |P|=5 ?")
    cap8, P8, restr = min_cap(14, 5)
    print(f"  min meas_GP(n=14,|P|=5) = {cap8} = {float(cap8):.6f}  bestP={P8}  restricted={restr}")
    print(f"  matches 2243/5880={float(Fraction(2243,5880)):.6f}? {cap8==Fraction(2243,5880)}")
    ms8 = meas_S(list(range(8)), 7)
    print(f"  meas_S(consec_8,s=7)={ms8}={float(ms8):.6f}  (canon 481/1470)")
    print()

    print("=" * 80)
    print("UNIFORM BINDING ROW |P|=5 : margin = cap - meas_S(consec_{n-6}, s=n/2)")
    print("=" * 80)
    claimed = {8: 0.073, 10: 0.204, 12: 0.192, 14: 0.054, 16: 0.128, 18: 0.100, 20: 0.136}
    print(f"{'n':>3} {'k=n-6':>6} {'meas_S(AP)':>12} {'cap(|P|=5)':>14} {'margin':>9} {'claimed':>8} {'restr':>6}")
    results = {}
    for n in [8, 10, 12, 14, 16, 18, 20]:
        s = n // 2
        k = n - 6
        ms = meas_S(list(range(k)), s) if k >= 1 else Fraction(0)
        cap, Pb, restr = min_cap(n, 5)
        margin = cap - ms
        results[n] = float(margin)
        print(f"{n:>3} {k:>6} {float(ms):>12.5f} {float(cap):>14.5f} {float(margin):>+9.4f} {claimed.get(n,0):>+8.3f} {str(restr):>6}  P={Pb}")
    print()
    tight = min(results, key=results.get)
    print(f"  TIGHTEST margin at n={tight} ({results[tight]:+.4f})")
    below = [n for n, m in results.items() if m < 0.10]
    print(f"  margins < 0.10: {below}")


if __name__ == "__main__":
    main()
