#!/usr/bin/env python3
"""
Regime-C arithmetic exploration (mac-mini-2026-07-03-S21).
Regime C = a covering 13-family with 7 near-equal SMALL far runners (w1..w1+6, w1 small) + 6 near.
The drifting floor fails; the problem is ARITHMETIC: at t=a/q the 7 far form an AP {base+j*a mod q}
on Z/q that must avoid the danger band [-q/14, q/14], SIMULTANEOUSLY with the near runners.

Questions:
 (1) Is every such family LONELY (exists t, all ||v_i t|| >= 1/14)?  [expect yes -- it's LRC14]
 (2) What is the arithmetic structure of the lonely time?  Rational a/q with small q?
 (3) The AP structure: at the lonely t, how do the 7 far AP terms sit in the safe band?
 (4) Is there a UNIFORM q that works (a single modulus covering all regime-C families)?
"""
from fractions import Fraction as F
from math import gcd

R = F(1, 14)  # danger radius

def min_dist_to_int(x):
    x = x % 1
    return min(x, 1 - x)

def M_numeric(speeds, N=200000):
    """max_t min_i ||v_i t|| via fine scan; returns (M, argmax t)."""
    best, bestt = F(0), F(0)
    for k in range(1, N):
        t = F(k, N)
        m = min(min_dist_to_int(v * t) for v in speeds)
        if m > best:
            best, bestt = m, t
    return best, bestt

def best_rational_lonely(speeds, qmax=400):
    """search t=a/q (q<=qmax) for the one maximizing min_i ||v_i t||; return (M, a, q)."""
    best, ba, bq = F(0), 0, 1
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            m = min(min_dist_to_int(F(v * a, q)) for v in speeds)
            if m > best:
                best, ba, bq = m, a, q
        if best >= R:  # found a lonely rational; keep searching for the SMALLEST q that works
            pass
    return best, ba, bq

def smallest_q_lonely(speeds, qmax=2000):
    """smallest q such that some a/q is lonely (min ||v_i a/q|| >= 1/14)."""
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            if all(min_dist_to_int(F(v * a, q)) >= R for v in speeds):
                # report the AP structure of the far block
                return q, a
    return None, None

if __name__ == "__main__":
    near = [1, 2, 3, 4, 5, 6]
    print("=" * 78)
    print("REGIME C: family = {1..6} + {w1..w1+6}, 13 speeds. Is it lonely? smallest lonely q?")
    print("=" * 78)
    print(f"{'w1':>6} {'M(scan)':>10} {'>=1/14?':>8} {'smallest q lonely':>18} {'a':>5} {'q/gcd hints'}")
    for w1 in [23, 24, 25, 29, 37, 50, 100, 200, 500, 1000, 2000, 5000, 7000]:
        far = list(range(w1, w1 + 7))
        speeds = near + far
        M, tstar = M_numeric(speeds, N=40000)
        q, a = smallest_q_lonely(speeds, qmax=1500)
        qstr = f"{q}" if q else ">1500"
        astr = f"{a}" if a else "-"
        # hint: is q related to 183=Phi6, or to the speeds?
        hint = ""
        if q:
            hint = f"q%7={q%7} q%14={q%14}"
        print(f"{w1:>6} {float(M):>10.5f} {str(M>=R):>8} {qstr:>18} {astr:>5} {hint}")

    print("\n" + "=" * 78)
    print("The AP structure at the lonely time (far block {w1..w1+6} at t=a/q):")
    print("=" * 78)
    for w1 in [23, 100, 1000, 5000]:
        far = list(range(w1, w1 + 7))
        speeds = near + far
        q, a = smallest_q_lonely(speeds, qmax=1500)
        if not q:
            print(f"w1={w1}: no lonely q<=1500")
            continue
        # the 7 far AP terms mod q
        ap = [(v * a) % q for v in far]
        band = q / F(14)  # danger half-width in Z/q units
        print(f"w1={w1}, t={a}/{q}: far AP mod q = {ap}")
        print(f"   danger band |x|<{float(band):.1f} (i.e. x in [0,{float(band):.1f}) or ({q}-{float(band):.1f},{q}))")
        print(f"   step a*1 mod q = {a% q} (consecutive far differ by a); near {[(v*a)%q for v in near]}")
