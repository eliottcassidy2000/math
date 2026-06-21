#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_sharp_fast_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

Fast integer-only verification of the sharp 1D bound and its two faces.

row 0 = histogram of pq points, bin assignment via pure integers.
Point (k,t):  7pv with v in band; the bin = floor(7 frac(pk/q) + (2t+1)/(2q))  mod 7.
On common denominator 2q:  value*2q = 14 (p k mod q) + (2t+1)  [the "7 frac(pk/q)" term is
7*(pk mod q)/q, times 2q = 14 (pk mod q)].  So
    bin(k,t) = floor( (14*(p*k % q) + (2t+1)) / (2q) )  mod 7.
We accumulate r_j, j=0..6, then e_j = 7 r_j - p q,  sum_j|e_j| = 7 * S.
Target faces:  sum|e| <= 12 p   and   sum|e| <= 20 q.
"""
from math import gcd

P = 7

def evec_fast(p, q):
    r = [0] * P
    for k in range(q):
        base = 14 * ((p * k) % q)         # = 2q * 7 frac(pk/q)
        for t in range(p):
            val = base + (2 * t + 1)        # numerator over 2q
            j = (val // (2 * q)) % P
            r[j] += 1
    return [7 * rj - p * q for rj in r]

def main():
    print("THREAD C fast: sharp faces  sum|e| <= 12p  and  sum|e| <= 20q")
    print("=" * 72)
    QMAX = 300  # saved .out is the q<=300 run (31506 ratios, 0 violations); ~80s
    viol_p = []
    viol_q = []
    achieve_p = []
    achieve_q = []
    best_p = (0, 1, None)   # ratio num/den as (num, den) for sum|e|/p
    best_q = (0, 1, None)
    n = 0
    for q in range(1, QMAX + 1):
        pmax = (43 * q) // 20
        for p in range(q + 1, pmax + 1):
            if 20 * p > 43 * q:
                continue
            if gcd(p, q) != 1:
                continue
            e = evec_fast(p, q)
            S7 = sum(abs(x) for x in e)
            n += 1
            if S7 > 12 * p:
                viol_p.append((p, q, S7))
            if S7 > 20 * q:
                viol_q.append((p, q, S7))
            if S7 == 12 * p:
                achieve_p.append((p, q))
            if S7 == 20 * q:
                achieve_q.append((p, q))
            # compare S7/p to best_p (num/den)
            if S7 * best_p[1] > best_p[0] * p:
                best_p = (S7, p, (p, q, S7))
            if S7 * best_q[1] > best_q[0] * q:
                best_q = (S7, q, (p, q, S7))
    from fractions import Fraction as Fr
    print(f"checked {n} window ratios, q<={QMAX}")
    print(f"  violations sum|e|<=12p : {len(viol_p)}   first: {viol_p[:3]}")
    print(f"  violations sum|e|<=20q : {len(viol_q)}   first: {viol_q[:3]}")
    print(f"  sup sum|e|/p = {Fr(best_p[0],best_p[1])} ({best_p[0]/best_p[1]:.6f})  at {best_p[2]}  (target 12)")
    print(f"  sup sum|e|/q = {Fr(best_q[0],best_q[1])} ({best_q[0]/best_q[1]:.6f})  at {best_q[2]}  (target 20)")
    print(f"  #achieve face-P (=12p): {len(achieve_p)}   sample: {achieve_p[:15]}")
    print(f"  #achieve face-Q (=20q): {len(achieve_q)}   sample: {achieve_q[:15]}")

    # which ratios are simultaneously large?  the binding face by region
    print("\nFace assignment (which of 12p, 20q is the binding/active envelope):")
    print("  12p < 20q  <=>  p/q < 20/12 = 5/3.  So:")
    print("   - for 1 < p/q < 5/3:  12p is the SMALLER target (so 'sharper face' is 12p? "
          "no -- D needs an UPPER bound, so the bound is min only if BOTH hold).")
    # Actually we want sum|e| <= min over the two? No: each is a separately-true upper bound.
    # The COMBINED sharp envelope is sum|e| <= min(12p, 20q) IF both always hold; test it.
    viol_min = []
    for q in range(1, 200):
        pmax = (43 * q) // 20
        for p in range(q + 1, pmax + 1):
            if 20 * p > 43 * q or gcd(p, q) != 1:
                continue
            e = evec_fast(p, q)
            S7 = sum(abs(x) for x in e)
            if S7 > min(12 * p, 20 * q):
                viol_min.append((p, q, S7, 12 * p, 20 * q))
    print(f"  violations of sum|e| <= min(12p,20q) for q<400: {len(viol_min)}  {viol_min[:5]}")

if __name__ == "__main__":
    main()
