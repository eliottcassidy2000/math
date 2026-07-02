#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S100 -- HYP-3856: the EXACT PER-PATTERN MIN-OVERLAP TABLE.
For coprime (P,Q), ov_{P,Q}(theta) = |{x : ||Px|| < r, ||Qx - theta|| < r}| is
piecewise linear in theta with rational breakpoints; we compute its EXACT min
(Fraction arithmetic) for all coprime P <= Q with PQ <= 64, at r = 1/14.
Purpose: replace THM-598(A)'s Fourier deficit bound (min >= (2r)^2 - 1/(3PQ))
with an exact, decide-checkable table -- the Lean-ready normal form.  Also
reports the TRUE dangerous list (exact min = 0 or below the THM-599 slack).
"""
from fractions import Fraction as F
from math import gcd

r = F(1, 14)

def overlap_at(P, Q, theta):
    """|{x in [0,1): ||Px||<r, ||Qx-theta||<r}| exactly."""
    iv = []
    for k in range(P):
        a, b = (F(k) - r) / P, (F(k) + r) / P
        iv.append((a, b))
    tot = F(0)
    for m in range(Q):
        c, d = (F(m) + theta - r) / Q, (F(m) + theta + r) / Q
        for (a, b) in iv:
            for sh in (-1, 0, 1):
                lo, hi = max(a, c + sh), min(b, d + sh)
                if hi > lo:
                    tot += hi - lo
    return tot

def exact_min(P, Q):
    # breakpoints of theta: where a Q-arc endpoint meets a P-arc endpoint:
    # (m + theta ± r)/Q = (k ± r)/P + integer  =>  theta in finite rational set
    bps = set()
    for k in range(P):
        for m in range(Q):
            for s1 in (-1, 1):
                for s2 in (-1, 1):
                    for j in (-1, 0, 1):
                        th = (F(k) + s1 * r) * Q / P - F(m) - s2 * r + j
                        th -= int(th)
                        if th < 0: th += 1
                        bps.add(th)
    bps = sorted(bps)
    best = None
    for i in range(len(bps)):
        a = bps[i]
        b = bps[(i + 1) % len(bps)] if i + 1 < len(bps) else bps[0] + 1
        mid = (a + b) / 2
        v = overlap_at(P, Q, mid)
        # also endpoints (piecewise linear: min at breakpoints or flat mids)
        ve = overlap_at(P, Q, a)
        for val in (v, ve):
            if best is None or val < best:
                best = val
    return best

indep = (2 * r) ** 2
print(f"r = 1/14, independence (2r)^2 = {indep} = {float(indep):.6f}")
print(f"{'(P,Q)':>8} {'PQ':>4} {'exact min':>12} {'float':>9} {'Fourier bd':>11} {'verdict':>22}")
dangerous = []
for Q in range(1, 33):
    for P in range(1, Q + 1):
        if gcd(P, Q) != 1 or P * Q > 64: continue
        mn = exact_min(P, Q)
        fb = indep - F(1, 3 * P * Q)
        verdict = "DANGEROUS (min=0)" if mn == 0 else ("safe (min>0)" if mn > 0 else "?")
        if mn == 0: dangerous.append((P, Q))
        if P * Q <= 20 or mn == 0:
            print(f"  ({P},{Q}) {P*Q:>4} {str(mn):>12} {float(mn):>9.5f} {float(fb):>11.5f} {verdict:>22}")
print(f"\nTRUE dangerous list (exact min = 0): {dangerous}")
print(f"Fourier-suspect list was PQ <= 16; exact table verdict above.")
