#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD4_span_robustness_macmini_S20.py  (THREAD 4, mac-mini-2026-06-21-S20)

ADVERSARIAL: the bounded-span exhaustive max in the cap audit only proves max measS7 <= cap
for E inside the box.  measS7 is scale-invariant, NOT span-bounded -- a primitive E can have
arbitrarily large span (e.g. {0,..,6,29}).  This script stresses whether the WORST measS7
escapes the box upward, especially at the BINDING row (k=8, margin only 0.054 -- the SMALLEST,
correcting the reframe's claim that k=10 binds).

For each k=8..11 we:
  (1) sweep span S = k-1 .. SMAX exhaustively (primitive E={0}+subset of {1..S}), tracking
      the running max measS7 and the argmax span;
  (2) report whether the max is ATTAINED at the consec span (=k-1) and STABLE as span grows
      (i.e. wider boxes never beat it) -- the empirical span-robustness of the maximizer.
A max that keeps RISING with span toward cap would threaten the cap; we check it plateaus.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

def measS7(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        ae = abs(e)
        for a in range(7 * ae + 1): bps.add(F(a, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo + hi) / 2
        hit = set()
        for e in E:
            v = e * xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += hi - lo
    return tot

def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1

CANON_CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91)}

print("#"*84)
print("# THREAD 4 -- SPAN-ROBUSTNESS of max measS7 (does the worst shape escape the box?)")
print("#"*84)

# SMAX chosen per k so the run finishes; k=8 widest (it's the binding row).
# (k=8 span<=18 = C(18,7) shapes per span; enough to show the plateau without OOM-time.)
SMAX = {8:18, 9:18, 10:16, 11:15}
for k in range(8, 12):
    cap = CANON_CAP[k]
    smax = SMAX[k]
    running_best = F(0); running_argE = None; running_argspan = None
    seen = set()
    print(f"\n=== k={k}, cap={cap}={float(cap):.6f}, span up to {smax} ===")
    for span in range(k - 1, smax + 1):
        # only NEW shapes with max(E)==span (avoid recount); enumerate E={0}+subset of {1..span} with span in it
        span_best = F(0); span_bestE = None
        for rest in itertools.combinations(range(1, span + 1), k - 1):
            if rest[-1] != span:  # ensure max==span (new shapes only)
                continue
            E = (0,) + rest
            if not primitive(E): continue
            m = measS7(list(E))
            if m > span_best: span_best, span_bestE = m, E
            if m > running_best: running_best, running_argE, running_argspan = m, E, span
        if span_bestE is not None:
            improved = "  <-- NEW GLOBAL MAX" if running_argspan == span and running_argE == span_bestE else ""
            print(f"   span={span:2d}: best-at-this-span measS7={float(span_best):.6f}  "
                  f"running max={float(running_best):.6f} (span {running_argspan}){improved}")
    margin = cap - running_best
    print(f"   --> running max measS7 = {running_best} = {float(running_best):.6f} at E={list(running_argE)} (span {running_argspan})")
    print(f"       cap - max = {float(margin):+.6f}  {'*** CAP VIOLATION ***' if margin<=0 else 'OK'}")
    print(f"       argmax span == consec span (k-1={k-1})? {running_argspan == k-1}")

print("\n" + "#"*84)
print("# If running max plateaus at the consec span and never approaches cap as span grows,")
print("# the bounded-span exhaustion is empirically the global max. (Not a proof of the")
print("# unbounded sup -- that is THM-557 single-block decorrelated extremality territory.)")
print("#"*84)
