#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) — the Erdos-Turan / single-frequency RATE for the far element, and the
resulting explicit B(k).
kind-pasteur-2026-06-19-S12.

We measured: max_B0 LIM(B0) is strictly below M_k (the bounded global max=consec):
   gap_k = M_k - max_B0 LIM(B0) = 0.13061 (k=8), 0.05406 (k=9).

We need:  |meas(S7)(B0 u {w}) - LIM(B0)| <= R(w),  with R(w) -> 0.
Then for w with R(w) < gap_k we have meas(S7)(B0 u {w}) < M_k for ALL bases B0,
giving an explicit B(k) = min w past which R(w) < gap_k.

The far element w hits sectors of [0,1) via frac(w x). Over the breakpoint
structure of the bounded orbit, the conditional law of frac(w x) on each
bounded-orbit cell is an arithmetic progression with step 1/w sampled across a
sub-interval; its sector-occupancy deviates from uniform by at most ~ (#walls)/w
= 6/w per cell, i.e. the single-frequency discrepancy is O(1/w) with an explicit
constant. We MEASURE max_B0 |meas(B0u{w}) - LIM(B0)| * w to bound the constant.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def dist_missed_count(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
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
        hit = set()
        for e in E:
            v = e*mid; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        t = sum(1 for j in range(1,7) if j not in hit)
        p[t] += (hi-lo)
    return p

def meas_S7(E):
    return dist_missed_count(E)[0]

def primitive(E):
    return reduce(gcd, [e for e in E if e != 0]) == 1

if __name__ == "__main__":
    print("=== Erdos-Turan rate: max_B0 |meas(B0u{w}) - LIM(B0)| * w ===\n")
    gap = {8: 0.13061, 9: 0.05406}
    # bases: consec (k-1) is the LIM-maximizer; test it plus a few others.
    for k in [8, 9]:
        km1 = k-1
        bases = [list(range(km1))]
        # add a couple of bounded non-consec bases
        if k == 8:
            bases += [[0,1,2,3,4,5,7],[0,1,2,3,4,6,7]]
        else:
            bases += [[0,1,2,3,4,5,6,8],[0,1,2,3,4,5,7,8]]
        print(f"--- k={k}  gap_k={gap[k]} ---")
        for B0 in bases:
            p = dist_missed_count(B0)
            lim = float(p[0] + p[1]*Fraction(1,7))
            maxscaled = 0.0; worstw = None; maxabs = 0.0
            for w in range(max(B0)+2, 200):
                E = B0 + [w]
                if not primitive(E): continue
                m = float(meas_S7(E))
                err = abs(m - lim)
                if err*w > maxscaled:
                    maxscaled = err*w; worstw = w
                # track when err first stays below gap
            # find B(k): smallest w0 s.t. for all w>=w0, err<gap
            firstgood = None
            errs = []
            for w in range(max(B0)+2, 200):
                E = B0 + [w]
                if not primitive(E): continue
                m = float(meas_S7(E)); err = abs(m-lim)
                errs.append((w, err))
            # last w where err >= gap
            bad = [w for (w,e) in errs if lim + e >= float({8:0.32721,9:0.41616}[k]) ]  # meas could reach M_k
            print(f"  B0={B0}: LIM={lim:.5f}, max|err|*w={maxscaled:.4f} (at w={worstw})")
            # report w past which meas < M_k guaranteed (empirically meas itself < M_k)
            Mk = {8:0.32721, 9:0.41616}[k]
            viol = [(w, round(lim+e,5)) for (w,e) in errs if (lim+e) > Mk]
            print(f"     LIM+max-tail-err exceeds M_k={Mk}?  worst meas-bound={lim+max(e for _,e in errs):.5f}  {'OK below M_k' if lim+max(e for _,e in errs)<=Mk else 'CHECK'}")
        print()
