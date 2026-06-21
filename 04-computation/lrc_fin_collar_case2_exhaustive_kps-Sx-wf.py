#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc_fin_collar_case2_exhaustive_kps-Sx-wf.py   (kps-2026-06-20, glue support)
=============================================================================
THM-547 case-2 EXHAUSTIVE finite check (NOT risk-ranked sampling).

Boundary collar: E = E' u {w}, E' subset {0,...,14} primitive-admissible with 0 in E',
|E'| = k-1, and 14 < w <= w*(k).  THM-547 case-1 (w > w*) is PROVED.  Case-2 was only
SPOT-checked upstream (k=8 risk-ranked sample, mac-mini).  Here we run case-2 EXHAUSTIVELY
and EXACTLY for k=8 (w* = 54): ALL C(14,7) choices of E' x ALL w in {15,...,54}.

If 0 violations: THM-547's case-2 is CLOSED EXACTLY for k=8 -> the entire boundary collar
(one far element) is PROVED for k=8 (case-1 + case-2 both done).  k=9,10 (w*=90/103) are the
same structure but heavier; we run a deterministic stride over E' for them and report.

EXACT rationals (fractions.Fraction).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}
WSTAR = {8: 54, 9: 90, 10: 103}
INNER = set(range(1, 7))

def sector_of(p):
    return int((p % 1) * 7)

def p0_inner(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    tot = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = {sector_of(e * xm) for e in E}
        if INNER <= secs:
            tot += x1 - x0
    return tot

def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, int(e))
    return g == 1

def run(k, exhaustive):
    base = list(range(1, 15))            # interior of {0..14}, 0 fixed in E'
    r = k - 1
    combos = list(itertools.combinations(base, r - 1))
    if not exhaustive:
        combos = combos[::max(1, len(combos) // 400)]
    worst = None
    viol = 0
    cnt = 0
    for combo in combos:
        Ep = (0,) + combo
        for w in range(15, WSTAR[k] + 1):
            E = Ep + (w,)
            if not primitive(E):
                continue
            cnt += 1
            p = p0_inner(E)
            m = CAPS[k] - p
            if m < 0:
                viol += 1
            if worst is None or m < worst[0]:
                worst = (m, E)
    tag = "EXHAUSTIVE" if exhaustive else "STRIDE-400"
    print(f"  k={k:2d} [{tag}] w in (14,{WSTAR[k]}]: configs={cnt:9d}  viol={viol}  "
          f"worst margin={float(worst[0]):+.6f} at E={worst[1]}")
    return viol

def main():
    print(__doc__)
    print("=" * 78)
    print("THM-547 CASE-2 EXHAUSTIVE FINITE CHECK")
    print("=" * 78)
    v8 = run(8, exhaustive=True)        # full: C(14,7)=3432 * 40 w ~ 137k configs, fast
    v9 = run(9, exhaustive=False)       # stride sample over E' (full is ~C(14,8)*76 heavy)
    v10 = run(10, exhaustive=False)
    print("\nVERDICT:")
    if v8 == 0:
        print("  k=8: THM-547 case-2 CLOSED EXACTLY (0 violations, EXHAUSTIVE).  Combined with")
        print("       case-1 (w>54, PROVED), the ONE-FAR boundary collar is fully PROVED for k=8.")
    else:
        print(f"  k=8: {v8} VIOLATIONS — recheck.")
    print(f"  k=9,10: stride samples, viol={v9},{v10} (full exhaustion is the remaining finite cost).")
    print("\nDONE.")

if __name__ == "__main__":
    main()
