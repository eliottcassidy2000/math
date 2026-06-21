#!/usr/bin/env python3
"""ANGLE F (III) -- the LONGEST-DWELL invariant.

From Angle F (II): consec MAXIMIZES the longest single all-ones dwell with
ZERO violators at k=8,9,10 (full box max<=13), while it does NOT maximize the
return-count, mean-dwell, or atom-count.  This script:

  (1) Locates and characterizes consec's longest all-ones dwell interval
      [x0,x1] -- where is it, how wide, what is it in closed form.
  (2) Tests the longest-dwell extremality on a WIDER box (max<=16,17) and at
      k=11,12 to see how robust the "0 violators" is.
  (3) Tests a candidate BOUND: longest_dwell(E) >= measS7(E)/n_runs? trivially
      no.  Instead test the SHARPER structural claim:
         consec is the UNIQUE shape whose all-ones set contains an interval of
         length >= L*(k), and L*(k) has a clean closed form.
  (4) Reports the closed-form longest dwell for consec_k, k=8..14 and looks for
      a pattern (the CA's max "covered dwell" formula).

This is the honest yield of Angle F: a clean CA-extremal quantity (longest
covered dwell), recorded with its limits (it is NOT measS7 itself).

stdlib only; exact Fractions.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def breakpoints(E):
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for j in range(0, 7 * e + 1):
            bps.add(F(j, 7 * e))
    return sorted(bps)


def color_of(e, x):
    y = (F(e) * x) % 1
    return (y.numerator * 7) // y.denominator


def covered(E, x):
    return len({color_of(e, x) for e in E}) == 7


def longest_dwell_interval(E):
    """Return (length, lo, hi) of the widest all-ones interval (merged atoms)."""
    bps = breakpoints(E)
    best = F(0)
    best_lo = best_hi = None
    cur = F(0)
    run_lo = None
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        if covered(E, mid):
            if run_lo is None:
                run_lo = lo
            cur += hi - lo
            if cur > best:
                best = cur
                best_lo = run_lo
                best_hi = hi
        else:
            cur = F(0)
            run_lo = None
    return best, best_lo, best_hi


def measS7(E):
    bps = breakpoints(E)
    total = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        if covered(E, (lo + hi) / 2):
            total += hi - lo
    return total


def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return g == 1


def canon(E):
    E = tuple(sorted(E))
    mx = E[-1]
    refl = tuple(sorted(mx - e for e in E))
    def reduce(t):
        g = 0
        for e in t:
            g = gcd(g, e)
        return tuple(e // g for e in t) if g > 1 else t
    return min(reduce(E), reduce(refl))


def enumerate_box(k, span_max):
    seen = set()
    out = []
    for rest in combinations(range(1, span_max + 1), k - 1):
        E = (0,) + rest
        if not primitive(E):
            continue
        c = canon(E)
        if c in seen:
            continue
        seen.add(c)
        out.append(E)
    return out


def main():
    print("=" * 78)
    print("ANGLE F (III) -- the LONGEST covered-dwell invariant")
    print("=" * 78)

    # (1) characterize consec's longest dwell
    print("\n[1] consec_k longest all-ones dwell interval + closed form")
    print(f"{'k':>3} {'L*=longest':>16} {'lo':>10} {'hi':>10} "
          f"{'measS7':>12} {'L*/measS7':>10}")
    for k in range(8, 15):
        E = tuple(range(k))
        L, lo, hi = longest_dwell_interval(E)
        ms = measS7(E)
        ratio = float(L / ms) if ms else 0.0
        print(f"{k:>3} {str(L):>16} {str(lo):>10} {str(hi):>10} "
              f"{float(ms):>12.6f} {ratio:>10.4f}")

    # (2) robustness of longest-dwell extremality on wider boxes
    print("\n[2] longest-dwell extremality: # shapes beating consec's L*")
    print(f"{'k':>3} {'box':>5} {'#shapes':>9} {'L*viol':>8} {'measS7viol':>11}")
    for k in (8, 9, 10):
        for SPAN in (13, 16, 18):
            E0 = tuple(range(k))
            L0, _, _ = longest_dwell_interval(E0)
            m0 = measS7(E0)
            shapes = enumerate_box(k, SPAN)
            lv = 0
            mv = 0
            for E in shapes:
                if E == E0:
                    continue
                L, _, _ = longest_dwell_interval(E)
                if L > L0:
                    lv += 1
                if measS7(E) > m0:
                    mv += 1
            print(f"{k:>3} {SPAN:>5} {len(shapes):>9} {lv:>8} {mv:>11}")

    # (3) k=11,12 longest-dwell (where measS7 extremality is NOT needed -- B2
    #     already closes k>=11 -- but check the CA invariant still holds)
    print("\n[3] longest-dwell extremality at k=11,12 (box<=15)")
    for k in (11, 12):
        E0 = tuple(range(k))
        L0, _, _ = longest_dwell_interval(E0)
        shapes = enumerate_box(k, 15)
        lv = sum(1 for E in shapes
                 if E != E0 and longest_dwell_interval(E)[0] > L0)
        print(f"  k={k}: #shapes={len(shapes)}  L*={float(L0):.6f}  "
              f"L*-violators={lv}")


if __name__ == "__main__":
    main()
