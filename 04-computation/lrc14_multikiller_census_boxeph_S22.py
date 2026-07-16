#!/usr/bin/env python3
"""
lrc14_multikiller_census_boxeph_S22.py — the MULTI-KILLER stratum census (HYP-6295)

boxeph-2026-07-12-S22.  Pure Python, integer-exact.

CONTEXT.  Single-killer covering-min is CLOSED (opus-S253/S254 + kps cont.57 +
mac-mini-S68): deep well {1..12,182} unique min 14/183.  The open frontier is
MULTI-KILLER: primitive covering 13-families with NO multiple of 182 (13-duty and
14-duty on different speeds).  klein's ILP certifies all speeds <= 182, so the live
sub-stratum has Vmax > 182.

METHOD.  Exact M via sum-ruler enumeration:
    M(S) = max over pairs (i<j), m in [1, q-1], q = v_i+v_j, of (1/q) min_u ||u m||_q
justified by THM-722(D1) (the max of f is attained at a sum-handoff, which lies on a
pair-sum ruler).  All-integer inner loop with early abort; no gcd skip (MISTAKE-114).
All-odd families short-circuit to M = 1/2 (t = 1/2).

BATTERY.
  (a) baselines + unit tests against S21 ledger values;
  (b) dropped-element shapes {1..12}\\{x} u {13a, 14b} (covering-filtered);
  (c) the 14-far ladder {1..11, 13, 14m} (6|m) and 13-far mirror {1..11, 14, 156k};
  (d) DESIGNED adversaries: CRT-placed killers imitating the deep well at base 183;
  (e) random multi-killer sweep (no 182-multiple), low-M hunt.

Question: does the multi-killer stratum dip below 14/183?  (Crux predicts NO.)
"""

from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
import sys

FLOOR = F(14, 183)

# ----------------------------------------------------------------------------
# exact M, integer inner loop
# ----------------------------------------------------------------------------

def exact_M(speeds):
    """Exact M(S) = max_t min_i ||v_i t||, via sum-ruler enumeration (THM-722 D1).
    Returns (M, witness_q, witness_m)."""
    n = len(speeds)
    assert len(set(speeds)) == n and all(v > 0 for v in speeds)
    if all(v % 2 for v in speeds):
        return F(1, 2), 2, 1                     # t = 1/2, all odd
    best_n, best_d = 0, 1                        # best margin best_n/best_d
    best_wit = (0, 0)
    qs = sorted({speeds[i] + speeds[j]
                 for i in range(n) for j in range(i + 1, n)}, reverse=True)
    for q in qs:
        rs = sorted({u % q for u in speeds})
        half = q // 2
        for m in range(1, q):
            d = q
            ok = True
            for r in rs:
                x = (r * m) % q
                if x > half:
                    x = q - x
                if x < d:
                    d = x
                    if d * best_d <= best_n * q:   # cannot beat current best
                        ok = False
                        break
            if ok and d * best_d > best_n * q:
                best_n, best_d = d, q
                best_wit = (q, m)
    g = gcd(best_n, best_d)
    return F(best_n, best_d), best_wit[0], best_wit[1]

# ----------------------------------------------------------------------------
# predicates
# ----------------------------------------------------------------------------

def is_primitive(s):
    return reduce(gcd, s) == 1

def is_covering(s):
    return all(any(v % d == 0 for v in s) for d in range(2, 15))

def is_multikiller(s):
    """covering, primitive, and NO multiple of 182 (13- and 14-duty split)."""
    return is_covering(s) and is_primitive(s) and not any(v % 182 == 0 for v in s)

# ----------------------------------------------------------------------------
# battery
# ----------------------------------------------------------------------------

def show(name, fam, expect=None):
    fam = sorted(fam)
    M, q, m = exact_M(fam)
    mk = is_multikiller(fam)
    flag = "BELOW-FLOOR!!" if M < FLOOR else ("=FLOOR" if M == FLOOR else "")
    exp = ""
    if expect is not None:
        exp = f" [expect {expect} {'OK' if M == expect else '** MISMATCH **'}]"
    print(f"  {name:44s} M = {str(M):>10s} = {float(M):.6f}  wit q={q} m={m}"
          f"  MK={mk} {flag}{exp}")
    return M, q, m, mk


def main():
    print("MULTI-KILLER STRATUM CENSUS (HYP-6295) — floor test vs 14/183 = "
          f"{float(FLOOR):.6f}\n")

    # (a) unit tests vs S21 ledger values
    print("(a) unit tests")
    show("AP {1..13}", list(range(1, 14)), F(1, 14))
    show("GW {1..11,13,24}", list(range(1, 12)) + [13, 24], F(1, 14))
    show("deep well {1..12,182}", list(range(1, 13)) + [182], F(14, 183))
    show("ladder k=7 {1..12,91}", list(range(1, 13)) + [91], F(7, 92))
    show("compressed 2*{1..12}u{13}", [2 * i for i in range(1, 13)] + [13],
         F(1, 13))
    show("baseline {2..14}", list(range(2, 15)), F(1, 8))
    print()

    # (c) ladders first (cheap, structured)
    print("(c1) the 14-far ladder {1..11, 13, 14m}, 6|m (12-duty on the killer)")
    lad = []
    for m in range(6, 121, 6):
        if m % 13 == 0:
            continue                              # would be a 182-multiple
        fam = list(range(1, 12)) + [13, 14 * m]
        if not is_multikiller(fam):
            continue
        lad.append((m, show(f"{{1..11,13,{14*m}}} (m={m})", fam)))
    print()

    print("(c2) the 13-far mirror {1..11, 14, 156k} (12-duty on the 13-killer)")
    for k in range(1, 11):
        if k % 7 == 0:
            continue                              # 156k would be 182-multiple
        fam = list(range(1, 12)) + [14, 156 * k]
        if not is_multikiller(fam):
            continue
        show(f"{{1..11,14,{156*k}}} (k={k})", fam)
    print()

    # (d) designed adversaries: CRT-placed killers imitating the deep well @183
    print("(d) designed adversaries (killers CRT-placed near the 183-band edge)")
    show("design {1..11, 7332, 2744}", list(range(1, 12)) + [7332, 2744])
    # variants: w13 == 12 mod 183 & 156|w13 ; w14 == 182 mod 183 & 14|w14, 13∤(w14/14)
    show("design {1..11, 7332, 2562+182=2744} alt", list(range(1, 12)) + [7332, 2744 + 183 * 14])
    # keep 12 in the core, drop a middle element instead (x=6), killers plug 6-duty
    # w13 = 13a with 6|13a -> 6|a ; place w13 == -15 mod 183 : a ...
    print()

    # (b) dropped-element scan {1..12}\{x} u {13a, 14b}
    print("(b) dropped-element scan {1..12}\\{x} u {13a,14b}, a,b <= 24, "
          "covering+primitive+no-182")
    results = []
    seen = set()
    for x in range(1, 13):
        core = [i for i in range(1, 13) if i != x]
        for a in range(1, 25):
            if a % 14 == 0:
                continue
            for b in range(1, 25):
                if b % 13 == 0:
                    continue
                w13, w14 = 13 * a, 14 * b
                if w13 == w14 or w13 in core or w14 in core:
                    continue
                fam = tuple(sorted(core + [w13, w14]))
                if fam in seen:
                    continue
                seen.add(fam)
                if len(set(fam)) < 13 or not is_multikiller(list(fam)):
                    continue
                M, q, m = exact_M(list(fam))
                results.append((M, x, a, b, q, m, fam))
    results.sort()
    print(f"   {len(results)} families; 10 LOWEST M:")
    for (M, x, a, b, q, m, fam) in results[:10]:
        flag = "BELOW-FLOOR!!" if M < FLOOR else ("=FLOOR" if M == FLOOR else "")
        print(f"   M = {str(M):>8s} = {float(M):.6f}  x={x:2d} a={a:2d} b={b:2d}"
              f"  wit q={q} m={m}  {flag}")
    below = [r for r in results if r[0] < FLOOR]
    print(f"   families below 14/183: {len(below)}")
    print(f"   stratum min (this scan): {results[0][0]} = "
          f"{float(results[0][0]):.6f}  margin over floor: "
          f"{float(results[0][0] - FLOOR):+.6f}")
    print()

    # (e) random multi-killer sweep
    print("(e) random multi-killer sweep, Vmax <= 500, >=1 speed > 182")
    rng = random.Random(20260712)
    found, tried, lows = 0, 0, []
    while found < 60 and tried < 200000:
        tried += 1
        fam = sorted(rng.sample(range(1, 501), 13))
        if max(fam) <= 182:
            continue
        if not is_multikiller(fam):
            continue
        found += 1
        M, q, m = exact_M(fam)
        lows.append((M, tuple(fam), q, m))
    lows.sort()
    print(f"   {found} random multi-killer families (of {tried} draws)")
    for (M, fam, q, m) in lows[:5]:
        flag = "BELOW-FLOOR!!" if M < FLOOR else ""
        print(f"   M = {str(M):>10s} = {float(M):.6f}  wit q={q} m={m} {flag}")
        print(f"        fam={list(fam)}")
    below = [r for r in lows if r[0] < FLOOR]
    print(f"   below 14/183: {len(below)}")
    print()

    print("CENSUS COMPLETE.")

if __name__ == '__main__':
    main()
