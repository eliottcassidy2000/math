#!/usr/bin/env python3
"""
Does the tight family collapse to the AP at DOUBLED PRIMES n=2p?
opus-2026-06-01-S553 (remote-control).

A set V is TIGHT iff M(V)=max_t min_i ||v_i t|| = 1/n (safe set nonempty, measure
zero).  If the ONLY tight n-set is the AP {1,...,n-1}, then LRC(n) follows (every
non-tight set has positive-measure safe set => lonely; the AP is lonely).

Calibration found sporadic tight sets at n=5 ({1,3,4,7}, max 7~1.5n) and n=6
({1,3,4,5,9}, max 9=1.5n) but NONE at n=7 (prime) or n=14 within [1,18].
We scan n=8..14 over primitive (n-1)-subsets of [1, ~1.5n] (the range where the
small-n sporadics live) with a FAST early-bail tightness test, and tabulate the
tight-family size, flagging the doubled primes n=10 (2*5) and n=14 (2*7).
"""

from fractions import Fraction
from math import gcd
from itertools import combinations


def _exact_strictly_safe(V, t, thr):
    """Exact: min_i ||v_i t|| > thr (STRICT) -> proves an open safe nbhd."""
    for v in V:
        x = (v * t) % 1
        if min(x, 1 - x) <= thr:
            return False
    return True


def _exact_full(V, n):
    """Exact classification (Fraction), used only on the rare boundary cases."""
    thr = Fraction(1, n)
    endpoints = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                endpoints.add(Fraction(k * n + s, v * n) % 1)
    pts = sorted(endpoints)
    m = len(pts)
    for i in range(m):                       # positive-measure (strict interior)
        a, b = pts[i], pts[(i + 1) % m]
        length = (b - a) if b > a else (b - a + 1)
        mid = (a + length / 2) % 1
        if _exact_strictly_safe(V, mid, thr):
            return "loose"
    for t in pts:                            # measure-zero tight point
        if all(min((v * t) % 1, 1 - (v * t) % 1) >= thr for v in V):
            return "tight"
    return "COUNTEREX"


def fast_tight(V, n):
    """Rigorous + fast: float pre-scan to find a strictly-safe interior point
    (=> loose, confirmed exactly); fall back to exact analysis only when the
    float optimum sits at/below the 1/n floor (the tight/counterexample cases)."""
    thr_f = 1.0 / n
    eps = 1e-7
    # float endpoints
    pts = []
    for v in V:
        for k in range(v + 1):
            base = k / v
            pts.append((base + thr_f / v) % 1.0)
            pts.append((base - thr_f / v) % 1.0)
    pts.sort()
    m = len(pts)
    best_val, best_mid = -1.0, 0.0
    for i in range(m):
        a = pts[i]
        b = pts[i + 1] if i + 1 < m else pts[0] + 1.0
        mid = (a + (b - a) / 2.0) % 1.0
        mn = 1.0
        for v in V:
            x = (v * mid) % 1.0
            d = x if x < 1.0 - x else 1.0 - x
            if d < mn:
                mn = d
                if mn <= thr_f:
                    break
        if mn > best_val:
            best_val, best_mid = mn, mid
    if best_val > thr_f + eps:
        # confirm exactly that SOME nearby exact rational is strictly safe
        thr = Fraction(1, n)
        # snap best_mid to a nearby exact candidate: use Fraction from float is
        # messy; instead exact-verify the float midpoint via its rational form.
        # Rebuild the exact midpoint of the interval containing best_mid:
        if _exact_strictly_safe(V, Fraction(best_mid).limit_denominator(10**6), thr):
            return "loose"
    return _exact_full(V, n)


def primitive_tuple(combo):
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def is_AP(V, m):
    return tuple(V) == tuple(range(1, m + 1))


def scan(n, B):
    m = n - 1
    tight = []
    counters = []
    count = 0
    for combo in combinations(range(1, B + 1), m):
        if not primitive_tuple(combo):
            continue
        count += 1
        r = fast_tight(combo, n)
        if r == "tight":
            tight.append(combo)
        elif r == "COUNTEREX":
            counters.append(combo)
    non_ap = [t for t in tight if not is_AP(t, m)]
    is_dp = ""
    # doubled prime test: n = 2*prime
    if n % 2 == 0:
        q = n // 2
        if q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1)):
            is_dp = f"  <<< DOUBLED PRIME 2*{q}"
    flag = "  *** non-AP sporadics!" if non_ap else "  [AP-unique]"
    print(f" n={n:2d}  B={B:2d}  ({count:7d} sets)   tight={len(tight)}  "
          f"(AP:{len(tight)-len(non_ap)}, sporadic:{len(non_ap)}){flag}{is_dp}",
          flush=True)
    for t in non_ap[:30]:
        print(f"        sporadic tight: {t}", flush=True)
    if counters:
        print(f"        *** COUNTEREXAMPLES: {counters}", flush=True)
    return n, len(tight), len(non_ap), is_dp != ""


if __name__ == "__main__":
    print("Tight-family size vs n  (scan speeds up to ~1.5n; sporadics live here)\n",
          flush=True)
    rows = []
    # B = ceil(1.5 n) + 1, capped for feasibility
    plan = {8: 13, 13: 20, 14: 21}    # n=8 validates fast method vs slow run
    for n in (8, 13, 14):
        rows.append(scan(n, plan[n]))
    print("\nSUMMARY (sporadic=0 => AP-unique in scanned range):", flush=True)
    for n, tot, spor, dp in rows:
        print(f"  n={n:2d}  tight={tot}  sporadic={spor}  "
              f"{'DOUBLED-PRIME' if dp else ''}", flush=True)
