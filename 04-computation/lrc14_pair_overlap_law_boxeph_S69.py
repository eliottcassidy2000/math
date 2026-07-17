#!/usr/bin/env python3
"""
THE PAIR-OVERLAP LAW (boxeph-2026-07-17-S69)
Extending the fleet's one open item (the 7-wall pair-crumb): the FULL-PERIOD
arithmetic core, exact.

Setting: danger sets D_v = {t in [0,1) : ||v t|| < 1/14} (the 1/14-danger
arcs: v arcs of width 1/(7v) at centers k/v).

(A) THE EXACT PAIR-OVERLAP FORMULA.  For speeds a <= b, g = gcd(a,b),
    delta = g/(ab) (the difference-grid spacing), W = (a+b)/(14ab) (the
    half-width sum), w_min = 1/(7b):
        mu(D_a cap D_b) = g * sum_{j in Z} min( (W - |j| delta)_+ , w_min ).
    Proof: arcs overlap iff their centers k/a, l/b are within W; center
    differences take exactly the values j*delta, each g times per period;
    the overlap at distance u is min((W-u)_+, w_min) (trapezoid).  5 lines.

(B) THE EQUIDISTRIBUTION VALUE: the capped-tent INTEGRAL is exactly
    (1/delta) int min((W-|u|)_+, w_min) du = ab * w_a w_b = 1/49 -- the
    independence value, for EVERY pair.  The overlap law is the discrete
    sampling of a trapezoid whose mean is exactly 1/49.

(C) THE NEAR-EQUAL FLOOR (the 7-wall face): for reduced consecutive ratios
    b/g = a/g + 1 with a' = a/g:
        mu = 1/(7(a'+1)) for a' <= 6;  mu >= 1/49 for ALL a' (verified to
    a' = 200; equality at a' = 6 and 7).  More generally the floor
    mu >= 1/49 holds on the near-equal domain (mapped exactly below) and
    FAILS for lopsided reduced ratios (mu(D_1 cap D_12) = 1/84 < 1/49):
    the pair-crumb needs the block's comparability, and the exact formula
    says precisely how much.

(D) THE WALL-CROSSING REFEREE: for c = 7 near-equal blocks, the hunter
    ledger good >= 1 - c/7 + sum_{consecutive} mu(D_i cap D_{i+1}) is
    positive using the EXACT pair values -- the union bound dies (1 - 7/7
    = 0) but the pair credits cross the wall; checked against the true
    good measure (exact sweep) on explicit blocks.

(E) THE WINDOW FACE (measured, for mac-mini): min over window positions of
    mu(D_a cap D_b cap I)/|I| at various |I| -- how fast the full-period
    floor localizes.

All exact (Fractions).  Independent of, and feeding, HYP-3874/HYP-4021.
"""

import sys
from fractions import Fraction as Fr
from math import gcd

FT = Fr(1, 14)


def mu_exact(a, b):
    """the (A)-formula, exact."""
    g = gcd(a, b)
    W = Fr(a + b, 14 * a * b)
    wmin = Fr(1, 7 * max(a, b))
    delta = Fr(g, a * b)
    tot = Fr(0)
    j = 0
    while True:
        u = j * delta
        if u >= W:
            break
        ol = min(W - u, wmin)
        tot += ol if j == 0 else 2 * ol
        j += 1
    return g * tot


def danger_breaks(v):
    """boundaries of D_v in [0,1): (m +- 1/14)/v."""
    out = []
    for m in range(v + 1):
        for s in (-1, 1):
            x = (Fr(m) + s * FT) / v
            if 0 <= x < 1:
                out.append(x)
    return out


def in_danger(v, t):
    f = (v * t) % 1
    return min(f, 1 - f) < FT


def mu_brute(speeds, want_all=True):
    """exact measure of {t : all/any speeds in danger} by breakpoint sweep."""
    bps = sorted(set([Fr(0), Fr(1)] + [x for v in speeds
                                       for x in danger_breaks(v)]))
    tot = Fr(0)
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        flags = [in_danger(v, mid) for v in speeds]
        if (all(flags) if want_all else any(flags)):
            tot += bps[i + 1] - bps[i]
    return tot


if __name__ == "__main__":
    print("THE PAIR-OVERLAP LAW -- referee (boxeph S69)")
    print("=" * 78)
    print("PART A -- exact formula vs brute sweep")
    nA = 0
    for a in range(1, 25):
        for b in range(a, 37, 3):
            m1 = mu_exact(a, b)
            m2 = mu_brute([a, b])
            assert m1 == m2, (a, b, m1, m2)
            nA += 1
    print(f"  formula == brute on {nA} pairs (exact rationals)")

    print()
    print("PART B/C -- the 1/49 equidistribution value and the domain map")
    # consecutive reduced ratios
    row = []
    for ap in range(1, 9):
        m = mu_exact(ap, ap + 1)
        row.append(f"({ap},{ap+1}): {m} ({'=1/%d' % (7*(ap+1)) if m == Fr(1, 7*(ap+1)) else '>1/49' if m > Fr(1,49) else '=1/49' if m == Fr(1,49) else 'BELOW'})")
    print("  consecutive: " + "; ".join(row[:6]))
    print("               " + "; ".join(row[6:]))
    worst_consec = min(mu_exact(ap, ap + 1) for ap in range(1, 201))
    print(f"  min over consecutive a' <= 200: {worst_consec} "
          f">= 1/49: {worst_consec >= Fr(1, 49)}")
    assert worst_consec >= Fr(1, 49)
    for ap in range(1, 7):
        assert mu_exact(ap, ap + 1) == Fr(1, 7 * (ap + 1))
    # the failure map: reduced pairs with mu < 1/49
    fails = []
    for ap in range(1, 15):
        for bp in range(ap + 1, 15):
            if gcd(ap, bp) == 1 and mu_exact(ap, bp) < Fr(1, 49):
                fails.append((ap, bp, float(bp / ap)))
    minratio = min(f[2] for f in fails)
    print(f"  reduced pairs a',b' <= 14 with mu < 1/49: {len(fails)}; "
          f"smallest failing ratio b'/a' = {minratio:.3f} "
          f"(e.g. {[f[:2] for f in fails if f[2] == minratio]})")
    print(f"  spot: mu(D_1, D_12) = {mu_exact(1, 12)} (= 1/84 < 1/49: "
          f"{mu_exact(1, 12) == Fr(1, 84)})")
    # the safe-ratio threshold
    ok_ratio = None
    for num, den in [(3, 2), (2, 1), (5, 2), (3, 1), (7, 2), (4, 1)]:
        allok = all(mu_exact(ap, bp) >= Fr(1, 49)
                    for ap in range(1, 30) for bp in range(ap + 1, 30)
                    if gcd(ap, bp) == 1 and Fr(bp, ap) <= Fr(num, den))
        print(f"  all reduced pairs with ratio <= {num}/{den} pass the 1/49 "
              f"floor: {allok}")
        if allok:
            ok_ratio = Fr(num, den)

    print()
    print("PART D -- crossing the 7-wall with exact pair credits")
    for v0, d in [(50, 1), (100, 1), (35, 7), (70, 14)]:
        block = [v0 + i * d for i in range(7)]
        union = mu_brute(block, want_all=False)
        good_true = 1 - union
        credits = sum(mu_exact(block[i], block[i + 1]) for i in range(6))
        ledger = 1 - Fr(7, 7) + credits
        print(f"  block {block}: TRUE good = {float(good_true):.5f}; hunter "
              f"ledger floor = 0 + credits {float(credits):.5f} = "
              f"{float(ledger):.5f} > 0: {ledger > 0} (<= true: "
              f"{ledger <= good_true})")
        assert ledger > 0 and ledger <= good_true

    print()
    print("PART E -- the window face (measured): min_I mu(I)/|I| vs |I|")
    a, b = 50, 51
    full = mu_exact(a, b)
    for Lden in (2, 5, 10, 25):
        L = Fr(1, Lden)
        worst = None
        for k in range(0, Lden * 2):
            t0 = Fr(k, 2 * Lden)
            bps = sorted(set([t0, t0 + L] +
                             [x for v in (a, b) for x in danger_breaks(v)
                              if t0 <= x <= t0 + L]))
            tot = Fr(0)
            for i in range(len(bps) - 1):
                mid = (bps[i] + bps[i + 1]) / 2
                if in_danger(a, mid) and in_danger(b, mid):
                    tot += bps[i + 1] - bps[i]
            r = tot / L
            worst = r if worst is None else min(worst, r)
        print(f"  (a,b)=(50,51): |I| = 1/{Lden}: min_I mu/|I| = "
              f"{float(worst):.5f}  (period value {float(full):.5f})")
    print("=" * 78)
    print("done")
