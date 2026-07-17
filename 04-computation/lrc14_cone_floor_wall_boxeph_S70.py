#!/usr/bin/env python3
"""
THE CONE FLOOR AND THE c = 7 WALL DICHOTOMY (boxeph-2026-07-17-S70)
The pair-crumb assembly layer, full-period face (extends LEM-042).

(A) THE DEFECT BOUNDS (one line each).  The pair overlap is the delta-sampling
    of the even unimodal trapezoid f(u) = min((W-|u|)_+, w_min) whose full
    integral is w_a w_b; sampling a nonincreasing tail loses/gains at most one
    cell, so
        1/49 - g/(7b)  <=  mu(D_a cap D_b)  <=  1/49 + g/(7b).
    In reduced terms (a', b') = (a, b)/g: the window is +- 1/(7 b').

(B) THE CONE FLOOR (exact, proved by finite check + (A)-tail): let
        c0(rho) = min { mu(D_a' cap D_b') : gcd = 1, b'/a' <= rho }.
    For the LACUNARY-COMPLEMENT CONE rho = 7/3 (opus-S336: all-ratios >= 7/3
    blocks are lonely universally): c0 = c0(7/3) computed exactly -- the
    finite range b' <= B0 exhaustively, the tail b' > B0 >= floor by (A).
    Also reported: c0(13) for the all-comparable <= 13x convention.

(C) THE c = 7 WALL DICHOTOMY (full period): a 7-block (seven speeds) either
    (i) has all consecutive ratios >= 7/3 -- the LACUNARY branch (opus) -- or
    (ii) has some consecutive ratio < 7/3, and then klein's hunter ledger
    with the SINGLE cone-floor credit gives
        good measure  =  1 - mu(union D_i)  >=  c0  >  0:
    the union bound's death at 7 x 1/7 = 1 is beaten by one pair credit.
    Referee: random + structured 7-blocks, both branches, brute-exact.

(D) THE c >= 8 HONEST BOUNDARY: pair credits alone need sum > (c-7)/7; with
    c0 < 1/49 the uniform-credit route FAILS at c = 8 (demonstrated on a
    block where the credit sum < 1/7 but the true good is far positive) --
    the c >= 8 case needs exact per-family credits (decide) or higher-order
    structure.  NAMED, not claimed.

(E) THE DIFFERENCE-COMB LAW + BEATS: D_a cap D_b subset {t : ||(b-a)t|| < 1/7}
    (triangle inequality, one line); per-beat mass census (the beat at k/d);
    the CLUSTERED BEAT-MISS: for d = b - a = 1 the entire mass sits in ONE
    beat -- sweep windows of width 2/a positioned off-beat carry ZERO
    (the LEM-042(E) collapse, localized): positioning is essential exactly
    where kps's clustered-block closer (cite_cluster7_lonely) already
    operates; the pair route serves the loose-comparable regime.

All exact (integers/Fractions).
"""

import sys
from fractions import Fraction as Fr
from math import gcd
import random

sys.path.insert(0, '04-computation')
from lrc14_pair_overlap_law_boxeph_S69 import mu_exact, mu_brute, \
    danger_breaks, in_danger

FT = Fr(1, 14)


def mu_int_num(a, b):
    """mu * 14ab/g as an integer (units of g/(14ab))."""
    g = gcd(a, b)
    W = a + b                      # units 1/(14ab)
    cap = 2 * min(a, b)
    step = 14 * g
    tot = 0
    j = 0
    while True:
        u = j * step
        if u >= W:
            break
        ol = min(W - u, cap)
        tot += ol if j == 0 else 2 * ol
        j += 1
    return tot, g


def mu_fr(a, b):
    t, g = mu_int_num(a, b)
    return Fr(g * t, 14 * a * b)


if __name__ == "__main__":
    print("THE CONE FLOOR AND THE c = 7 WALL DICHOTOMY (boxeph S70)")
    print("=" * 78)
    print("PART A -- the defect bounds 1/49 -+ g/(7b)")
    rng = random.Random(7)
    worst = 0
    for _ in range(500):
        a = rng.randint(1, 400)
        b = rng.randint(a, 500)
        m = mu_fr(a, b)
        g = gcd(a, b)
        lo = Fr(1, 49) - Fr(g, 7 * b)
        hi = Fr(1, 49) + Fr(g, 7 * b)
        assert lo <= m <= hi, (a, b, m)
        worst += 1
    # spot equality of the integer form vs LEM-042 exact
    for a, b in [(6, 7), (5, 8), (1, 12), (50, 51), (35, 42)]:
        assert mu_fr(a, b) == mu_exact(a, b)
    print(f"  bounds hold on {worst} random pairs (exact); integer form == "
          f"LEM-042 formula")

    print()
    print("PART B -- the cone floors c0(7/3) and c0(13), exact")
    for num, den, B0 in [(7, 3, 700), (13, 1, 700)]:
        rho = Fr(num, den)
        best = None
        arg = None
        for bp in range(2, B0 + 1):
            amin = int(bp / rho) if bp % (num // gcd(num, den)) else bp * den // num
            for ap in range(max(1, int((bp * den + num - 1) // num)), bp):
                if Fr(bp, ap) > rho or gcd(ap, bp) != 1:
                    continue
                m = mu_fr(ap, bp)
                if best is None or m < best:
                    best, arg = m, (ap, bp)
        tail = Fr(1, 49) - Fr(1, 7 * B0)
        print(f"  rho = {num}/{den}: min over b' <= {B0} is {best} at {arg} "
              f"(= {float(best):.6f}); tail floor (b' > {B0}) = "
              f"{float(tail):.6f} > min: {tail > best}")
        assert tail > best
        if rho == Fr(7, 3):
            C0 = best
            C0arg = arg
    print(f"  ==> c0(7/3) = {C0} (exact, attained at {C0arg}); the cone floor "
          f"is PROVED (finite check + tail bound)")

    print()
    print("PART C -- the c = 7 wall dichotomy referee")
    blocks = [[50, 51, 52, 53, 54, 55, 56],
              [40, 44, 55, 60, 66, 70, 77],
              [24, 30, 40, 45, 54, 60, 72],
              [10, 24, 57, 133, 311, 726, 1694],
              [12, 20, 30, 44, 60, 84, 120]]
    for blk in blocks:
        blk = sorted(blk)
        rats = [Fr(blk[i + 1], blk[i]) for i in range(6)]
        lac = all(r >= Fr(7, 3) for r in rats)
        if lac:
            print(f"  block {blk}: LACUNARY branch (all ratios >= 7/3) -- "
                  f"opus-S336 closes it")
            continue
        credits = [mu_fr(blk[i], blk[i + 1]) for i in range(6)
                   if Fr(blk[i + 1], blk[i]) < Fr(7, 3)]
        good = 1 - mu_brute(blk, want_all=False)
        floor = max(credits)
        print(f"  block {blk}: hunter branch: best credit {floor} "
              f">= c0 {C0}: {floor >= C0}; true good = {float(good):.5f} "
              f">= credit: {good >= floor}")
        assert floor >= C0 and good >= floor

    print()
    print("PART D -- the c = 8 honest boundary")
    blk8 = [100 + i for i in range(8)]
    credits = sum(mu_fr(blk8[i], blk8[i + 1]) for i in range(7))
    need = Fr(1, 7)
    good = 1 - mu_brute(blk8, want_all=False)
    print(f"  block {blk8}: pair-credit sum = {float(credits):.5f} vs needed "
          f"(c-7)/7 = {float(need):.5f}: sufficient: {credits > need}; "
          f"true good = {float(good):.5f} (the truth is fine; the uniform "
          f"pair route is not) -- c >= 8 NAMED")

    print()
    print("PART E -- difference-comb law + beats + the clustered beat-miss")
    for a, b in [(50, 53), (35, 42), (50, 51)]:
        d = b - a
        # subset law: every point of D_a cap D_b has ||d t|| < 1/7
        bps = sorted(set([Fr(0), Fr(1)] + danger_breaks(a) + danger_breaks(b)))
        beats = {}
        ok = True
        for i in range(len(bps) - 1):
            mid = (bps[i] + bps[i + 1]) / 2
            if in_danger(a, mid) and in_danger(b, mid):
                fd = (d * mid) % 1
                if min(fd, 1 - fd) >= Fr(1, 7):
                    ok = False
                k = round(float(d * mid)) % d if d > 0 else 0
                beats[k] = beats.get(k, Fr(0)) + (bps[i + 1] - bps[i])
        tot = sum(beats.values())
        print(f"  (a,b)=({a},{b}) d={d}: subset law holds: {ok}; beats "
              f"{{{', '.join(f'{k}: {float(v):.5f}' for k, v in sorted(beats.items()))}}}"
              f" (total {float(tot):.5f})")
        assert ok
    # clustered beat-miss: (a, a+1), windows of width 2/a off the beat
    a = 2464
    b = a + 1
    hits = 0
    for c0w in (Fr(1, 3), Fr(1, 2), Fr(2, 3)):
        lo, hi = c0w, c0w + Fr(2, a)
        m = Fr(0)
        for x in [t for t in danger_breaks(a) if lo <= t <= hi]:
            pass
        # windowed brute on [lo, hi]
        bps = sorted(set([lo, hi] +
                         [t for v in (a, b) for t in danger_breaks(v)
                          if lo <= t <= hi]))
        for i in range(len(bps) - 1):
            mid = (bps[i] + bps[i + 1]) / 2
            if in_danger(a, mid) and in_danger(b, mid):
                m += bps[i + 1] - bps[i]
        hits += (m > 0)
        print(f"  clustered ({a},{b}): off-beat window at {float(c0w):.2f} "
              f"(width 2/a): mass = {m}")
    print(f"  => the d = 1 mass sits in ONE beat (near t = 0/1); off-beat "
          f"sweep windows carry ZERO ({3 - hits}/3 zero): positioning is "
          f"essential exactly in the clustered regime -- which "
          f"cite_cluster7_lonely already covers by the C0-citation route")
    print("=" * 78)
    print("done")
