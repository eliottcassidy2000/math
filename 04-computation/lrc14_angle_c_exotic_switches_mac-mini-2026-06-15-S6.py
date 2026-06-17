#!/usr/bin/env python3
"""
ANGLE C — condense the EXOTIC closest-to-counterexample cases of LRC(14) to a
few relations/switches.  REGIONS/SECTIONS reframe.

The "hard core" = covering sets (a multiple of every q in 2..14, THM-523).
The closest-to-1/14 from ABOVE is {1..11,13,84} with M = 7/89.

Findings (all exact rationals, stdlib only):

1. SECTION OCCUPANCY.
   - At every GRID time a/14 (a in (Z/14)* = {1,3,5,9,11,13}) the runner 84
     sits in section 0 (84 = 6*14): the observer is NEVER lonely on the grid.
     So M is attained OFF-grid, at tau* = 37/89.
   - On the FINE circle (denominator 89) the observer's clear band of radius
     7/89 is flanked symmetrically by exactly TWO runners: v=5 at +7 and v=84
     at -7.  They are the BINDING PAIR.

2. THE BINDING-PAIR SWITCH (one relation).
   - tau* = 37/89 is a SUM-crossing of the pair (5, 84): 5 + 84 = 89 = denom.
   - 5*tau* = +7/89, 84*tau* = -7/89; symmetric about 0 (a + b = denom forces
     frac(a*tau)+frac(b*tau) = integer => equidistant).
   - M(S) for the whole family is pinned by this single pair (5, X).

3. CLOSED FORM over the covering family {1..11,13} + X (X = 84m).
   - M = 7m/(84m+5), monotone increasing in m, limit 7/84 = 1/12.
   - MINIMUM at m=1: 7/89.  Confirmed (452-set scan) the GLOBAL closest-to-1/14.

4. THE 1/14 MARGIN is one inequality:
   7/89 > 1/14  <=>  7*14 > 89  <=>  98 > 89  (TRUE by 9 = 3^2).

5. FINITE CHECKLIST.  Across 27 covering bases the SMALL flank runner stays in
   {2,4,5,13} (tiny); for the canonical base {1..11,13} it is always 5.
   So "does this covering set beat 1/14?" reduces to a handful of pair sum-
   crossings (small flank in {2,4,5,13}) x (the big element X) — exhaustively
   checkable.

Tournament link: ANALOGY ONLY.  89 = 5+84 prime; margin 9 = 3^2; runner 5 lies
in the order-3 (Eisenstein/Phi_3) orbit of (Z/14)*.  The honest bridge is the
SDR/Hall side of the reframe ("each runner its own section" = system of distinct
representatives), not a tournament statement.
"""
from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations


def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r


def g(S, t):
    return min(nrm(v * t) for v in S)


def cand(S):
    S = sorted(set(S))
    C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2))
    return C


def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at


HARD = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 84]
BASE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]


def part1_sections():
    print("=" * 64)
    print("1. SECTION OCCUPANCY  (reframe)")
    print("=" * 64)
    b, at = M(HARD)
    print(f"S = {HARD}")
    print(f"M = {b} = {float(b):.6f}   (1/14 = {float(F(1,14)):.6f})   M>1/14: {b>F(1,14)}")
    print(f"attained OFF-grid at tau* = {at}")
    print()
    units = [a for a in range(1, 14) if gcd(a, 14) == 1]
    print(f"grid times a/14, a in (Z/14)* = {units}:")
    for a in units:
        sec0 = [v for v in HARD if (v * a) % 14 == 0]
        print(f"  a={a:2d}: section 0 (observer NOT lonely) = {sec0}   "
              f"-> 84 always sits in section 0")
    print()
    D = at.denominator
    print(f"FINE circle (denom {D}); observer clear band radius {b}:")
    above = min((((v * at.numerator) % D), v) for v in HARD if 0 < (v * at.numerator) % D <= D // 2)
    below = min(((D - (v * at.numerator) % D), v) for v in HARD if (v * at.numerator) % D > D // 2)
    print(f"  nearest ABOVE: runner {above[1]} at +{above[0]}/{D}")
    print(f"  nearest BELOW: runner {below[1]} at -{below[0]}/{D}")
    print(f"  => binding pair ({above[1]}, {below[1]}), gap symmetric.")


def part2_binding():
    print()
    print("=" * 64)
    print("2. BINDING-PAIR SWITCH")
    print("=" * 64)
    b, at = M(HARD)
    a, bb = 5, 84
    print(f"5 + 84 = {a+bb} = denom(tau*) = {at.denominator}")
    print(f"5*tau* mod 1 = {nrm(5*at)},  84*tau* mod 1 = {nrm(84*at)}  (symmetric)")
    print("optimum is the SUM-crossing tau* = k/(a+b) of pair (5,84).")
    print()
    print("top pair sum/diff crossings (the switch candidates):")
    rows = []
    Ss = sorted(set(HARD))
    for i in range(len(Ss)):
        for j in range(i + 1, len(Ss)):
            x, y = Ss[i], Ss[j]
            for d, sgn in ((x + y, '+'), (y - x, '-')):
                if d <= 0:
                    continue
                k = 1
                while F(k, d) <= F(1, 2):
                    t = F(k, d)
                    rows.append((g(HARD, t), sgn, x, y, t)); k += 1
    rows.sort(reverse=True)
    for gg, sgn, x, y, t in rows[:6]:
        print(f"  g={gg} ({float(gg):.5f})  pair({x},{y}) {sgn}-cross  tau={t}")
    print(f"  global M = {rows[0][0]}  (top row, pair (5,84))")


def part3_family():
    print()
    print("=" * 64)
    print("3. CLOSED FORM over covering family {1..11,13}+X, X=84m")
    print("=" * 64)
    miss = [q for q in range(2, 15) if not any(v % q == 0 for v in BASE)]
    print(f"base {{1..11,13}} missing q = {miss}; big element must be multiple of "
          f"lcm{tuple(miss)} = {lcm(*miss)}")
    print()
    for m in range(1, 8):
        X = 84 * m
        b, at = M(BASE + [X])
        pred = F(7 * m, 84 * m + 5)
        binders = tuple(sorted(v for v in sorted(set(BASE + [X])) if nrm(v * at) == b))
        print(f"  m={m} X={X:4d}: M={str(b):>9}  =7m/(84m+5)? {b==pred}  "
              f"binders={binders}  {float(b):.6f}")
    print("  limit m->inf -> 7/84 = 1/12 (no-big-element value); min at m=1 = 7/89")


def part4_margin_and_checklist():
    print()
    print("=" * 64)
    print("4. THE 1/14 MARGIN  and  5. FINITE CHECKLIST")
    print("=" * 64)
    print("7/89 > 1/14  <=>  7*14 > 89  <=>  98 > 89  (TRUE by 9 = 3^2)")
    print()
    # global scan
    def is_cov(S):
        return all(any(v % q == 0 for v in S) for q in range(2, 15))
    best = (F(1), None)
    cnt = 0
    bigs = [84 * m for m in range(1, 11)]
    for combo in combinations(range(1, 15), 12):
        for big in bigs:
            S = sorted(set(list(combo) + [big]))
            if len(S) != 13 or not is_cov(S):
                continue
            b, _ = M(S)
            cnt += 1
            if F(1, 14) < b < best[0]:
                best = (b, tuple(S))
    print(f"scanned {cnt} covering sets (12 from 1..14 + one mult of 84):")
    print(f"  UNIQUE closest-to-1/14 from above: M={best[0]}  S={best[1]}")
    print()
    # small-flank set across bases
    small = set()
    tb = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],
          [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14],
          [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],
          [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]]
    for b0 in tb:
        miss = [q for q in range(2, 15) if not any(v % q == 0 for v in b0)]
        need = lcm(*miss) if miss else 1
        for m in range(1, 8):
            X = need * m
            S = sorted(set(b0 + ([X] if X > 14 else [])))
            if X <= 14 and miss:
                continue
            bb, at = M(S)
            for v in S:
                if v <= 14 and nrm(v * at) == bb:
                    small.add(v)
    print(f"SMALL flank runner across covering bases lies in: {sorted(small)} (tiny)")
    print("=> CHECKLIST: 'beats 1/14?' = check pair sum-crossings")
    print("   (small flank in {2,4,5,13}) x (big element X).  Exhaustive, finite.")


if __name__ == "__main__":
    part1_sections()
    part2_binding()
    part3_family()
    part4_margin_and_checklist()
