#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION of the 'subtorus-relation-lattice' angle on the LRC(14)
uniform floor B(k).

Author: kind-pasteur-2026-06-18-S5-wf (verifier)

The claim under test is marked PARTIAL. It explicitly states:
  - It PROVES a structural backbone (scale-invariance L1, single-circle 1D
    integral, subtorus LIMIT object, finitely many subtorus types).
  - It does NOT close the 2/7 uniform floor B(k).
  - It claims the 2/7 object is REFUTED as a floor mechanism (rho*_{2/7}=0).
  - The live target is the 1/7 weak spread bound; the angle grounds
    'consecutive minimizes' threshold-agnostically via scale-invariance.

GOAL: Re-derive every quoted EXACT number; stress-test each sub-claim; name the
precise surviving gap. We DO NOT take the canon mu_exact on faith -- we build a
RIGOROUS exact engine (adds gap=threshold crossing breakpoints) and cross-check.
"""

from fractions import Fraction as F
from itertools import combinations
from math import comb

# ----------------------------------------------------------------------------
# RIGOROUS exact mu engine.
#
# mu(E) = meas{ x in [0,1) : circular max-gap of { frac(e_i x) } > thr }.
#
# Subtlety: the *order* of the k phase-points can only change at x = t/(e_i-e_j).
# Between two consecutive order-change breakpoints, each point frac(e_i x) is an
# affine function of x with integer slope e_i (mod the wrap), the sorted order is
# fixed, so every circular gap g_p(x) is affine in x. Therefore max-gap is a max
# of affine functions = piecewise-linear CONVEX, and the set {max-gap > thr} on
# the sub-interval is the complement of finitely many points where each affine
# gap crosses thr. We must add those crossing points as breakpoints, then on each
# atomic interval the truth of (max-gap > thr) is constant -> sample the midpoint.
#
# Crossing points: for each ordered pair giving a gap that is affine a*x+b, solve
# a*x+b = thr. But it is cleaner & fully rigorous to just collect ALL candidate
# breakpoints from BOTH sources and sample midpoints of atomic cells:
#   (A) order changes: x = t/(e_i-e_j), t integer
#   (B) gap = thr crossings.
# For (B): on a cell the sorted points are p_(0)<...<p_(k-1) with p_(r)=frac(e_{s(r)} x).
# Each is affine c_r + e_{s(r)} x on the cell. A gap = p_(r+1)-p_(r) (or wrap) is
# affine; set = thr. We don't know s(r) a priori per cell, so instead we add ALL
# x where ANY (e_i x - e_j x) hits an integer +/- thr, i.e. (e_i-e_j) x = t +/- thr
# AND e_i x = t +/- thr (wrap gap involves a single point distance to 0/1). To be
# safe and fully general we collect:
#    for all i<j, all integers t: x=(t)/(ei-ej)            [order change / collision]
#    for all i<j, all integers t: x=(t+thr)/(ei-ej), (t-thr)/(ei-ej)  [pairwise gap=thr]
#    for all i, all integers t:   x=(t+thr)/ei, (t-thr)/ei            [gap to wrap=thr]
# Within [0,1). This OVER-collects (harmless: extra breakpoints only refine cells).
# ----------------------------------------------------------------------------

def mu_rig(E, thr=F(2,7)):
    E = sorted(set(int(e) for e in E))
    assert E[0] == 0 or 0 in E, "0 must be in E"
    bps = set([F(0), F(1)])
    diffs = []
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            d = E[j] - E[i]
            diffs.append(d)
    # pairwise collisions and pairwise gap=thr
    for d in set(diffs):
        if d == 0:
            continue
        # order changes: t/d for t in 0..d
        for t in range(0, d+1):
            bps.add(F(t, d))
        # gap=thr crossings for pairwise differences
        # x where frac((d)x) = thr or 1-thr  => d*x = t+thr or t-thr
        # range of d*x over [0,1) is [0,d)
        t = 0
        # cover t from -1 to d to be safe
        for t in range(-1, d+1):
            for s in (thr, -thr):
                x = (F(t) + s) / d
                if 0 <= x < 1:
                    bps.add(x)
    # gap-to-wrap = thr: frac(e x) = thr or 1-thr for each e
    for e in E:
        if e == 0:
            continue
        for t in range(-1, e+1):
            for s in (thr, -thr):
                x = (F(t) + s) / e
                if 0 <= x < 1:
                    bps.add(x)
    bps = sorted(b for b in bps if 0 <= b < 1)
    tot = F(0)
    for a, b in zip(bps, bps[1:] + [F(1)]):
        if b <= a:
            continue
        mid = (a + b) / 2
        pts = sorted(set((F(e) * mid) % 1 for e in E))
        if len(pts) == 1:
            mg = F(1)
        else:
            gaps = [pts[t+1] - pts[t] for t in range(len(pts)-1)]
            gaps.append(pts[0] + 1 - pts[-1])
            mg = max(gaps)
        if mg > thr:
            tot += (b - a)
    return tot


def mu_consec(k, thr=F(2,7)):
    return mu_rig(list(range(k)), thr)


# ----------------------------------------------------------------------------
# F(k): iid uniform max-gap > 2/7 ceiling, exact.
# P(maxgap > g) = sum_{j>=1} (-1)^{j+1} C(k,j) (1-j*g)_+^{k-1}
# with g = 2/7 (gap threshold). (Standard order-statistics inclusion-exclusion.)
# ----------------------------------------------------------------------------

def Fk(k, g=F(2,7)):
    tot = F(0)
    for j in range(1, k+1):
        base = 1 - j*g
        if base <= 0:
            continue
        tot += F((-1)**(j+1)) * comb(k, j) * base**(k-1)
    return tot


def report(label, val):
    print(f"{label:55s} = {str(val):>22s}  ~= {float(val):.6f}")


if __name__ == "__main__":
    print("="*90)
    print("PART 1: cross-check RIGOROUS engine vs canon a(k) values (2/7 threshold)")
    print("="*90)
    canon_a = {3:F(1), 4:F(19,21), 5:F(9,14), 6:F(4,7), 7:F(83,210),
               8:F(163,490), 9:F(277,980), 10:F(557,2205), 11:F(127,588),
               12:F(1313,6468), 13:F(829,4620)}
    allmatch = True
    for k in range(3, 14):
        v = mu_consec(k)
        ok = (v == canon_a[k])
        allmatch &= ok
        print(f"  a({k:2d})  canon={str(canon_a[k]):>14s}  rig={str(v):>14s}  "
              f"{'MATCH' if ok else '*** MISMATCH ***'}")
    print(f"\n  ALL consecutive a(k) match canon: {allmatch}")

    print("\n" + "="*90)
    print("PART 2: canon mu_min(k) bounded-spread (consecutive is NOT global min for k>=7)")
    print("="*90)
    # canon claims: mu_7=13/35 at (0,2,3,4,5,6,8)
    tests = {
        (0,2,3,4,5,6,8): F(13,35),
        (0,2,3,4,5,6,8,11): F(71,220),
    }
    for E, claimed in tests.items():
        v = mu_rig(list(E))
        print(f"  mu{E} = {str(v):>14s}  claimed {str(claimed):>10s}  "
              f"{'MATCH' if v==claimed else '*** MISMATCH ***'}")

    print("\n" + "="*90)
    print("PART 3: F(k) iid CEILING exact recompute")
    print("="*90)
    for k in (7, 11, 13):
        report(f"F({k})", Fk(k))
