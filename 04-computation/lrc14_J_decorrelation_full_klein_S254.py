# -*- coding: utf-8 -*-
# klein-2026-07-11-S254: THE J-FUNCTIONAL DECORRELATION + EXHAUSTIVE ARGMIN (exact).
#
# GOAL: reduce the k=9 base row (J = 6m1 - m2 >= 432/91 for all 9-cores) to a FINITE check
# by showing large-spread cores are safe (J -> J_iid ~ 7.2 >> floor) and consec is the
# global J-argmin.  This is the moment-ladder route's genuinely-hard remaining item
# (kps cont.36 showed the clean-ruler route is decoupled/finite).
#
#  (A) EXHAUSTIVE ARGMIN: global min J over ALL k-cores in [1..W] -- is it consec?  (Unlike
#      Var(N), which mod-7 beats -- mac-mini cont.37.)  Report argmin + the runner-up shape.
#  (B) DECORRELATION LAW: worst (min) J over cores of spread >= s, as a function of s.
#      If this RISES above the floor for s beyond a threshold s0, then [box: spread < s0]
#      + [decorrelation: spread >= s0] closes the row.  Find s0 and the box size needed.
#  (C) IID LIMIT: J_iid = 6*7*(6/7)^k - 42*(5/7)^k exact -- the dissociated value; the gap
#      J_iid - floor is the decorrelation budget.
#  (D) THREE-GAP MECHANISM at consec: J(consec_k) closed-form check + why it is the min.
#
# Conventions == mac-mini/kps (sectors [c/7,(c+1)/7); N = 7 - #occupied; J = 6m1 - m2).

import sys
from math import gcd
from fractions import Fraction as F
from itertools import combinations

def lcm(a, b):
    return a // gcd(a, b) * b

def J_exact(E):
    nz = [abs(e) for e in E if e != 0]
    has0 = len(nz) < len(E)
    L = 1
    for e in nz:
        L = lcm(L, e)
    D = 7 * L
    pts = set([0, D])
    for e in nz:
        step = L // e
        pts.update(range(0, D + 1, step))
    pts = sorted(pts)
    pn = [0] * 8
    base_hit = 1 if has0 else 0
    for t1, t2 in zip(pts, pts[1:]):
        s = t1 + t2
        hit = base_hit
        for e in nz:
            c = (7 * e * s // (2 * D)) % 7
            hit |= 1 << c
        n = 7 - bin(hit).count("1")
        pn[n] += t2 - t1
    p = [F(x, D) for x in pn]
    m1 = sum(F(j) * p[j] for j in range(8))
    m2 = sum(F(j * (j - 1)) * p[j] for j in range(8))
    return 6 * m1 - m2

def norm(E):
    g = 0
    for e in E:
        if e:
            g = gcd(g, e)
    return tuple(e // g for e in E) if g > 1 else tuple(E)

THRESH = F(432, 91)

def J_iid(k):
    return 6 * 7 * F(6, 7) ** k - 42 * F(5, 7) ** k

def exhaustive_argmin(k, W):
    print(f"=== (A)+(B) k={k}, box [1..{W}], floor 432/91 = {float(THRESH):.5f} ===")
    print(f"    J_iid({k}) = {J_iid(k)} ~ {float(J_iid(k)):.5f}  (decorrelation budget "
          f"{float(J_iid(k) - THRESH):+.4f})")
    seen = set()
    gmin = None; gmin_set = None
    # decorrelation: min J among cores with spread >= s, tabulated
    by_minspread = {}   # spread threshold -> running min J over spread>=that (filled after)
    all_by_spread = {}  # spread -> min J at exactly-ish that spread (bucketed)
    n = 0
    recs = []
    for combo in combinations(range(1, W + 1), k):
        key = norm(combo)
        if key in seen:
            continue
        seen.add(key)
        n += 1
        J = J_exact(combo)
        spread = combo[-1] - combo[0]
        recs.append((J, spread, combo))
        if gmin is None or J < gmin:
            gmin = J; gmin_set = combo
    recs.sort()
    consec = tuple(range(1, k + 1))
    print(f"  {n} normalized cores. GLOBAL min J = {gmin} ~ {float(gmin):.5f} at {gmin_set}"
          f"  {'== CONSEC' if norm(gmin_set)==norm(consec) else '*** NOT CONSEC ***'}")
    print(f"  five lowest-J cores:")
    for J, spread, combo in recs[:5]:
        print(f"    J={float(J):.5f} margin{float(J-THRESH):+.5f} spread={spread}  {combo}")
    # decorrelation: min J over cores with spread >= s
    print(f"  DECORRELATION -- min J among cores with spread >= s:")
    spreads = sorted(set(r[1] for r in recs))
    for s0 in spreads:
        mJ = min((r[0] for r in recs if r[1] >= s0), default=None)
        if mJ is not None:
            flag = "  <-- floor cleared for all spread>=s" if mJ > THRESH else ""
            print(f"    s>={s0:3d}: min J = {float(mJ):.5f} (margin {float(mJ-THRESH):+.5f}){flag}")

def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 22
    exhaustive_argmin(k, W)

if __name__ == "__main__":
    main()
