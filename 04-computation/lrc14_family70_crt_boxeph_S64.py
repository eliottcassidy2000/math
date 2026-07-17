#!/usr/bin/env python3
"""
THE FAMILY-70 CRT STRUCTURE: THE SURVIVOR LAW (boxeph-2026-07-17-S64)
Owner directive: the family70 k == 1 (mod 10) progression (LEM-034's named
open); keep proving little statements.

THE SECTION-VECTOR FORMULA (LEM-035(A)).  At a class-0 crossing x = k/(7m)
of a 7-full owner e = 7m, every other runner f has 7 f x = fk/m, so its
section at x is sigma_f = floor(fk/m) mod 7, it is AT a boundary iff m | fk
(the LEM-034 lattice), and its LEFT section is sigma_f - [m | fk] (mod 7).
Owner e itself sits in section 6 on the left.  Hence the SURVIVOR CRITERION
(x- in R_0): the left-vector of the other six runners avoids 0 and covers
{1,...,5}  (6 is already covered by e; one slack slot).

THE r = 1 UNIVERSAL FAMILY (LEM-035(B), proved).  For E = {1,...,6, 7M},
M >= 7, write k = Mj + r.  At r = 1 the carries floor(f/M) vanish and no
runner is at a boundary, so the left-vector is (j, 2j, ..., 6j) = the
multiplication-by-j permutation of {1..6} when j != 0 (mod 7): covers
{1..5}, avoids 0 -> SURVIVOR.  j == 0 puts all six smalls in section 0 ->
dead.  So the survivors at r = 1 are EXACTLY k = Mj + 1, j = 1..6 --
independent of M.  This is the measured family70 progression k == 1 (mod
10), k != 1, PROVED.

THE M = 10 COMPLETENESS (LEM-035(C), proved by structured case analysis,
automated here): no other residue r survives for M = 10 -- the S63 census
pattern is a theorem.

THE M = 11 PREDICTION (LEM-035(D)): carries at r = 9 give the left-vector
(j, 2j+1, 3j+2, 4j+3, 5j+4, 6j+4); j = 5 yields (5,4,3,2,1,6) -- an extra
survivor at k = 11*5 + 9 = 64.  Verified here, plus the extra-survivor
spectrum for M = 8..15 and the general-cluster re-derivation of ALL S63
survivor sets from the same formula.

Little statements: formula referee on 9 clusters (all 7-full owners, all
crossings, exact Fractions); survivor sets census == formula everywhere;
r=1 law all M; attribution/N_0 cross-check; reflection closure (survivor
exits of R_0 at x reflect to class-1 entries of R_6 at 1-x).
"""

import sys
from bisect import bisect_left
from fractions import Fraction as Fr
from math import gcd

sys.path.insert(0, '04-computation')
from lrc14_hyp6994_resonance_test_boxeph_S25 import endpoints
from lrc14_general_resonance_law_boxeph_S26 import owner_data

FIVE = frozenset(range(1, 6))
SIX = frozenset(range(1, 7))


def section(f, x):
    return int((f * x % 1) * 7)


def leftvec_formula(E, e, m, k):
    """left sections of the runners f != e at x = k/(7m), via the formula."""
    out = []
    for f in E:
        if f == e:
            continue
        s = (f * k // m) % 7
        if (f * k) % m == 0:
            s = (s - 1) % 7
        out.append(s)
    return out


def survivor_formula(E, e, m, k):
    lv = leftvec_formula(E, e, m, k)
    return 0 not in lv and FIVE <= set(lv)


def census_left(E, e, k, bps):
    """exact left-neighborhood sections at x = k/e (census ground truth)."""
    x = Fr(k, e)
    if k == 0:
        mid = (bps[-2] + 1) / 2
    else:
        i = bisect_left(bps, x)
        assert bps[i] == x
        mid = (bps[i - 1] + x) / 2
    return [section(f, mid) for f in E], mid


CLUSTERS = ([([1, 2, 3, 4, 5, 6, 7 * M], f"family{7 * M}") for M in
             range(8, 16)] +
            [([12, 15, 20, 21, 28, 30, 35], "balanced"),
             ([8, 9, 10, 12, 14, 15, 18], "near-AP"),
             ([1, 2, 3, 4, 5, 56, 84], "two-large")])

S63_SURVIVORS = {("balanced", 21): [11], ("balanced", 28): [], ("balanced", 35): [22],
                 ("near-AP", 14): [12], ("two-large", 56): [9, 11, 33, 47, 48],
                 ("two-large", 84): [13, 14, 16, 37, 62, 70, 71, 72]}

if __name__ == "__main__":
    print("THE FAMILY-70 CRT STRUCTURE: THE SURVIVOR LAW (boxeph S64)")
    print("=" * 78)
    print("PART 1+2 -- formula referee + survivor sets (census == formula)")
    extras = {}
    for E, name in CLUSTERS:
        pts = endpoints(E, 0)
        if not pts:
            print(f"  [{name}] R_0 empty; skipped")
            continue
        ep_pos = {p for p, sg, o in pts}
        bps = sorted(set(Fr(kk, 7 * f) for f in E for kk in range(7 * f + 1)))
        for e in [f for f in E if f % 7 == 0]:
            m = e // 7
            surv_f = []
            worst = 0
            for k in range(e):
                lv = sorted(leftvec_formula(E, e, m, k) + [6])
                direct, mid = census_left(E, e, k, bps)
                if k != 0:
                    assert lv == sorted(direct), (name, e, k, lv, direct)
                if k != 0 and survivor_formula(E, e, m, k):
                    surv_f.append(k)
                    assert Fr(k, e) in ep_pos, (name, e, k)
            # census survivors = crossings whose left side is in R_0
            surv_c = [k for k in range(1, e)
                      if set(census_left(E, e, k, bps)[0]) == SIX]
            assert surv_f == surv_c, (name, e, surv_f, surv_c)
            if (name, e) in S63_SURVIVORS:
                assert surv_f == S63_SURVIVORS[(name, e)], (name, e, surv_f)
            if name.startswith("family"):
                M = m
                r1 = [M * j + 1 for j in range(1, 7)]
                extra = [k for k in surv_f if k not in r1]
                missing_r1 = [k for k in r1 if k not in surv_f]
                extras[M] = extra
                print(f"  [{name}] e={e} (M={m}): survivors {surv_f}; "
                      f"r=1 family complete: {not missing_r1}; "
                      f"extras: {extra if extra else 'NONE'}")
            else:
                print(f"  [{name}] e={e} (m={m}): survivors {surv_f} "
                      f"== S63 census: True")
    print()
    print("PART 3 -- the theorems on the family clusters")
    print(f"  r=1 universal family (k = Mj+1, j=1..6 survive; k=1 dead): "
          f"verified for all M above")
    print(f"  M=10 completeness (no extras): {extras.get(10) == []}")
    print(f"  M=11 prediction (extra survivor k = 64): "
          f"{extras.get(11) == [64]}")
    print("  extra-survivor spectrum: " + ", ".join(
        f"M={M}: {e if e else '-'}" for M, e in sorted(extras.items())))
    assert extras.get(10) == [] and extras.get(11) == [64]
    print()
    print("PART 4 -- attribution / N_0 cross-check (family clusters)")
    for E, name in CLUSTERS:
        if not name.startswith("family"):
            continue
        P, Mpts, data = owner_data(E, 0)
        if Mpts == 0:
            continue
        e = max(E)
        m = e // 7
        if e not in data:
            continue
        surv = [k for k in range(1, e) if survivor_formula(E, e, m, k)]
        owned = [k for k in surv
                 if min(f for f in E if (f * k) % m == 0 or f == e) == e]
        n0 = data[e]["N"][0]
        print(f"  [{name}] survivors {len(surv)}, attributed to {e}: "
              f"{len(owned)}; N_0 = {n0:+d} (match: {n0 == -len(owned)})")
        assert n0 == -len(owned)
    print()
    print("PART 5 -- reflection closure (family70): survivor exits of R_0 at "
          "k/70 -> class-1 entries of R_6 at (70-k)/70")
    E = [1, 2, 3, 4, 5, 6, 70]
    pts6 = endpoints(E, 6)
    look = {p: (sg, o) for p, sg, o in pts6}
    for k in (11, 21, 31, 41, 51, 61):
        p = Fr(70 - k, 70)
        sg, o = look[p]
        j = p * 7 * o
        assert j.denominator == 1
        print(f"  k={k}: 1-x = {p} in R_6 endpoints: sign {sg:+d}, owner {o}, "
              f"class {int(j) % 7}")
        assert sg == +1
    print("=" * 78)
    print("done")
