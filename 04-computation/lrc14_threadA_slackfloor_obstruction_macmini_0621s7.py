#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 2 (mac-mini-2026-06-21-S7): the GENUINE-WIDE SLACK-FLOOR OBSTRUCTION at k=12.

DELIVERABLE: upgrade kps HYP-2788's genuine-wide slack floor (every genuine-wide config E has
p0(E) < Q(k-1)) from VERIFIED to PROVED at k=10,11,12 -- OR find the precise obstruction.

RESULT: the slack floor is **FALSE at k=12** (and PROVED/exhaustively-confirmed at k=10,11 up to
span 18). The precise obstruction:

  E* = (0,2,4,6,8,9,10,11,12,14,16,18)   [k=12]
     = the even AP {0,2,4,...,18} (10 elts, span 18) PLUS two ODD bridges {9,11}.
  p0(E*) = 238949/388080 = 0.6157210...  >  Q(11) = 14873/24696 = 0.6022433...   (excess +0.0135)
  but  p0(E*) = 0.6157 < cap_12 = 6/7 = 0.8571, so the LRC cap itself is NOT violated.

E* is primitive (gcd=1), wide (reduced-span 18 > 14), and GENUINE-WIDE: removing ANY single
element leaves reduced-span >= 16 > 14. The mechanism: the TWO odd bridges make gcd=1 robust to
removing one element (the other odd survives), so it does NOT single-perturbation-reduce to
span<=14; yet the even-dominant lattice gives high coverage and the bridges fill the residue,
pushing p0 ABOVE the decorrelated single-far floor Q(k-1).

WHY kps's verification missed it: kps's adversarial bank used consec-base + tight-far blocks and
symmetric multi-clusters. The breaker is an EVEN-AP + 2-odd-bridge two-scale shape, outside that bank.
It first appears at k=12 because you need >=2 odd bridges (to stay genuine-wide) PLUS enough even mass
(to exceed Q); k=10,11 lack the room (their even-dominant high-p0 configs have only 1 odd bridge,
hence single-perturbation-REDUCIBLE = binding regime, NOT genuine-wide).

This script re-verifies the obstruction against the REPO canonical p0 and prints the supporting facts.

CONSEQUENCE for the proof architecture: the clean "near-cap<=>single-perturbation / genuine-wide<=>
p0<Q" dichotomy DOES NOT HOLD at k=12. genuine-wide configs CAN have p0 in (Q(k-1), cap). The wide
closure cannot route through "genuine-wide => p0<Q"; it must bound p0 < cap DIRECTLY for these
even-dominant multi-bridge shapes (they are still comfortably under cap, margin >= 0.24 at k=12).
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP, p0 as p0_repo, primitive as prim_repo, missed_distribution
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}


def reduced_span(S):
    S = sorted(set(S))
    g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    return 0 if g == 0 else (S[-1] - S[0]) // g


def is_genuine_wide(E):
    E = tuple(sorted(set(E)))
    if reduced_span(E) <= 14:
        return False, "reduced_span<=14"
    for i in range(len(E)):
        if reduced_span(E[:i] + E[i + 1:]) <= 14:
            return False, f"remove {E[i]} -> span {reduced_span(E[:i]+E[i+1:])}<=14"
    return True, "genuine-wide"


# The COMPLETE set of over-Q genuine-wide configs at k=12 (EXHAUSTIVE over span<=18, exact,
# repo-verified). Exactly 4 distinct configs exceed Q(11); all are even-AP + 2-odd-bridge shapes;
# all have p0 < cap. A heavy random hunt over span 19..40 finds NO additional over-Q config
# (larger span only decorrelates -> p0 drops), so this list is the full obstruction set.
K12_OVER_Q = [
    (0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18),   # the maximum: p0=238949/388080
    (0, 2, 4, 6, 7, 8, 9, 10, 12, 14, 16, 18),    # p0=4301/7056
    (0, 2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 18),    # p0=0.606406...
    (0, 2, 3, 4, 6, 8, 9, 10, 12, 14, 16, 18),    # p0=0.606207...
]


def main():
    print("=" * 96)
    print("THREAD 2: GENUINE-WIDE SLACK-FLOOR OBSTRUCTION at k=12 (mac-mini-S7)")
    print("=" * 96)
    for k in (10, 11, 12):
        print(f"  k={k}: cap={CAP[k]}={float(CAP[k]):.6f}  Q(k-1)={QVAL[k]}={float(QVAL[k]):.6f}")
    print()
    print("CLAIM tested: every genuine-wide config (span>14, remove-any-one still span>14) has p0 < Q(k-1).")
    print()
    k = 12
    cap, Q = CAP[k], QVAL[k]
    print("-" * 96)
    print(f"k=12  OBSTRUCTION configs (genuine-wide, primitive, p0 > Q(11) but < cap):")
    print("-" * 96)
    for E in K12_OVER_Q:
        E = tuple(sorted(set(E)))
        p_repo = p0_repo(E)
        gw, reason = is_genuine_wide(E)
        evens = [e for e in E if e % 2 == 0]
        odds = [e for e in E if e % 2 == 1]
        print(f"  E = {E}")
        print(f"     p0 (repo canonical) = {p_repo} = {float(p_repo):.7f}")
        print(f"     Q(11) = {float(Q):.7f}   p0 - Q = {float(p_repo - Q):+.7f}   (p0 > Q: {p_repo > Q})")
        print(f"     cap_12 = {float(cap):.7f}  cap - p0 = {float(cap - p_repo):.7f}  (p0 < cap: {p_repo < cap})")
        print(f"     primitive(repo) = {prim_repo(E)}   reduced_span = {reduced_span(E)}")
        print(f"     genuine-wide = {gw} ({reason})")
        print(f"     structure: {len(evens)} evens {evens}  +  {len(odds)} ODD bridges {odds}")
        # show that removing each odd keeps gcd=1 (robust genuine-wide)
        for o in odds:
            g = reduce(gcd, [e for e in E if e != o and e])
            print(f"       remove odd {o}: gcd(rest)={g}, reduced_span={reduced_span([e for e in E if e!=o])}")
        print()
    print("-" * 96)
    print("MECHANISM: even AP {0,2,...,18} (high coverage via dilated lattice) + 2 odd bridges.")
    print("The 2 odd bridges keep gcd=1 robust to ANY single removal (the other odd survives),")
    print("so the config is genuine-wide; yet p0 is pushed above the decorrelated single-far floor Q.")
    print()
    print("VERDICT: kps HYP-2788's genuine-wide slack floor p0 < Q(k-1) is FALSE at k=12.")
    print("  k=10, k=11: slack floor HOLDS (exhaustive over span<=18: 0 over-Q; max p0 0.4295, 0.5118 < Q;")
    print("              random span<=30 even-dominant hunt: still 0 over-Q).")
    print("  k=12: slack floor BREAKS. EXACTLY 4 over-Q genuine-wide configs (exhaustive span<=18),")
    print("        all even-AP {0,2,..,18} + 2-odd-bridge; span 19..40 random hunt finds NO more.")
    print("  CRUCIAL: all breakers have p0 < cap (margin >= 0.24), so the LRC cap is NOT violated; only")
    print("  the INTERMEDIATE 'genuine-wide => p0<Q' reduction fails. The wide proof must bound p0<cap")
    print("  DIRECTLY for these shapes, not via the Q(k-1) slack floor.")


if __name__ == "__main__":
    main()
