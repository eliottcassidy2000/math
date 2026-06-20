#!/usr/bin/env python3
"""INDEPENDENT adversarial verification of the finite Newton/Mobius peel identity.

CLAIM under test (attempt to REFUTE):
  Let B be a "base" set of runners and F a set of "far" runners, all positive
  integers (speeds). Define p0(A) = the survival measure described below.

  Finite Newton expansion:
     p0(B u F) = sum_{S subset F} Delta_S(B),
  where the finite Mobius / forward-difference operator is
     Delta_S(B) = sum_{T subset S} (-1)^{|S|-|T|} p0(B u T).
  AND Delta_S(B) = 0 whenever |S| > 6.

  Rationale offered: B misses at most 6 sectors, so more than 6 far runners
  cannot all be *jointly* needed -> the |S|>6 differences vanish.

This script writes its OWN exact p0 from scratch (Fraction arithmetic, no import
of the project's existing implementation) and tests the identity term by term.

----------------------------------------------------------------------------
THE MODEL (seven-sector LRC(14) survival, matching the project canon):

  Circle [0,1) of phases t. Sectors are the 7 half-open intervals
  [j/7, (j+1)/7) for j=0..6. The CENTRAL forbidden sector for a runner is the
  one straddling 0, i.e. the runner c "dies" at phase t iff {c t} lands in the
  forbidden window of width 2/14 = 1/7 centered at 0, i.e. ||c t|| < 1/14
  (equivalently {c t} in [0,1/14) U (13/14,1)).

  A runner c SURVIVES the phase t iff ||c t|| >= 1/14.

  p0(A) = meas{ t in [0,1) : every c in A survives },  exact Fraction.

  This is EXACTLY the safe-measure object meas(G_A) of
  04-computation/lrc14_core_gap_survival_bridge_codex_s36.py (danger arcs of
  half-width 1/14 around each multiple of 1/c). It is implemented here
  independently so the verification does not borrow that code.

  Runner c=0 is degenerate (constant phase), so for the model to be the genuine
  seven-sector object we use POSITIVE runners. To honor the requested literal
  inputs B=(0,1,2,3,4,5) we ALSO test by shifting to positive runners
  (1,2,3,4,5,6) -- a runner labelled 0 has {0*t}=0 < 1/14 for ALL t, killing
  every phase and giving p0==0 identically, which is a degenerate (but still
  valid) instance of the identity. Both are reported.
----------------------------------------------------------------------------
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


# ---- exact survival measure p0 (written from scratch) ----------------------

HALF = Fraction(1, 14)  # danger half-width: ||c t|| < 1/14 is fatal


def _danger_arcs(c: int):
    """Arcs of t in [0,1) where runner c dies: ||c t|| < 1/14.

    {c t} in (a - 1/14, a + 1/14) for each integer a, i.e.
    t in (a/c - 1/(14c), a/c + 1/(14c)) for a = 0..c-1, wrapped to [0,1).
    Returns a list of (lo, hi) sub-intervals already clipped/wrapped to [0,1).
    """
    if c == 0:
        # {0*t} = 0 < 1/14 for every t -> entire circle is fatal.
        return [(Fraction(0), Fraction(1))]
    arcs = []
    d = c  # scale
    hw = Fraction(1, 14 * d)  # half width in t around each a/c
    for a in range(c):
        center = Fraction(a, c)
        lo = center - hw
        hi = center + hw
        # wrap into [0,1)
        if lo < 0:
            arcs.append((Fraction(0), hi))
            arcs.append((lo + 1, Fraction(1)))
        elif hi > 1:
            arcs.append((lo, Fraction(1)))
            arcs.append((Fraction(0), hi - 1))
        else:
            arcs.append((lo, hi))
    return arcs


def _merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append([a, b])
    return out


def p0(runners) -> Fraction:
    """Exact measure of phases where EVERY runner survives (||c t||>=1/14)."""
    danger = []
    for c in runners:
        danger.extend(_danger_arcs(c))
    merged = _merge(danger)
    bad = sum((b - a for a, b in merged), Fraction(0))
    return Fraction(1) - bad


# ---- finite Mobius difference operator -------------------------------------

def delta_S(B, S) -> Fraction:
    """Delta_S(B) = sum_{T subset S} (-1)^{|S|-|T|} p0(B u T)."""
    S = list(S)
    total = Fraction(0)
    for k in range(len(S) + 1):
        sign = (-1) ** (len(S) - k)
        for T in combinations(S, k):
            total += sign * p0(tuple(B) + T)
    return total


def newton_sum(B, F) -> Fraction:
    """sum_{S subset F} Delta_S(B), over ALL subsets S of F."""
    F = list(F)
    total = Fraction(0)
    for k in range(len(F) + 1):
        for S in combinations(F, k):
            total += delta_S(B, S)
    return total


# ---- experiments -----------------------------------------------------------

def report_pair(B, F, label):
    print(f"=== {label} ===")
    print(f"B = {B}, F = {F}  (r = |F| = {len(F)})")
    direct = p0(tuple(B) + tuple(F))
    via_newton = newton_sum(B, F)
    print(f"  p0(B u F) direct        = {direct} = {float(direct):.12f}")
    print(f"  sum_S Delta_S(B)        = {via_newton} = {float(via_newton):.12f}")
    equal = direct == via_newton
    print(f"  EXACTLY EQUAL?          = {equal}")
    # per-size breakdown of Delta contributions
    F = list(F)
    print("  Delta breakdown by |S|:")
    for k in range(len(F) + 1):
        ssum = Fraction(0)
        worst = Fraction(0)
        for S in combinations(F, k):
            d = delta_S(B, S)
            ssum += d
            if abs(d) > abs(worst):
                worst = d
        print(f"    |S|={k}: sum over subsets = {ssum}"
              f"   (max |Delta_S| at this size = {abs(worst)})")
    print()
    return equal, direct, via_newton


def report_top_vanishing(B, F):
    """Confirm Delta_S = 0 for the TOP-size subset (|S|=|F|), the |S|>6 claim."""
    r = len(F)
    print(f"--- top-order difference test: B={B}, F={F}, |S|=r={r} ---")
    d_full = delta_S(B, tuple(F))
    print(f"  Delta_F(B) (|S|={r}) = {d_full}  ->  zero? {d_full == 0}")
    print()
    return d_full


def main():
    print("INDEPENDENT Newton-peel verification (exact Fraction arithmetic)\n")

    # ---- requested case 1: literal B=(0,1,2,3,4,5), F=(17,19,23) ----
    # Note: runner 0 is degenerate (kills all phases), so p0 == 0 identically.
    report_pair((0, 1, 2, 3, 4, 5), (17, 19, 23),
                "requested literal: B=(0,1,2,3,4,5), F=(17,19,23)")

    # ---- non-degenerate analog: shift base to positive runners ----
    report_pair((1, 2, 3, 4, 5, 6), (17, 19, 23),
                "non-degenerate base: B=(1,2,3,4,5,6), F=(17,19,23)")

    # A genuinely "6-missed-sector" base with a real positive zero-survival
    # structure: use a base of 6 positive runners and 3 dissociated far ones.
    report_pair((1, 2, 3, 4, 5, 7), (17, 19, 23),
                "alt base with a hole: B=(1,2,3,4,5,7), F=(17,19,23)")

    # ---- r = 7 far set: confirm the |S|=7 term Delta_F is 0 ----
    print("################ r=7 top-order vanishing tests ################\n")
    # 7 dissociated far values (well-separated, mutually coprime-ish, far from base)
    F7 = (17, 19, 23, 29, 31, 37, 41)
    report_top_vanishing((1, 2, 3, 4, 5, 6), F7)
    report_top_vanishing((0, 1, 2, 3, 4, 5), F7)
    report_top_vanishing((1, 2, 3, 4, 5, 7), F7)

    # Also confirm full identity still holds for an r=7 case (sum over ALL S).
    report_pair((1, 2, 3, 4, 5, 6), F7,
                "FULL identity at r=7: B=(1,2,3,4,5,6), F=7 dissociated far")

    # ---- ADVERSARIAL: does |S|>6 REALLY always vanish? Try to break it. ----
    print("################ adversarial |S|>6 vanishing probes ################\n")
    # Make far runners NOT dissociated / small, base with FEWER than 6 elements,
    # base with many holes -- try to manufacture a nonzero |S|=7 difference.
    probes = [
        ((1, 2), (3, 5, 7, 9, 11, 13, 15)),          # tiny base, structured far
        ((1,), (2, 3, 4, 5, 6, 7, 8)),                # base size 1
        ((1, 2, 3, 4, 5, 6, 7, 8), (9, 10, 11, 12, 13, 14, 15)),  # big base
        ((2, 3, 5, 7, 11), (13, 17, 19, 23, 29, 31, 37)),  # primes
    ]
    any_nonzero = False
    for B, F in probes:
        d = delta_S(B, tuple(F))
        flag = "" if d == 0 else "  <<< NONZERO!"
        print(f"  B={B}, |F|={len(F)} F={F}: Delta_F = {d}{flag}")
        if d != 0:
            any_nonzero = True
    print()
    print(f"Any |S|=7 top-difference nonzero among adversarial probes? {any_nonzero}")

    # ---- ADVERSARIAL: |S|=7 vs |S|=6 -- is 6 really the threshold? ----
    print("\n################ is 6 the exact threshold? |S|=6 nonzero? ################\n")
    # If the claim is right we expect SOME |S|=6 difference to be nonzero
    # (otherwise the threshold would be lower than 6).
    B6 = (1, 2, 3, 4, 5, 6)
    F6 = (17, 19, 23, 29, 31, 37)
    d6 = delta_S(B6, F6)
    print(f"  B={B6}, |S|=6 F={F6}: Delta_S = {d6}  -> nonzero? {d6 != 0}")
    # and a base SMALLER than 6 -- threshold may drop with base size
    for bsize in range(1, 8):
        B = tuple(range(1, bsize + 1))
        # use |S| = (number of missed sectors could be up to 6)
        for ssize in (6, 7):
            F = tuple(range(101, 101 + ssize))  # large dissociated far runners
            d = delta_S(B, F)
            print(f"  |B|={bsize} B={B}, |S|={ssize}: Delta_S nonzero? {d != 0}"
                  f"  value={d}")


if __name__ == "__main__":
    main()
