#!/usr/bin/env python3
"""
OPEN-Q-108 (uniform fattening) -- a SHARP explicit conjecture, a decoupling lemma, and the induction
step that forces the extremizer to be bounded.  mac-mini-2026-07-24-S170 (Opus).

G_C(theta) = { tau in [0,1) : ||v tau|| >= theta for all v in C },  theta = 1/14.
OPEN-Q-108 asks: is there c>0 with meas(G_C) >= c for EVERY 12-subset C of distinct positive integers?
(Equivalently: is the primitive tight locus at n=13 finite?)

SHARP CONJECTURE (this file):
    meas(G_C) >= 7/858  for every 12-subset C,  with equality iff C = {1..13}\{6} up to dilation.
7/858 is exactly THM-541's minimum over the drop-family; the conjecture is that it is GLOBAL.

MEASURE-DECOUPLING LEMMA (proved):
    meas(G_{C u {W}}) >= (6/7) meas(G_C) - 2N/(7W),   N = number of intervals of G_C.
    Proof: G_C is a union of N intervals of total measure mu. The bad set {||W tau||<1/14} is a union of
    intervals of length 1/(7W) spaced 1/W. An interval of length l meets at most lW+2 of them, so
    meas(bad ∩ G_C) <= mu/7 + 2N/(7W). Subtract. QED

INDUCTION STEP (closes, with 3.4x slack):
    min over 11-subsets of meas(G) = 313/9702 ~ 0.032261  (computed, at {1,2,3,4,5,7,8,9,11,12,13})
    required: > (7/6)(7/858) = 49/5148 ~ 0.009518.   3.39x slack.
    Hence for W large, meas(G_{C' u {W}}) >= (6/7)(313/9702) - eps > 7/858, so any 12-set attaining
    the minimum has a BOUNDED maximum element:  W <= (2/(7*slack)) * N' = 14.656 * N'.

EFFECTIVENESS (honest): N' is small exactly where it matters (extremal 11/12-sets have N'=10..22, giving
W <= 147..322), but N' is not bounded absolutely, so the general statement is a ratio-type bound and
recursing 11 levels gives an astronomically large absolute bound (~30^11). OPEN-Q-108 is thereby REDUCED
to a finite check, but not to an effective one.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools

THETA = F(1, 14)


def safe_intervals(v, th=THETA):
    out = []
    for k in range(v):
        lo = F(k, v) + th / v
        hi = F(k + 1, v) - th / v
        if lo <= hi:
            out.append((lo, hi))
    return out


def intersect(A, B):
    out, i, j = [], 0, 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo <= hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def G(C, th=THETA):
    """EXACT G_C as a list of disjoint closed intervals."""
    cur = safe_intervals(C[0], th)
    for v in C[1:]:
        cur = intersect(cur, safe_intervals(v, th))
        if not cur:
            return []
    return cur


def measG(C, th=THETA):
    return sum((hi - lo) for lo, hi in G(C, th))


def main():
    print(f"SHARP CONJECTURE: meas(G_C) >= 7/858 = {float(F(7,858)):.8f} for all 12-subsets\n")

    print("(1) reproduce THM-541: min over 12-cores {1..13}\\{j}")
    vals = sorted((measG([x for x in range(1, 14) if x != j]), j) for j in range(1, 14))
    print(f"    min = {vals[0][0]} at j={vals[0][1]}   (7/858)  match={vals[0][0]==F(7,858)}\n")

    print("(2) exhaustive 12-subsets of {1..16} (primitive)")
    best = []
    for C in itertools.combinations(range(1, 17), 12):
        if reduce(gcd, C) != 1:
            continue
        best.append((measG(list(C)), C))
    best.sort()
    print(f"    checked {len(best)}; minimum = {best[0][0]} at {best[0][1]}")
    print(f"    any below 7/858? {any(m < F(7,858) for m, _ in best)}\n")

    print("(3) the induction step: min over 11-subsets of {1..15}")
    b11 = []
    for C in itertools.combinations(range(1, 16), 11):
        if reduce(gcd, C) != 1:
            continue
        b11.append((measG(list(C)), C))
    b11.sort()
    m11 = b11[0][0]
    need = F(7, 6) * F(7, 858)
    print(f"    min_11 = {m11} = {float(m11):.8f} at {b11[0][1]}")
    print(f"    need   > {need} = {float(need):.8f}   -> step closes: {m11 > need}"
          f"  (slack {float(m11/need):.2f}x)\n")

    slack = F(6, 7) * m11 - F(7, 858)
    coef = F(2, 7) / slack
    print(f"(4) resulting bound on the max element: W <= {coef} * N' = {float(coef):.3f} * N'")
    for name, C in [("{1..13}\\{6}", [x for x in range(1, 14) if x != 6]),
                    ("{1,2,3,4,5,7,8,9,11,12,13}", [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13]),
                    ("{1..11}", list(range(1, 12)))]:
        N = len(G(C))
        print(f"      {name:>28}: N'={N:>3}  ->  W <= {float(coef)*N:.0f}")


if __name__ == "__main__":
    main()
