#!/usr/bin/env python3
"""THREAD F (mac-mini-2026-06-20-SF): verify the LRC(14) S3 reduction target
and pin down the EXACT S3-closing condition.

Questions audited (all exact Fractions, reusing the canonical
missed_distribution / CAP from lrc14_wide_branch_ridge_codex_s47):

  Q1. Is p0(E)=meas(S7(E)) the right object, and does it equal the full-cover
      residual capacity q_E(R={1..6})?  (HYP-2700: |R|=6 layer of the cone.)
  Q2. cap_k values; relation of cap_k to the actual loneliness threshold.
  Q3. HYP-2700b: U4 = p0+p5+5*p6 vs p0.  Does "consec maximizes U4 and U4<=cap"
      close k=8?  (NO, per HYP-2700b.)  Confirm exactly.
  Q4. The S3-closing condition: is it  p0(E)<=cap_k for all primitive E (k=8..12)?
      Or the cone  sum_R w_R q_C(R) <= ... ?  Test: does p0<=cap suffice, i.e.
      is the cone's |R|=6 layer the ONLY layer the LRC predicate needs?
  Q5. iid_k floor = 7!*S(k,7)/7^k and the wide-budget cap_k - iid_k.
  Q6. consec_k = {0..k-1}: p0 (=measS7) <= cap_k with explicit margin.
  Q7. Cone small-|R| violations (HYP-2700a): do they matter for the LRC
      predicate, given the predicate only reads |R|=6?
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import comb, factorial

from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    missed_distribution,
    p0 as p0_func,
    primitive,
    wall_breakpoints,
)

ALL_INNER = 0b1111110


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def stirling2(n: int, k: int) -> int:
    # surjections of n labels onto k boxes / k!  -> S(n,k)
    s = 0
    for j in range(k + 1):
        s += (-1) ** j * comb(k, j) * (k - j) ** n
    return s // factorial(k)


def iid_k(k: int) -> Fraction:
    # probability that k iid uniform Z/7 colors hit ALL of Z/7 = 7! S(k,7)/7^k
    return Fraction(factorial(7) * stirling2(k, 7), 7 ** k)


def residual_capacity(row, R) -> Fraction:
    """q_E(R) = meas{x: image of coloring c(e,x)=floor(7 frac(ex)) over e in E
    contains every inner sector j in R}.  R subset {1..6}.  This is the
    cone coordinate of HYP-2697/2700.  For R={1..6} this is p0=meas(S7)."""
    nonzero = [x for x in row if x]
    if not nonzero:
        return Fraction(0)
    from functools import reduce
    from math import gcd

    def lcm(a, b):
        return a // gcd(a, b) * b

    l = reduce(lcm, nonzero)
    d, bps = wall_breakpoints(row)
    den2 = 2 * l
    Rmask = 0
    for j in R:
        Rmask |= 1 << j
    num = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nonzero:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & Rmask) == Rmask:
            num += hi - lo
    return Fraction(num, d)


def main():
    print("=" * 78)
    print("THREAD F — LRC(14) S3-closing condition audit (exact)")
    print("=" * 78)

    # ---- Q1: p0 = q_E({1..6}) and definition check ----
    print("\n[Q1] p0(E)=meas(S7(E)) equals the full-cover residual capacity q_E({1..6})")
    for E in [tuple(range(8)), tuple(range(9)), (0, 1, 2, 30, 31, 32, 60, 61, 62)]:
        p0v = p0_func(E)
        qfull = residual_capacity(E, (1, 2, 3, 4, 5, 6))
        print(f"  E={E}")
        print(f"    p0={fmt(p0v)}  q_E(R=1..6)={fmt(qfull)}  equal={p0v==qfull}")

    # ---- Q2/Q5: cap_k, iid_k, wide budget ----
    print("\n[Q2/Q5] cap_k, iid_k=7!S(k,7)/7^k, wide budget cap_k - iid_k")
    for k in range(8, 13):
        cap = CAP[k]
        ii = iid_k(k)
        print(f"  k={k}: cap={fmt(cap)}  iid={fmt(ii)}  cap-iid={fmt(cap-ii)}")

    # ---- Q6: consec_k p0 vs cap ----
    print("\n[Q6] consec_k = {0..k-1}: p0=measS7 vs cap_k (the bounded-span extreme)")
    for k in range(8, 13):
        E = tuple(range(k))
        p0v = p0_func(E)
        cap = CAP[k]
        print(f"  k={k}: consec p0={fmt(p0v)}  cap={fmt(cap)}  cap-p0={fmt(cap-p0v)}  "
              f"p0<=cap={p0v<=cap}")

    # ---- Q3: the HYP-2700b U4-vs-p0 discrepancy ----
    print("\n[Q3] HYP-2700b: U4=p0+p5+5p6 OVERSHOOTS cap for consec_8 (U4 cannot close k=8)")
    for k in range(8, 13):
        E = tuple(range(k))
        dist = missed_distribution(E)
        p0v = dist[0]
        U4 = dist[0] + dist[5] + 5 * dist[6]
        cap = CAP[k]
        print(f"  k={k}: consec p0={fmt(p0v)}  U4={fmt(U4)}  cap={fmt(cap)}  "
              f"p0<=cap={p0v<=cap}  U4<=cap={U4<=cap}")

    # ---- Q4: which functional is the LRC predicate?  p0<=cap is the predicate.
    #          The cone sum_R w_R q_C(R) is a PROOF DEVICE for p0, not the predicate.
    print("\n[Q4] The LRC S3 predicate IS p0(E)<=cap_k.  The cone is a route to it.")
    print("     Test: exhaustive bounded span check that consec maximizes p0,")
    print("     so p0(E)<=p0(consec_k)<=cap_k on the bounded box (k=8,9,10).")
    for k, bnd in ((8, 13), (9, 13), (10, 13)):
        capk = CAP[k]
        consec_p0 = p0_func(tuple(range(k)))
        maxp0 = Fraction(0)
        argmax = None
        viol_cap = 0
        viol_consec = 0
        nrows = 0
        for combo in combinations(range(1, bnd + 1), k - 1):
            E = (0,) + combo
            if not primitive(E):
                continue
            nrows += 1
            pv = p0_func(E)
            if pv > maxp0:
                maxp0 = pv
                argmax = E
            if pv > capk:
                viol_cap += 1
            if pv > consec_p0:
                viol_consec += 1
        print(f"  k={k} box max(E)<={bnd}: rows={nrows}")
        print(f"    consec p0={fmt(consec_p0)}  global max p0={fmt(maxp0)} at {argmax}")
        print(f"    rows with p0>cap_k: {viol_cap}   rows with p0>consec_p0: {viol_consec}")

    # ---- Q7: cone small-|R| violations exist but do NOT touch the LRC predicate
    print("\n[Q7] HYP-2700a small-|R| cone violations vs the LRC predicate (|R|=6 only)")
    K8 = tuple(range(8))   # consec
    # a known small-|R| challenger (HYP-2697 style, nonconsecutive)
    C8 = (0, 1, 2, 3, 4, 5, 6, 8)
    print(f"  K=consec_8={K8}   C={C8}")
    for R in [(1, 2, 3, 4, 5, 6), (1, 3, 5), (1, 6), (1,)]:
        qK = residual_capacity(K8, R)
        qC = residual_capacity(C8, R)
        flag = "  <-- C beats K (small-R)" if qC > qK else ""
        print(f"    R={R}: q_K={fmt(qK)}  q_C={fmt(qC)}  K>=C={qK>=qC}{flag}")
    print("  NOTE: the LRC predicate reads ONLY R={1..6} (full cover = p0).")
    print("        Small-|R| violations are internal to the cone PROOF of consec-max,")
    print("        not failures of p0(E)<=cap_k.")


if __name__ == "__main__":
    main()
