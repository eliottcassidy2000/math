#!/usr/bin/env python3
"""
HYP-5853: chain-Bonferroni for the 2-adic E3 diagonal.
boxeph-2026-07-09-S11.  Exact rationals throughout.

A maximal doubling chain {a, 2a, ..., 2^(L-1)·a} ⊆ S is jointly safe at t iff
y = a·t satisfies  ||2^j y|| ≥ 1/14  for j = 0..L-1  — the run-length-limited
(RLL) shift set of the doubling map T(y) = 2y mod 1 relative to the band
B = [1/14, 13/14].  The dangers NEST 2-adically, so

    mu_L := meas( ∩_{j<L} T^{-j}(B) )

is far larger than the union bound 1 - L/7.  Exact computation: T^{-j}(B) is
a union of 2^j intervals with denominators 14·2^j; intersect interval lists
with Fractions.

CHAIN-BONFERRONI: S (13 speeds) decomposes into m maximal chains of lengths
L_1..L_m (Σ L_i = 13, W0 = doublingCount = 13 - m).  Since each chain's safe
set is a dilate (scale a_i) of the SAME RLL set, Bonferroni gives

    meas{t : ALL 13 runners clear 1/14} ≥ 1 - Σ_i (1 - mu_{L_i}),

and positivity ⟹ S is LONELY OUTRIGHT (no realization machinery: measure > 0
of the non-strict full safe set is a witness).  Output: mu_L table (L ≤ 13),
the full partition table of 13, and the closing frontier in W0.
"""
from fractions import Fraction as F
from functools import lru_cache

LO, HI = F(1, 14), F(13, 14)


def intersect(A, B):
    out, i, j = [], 0, 0
    while i < len(A) and j < len(B):
        a = max(A[i][0], B[j][0])
        b = min(A[i][1], B[j][1])
        if a < b:
            out.append((a, b))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def preimage_band(k):
    """T^{-k}(B) = union over i < 2^k of [(i+1/14)/2^k, (i+13/14)/2^k]."""
    d = 2 ** k
    return [((i + LO) / d, (i + HI) / d) for i in range(d)]


def chain_safe_measure(L):
    cur = [(F(0), F(1))]
    for k in range(L):
        cur = intersect(cur, preimage_band(k))
    meas = sum(b - a for a, b in cur)
    return meas, len(cur)


def partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield []
        return
    for first in range(min(n, mx), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest


def main():
    print("HYP-5853: chain-Bonferroni -- exact RLL measures of doubling chains")
    print("\n== per-chain safe measures mu_L (band [1/14, 13/14]) ==")
    mu = {}
    for L in range(1, 14):
        m, nint = chain_safe_measure(L)
        mu[L] = m
        union = 1 - F(L, 7)
        print(f"  L={L:2d}: mu_L = {m} = {float(m):.6f}   "
              f"[union bound {float(union):+.4f}]  "
              f"cost 1-mu = {float(1-m):.4f} vs L/7 = {float(F(L,7)):.4f}  "
              f"({nint} intervals)")

    print("\n== chain-Bonferroni partition table (partitions of 13) ==")
    print("   budget = Sum_i (1 - mu_{L_i});  PASS iff budget < 1 "
          "(=> lonely outright)")
    results = []
    for part in partitions(13):
        budget = sum(1 - mu[L] for L in part)
        W0 = 13 - len(part)
        results.append((W0, part, budget))
    results.sort(key=lambda r: (r[0], [-x for x in r[1]]))
    # summarize per W0: worst (max) budget over partitions with that W0
    from collections import defaultdict
    byW0 = defaultdict(list)
    for W0, part, budget in results:
        byW0[W0].append((budget, part))
    frontier = None
    for W0 in sorted(byW0):
        worst = max(byW0[W0], key=lambda x: x[0])
        best = min(byW0[W0], key=lambda x: x[0])
        all_pass = worst[0] < 1
        any_pass = best[0] < 1
        tag = "ALL PASS" if all_pass else ("some pass" if any_pass else "none")
        print(f"  W0={W0:2d} (m={13-W0:2d} chains): worst budget "
              f"{float(worst[0]):.4f} {worst[1]}; best {float(best[0]):.4f} "
              f"{best[1]}  -> {tag}")
        if all_pass and frontier is None:
            frontier = W0
    print(f"\n  CLOSING FRONTIER: every 13-set with doublingCount W0 >= "
          f"{frontier} is LONELY OUTRIGHT by chain-Bonferroni"
          if frontier is not None else
          "\n  no all-pass W0 level (chain gains alone insufficient)")

    # sanity: the 1/3 anchor is in every safe set
    print("\n== sanity: 1/3-anchor ==")
    for L in (4, 8, 13):
        cur = [(F(0), F(1))]
        for k in range(L):
            cur = intersect(cur, preimage_band(k))
        inside = any(a <= F(1, 3) <= b for a, b in cur)
        print(f"  L={L}: 1/3 in safe set: {inside}")


if __name__ == '__main__':
    main()
