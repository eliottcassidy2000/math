#!/usr/bin/env python3
"""
THE CONSECUTIVE CLOSED FORM, THE c = 8 THEOREM, THE EDGE CAP, AND THE
TREE-HUNTER CEILING (boxeph-2026-07-17-S71)

(A) THE CONSECUTIVE CLOSED FORM:  for r = k mod 7,
        mu(D_k cap D_{k+1}) = 1/49 + r(6-r) / (49 k (k+1)).
    Derivation: the LEM-042 integer sum for (k, k+1) gives
    49*Sigma - 14k(k+1) = 14 r(6-r) exactly (the quasi-polynomial
    collapses).  Zero excess exactly at k == 0 (mod 7); max at r = 3.
    Specializes to mu = 1/(7(k+1)) for k <= 6 and 11/504 at k = 8.

(B) THE c = 8 CONSECUTIVE THEOREM: for every consecutive block
    [v, v+1, ..., v+7], the seven path credits sum to
        1/7 + sum_{i: r_i != 0} r_i(6-r_i)/(49(v+i)(v+i+1)) ,
    and a window of seven consecutive k's contains EXACTLY ONE r = 0, so
    the excess is STRICTLY positive:  good >= excess > 0.  EVERY
    consecutive 8-block crosses the wall (uniformly in v).
    The c = 9 consecutive boundary: needs excess > 6/49 -- fails for all
    v >= 2 (measured; the excesses decay like 1/v^2).

(C) THE EDGE CAP: mu(D_a cap D_b) <= 1/14 for ALL distinct pairs, with
    equality iff the reduced pair is (1,2).  Proof architecture: each term
    of the LEM-042 sum is <= the cap 2a, and there are 2n+1 <= (a+b)/7 + 1
    terms, giving mu <= (a+b+7)/(49b) <= 1/14 whenever 35b >= 14a + 98,
    i.e. outside a FINITE set (a <= 4, b <= 4), checked directly.

(D) THE TREE-HUNTER (Hunter 1976, arbitrary spanning tree; klein's Lean
    lemma is the path case; the tree case = leaf-plucking induction,
    named):  mu(union A_i) <= sum mu(A_i) - sum_{(i,j) in T} mu(A_i cap A_j).
    Referee: 60 random families x random spanning trees, exact.
    THE CEILING: credits <= (c-1)/14 by (C), and crossing needs
    > (c-7)/7 = (2c-14)/14, so the tree-hunter can cross ONLY for c <= 12
    (at c = 13: 12/14 = 12/14, not strict -- IMPOSSIBLE for every tree).
    Feasibility demos: doubling families cross at c = 9..12 (credits
    (c-1)/14 vs needed; brute-exact good verified).

All exact (Fractions / integers).
"""

import sys
import random
from fractions import Fraction as Fr
from math import gcd

sys.path.insert(0, '04-computation')
from lrc14_pair_overlap_law_boxeph_S69 import mu_brute
from lrc14_cone_floor_wall_boxeph_S70 import mu_fr


def consecutive_closed(k):
    r = k % 7
    return Fr(1, 49) + Fr(r * (6 - r), 49 * k * (k + 1))


if __name__ == "__main__":
    print("CONSECUTIVE FORM / c=8 THEOREM / EDGE CAP / TREE-HUNTER (boxeph S71)")
    print("=" * 78)
    print("PART A -- the consecutive closed form")
    for k in range(1, 401):
        assert mu_fr(k, k + 1) == consecutive_closed(k), k
    print("  mu(k, k+1) == 1/49 + r(6-r)/(49k(k+1)) exact for k = 1..400")

    print()
    print("PART B -- the c = 8 consecutive theorem; the c = 9 boundary")
    worst_excess = None
    for v in range(1, 201):
        ex = sum(consecutive_closed(v + i) - Fr(1, 49) for i in range(7))
        zeros = sum(1 for i in range(7) if (v + i) % 7 == 0)
        assert zeros == 1 and ex > 0, v
        if worst_excess is None or ex < worst_excess[0]:
            worst_excess = (ex, v)
    print(f"  excess > 0 for all v <= 200 (exactly one r = 0 per window); "
          f"min excess {float(worst_excess[0]):.2e} at v = {worst_excess[1]}")
    for v in (1, 13, 50):
        blk = [v + i for i in range(8)]
        credits = sum(mu_fr(blk[i], blk[i + 1]) for i in range(7))
        good = 1 - mu_brute(blk, want_all=False)
        ledger = credits - Fr(1, 7)
        print(f"  block [{v}..{v+7}]: ledger = {float(ledger):.6f} > 0: "
              f"{ledger > 0}; true good = {float(good):.5f} >= ledger: "
              f"{good >= ledger}")
        assert 0 < ledger <= good
    c9 = [v for v in range(1, 100)
          if sum(consecutive_closed(v + i) - Fr(1, 49) for i in range(8))
          > Fr(6, 49)]
    print(f"  c = 9 consecutive route: crossing v's = {c9 if c9 else 'NONE'}"
          f" (the pure consecutive-pair route dies at c = 9 for v >= 2)")

    print()
    print("PART C -- the edge cap mu <= 1/14, equality iff reduced (1,2)")
    mx = None
    eq = []
    for a in range(1, 301):
        for b in range(a + 1, 302):
            if gcd(a, b) != 1:
                continue
            m = mu_fr(a, b)
            assert m <= Fr(1, 14), (a, b)
            if m == Fr(1, 14):
                eq.append((a, b))
            if mx is None or m > mx[0]:
                mx = (m, (a, b))
    print(f"  max over reduced pairs a < b <= 300: {mx[0]} at {mx[1]}; "
          f"equality cases: {eq}")
    assert eq == [(1, 2)]

    print()
    print("PART D -- tree-hunter validity + the ceiling + doubling demos")
    rng = random.Random(14)
    for trial in range(60):
        c = rng.randint(3, 7)
        fam = sorted(rng.sample(range(2, 60), c))
        # random spanning tree on 0..c-1
        edges = []
        nodes = list(range(c))
        rng.shuffle(nodes)
        for i in range(1, c):
            edges.append((nodes[i], nodes[rng.randrange(i)]))
        union = mu_brute(fam, want_all=False)
        rhs = Fr(c, 7) - sum(mu_brute([fam[i], fam[j]]) for i, j in edges)
        assert union <= rhs, (fam, edges)
    print("  tree-hunter inequality exact on 60 random family/tree instances")
    for c in (9, 10, 12):
        fam = [3 * 2 ** i for i in range(c)]
        credits = sum(mu_fr(fam[i], fam[i + 1]) for i in range(c - 1))
        need = Fr(c - 7, 7)
        good = 1 - mu_brute(fam, want_all=False)
        ledger = credits - need
        print(f"  doubling c={c} {fam[:4]}...: credits {credits} "
              f"(= (c-1)/14: {credits == Fr(c - 1, 14)}) vs needed {need}; "
              f"ledger = {float(ledger):.5f} > 0: {ledger > 0}; true good = "
              f"{float(good):.5f} >= ledger: {good >= ledger}")
        assert credits == Fr(c - 1, 14) and ledger > 0 and good >= ledger
    print(f"  THE CEILING: at c = 13 max credits 12/14 vs needed 12/14 -- "
          f"NOT strict: the tree-hunter cannot cross c = 13 (any tree, any "
          f"family): {Fr(12, 14) > Fr(12, 14)} strict fails, PROVED by (C)")
    print("=" * 78)
    print("done")
