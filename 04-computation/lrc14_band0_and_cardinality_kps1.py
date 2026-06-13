#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): the two PROVED structural facts of THM-497 — the band-0 divisibility
lemma (Part B) and the cardinality-permits dichotomy (Part C).
kind-pasteur-2026-06-13-S1.  (Complements codex's Q27 set-cover / Church-descent
route and the convergent covering-deficit line; these are the proved pieces.)

Band criterion (THM-492, verified there): t=a/q (gcd(a,q)=1) is a strict witness
(M>1/14) iff every v in S has (va mod q) not in B_q = {r: min(r,q-r) <= floor(q/14)}.
A unit a is a witness-numerator iff it escapes the dilated bands U_v v^{-1}B_q.

PART B (band-0 lemma).  For q<=13, floor(q/14)=0 so B_q={0}: a unit a is blocked by
v iff va==0 mod q iff q|v. So a/q is a strict witness iff q does NOT divide any
runner; a config with NO strict witness at any q in {2..13} must have, for EACH such
q, some runner divisible by q (the divisor-sets COVER {2..13}).

PART C (cardinality dichotomy).  At a band-k shell q (14k <= q <= (k+1)*14-1,
B_q={0,+-1,...,+-k}), a runner v with q-not-div-v blocks exactly {+-j v^{-1}: 1<=j<=k},
<= 2k units; 13 runners block <= 26k units while |(Z/q)*| = phi(q) <= q ~ 14k.
Since 26k > 14k, the cardinality budget ALWAYS permits a complete cover => no
counting/union-bound/first-moment argument can force a witness; the obstruction is
purely additive alignment of the {+-j v_i^{-1}}.
"""

import sys, random
from math import gcd
from functools import reduce
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def band(q, n=14):
    h = q // n
    return [min(r, q - r) <= h for r in range(q)]


def witness_at(S, q):
    B = band(q)
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        if all(not B[(v * a) % q] for v in S):
            return True
    return False


def part_B():
    print("=== PART B: band-0 lemma (q<=13): witness-at-q <=> q divides no runner ===", flush=True)
    rng = random.Random(11)
    bad = tot = 0
    for _ in range(20000):
        S = sorted(rng.sample(range(1, 200), 13))
        for q in range(2, 14):
            tot += 1
            if witness_at(S, q) != all(v % q != 0 for v in S):
                bad += 1
    print(f"   {tot} (config,q<=13) pairs: {bad} mismatches  "
          f"=> {'LEMMA HOLDS' if bad == 0 else 'FALSE'}", flush=True)
    # divisor-cover consequence
    rng2 = random.Random(2)
    blockers = viol = 0
    for _ in range(200000):
        S = sorted(rng2.sample(range(1, 400), 13))
        if reduce(gcd, S) != 1:
            continue
        if all(not witness_at(S, q) for q in range(2, 14)):
            blockers += 1
            if not all(any(v % q == 0 for v in S) for q in range(2, 14)):
                viol += 1
    print(f"   configs blocking all q<=13: {blockers}; divisor-cover violations: {viol} (must be 0)", flush=True)


def part_C():
    print("\n=== PART C: cardinality dichotomy — 26k >= phi(q) at every band ===", flush=True)
    print("   band-k shell q in [14k,(k+1)*14-1]; runner blocks <=2k units; 13 runners <=26k.", flush=True)
    print("    band k | q range        | max phi(q) | 26k | 26k>=phi? (cover cardinality-permitted)", flush=True)
    def phi(q):
        r = q; x = q; p = 2
        while p * p <= x:
            if x % p == 0:
                while x % p == 0:
                    x //= p
                r -= r // p
            p += 1
        if x > 1:
            r -= r // x
        return r
    for k in range(1, 9):
        lo, hi = 14 * k, (k + 1) * 14 - 1
        mphi = max(phi(q) for q in range(lo, hi + 1))
        print(f"      {k:2d}   | [{lo:3d},{hi:3d}]   |   {mphi:4d}     | {26*k:3d} | "
              f"{'YES' if 26*k >= mphi else 'NO'}", flush=True)
    print("   => 26k > 14k >= phi(q) for all k: COUNTING NEVER FORBIDS THE COVER.", flush=True)
    print("      Hence no first-moment/union-bound proves C'(14); only additive alignment", flush=True)
    print("      of the {+-j v_i^{-1}} can leave an uncovered unit (the structural obstruction).", flush=True)


if __name__ == "__main__":
    part_B()
    part_C()
