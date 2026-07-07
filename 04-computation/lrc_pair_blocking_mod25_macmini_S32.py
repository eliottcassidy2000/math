#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S32 -- sharpening kps-S41's mod-25 covering core via the
PAIR-BLOCKING reformulation.

kps-S41 (HYP-4567, LRCMod25Floor.lean GREEN): M >= 2/25 is witnessed at t=c/25 if
some unit c has all v_i*c mod 25 in [2,23].  Open core: prove near-tight
no-mult-of-25 families ARE mod-25-clearable.

MY REFORMULATION (verified): at denominator 25, clearance-2 forbidden zone =
residues {0,1,24}.  For a UNIT c, v*c in {0,1,24} <=> v ≡ 0, +-c^{-1} mod 25.
Non-units mod 25 (0,5,10,15,20) are ALWAYS safe.  Units come in 10 +-pairs
{1,24},{2,23},...,{12,13}.  So:
    mod-25 covering FAILS (no clearing c)  <=>  V blocks ALL 10 unit +-pairs.
The AP {1..12} blocks all 10 (M=1/13, correctly not clearable).

THE SHARP CORE (this script): which 12-speed families (no mult of 25) block all
10 pairs, and are the pair-blockers with M < 2/25 EXACTLY the AP (+ dilates)?
If yes, then: non-blocker => M>=2/25 (kps floor); blocker => AP-like (this) =>
swept/1/13; mult-of-25 => loose.  That closes the lower bound = (G).
"""
import itertools
import random
from fractions import Fraction as F
from math import gcd

LO, HI = F(1, 13), F(2, 25)
UNIT_PAIRS = [(1, 24), (2, 23), (3, 22), (4, 21), (6, 19),
              (7, 18), (8, 17), (9, 16), (11, 14), (12, 13)]


def blocks_all_pairs(V):
    """True iff V's residues mod 25 hit every unit +-pair (mod-25 covering fails)."""
    res = {v % 25 for v in V}
    if 0 in res:
        return False  # a multiple of 25 -> handled separately (loose)
    for a, b in UNIT_PAIRS:
        if a not in res and b not in res:
            return False  # this pair is free -> c=inv clears -> M>=2/25
    return True


def _dens(W):
    d = set()
    for v, w in itertools.combinations(W, 2):
        d.add(v + w)
        if v != w:
            d.add(abs(v - w))
    for v in W:
        d.add(2 * v)
    d.discard(0)
    return d


def exact_M(W):
    best = F(0)
    seen = set()
    for s in _dens(W):
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best


def longest_dilated_AP(W):
    S = set(W)
    Ws = sorted(W)
    best = 1
    for i in range(len(Ws)):
        for j in range(i + 1, len(Ws)):
            a, e = Ws[i], Ws[j] - Ws[i]
            L = 2
            nxt = Ws[j] + e
            while nxt in S:
                L += 1
                nxt += e
            if L > best:
                best = L
    return best


def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(set(v // g for v in W))) if g > 1 else tuple(sorted(set(W)))


def gen(rng, hmax=50, n_random=120000):
    fams = set()
    # structured: dilated AP length L + defects (near-AP, likely pair-blockers)
    for e in range(1, 6):
        for a in range(1, e + 2):
            for L in range(8, 13):
                ap = [a + i * e for i in range(L)]
                if ap[-1] > hmax:
                    continue
                need = 12 - L
                pool = [x for x in range(1, hmax + 1) if x not in ap]
                if need == 0:
                    W = primitive(tuple(ap))
                    if len(W) == 12:
                        fams.add(W)
                    continue
                for _ in range(600):
                    W = primitive(tuple(sorted(set(ap) | set(rng.sample(pool, need)))))
                    if len(W) == 12:
                        fams.add(W)
    # random
    for _ in range(n_random):
        W = primitive(tuple(sorted(rng.sample(range(1, hmax + 1), 12))))
        if len(W) == 12:
            fams.add(W)
    return fams


def main():
    rng = random.Random(32)
    print("=" * 88)
    print("PAIR-BLOCKING mod 25: which families block all 10 unit +-pairs, and their M?")
    print("=" * 88)
    fams = gen(rng)
    print(f"  generated {len(fams)} 12-speed families (height<=50)")
    blockers = [W for W in fams if blocks_all_pairs(W)]
    print(f"  families BLOCKING all 10 pairs (mod-25 covering fails): {len(blockers)}")
    print(f"    (the rest -> some pair free -> M>=2/25 by kps's LRCMod25Floor)")
    print()
    # for pair-blockers, compute M and longest-AP
    below = []       # blockers with M < 2/25 (the residual to characterize)
    ap_like = 0
    minM_by_d = {}
    for W in blockers:
        M = exact_M(W)
        d = 12 - longest_dilated_AP(W)
        if d not in minM_by_d or M < minM_by_d[d][0]:
            minM_by_d[d] = (M, W)
        if M < HI:
            below.append((M, d, W))
    print(f"  pair-blockers with M < 2/25 (the sharp residual): {len(below)}")
    for M, d, W in sorted(below)[:15]:
        tag = "  <== M<1/13?!" if M < LO else ("  = AP/tight" if M == LO else "")
        print(f"    M={M} ({float(M):.5f}) defects={d}  W={list(W)}{tag}")
    print()
    print("  MIN M among pair-blockers, by defect count d=12-longestAP:")
    print(f"    {'d':>3} {'min M':>14} {'<2/25?':>7}  min-family")
    for d in sorted(minM_by_d):
        M, W = minM_by_d[d]
        Ws = str(list(W)) if len(str(list(W))) < 44 else str(list(W)[:10]) + "..]"
        print(f"    {d:>3} {str(M)+' ('+format(float(M),'.4f')+')':>20} "
              f"{'YES' if M < HI else 'no':>7}  {Ws}")
    print()
    print("  INTERPRETATION: if the ONLY pair-blockers with M<2/25 are the AP (d<=?)")
    print("  and near-AP families that are ALSO caught by the swept <=2-defect world,")
    print("  then kps's mod-25 core reduces to a clean finite statement: 'blocking all")
    print("  10 unit-pairs mod 25 + M<2/25 => AP-like'.  A d>=3 pair-blocker with")
    print("  M<2/25 would be the crux counterexample (report loudly).")


if __name__ == "__main__":
    main()
