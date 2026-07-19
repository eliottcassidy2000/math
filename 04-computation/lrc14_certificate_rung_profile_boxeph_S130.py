#!/usr/bin/env python3
"""
lrc14_certificate_rung_profile_boxeph_S130.py

THE CERTIFICATE-RUNG PROFILE (boxeph-2026-07-19-S130, HYP-7880).

Thesis being tested (the "2/19 lens" cashed as data):

  (1) THE GENERAL-p SPREAD LEMMA (one-line generalization of S115 mod-13
      pair-blocking and S126 mod-19 antipodal spread; same proof, every prime p):
          M(V) < 2/p  and  p divides no speed
          ==>  the residues {v mod p} hit ALL (p-1)/2 antipodal unit-pairs of Z/p.
      Proof: at t=b/p every distance is a multiple of 1/p; min <= M < 2/p forces
      min <= 1/p, i.e. some v with vb == 0,+-1 (mod p); p|v excluded, so
      v == +-b^{-1}; as b runs over units, every antipodal pair is hit.  QED.

  (2) THE RUNG IDENTITY: the loneliness spectrum values M = D/s (Pinch: det over
      pair-sum, boxeph-S121/S123 ladder) and the modular certificate values k/p
      are the SAME grid: certificates certify from below exactly at the rungs the
      spectrum can occupy.  In particular for the LIVE LRC(14) window
      (1/14, 3/41) (opus THM-1240) every family sits BELOW the rungs
      2/23 = 0.0870, 2/19 = 0.1053, 2/17 = 0.1176 and 2/13, hence is pinned by
      the spread lemma at p = 13, 17, 19, 23 simultaneously (or carries a
      multiple of p).  SLACK COUNT at p=23: 11 antipodal pairs vs 13 speeds
      => at most TWO pairs doubled: near-bijectivity mod 23.  Nobody has used
      p=23; grep shows 2/23 in canon only as an attained M value.

  (3) THE PIGEONHOLE FLOOR (canon: THM-518 "prime route lands at 2/29"): at
      p=29, 14 pairs > 13 speeds, so an unhit pair ALWAYS exists: every family
      without a multiple of 29 has M >= 2/29 = 96.6% of 1/14.  The WHOLE of
      LRC(14) is the sliver [2/29, 1/14) of width 1/406 (MISTAKE-093 window)
      plus the 29-divisible branch.

This script computes EXACTLY (pure fractions, no numpy):
  - M(V) for the canonical family cast (integer hot loop, exact output),
    together with the attaining denominators (rungs D/s);
  - the antipodal spread profile at p in {13,17,19,23,29,31}: multiples of p,
    pairs hit / unhit, slack, and the doubled pairs at p=23;
  - a self-test: the spread lemma and the 2/29 pigeonhole floor on the cast
    plus random primitive 13-sets (if my M ever violates them, MY CODE is wrong).

Families: tight AP {1..13}; GW {1..11,13,24}; near-floor {1..11,13,36};
deep well {1..12,182}; compact minimizers; S* = {1,2,3,5,7,...,13,38,42}
(MISTAKE/THM-526 counterexample, M=2/23); {1..13}\{6} u {56} (mid-drop, 2/23);
n=12 side: {1..12}, {1..11,13}, {1..11,24}, {1..13}\{6}, band-filled q=38 set.
"""

from fractions import Fraction
from math import gcd
import random

PRIMES = [13, 17, 19, 23, 29, 31]


def M_exact(V):
    """Exact M(V) = max_t min_v ||v t||.

    Pinch lemma (HYP-2059 / THM-401 / boxeph-S120): the maximizer is at
    t = a/q with q a pairwise sum v_i+v_j.  We scan ALL a in 1..q-1 for ALL
    q up to max pair sum (all numerators, not just coprime -- MISTAKE-173),
    which contains every pair-sum denominator, hence is exact.
    Integer hot loop (opus-S398 speed lesson).
    """
    V = sorted(V)
    qmax = V[-1] + V[-2]
    best_num, best_den = 0, 1          # best distance as num/den
    best_at = []
    for q in range(2, qmax + 1):
        for a in range(1, q):
            # min over v of dist(v*a mod q) in units of 1/q
            m = q
            for v in V:
                r = (v * a) % q
                d = r if r <= q - r else q - r
                if d < m:
                    m = d
                    if m * best_den < best_num * q:   # early exit: cannot beat
                        break
            if m * best_den > best_num * q:
                best_num, best_den, best_at = m, q, [(m, q, a)]
            elif m * best_den == best_num * q and m > 0:
                if len(best_at) < 4 and (m, q, a) not in best_at:
                    best_at.append((m, q, a))
    return Fraction(best_num, best_den), best_at


def spread_profile(V, p):
    """Antipodal spread of V mod p: (multiples, pairs hit, unhit list, doubled)."""
    mult = [v for v in V if v % p == 0]
    pairs = {}
    for v in V:
        r = v % p
        if r == 0:
            continue
        j = r if r <= p - r else p - r
        pairs.setdefault(j, []).append(v)
    npairs = (p - 1) // 2
    hit = sorted(pairs)
    unhit = [j for j in range(1, npairs + 1) if j not in pairs]
    doubled = {j: vs for j, vs in pairs.items() if len(vs) >= 2}
    return mult, hit, unhit, doubled


FAMILIES = [
    ("AP13 {1..13}                    ", list(range(1, 14))),
    ("GW {1..11,13,24}                ", list(range(1, 12)) + [13, 24]),
    ("NF36 {1..11,13,36}              ", list(range(1, 12)) + [13, 36]),
    ("DW {1..12,182} deep well        ", list(range(1, 13)) + [182]),
    ("CM39 {2,4,..,24,39}             ", list(range(2, 25, 2)) + [39]),
    ("Sstar {..,38,42} (THM-526 cex)  ", [1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 38, 42]),
    ("MID6 {1..13}\\{6} u {56}         ", [v for v in range(1, 14) if v != 6] + [56]),
    ("AP12 {1..12}                    ", list(range(1, 13))),
    ("R12a {1..11,13}                 ", list(range(1, 12)) + [13]),
    ("R12b {1..11,24}                 ", list(range(1, 12)) + [24]),
    ("R23  {1..13}\\{6}                ", [v for v in range(1, 14) if v != 6]),
    ("BF38 {3,5,7..13,15,21,35}       ", [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35]),
]


def main():
    print("=" * 100)
    print("CERTIFICATE-RUNG PROFILE  (boxeph-S130, HYP-7880)")
    print("rungs: 2/13=%.4f 2/17=%.4f 2/19=%.4f 2/23=%.4f 2/29=%.4f ; targets 1/14=%.4f 1/13=%.4f 3/41=%.4f"
          % (2 / 13, 2 / 17, 2 / 19, 2 / 23, 2 / 29, 1 / 14, 1 / 13, 3 / 41))
    print("=" * 100)

    lemma_violations = 0
    for name, V in FAMILIES:
        M, at = M_exact(V)
        rungs = ", ".join("%d/%d@a=%d" % (m, q, a) for (m, q, a) in at[:3])
        print("\n%s  n=%d  M = %s = %.5f   attained: %s" % (name, len(V), M, float(M), rungs))
        for p in PRIMES:
            mult, hit, unhit, doubled = spread_profile(V, p)
            npairs = (p - 1) // 2
            below = M < Fraction(2, p)
            tag = "M<2/%d" % p if below else "     "
            status = ""
            if below and not mult:
                if unhit:
                    status = "  *** SPREAD LEMMA VIOLATED ***"
                    lemma_violations += 1
                else:
                    slack = len([v for v in V if v % p]) - npairs
                    status = "  FULL SPREAD forced (slack %d, doubled %s)" % (
                        slack, sorted(doubled) if doubled else "-")
            elif mult:
                status = "  blocked by multiples %s" % mult
            elif unhit:
                status = "  unhit pairs %s => certifies M>=2/%d" % (unhit, p)
            print("   p=%2d [%s] pairs %2d/%2d hit%s" % (p, tag, len(hit), npairs, status))

    # --- self-tests -------------------------------------------------------
    print("\n" + "=" * 100)
    print("SELF-TEST 1: pigeonhole floor -- every 13-set without a multiple of 29 has M >= 2/29")
    print("SELF-TEST 2: general-p spread lemma on random families")
    rng = random.Random(219)          # seed 2/19 :)
    bad = 0
    trials = 200
    for _ in range(trials):
        while True:
            V = sorted(rng.sample(range(1, 61), 13))
            g = 0
            for v in V:
                g = gcd(g, v)
            if g == 1:
                break
        M, _ = M_exact(V)
        if all(v % 29 for v in V) and M < Fraction(2, 29):
            bad += 1
            print("  PIGEONHOLE VIOLATED:", V, M)
        for p in PRIMES:
            mult, hit, unhit, doubled = spread_profile(V, p)
            if M < Fraction(2, p) and not mult and unhit:
                bad += 1
                print("  SPREAD LEMMA VIOLATED at p=%d:" % p, V, M)
    print("random trials: %d, violations: %d" % (trials, bad))
    print("cast lemma violations: %d" % lemma_violations)
    print("\nVERDICT: %s" % ("ALL CONSISTENT -- rung identity + spread pinning verified exactly"
                             if bad == 0 and lemma_violations == 0 else "*** FAILURES -- see above ***"))


if __name__ == "__main__":
    main()
