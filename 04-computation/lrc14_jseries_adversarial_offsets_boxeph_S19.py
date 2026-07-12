#!/usr/bin/env python3
"""
lrc14_jseries_adversarial_offsets_boxeph_S19.py  (boxeph-2026-07-12-S19)

death-star-S14's named "cheap next probe": the j-series with ADVERSARIAL offsets.
Family class (compressed stratum, j off-lattice runners):

    v(L, j, r) = L*{1..13-j}  u  {(13-j+i)*L + r_i : i = 1..j},   1 <= r_i <= L-1.

death-star's uniform +1 offsets gave EXACT M = 1/(14-j) (j = 1..7). Question: can
adversarial offsets push the stratum floor BELOW 1/13 (toward 1/14), or does the
floor stay at ~1/13 for L >= 2?

THE FRAME (new observation): at L = 1 the offset constraint r_i < L is vacuous and
the class boundary contains the KNOWN tight points -- GW {1..11, 13, 24} IS the
(L=1, j=2) member with r = (1, 11), M = 1/14 exactly (THM-708). THM-709 (the
doubling-tight locus is a singleton) suggests no L >= 2 analog exists. So the probe
tests: [L >= 2 compressed strata floor = 1/13 (conjecture)] vs [tight points leak
into L >= 2]. Either outcome shapes death-star's "extend the proved stratum to
j <= 12" plan.

Method: coordinate descent on r (exact capped evaluator; loose candidates abort on
the small-q grid). DC verified programmatically; primitivity forced by odd/coprime
offsets. L chosen so the LATTICE part carries all divisor duties:
  j <= 6: L = 858 = 2*3*11*13  (kept {1..7}L covers d=2..14: 8|4L, 9|3L, 5|5L, 7|7L...)
  j = 7:  L = 6006 = 858*7     (kept {1..6}L needs 7 | L)
"""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce

def is_dc(v):
    return all(any(x % d == 0 for x in v) for d in range(2, 15))

def M_capped(v, capn, capd):
    """Exact M if M < cap else None. Small q first for fast aborts."""
    qs = list(range(2, 61))
    big = set()
    n = len(v)
    for i in range(n):
        for j in range(i, n):
            big.add(v[i] + v[j])
            if v[i] != v[j]: big.add(abs(v[i] - v[j]))
    qs += sorted(q for q in big if q > 60)
    mb, qb = 0, 1
    for q in qs:
        for p in range(1, q // 2 + 1):
            m = q
            for x in v:
                r = (x * p) % q
                if r > q - r: r = q - r
                if r < m:
                    m = r
                    if m * qb <= mb * q: break
            if m * qb > mb * q:
                mb, qb = m, q
                if mb * capd >= capn * qb:
                    return None
    return (mb, qb)

def family(L, j, rs):
    v = [L * k for k in range(1, 14 - j)] + [(13 - j + i + 1) * L + rs[i] for i in range(j)]
    # note: i runs 0..j-1 => multipliers (14-j)..13
    return sorted(v)

def probe(L, j, iters=140, seed=7):
    random.seed(seed * 1000 + j)
    # start at death-star's +1 point
    cur = [1] * j
    v0 = family(L, j, cur)
    assert len(set(v0)) == 13 and is_dc(v0) and reduce(gcd, v0) == 1, "bad base family"
    best = M_capped(v0, 1, 1)  # cap 1 => always exact
    bestr = list(cur)
    moves = [1, -1, 2, -2, 3, -3, 5, -5, L // 13, -(L // 13), L // 4, -(L // 4)]
    for it in range(iters):
        cand = list(bestr if random.random() < 0.7 else cur)
        i = random.randrange(j)
        cand[i] = cand[i] + random.choice(moves)
        if not (1 <= cand[i] <= L - 1): continue
        v = family(L, j, cand)
        if len(set(v)) != 13: continue
        if reduce(gcd, v) != 1 or not is_dc(v): continue
        r = M_capped(v, best[0], best[1])
        if r is not None:
            best = r; bestr = list(cand)
        cur = cand
    return best, bestr

def main():
    fast = '--fast' in sys.argv
    print("j-series adversarial offsets probe (death-star-S14 handoff)")
    print(f"{'j':>2s} {'L':>5s} {'+1 baseline 1/(14-j)':>21s} {'adversarial min M':>20s} {'offsets':>28s}")
    for j in range(1, 8):
        L = 858 if j <= 6 else 6006
        iters = 60 if (fast or j == 7) else 140
        (m, q), rs = probe(L, j, iters=iters)
        base = F(1, 14 - j)
        print(f"{j:2d} {L:5d} {str(base):>21s} {m}/{q} = {m/q:20.6f} {str(rs):>28s}")
    print("\n1/13 = 0.076923, 1/14 = 0.071429.  L>=2 floor staying at ~1/13 supports extending the")
    print("proved compressed stratum (death-star); a dip below 1/13 would locate an L>=2 GW-analog")
    print("(contradicting the THM-709 singleton reading) and mark the stratum genuinely tight-adjacent.")

if __name__ == '__main__':
    main()
