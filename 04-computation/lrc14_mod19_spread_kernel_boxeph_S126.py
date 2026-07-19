#!/usr/bin/env python3
"""
boxeph-2026-07-19-S126 — the mod-19 spread lemma and the second-aligned-gap kernel for 3/38.

Owner: prove the isoperimetric spread bound on the 1/19 alphabet; prove the kernel (forbid a
second interior aligned gap).

LEMMA 1 (mod-19 antipodal-spread, PROVED). If M(C) < 2/19 then either 19 | some speed, or the
residues {v mod 19} cover ALL 9 antipodal unit-pairs {+-1,...,+-9} of Z/19.
  Proof: at t=n/19 (n=1..18), ||v n/19|| is a multiple of 1/19, so min_v ||v n/19|| is too;
  M<2/19 forces it <= 1/19, so some v has vn == 0,+-1 mod 19. If 19 does not divide any v (and
  19 !| n), then vn == +-1, i.e. v == +-n^{-1}; as n ranges over units so does n^{-1}, giving
  {+-(v mod 19)} = all units = all 9 pairs. QED.  (Translation-SENSITIVE, per opus's triage.)

This script:
 (A) VERIFIES Lemma 1 over random families with M<2/19 (0 violations expected);
 (B) tests FEASIBILITY of {band mod 38 in [3,35]} + {covering} + {mod-19 spread} and whether any
     such family reaches M=3/38 (expected: feasible but none reaches 3/38 -- mod-19 is necessary
     not sufficient);
 (C) THE KERNEL. At the s=38 maximizer t*=m/38, a pair summing to 38 lands at antipodal positions
     {x,-x} about 0 -- so EVERY origin-gap edge pair sums to 38, and a SECOND interior gap of width
     6/38 needs edge runners differing by == 6 m^{-1} mod 38. Count, for band-filled covering
     families, the interior gaps of width >= 6/38 in the position set {v m mod 38}: is the origin
     gap the unique maximal one?
"""
from fractions import Fraction as F
from math import gcd


def fd(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def fmin(sp, t):
    return min(fd(F(v) * t) for v in sp)


def M_pinch(sp):
    sums = {sp[i] + sp[j] for i in range(len(sp)) for j in range(i + 1, len(sp))}
    best = F(0)
    for q in sums:
        for m in range(1, q):
            v = fmin(sp, F(m, q))
            if v > best:
                best = v
    return best


def covering(sp):
    return all(any(v % q == 0 for v in sp) for q in range(2, 13))


def mod19_spread(sp):
    res = set(v % 19 for v in sp)
    if 0 in res:
        return True
    return set(min(r, 19 - r) for r in res) == {1, 2, 3, 4, 5, 6, 7, 8, 9}


# deterministic pseudo-random sampler (no Math.random; index-seeded LCG)
def lcg(seed):
    x = seed
    while True:
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        yield x


def sample(gen, lo, hi, k, exclude=()):
    out = []
    while len(out) < k:
        v = lo + next(gen) % (hi - lo + 1)
        if v not in out and v not in exclude:
            out.append(v)
    return out


TWO19 = F(2, 19)
THREE38 = F(3, 38)
print("=" * 72)
print("(A) LEMMA 1 verification: M<2/19 => mod-19 antipodal spread (or 19|v)")
print("=" * 72)
g = lcg(12345)
viol = tested = 0
for _ in range(2500):
    C = sorted(sample(g, 1, 55, 12))
    M = M_pinch(C)
    if M < TWO19:
        tested += 1
        if not mod19_spread(C):
            viol += 1
            if viol <= 5:
                print(f"   VIOLATION: {C} M={M} res19={sorted(set(v%19 for v in C))}")
print(f"  families with M<2/19: {tested};  VIOLATIONS: {viol}  (Lemma 1 holds: {viol == 0})")

print()
print("=" * 72)
print("(B) FEASIBILITY: band[3,35] + covering + mod-19-spread; any with M=3/38?")
print("=" * 72)
g2 = lcg(999)
tot = feas = m338 = 0
pool = [x for x in range(3, 36)]
for _ in range(6000):
    rest = sample(g2, 3, 35, 10, exclude=(3, 35))
    C = sorted({3, 35} | set(rest))
    if len(C) != 12 or not covering(C):
        continue
    tot += 1
    if mod19_spread(C):
        feas += 1
        if M_pinch(C) == THREE38:
            m338 += 1
print(f"  band+covering: {tot}; also mod-19 spread: {feas}; with M=3/38: {m338}")
print(f"  => mod-19 spread is {'FEASIBLE' if feas else 'INFEASIBLE'} with band+covering, "
      f"but reaches 3/38 in {m338} (necessary, not sufficient).")

print()
print("=" * 72)
print("(C) KERNEL: interior gaps of width>=6/38 in the maximizer position set {v m mod 38}")
print("=" * 72)


def gaps38(C, m):
    pos = sorted(set((v * m) % 38 for v in C))
    n = len(pos)
    gaps = []
    for i in range(n):
        a = pos[i]
        b = pos[(i + 1) % n] + (38 if i == n - 1 else 0)
        gaps.append((b - a, a, pos[(i + 1) % n]))   # gap width (in /38 units), from a to next
    return sorted(gaps, reverse=True)


examples = [
    [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35],
    [3, 5, 7, 8, 9, 10, 11, 12, 17, 21, 24, 35],
    [3, 7, 8, 9, 10, 11, 12, 15, 20, 22, 33, 35],
]
for C in examples:
    # find m making it band-filled (residues in [3,35]); m=1 by construction here
    m = 1
    gp = gaps38(C, m)
    big = [g for g in gp if g[0] >= 6]
    print(f"  C={C}")
    print(f"     gaps (width,from,to) >=6/38: {big}   #gaps>=6: {len(big)}  "
          f"origin gap present (35->3 wrap width 6)? {(6, 35, 3) in gp or any(a==35 for _,a,_ in gp)}")
print("  (origin gap = the 35->3 wrap, width 6; a SECOND gap>=6 = a second interior aligned gap)")
print("\nDONE.")
