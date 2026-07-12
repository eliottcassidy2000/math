# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont49: DC-completeness is combinatorially CHEAP -> the "<=6 coprime-to-30030 core"
# (opus-S243/S244, mac-mini cont.50) is NOT forced by divisor-completeness.
#
# CLAIM: {2,...,14} is coverable by ONE integer (360360 = lcm(2..14)) or TWO (2520, 143). So a DC family
# needs only 1 non-coprime-to-30030 runner; the other 12 can be coprime to 30030 = 2*3*5*7*11*13. Hence a
# spread primitive DC family can have up to 12 coprime-to-30030 runners -- the empirical "<=6 core" is a
# property of the fleet's bounded-diameter/adversarial SAMPLE, not a theorem. The correct universal invariant
# is: the coprime core is <=12 (>=1 runner covers DC), hence a <=12-speed family, LRC(<=13)-protected.
from math import gcd, lcm
from functools import reduce
from itertools import combinations
from fractions import Fraction as F

P30030 = 2*3*5*7*11*13   # 30030
def is_DC(v):  return all(any(x % d == 0 for x in v) for d in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1
def lrun(v):
    v = sorted(set(v)); b = m = 1
    for i in range(1, len(v)):
        if v[i] == v[i-1] + 1: m += 1; b = max(b, m)
        else: m = 1
    return b
def ncoprime(v): return sum(1 for x in v if gcd(x, P30030) == 1)
def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def Mlb(v, qcap=200):
    best = F(0)
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best = m
    return best
def min_clear_q(v, hi=200):
    for q in range(15, hi):
        if q % 14 and any(all(-(-q//14) <= (vi*p) % q <= q - -(-q//14) for vi in v) for p in range(1, q)):
            return q
    return None

def min_cover_of_reqs():
    # minimum number of integers whose divisor-sets cover {2,...,14}
    L = lcm(*range(2, 15))
    print(f"  lcm(2..14) = {L}; ONE runner = {L} covers all of {{2..14}}: {is_DC([L])}")
    # a nice 2-cover
    print(f"  two-cover {{2520,143}}: covers {{2..14}}? {is_DC([2520,143])}  (2520=2^3*3^2*5*7 covers 2..10,12,14; 143=11*13 covers 11,13)")

def main():
    print("(1) DC-completeness is CHEAP -- minimum integers covering the requirements {2,...,14}:")
    min_cover_of_reqs()

    print("\n(2) MAX coprime-to-30030 runners in a spread primitive DC family (constructive):")
    primes_gt13 = [17,19,23,29,31,37,41,43,47,53,59,61,67,71]
    # 1 DC-covering runner (360360) + 12 primes>13
    fam1 = [360360] + primes_gt13[:12]
    print(f"  A = {{360360}} + 12 primes>13: {sorted(fam1)}")
    print(f"     DC? {is_DC(fam1)}  primitive? {prim(fam1)}  longest-run={lrun(fam1)}  #coprime-to-30030 = {ncoprime(fam1)}")
    print(f"     M (lower bd, q<200) = {Mlb(fam1)} = {float(Mlb(fam1)):.4f}  (LRC needs >=1/14={1/14:.4f}); min-clear-q={min_clear_q(fam1)}")
    # 2-cover version (smaller diameter): {2520,143} + 11 primes>13
    fam2 = [2520, 143] + primes_gt13[:11]
    print(f"  B = {{2520,143}} + 11 primes>13: {sorted(fam2)}")
    print(f"     DC? {is_DC(fam2)}  primitive? {prim(fam2)}  longest-run={lrun(fam2)}  #coprime-to-30030 = {ncoprime(fam2)}")
    print(f"     M (lower bd) = {float(Mlb(fam2)):.4f}; min-clear-q={min_clear_q(fam2)}")

    print("\n(3) => the '<=6 coprime-to-30030 core' is NOT forced by DC (found families with 11-12 coprime).")
    print("    The UNIVERSAL invariant: coprime core <= 13 - (#runners needed to cover {2..14}) <= 12,")
    print("    a <=12-speed family, so LRC(<=13) gives it reach >= 1/13 > 1/14 REGARDLESS of size.")

    print("\n(4) is the '<=6 core' a BOUNDED-DIAMETER artifact? max #coprime-to-30030 over SMALL-diameter DC:")
    import random
    rng = random.Random(11)
    for vmax in (25, 40, 80, 200):
        mx = 0; ex = None
        for _ in range(30000):
            v = sorted(rng.sample(range(1, vmax+1), 13))
            if prim(v) and is_DC(v) and lrun(v) <= 7:
                c = ncoprime(v)
                if c > mx: mx, ex = c, v
        print(f"    Vmax<={vmax}: max #coprime-to-30030 = {mx}  (e.g. {ex})")
    print("    => if max grows with Vmax, the '<=6' is a small-diameter artifact; the primes-heavy large-diameter")
    print("       families (A,B above) reach 11-12 -- the core SIZE is unbounded, but always <=12 (LRC-protected).")

if __name__ == "__main__":
    main()
