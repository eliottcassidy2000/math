#!/usr/bin/env python3
"""
CRITICAL adversarial CONSTRUCTION (mac-mini-2026-07-03-S21, HYP-3877).
Can a covering gcd=1 hge7 family BLOCK the whole band {15..33} by DIVISION (every band q | some speed)?
If yes, the band {15..33} FAILS and Q* > 33 (or the route needs rethinking).
Strategy: pack lcm-multiples covering all 19 band moduli using <=13 speeds, ensure covering + gcd=1 + >=7 far.
Band primes 17,19,23,29,31 (lcm 17.19.23=7429<=8900, 29.31=899). Then find the smallest lonely q anyway.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce

def gcd_all(xs): return reduce(gcd, xs)

def danger_residues(q):
    return {r for r in range(q) if min(r, q - r) * 14 < q}

def lonely_at_q(speeds, q):
    dang = danger_residues(q)
    for a in range(1, q):
        if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
            return a
    return None

def smallest_lonely_q(speeds, qmax=2000):
    """smallest q (ANY q, not just band) with a lonely a/q -- uses correct danger test min-dist>=1/14."""
    for q in range(2, qmax + 1):
        dang = {r for r in range(q) if min(r, q - r) * 14 < q}
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
                return q, a
    return None, None

def smallest_band_q(speeds, qmax=2000):
    for q in range(15, qmax + 1):
        if lonely_at_q(speeds, q) is not None:
            return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def band_blocked_by_division(speeds):
    """which band moduli q in {15..33} divide some speed (=blocked by division / residue 0)."""
    return {q for q in range(15, 34) if any(v % q == 0 for v in speeds)}

if __name__ == "__main__":
    print("Constructing band-blockers: cover {15..33} by division with <=13 covering gcd=1 speeds.")
    print("=" * 88)
    # explicit packers: each speed divisible by a group whose lcm<=8900, together covering 15..33
    candidates = [
        # (label, speeds) -- try to divide all of 15..33
        ("primes+pp packing",
         [7429, 899, 800, 432, 462, 360, 26*7, 28*11, 33*13, 19*17, 25*23, 29*20, 31*18]),
        ("lcm-dense",
         [2*3*5*7*11, 17*19, 23*29, 31*16, 25*27, 2*3*5*7*13, 32*15, 18*22, 20*26, 21*24, 28*30, 33*23, 19*29]),
        ("maximal band cover",
         [17*19*23, 29*31, 25*32, 27*16, 21*22, 15*18, 20*26, 24*28, 30*33, 2*3*7*11, 5*13*2, 4*9*5, 8*3*7]),
    ]
    for label, speeds in candidates:
        speeds = sorted(set(v for v in speeds if v > 0))
        g = gcd_all(speeds)
        speeds1 = [v // g for v in speeds]  # gcd=1 reduce
        cov = is_covering(speeds1)
        blocked = band_blocked_by_division(speeds1)
        nfar = sum(1 for v in speeds1 if v > 22)
        bq = smallest_band_q(speeds1, qmax=400)
        sq, sa = smallest_lonely_q(speeds1, qmax=400)
        print(f"\n[{label}] {len(speeds1)} speeds, gcd(reduced)=1, covering={cov}, far={nfar}")
        print(f"   speeds (gcd-reduced) = {speeds1}")
        print(f"   band moduli blocked by division: {sorted(blocked)}  ({len(blocked)}/19)")
        print(f"   smallest BAND q (15..) lonely = {bq}")
        print(f"   smallest ANY q lonely = {sq} (a={sa})   [must exist -- LRC14]")

    print("\n" + "=" * 88)
    print("TARGETED: maximize band moduli blocked by division among COVERING gcd=1 hge7 families.")
    print("=" * 88)
    import random
    rng = random.Random(999)
    best_block = (0, None)
    for _ in range(400000):
        # bias speeds toward band multiples
        pool = []
        for q in [17,19,23,29,31,25,27,32,16,21,22,33,15,18,20,24,26,28,30]:
            pool.append(q * rng.randint(1, 8900 // q))
        speeds = rng.sample(pool, min(13, len(pool)))
        while len(speeds) < 13:
            speeds.append(rng.randint(1, 22))
        speeds = sorted(set(speeds))
        if len(speeds) != 13: continue
        if gcd_all(speeds) != 1: continue
        if sum(1 for v in speeds if v > 22) < 7: continue
        if not is_covering(speeds): continue
        blocked = band_blocked_by_division(speeds)
        if len(blocked) > best_block[0]:
            best_block = (len(blocked), speeds)
    nb, sp = best_block
    print(f"MAX band moduli blocked by division (over covering gcd=1 hge7) = {nb}/19")
    if sp:
        print(f"   family: {sp}")
        print(f"   blocked: {sorted(band_blocked_by_division(sp))}")
        print(f"   smallest band-q lonely = {smallest_band_q(sp, qmax=600)}")
    print("\n=> if MAX blocked < 19, the band {15..33} is NEVER fully division-blocked -> a free q always exists")
    print("   (and the +-1 refinement only needs a generic residue among the ~unblocked moduli).")
