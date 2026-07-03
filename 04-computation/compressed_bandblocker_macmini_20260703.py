#!/usr/bin/env python3
"""
Do COMPRESSED hge7 band-blockers need q -> infinity? (mac-mini-2026-07-03-S22, HYP-4040 refinement)
The hge7 OBLIGATION requires COMPRESSED families: no ratio-13 dominant, i.e. every runner i has some j with
|v_i| < 13|v_j|. The clean lcm family {1..11,13,lcm(2..X)} is DOMINANT (lcm is >13x all others) -> handled by
the dominant-runner route, NOT hge7. So HYP-4040's lower bound (q>X) is about an ALREADY-CLOSED case.
The REAL question: can a COMPRESSED family (>=7 far, no dominant) also be a band-blocker with q -> infinity?

Counting: to block band {15..Q} a family needs a runner divisible by each band prime; with 13 runners each
<= M, the product of covered primes <= M^13, so (primorial ~ e^Q) <= M^13 => Q <= 13 log M. So compressed
band-blockers can push q up to ~13 log(max-speed). We CONSTRUCT compressed geometric-chain band-blockers
(runners at ratios ~12 < 13, each highly composite) and check the witness denominator vs magnitude.
"""
from math import gcd, log
from functools import reduce
from sympy import primerange, nextprime
import random

def gcd_all(xs): return reduce(gcd, xs)
def danger(q): return {r for r in range(q) if min(r, q - r) * 14 < q}

def witness_denominator(speeds, qmax=2000):
    for q in range(15, qmax + 1):
        d = danger(q)
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in d for v in speeds):
                return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))
def far_count(speeds): return sum(1 for v in speeds if abs(v) > 22)
def is_compressed(speeds):
    return all(any(j != i and abs(speeds[i]) < 13 * abs(speeds[j]) for j in range(len(speeds)))
               for i in range(len(speeds)))

def build_compressed_blocker(scale, band_primes, rng):
    """13 runners in a geometric chain (ratio ~12), collectively divisible by band_primes, compressed+covering."""
    # anchors along a chain scale, 12*scale, 144*scale, ... within ratio 13
    anchors = [scale * (12 ** i) for i in range(4)]  # 4 scale rungs
    # distribute band primes across runners near anchors; each runner = anchor rounded to a multiple of its prime(s)
    speeds = []
    pi = 0
    for a in anchors:
        # put ~ (len(band_primes)/4) primes on a runner near this anchor as one composite <= 13*a
        prod = 1; used = []
        while pi < len(band_primes) and prod * band_primes[pi] * (a // max(prod,1) or 1) <= 13 * a:
            if prod * band_primes[pi] <= 13:  # keep building
                prod *= band_primes[pi]; used.append(band_primes[pi]); pi += 1
            else:
                if prod * band_primes[pi] * (a // (prod) if prod else a) <= 13 * a:
                    prod *= band_primes[pi]; used.append(band_primes[pi]); pi += 1
                else: break
        # runner = smallest multiple of prod in [a, 13a)
        if prod == 0: prod = 1
        m = ((a + prod - 1) // prod) * prod
        if m < a: m = a
        speeds.append(m)
    # remaining band primes: attach as multiples near the top anchor
    while pi < len(band_primes):
        p = band_primes[pi]; pi += 1
        top = anchors[-1]
        m = ((top + p - 1) // p) * p
        speeds.append(m)
    # pad with small covering runners to 13
    for q in [8,9,5,7,11,13,2,3]:
        if len(speeds) >= 13: break
        speeds.append(q)
    speeds = sorted(set(speeds))
    while len(speeds) < 13:
        speeds.append(rng.randint(2, 22))
        speeds = sorted(set(speeds))
    return speeds[:13]

if __name__ == "__main__":
    rng = random.Random(2026)
    print("COMPRESSED hge7 band-blockers: witness denominator vs magnitude scale.")
    print("=" * 84)
    print(f"{'scale':>8} {'#band primes':>13} {'max-speed':>12} {'compressed':>11} {'covering':>9} {'far':>4} {'q(witness)':>11}")
    for scale in [1, 5, 20, 60, 200, 600, 2000]:
        # band primes to try to block: grows with scale (more room)
        Q = 15 + int(10 * log(scale + 2))
        bps = list(primerange(15, Q))
        best = (0, None)
        for _try in range(4000):
            # randomize a bit: shuffle band primes and anchor multipliers
            rng.shuffle(bps)
            sp = build_compressed_blocker(scale, bps, rng)
            if len(sp) != 13 or gcd_all(sp) != 1: continue
            if far_count(sp) < 7 or not is_covering(sp) or not is_compressed(sp): continue
            q = witness_denominator(sp, qmax=1500)
            if q and q > best[0]: best = (q, sp)
        if best[1]:
            sp = best[1]
            print(f"{scale:>8} {len(bps):>13} {max(sp):>12} {str(is_compressed(sp)):>11} {str(is_covering(sp)):>9} {far_count(sp):>4} {best[0]:>11}")
        else:
            print(f"{scale:>8} {len(bps):>13} {'-':>12} {'(no valid compressed blocker found)':>40}")
    print("\n=> if q(witness) GROWS with scale for COMPRESSED families: the REAL hge7 obligation has")
    print("   unbounded-witness families -> arithmetic band CANNOT close hge7 alone (magnitude split forced).")
    print("   if q stays bounded: compressed forces near-equal-ish structure -> band {15..~50} may suffice.")
