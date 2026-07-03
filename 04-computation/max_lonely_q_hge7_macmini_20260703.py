#!/usr/bin/env python3
"""
Is the smallest lonely q UNIFORMLY BOUNDED over covering gcd=1 hge7 families? (mac-mini-S21, HYP-3877)
Band-blockers (speeds divisible by all small moduli) push the smallest lonely q above 33. To need a large
lonely q, a family must block {2..q-1} by division/+-1 -> needs speeds covering primes up to q-1 -> with 13
speeds and bounded magnitude, how large can q get? If q grows with the speed bound, NO uniform band.

We (a) greedily CONSTRUCT band-blockers up to a magnitude bound, (b) random-search maximize smallest-lonely-q,
across speed magnitude bounds 10^3..10^7, to see if max-q plateaus (bounded route) or grows (unbounded).
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)

def smallest_lonely_q(speeds, qmax=4000):
    for q in range(2, qmax + 1):
        dang = {r for r in range(q) if min(r, q - r) * 14 < q}
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
                return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def far_count(speeds): return sum(1 for v in speeds if abs(v) > 22)

def greedy_blocker(prime_bound, mag_bound):
    """13 speeds each = product of primes to cover as many small moduli as possible under mag_bound."""
    primes = [p for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71] if p <= prime_bound]
    speeds = []
    # assign each speed a bundle of primes whose product <= mag_bound, covering distinct primes
    i = 0
    while i < len(primes) and len(speeds) < 13:
        v = 1; used = []
        while i < len(primes) and v * primes[i] <= mag_bound:
            v *= primes[i]; used.append(primes[i]); i += 1
        if v > 1: speeds.append(v)
        else: i += 1
    # pad with small covering speeds
    for q in range(2, 15):
        if len(speeds) >= 13: break
        if not any(s % q == 0 for s in speeds): speeds.append(q)
    while len(speeds) < 13: speeds.append(random.randint(2, 22))
    speeds = sorted(set(speeds))[:13]
    return speeds

if __name__ == "__main__":
    print("Max smallest-lonely-q over covering gcd=1 hge7 families, vs speed magnitude bound.")
    print("=" * 88)
    print(f"{'mag_bound':>10} {'method':>14} {'#families':>10} {'MAX smallest-lonely-q':>22} {'example q':>10}")
    for mag in [10**3, 10**4, 10**5, 10**6, 10**7]:
        rng = random.Random(mag)
        best = (0, None)
        # (a) greedy blockers with various prime bounds
        for pb in range(15, 75):
            sp = greedy_blocker(pb, mag)
            if len(sp) == 13 and gcd_all(sp) == 1 and far_count(sp) >= 7 and is_covering(sp):
                q = smallest_lonely_q(sp, qmax=4000)
                if q and q > best[0]: best = (q, sp)
        # (b) random band-multiple families
        band = list(range(15, min(mag, 120)))
        for _ in range(120000):
            pool = [q * rng.randint(1, max(1, mag // q)) for q in band]
            speeds = sorted(set(rng.sample(pool, min(11, len(pool))) + rng.sample(range(2,23), 2)))
            if len(speeds) < 13:
                speeds = sorted(set(speeds) | {rng.randint(2, mag) for _ in range(13 - len(speeds))})
            speeds = sorted(set(speeds))[:13]
            if len(speeds) != 13 or gcd_all(speeds) != 1 or far_count(speeds) < 7 or not is_covering(speeds):
                continue
            q = smallest_lonely_q(speeds, qmax=1500)
            if q and q > best[0]: best = (q, speeds)
        print(f"{mag:>10} {'greedy+random':>14} {'~120k':>10} {best[0]:>22} {'':>10}")
        if best[1]:
            print(f"           example: {best[1]}  -> smallest lonely q = {best[0]}")
    print("\n=> if MAX-q PLATEAUS as mag grows: uniform band exists (route closes hge7 with a bigger band).")
    print("   if MAX-q GROWS with mag: no uniform band -> band route closes near-equal/generic only, not all hge7.")
