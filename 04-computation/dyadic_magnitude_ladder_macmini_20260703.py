#!/usr/bin/env python3
"""
THE DYADIC MAGNITUDE LADDER (mac-mini-2026-07-03-S22, HYP-4040).
Analogy to arXiv:2607.00876 (Bairaktari-Larsen, binary-tree continual counting, Omega(log^{3/2} n)).
That paper: controlling prefix sums across ALL time scales costs ~log (via a dyadic/binary-tree decomposition),
and the discrepancy lower bound is tight. LRC analogue: a covering family must dodge the danger residues at
ALL small moduli (scales) at once; the band-blockers (lcm families) are engineered to block every small scale,
pushing the smallest lonely denominator up. QUESTION: does the worst smallest-lonely-q Q(k) over covering
gcd=1 hge7 families with max-speed in [2^k, 2^{k+1}) grow LINEARLY in the dyadic level k = log2(max-speed)?
(linear in k  <=>  q ~ log(max-speed) -- the LRC "all-scales" cost, the discrepancy exponent.)
"""
from math import gcd, log2
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)

def smallest_lonely_q(speeds, qmax=4000):
    for q in range(15, qmax + 1):
        dang = {r for r in range(q) if min(r, q - r) * 14 < q}
        for a in range(1, q):
            if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
                return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def far_count(speeds): return sum(1 for v in speeds if abs(v) > 22)

def greedy_lcm_blocker(mag_lo, mag_hi, rng):
    """band-blocker: speeds = prime-power bundles, each in [mag_lo, mag_hi], covering many small moduli."""
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
    speeds = []; i = 0
    while i < len(primes) and len(speeds) < 11:
        v = 1;
        while i < len(primes) and v * primes[i] <= mag_hi:
            if v * primes[i] >= mag_lo or v > 1: v *= primes[i]
            i += 1
        if mag_lo <= v <= mag_hi: speeds.append(v)
        elif v > 1 and v < mag_lo:
            v *= rng.randint(max(1, mag_lo // v), max(1, mag_hi // v))
            if mag_lo <= v <= mag_hi: speeds.append(v)
    for q in range(2, 15):
        if len(speeds) >= 13: break
        if not any(s % q == 0 for s in speeds): speeds.append(q)
    while len(speeds) < 13: speeds.append(rng.randint(2, 22))
    return sorted(set(speeds))[:13]

if __name__ == "__main__":
    print("Dyadic magnitude ladder: worst smallest-lonely-q Q(k) vs level k=log2(max-speed).")
    print("=" * 84)
    print(f"{'k':>4} {'[2^k, 2^k+1)':>18} {'#families':>10} {'Q(k) worst':>11} {'Q(k)/k':>8} {'example (worst)':>15}")
    rows = []
    for k in range(5, 25):
        lo, hi = 2**k, 2**(k+1)
        if lo < 23: lo = 23
        rng = random.Random(1000 + k)
        worst = (0, None); ncov = 0
        # (a) greedy lcm band-blockers (discrepancy-extremal)
        for _ in range(3000):
            sp = greedy_lcm_blocker(lo, hi - 1, rng)
            if len(sp) == 13 and gcd_all(sp) == 1 and far_count(sp) >= 7 and is_covering(sp) and lo <= max(sp) < hi:
                q = smallest_lonely_q(sp, qmax=3000)
                ncov += 1
                if q and q > worst[0]: worst = (q, sp)
        # (b) random band-multiple families in the dyadic band
        bandmods = list(range(15, 90))
        for _ in range(40000):
            speeds = []
            for q in rng.sample(bandmods, rng.randint(6, 11)):
                m = q * rng.randint(max(1, lo // q), max(1, (hi - 1) // q))
                if lo <= m < hi: speeds.append(m)
            while len(set(speeds)) < 13:
                speeds.append(rng.randint(2, 22) if rng.random() < 0.35 else rng.randint(lo, hi - 1))
            speeds = sorted(set(speeds))[:13]
            if len(speeds) != 13 or gcd_all(speeds) != 1 or far_count(speeds) < 7 or not is_covering(speeds):
                continue
            if not (lo <= max(speeds) < hi): continue
            q = smallest_lonely_q(speeds, qmax=3000); ncov += 1
            if q and q > worst[0]: worst = (q, speeds)
        ratio = worst[0] / k if k else 0
        ex = str(worst[1][:4]) + "..." if worst[1] else "-"
        print(f"{k:>4} {'[%d,%d)'%(lo,hi):>18} {ncov:>10} {worst[0]:>11} {ratio:>8.2f} {ex:>15}")
        rows.append((k, worst[0]))
    print("\nQ(k) sequence:", [q for _, q in rows])
    # linear fit Q(k) ~ a*k + b
    ks = [k for k, _ in rows]; qs = [q for _, q in rows]
    n = len(ks); sk = sum(ks); sq = sum(qs); skk = sum(k*k for k in ks); skq = sum(k*q for k, q in rows)
    a = (n*skq - sk*sq) / (n*skk - sk*sk); b = (sq - a*sk) / n
    print(f"linear fit: Q(k) ~ {a:.2f}*k + {b:.1f}   (k = log2 max-speed; slope {a:.2f} per doubling)")
    print(f"=> if slope ~ const: q ~ {a:.1f}*log2(max-speed) -- the LRC all-scales/discrepancy cost is LOGARITHMIC")
    print("   (linear-in-k = the dyadic ladder; each level a finite band; magnitude split = 2-level coarsening)")
