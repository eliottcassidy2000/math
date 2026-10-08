#!/usr/bin/env python3
"""q_p = P(y and y+1 merge at equal Terras time) for the px+1 map on Z_2 (Haar y), p = 5, 7, 9, 11, 13.
Integer chain: state (k, E) with e = E / p^max(0,-k) (E integer). Merges for p >= 5 happen early (the
normalized offset grows geometrically afterwards), so T = 200 is ample; we also record the latest merge time.
Usage: NSAMP SEED"""
import sys, random
from fractions import Fraction as Fr
def run(p, nsamp, T, rnd):
    merged = 0; latest = 0; times = []
    for _ in range(nsamp):
        k, e = 0, Fr(1)
        for t in range(1, T + 1):
            sig = (e.numerator * pow(e.denominator, -1, 2)) % 2
            b = rnd.getrandbits(1)
            if sig == 0:
                e = e / 2 if b == 0 else (p * e + 1 - Fr(p) ** k) / 2
            elif b == 0:
                k, e = k + 1, (p * e + 1) / 2
            else:
                k, e = k - 1, (e - Fr(p) ** (k - 1)) / 2
            if k == 0 and e == 0:
                merged += 1; latest = max(latest, t); times.append(t); break
            if abs(e) > Fr(p) ** (max(k, 0)) * 10**30: break   # |f| > 1e30: unreachable back (Lemma), stop
    return merged, latest, times
if __name__ == '__main__':
    nsamp, seed = int(sys.argv[1]), int(sys.argv[2])
    rnd = random.Random(seed)
    for p in (5, 7, 9, 11, 13):
        m, latest, times = run(p, nsamp, 400, rnd)
        q = m / nsamp
        se = (q * (1 - q) / nsamp) ** 0.5
        print(f"p={p:2d}: merged {m}/{nsamp}  q_p = {q:.5f} +- {se:.5f}; latest merge time {latest}; merge-time histogram (first 12): {sorted(times)[:12]}", flush=True)
