#!/usr/bin/env python3
"""Rank-one prediction: the equal-time merge time of Haar y and y + E (3x+1, pair chain from (0, E)) is a sum of about
c log E debt excursions with tail t^(-1/2), hence of order (log E)^2 (a 1/2-stable sum), and
P(no merge by T) ~ (a + b log E) T^(-1/2).  Usage: NSAMP TMAX SEED"""
import sys, random, math
from fractions import Fraction as Fr
def par(e): return (e.numerator * pow(e.denominator, -1, 2)) % 2
def merge_time(E, tmax, rnd):
    k, e = 0, Fr(E); t = 0
    while t < tmax:
        s = par(e); b = rnd.getrandbits(1)
        if s == 0:
            e = e / 2 if b == 0 else (3 * e + 1 - Fr(3) ** k) / 2
        elif b == 0:
            k, e = k + 1, (3 * e + 1) / 2
        else:
            k, e = k - 1, (e - Fr(3) ** (k - 1)) / 2
        t += 1
        if k == 0 and e == 0: return t
    return None
if __name__ == '__main__':
    nsamp, tmax, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rnd = random.Random(seed)
    print(f"3x+1 pair chain from (0, E): {nsamp} samples, cap {tmax}")
    for j in (1, 2, 4, 8, 16, 32, 64, 128):
        E = 2**j + 1 if j > 1 else 1   # odd offsets (E = 2^j + 1), plus E = 1
        ts = [merge_time(E, tmax, rnd) for _ in range(nsamp)]
        done = sorted(x for x in ts if x is not None)
        q = {T: sum(1 for x in ts if x is None or x > T) / nsamp for T in (1000, 4000, 16000) if T <= tmax}
        med = done[len(ts) // 2] if len(done) > len(ts) // 2 else None
        L = math.log(E) if E > 1 else 1.0
        qs = "  ".join(f"sqrt(T)q({T})={math.sqrt(T)*v:6.2f}" for T, v in q.items())
        print(f"  E=2^{j}+1: median merge time {med}  median/(ln E)^2 = {med / L**2 if med else float('nan'):.2f}  {qs}", flush=True)
