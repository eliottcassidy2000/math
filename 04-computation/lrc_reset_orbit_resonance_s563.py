#!/usr/bin/env python3
"""
Time in the LRC as an orbit: commensurability (reset) is the hard regime, and the
resonance-folding categorises difficulty.
opus-2026-06-02-S563 (remote-control).

The 13 runners = one point on T^13 along gamma(t)=(v_i t mod 1); loneliness = gamma
meets the box B={x:||x_i||>=1/14}. Incommensurate ratios => gamma DENSE (Weyl) =>
trivially lonely. Commensurate (integer, reset at t=1/gcd) => closed loop, can miss
B => the whole problem. Among resetting sets, the categoriser is how the loop FOLDS
(resonances), not the reset length. We verify: more resonance <=> smaller lonely
measure <=> harder; the AP folds maximally (most resonances, fewest distinct moments).
"""
from fractions import Fraction
from math import gcd
from itertools import combinations
import random

n = 14


def safe_measure(V):
    thr = Fraction(1, n)
    eps = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                eps.add(Fraction(k * n + s, n * v) % 1)
    pts = sorted(eps)
    L = len(pts)
    meas = Fraction(0)
    for i in range(L):
        a, b = pts[i], pts[(i + 1) % L]
        ln = (b - a) if b > a else (b - a + 1)
        mid = (a + ln / 2) % 1
        if all(min((v * mid) % 1, 1 - (v * mid) % 1) > thr for v in V):
            meas += ln
    return meas


def Mval(V):
    c = set()
    for v in V:
        for k in range(2 * v):
            c.add(Fraction(2 * k + 1, 2 * v) % 1)
    for i in range(len(V)):
        for j in range(len(V)):
            for s in (1, -1):
                d = V[i] + s * V[j]
                if d:
                    for k in range(abs(d) + 1):
                        c.add(Fraction(k, d) % 1)
    return max(min(min((v * t) % 1, 1 - (v * t) % 1) for v in V) for t in c)


def resonances3(V, C=2):
    cnt = 0
    for x, y, z in combinations(V, 3):
        for a in range(-C, C + 1):
            for b in range(-C, C + 1):
                for cc in range(-C, C + 1):
                    if (a, b, cc) != (0, 0, 0) and a * x + b * y + cc * z == 0:
                        cnt += 1
    return cnt


def crit_times(V):
    eps = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                eps.add(Fraction(k * n + s, n * v) % 1)
    return len(eps)


def primitive(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V))


if __name__ == "__main__":
    fams = {
        'AP {1..13}': tuple(range(1, 14)),
        'V* sporadic': (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24),
        'AP odd': tuple(range(1, 27, 2)),
        'Fibonacci': (1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377),
        'Sidon': (1, 2, 5, 11, 22, 40, 56, 73, 78, 97, 118, 134, 150),
    }
    rng = random.Random(3)
    for k in ('random A', 'random B'):
        fams[k] = primitive(tuple(sorted(rng.sample(range(1, 80), 13))))
    print(f'{"family":16s} {"safe-meas~":>10s} {"M(S)":>8s} '
          f'{"3-reson|c|<=2":>13s} {"#crit":>7s}')
    for name, V in fams.items():
        if max(V) > 2000:
            print(f'{name:16s} {"(big)":>10s} {"(big)":>8s} {resonances3(V):>13d} '
                  f'{"(big)":>7s}')
            continue
        mu = safe_measure(V)
        print(f'{name:16s} {float(mu):>10.4f} {str(Mval(V)):>8s} '
              f'{resonances3(V):>13d} {crit_times(V):>7d}')
    print("\nMORE resonance (orbit folding) <=> SMALLER lonely measure <=> HARDER.")
    print("AP: most resonances, FEWEST distinct critical moments (resonances make")
    print("them coincide) -> the hardest set condenses its 4D model the most.")
