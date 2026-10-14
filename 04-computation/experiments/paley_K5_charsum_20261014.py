#!/usr/bin/env python3
"""The oriented K_5 character sum T(q) = sum_{x,y,z} chi(x(x-1) y(y-1) z(z-1) (y-x)(z-x)(z-y)) for primes q = 3 mod 4.

Vertices 0, 1, x, y, z with all ten arcs oriented lower -> higher: the sum that controls 5-vertex subtournament
counts of the Paley tournament beyond what double regularity fixes.  If T(q) is not a polynomial in q, a genuine
Hasse-Weil (Deligne) term enters the Paley tournament at five vertices; we print T(q), T(q)/q^(3/2) and a polynomial
test.  Session opus-2026-10-14-S34.
"""
import numpy as np
import numba as nb


@nb.njit(cache=True)
def T(q, chi):
    tot = 0
    for x in range(2, q):
        cx = chi[(x * (x - 1)) % q]
        for y in range(2, q):
            if y == x:
                continue
            cxy = cx * chi[(y * (y - 1)) % q] * chi[(y - x) % q]
            for z in range(2, q):
                if z == x or z == y:
                    continue
                tot += cxy * chi[(z * (z - 1)) % q] * chi[(z - x) % q] * chi[(z - y) % q]
    return tot


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, n + 1) if s[i]]


vals = {}
for q in primes_upto(160):
    if q % 4 != 3:
        continue
    chi = np.zeros(q, dtype=np.int64)
    sq = set((x * x) % q for x in range(1, q))
    for x in range(1, q):
        chi[x] = 1 if x in sq else -1
    t = T(q, chi)
    vals[q] = t
    print(f'   q = {q:3d}: T = {t:9d},  T/q = {t / q:10.3f},  T/q^2 = {t / q ** 2:8.4f},  (T - q(q-?)...)  T mod q = {t % q}', flush=True)
