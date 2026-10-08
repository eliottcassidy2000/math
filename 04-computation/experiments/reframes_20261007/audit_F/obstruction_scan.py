#!/usr/bin/env python3
"""audit_F (C1b): scan maps for the affine congruence obstruction of obstruction_check.py:
a prime p not dividing d*prod(m_i) and c mod p with (m_i - d) c + r_i = 0 mod p for all i (then y, y+e never meet
when p does not divide e).  Scans p < 2000 (and p^2 for p < 50).  Light."""
import math
from maps_check import MAPS


def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def obstructions(d, br, P=2000):
    out = []
    prod = d
    for m, r in br:
        prod *= m
    for p in primes(P):
        if prod % p == 0:
            continue
        for q in ([p, p * p] if p < 50 else [p]):
            for c in range(q):
                if all(((m - d) * c + r) % q == 0 for m, r in br):
                    out.append((q, c))
                    break
    return out


if __name__ == '__main__':
    extra = {
        '3x+5':            (2, [(1, 0), (3, 5)]),
        'x+5 (rank 0)':    (2, [(1, 0), (1, 5)]),
        'Z3_115 r=0,2,2':  (3, [(1, 0), (1, 2), (5, 2)]),
        'Z3_111 r=0,2,4':  (3, [(1, 0), (1, 2), (1, 4)]),
        'Z3_125 r=0,7,14': (3, [(1, 0), (2, 7), (5, 14)]),
        'Z3_125 r=9,1,5':  (3, [(1, 9), (2, 1), (5, 5)]),
    }
    for name, val in list(MAPS.items()) + list(extra.items()):
        d, br = val[0], val[1]
        print(f"{name:18s} d={d} br={br}: obstructions (modulus, c) with p < 2000: {obstructions(d, br)}")
