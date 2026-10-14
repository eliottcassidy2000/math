#!/usr/bin/env python3
"""The K_4 character sum S(q) = sum_{x,y in F_q} chi(x (x-1) y (y-1) (y-x)) for primes q, by residue class mod 4.

S(q) is the oriented K_4 sum (vertices 0, 1, x, y; all six arcs lower -> higher) that first appears in the Paley
cluster expansion when glued path quotients acquire a K_4 minor.  For q = 3 mod 4 (Paley tournaments exist) we test
whether S(q) is a polynomial in q; for q = 1 mod 4 we compare with q = a^2 + b^2 (a odd), the Gaussian (CM by Q(i))
signature.  Session opus-2026-10-14-S34.
"""
import numpy as np


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, n + 1) if s[i]]


def S(q):
    chi = np.zeros(q, dtype=np.int64)
    sq = set((x * x) % q for x in range(1, q))
    for x in range(1, q):
        chi[x] = 1 if x in sq else -1
    xs = np.arange(q)
    fx = chi[(xs * (xs - 1)) % q]              # chi(x(x-1))
    tot = 0
    for y in range(q):
        cy = chi[(y * (y - 1)) % q]
        if cy == 0:
            continue
        tot += cy * int(np.dot(fx, chi[(y - xs) % q]))
    return tot


def two_squares(q):
    for a in range(1, int(q ** 0.5) + 1, 2):
        b2 = q - a * a
        b = int(round(b2 ** 0.5))
        if b * b == b2:
            return a, b
    return None


for cls in (3, 1):
    print(f'primes q = {cls} mod 4:')
    for q in primes_upto(200):
        if q % 4 != cls:
            continue
        s = S(q)
        extra = ''
        if cls == 1:
            a, b = two_squares(q)
            extra = f'  q = {a}^2 + {b}^2;  (S + ...)/a^2 checks: S/(4a^2 - 2q) = {s / (4 * a * a - 2 * q):.4f}' if 4 * a * a != 2 * q else ''
        print(f'   q = {q:3d}: S = {s:7d},  S/q = {s / q:9.4f},  S + 3q = {s + 3 * q}{extra}')
