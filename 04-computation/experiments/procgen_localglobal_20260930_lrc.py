#!/usr/bin/env python3
"""procgen_localglobal_20260930_lrc.py -- LRC side of the local/global lane (helpers + experiments).

Convention: a speed set V has k distinct positive integer speeds; the LRC threshold is
delta = 1/(k+1) (so LRC(14) = 13 speeds, delta = 1/14).  M(V) = max_t min_v ||t v||.
All loneliness tests are exact integer arithmetic.

  M_exact(V)            exact M(V) (maximum over the candidate times c/(v_i+v_j), which contain
                        every local maximum of t -> min_v ||t v||), as a Fraction, plus maximisers
  lonely_count(V,q)     #{m in [1,q-1] : (k+1) * |m v mod q|_centered >= q for all v}
  least_lonely_den      least q with a lonely time m/q (gcd(m,q)=1 is automatic for the least q)
  lonely_primes         primes q <= Q with a lonely time m/q
  covers(V,q)           some speed divisible by q   (the length-1 relation in Lambda_q(V))
stdout only; checks raise on failure.
"""
import math
import random
from fractions import Fraction

import numpy as np


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(x) for x in np.nonzero(s)[0]]


PRIMES = primes_upto(20000)


def M_exact(V):
    """exact max_t min_v ||t v|| for integer speeds V (Fraction) and the list of maximisers in (0,1/2]."""
    V = sorted(set(int(v) for v in V))
    Va = np.array(V, dtype=np.int64)
    best = Fraction(0)
    arg = []
    dens = sorted(set(V[i] + V[j] for i in range(len(V)) for j in range(i, len(V))))
    for d in dens:
        c = np.arange(1, d // 2 + 1, dtype=np.int64)  # t = c/d in (0, 1/2]
        r = (c[:, None] * Va[None, :]) % d
        num = np.minimum(r, d - r).min(axis=1)  # min_v ||c v / d|| * d
        i = int(np.argmax(num))
        val = Fraction(int(num[i]), d)
        if val > best:
            best = val
            arg = [Fraction(int(cc), d) for cc in c[num * best.denominator == best.numerator * d]]
        elif val == best:
            arg += [Fraction(int(cc), d) for cc in c[num * best.denominator == best.numerator * d]]
    arg = sorted(set(arg))
    return best, arg


def lonely_count(V, q, k=None):
    """number of m in [1, q-1] with ||m v/q|| >= 1/(k+1) for all v (exact)."""
    if k is None:
        k = len(V)
    Va = np.array(V, dtype=np.int64) % q
    if np.any(Va == 0):
        return 0
    m = np.arange(1, q, dtype=np.int64)
    r = (m[:, None] * Va[None, :]) % q
    d = np.minimum(r, q - r)
    return int(np.count_nonzero(((k + 1) * d >= q).all(axis=1)))


def least_lonely_den(V, Qmax=5000, k=None):
    for q in range(2, Qmax + 1):
        if lonely_count(V, q, k) > 0:
            return q
    return None


def lonely_primes(V, Qmax=500, k=None):
    return [q for q in PRIMES if q <= Qmax and lonely_count(V, q, k) > 0]


def least_lonely_prime(V, Qmax=20000, k=None):
    for q in PRIMES:
        if q > Qmax:
            return None
        if lonely_count(V, q, k) > 0:
            return q
    return None


def covers(V, q):
    return any(v % q == 0 for v in V)


def divisor_complete(V, k=None):
    """every q in {2..k+1} divides some speed (the repo's 'covering' class)."""
    if k is None:
        k = len(V)
    return all(covers(V, q) for q in range(2, k + 2))


def prime_cover_depth(V):
    """largest P such that every prime <= P divides some speed (0 if 2 is uncovered)."""
    P = 0
    for q in PRIMES:
        if covers(V, q):
            P = q
        else:
            return P
    return P


def smooth_upto(N, S):
    out = [1]
    for p in S:
        new = []
        for x in out:
            y = x
            while y <= N:
                new.append(y)
                y *= p
        out = sorted(set(new))
    return out


if __name__ == "__main__":
    # quick self-test
    for k in range(3, 9):
        M, arg = M_exact(range(1, k + 1))
        assert M == Fraction(1, k + 1), (k, M)
    M, arg = M_exact(list(range(1, 13)) + [182])
    assert M == Fraction(14, 183), M
    print("lrc helpers OK", M, arg)
