#!/usr/bin/env python3
"""A 'log m / 6' law for the harmonic twin of sum C(2k,k)^3 / m^k at CM points (S43, after (H)).

Twin census -> (H): sum_{k<=n} C(2k,k)^3 (H_2k - H_k) == 0 (mod p), n = (p-1)/2, verified p < 3000.  Reflection
k -> n-k turns m = 1 into m = 4096 and gives sum_{k<=n} C(2k,k)^3/4096^k (H_2k - H_k) == 2 q_p(2) sum C(2k,k)^3/4096^k.
Both are the case m = 1, 4096 of the guessed law
    (L_p)  sum_{k<=n} C(2k,k)^3 m^-k (H_2k - H_k) == (q_p(m)/6) sum_{k<=n} C(2k,k)^3 m^-k   (mod p),
q_p(m) = (m^(p-1) - 1)/p.  Pre-registered prediction (written before this run): (L_p) holds at every one of Sun's
CM values m in {1, -8, 16, -64, 64, 256, -512, 4096} (j = (m-16)^3/m CM) for all p not dividing 6m, and FAILS for
generic m (3, 5, 7, 10, 100, -3).  Archimedean analogue tested numerically:
    (L_inf) sum_k C(2k,k)^3 m^-k (H_2k - H_k) = (log|m| / 6) sum_k C(2k,k)^3 m^-k    for CM |m| > 64.
Usage: python cm_harmonic_log_law_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json
from math import comb, log
from fractions import Fraction


def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0] = s[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]


def test_m(m, primes):
    """returns (#primes tested, #fail, #primes with F == 0 mod p, list of first failures)."""
    fails = []
    zero = 0
    tested = 0
    for p in primes:
        if (6 * m) % p == 0:
            continue
        n = (p - 1) // 2
        invm = pow(m % p, p - 2, p)
        c = 1          # C(2k,k) mod p
        H = [0] * (2 * n + 1)
        for i in range(1, 2 * n + 1):
            H[i] = (H[i - 1] + pow(i, p - 2, p)) % p
        F = 0
        h = 0
        w = 1
        for k in range(n + 1):
            if k:
                c = c * (2 * k) * (2 * k - 1) % p * pow(k * k, p - 2, p) % p
                w = w * invm % p
            t = pow(c, 3, p) * w % p
            F = (F + t) % p
            h = (h + t * (H[2 * k] - H[k])) % p
        q = (pow(m % (p * p), p - 1, p * p) - 1) // p % p
        tested += 1
        if F == 0:
            zero += 1
        if (6 * h - q * F) % p:
            fails.append(p)
    return tested, len(fails), zero, fails[:8]


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    primes = primes_upto(PMAX)[2:]
    out = {}
    CM = [1, -8, 16, -64, 64, 256, -512, 4096]
    GEN = [3, 5, 7, 10, 100, -3, 2, -1, 32, 128]
    print('(L_p): 6 sum_{k<=n} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=n} C(2k,k)^3 m^-k (mod p), 5 <= p <= %d' % PMAX)
    for m in CM + GEN:
        t, nf, z, ff = test_m(m, primes)
        tag = 'CM ' if m in CM else 'gen'
        print('  %s m = %6d : tested %4d primes, failures %4d (F == 0 mod p at %4d primes)%s'
              % (tag, m, t, nf, z, '' if not ff else ', first failures %s' % ff))
        out[str(m)] = dict(cm=m in CM, tested=t, failures=nf, F_zero=z, first_failures=ff)
    # archimedean analogue
    print('(L_inf): sum_k C(2k,k)^3 m^-k (H_2k - H_k)  vs  (log|m|/6) sum_k C(2k,k)^3 m^-k')
    arch = {}
    from decimal import Decimal, getcontext
    getcontext().prec = 50
    for m in [-512, 256, 4096, -4096, 1000, 512]:
        F = Decimal(0)
        h = Decimal(0)
        c = Decimal(1)
        Hk = Decimal(0)
        H2k = Decimal(0)
        term = Decimal(1)
        K = 400
        for k in range(K):
            if k:
                c = c * (2 * k) * (2 * k - 1) / (k * k)
                Hk += Decimal(1) / k
                H2k += Decimal(1) / (2 * k - 1) + Decimal(1) / (2 * k)
                term = term / m
            t = c ** 3 * term
            F += t
            h += t * (H2k - Hk)
        ratio = h / F
        target = Decimal(log(abs(m))) / 6
        print('  m = %6d: h/F = %.15f ; log|m|/6 = %.15f ; diff %.3e' % (m, ratio, target, ratio - target))
        arch[str(m)] = dict(ratio=str(ratio), log_over_6=str(target))
    out['archimedean'] = arch
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_harmonic_log_law_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
