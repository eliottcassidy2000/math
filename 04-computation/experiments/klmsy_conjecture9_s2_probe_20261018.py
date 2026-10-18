#!/usr/bin/env python3
"""Exploratory probe of KLMSY Conjecture 9 at s = 2 (S43).  No prediction registered; small primes only.

Conjecture 9 (s = 2, M = 1): F_{(p^2-1)/2}(lambda) == ((lambda-1)/p) alpha^2 F_{(p-1)/2}(lambda)  (mod p^4).
Without alpha: at ordinary p, F_{(p-1)/2} == eps alpha^2 (mod p^2) only, so we test the alpha-free consequences
    (a) F_{(p^2-1)/2} == F_{(p-1)/2}^2 (mod p^2)  [implied by s = 1 block factorization; sanity check]
    (b) F_{(p^2-1)/2} * F_0 == F_{(p-1)/2} * F_{(p-1)/2} (mod p^3)?  (a guess at the next order)
and record the p-adic valuation of F_{(p^2-1)/2} - F_{(p-1)/2}^2 at CM vs generic lambda.
Usage: python klmsy_conjecture9_s2_probe_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, json, os
from fractions import Fraction
from math import comb


def vp(x, p):
    if x == 0:
        return 99
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    lams = {('CM', m): Fraction(64, m) for m in (1, -8, 16, -64, 64, 256, -512, 4096)}
    for g in (3, 5, 7):
        lams[('gen', g)] = Fraction(64, g)
    Kmax = (PMAX * PMAX - 1) // 2 + 1
    c = [1] * (Kmax + 1)
    for k in range(1, Kmax + 1):
        c[k] = c[k - 1] * (2 * k) * (2 * k - 1) // (k * k)
    out = {}
    for key, lam in lams.items():
        vals = {}
        for p in primes:
            if lam.numerator % p == 0 or lam.denominator % p == 0:
                continue
            P6 = p ** 6
            l = lam.numerator * pow(lam.denominator, -1, P6) % P6
            w = l * pow(64, -1, P6) % P6
            def F(N):
                s = 0
                t = 1
                for k in range(N + 1):
                    s = (s + pow(c[k] % P6, 3, P6) * t) % P6
                    t = t * w % P6
                return s
            n = (p - 1) // 2
            F1 = F(n)
            F2 = F((p * p - 1) // 2)
            d = (F2 - F1 * F1) % P6
            vals[p] = dict(v_diff=min(vp(d, p), 6), v_F1=min(vp(F1, p), 6))
        out['%s %s' % key] = vals
        print('%s m = %5d: v_p(F_{(p^2-1)/2} - F_{(p-1)/2}^2) = %s' % (key[0], key[1],
              ', '.join('%d:%d%s' % (p, v['v_diff'], '' if v['v_F1'] == 0 else '(ss)') for p, v in vals.items())), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'klmsy_conjecture9_s2_probe_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
