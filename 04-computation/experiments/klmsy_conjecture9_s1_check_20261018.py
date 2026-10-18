#!/usr/bin/env python3
"""KLMSY Conjecture 9 at s = 1 from THM-4686 (S43).

Kibelbek-Long-Moss-Sheller-Yuan (J. Number Theory 2016) Conjecture 9: for CM lambda, odd M and s >= 1,
    F_{(M p^s - 1)/2}(lambda) == ((lambda-1)/p) alpha^2 F_{(M p^(s-1) - 1)/2}(lambda)   (mod p^(2s)),
F_N(lambda) = sum_{k<=N} ((1/2)_k/k!)^3 lambda^k.  They prove M = s = 1 (Theorem 2).  THM-4686's law plus the p-shift lemma
give block factorization for every base-p digit, hence for every odd M:
    (D)  F_{(Mp-1)/2} == F_{(p-1)/2} * F_{(M-1)/2}   (mod p^2),
which with Theorem 2 is Conjecture 9 at s = 1.  This script checks (D) at the eight CM values lambda = 64/m,
m in {1,-8,16,-64,64,256,-512,4096}, for odd M <= 11 and primes 5 <= p <= PMAX, and at generic lambda as a control.
Usage: python klmsy_conjecture9_s1_check_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, json, os
from fractions import Fraction
from math import comb


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    lambdas = {('CM', m): Fraction(64, m) for m in (1, -8, 16, -64, 64, 256, -512, 4096)}
    for g in (3, 5, 7, 10):
        lambdas[('gen', g)] = Fraction(64, g)
    Ms = [3, 5, 7, 9, 11]
    out = {}
    Kmax = (11 * PMAX - 1) // 2 + 1
    c = [1] * (Kmax + 1)
    for k in range(1, Kmax + 1):
        c[k] = c[k - 1] * (2 * k) * (2 * k - 1) // (k * k)
    for key, lam in lambdas.items():
        res = {}
        for M in Ms:
            fails = []
            for p in primes:
                if lam.numerator % p == 0 or lam.denominator % p == 0:
                    continue
                pp = p * p
                l = lam.numerator * pow(lam.denominator, -1, pp) % pp
                inv64 = pow(64, -1, pp)
                w = l * inv64 % pp          # a_k^3 lambda^k = C(2k,k)^3 (lambda/64)^k
                def F(N):
                    s = 0
                    t = 1
                    for k in range(N + 1):
                        s = (s + (c[k] % pp) ** 3 % pp * t) % pp
                        t = t * w % pp
                    return s
                n = (p - 1) // 2
                lhs = F((M * p - 1) // 2)
                rhs = F(n) * F((M - 1) // 2) % pp
                if (lhs - rhs) % pp:
                    fails.append(p)
            res[M] = fails
        out['%s %s' % key] = {str(M): v for M, v in res.items()}
        print('%s m = %5d (lambda = %s): ' % (key[0], key[1], lam) + '; '.join(
            'M=%d %s' % (M, 'ok' if not res[M] else 'fails %d (first %s)' % (len(res[M]), res[M][:4])) for M in Ms), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'klmsy_conjecture9_s1_check_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
