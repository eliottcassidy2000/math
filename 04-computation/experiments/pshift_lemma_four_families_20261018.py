#!/usr/bin/env python3
"""p-shift lemma for four hypergeometric families (S43, HYP-9311).

Claim: for A(k) = prod of factorial ratios, w_k = d/dk log A(k) (Gamma continuation), a carry-free index j (every
factorial argument c j < p) and a digit a with c a < p for every c:
    A(ap + j) == A(a) A(j) (1 + a p w_j)   (mod p^2).
Families: A1 = (2k)!^3/k!^6, A2 = (2k)!(3k)!/k!^5, A3 = (2k)!(4k)!/(k!^4 (2k)!) * ... (C(2k,k)^2 C(4k,2k) = (4k)!/(k!^4)),
A4 = C(2k,k)C(3k,k)C(6k,3k) = (6k)!/(k!^3 (3k)!).  Checked exactly for primes 7 <= p <= 61.
Usage: python pshift_lemma_four_families_20261018.py
Session opus-2026-10-17-S43.
"""
from fractions import Fraction
from math import comb, factorial
import json, os

FAM = {
    'A1': (lambda k: comb(2 * k, k) ** 3, [(6, 2), (-6, 1)], 2),
    'A2': (lambda k: comb(2 * k, k) ** 2 * comb(3 * k, k), [(3, 3), (2, 2), (-5, 1)], 3),
    'A3': (lambda k: comb(2 * k, k) ** 2 * comb(4 * k, 2 * k), [(4, 4), (-4, 1)], 4),
    'A4': (lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k), [(6, 6), (-3, 3), (-3, 1)], 6),
}


def main():
    primes = [p for p in range(7, 62) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    H = [Fraction(0)]
    for i in range(1, 6 * 62 * 6 + 10):
        H.append(H[-1] + Fraction(1, i))
    res = {}
    for name, (Af, wt, cmax) in FAM.items():
        bad = 0
        checks = 0
        for p in primes:
            pp = p * p
            for a in range(0, (p - 1) // cmax + 1):
                for j in range(0, (p - 1) // cmax + 1):
                    w = sum(c * H[s * j] for c, s in wt)
                    lhs = Af(a * p + j) % pp
                    rhs = Fraction(Af(a) * Af(j)) * (1 + a * p * w)
                    r = rhs.numerator * pow(rhs.denominator, -1, pp) % pp
                    checks += 1
                    if (lhs - r) % pp:
                        bad += 1
        res[name] = dict(checks=checks, failures=bad)
        print('%s: A(ap+j) == A(a)A(j)(1 + ap w_j) mod p^2 for carry-free a, j, 7 <= p <= 61: %d checks, %d failures'
              % (name, checks, bad))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(res, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'pshift_lemma_four_families_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
