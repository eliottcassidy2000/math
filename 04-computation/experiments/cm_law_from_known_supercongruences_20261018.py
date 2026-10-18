#!/usr/bin/env python3
"""The CM harmonic law (HYP-9310) as a corollary of the known mod-p^2 CM supercongruences at both ends of each
reflection (2-isogeny) pair (S43).

With n = (p-1)/2, a_k = C(2k,k)/4^k, F(z) = sum_{k<=n} a_k^3 z^k, G(w) = sum_{k<=n} a_k^3 w^k Delta_k,
Delta_k = O_(n-k) - O_k, the shadow a_(n-k) == (-1)^n a_k (1 + p Delta_k) (mod p^2) gives the exact-mod-p^2 reflection
    (R)  F(z) == (-1)^n z^n [ F(1/z) + 3p G(1/z) ]   (mod p^2).
If F(z) == eps1 a(p) and F(1/z) == eps2 a(p) (mod p^2) are known (the CM supercongruences), (R) determines G(1/z) mod p,
and h = -G + q_p(2) F (mod p) then gives the harmonic law 6h == q_p(m) F at m = 64/z.  This script checks:
(1) (R) at z = 64, -8, 4 (m = 1, -8, 16) for 5 <= p <= PMAX;
(2) the derived value of G(1/z) mod p, and hence the law at m' = 4096/m, from the two supercongruence values alone.
Usage: python cm_law_from_known_supercongruences_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, json, os
from fractions import Fraction
from math import comb


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    out = {}
    for m, z in ((1, Fraction(64)), (-8, Fraction(-8)), (16, Fraction(4))):
        badR = []
        badLaw = []
        for p in primes:
            if (m * 6) % p == 0 or z.numerator % p == 0:
                continue
            n = (p - 1) // 2
            pp = p * p
            def md(x):
                x = Fraction(x)
                return x.numerator * pow(x.denominator, -1, pp) % pp
            a3 = [md(Fraction(comb(2 * k, k), 4 ** k) ** 3) for k in range(n + 1)]
            O = [Fraction(0)]
            for j in range(1, n + 1):
                O.append(O[-1] + Fraction(1, 2 * j - 1))
            Dl = [md(O[n - k] - O[k]) for k in range(n + 1)]
            zz = md(z)
            zi = md(1 / z)
            F = sum(a3[k] * pow(zz, k, pp) for k in range(n + 1)) % pp
            Fi = sum(a3[k] * pow(zi, k, pp) for k in range(n + 1)) % pp
            Gi = sum(a3[k] * pow(zi, k, pp) * Dl[k] for k in range(n + 1)) % pp
            sgn = (-1) ** n
            if (F - sgn * pow(zz, n, pp) * (Fi + 3 * p * Gi)) % pp:
                badR.append(p)
            # harmonic law at m' = 4096/m (i.e. at 1/z) computed directly: 6 h(1/z) == q_p(m') F(1/z) mod p
            mp = Fraction(4096, m)
            H = [Fraction(0)]
            for i in range(1, 2 * n + 1):
                H.append(H[-1] + Fraction(1, i))
            h = sum(a3[k] * pow(zi, k, pp) * md(H[2 * k] - H[k]) for k in range(n + 1)) % p
            qm = (pow(md(mp) % pp, p - 1, pp) - 1) // p % p
            if (6 * h - qm * Fi) % p:
                badLaw.append(p)
        out[str(m)] = dict(reflection_failures=badR, law_failures_at_4096_over_m=badLaw)
        print('m = %3d (z = %s): reflection (R) mod p^2 failures %s; law at m\' = 4096/m failures %s'
              % (m, z, badR if badR else 'none', badLaw if badLaw else 'none'), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_law_from_known_supercongruences_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
