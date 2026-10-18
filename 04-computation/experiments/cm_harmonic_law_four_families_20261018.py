#!/usr/bin/env python3
"""The CM harmonic log law for four hypergeometric families, with Sun's weights (S43; generalizes HYP-9310).

For A(k) in {C(2k,k)^3, C(2k,k)^2 C(3k,k), C(2k,k)^2 C(4k,2k), C(2k,k) C(3k,k) C(6k,3k)} let
    w_k = d/dk log A(k)  (Gamma continuation) = 6(H_2k - H_k), 3H_3k + 2H_2k - 5H_k, 4(H_4k - H_k), 3(2H_6k - H_3k - H_k),
which are exactly the harmonic weights of Z.-W. Sun's Conjecture 29 (arXiv:2210.07238) and exactly the p-shift weights:
A(ap + j) == A(a) A(j) (1 + ap w_j) (mod p^2) for carry-free j.  Over the full period the poles of w_k at carry indices
are cancelled by the zeros of A(k), so t_k w_k is p-integral for every k < p (checked below).
PRE-REGISTERED (written before this run): at every CM parameter m found by cm_detector_lucas_block_scan_20261018.py
(33 values over four families) and every prime 7 <= p <= PMAX not dividing 6m,
    (L_p^A)   sum_{k<p} A(k) m^-k w_k == q_p(m) sum_{k<p} A(k) m^-k   (mod p),
and it fails at generic m.  (For A = C(2k,k)^3 this is HYP-9310, since the terms with (p-1)/2 < k < p vanish mod p.)
Exact rational arithmetic (fractions), so PMAX is modest.
Usage: python cm_harmonic_law_four_families_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json
from fractions import Fraction
from math import comb


def primes_between(a, b):
    return [p for p in range(a, b + 1) if p > 1 and all(p % d for d in range(2, int(p ** 0.5) + 1))]


FAM = {
    'A1': (lambda k: comb(2 * k, k) ** 3, [(6, 2), (-6, 1)]),
    'A2': (lambda k: comb(2 * k, k) ** 2 * comb(3 * k, k), [(3, 3), (2, 2), (-5, 1)]),
    'A3': (lambda k: comb(2 * k, k) ** 2 * comb(4 * k, 2 * k), [(4, 4), (-4, 1)]),
    'A4': (lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k), [(6, 6), (-3, 3), (-3, 1)]),
}
SURV = {
    'A1': [-512, -64, -8, 1, 16, 64, 256, 4096],
    'A2': [-8640, -1728, -192, -27, 8, 64, 108, 216, 1458, 3375],
    'A3': [-82944, -12288, -3969, -1024, -144, 81, 256, 648, 2304, 20736],
    'A4': [-32768, -3375, 1728, 8000, 54000],
}
GEN = [2, 3, 5, 7, 10, 100]


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    primes = primes_between(7, PMAX)
    Hmax = 6 * PMAX + 6
    H = [Fraction(0)]
    for i in range(1, Hmax + 1):
        H.append(H[-1] + Fraction(1, i))
    out = {}
    for fam, (Af, wt) in FAM.items():
        Avals = [Af(k) for k in range(PMAX)]
        TW = [Fraction(Avals[k]) * sum(c * H[a * k] for c, a in wt) for k in range(PMAX)]
        res = {}
        for m in SURV[fam] + GEN:
            fails = []
            nonint = 0
            tested = 0
            for p in primes:
                if (6 * m) % p == 0:
                    continue
                tested += 1
                S = 0
                Sw = 0
                pp = p
                minv = pow(m % pp, -1, pp)
                for k in range(p):
                    term = Avals[k]
                    tw = TW[k]
                    if tw.denominator % p == 0:
                        nonint += 1
                        continue
                    mk = pow(minv, k, pp)
                    S = (S + term * mk) % pp
                    Sw = (Sw + tw.numerator * pow(tw.denominator, -1, pp) * mk) % pp
                q = (pow(m % (p * p), p - 1, p * p) - 1) // p % p
                if (Sw - q * S) % p:
                    fails.append(p)
            res[str(m)] = dict(cm=m in SURV[fam], tested=tested, fails=fails, nonintegral_terms=nonint)
            print('%s m = %7d (%s): tested %d primes, failures %s%s' % (
                fam, m, 'CM ' if m in SURV[fam] else 'gen', tested,
                'none' if not fails else '%d, first %s' % (len(fails), fails[:6]),
                '' if not nonint else '  [non-p-integral terms: %d]' % nonint), flush=True)
        out[fam] = res
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_harmonic_law_four_families_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
