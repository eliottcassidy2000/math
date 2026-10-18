#!/usr/bin/env python3
"""Lucas-block form of the CM harmonic law (S43, HYP-9310).

For 0 <= j <= n = (p-1)/2 and every digit a with 2a < p:  C(2(ap+j), ap+j) == C(2a,a) C(2j,j) (1 + 2ap(H_2j - H_j))
(mod p^2).  Hence with t_k = C(2k,k)^3 m^-k and F_m = sum_{k<=n} t_k:
    block_a := sum_{j<=n} t_(ap+j) == C(2a,a)^3 m^-a m^(-a(p-1)) (F_m + 6 a p h_m)   (mod p^2),
so the law 6 h_m == q_p(m) F_m (mod p) is EQUIVALENT to the second-order Lucas factorization of the first block
    (B1)  sum_{k=p}^{p+n} C(2k,k)^3 m^-k == (8/m) sum_{k=0}^{n} C(2k,k)^3 m^-k   (mod p^2),
and implies block_a == C(2a,a)^3 m^-a F_m for every digit a, hence F_(p^2) == F_p^2 (mod p^2) with F_N = sum_{k<N} t_k.
Checks: (1) the binomial lemma; (2) (B1) at the 8 CM values and at generic m; (3) F_(p^2) == F_p^2 for small p.
Usage: python cm_law_lucas_block_form_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json
from math import comb


def primes_upto(N):
    return [p for p in range(5, N + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    primes = primes_upto(PMAX)
    out = {}
    # (1) binomial lemma
    bad = 0
    for p in primes[:12]:
        n = (p - 1) // 2
        pp = p * p
        H = [0] * (2 * n + 1)
        for i in range(1, 2 * n + 1):
            H[i] = (H[i - 1] + pow(i, -1, pp)) % pp
        for a in range(0, (p + 1) // 2):
            for j in range(n + 1):
                lhs = comb(2 * (a * p + j), a * p + j) % pp
                rhs = comb(2 * a, a) * comb(2 * j, j) * (1 + 2 * a * p * (H[2 * j] - H[j])) % pp
                if (lhs - rhs) % pp:
                    bad += 1
    print('(1) C(2(ap+j),ap+j) == C(2a,a)C(2j,j)(1 + 2ap(H_2j - H_j)) mod p^2, p <= %d, 2a < p, j <= n: %s'
          % (primes[11], 'holds' if not bad else '%d failures' % bad))
    out['lemma_failures'] = bad
    # (2) (B1)
    CM = [1, -8, 16, -64, 64, 256, -512, 4096]
    GEN = [3, 5, 7, -3, 10]
    res = {}
    for m in CM + GEN:
        fails = []
        for p in primes:
            if m % p == 0:
                continue
            n = (p - 1) // 2
            pp = p * p
            im = pow(m % pp, -1, pp)
            F = sum(comb(2 * k, k) ** 3 * pow(im, k, pp) for k in range(n + 1)) % pp
            B1 = sum(comb(2 * k, k) ** 3 * pow(im, k, pp) for k in range(p, p + n + 1)) % pp
            if (B1 - 8 * im * F) % pp:
                fails.append(p)
        res[str(m)] = fails
        print('(2) (B1) at m = %5d (%s): %s' % (m, 'CM ' if m in CM else 'gen', 'holds for all p <= %d' % PMAX if not fails else 'fails at %d primes, first %s' % (len(fails), fails[:6])))
    out['B1'] = res
    # (3) F_(p^2) == F_p^2 mod p^2 for small p
    res3 = {}
    for m in CM + GEN:
        fails = []
        for p in primes_upto(40):
            if m % p == 0:
                continue
            pp = p * p
            im = pow(m % pp, -1, pp)
            Fp = sum(comb(2 * k, k) ** 3 * pow(im, k, pp) for k in range(p)) % pp
            Fp2 = 0
            c = 1
            w = 1
            for k in range(p * p):
                if k:
                    c = c * (2 * k) * (2 * k - 1) // (k * k)
                    w = w * im % pp
                Fp2 = (Fp2 + (c % pp) ** 3 * w) % pp
            if (Fp2 - Fp * Fp) % pp:
                fails.append(p)
        res3[str(m)] = fails
        print('(3) F_(p^2) == F_p^2 mod p^2 at m = %5d: %s' % (m, 'holds for p <= 37' if not fails else 'fails at %s' % fails[:6]))
    out['Fp2'] = res3
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_law_lucas_block_form_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
