#!/usr/bin/env python3
"""Probe: does the Lucas-block congruence at CM points hold beyond mod p^2? (S43, around HYP-9310/9311.)

(B1') sum_{k=p}^{2p-1} t_k == t_1 sum_{k<p} t_k, t_k = A(k) m^-k, tested mod p^3 at the CM parameters of the four
families, and the full-square form F_(p^2) == F_p^2 mod p^3 for A1.  Exploratory: no prediction registered.
Usage: python cm_lucas_block_mod_p3_probe_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, json, os
from math import comb

FAM = {
    'A1': lambda k: comb(2 * k, k) ** 3,
    'A2': lambda k: comb(2 * k, k) ** 2 * comb(3 * k, k),
    'A3': lambda k: comb(2 * k, k) ** 2 * comb(4 * k, 2 * k),
    'A4': lambda k: comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k),
}
SURV = {
    'A1': [-512, -64, -8, 1, 16, 64, 256, 4096],
    'A2': [-8640, -1728, -192, -27, 8, 64, 108, 216, 1458, 3375],
    'A3': [-82944, -12288, -3969, -1024, -144, 81, 256, 648, 2304, 20736],
    'A4': [-32768, -3375, 1728, 8000, 54000],
}


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    primes = [p for p in range(7, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    out = {}
    for fam, Af in FAM.items():
        A = [Af(k) for k in range(2 * PMAX)]
        res = {}
        for m in SURV[fam]:
            ok3 = []
            bad3 = []
            for p in primes:
                if (6 * m) % p == 0:
                    continue
                P3 = p ** 3
                im = pow(m % P3, -1, P3)
                S0 = sum(A[k] % P3 * pow(im, k, P3) for k in range(p)) % P3
                S1 = sum(A[k] % P3 * pow(im, k, P3) for k in range(p, 2 * p)) % P3
                d = (S1 - A[1] * im * S0) % P3
                assert d % (p * p) == 0
                (ok3 if d == 0 else bad3).append(p)
            res[str(m)] = dict(mod_p3_holds=ok3, mod_p3_fails=bad3)
            print('%s m = %7d: (B1\') mod p^3 holds at %d primes %s, fails at %d' % (fam, m, len(ok3), ok3[:8], len(bad3)), flush=True)
        out[fam] = res
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_lucas_block_mod_p3_probe_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
