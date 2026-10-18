#!/usr/bin/env python3
"""P9 of PREREGISTRATION.md (fifth pre-registration, pushed in 1846da553b before this script was written).
kind-pasteur-2026-10-18-S2.  Lucas-block detector (B1') on weight-3 (Calabi-Yau threefold) families, 0 < |m| <= 10^4.
Also: defect valuations at the survivors (split by F_p(1/m) mod p) and the non-ordinary density.
"""
import json
from math import comb, factorial
from collections import Counter
from sporadic_cm_detector import FAMILIES, primes_between
from ordinary_detector_and_mechanism import scan
from defect_valuation_probe import val

N = 520


def hulek_verrill():
    domb = FAMILIES['Domb']()               # Domb(n) = sum over 4-part compositions of multinomial^2
    return [sum(comb(n, a) ** 2 * domb[n - a] for a in range(n + 1)) for n in range(N + 1)]


def a_3k_squared():
    return [(factorial(3 * k) // factorial(k) ** 3) ** 2 for k in range(N + 1)]


def c2k_fourth():
    return [comb(2 * k, k) ** 4 for k in range(N + 1)]


def quintic():
    return [factorial(5 * k) // factorial(k) ** 5 for k in range(N + 1)]


FAMS = {'Hulek-Verrill (AESZ 34)': (hulek_verrill, [-7], [1, 9, 25]),
        '((3k)!/k!^3)^2 = 4F3(1/3,1/3,2/3,2/3)': (a_3k_squared, [-5832], [729]),
        'C(2k,k)^4 (control)': (c2k_fourth, [], [256]),
        '(5k)!/k!^5 quintic (control)': (quintic, [], [3125])}


def main():
    out = {}
    for name, (gen, attractors, singular) in FAMS.items():
        u = gen()
        assert u[:3] == [1, u[1], u[2]]
        surv = scan(u, 10000, False)
        surv_ord = scan(u, 10000, True)
        vals = {}
        for m in sorted(set(surv) | set(attractors)) + [5]:
            c = {'ord': Counter(), 'nonord': Counter()}
            for p in primes_between(7, 199):
                if m % p == 0 or u[1] % p == 0:
                    continue
                typ, v = val(u, m, p)
                c['ord' if typ == 'ord' else 'nonord'][v] += 1
            vals[m] = {k: dict(sorted(v.items())) for k, v in c.items()}
        p9a = all(a in surv for a in attractors)
        p9b = all(s in set(attractors) | set(singular) for s in surv)
        out[name] = {'u_head': u[:6], 'survivors': surv, 'survivors_ord': surv_ord, 'P9a': p9a, 'P9b': p9b,
                     'valuations': {str(k): v for k, v in vals.items()}}
        print(f"{name}: u = {u[:6]}; (B1') survivors {surv}; (B1'_ord) {surv_ord}; P9a {'PASS' if p9a else 'FAIL'}; "
              f"P9b {'PASS' if p9b else 'FAIL'}", flush=True)
        for m, v in vals.items():
            print(f"      m = {m}: valuation histograms {v}", flush=True)
    json.dump(out, open('weight3_attractor_scan.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
