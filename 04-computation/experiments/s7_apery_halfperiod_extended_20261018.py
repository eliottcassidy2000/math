#!/usr/bin/env python3
"""THM-4669 controls: s_7(n) == (-1)^n A(n) (mod p^2) at n = (p^r - 1)/2, with mod-p^3 status (S43; rewritten after the
S43 audit, finding m10, so that every line of the committed .out is produced by this script).

(1) r = 1: all primes 5 <= p < 1500 (mod p^2), and the mod-p^3 status for 5 <= p < 250.
(2) Frobenius points r >= 2 for p in {3, 5, 7, 11, 13, 17, 19, 23, 29, 31} with n <= 1300: mod p^2 and mod p^3.
    At the non-ordinary prime p = 11 (a(11) = -44 for eta(2t)^4 eta(4t)^4), v_11(A) = v_11(s_7) = r, so the congruence
    holds to higher order there (audit m9).
(3) The control points n = p - 1 and n = (3p - 1)/2 for 5 <= p <= 300 (mod p and mod p^2).
Usage: python s7_apery_halfperiod_extended_20261018.py OUT.json
Session opus-2026-10-17-S43.
"""
import sys, json
from math import comb


def A(n):
    return sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))


def s7(n):
    return sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range((n + 1) // 2, n + 1))


def vp(x, p):
    if x == 0:
        return 99
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def main():
    res = {}
    primes = [p for p in range(5, 1500) if all(p % d for d in range(2, int(p ** .5) + 1))]
    fails = []
    p3 = []
    for p in primes:
        n = (p - 1) // 2
        d = s7(n) - (-1) ** n * A(n)
        if d % (p * p):
            fails.append(p)
        if p < 250 and d % (p ** 3) == 0:
            p3.append(p)
    print('(1) r = 1: primes tested %d (5..%d), mod p^2 failures %s; primes p < 250 where it holds mod p^3: %s'
          % (len(primes), primes[-1], fails, p3), flush=True)
    res['half'] = dict(tested=len(primes), pmax=primes[-1], fails=fails, mod_p3_holds_below_250=p3)
    hi = []
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        for r in (2, 3, 4):
            n = (p ** r - 1) // 2
            if n > 1300:
                continue
            d = s7(n) - (-1) ** n * A(n)
            ok2 = d % (p * p) == 0
            ok3 = d % (p ** 3) == 0
            hi.append(dict(p=p, r=r, n=n, mod_p2=ok2, mod_p3=ok3, v_p_diff=vp(d, p), v_p_A=vp(A(n), p)))
            print('(2) p=%d r=%d n=%d: mod p^2 %s, mod p^3 %s (v_p(diff) = %d, v_p(A) = %d)'
                  % (p, r, n, ok2, ok3, vp(d, p), vp(A(n), p)), flush=True)
    res['higher'] = hi
    ctl = []
    for p in [q for q in range(5, 301) if all(q % d for d in range(2, int(q ** .5) + 1))]:
        for name, n in (('p-1', p - 1), ('(3p-1)/2', (3 * p - 1) // 2)):
            d = s7(n) - (-1) ** n * A(n)
            ctl.append(dict(p=p, point=name, mod_p=d % p == 0, mod_p2=d % (p * p) == 0))
    for name in ('p-1', '(3p-1)/2'):
        m1 = [c['p'] for c in ctl if c['point'] == name and c['mod_p']]
        m2 = [c['p'] for c in ctl if c['point'] == name and c['mod_p2']]
        print('(3) n = %s, 5 <= p <= 300: holds mod p at %s; mod p^2 at %s' % (name, m1, m2), flush=True)
    res['controls'] = ctl
    json.dump(res, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
