"""Camion (strong-tournament rotation) and tournament-count tests (THM-4746) on products of 19, 79, 101, 169.

Korselt form (THM-4746 Thm 4, with the mod p^2 law (7) for the factor 13^2): s(n) = c(n) (mod n) iff for each odd
prime p || n, p | D(r_p) with r_p = n mod (p-1) in 1..p-1; for p^2 = 169 || n, 169 | D(r) with r = n mod p(p-1)
taken >= 3 (period p(p-1) from m >= 3).  D = s - c, s = strong labelled tournaments (Wright), c = 2^((m-1)/2).
The tournament-count test: T(n) = c(n) (mod n) iff ord_n(2) | (n-1)^2/2."""
from math import comb
from sympy import n_order, factorint

def strong_mod(rmax, mod):
    T = [pow(2, comb(m, 2), mod) for m in range(rmax + 1)]
    s = [0] * (rmax + 1)
    for m in range(1, rmax + 1):
        acc = T[m]
        for k in range(1, m):
            acc -= comb(m, k) % mod * s[k] % mod * T[m - k]
        s[m] = acc % mod
    return s

def D_mod(r, mod):
    s = strong_mod(r, mod)
    c = pow(2, (r - 1) // 2, mod) if r % 2 else 0
    return (s[r] - c) % mod

cases = {'19*79': 19 * 79, '19*101': 19 * 101, '79*101': 79 * 101, '19*79*101': 19 * 79 * 101,
         '19*169': 19 * 169, '79*169': 79 * 169, '101*169': 101 * 169, '19*79*169 (cycle product)': 19 * 79 * 169,
         '19*79*101*169': 19 * 79 * 101 * 169, '169': 169, '13*101': 13 * 101}
for name, n in cases.items():
    f = factorint(n)
    verdict = True
    detail = []
    for p, e in f.items():
        if e == 1:
            r = n % (p - 1) or (p - 1)
            ok = D_mod(r, p) == 0
            detail.append('p=%d r=%d %s' % (p, r, 'zero' if ok else 'nonzero'))
        else:
            assert p == 13 and e == 2
            per = p * (p - 1)
            r = n % per
            while r < 3:
                r += per
            ok = D_mod(r, p * p) == 0
            detail.append('169: r=%d %s' % (r, 'zero' if ok else 'nonzero'))
        verdict &= ok
    count_test = ((n - 1) ** 2 // 2) % n_order(2, n) == 0
    print('%-26s n=%-9d strong Camion test: %-5s [%s]; tournament-count test: %s' % (name, n, verdict, '; '.join(detail), count_test))
