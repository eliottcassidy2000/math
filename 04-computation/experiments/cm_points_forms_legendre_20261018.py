#!/usr/bin/env python3
"""CM points of sum C(2k,k)^3/m^k: governing newforms, Legendre squares, and the general twin (S43).

(1) Exact identity check: P_N(x)^2 = sum_k C(N+k,2k) C(2k,k)^2 ((x^2-1)/4)^k for N <= 25 (polynomial identity in x).
    At N = n = (p-1)/2, C(n+k,2k) == (-1)^k a_k / 4^k (mod p^2), so sum_{k<=n} C(2k,k)^3 m^-k == P_n(x)^2 (mod p^2) with
    x^2 = 1 - 64/m  (checked for the CM m below, p <= PMAX (the committed run used PMAX = 200; corrected after S43b m6), x^2 taken mod p^2).
(2) For each of Sun's CM values m, find the eta-quotient f, character eps and scale s with
    sum_{k<=n} C(2k,k)^3 m^-k == eps(p) s(p) a_f(p) (mod p^2), 5 <= p <= PMAX (committed run: 200), p not dividing m
    (library: 4C = eta(4t)^6 [Q(i)], 8A = eta(t)^2 eta(2t) eta(4t) eta(8t)^2 [Q(sqrt-2)], 7AB = eta(t)^3 eta(7t)^3
    [Q(sqrt-7)], 12 = eta(2t)^3 eta(6t)^3 [Q(sqrt-3), not an M24 cycle shape]).
(3) The general twin: T_m = sum_j C(2j,j) C(p-1-j,j)^2 m^-j  ==  (1 - p q_p(m)/3) F_m  (mod p^2)  <=>  (L_p) at m.
Usage: python cm_points_forms_legendre_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json
from math import comb
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from halfperiod_twin_class_forms_20261018 import eta_quotient, legendre


def legendreP_coeffs(N):
    """coefficients of P_N(x) as Fractions, via the three-term recurrence."""
    P0 = [Fraction(1)]
    if N == 0:
        return P0
    P1 = [Fraction(0), Fraction(1)]
    for n in range(1, N):
        # (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
        xP = [Fraction(0)] + P1
        Pn1 = [Fraction(0)] * (n + 2)
        for i, c in enumerate(xP):
            Pn1[i] += Fraction(2 * n + 1, n + 1) * c
        for i, c in enumerate(P0):
            Pn1[i] -= Fraction(n, n + 1) * c
        P0, P1 = P1, Pn1
    return P1


def polymul(a, b):
    c = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    out = {}
    # (1) identity
    ok = True
    for N in range(0, 26):
        P = legendreP_coeffs(N)
        lhs = polymul(P, P)
        # rhs in powers of x: ((x^2-1)/4)^k
        rhs = [Fraction(0)] * (2 * N + 1)
        base = [Fraction(-1, 4), Fraction(0), Fraction(1, 4)]
        pw = [Fraction(1)]
        for k in range(N + 1):
            c = comb(N + k, 2 * k) * comb(2 * k, k) ** 2
            for i, v in enumerate(pw):
                rhs[i] += c * v
            pw = polymul(pw, base)
        lhs += [Fraction(0)] * (len(rhs) - len(lhs))
        if lhs != rhs[:len(lhs)] or any(rhs[len(lhs):]):
            ok = False
    print('(1) P_N(x)^2 = sum_k C(N+k,2k) C(2k,k)^2 ((x^2-1)/4)^k for 0 <= N <= 25:', ok)
    out['legendre_square_identity'] = ok
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    CM = [1, -8, 16, -64, 64, 256, -512, 4096]
    forms = {'4C': {4: 6}, '8A': {1: 2, 2: 1, 4: 1, 8: 2}, '7AB': {1: 3, 7: 3}, '12(-3)': {2: 3, 6: 3}}
    coef = {f: eta_quotient(s, PMAX) for f, s in forms.items()}
    chars = {'1': lambda p: 1, '(-1/p)': lambda p: legendre(-1, p), '(2/p)': lambda p: legendre(2, p),
             '(-2/p)': lambda p: legendre(-2, p), '(3/p)': lambda p: legendre(3, p), '(-3/p)': lambda p: legendre(-3, p),
             '(7/p)': lambda p: legendre(7, p), '(-7/p)': lambda p: legendre(-7, p)}
    legok = True
    twin_ok = {}
    ident = {}
    for m in CM:
        found = []
        vals = {}
        Tvals = {}
        for p in primes:
            if m % p == 0:
                continue
            n = (p - 1) // 2
            pp = p * p
            im = pow(m % pp, -1, pp)
            F = sum(comb(2 * k, k) ** 3 * pow(im, k, pp) for k in range(n + 1)) % pp
            T = sum(comb(2 * j, j) * comb(p - 1 - j, j) ** 2 * pow(im, j, pp) for j in range(n + 1)) % pp
            vals[p] = F
            Tvals[p] = T
            # Legendre square: P_n(x)^2 is a polynomial in y = x^2; evaluate at y = 1 - 64/m mod p^2
            y = (1 - 64 * im) % pp
            Pn = legendreP_coeffs(n)
            Pn2 = polymul(Pn, Pn)
            s = 0
            for i in range(0, len(Pn2), 2):
                c = Pn2[i]
                if c:
                    s += c.numerator * pow(c.denominator, -1, pp) * pow(y, i // 2, pp)
            if (s - F) % pp:
                legok = False
        for fname in forms:
            for cname, ch in chars.items():
                for i in range(-8, 9):
                    good = True
                    for p in vals:
                        n = (p - 1) // 2
                        pp = p * p
                        sc = pow(2, i * n, pp) if i >= 0 else pow(pow(2, -i * n, pp), -1, pp)
                        if (vals[p] - ch(p) * sc * coef[fname][p - 1]) % pp:
                            good = False
                            break
                    if good:
                        found.append('%s * 2^(%d n) * a_%s(p)' % (cname, i, fname))
        # twin
        tok = True
        for p in vals:
            pp = p * p
            q = (pow(m % pp, p - 1, pp) - 1) // p % p
            rhs = vals[p] * (1 - p * q * pow(3, -1, pp)) % pp
            if (Tvals[p] - rhs) % pp:
                tok = False
        twin_ok[str(m)] = tok
        ident[str(m)] = found
        jinv = Fraction((m - 16) ** 3, m)
        print('(2) m = %5d, j = %s: sum C(2k,k)^3 m^-k == %s (mod p^2), p <= %d' % (m, jinv, found if found else 'NO MATCH', PMAX))
        print('    (3) twin T_m == (1 - p q_p(m)/3) F_m (mod p^2):', tok)
    print('(1b) sum_{k<=n} C(2k,k)^3 m^-k == P_n(x)^2, x^2 = 1 - 64/m (mod p^2), all CM m, p <= %d: %s' % (PMAX, legok))
    out['legendre_mod_p2'] = legok
    out['forms'] = ident
    out['twin'] = twin_ok
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_points_forms_legendre_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
