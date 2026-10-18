#!/usr/bin/env python3
"""S43b referee: modular identifications (written from scratch).

(1) eta-products 4C = eta(4t)^6, 8A = eta(t)^2 eta(2t) eta(4t) eta(8t)^2, 7AB = eta(t)^3 eta(7t)^3,
    E12 = eta(2t)^3 eta(6t)^3, 4A = eta(2t)^4 eta(4t)^4: CM check a(p) = 4x^2 - 2p (p = x^2 + d y^2, split) / 0 (inert).
(2) F_m = sum_{k<=n} C(2k,k)^3 m^-k mod p^2 at the 8 CM points against chi(p) a_f(p), 5 <= p <= PMAX.
(3) census class values mod p^2 (F, 4^n F, 16^n F, 7AB, 4A) and the mod-p-only classes 6, 7.
(4) j = (m-16)^3/m, CM discriminants; the M24 cycle-shape list and 2^3 6^3.
Usage: python audit_forms.py [PMAX]
"""
import sys
from math import comb, isqrt


def eta_prod(shape, N):
    """coefficients of q^{sum a e/24} * prod_a prod_n (1-q^{an})^{e_a}, shifted so index 0 = leading power."""
    b = [0] * (N + 1)
    b[0] = 1
    for a, e in shape.items():
        for _ in range(e):
            for n in range(1, N // a + 1):
                s = a * n
                for i in range(N, s - 1, -1):
                    b[i] -= b[i - s]
    return b


def primes_in(a, b):
    return [p for p in range(max(a, 2), b + 1) if all(p % d for d in range(2, isqrt(p) + 1))]


def leg(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def rep(p, d):
    """x with p = x^2 + d y^2 (or None)."""
    for y in range(0, isqrt(p // d) + 1):
        r = p - d * y * y
        x = isqrt(r)
        if x * x == r:
            return x, y
    return None


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    N = PMAX + 2
    shapes = {'4C': {4: 6}, '8A': {1: 2, 2: 1, 4: 1, 8: 2}, '7AB': {1: 3, 7: 3}, 'E12': {2: 3, 6: 3}, '4A': {2: 4, 4: 4}}
    coef = {}
    for nm, sh in shapes.items():
        assert sum(a * e for a, e in sh.items()) == 24
        b = eta_prod(sh, N)
        coef[nm] = lambda m, b=b: b[m - 1]          # q^1 * prod -> a(m) = b[m-1]
    primes = primes_in(5, PMAX)
    # (1) CM structure of the weight-3 forms
    cmd = {'4C': 1, '8A': 2, '7AB': 7, 'E12': 3}
    bad1 = []
    for nm, d in cmd.items():
        for p in primes:
            if p == d or (nm == '7AB' and p == 7) or (nm == 'E12' and p == 3):
                continue
            r = rep(p, d)
            a = coef[nm](p)
            if r is None:
                if a != 0:
                    bad1.append((nm, p, 'inert but a != 0'))
            else:
                x, y = r
                if nm == '4C' and x % 2 == 0:
                    x = y
                if abs(a) != abs(4 * x * x - 2 * p):
                    bad1.append((nm, p, a, 4 * x * x - 2 * p))
    print('[C1] CM weight-3 eta-products: |a(p)| = |4x^2 - 2p| at split p, 0 at inert p, 5<=p<=%d: %s'
          % (PMAX, 'all hold' if not bad1 else bad1[:6]))
    # signs: report a(p) vs 4x^2-2p sign agreement
    sg = {}
    for nm, d in cmd.items():
        agree = 0
        tot = 0
        for p in primes:
            if p == 7 and nm == '7AB':
                continue
            r = rep(p, d)
            if r and coef[nm](p):
                x, y = r
                if nm == '4C' and x % 2 == 0:
                    x = y
                tot += 1
                agree += coef[nm](p) == 4 * x * x - 2 * p
        sg[nm] = (agree, tot)
    print('     a(p) == +(4x^2 - 2p) (x odd for 4C): %s' % sg)
    # (2) F_m against chi(p) a_f(p)
    assign = {1: ('7AB', lambda p: 1), 4096: ('7AB', lambda p: leg(-1, p)), -8: ('4C', lambda p: 1),
              64: ('4C', lambda p: 1), -512: ('4C', lambda p: leg(2, p)), -64: ('8A', lambda p: leg(-1, p)),
              16: ('E12', lambda p: 1), 256: ('E12', lambda p: leg(-1, p))}
    for m, (f, ch) in assign.items():
        bad = []
        for p in primes:
            pp = p * p
            n = (p - 1) // 2
            w = pow(m % pp, -1, pp)
            F = 0
            c = 1
            pw = 1
            for k in range(n + 1):
                if k:
                    c = c * (2 * k) * (2 * k - 1) // (k * k)
                F = (F + (c % pp) ** 3 * pw) % pp
                pw = pw * w % pp
            if (F - ch(p) * coef[f](p)) % pp:
                bad.append(p)
        jm = (m - 16) ** 3 / m
        print('[C2] m = %5d (j = %s): F_m == chi(p) a_%s(p) mod p^2 for all 5<=p<=%d: %s'
              % (m, int(jm) if jm == int(jm) else jm, f, PMAX, 'yes' if not bad else 'FAILS at %s' % bad[:6]))
    # (3) census class values mod p^2 / mod p
    bad3 = {}
    for p in primes:
        pp = p * p
        n = (p - 1) // 2
        a4C = coef['4C'](p)
        a7 = coef['7AB'](p)
        a4A = coef['4A'](p)
        sgn = -1 if n % 2 else 1
        vals = {
            'b(n)': sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1)),
            'sum C(n+k,k)^3': sum(comb(n + k, k) ** 3 for k in range(n + 1)),
            '(i-a)': sum(comb(n + k, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(n + 1)),
            '(i-b)': sum(comb(n, k) ** 3 * comb(2 * n, n + k) for k in range(n + 1)),
            '(i-c)': sum(comb(n, k) * comb(n + k, k) ** 2 * comb(2 * n, n + k) for k in range(n + 1)),
            'sum C(2k,k)^3': sum(comb(2 * k, k) ** 3 for k in range(n + 1)),
            'T_1': sum(comb(2 * j, j) * comb(2 * n - j, j) ** 2 for j in range(n + 1)),
            'Apery A': sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1)),
        }
        targ = {'b(n)': a4C, 'sum C(n+k,k)^3': a4C, '(i-a)': pow(4, n, pp) * a4C, '(i-b)': pow(16, n, pp) * a4C,
                '(i-c)': pow(16, n, pp) * a4C, 'sum C(2k,k)^3': a7, 'T_1': a7, 'Apery A': a4A}
        for k_, v in vals.items():
            if (v - targ[k_]) % pp:
                bad3.setdefault(k_, []).append(p)
        u6a = sum(comb(n + k, k) ** 3 * comb(n, 2 * k) for k in range(n // 2 + 1))
        u6b = sum(comb(n + k, k) * comb(2 * k, k) * comb(2 * k, n) * comb(2 * n - 2 * k, n - k) for k in range(n + 1))
        c7a = sum(comb(2 * n - 2 * k, n - k) ** 2 * comb(n + k, 2 * k) * comb(2 * n, n + k) for k in range(n + 1))
        c7b = sum(comb(n + k, 2 * k) ** 3 * comb(2 * n, n + k) for k in range(n + 1))
        if (u6a - leg(-1, p) * a4A) % p or (u6b - a4A) % p:
            bad3.setdefault('class 6 mod p', []).append(p)
        if (u6a - leg(-1, p) * a4A) % pp == 0:
            bad3.setdefault('class 6 holds mod p^2 (unexpected)', []).append(p)
        if (c7a - a7) % p or (c7b - a7) % p:
            bad3.setdefault('class 7 mod p', []).append(p)
        if (c7a - a7) % pp == 0:
            bad3.setdefault('class 7a == a7 mod p^2 (sporadic)', []).append(p)
    print('[C3] census class values (b, sum C(n+k,k)^3 == a_4C; (i-a) == 4^n a_4C; (i-b),(i-c) == 16^n a_4C;'
          ' sum C(2k,k)^3, T_1 == a_7AB; Apery A == a_4A, all mod p^2; classes 6, 7 mod p), 5 <= p <= %d:' % PMAX)
    print('     exceptions: %s' % ({k: v[:8] for k, v in bad3.items()} if bad3 else 'none'))
    # (4) j-invariants and discriminants; M24 shapes
    rational_cm = {0: -3, 1728: -4, -3375: -7, 8000: -8, -32768: -11, 54000: -12, 287496: -16, -884736: -19,
                   -12288000: -27, 16581375: -28, -884736000: -43, -147197952000: -67, -262537412640768000: -163}
    for m in (1, 4096, -8, 64, -512, -64, 16, 256):
        j = (m - 16) ** 3 // m if (m - 16) ** 3 % m == 0 else None
        print('[C4] m = %5d: j = %s, D = %s' % (m, j, rational_cm.get(j)))
    m24 = [{1: 24}, {1: 8, 2: 8}, {2: 12}, {1: 6, 3: 6}, {3: 8}, {2: 4, 4: 4}, {1: 4, 2: 2, 4: 4}, {4: 6}, {1: 4, 5: 4},
           {1: 2, 2: 2, 3: 2, 6: 2}, {6: 4}, {1: 3, 7: 3}, {1: 2, 2: 1, 4: 1, 8: 2}, {2: 2, 10: 2}, {1: 2, 11: 2},
           {2: 1, 4: 1, 6: 1, 12: 1}, {12: 2}, {1: 1, 2: 1, 7: 1, 14: 1}, {1: 1, 3: 1, 5: 1, 15: 1}, {3: 1, 21: 1},
           {1: 1, 23: 1}]
    print('[C4] 21 M24 cycle shapes; 2^3 6^3 in list: %s' % ({2: 3, 6: 3} in m24))


if __name__ == '__main__':
    main()
