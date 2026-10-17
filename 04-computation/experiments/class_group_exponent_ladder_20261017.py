#!/usr/bin/env python3
"""The exponent ladder of imaginary quadratic class groups: for each e, every discriminant D with -X <= D < 0 whose
form class group Cl(D) has exponent dividing e (S40).

e = 1 is the class number one problem (largest |D| = 163, Heegner-Baker-Stark), e = 2 is "one class per genus"
(the idoneal numbers are the n with D = -4n; largest 1848, complete under GRH and up to at most one or two unknown
further values unconditionally, Weinberger 1973). Larger e give the analogous finite lists (conditionally complete; see
the note). Method: (1) count reduced primitive forms of every discriminant down to -X in one pass over (a, b, c);
(2) genus count 2^(mu - 1) from the factorisation (Gauss; Cox, Primes of the form x^2 + ny^2, Thm 3.15 / Prop 3.11);
(3) for every D whose class number has only prime factors dividing e (and odd class number if e is odd), list the
reduced forms and compute the exponent of Cl(D) exactly by Dirichlet composition (three-term gcd form, Buchmann-Vollmer
section 5.2) and reduction. Validation: the e = 2 list restricted to D = -4n reproduces Euler's 65 idoneal numbers.
Usage: python class_group_exponent_ladder_20261017.py X
Session opus-2026-10-17-S40.
"""
import sys, os, json, math
import numpy as np
from numba import njit


@njit(cache=True)
def gcd(a, b):
    a = abs(a); b = abs(b)
    while b:
        a, b = b, a % b
    return a


@njit(cache=True)
def count_forms(X):
    """h[n] = number of reduced primitive forms of discriminant -n (n = 0, 3 mod 4), 0 < n <= X."""
    h = np.zeros(X + 1, dtype=np.int32)
    a = 1
    while 3 * a * a <= X:
        for b in range(-a + 1, a + 1):
            c = a
            if b < 0:
                c = a + 1          # a = c needs b >= 0
            while True:
                n = 4 * a * c - b * b
                if n > X:
                    break
                if gcd(gcd(a, b), c) == 1:
                    h[n] += 1
                c += 1
        a += 1
    return h


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def genus_mu(n):
    """mu for the discriminant D = -n (n = 0 or 3 mod 4): number of genera is 2^(mu - 1)."""
    if n % 4 == 3:
        return len([p for p in factor(n) if p != 2])
    m = n // 4
    r = len([p for p in factor(m) if p != 2]) if m > 1 else 0
    if m % 4 == 3:
        return r
    if m % 4 in (1, 2):
        return r + 1
    if m % 8 == 4:
        return r + 1
    return r + 2                                  # m = 0 mod 8


@njit(cache=True)
def reduced_forms_nb(D, hD):
    out = np.zeros((hD, 3), dtype=np.int64)
    k = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (b < 0 and a == c):
                continue
            if gcd(gcd(a, b), c) != 1:
                continue
            out[k, 0] = a; out[k, 1] = b; out[k, 2] = c
            k += 1
        a += 1
    return out


def reduced_forms(D, hD):
    return [tuple(int(x) for x in r) for r in reduced_forms_nb(D, hD)]


def power(f, k, D):
    b0 = D % 2
    r = (1, b0, (b0 * b0 - D) // 4)
    while k:
        if k & 1:
            r = compose(r, f, D)
        f = compose(f, f, D)
        k >>= 1
    return r


def exponent_divides(D, forms, e, rng):
    """Exact test that f^e is the identity for every reduced form; random forms first (a proper subgroup Cl[e]
    has index >= 2, so most failures are found at once)."""
    b0 = D % 2
    ident = (1, b0, (b0 * b0 - D) // 4)
    idx = rng.permutation(len(forms))
    for i in idx:
        if power(forms[i], e, D) != ident:
            return False
    return True


def reduce_form(a, b, c):
    while True:
        if c < a:
            a, b, c = c, -b, a
            continue
        if b > a or b <= -a:
            k = (a - b) // (2 * a)                 # b + 2ak in (-a, a]
            c = c + k * b + k * k * a
            b = b + 2 * a * k
            continue
        if a == c and b < 0:
            b = -b
        return (a, b, c)


def xgcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def compose(f, g, D):
    """Dirichlet composition of primitive forms of discriminant D (Buchmann-Vollmer, three-term extended gcd)."""
    a1, b1, c1 = f
    a2, b2, c2 = g
    s = (b1 + b2) // 2
    d0, u0, v0 = xgcd(a1, a2)                     # u0 a1 + v0 a2 = d0
    e, x, w = xgcd(d0, s)                         # x d0 + w s = e
    u, v = x * u0, x * v0                         # u a1 + v a2 + w s = e
    a3 = a1 * a2 // (e * e)
    b3 = (u * a1 * b2 + v * a2 * b1 + w * (b1 * b2 + D) // 2) // e
    b3 %= 2 * a3
    c3 = (b3 * b3 - D) // (4 * a3)
    return reduce_form(a3, b3, c3)


def exponent(D, forms):
    b0 = D % 2
    ident = (1, b0, (b0 * b0 - D) // 4)
    L = 1
    for f in forms:
        g, k = f, 1
        while g != ident:
            g = compose(g, f, D)
            k += 1
            if k > len(forms):
                raise RuntimeError('order exceeds class number at D = %d' % D)
        L = L * k // math.gcd(L, k)
    return L


def is_fundamental(n):
    if n % 4 == 3:
        return all(e == 1 for e in factor(n).values())
    if n % 4 == 0:
        m = n // 4
        return m % 4 in (1, 2) and all(e == 1 for e in factor(m).values())
    return False


def main():
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 6
    h = count_forms(X)
    EXPS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16]
    lists = {e: [] for e in EXPS}
    rng = np.random.default_rng(20261017)
    idoneal = []
    ncand = 0
    for n in range(3, X + 1):
        if n % 4 not in (0, 3):
            continue
        hn = int(h[n])
        if hn == 0:
            continue
        mu = genus_mu(n)
        if hn == 2 ** (mu - 1) and n % 4 == 0:
            idoneal.append(n // 4)
        pf = set(factor(hn)) if hn > 1 else set()
        need = [e for e in EXPS if pf <= set(factor(e)) and (e % 2 == 0 or mu == 1)]
        if not need:
            continue
        if hn == 2 ** (mu - 1):
            ok = [e for e in need if e % (1 if hn == 1 else 2) == 0]
        else:
            ncand += 1
            forms = reduced_forms(-n, hn)
            ok = []
            for e in sorted(need):
                if any(e % d == 0 for d in ok) or exponent_divides(-n, forms, e, rng):
                    ok.append(e)
            if n <= 20000 and ok:
                ex = exponent(-n, forms)
                assert all(e % ex == 0 for e in ok) and all(e % ex for e in need if e not in ok), n
        for e in ok:
            lists[e].append((n, hn, is_fundamental(n)))
    print('X = %d; exponent computed by composition for %d candidate discriminants' % (X, ncand))
    print('idoneal numbers recovered (D = -4n with one class per genus): %d, largest %d' % (len(idoneal), max(idoneal)))
    known65 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25, 28, 30, 33, 37, 40, 42, 45, 48, 57, 58,
               60, 70, 72, 78, 85, 88, 93, 102, 105, 112, 120, 130, 133, 165, 168, 177, 190, 210, 232, 240, 253, 273, 280,
               312, 330, 345, 357, 385, 408, 462, 520, 760, 840, 1320, 1365, 1848]
    print('equals Euler\'s list of 65:', sorted(idoneal) == known65)
    hist = {}
    for m in idoneal:
        hist[int(h[4 * m])] = hist.get(int(h[4 * m]), 0) + 1
    print('class numbers of the idoneal discriminants -4n:', dict(sorted(hist.items())))
    summary = {}
    print('  e   #D   #fund   largest |D|   largest fundamental |D|   class group order at the largest')
    for e in EXPS:
        L = lists[e]
        fund = [t for t in L if t[2]]
        big = max(L) if L else None
        bigf = max(fund) if fund else None
        summary[e] = dict(count=len(L), fundamental=len(fund), largest=big, largest_fundamental=bigf,
                          top=sorted(L)[-6:])
        print('%3d %5d %6d   %10s   %10s   %s' % (e, len(L), len(fund), big[0] if big else '-', bigf[0] if bigf else '-',
                                               (big[1], bigf[1] if bigf else None)))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(X=X, idoneal=sorted(idoneal), idoneal_class_numbers=hist, ladder=summary),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'class_group_exponent_ladder_20261017.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
