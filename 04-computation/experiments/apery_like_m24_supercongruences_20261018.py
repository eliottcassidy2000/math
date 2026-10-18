#!/usr/bin/env python3
"""Apery-like sequences against the 26 Mason eta-products of M24: a supercongruence battery (S43).

For each sequence u(n) and each eta-product eta_g = sum a(n) q^n (all 26 M24 cycle shapes, computed as in
m24_paley_golay_eta_20261018.py), and each quadratic twist by the Kronecker symbol (d/p), d in D, we test
    u((p-1)/2)  ==  (d/p) a(p)   (mod p^2)        and the weaker congruence mod p,
for all primes 5 <= p <= PMAX not dividing 6 * level(g) * d. A pair is reported as a mod-p^2 match if every tested
prime satisfies the congruence. Sequences: zeta(3)-type (Apery A_n, Domb, Almkvist-Zudilin, Franel order 4, Cooper's
s_7, s_10, s_18, Yang-Zudilin) and zeta(2)-type (Apery b_n, Franel, Zagier's sporadic C, D, E, F, and the Legendre /
hypergeometric controls binom(2n,n)^2, binom(2n,n)^3).
Usage: python apery_like_m24_supercongruences_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, os, json, math
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m24_paley_golay_eta_20261018 import SHAPES, eta_product


def kron(d, p):
    """Kronecker symbol (d/p) for an odd prime p."""
    d %= p
    if d == 0:
        return 0
    return 1 if pow(d, (p - 1) // 2, p) == 1 else -1


SEQ = {
    # zeta(3)-type (second-order recurrences of degree 3)
    'Apery A (zeta3)': lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1)),
    'Domb': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(n + 1)),
    'Almkvist-Zudilin': lambda n: sum((-1) ** k * 3 ** (n - 3 * k) * comb(n, 3 * k) * comb(n + k, n)
                                      * math.factorial(3 * k) // math.factorial(k) ** 3 for k in range(n // 3 + 1)),
    'Franel4 (=s10)': lambda n: sum(comb(n, k) ** 4 for k in range(n + 1)),
    'Cooper s7': lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n) for k in range(n + 1)),
    'Cooper s18': lambda n: sum((-1) ** k * comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k)
                                * (comb(2 * k, n) + comb(2 * n - 2 * k, n)) for k in range(n + 1)),
    'Yang-Zudilin': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) ** 2 for k in range(n + 1)),
    # zeta(2)-type (second-order recurrences of degree 2)
    'Apery b (zeta2)': lambda n: sum(comb(n, k) ** 2 * comb(n + k, k) for k in range(n + 1)),
    'Franel': lambda n: sum(comb(n, k) ** 3 for k in range(n + 1)),
    'Zagier C': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) for k in range(n + 1)),
    'Zagier D': lambda n: sum(comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(n + 1)),
    'Zagier E': lambda n: sum(comb(n, 2 * k) * comb(2 * k, k) ** 2 * 4 ** (n - 2 * k) for k in range(n // 2 + 1)),
    'Zagier F': lambda n: sum((-1) ** k * 8 ** (n - k) * comb(n, k) * sum(comb(k, j) ** 3 for j in range(k + 1))
                              for k in range(n + 1)),
    'Zagier A (Franel-type #A)': lambda n: sum(comb(n, k) ** 3 for k in range(n + 1)),
    'Zagier B': lambda n: sum((-1) ** k * 3 ** (n - 3 * k) * comb(n, 3 * k) * math.factorial(3 * k) // math.factorial(k) ** 3
                              for k in range(n // 3 + 1)),
    'binom(2n,n)^2 (control)': lambda n: comb(2 * n, n) ** 2,
    'binom(2n,n)^3 (control)': lambda n: comb(2 * n, n) ** 3,
}

D = [1, -1, 2, -2, 3, -3, 5, -5, 6, -6, 7, -7]


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    N = PMAX + 2
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    forms = {}
    for name, (sh, cen) in SHAPES.items():
        b = eta_product(sh, N)
        a = [0] * (N + 1)
        for n in range(N):
            a[n + 1] = b[n]
        level = math.lcm(*sh.keys()) * min(sh.keys())
        forms[name] = (a, level, sum(sh.values()) // 2)
    vals = {}
    for sname, fn in SEQ.items():
        vals[sname] = {p: fn((p - 1) // 2) for p in primes}
    results = []
    print('mod p^2 matches  u((p-1)/2) == (d/p) a_g(p), all primes 5 <= p <= %d with p not dividing 6*level*d:' % PMAX)
    for sname in SEQ:
        found = []
        for gname, (a, level, wt) in forms.items():
            for d in D:
                tested = [p for p in primes if level % p and 6 % p and d % p]
                if len(tested) < 10:
                    continue
                ok2 = all((vals[sname][p] - kron(d, p) * a[p]) % (p * p) == 0 for p in tested)
                ok1 = all((vals[sname][p] - kron(d, p) * a[p]) % p == 0 for p in tested)
                if ok2 or ok1:
                    found.append(dict(seq=sname, cls=gname, weight=wt, level=level, d=d, mod='p^2' if ok2 else 'p',
                                      primes=len(tested)))
        results += found
        best = [f for f in found if f['mod'] == 'p^2']
        weak = [f for f in found if f['mod'] == 'p']
        print('  %-26s p^2: %s' % (sname, ', '.join('%s(w%d,N%d)*(%d/p)' % (f['cls'], f['weight'], f['level'], f['d']) for f in best) or '-'))
        if weak:
            print('  %-26s  p only: %s' % ('', ', '.join('%s*(%d/p)' % (f['cls'], f['d']) for f in weak)))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(PMAX=PMAX, results=results), open(os.path.join(here, '..', '..', '05-knowledge', 'results',
              'apery_like_m24_supercongruences_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
