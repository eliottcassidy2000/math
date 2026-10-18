#!/usr/bin/env python3
"""Closed forms in N for the Witt-forced degree-N generators G_5, G_7, G_9 (S44).

For each N in a range, compute the unique admissible lowest-weight vector G_L(N) (witt_closure_forces_generators:
generators at Li 3, 5, 7, 9; mod six primes, CRT, rational reconstruction), in Cayley coordinates.  Then reconstruct
each coefficient as a rational function of N: the smallest (deg num, deg den) fitting the first points, accepted only
if it reproduces every remaining point (the number of confirming points is reported).  dim S_L is recorded for every N.
Usage: python witt_rigidity_interpolation_20261018.py [Nmin Nmax]
"""
import os, sys, json
from fractions import Fraction
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
from witt_closure_forces_generators_20261018 import closure, cayley_of
from semiinvariant_generators_20261018 import crt, ratrec_big

PRIMES = (2147483647, 2147483629, 2147483587, 2147483579, 2147483563, 2147483549)


def forced_generators(N, Lmax=9):
    per = [closure(N, Lmax, P, gen_pattern={3, 5, 7, 9})[0] for P in PRIMES]
    dims = {L: len(per[0][L]) for L in per[0]}
    out = {}
    for L in (5, 7, 9):
        if L <= Lmax and dims.get(L) == 1 and N >= L:
            co = [cayley_of(s[L][0], L, P) for s, P in zip(per, PRIMES)]
            out[L] = {nu: ratrec_big(*crt([c[nu] for c in co], PRIMES)) for nu in co[0]}
    return dims, out


def fit_rational(points, maxdeg=8):
    Nn = sp.Symbol('N')
    pts = [(sp.Integer(n), sp.Rational(v.numerator, v.denominator)) for n, v in points]
    for tot in range(0, 2 * maxdeg + 1):
        for dn in range(0, tot + 1):
            d = tot - dn
            nunk = (d + 1) + dn
            if nunk + 2 > len(pts):
                continue
            a = sp.symbols(f'a0:{d+1}')
            b = sp.symbols(f'b0:{dn}') if dn else ()
            eqs = []
            for n, v in pts[:nunk]:
                num = sum(a[i] * n ** i for i in range(d + 1))
                den = n ** dn + sum(b[i] * n ** i for i in range(dn))
                eqs.append(num - v * den)
            sol = sp.solve(eqs, list(a) + list(b), dict=True)
            if not sol:
                continue
            s = sol[0]
            if any(sym not in s for sym in list(a) + list(b)):
                continue
            num = sum(s[a[i]] * Nn ** i for i in range(d + 1))
            den = Nn ** dn + sum(s[b[i]] * Nn ** i for i in range(dn))
            expr = sp.factor(num / den)
            if all(sp.simplify(expr.subs(Nn, n) - v) == 0 for n, v in pts):
                return expr, len(pts) - nunk
    return None, 0


if __name__ == '__main__':
    Nmin = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    Nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    data = {}
    dimsall = {}
    for N in range(Nmin, Nmax + 1):
        dims, g = forced_generators(N)
        dimsall[N] = dims
        data[N] = g
        print(f'N={N}: dim S_L {dims}', flush=True)
    res = {'dims': {N: d for N, d in dimsall.items()}, 'fits': {}}
    for L in (5, 7, 9):
        nus = sorted(set().union(*[set(data[N].get(L, {})) for N in data]), key=lambda t: (-t[0], t))
        for nu in nus:
            if nu == (L,):
                continue
            pts = [(N, data[N][L][nu]) for N in sorted(data) if L in data[N] and data[N][L].get(nu) is not None]
            expr, conf = fit_rational(pts)
            res['fits'][f'G_{L} {nu}'] = {'expr': str(expr), 'confirming_points': conf, 'n_points': len(pts)}
            print(f'G_{L} coefficient of A{list(nu)}: {expr}   ({conf} confirming points of {len(pts)})', flush=True)
    json.dump(res, open(os.path.join(HERE, f'witt_rigidity_interpolation_{Nmin}_{Nmax}.json'), 'w'), indent=1, default=str)
