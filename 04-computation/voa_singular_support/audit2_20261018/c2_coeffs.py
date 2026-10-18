#!/usr/bin/env python3
"""Li-degree-0 symbols (C_2 relations) of W_3(p,q) at given weights, normalised, compared with given rationals at two
primes (robust: compares mod P instead of rational reconstruction).  usage: python c2_coeffs.py p q w1 [w2 ...]"""
import sys
from fractions import Fraction
from aud_algebras import W3, LW, WW
from vamod import rowspace_basis
p, q = int(sys.argv[1]), int(sys.argv[2]); ws = [int(a) for a in sys.argv[3:]]
c = 2 * (1 - Fraction(12 * (p - q) ** 2, p * q))
for P in (2147483647, 2147483629):
    va = W3(c, P)
    G = va.gram_all(2 * max(ws))
    I = va.radical(G)
    for w in ws:
        B = va.basis(2 * w)
        order = sorted(B, key=lambda m: (va.mono_li(m), m))
        basisI = rowspace_basis(I[(2 * w, 0)], order, P)
        li0 = [v for v in basisI if any(va.mono_li(m) == 0 for m in v)]
        for v in li0:
            sym = {m: cc for m, cc in v.items() if va.mono_li(m) == 0}
            terms = {}
            for m, cc in sym.items():
                nx = sum(1 for X in m if X[0] == LW); ny = sum(1 for X in m if X[0] == WW)
                terms[(nx, ny)] = cc
            print(f'P={P} weight {w}: raw Li-0 symbol (x^a y^b: coeff mod P) =', terms)
            # normalise to the y-highest monomial
            yb = max(terms, key=lambda t: t[1])
            inv = pow(terms[yb], P - 2, P)
            print('      normalised to', yb, ':', {t: cc * inv % P for t, cc in terms.items()})
    # author's q=13 rationals mod P for comparison
    def F(x): x = Fraction(x); return x.numerator % P * pow(x.denominator % P, P - 2, P) % P
    if (p, q) == (3, 13):
        print('   author y^4 - (13/56)x^3y^2 - (28529/3278)x^6 mod P: coeff x^3y^2 =', F(Fraction(-13, 56)), ' x^6 =', F(Fraction(-28529, 3278)))
        print('   author xy^3 - (13/448)x^4y: x^4y =', F(Fraction(-13, 448)))
