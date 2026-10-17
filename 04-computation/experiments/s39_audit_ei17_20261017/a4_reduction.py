#!/usr/bin/env python3
"""Audit of THM-4665, item 4: the residue fields and the action of complex conjugation at the 46 primes.

For each p in the two lists: the order of p in (Z/55)^x, the factorization pattern of Phi_55 mod p (sympy), and the
decomposition group <p mod 55>; whether sigma_(-1) lies in it (conjugation fixes P) or not (conjugation moves P).
Concrete check of the two reduction models with a sample unit direction u = w / conj(w), w = 3 + z - 2 z^7 + z^20
(z = zeta_55), which has norm one in L:
  * p = 1 mod 55: P = (p, z - r) for each of the 40 roots r of Phi_55 in F_p; conj(z) = z^54, so
    (u mod P, conj(u) mod P) = (w(r)/w(1/r), w(1/r)/w(r)): product 1, i.e. a point of the split model, and
    conj(P) = (p, z - 1/r) != P;
  * p = 54 mod 55: P = (p, g(z)) with g an irreducible quadratic factor; in F_p[y]/(g), y^p = y^(-1) (Frobenius =
    sigma_p = sigma_(-1) on the image of z, hence on O_L / P), and u mod P = w(y)^(1-p) has norm u^(p+1) = 1.
Also: complex conjugation on Q(zeta_55) is sigma_(-1): conj(exp(2 pi i/55)) = exp(-2 pi i/55) = zeta^54 (numerical).
"""
import cmath

import sympy as sp

Z = sp.Symbol('z')
PHI = sp.Poly(sp.cyclotomic_poly(55, Z), Z)
S1 = [331, 661, 881, 991, 1321, 1871, 2311, 2531, 2971, 3301, 3631, 3851, 5171, 5281, 5501, 6271, 6491, 7151, 7591, 8581,
      9241, 9461, 9791, 9901, 10891, 11551]
S1_POS = [2861, 3191, 4621, 4951, 7481, 10781]
T54 = [109, 439, 659, 769, 1319, 1429, 1759, 1979, 2089, 2309, 2749, 2969, 3079, 3299, 3739, 4289, 4729, 5059, 5279, 5939]


def order_mod55(p):
    k, x = 1, p % 55
    while x != 1:
        x = x * p % 55; k += 1
    return k


def factor_degrees(p):
    fl = sp.Poly(PHI.as_expr(), Z, modulus=p).factor_list()[1]
    return sorted(f.degree() for f, e in fl), all(e == 1 for f, e in fl), fl


z = cmath.exp(2j * cmath.pi / 55)
print('complex conjugation: |conj(zeta) - zeta^54| = %.2e (so conj = sigma_(-1) on Q(zeta_55))' % abs(z.conjugate() - z ** 54))
print('phi(55) = %d = [L:Q]' % sp.totient(55))
ok_all = True
for p in S1 + S1_POS + T54:
    o = order_mod55(p)
    degs, sqf, fl = factor_degrees(p)
    dec = sorted({pow(p, k, 55) for k in range(o)})
    fixed = (54 in dec)
    want = (o == 1 and degs == [1] * 40 and not fixed) or (o == 2 and degs == [2] * 20 and fixed)
    ok_all &= want and sqf
    print('p = %5d  p mod 55 = %2d  order %d  Phi_55 mod p: %2d factors of degree %d (squarefree %s)  decomposition group %s  '
          'conjugation fixes P: %s' % (p, p % 55, o, len(degs), degs[0], sqf, dec, fixed))
print('all 52 primes behave as stated (split: 40 primes of degree 1, conjugation moves P; 54 mod 55: 20 primes of degree 2,'
      ' conjugation = Frobenius): %s' % ok_all)

# sample unit direction
W = 3 + Z - 2 * Z ** 7 + Z ** 20
for p in (331, 2861):
    roots = [r for r in range(2, p) if PHI.eval(r) % p == 0]
    assert len(roots) == 40
    good = 0
    for r in roots:
        ri = pow(r, -1, p)
        a = int(W.subs(Z, r)) % p; ab = int(W.subs(Z, ri)) % p
        if a == 0 or ab == 0:
            continue
        x = a * pow(ab, -1, p) % p; xb = ab * pow(a, -1, p) % p
        assert x * xb % p == 1 and ri in roots and ri != r
        good += 1
    print('p = %d: 40 roots r of Phi_55; for the %d primes P = (p, z - r) where w is a P-unit, (u, conj u) mod P = '
          '(x, 1/x) and conj(P) = (p, z - 1/r) != P' % (p, good))

for p in (109, 5939):
    degs, sqf, fl = factor_degrees(p)
    g = fl[0][0]
    y = sp.Poly(Z, Z, modulus=p)
    G = sp.Poly(g.as_expr(), Z, modulus=p)
    yp = pow(y, p, G) if False else sp.Poly(sp.rem(sp.Poly(Z ** p, Z, modulus=p).as_expr(), G.as_expr(), modulus=p), Z, modulus=p)
    yinv = sp.Poly(sp.invert(Z, G.as_expr(), modulus=p), Z, modulus=p)
    print('p = %d: g = %s; y^p = %s, y^(-1) = %s, equal: %s' % (p, G.as_expr(), yp.as_expr(), yinv.as_expr(),
                                                                 sp.rem((yp - yinv).as_expr(), G.as_expr(), modulus=p) == 0))
    # u mod P = w(y) / w(y^-1) = w(y) / w(y)^p, norm u^(p+1) = 1
    wy = sp.rem(W, G.as_expr(), modulus=p)
    wyp = sp.rem(sp.Poly(W.subs(Z, Z ** (p % 55)), Z).as_expr(), G.as_expr(), modulus=p)   # w(y^p) = w(y^(p mod 55))
    wbar = sp.rem(sp.Poly(W.subs(Z, Z ** 54), Z).as_expr(), G.as_expr(), modulus=p)          # w(y^-1)
    print('        w(y^p) = w(y^-1) in F_p[y]/(g): %s' % (sp.rem(wyp - wbar, G.as_expr(), modulus=p) == 0))
    u = sp.rem(wy * sp.invert(wbar, G.as_expr(), modulus=p), G.as_expr(), modulus=p)
    # norm: u * u^p = u * u(y^-1)
    up = sp.rem(sp.Poly(sp.expand(u.subs(Z, Z ** 54)), Z).as_expr(), G.as_expr(), modulus=p)
    print('        u mod P has norm u * u^p = %s' % sp.rem(sp.expand(u * up), G.as_expr(), modulus=p))
