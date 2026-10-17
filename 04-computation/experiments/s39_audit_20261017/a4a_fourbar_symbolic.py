#!/usr/bin/env python3
"""S39 audit, item 4 (symbolic part): the four-bar curve, 30a1, and the two Gamma_1(6) models of THM-4663 (8)-(10).

  * quartic invariants I, J, 4I^3 - J^2, j for E'_c : y^2 = (t^2 - 4c)(t - c - 1)(t - c + 3) and for the linkage model
    L_c : w^2 = (t^2 + (1+c)t + 1)^2 - 4ct(t+1)^2; fibre types at c = 0, 1, 9, oo (orders of Delta, I != 0, degrees);
  * the relation between the two families: Phi_2(j(E'_c), j(L_c)) = 0 identically (2-isogeny over Qbar(c)), and the
    abstract isomorphism j(L_c) = j(E'_(sigma(c))), sigma(c) = (c - 9)/(c - 1) (base reparametrisation permuting cusps);
  * the Mordell-Weil torsion: a Q(c)-rational point of order 6 on E'_c (proof of 6 | #E'_c(F_p)) and on L_c;
  * 30a1 = [1,0,1,1,2]: Delta, j; the four-bar quartic t^4 - 6t^3 - 5t^2 - 6t + 1: I, J, disc, j;
  * smoothness of the plane cubic 2u^2 v + 2uv^2 - u^2 - 5uv - v^2 + 2u + 2v (bad primes), the fixed points of iota at p = 5;
  * twist witnesses: for each squarefree D | 30 (both signs) a prime 7 <= p <= 139 with (D/p) = -1 and a_p(30a1) != 0.
Usage: python -B a4a_fourbar_symbolic.py
"""
import sympy as sp
from sympy import symbols, expand, factor, Poly, cancel, together, simplify

c, t, u, v, w = symbols('c t u v w')


def invariants(Q):
    P = Poly(expand(Q), t)
    a4, a3, a2, a1, a0 = [P.coeff_monomial(t ** k) for k in (4, 3, 2, 1, 0)]
    I = expand(12 * a4 * a0 - 3 * a3 * a1 + a2 ** 2)
    J = expand(72 * a4 * a2 * a0 + 9 * a3 * a2 * a1 - 27 * a4 * a1 ** 2 - 27 * a0 * a3 ** 2 - 2 * a2 ** 3)
    D = factor(4 * I ** 3 - J ** 2)
    j = factor(cancel(6912 * I ** 3 / (4 * I ** 3 - J ** 2)))
    return factor(I), J, D, j


QE = (t ** 2 - 4 * c) * (t - c - 1) * (t - c + 3)
QL = (t ** 2 + (1 + c) * t + 1) ** 2 - 4 * c * t * (t + 1) ** 2
out = {}
for name, Q in (('E\'_c', QE), ('L_c (linkage)', QL)):
    I, J, D, j = invariants(Q)
    out[name] = (I, J, D, j)
    degI, degJ, degD = sp.degree(I, c), sp.degree(J, c), sp.degree(sp.expand(D), c)
    print('%s: I = %s, deg I = %d, deg J = %d;  4I^3 - J^2 = %s (degree %d, so ord_oo Delta = %d)'
          % (name, I, degI, degJ, D, degD, 12 - degD))
    print('    j = %s' % j)
    for c0 in (0, 1, 9):
        ordD = sp.Poly(sp.expand(D), c).as_expr()
        k = 0
        while sp.expand(ordD).subs(c, c0) == 0:
            ordD = sp.cancel(ordD / (c - c0)); k += 1
        print('    c = %d: ord Delta = %d, I(c) = %s -> %s' % (c0, k, I.subs(c, c0), ('I_%d' % k) if I.subs(c, c0) != 0 else 'additive'))
    print('    c = oo: ord Delta = %d, deg I = %d -> %s' % (12 - degD, degI, ('I_%d' % (12 - degD)) if degI == 4 else 'check'))

jE = out['E\'_c'][3]; jL = out['L_c (linkage)'][3]
X, Y = symbols('X Y')
Phi2 = (X ** 3 + Y ** 3 - X ** 2 * Y ** 2 + 1488 * (X ** 2 * Y + X * Y ** 2) - 162000 * (X ** 2 + Y ** 2)
        + 40773375 * X * Y + 8748000000 * (X + Y) - 157464000000000)
val = sp.cancel(together(Phi2.subs({X: jE, Y: jL})))
print('\nPhi_2(j(E\'_c), j(L_c)) =', val, ' (0: the two families are 2-isogenous over Qbar(c))')
sig = (c - 9) / (c - 1)
print('j(L_c) - j(E\'_(sigma(c))), sigma(c) = (c-9)/(c-1):', sp.cancel(jL - jE.subs(c, sig)))
print('j(L_c) - j(E\'_c):', sp.factor(sp.cancel(jL - jE)) != 0 and 'nonzero (not isomorphic over Qbar(c) at the same c)')


# ---------------------------------------------------------------- Weierstrass models and torsion
def weierstrass_from_quartic(Q, t0):
    """y^2 = Q(t), Q monic quartic with rational root t0: t = t0 + 1/s, y = Y/s^2 gives Y^2 = A s^3 + B s^2 + C s + 1;
    X = A s, Y' = A Y gives Y'^2 = X^3 + B X^2 + A C X + A^2."""
    s = symbols('s')
    G = sp.expand(sp.cancel(s ** 4 * Q.subs(t, t0 + 1 / s)))
    P = Poly(G, s)
    A, B, C, D = [P.coeff_monomial(s ** k) for k in (3, 2, 1, 0)]
    assert sp.simplify(D - 1) == 0
    return sp.simplify(A), sp.simplify(B), sp.simplify(A * C), sp.simplify(A ** 2)


def ec_add(P1, P2, a2, a4):
    if P1 is None: return P2
    if P2 is None: return P1
    x1, y1 = P1; x2, y2 = P2
    if sp.simplify(x1 - x2) == 0:
        if sp.simplify(y1 + y2) == 0:
            return None
        lam = sp.cancel((3 * x1 ** 2 + 2 * a2 * x1 + a4) / (2 * y1))
    else:
        lam = sp.cancel((y2 - y1) / (x2 - x1))
    x3 = sp.cancel(lam ** 2 - a2 - x1 - x2)
    y3 = sp.cancel(-(y1 + lam * (x3 - x1)))
    return (x3, y3)


for name, Q, t0 in (('E\'_c', QE, c + 1), ('L_c', QL, None)):
    if t0 is None:
        # L_c: rational points with w = 0? Q_L(t) = 0 has the root t = -1? Q_L(-1) = (1 - 1 - c + 1)^2 - 0 = (1 - c)^2 != 0.
        # use the point at t = 0: Q_L(0) = 1, so (0, +-1) are rational; move it to infinity by t = 1/s' (then quartic in s'
        # with leading coefficient Q_L(0) = 1) -- equivalent: quartic t^4 Q_L(1/t) is monic with constant 1.
        Q2 = sp.expand(t ** 4 * QL.subs(t, 1 / t))
        # find a rational root of Q2 or of QL over Q(c): none in general; use the two-points-at-infinity model directly
        print('\n%s: no rational root over Q(c) used; torsion read from point counts in a4b' % name)
        continue
    A, B, AC, A2 = weierstrass_from_quartic(Q, t0)
    a2, a4, a6 = B, AC, A2
    # points: P = (0, A) (image of a point at infinity of the quartic), T2 = (-A/4, 0) (image of t = c - 3)
    P = (sp.Integer(0), A)
    T2 = (sp.cancel(-A / 4), sp.Integer(0))
    chk = lambda Pt: sp.simplify(Pt[1] ** 2 - (Pt[0] ** 3 + a2 * Pt[0] ** 2 + a4 * Pt[0] + a6)) == 0
    assert chk(P) and chk(T2)
    mult = [None, P]
    for k in range(2, 7):
        mult.append(ec_add(mult[-1], P, a2, a4))
    order = next((k for k in range(1, 7) if mult[k] is None), None)
    print('\n%s: Weierstrass Y^2 = X^3 + (%s) X^2 + (%s) X + (%s)' % (name, sp.factor(a2), sp.factor(a4), sp.factor(a6)))
    print('    the image of a point at infinity of the quartic, P = (0, A), has order %s over Q(c); 3P = %s, T2 = %s'
          % (order, mult[3], T2))

# ---------------------------------------------------------------- 30a1 and the four-bar quartic
a1, a2_, a3, a4_, a6_ = 1, 0, 1, 1, 2
b2 = a1 ** 2 + 4 * a2_; b4 = 2 * a4_ + a1 * a3; b6 = a3 ** 2 + 4 * a6_
b8 = a1 ** 2 * a6_ + 4 * a2_ * a6_ - a1 * a3 * a4_ + a2_ * a3 ** 2 - a4_ ** 2
c4 = b2 ** 2 - 24 * b4; c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
Dl = -b2 ** 2 * b8 - 8 * b4 ** 3 - 27 * b6 ** 2 + 9 * b2 * b4 * b6
print('\n30a1: c4 = %d, c6 = %d, Delta = %d = %s, j = %s = %s' % (c4, c6, Dl, sp.factorint(Dl), sp.Rational(c4 ** 3, Dl), sp.factorint(c4 ** 3)))
Iq, Jq, Dq, jq = invariants(t ** 4 - 6 * t ** 3 - 5 * t ** 2 - 6 * t + 1)
print('four-bar quartic t^4 - 6t^3 - 5t^2 - 6t + 1: I = %s, J = %s, (4I^3 - J^2)/27 = %s = %s, j = %s'
      % (Iq, Jq, sp.Rational(Dq, 27), sp.factorint(sp.Rational(Dq, 27)), jq))
print('L_4 is that quartic:', sp.expand(QL.subs(c, 4) - (t ** 4 - 6 * t ** 3 - 5 * t ** 2 - 6 * t + 1)) == 0)

# twist witnesses
def ap30(p):
    cnt = 1
    for x in range(p):
        for y in range(p):
            if (y * y + x * y + y - x ** 3 - x - 2) % p == 0:
                cnt += 1
    return p + 1 - cnt
primes = list(sp.primerange(7, 140))
ap = {p: ap30(p) for p in primes}
wit = {}
for D in [s * m for m in (1, 2, 3, 5, 6, 10, 15, 30) for s in (1, -1)]:
    if D == 1:
        continue
    wit[D] = next((p for p in primes if sp.jacobi_symbol(D % p, p) == -1 and ap[p] != 0), None)
print('twist witnesses (D: least prime with (D/p) = -1, a_p != 0):', wit)
assert all(v is not None for v in wit.values())

# ---------------------------------------------------------------- the plane cubic: singular points
Fc = 2 * u ** 2 * v + 2 * u * v ** 2 - u ** 2 - 5 * u * v - v ** 2 + 2 * u + 2 * v
Fu, Fv = sp.diff(Fc, u), sp.diff(Fc, v)
G = sp.groebner([Fc, Fu, Fv], u, v, order='lex')
print('\nplane cubic: Groebner basis of (F, F_u, F_v) over Q =', list(G))
R1 = sp.resultant(Fu, Fv, u); R2 = sp.resultant(Fc, Fu, u)
N = sp.resultant(sp.Poly(R1, v), sp.Poly(R2, v))
print('  integer resultant (affine singular points possible only mod its primes):', sp.factorint(N))
# at infinity: homogenize
Fh = sp.expand(w ** 3 * Fc.subs({u: u / w, v: v / w}))
for (uu, vv) in ((1, 0), (0, 1), (1, -1)):
    grad = [sp.diff(Fh, z).subs({u: uu, v: vv, w: 0}) for z in (u, v, w)]
    print('  point at infinity (%d:%d:0): gradient %s' % (uu, vv, grad))
# iota fixed points: u_i = +-1 with sum 2
for p in (3, 5, 7, 11):
    fix = [(a, b, d) for a in (1, -1) for b in (1, -1) for d in (1, -1) if (a + b + d - 2) % p == 0]
    print('  p = %d: iota-fixed points (u_i = +-1, sum 2):' % p, fix)
print('DONE')
