#!/usr/bin/env python3
"""Referee check B1 for THM-4668 (second S41 audit): every algebraic step of the proof, recomputed independently.

1. Norm components of x = c0 + c1 s + c2 h + c3 sh under conj (c_i -> conj c_i, s -> -s, h -> h), s^2 = D, h^2 = R.
2. Coordinates of x = a + u b + v c + e g, v = Z(-1+s)/2, e = Z s/h.
3. Case c != 0 and case c = 0 identities (with D = 1 - 4/N, R = 4 - N, N = Z Zbar).
4. The coefficient equations of (a + u b)^2 (1 + dbar u) - k u (d + u) modulo u^3 = m2 u^2 + m1 u + m0.
5. The minimal polynomial of u over Q(j), recomputed from h^3 - 7h^2 + 11h - 1 (h = w + 3), and irreducibility.
6. Field degrees: [Q(zeta15, j):Q] = 16 and [F:K] = 3 (prime certificate), so [F:Q] = 48.
7. High-precision values of N, D, R at the intended embedding and at tau (j -> -j).
8. Own implementation of Q(zeta5, j) arithmetic: lambda, lambda^2 + dbar, its absolute norm, the resultant form, and
   comparison with the numerator printed in golden_e_direction_proof_check_20261017.out.
Usage: python b1_proof_algebra.py
"""
from fractions import Fraction as Fr
import itertools
import sympy as sp
import mpmath as mp

print('=== 1. norm components ===')
c0, c1, c2, c3, b0, b1, b2, b3, S, H, D, R = sp.symbols('c0 c1 c2 c3 cb0 cb1 cb2 cb3 s h D R')
x = c0 + c1 * S + c2 * H + c3 * S * H
xb = b0 - b1 * S + b2 * H - b3 * S * H            # conj: c_i -> conj c_i, s -> -s, h -> h
prod = sp.expand(x * xb)
prod = sp.expand(prod.subs(S ** 2, D).subs(H ** 2, R))
poly = sp.Poly(prod, S, H)
comp = {m: sp.factor(poly.coeff_monomial(m)) for m in [1, S, H, S * H]}
for m, v in comp.items():
    print('  coeff of', m, ':', sp.expand(v))
c2zero = {c2: 0, b2: 0}
eq_i = sp.expand(comp[1].subs(c2zero) - (c0 * b0 - D * c1 * b1 - D * R * c3 * b3))
eq_ii = sp.expand(comp[S].subs(c2zero) - (c1 * b0 - c0 * b1))
eq_iii = sp.expand(comp[H].subs(c2zero) + D * (c1 * b3 + c3 * b1))
eq_iv = sp.expand(comp[S * H].subs(c2zero) - (c3 * b0 - c0 * b3))
print('  with c2 = 0: (i), (ii), (iii)*(-D), (iv) residuals:', eq_i, eq_ii, eq_iii, eq_iv)

print('=== 2. coordinates of x = a + u b + v c + e g ===')
a, b, c, g, u, Z = sp.symbols('a b c g u Z')
v = Z * (-1 + S) / 2
e = Z * S / H
xx = sp.expand((a + u * b + v * c + e * g) * H * R)   # multiply by h R to clear 1/h; h^2 = R
xx = sp.expand(xx.subs(H ** 2, R))
# x = (a + u b - Z c/2) + (Z c/2) s + 0 h + (Z g/R) s h ; check x*h*R = h R c0 + h R c1 s + R c3 s h^2 ...
target = sp.expand(((a + u * b - Z * c / 2) + (Z * c / 2) * S + (Z * g / R) * S * H) * H * R)
target = sp.expand(target.subs(H ** 2, R))
print('  x - [c0 + c1 s + c3 sh] (times hR):', sp.simplify(xx - target))

print('=== 3. case identities ===')
N, Zb, cb, gb, t = sp.symbols('N Zb cb gb t')
Dv = 1 - 4 / N; Rv = 4 - N
C1, C1b = Z * c / 2, Zb * cb / 2
C3, C3b = Z * g / Rv, Zb * gb / Rv
iii = sp.simplify((C1 * C3b + C3 * C1b).subs(Z * Zb, N))
iii = sp.simplify(iii.subs(Zb, N / Z))
print('  (iii) c1 c3b + c3 c1b =', sp.factor(iii), ' [claimed (N/(2R))(c gb + g cb)]')
i_c0zero = sp.simplify((-Dv * C1 * C1b - Dv * Rv * C3 * C3b).subs(Zb, N / Z))
print('  (i) with c0 = 0:', sp.factor(i_c0zero), ' minus claim:', sp.simplify(i_c0zero - (Rv * c * cb / 4 + g * gb)))
C0, C0b = t * Z * g, t * Zb * gb
i_c0 = sp.simplify((C0 * C0b - Dv * Rv * C3 * C3b).subs(Zb, N / Z))
print('  (i) with c = 0, c0 = t Z g:', sp.factor(i_c0), ' minus claim:', sp.simplify(i_c0 - g * gb * (1 + t ** 2 * N)))
# c0^2 N = k Z^2 g^2, N = Z Zb, Zb = dbar + 1/u  =>  c0^2 (1 + dbar u) = k g^2 u (d + u)
dd, db, k, A0 = sp.symbols('d dbar kappa c0val')
lhs = A0 ** 2 * Z * (db + 1 / u) - k * Z ** 2 * g ** 2        # c0^2 N - k Z^2 g^2 with N = Z (dbar + 1/u)
lhs = sp.expand(lhs.subs(Z, dd + u) * u / (dd + u))
rhs = sp.expand(A0 ** 2 * (1 + db * u) - k * g ** 2 * u * (dd + u))
print('  u (c0^2 N - k Z^2 g^2)/Z - [c0^2 (1 + dbar u) - k g^2 u (d + u)] =', sp.simplify(lhs - rhs))

print('=== 4. coefficient equations ===')
m0, m1, m2, kp = sp.symbols('m0 m1 m2 kp')
expr = sp.expand((a + u * b) ** 2 * (1 + db * u) - kp * u * (dd + u))
rem = sp.rem(sp.Poly(expr, u), sp.Poly(u ** 3 - m2 * u ** 2 - m1 * u - m0, u))
co = rem.all_coeffs()[::-1]
claimed = [a ** 2 + db * b ** 2 * m0, 2 * a * b + db * a ** 2 + db * m1 * b ** 2 - kp * dd, b ** 2 + 2 * db * a * b + db * m2 * b ** 2 - kp]
for i, (x1, x2) in enumerate(zip(co, claimed)):
    print('  [u^%d] coefficient:' % i, sp.expand(x1), ' | minus claimed:', sp.expand(x1 - x2))
# lambda derivation: a = lam b, b != 0; [1] lam^2 = -dbar m0; [u^2] kp = b^2 (1 + 2 dbar lam + dbar m2); [u] then
lam = sp.symbols('lam')
eq_u = sp.expand((co[1].subs(a, lam * b).subs(kp, b ** 2 * (1 + 2 * db * lam + db * m2))) / b ** 2)
print('  [u] after a = lam b and kp from [u^2]:', sp.collect(eq_u, lam), '  (quadratic in lam: needs [1] to linearise)')
eq_u_lin = sp.expand(eq_u.subs(lam ** 2, -db * m0))
lam_sol = sp.solve(eq_u_lin, lam)[0]
lam_claim = (db ** 2 * m0 - db * m1 + dd + dd * db * m2) / (2 - 2 * dd * db)
print('  lambda from [u] with lam^2 = -dbar m0:', sp.simplify(lam_sol), '; minus THM/script formula:', sp.simplify(lam_sol - lam_claim))

print('=== 5. minimal polynomial of u over Q(j) ===')
J, Hh, U = sp.symbols('j h U')
cub_h = sp.expand((Hh - 3) ** 3 + 2 * (Hh - 3) ** 2 - 4 * (Hh - 3) - 4)
print('  cubic in h = w + 3:', cub_h)
h_of_u = J * (1 + U) / (1 - U)
num = sp.expand(sp.numer(sp.together(cub_h.subs(Hh, h_of_u))))
num = sp.expand(num.subs(J ** 3, -11 * J).subs(J ** 2, -11))
P = sp.Poly(num, U)
cs = P.all_coeffs()
mon = [sp.nsimplify(sp.radsimp(ci / cs[0])) for ci in cs]
print('  monic cubic: U^3 + (%s) U^2 + (%s) U + (%s)' % tuple(sp.simplify(m.subs(J, sp.sqrt(-11)).subs(sp.sqrt(-11), J)) for m in mon[1:]))
mm2 = sp.simplify(-mon[1]); mm1 = sp.simplify(-mon[2]); mm0 = sp.simplify(-mon[3])
print('  u^3 = m2 u^2 + m1 u + m0 with m2 =', mm2, ', m1 =', mm1, ', m0 =', mm0)
print('  equals THM values (m2 = -(20+11j)/19, m1 = (20-11j)/19, m0 = 1):',
      sp.simplify(mm2 + (20 + 11 * J) / 19) == 0, sp.simplify(mm1 - (20 - 11 * J) / 19) == 0, sp.simplify(mm0 - 1) == 0)
X = sp.symbols('X')
cubU = sp.Poly(X ** 3 - mm2.subs(J, sp.sqrt(-11)) * X ** 2 - mm1.subs(J, sp.sqrt(-11)) * X - mm0, X, extension=sp.sqrt(-11))
fl = sp.factor_list(cubU.as_expr(), X, extension=sp.sqrt(-11))
print('  factorisation over Q(sqrt(-11)):', [(sp.degree(f, X), m) for f, m in fl[1]])
mp.mp.dps = 50
wr = [r for r in mp.polyroots([1, 2, -4, -4]) if abs(mp.im(r)) < mp.mpf(10) ** -40 and 1 < mp.re(r) < 2][0]
wr = mp.re(wr)
jn = mp.mpc(0, mp.sqrt(11))
un = (wr + 3 - jn) / (wr + 3 + jn)
m2n = -(20 + 11 * jn) / 19; m1n = (20 - 11 * jn) / 19
print('  w =', mp.nstr(wr, 30), '; |u^3 - m2 u^2 - m1 u - 1| at the embedding (50 digits):', mp.nstr(abs(un ** 3 - m2n * un ** 2 - m1n * un - 1), 5))

print('=== 6. field degrees ===')
x5 = sp.symbols('x5')
fl15 = sp.factor_list(sp.cyclotomic_poly(15, x5), x5, extension=sp.sqrt(-11))
print('  Phi_15 over Q(sqrt(-11)) factor degrees:', [sp.degree(f, x5) for f, m in fl15[1]], '=> [K:Q] = 16' if len(fl15[1]) == 1 else '')
cert = []
for p in sp.primerange(17, 5000):
    if p % 15 != 1 or sp.legendre_symbol(p - 11 % p, p) != 1:
        continue
    roots = [y for y in range(p) if (y ** 3 + 2 * y * y - 4 * y - 4) % p == 0]
    if not roots:
        cert.append(p)
    if len(cert) == 3:
        break
print('  primes p = 1 mod 15 with (-11/p) = 1 (split completely in K) where w^3+2w^2-4w-4 has NO root mod p:', cert)
print('  => the monic integral cubic has no root in K (a root would be integral in K and reduce to F_p), so [F:K] = 3, [F:Q] = 48')
disc = sp.discriminant(X ** 3 + 2 * X ** 2 - 4 * X - 4, X)
print('  disc of the cubic over Q:', disc, sp.factorint(disc), '(not a square: S3 cubic, three real roots)')

print('=== 7. N, D, R at the intended embedding and at tau ===')
z5 = mp.exp(2j * mp.pi / 5)
for name, jj in (('intended (j = +i sqrt11)', jn), ('tau (j = -i sqrt11)', -jn)):
    uu = (wr + 3 - jj) / (wr + 3 + jj)
    Zn = 1 + z5 + uu
    Nn = abs(Zn) ** 2
    print('  %s: u = %s, N = %s, 37/15 < N < 59/20: %s, R = 4 - N = %s, D = 1 - 4/N = %s' % (
        name, mp.nstr(uu, 12), mp.nstr(Nn, 20), 37 / 15 < Nn < 59 / 20, mp.nstr(4 - Nn, 15), mp.nstr(1 - 4 / Nn, 15)))
print('  (at the intended embedding R > 0, so h = sqrt(R) is real and conj fixes h; D < 0 so s is imaginary and conj s = -s)')
# other embeddings of zeta5 at the same u: for the record (R changes sign for some)
for kk in (1, 2, 3, 4):
    for jj in (jn, -jn):
        for wv in [mp.re(r) for r in mp.polyroots([1, 2, -4, -4])]:
            uu = (wv + 3 - jj) / (wv + 3 + jj)
            Nn = abs(1 + mp.exp(2j * mp.pi * kk / 5) + uu) ** 2
print('  (other embeddings are irrelevant: the theorem is about |x| = 1 at the intended embedding)')


print('=== 8. exact lambda^2 + dbar in Q(zeta5, j), own arithmetic ===')
# element: dict {(i, k): Fraction} with i < 4 (zeta5 power), k < 2 (j power)
BASIS = [(i, k) for k in range(2) for i in range(4)]


def norm_elt(dct):
    out = {}
    for (i, k), cval in dct.items():
        if cval == 0:
            continue
        # reduce j^k: j^2 = -11
        coef = Fr(cval)
        while k >= 2:
            k -= 2; coef *= -11
        # reduce zeta^i: zeta^5 = 1; zeta^4 = -1 - z - z^2 - z^3
        i %= 5
        if i == 4:
            for i2 in range(4):
                out[(i2, k)] = out.get((i2, k), Fr(0)) - coef
        else:
            out[(i, k)] = out.get((i, k), Fr(0)) + coef
    return {kk: vv for kk, vv in out.items() if vv != 0}


def add(x1, x2):
    out = dict(x1)
    for kk, vv in x2.items():
        out[kk] = out.get(kk, Fr(0)) + vv
    return {kk: vv for kk, vv in out.items() if vv != 0}


def neg(x1):
    return {kk: -vv for kk, vv in x1.items()}


def mul(x1, x2):
    out = {}
    for (i1, k1), v1 in x1.items():
        for (i2, k2), v2 in x2.items():
            kk = (i1 + i2, k1 + k2)
            out[kk] = out.get(kk, Fr(0)) + v1 * v2
    return norm_elt(out)


def sc(x1, q):
    return {kk: vv * Fr(q) for kk, vv in x1.items() if vv * Fr(q) != 0}


def const(q):
    return {(0, 0): Fr(q)} if q != 0 else {}


def vec(x1):
    return [x1.get(bb, Fr(0)) for bb in BASIS]


def mulmatrix(x1):
    cols = []
    for bb in BASIS:
        cols.append(vec(mul(x1, {bb: Fr(1)})))
    return sp.Matrix(8, 8, lambda r, cidx: sp.Rational(cols[cidx][r].numerator, cols[cidx][r].denominator))


def inv(x1):
    Mx = mulmatrix(x1)
    sol = Mx.LUsolve(sp.Matrix([1, 0, 0, 0, 0, 0, 0, 0]))
    out = {bb: Fr(int(sp.numer(sol[r])), int(sp.denom(sol[r]))) for r, bb in enumerate(BASIS)}
    out = {kk: vv for kk, vv in out.items() if vv != 0}
    assert mul(out, x1) == const(1)
    return out


def norm_Q(x1):
    return mulmatrix(x1).det()


z = {(1, 0): Fr(1)}
jj = {(0, 1): Fr(1)}
one = const(1)
d = add(one, z)
dbar = add(one, norm_elt({(4, 0): Fr(1)}))
dd = mul(d, dbar)
m1e = sc(add(const(20), sc(jj, -11)), Fr(1, 19))
m2e = sc(add(const(20), sc(jj, 11)), Fr(-1, 19))
print('  |d|^2 =', dd, ' (= 2 + z + z^4; equals 1?', dd == one, ')')
q1 = add(const(2), sc(dd, -2))
q0 = add(add(mul(dbar, m1e), neg(d)), neg(mul(dd, m2e)))
lam = mul(add(mul(dbar, dbar), neg(q0)), inv(q1))       # (dbar^2 - dbar m1 + d + |d|^2 m2)/(2 - 2|d|^2)
Eval = add(mul(lam, lam), dbar)
print('  lambda =', {kk: str(vv) for kk, vv in lam.items()})
print('  lambda^2 + dbar =', {kk: str(vv) for kk, vv in Eval.items()})
print('  lambda^2 + dbar is zero?', Eval == {})
nE = norm_Q(Eval)
nEr = sp.Rational(nE)
print('  absolute norm N_{Q(zeta5,j)/Q}(lambda^2 + dbar) =', nEr, '= num', sp.factorint(sp.numer(nEr)), '/ den', sp.factorint(sp.denom(nEr)))
# resultant form (no linearisation): Res_lam(lam^2 + dbar, dbar lam^2 + q1 lam + q0) = (q0 - dbar^2)^2 + q1^2 dbar
res = add(mul(add(q0, neg(mul(dbar, dbar))), add(q0, neg(mul(dbar, dbar)))), mul(mul(q1, q1), dbar))
print('  Res = (q0 - dbar^2)^2 + q1^2 dbar equals q1^2 (lambda^2 + dbar):', res == mul(mul(q1, q1), Eval), '; Res nonzero:', res != {})
# their printed numerator / denominator
num_theirs = {(3, 1): Fr(-44), (2, 1): Fr(396), (1, 1): Fr(880), (0, 1): Fr(902), (3, 0): Fr(-3770), (2, 0): Fr(-5788), (1, 0): Fr(-5195), (0, 0): Fr(-893)}
den_theirs = {(3, 0): Fr(-1444), (2, 0): Fr(-1444), (0, 0): Fr(1444)}
print('  printed numerator / printed den^2 equals my lambda^2 + dbar:', mul(Eval, den_theirs) == num_theirs)
# numeric check at all 8 embeddings
vals = []
for kk in (1, 2, 3, 4):
    for sgn in (1, -1):
        zz = mp.exp(2j * mp.pi * kk / 5); jv = sgn * mp.mpc(0, mp.sqrt(11))
        val = sum(mp.mpf(vv.numerator) / vv.denominator * zz ** i * jv ** k2 for (i, k2), vv in Eval.items())
        vals.append(abs(val))
print('  |lambda^2 + dbar| over the 8 embeddings: min %s, max %s' % (mp.nstr(min(vals), 8), mp.nstr(max(vals), 8)))
# also the alternative: what if one (wrongly) used only [u] and [u^2]: both roots of the quadratic, check lam^2 + dbar at each
print('  (cross-check) roots of dbar lam^2 + q1 lam + q0 at the intended embedding and lam^2 + dbar there:')
zz = mp.exp(2j * mp.pi / 5); jv = mp.mpc(0, mp.sqrt(11))
ev = lambda X1: sum(mp.mpf(vv.numerator) / vv.denominator * zz ** i * jv ** k2 for (i, k2), vv in X1.items())
A2, A1, A0n = ev(dbar), ev(q1), ev(q0)
disc_ = mp.sqrt(A1 ** 2 - 4 * A2 * A0n)
for r in ((-A1 + disc_) / (2 * A2), (-A1 - disc_) / (2 * A2)):
    print('     root', mp.nstr(r, 10), ' |root^2 + dbar| =', mp.nstr(abs(r ** 2 + A2), 10))
