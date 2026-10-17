#!/usr/bin/env python3
"""S39 audit, item 8 and the THM-4664 addendum.

A. Gao's side conditions for the tangent-sweep family (derived here from x | (p(w) + 2g)/g and x^2 | q(w)/g^2 + u at
   x = 0): p(g0) = -2 g0, q(g0) = -g0^2, alpha (p'(g0) + 4) = -g0 (p'(g0) + 2).  Checked on THM-1300 and Gao's G; tested
   on the d = 4, 5 polynomials used in keller_sweep_value_set_law_20261017.py.
B. The actual Keller maps (THM-1300, d = 2; Gao's G, d = 3) evaluated on all of F_p^3: fibre-size histogram over
   F_3 != 0, the S_(d+1) fixed-point law, no fibre of size d (an independent check of the value-set reduction).
C. THM-4664: (i) the map formula = THM-4625 data, det JF = -2; (ii) the (F_2, F_3)-fibre parametrisation is a bijection
   {(b, c, w) : c != 0, p(w) != cb} -> {x != 0, Gamma != 0} for EVERY b, with F = (H/c^2, b, c) (exhaustive, p <= 13);
   (iii) the decomposition S = M - E + S_Gamma verified EXACTLY in Z[zeta_p] (integer exponent histograms, modulo the
   all-ones vector) at every frequency for p = 5, 7, 11; (iv) the large-coefficient set = the cone, with margins, by an
   independent evaluation of F and numpy FFT, 11 <= p <= 41; (v) the Jelonek quartic as the discriminant of the fibre
   cubic, its top-degree part, and the dual curve of 27 a^2 c + b^3 = 0 by elimination.
Usage: python -B a8_keller_and_thm4664.py
"""
import numpy as np
import sympy as sp
from numba import njit
from fractions import Fraction
from math import factorial

x, y, z, w, g0, al = sp.symbols('x y z w g0 alpha')

# ------------------------------------------------------------------ A
def sweep(pw):
    qw = sp.integrate(sp.Symbol('s') / 2 * sp.diff(pw.subs(w, sp.Symbol('s')), sp.Symbol('s')), (sp.Symbol('s'), 0, w))
    return sp.expand(pw), sp.expand(qw)


def side_conditions(pw):
    pw, qw = sweep(pw)
    c1 = sp.expand(pw.subs(w, g0) + 2 * g0)
    c2 = sp.expand(qw.subs(w, g0) + g0 ** 2)
    res = sp.resultant(sp.Poly(sp.cancel(c1 / g0), g0), sp.Poly(sp.cancel(c2 / g0 ** 2), g0))
    return c1, c2, res


print('A. side conditions p(g0) = -2 g0, q(g0) = -g0^2 (common nonzero root needed):')
for name, pw in (('THM-1300 d=2', -3 * w ** 2 + 4 * w), ('Gao G d=3', w ** 3 - 6 * w ** 2 + 6 * w),
                 ('script d=4', w ** 4 - 2 * w ** 3 + 3 * w ** 2 + w), ('script d=5', 2 * w ** 5 + w ** 4 - 3 * w ** 3 + w ** 2 - w)):
    c1, c2, res = side_conditions(pw)
    common = sp.gcd(sp.Poly(sp.cancel(c1 / g0), g0), sp.Poly(sp.cancel(c2 / g0 ** 2), g0))
    print('   %-13s resultant of the two conditions = %s; common factor %s' % (name, res, common.as_expr()))
# scaling: mu p is a member iff P(g0) = int_0^g0 p = 0 for some g0 != 0 (then mu = -2 g0/p(g0))
for name, pw in (('script d=4', w ** 4 - 2 * w ** 3 + 3 * w ** 2 + w), ('script d=5', 2 * w ** 5 + w ** 4 - 3 * w ** 3 + w ** 2 - w)):
    P = sp.integrate(pw, (w, 0, w))
    print('   %s: P(w)/w^2 = %s, factorization over Q: %s' % (name, sp.expand(P / w ** 2), sp.factor(P / w ** 2)))

# ------------------------------------------------------------------ B
def keller_map(pw, gam):
    pw_, qw_ = sweep(pw)
    u = 1 + x * y
    W = gam * u
    F1 = sp.cancel((qw_.subs(w, W) + gam * W) / (gam * x) ** 2)
    F2 = sp.cancel((pw_.subs(w, W) + 2 * gam) / (gam * x))
    F3 = sp.expand(gam * x)
    F = [sp.expand(F1), sp.expand(F2), F3]
    for f in F:
        assert sp.denom(sp.together(f)).free_symbols == set()     # polynomial (a constant denominator is allowed)
    J = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in F]).det()
    return F, sp.simplify(J)


def law(n):
    return [Fraction(1, factorial(k)) * sum(Fraction((-1) ** j, factorial(j)) for j in range(n - k + 1)) for k in range(n + 1)]


print('\nB. actual Keller maps on F_p^3 (fibre sizes over F_3 != 0):')
maps = {'THM-1300': (-3 * w ** 2 + 4 * w, 2 - 3 * x * y - x ** 2 * z, 2),
        'Gao G': (w ** 3 - 6 * w ** 2 + 6 * w, 2 - 4 * x * y - x ** 2 * z, 3)}
FMAP = {}
for name, (pw, gam, d) in maps.items():
    F, J = keller_map(pw, gam)
    FMAP[name] = F
    print('   %s: F = %s; det JF = %s' % (name, [str(f)[:60] + '...' for f in F], J))
    dens = [int(sp.denom(sp.together(f))) for f in F]
    fl = [sp.lambdify((x, y, z), sp.expand(f * D), 'numpy') for f, D in zip(F, dens)]
    for p in (31, 61):
        X, Y, Z = np.meshgrid(np.arange(p, dtype=object), np.arange(p, dtype=object), np.arange(p, dtype=object), indexing='ij')
        vals = [np.vectorize(lambda t: int(t) % p, otypes=[np.int64])(f(X, Y, Z)) * pow(D, -1, p) % p
                for f, D in zip(fl, dens)]
        idx = (vals[0] * p + vals[1]) * p + vals[2]
        N = np.bincount(idx.ravel(), minlength=p ** 3).reshape(p, p, p)
        Nc = N[:, :, 1:]                                       # F_3 = c != 0
        hist = np.bincount(Nc.ravel(), minlength=d + 2)
        L = law(d + 1)
        emp = hist / Nc.size
        print('     p = %d: frequencies %s | S_%d law %s | size-%d fibres over c != 0: %d | image density %.4f vs %.4f'
              % (p, ' '.join('%d:%.4f' % (k, emp[k]) for k in range(d + 2)), d + 1,
                 ' '.join('%.4f' % float(L[k]) for k in range(d + 2)), d, hist[d], 1 - emp[0], 1 - float(L[0])))
        assert hist[d] == 0

# ------------------------------------------------------------------ C
print('\nC. THM-4664')
u = 1 + x * y
G = 2 - 3 * x * y - x ** 2 * z
F1 = u ** 3 * z + y ** 2 * u * (4 + 3 * x * y); F2 = y + 3 * x * u ** 2 * z + 3 * x * y ** 2 * (4 + 3 * x * y); F3 = 2 * x - 3 * x ** 2 * y - x ** 3 * z
print('   (i) THM-4664 / THM-4621 map equals the THM-4625 construction:',
      all(sp.expand(a - b) == 0 for a, b in zip([F1, F2, F3], FMAP['THM-1300'])), '; F_3 = x Gamma:', sp.expand(F3 - x * G) == 0)
f1n = sp.lambdify((x, y, z), F1); f2n = sp.lambdify((x, y, z), F2); f3n = sp.lambdify((x, y, z), F3)


def Fp(v, p):
    a, b, c = v
    return (f1n(a, b, c) % p, f2n(a, b, c) % p, f3n(a, b, c) % p)


for p in (5, 7, 11, 13):
    inv = lambda t: pow(t % p, -1, p)
    seen = {}
    ok = True
    for bb in range(p):
        for cc in range(1, p):
            for ww in range(p):
                pwv = (-3 * ww * ww + 4 * ww) % p
                if pwv == cc * bb % p:
                    continue
                xx = 2 * cc * inv(cc * bb - pwv) % p
                yy = (ww * inv(cc) - inv(xx)) % p
                zz = (2 - 3 * xx * yy - cc * inv(xx)) * inv(xx * xx) % p
                H = (ww ** 3 * inv(2) - ww * ww + cc * bb * inv(2) * ww) % p
                if Fp((xx, yy, zz), p) != (H * inv(cc * cc) % p, bb, cc):
                    ok = False
                key = (xx, yy, zz)
                if key in seen:
                    ok = False
                seen[key] = 1
    target = {(a, b, c) for a in range(1, p) for b in range(p) for c in range(p) if (2 - 3 * a * b - a * a * c) % p}
    print('   (ii) p = %2d: parametrisation injective, F = (H/c^2, b, c) for every b: %s; image = {x != 0, Gamma != 0}: %s'
          % (p, ok, set(seen) == target))
    assert ok and set(seen) == target


@njit
def hist_counts(p, al, be, ga, F1a, F2a, F3a):
    h = np.zeros(p, dtype=np.int64)
    for i in range(F1a.shape[0]):
        h[(al * F1a[i] + be * F2a[i] + ga * F3a[i]) % p] += 1
    return h


def exact_check(p):
    pts = [(a, b, c) for a in range(p) for b in range(p) for c in range(p)]
    Fv = np.array([Fp(v, p) for v in pts], dtype=np.int64)
    F1a, F2a, F3a = Fv[:, 0].copy(), Fv[:, 1].copy(), Fv[:, 2].copy()
    inv = [0] + [pow(t, -1, p) for t in range(1, p)]
    # E terms and S_Gamma terms as point lists
    Epts = [(c, ww) for c in range(1, p) for ww in range(p)]
    Gpts = []
    for xx in range(1, p):
        for yy in range(p):
            zz = (2 - 3 * xx * yy) * inv[xx] * inv[xx] % p
            Gpts.append(Fp((xx, yy, zz), p))
    Gp = np.array(Gpts, dtype=np.int64)
    worst = 0
    for a_ in range(p):
        for b_ in range(p):
            for c_ in range(p):
                if a_ == 0:
                    continue
                hS = hist_counts(p, a_, b_, c_, F1a, F2a, F3a)
                hE = np.zeros(p, dtype=np.int64)
                for (c, ww) in Epts:
                    qw = (-ww ** 3 + ww * ww) % p; pw = (-3 * ww * ww + 4 * ww) % p
                    hE[(a_ * qw * inv[c] * inv[c] + b_ * pw * inv[c] + c_ * c) % p] += 1
                hG = hist_counts(p, a_, b_, 0, Gp[:, 0].copy(), Gp[:, 1].copy(), Gp[:, 2].copy())
                hM = np.zeros(p, dtype=np.int64)
                ph = (-4 * b_ * b_ * inv[a_]) % p
                onc = (a_ * a_ * c_ - 4 * b_ ** 3) % p == 0
                hM[ph] += p * (p * (1 if onc else 0) - 1)
                diff = hS - (hM - hE + hG)
                if not np.all(diff == diff[0]):          # zero in Z[zeta_p] iff constant
                    worst += 1
    return worst


for p in (5, 7, 11):
    print('   (iii) p = %2d: exact Z[zeta_p] check of S = M - E + S_Gamma at every frequency with alpha != 0: failures %d'
          % (p, exact_check(p)))


def N_counts(p):
    a = np.arange(p, dtype=np.int64)
    X, Y, Z = np.meshgrid(a, a, a, indexing='ij')
    U = (1 + X * Y) % p
    f1 = (U * U % p * U % p * Z + Y * Y % p * U % p * ((4 + 3 * X * Y) % p)) % p
    f2 = (Y + 3 * X % p * U % p * U % p * Z + 3 * X % p * (Y * Y % p) % p * ((4 + 3 * X * Y) % p)) % p
    f3 = (2 * X - 3 * X * X % p * Y - X * X % p * X % p * Z) % p
    return np.bincount(((f1 * p + f2) * p + f3).ravel(), minlength=p ** 3).reshape(p, p, p)


print('   (iv) large coefficients vs the cone alpha^2 gamma = 4 beta^3 (independent evaluation, numpy FFT):')
for p in (11, 13, 17, 19, 23, 29, 31, 37, 41):
    N = N_counts(p).astype(float)
    S = np.abs(np.fft.fftn(N)); S[0, 0, 0] = 0
    a = np.arange(p)
    A_, B_, C_ = np.meshgrid(a, a, a, indexing='ij')
    cone = ((A_ * A_ % p) * C_ - 4 * (B_ ** 3 % p)) % p == 0
    cone[0, 0, 0] = False
    big = S > p * p / 2
    print('        p = %2d: sets equal %s; min on cone %.4f p^2, max off cone %.4f p^2 (margin to p^2/2: %.3f, %.3f)'
          % (p, bool((big == cone).all()), S[cone].min() / p ** 2, S[~cone].max() / p ** 2,
             S[cone].min() / p ** 2 - 0.5, 0.5 - S[~cone].max() / p ** 2))

a, b, c = sp.symbols('a b c')
H = w ** 3 / 2 - w ** 2 + c * b * w / 2 - c ** 2 * a
disc = sp.factor(sp.discriminant(H, w))
c3 = 27 * a ** 2 * c ** 2 - 18 * a * b * c + 16 * a + b ** 3 * c - b ** 2
print('   (v) disc_w(H_(b,c)(w) - c^2 a) =', disc, '; ratio to c_3 =', sp.simplify(disc / c3))
top = sp.Poly(c3, a, b, c)
deg = top.total_degree()
toppart = sum(cf * a ** m[0] * b ** m[1] * c ** m[2] for m, cf in zip(top.monoms(), top.coeffs()) if sum(m) == deg)
print('       top-degree part of c_3:', sp.factor(toppart))
A_, B_, C_, t = sp.symbols('A B C t')
# dual of 27 a^2 c + b^3 = 0: tangent line coordinates (A:B:C) = gradient
Fcub = 27 * a ** 2 * c + b ** 3
grad = [sp.diff(Fcub, v) for v in (a, b, c)]
par = {a: 1, b: -3 * t, c: t ** 3}
gl = [sp.factor(gg.subs(par)) for gg in grad]
print('       gradient along (1 : -3t : t^3):', gl, '; satisfies A^2 C = 4 B^3:',
      sp.simplify(gl[0] ** 2 * gl[2] - 4 * gl[1] ** 3) == 0)
print('DONE')
