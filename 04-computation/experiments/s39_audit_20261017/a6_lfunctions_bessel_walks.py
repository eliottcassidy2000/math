#!/usr/bin/env python3
"""S39 audit, item 5b and the p_3 addendum: THM-4663 statement 7 and the real-place identity of statement 9.

L-functions: Lambda(s) = N^(s/2) (2 pi)^(-s) Gamma(s) L(f, s) evaluated by the theta splitting at an ARBITRARY point
y0 = A/sqrt(N):  Lambda(s) = sum_n a_n [ N^(s/2) (2 pi n)^(-s) Gamma(s, 2 pi n A/sqrt N)
                                         + eps N^((k-s)/2) (2 pi n)^(-(k-s)) Gamma(k-s, 2 pi n/(A sqrt N)) ];
independence of A (A = 1, 1.37, 0.71) certifies the functional equation with root number eps; a_n are taken from the
q-expansions of the eta quotients (a2), not from a Hecke recursion.
Bessel moments by mpmath quadrature (40 digits).
Walk densities without Bessel functions: p_2(s) = 2/(pi sqrt(4 - s^2)); adding a unit step to a vector of length s gives
radius r with density 2r/(pi sqrt(4 s^2 - (r^2 - s^2 - 1)^2)).  Hence (substituting t = r^2 + 1 - s^2)
    p_3(r) = (2r/pi^2) Omega(r^2),  Omega(c) = int dt/sqrt(Q_c(t)) over the bounded interval where Q_c > 0,
an exact identity; and p_n'(0) = p_(n-1)(1), so p_5'(0) = p_4(1) = int_0^2 p_3(s) 2/(pi s sqrt(4 - s^2)) ds
    = (4/pi^3) int_0^2 Omega(s^2) / sqrt(4 - s^2) ds.
Omega by the Legendre-form/AGM formula, cross-checked by direct quadrature.
Usage: python -B a6_lfunctions_bessel_walks.py
"""
import json
import mpmath as mp

mp.mp.dps = 40
F = json.load(open('a2_forms_to_3000.json'))
f15 = F['f15']; b46 = F['b46']


def Lam(a, N, k, s, eps, A, nmax=400):
    sN = mp.sqrt(N)
    tot = mp.mpf(0)
    for n in range(1, nmax + 1):
        if a[n] == 0:
            continue
        x1 = 2 * mp.pi * n * A / sN
        x2 = 2 * mp.pi * n / (A * sN)
        tot += a[n] * (mp.power(N, s / 2) * mp.power(2 * mp.pi * n, -s) * mp.gammainc(s, x1)
                       + eps * mp.power(N, (k - s) / 2) * mp.power(2 * mp.pi * n, -(k - s)) * mp.gammainc(k - s, x2))
    return tot


def Lval(a, N, k, s, eps, A=1):
    s = mp.mpf(s)
    return Lam(a, N, k, s, eps, A) / (mp.power(N, s / 2) * mp.power(2 * mp.pi, -s) * mp.gamma(s))


print('functional-equation test (A-independence) and root numbers:')
for name, a, N, k, pts in (('f15 (weight 3, level 15)', f15, 15, 3, [1, 2, mp.mpf('1.3')]),
                           ('b = (eta1 eta2 eta3 eta6)^2 (weight 4, level 6)', b46, 6, 4, [1, 2, 3, mp.mpf('2.6')])):
    for eps in (1, -1):
        spread = max(abs(Lval(a, N, k, s, eps, 1) - Lval(a, N, k, s, eps, A)) for s in pts for A in (mp.mpf('1.37'), mp.mpf('0.71')))
        print('   %s, eps = %+d: max |L(A=1) - L(A)| over test points = %s' % (name, eps, mp.nstr(spread, 3)))
L15 = {s: Lval(f15, 15, 3, s, 1) for s in (1, 2)}
L6 = {s: Lval(b46, 6, 4, s, 1) for s in (1, 2, 3)}
print('   L(f15,1) = %s, L(f15,2) = %s' % (mp.nstr(L15[1], 30), mp.nstr(L15[2], 30)))
print('   L(b,1) = %s, L(b,2) = %s, L(b,3) = %s' % (mp.nstr(L6[1], 30), mp.nstr(L6[2], 30), mp.nstr(L6[3], 30)))

I0, K0 = mp.besseli, mp.besselk
pts = [0, 1, 5, 20, 60, mp.inf]
IK4 = mp.quad(lambda t: t * I0(0, t) * K0(0, t) ** 4, pts)
I2K3 = mp.quad(lambda t: t * I0(0, t) ** 2 * K0(0, t) ** 3, pts)
IK5 = mp.quad(lambda t: t * I0(0, t) * K0(0, t) ** 5, pts)
I2K4 = mp.quad(lambda t: t * I0(0, t) ** 2 * K0(0, t) ** 4, pts)
pi = mp.pi
rows = [('int t I0 K0^4', IK4, pi ** 2 * L15[1] / 5), ('int t I0^2 K0^3', I2K3, 3 * L15[2] / 4),
        ('int t I0 K0^5', IK5, pi ** 2 * L6[2] / 2), ('int t I0^2 K0^4 vs pi^2 L6(1)/2', I2K4, pi ** 2 * L6[1] / 2),
        ('int t I0^2 K0^4 vs 3 L6(3)/2', I2K4, 3 * L6[3] / 2), ('L15(1)/L15(2) vs sqrt15/(2pi)', L15[1] / L15[2], mp.sqrt(15) / (2 * pi))]
print('\nBessel moments vs L-values (40-digit arithmetic):')
for nm, x, y in rows:
    print('   %-34s %s  vs  %s   |diff| = %s' % (nm, mp.nstr(x, 25), mp.nstr(y, 25), mp.nstr(abs(x - y), 3)))


# ---------------------------------------------------------------- walk densities
def Omega(c):
    """int dt/sqrt(Q_c(t)) over the bounded interval with Q_c > 0, Q_c = (t^2 - 4c)(t - c - 1)(t - c + 3), 0 < c < 9, c != 1.
    Legendre form with the root differences factored (no cancellation near c = 1):
      c < 1: (e4-e2)(e3-e1) = (1+r)^3 (3-r),  k'^2 = (1-r)^3 (3+r) / ((1+r)^3 (3-r));
      c > 1: (e4-e2)(e3-e1) = 16 r,           k'^2 = (r-1)^3 (r+3) / (16 r);       r = sqrt(c)."""
    r = mp.sqrt(mp.mpf(c))
    if r < 1:
        den = (1 + r) ** 3 * (3 - r)
        kp2 = (1 - r) ** 3 * (3 + r) / den
    else:
        den = 16 * r
        kp2 = (r - 1) ** 3 * (r + 3) / den
    if kp2 == 0:
        kp2 = mp.mpf(10) ** (-(3 * mp.mp.dps))      # measure-zero node at the log singularity
    return mp.pi / (mp.agm(1, mp.sqrt(kp2)) * mp.sqrt(den))


def Omega_quad(c):
    c = mp.mpf(c)
    rts = sorted([-2 * mp.sqrt(c), 2 * mp.sqrt(c), c + 1, c - 3])
    e2, e3 = rts[1], rts[2]
    Q = lambda t: (t * t - 4 * c) * (t - c - 1) * (t - c + 3)
    return mp.quad(lambda t: 1 / mp.sqrt(Q(t)), [e2, (e2 + e3) / 2, e3])


def p3_bswz(r):
    r = mp.mpf(r)
    z = r ** 2 * (9 - r ** 2) ** 2 / (3 + r ** 2) ** 3
    return 2 * mp.sqrt(3) / pi * r / (3 + r ** 2) * mp.hyp2f1(mp.mpf(1) / 3, mp.mpf(2) / 3, 1, z)


def p3_conv(r):
    """direct convolution p_3(r) = int p_2(s) * 2r/(pi sqrt(4s^2 - (r^2 - s^2 - 1)^2)) ds over |r-1| < s < min(2, r+1)."""
    r = mp.mpf(r)
    lo, hi = abs(r - 1), min(mp.mpf(2), r + 1)
    g = lambda s: (2 / (pi * mp.sqrt(4 - s * s))) * 2 * r / (pi * mp.sqrt(4 * s * s - (r * r - s * s - 1) ** 2))
    return mp.quad(g, [lo, (lo + hi) / 2, hi])


print('\np_3(r) = (2/pi^2) r Omega(r^2) (proved by the convolution substitution; numerical confirmation):')
mp.mp.dps = 30
for r in ('0.2', '0.5', '0.9', '1.1', '1.7', '2.5', '2.9'):
    a = 2 / pi ** 2 * mp.mpf(r) * Omega(mp.mpf(r) ** 2)
    b = p3_bswz(r)
    cq = p3_conv(r)
    oq = 2 / pi ** 2 * mp.mpf(r) * Omega_quad(mp.mpf(r) ** 2)
    print('   r = %s: (2/pi^2) r Omega = %s;  direct convolution %s;  Omega by quadrature %s;  BSWZ 2F1 %s;  max|diff| %s'
          % (r, mp.nstr(a, 22), mp.nstr(cq, 22), mp.nstr(oq, 22), mp.nstr(b, 22), mp.nstr(max(abs(a - b), abs(a - cq), abs(a - oq)), 3)))
mp.mp.dps = 40
p5 = 4 / pi ** 3 * mp.quad(lambda s: Omega(s * s) / mp.sqrt(4 - s * s), [0, 1, 2])
print('\np_5\'(0) = p_4(1) = (4/pi^3) int_0^2 Omega(s^2)/sqrt(4 - s^2) ds = %s' % mp.nstr(p5, 30))
print('6 L(f15,1)/pi^2                                   = %s' % mp.nstr(6 * L15[1] / pi ** 2, 30))
print('(30/pi^4) int t I0 K0^4                           = %s' % mp.nstr(30 * IK4 / pi ** 4, 30))
print('|diff| = %s' % mp.nstr(abs(p5 - 6 * L15[1] / pi ** 2), 3))
mp.mp.dps = 25
J5 = mp.quadosc(lambda x: x * mp.besselj(0, x) ** 5, [0, mp.inf], period=2 * mp.pi)
print('int_0^oo x J0(x)^5 dx (mpmath quadosc, period 2 pi)  = %s' % mp.nstr(J5, 22))
