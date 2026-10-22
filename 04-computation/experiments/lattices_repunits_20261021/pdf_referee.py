"""Referee checks for three owner-supplied PDFs (session opus-2026-10-21-S52; note
05-knowledge/results/lattices_modular_repunits_20261021.md).  The PDFs are data; nothing here reads them.

(D) Modular identities quoted by the owner: Delta = (E4^3 - E6^2)/1728, E8 = E4^2, E10 = E4 E6, E14 = E4^2 E6, and
    Delta = q prod (1 - q^n)^24 = eta^24 (the q^(1/24) of eta), exact integer q-expansions to q^40.
(N) 'Sub-double-norm Voronoi relevance ... dimension 48' (nebe_voronoi_rigidity.pdf): the extremal weight-24 theta
    series E4^6 - 1440 E4^3 Delta + 125280 Delta^2 (coefficients forced by a_1 = a_2 = 0), its shells 52,416,000,
    39,007,332,000, 6,609,020,221,440, and the density (3 pi/2)^24/24!.  Also the Leech theta series E4^3 - 720 Delta.
(M) 'Arithmetic of the McLaughlin slice' (mclaughlin_slice_arithmetic.pdf): the even rank-2 forms of determinant 15
    with minimum >= 4 (Gauss reduction), the norm form of p_2 = (2, omega) in Q(sqrt(-15)), and Proposition 6.1
    (reduced forms with a = c = 2 and fundamental discriminant).
(E) 'E8 mass kernel' (e8.pdf): zeta_E8(10) = 240 * 2^-5 zeta(5) zeta(2) = (5/4) pi^2 zeta(5) (classical), checked by a
    direct convergent sum; and the correct Leech analogue at s = 10 (only by analytic continuation, the series
    converges for s > 24): zeta_Leech(s) = (65520/691) 2^(-s/2) [zeta(s/2) zeta(s/2 - 11) - L(Delta, s/2)], so
    zeta_Leech(10) = -(65520/691) 2^-5 L(Delta, 5) because zeta(-6) = 0.  Inside the convergent range the 'odd zeta'
    Leech value the PDF anticipates carries an irremovable cusp term: zeta_Leech(26) = (65520/691) 2^-13 [zeta(13)
    zeta(2) - L(Delta, 13)], checked against the direct shell sum.  L(Delta, s) by the rapidly convergent
    incomplete-gamma formula, validated against the absolutely convergent Dirichlet series at s = 10 and s = 13.
The 'Bit Depth' trace claim of e8.pdf is tested separately in bitdepth_trace_test.py.  Prints ALL CHECKS PASSED."""
import math
from fractions import Fraction
import mpmath

mpmath.mp.dps = 40
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


N = 40


def sigma(k, n):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)


def eis(coef, k):
    return [1] + [coef * sigma(k - 1, n) for n in range(1, N + 1)]


def mul(a, b):
    c = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                c[i + j] += x * b[j]
    return c


E4, E6, E8, E10, E14 = eis(240, 4), eis(-504, 6), eis(480, 8), eis(-264, 10), eis(-24, 14)
prod = [1] + [0] * N
for n in range(1, N + 1):                       # prod (1 - q^n)^24
    for _ in range(24):
        prod = [prod[i] - (prod[i - n] if i >= n else 0) for i in range(N + 1)]
Delta = [0] + prod[:N]                          # q * prod
E4c = mul(mul(E4, E4), E4); E6s = mul(E6, E6)
print('(D) modular identities')
check(all((E4c[i] - E6s[i]) == 1728 * Delta[i] for i in range(N + 1)), 'Delta = (E4^3 - E6^2)/1728 = q prod (1 - q^n)^24 to q^40 (tau(1..6) = %s)' % Delta[1:7])
check(mul(E4, E4) == E8 and mul(E4, E6) == E10 and mul(mul(E4, E4), E6) == E14, 'E8 = E4^2, E10 = E4 E6, E14 = E4^2 E6 to q^40')
check(248 == 8 + 240 == 8 + 3 * 80 and E4[1] == 240, 'dim E8 = 248 = rank 8 + 240 roots (240 = 3 * 80 = the q-coefficient of E4 = the E8 kissing number)')

print('(N) the extremal even unimodular theta series in rank 48, and the Leech theta series')
E4p = [1] + [0] * N
for _ in range(6):
    E4p = mul(E4p, E4)
E43D = mul(E4c, Delta)
D2 = mul(Delta, Delta)
# theta = E4^6 + x E4^3 Delta + y Delta^2 with a_1 = a_2 = 0
x = -Fraction(E4p[1], E43D[1])
y = -Fraction(E4p[2] + x * E43D[2], D2[2])
theta48 = [E4p[i] + x * E43D[i] + y * D2[i] for i in range(N + 1)]
check(x == -1440 and y == 125280, f'a_1 = a_2 = 0 forces theta = E4^6 {x} E4^3 Delta + {y} Delta^2')
check(theta48[1:6] == [0, 0, 52416000, 39007332000, 6609020221440], f'shell sizes 0, 0, {theta48[3]}, {theta48[4]}, {theta48[5]} (as in the PDF)')
check(sum(theta48[3:6]) == 6648079969440, 'sum of the first three shells = 6,648,079,969,440 (the PDF facet count)')
dens = mpmath.mpf(3) * mpmath.pi / 2
dens = dens ** 24 / mpmath.factorial(24)
check(abs(dens - mpmath.mpf('2.31782953e-8')) < 1e-15, f'density (3 pi/2)^24/24! = {mpmath.nstr(dens, 10)}')
leech = [E4c[i] - 720 * Delta[i] for i in range(N + 1)]
check(leech[:5] == [1, 0, 196560, 16773120, 398034000], f'Leech theta E4^3 - 720 Delta: {leech[:5]}')
check(all(691 * leech[n] == 65520 * (sigma(11, n) - Delta[n]) for n in range(1, N + 1)), 'Leech shells a(n) = (65520/691)(sigma_11(n) - tau(n)) for n <= 40')

print('(M) the McLaughlin slice')
forms = [(a, b, c) for a in range(1, 10) for c in range(a, 20) for b in range(-a, a + 1) if 4 * a * c - b * b == 15]
reduced_min4 = [f for f in forms if 2 * f[0] >= 4]
check(sorted(set((2 * a, abs(b), 2 * c) for a, b, c in reduced_min4)) == [(4, 1, 4)], f'reduced even forms [2a, b, 2c] of determinant 15: {forms}; with minimum >= 4 only [4, +-1, 4]')
# norm form of p_2 = (2, omega), omega = (1 + sqrt(-15))/2: N(2a + b omega) = 4a^2 + 2ab + 4b^2
check(all(4 * a * a + 2 * a * b + 4 * b * b == (2 * a + b / 2) ** 2 + 15 * (b / 2) ** 2 for a in range(-4, 5) for b in range(-4, 5)), 'N(2a + b omega) = 4a^2 + 2ab + 4b^2 (Gram [[4,1],[1,4]], det 15)')
# Proposition 6.1: reduced forms (a, b, c) with a = c = 2, b in [-2, 2], discriminant b^2 - 16 fundamental
def fundamental(D):
    if D % 4 == 1:
        return all(D % (p * p) for p in range(3, int(abs(D) ** 0.5) + 2, 2))
    if D % 4 == 0:
        m = D // 4
        return m % 4 in (2, 3) and all(m % (p * p) for p in range(3, int(abs(m) ** 0.5) + 2, 2))
    return False
props = [(b, b * b - 16, fundamental(b * b - 16)) for b in range(-2, 3)]
check([D for b, D, f in props if f] == [-15, -15], f'reduced forms with a = c = 2: (b, D, fundamental) = {props}; only D = -15 is fundamental (Prop. 6.1, elementary)')

print('(E) Epstein zeta values')
z_closed = 240 * mpmath.mpf(2) ** -5 * mpmath.zeta(5) * mpmath.zeta(2)
check(abs(z_closed - mpmath.mpf(5) / 4 * mpmath.pi ** 2 * mpmath.zeta(5)) < mpmath.mpf(10) ** -35, 'zeta_E8(10) = 240 * 2^-5 zeta(5) zeta(2) = (5/4) pi^2 zeta(5) (exact)')
# direct sum over shells: sum_n 240 sigma_3(n) (2n)^-5, tail estimated by the mean value sigma_3(n) ~ zeta(4) n^3
M = 200000
sig3 = [0] * (M + 1)
for d in range(1, M + 1):
    d3 = d ** 3
    for m in range(d, M + 1, d):
        sig3[m] += d3
part = math.fsum(240 * sig3[n] / (2.0 * n) ** 5 for n in range(1, M + 1))
tail = 240 / 32 * float(mpmath.zeta(4)) / M
check(abs(part + tail - float(z_closed)) < 1e-8, f'direct shell sum {part + tail:.10f} (N = 2e5, with tail) = closed form {float(z_closed):.10f}')
taus = Delta  # tau(n) for n <= 40
def L_delta(s):
    s = mpmath.mpf(s)
    tot = mpmath.mpf(0)
    for n in range(1, N + 1):
        x = 2 * mpmath.pi * n
        tot += taus[n] * (mpmath.gammainc(s, x) / x ** s + mpmath.gammainc(12 - s, x) / x ** (12 - s))
    return (2 * mpmath.pi) ** s / mpmath.gamma(s) * tot
# tau(n) for n <= 3000 from the q-product, for the absolutely convergent check at s = 10
pb = [0] * 3001; pb[0] = 1
for n in range(1, 3001):
    for _ in range(24):
        for i in range(3000, n - 1, -1):
            pb[i] -= pb[i - n]
tau3000 = [0] + pb[:3000]
direct10 = mpmath.fsum(mpmath.mpf(tau3000[n]) / mpmath.mpf(n) ** 10 for n in range(1, 3001))
check(abs(L_delta(10) - direct10) < mpmath.mpf(10) ** -10, f'L(Delta, 10) by the incomplete-gamma formula = direct Dirichlet sum = {mpmath.nstr(direct10, 15)}')
L5 = L_delta(5)
zl = -mpmath.mpf(65520) / 691 * mpmath.mpf(2) ** -5 * L5
check(abs(mpmath.zeta(-6)) < mpmath.mpf(10) ** -30 and L5 > 0, f'L(Delta, 5) = {mpmath.nstr(L5, 20)}; zeta(-6) = 0, so zeta_Leech(10) = -(65520/691) 2^-5 L(Delta, 5) = {mpmath.nstr(zl, 20)} (analytic continuation)')
# inside the convergent range the 'odd zeta' Leech value carries an irremovable cusp term: s = 26
L13 = L_delta(13)
direct13 = mpmath.fsum(mpmath.mpf(tau3000[n]) / mpmath.mpf(n) ** 13 for n in range(1, 3001))
z26 = mpmath.mpf(65520) / 691 * mpmath.mpf(2) ** -13 * (mpmath.zeta(13) * mpmath.zeta(2) - L13)
# direct shell sum sum_n a(n) (2n)^-13 with a(n) = (65520/691)(sigma_11(n) - tau(n)), plus the mean-value tail
M2 = 3000
sig11 = [0] * (M2 + 1)
for d in range(1, M2 + 1):
    for m in range(d, M2 + 1, d):
        sig11[m] += d ** 11
shells = mpmath.fsum(mpmath.mpf(65520 * (sig11[k] - tau3000[k])) / 691 / (2 * mpmath.mpf(k)) ** 13 for k in range(1, M2 + 1))
tail26 = mpmath.mpf(65520) / 691 * mpmath.mpf(2) ** -13 * mpmath.zeta(12) / M2   # mean sigma_11(n) = zeta(12) n^11; sum_{n > M} n^-2 ~ 1/M
check(abs(L13 - direct13) < mpmath.mpf(10) ** -12 and abs(shells + tail26 - z26) < 1e-5 * z26,
      f'zeta_Leech(26) = (65520/691) 2^-13 [zeta(13) zeta(2) - L(Delta, 13)] = {mpmath.nstr(z26, 15)} (shell sum {mpmath.nstr(shells + tail26, 10)}); L(Delta, 13) = {mpmath.nstr(L13, 15)} is not removable')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
