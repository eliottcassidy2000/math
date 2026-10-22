"""Audit G (synthesis note section 4): Rogers-Ramanujan-type identities, c_eff, q -> 1 ratios, the claimed relations,
the fifth-order mock theta Hecke-type identity, and the actual coefficient growth of the signed families.
Auditor for opus-2026-10-22-S55.  Own q-series code (exact integers), written from the literature forms
(Slater 1952 numbering as in McLaughlin-Sills-Zimmer, Electron. J. Combin. DS15).
Usage: python g_rr.py [M]"""
import sys, math, cmath
from fractions import Fraction

M = int(sys.argv[1]) if len(sys.argv) > 1 else 250


class PS:
    """truncated power series with integer (or Fraction) coefficients mod q^(M+1)."""
    def __init__(self, c=None):
        self.c = (list(c) + [0] * (M + 1))[:M + 1] if c is not None else [0] * (M + 1)

    @staticmethod
    def one():
        return PS([1])

    def __add__(self, o):
        return PS([a + b for a, b in zip(self.c, o.c)])

    def __sub__(self, o):
        return PS([a - b for a, b in zip(self.c, o.c)])

    def __mul__(self, o):
        if isinstance(o, (int, Fraction)):
            return PS([a * o for a in self.c])
        out = [0] * (M + 1)
        nz = [(j, b) for j, b in enumerate(o.c) if b]
        for i, a in enumerate(self.c):
            if a:
                for j, b in nz:
                    if i + j > M:
                        break
                    out[i + j] += a * b
        return PS(out)
    __rmul__ = __mul__

    def shift(self, k):
        if k > M:
            return PS()
        return PS([0] * k + self.c[:M + 1 - k])

    def inv(self):
        assert self.c[0] in (1, -1)
        out = [0] * (M + 1); out[0] = Fraction(1, self.c[0]) if self.c[0] != 1 else 1
        nz = [(k, a) for k, a in enumerate(self.c) if a and k]
        for n in range(1, M + 1):
            s = 0
            for k, a in nz:
                if k > n:
                    break
                s += a * out[n - k]
            out[n] = -s * out[0]
        return PS(out)

    def __truediv__(self, o):
        return self * o.inv()

    def subs_power(self, k):                      # f(q^k)
        out = [0] * (M + 1)
        for i, a in enumerate(self.c):
            if i * k > M:
                break
            out[i * k] = a
        return PS(out)


def binom1(e, sign):                              # 1 + sign*q^e
    c = [0] * (M + 1); c[0] = 1
    if e <= M:
        c[e] += sign
    return PS(c)


def poch(a, step, count=None, sign=-1):           # (sign... ) prod_{k<count} (1 + sign q^(a+k step))
    out = PS.one(); k = 0
    while count is None or k < count:
        e = a + k * step
        if e > M:
            break
        if e == 0:
            out = out * PS([1 + sign])
        else:
            out = out * binom1(e, sign)
        k += 1
    return out


def qp(a, step, count=None):
    return poch(a, step, count, -1)


def qpn(a, step, count=None):                     # (-q^a; q^step)_count
    return poch(a, step, count, +1)


def jtp(a, m):                                    # (q^a, q^(m-a), q^m; q^m)_inf  via the triple product SERIES (independent of products)
    c = [0] * (M + 1)
    j = 0
    # theta_{a,m} = sum_j (-1)^j q^(m j(j-1)/2 + a j)
    for j in range(-200, 201):
        e = m * j * (j - 1) // 2 + a * j
        if 0 <= e <= M:
            c[e] += (-1) ** (j % 2)
    return PS(c)


def ssum(term, nmax=200):
    tot = PS()
    for n in range(nmax):
        t = term(n)
        if t is None:
            break
        tot = tot + t
    return tot


def T(e, num, den):
    return None if e > M else (num / den).shift(e)


Q1 = qp(1, 1)
ids = {}
ids['RR G'] = (ssum(lambda n: T(n * n, PS.one(), qp(1, 1, n))), (qp(1, 5) * qp(4, 5)).inv())
ids['RR H'] = (ssum(lambda n: T(n * n + n, PS.one(), qp(1, 1, n))), (qp(2, 5) * qp(3, 5)).inv())
ids['RS A'] = (ssum(lambda n: T(2 * n * n, PS.one(), qp(2, 2, n) * qpn(1, 1, 2 * n))), jtp(3, 7) / qp(2, 2))
ids['RS B'] = (ssum(lambda n: T(2 * n * n + 2 * n, PS.one(), qp(2, 2, n) * qpn(1, 1, 2 * n))), jtp(2, 7) / qp(2, 2))
ids['RS C'] = (ssum(lambda n: T(2 * n * n + 2 * n, PS.one(), qp(2, 2, n) * qpn(1, 1, 2 * n + 1))), jtp(1, 7) / qp(2, 2))
ids['R14 (61)'] = (ssum(lambda n: T(n * n, PS.one(), qp(1, 2, n) * qp(1, 1, n))), jtp(6, 14) / Q1)
ids['R14 (60)'] = (ssum(lambda n: T(n * n + n, PS.one(), qp(1, 2, n + 1) * qp(1, 1, n))), jtp(4, 14) / Q1)
ids['R14 (59)'] = (ssum(lambda n: T(n * n + 2 * n, PS.one(), qp(1, 2, n + 1) * qp(1, 1, n))), jtp(2, 14) / Q1)
ids['GG 1'] = (ssum(lambda n: T(n * n, qpn(1, 2, n), qp(2, 2, n))), (qp(1, 8) * qp(4, 8) * qp(7, 8)).inv())
ids['GG 2'] = (ssum(lambda n: T(n * n + 2 * n, qpn(1, 2, n), qp(2, 2, n))), (qp(3, 8) * qp(4, 8) * qp(5, 8)).inv())
ids['B9 A (S.42)'] = (ssum(lambda n: T(3 * n * n, qp(1, 1, 3 * n), qp(3, 3, n) * qp(3, 3, 2 * n))), jtp(4, 9) / qp(3, 3))
ids['B9 B (S.41)'] = (ssum(lambda n: T(3 * n * n + 3 * n, qp(1, 1, 3 * n) * binom1(3 * n + 2, -1), qp(3, 3, n) * qp(3, 3, 2 * n + 1))), jtp(2, 9) / qp(3, 3))
ids['B9 C (S.40)'] = (ssum(lambda n: T(3 * n * n + 3 * n, qp(1, 1, 3 * n + 1), qp(3, 3, n) * qp(3, 3, 2 * n + 1))), jtp(1, 9) / qp(3, 3))
ids['D27 A (S.93)'] = (ssum(lambda n: PS.one() if n == 0 else T(n * n, qp(3, 3, n - 1), qp(1, 1, n) * qp(1, 1, 2 * n - 1))), jtp(12, 27) / Q1)
ids['D27 B (S.92)'] = (ssum(lambda n: T(n * n + n, qp(3, 3, n), qp(1, 1, n) * qp(1, 1, 2 * n + 1))), qp(9, 9) / Q1)
ids['D27 C (S.91)'] = (ssum(lambda n: T(n * n + 2 * n, qp(3, 3, n), qp(1, 1, n) * qp(1, 1, 2 * n + 2))), jtp(6, 27) / Q1)
ids['D27 D (S.90)'] = (ssum(lambda n: T(n * n + 3 * n, qp(3, 3, n), qp(1, 1, n) * qp(1, 1, 2 * n + 2))), jtp(3, 27) / Q1)
JS_sum = ssum(lambda n: T(2 * n * n, PS.one(), qp(1, 1, 2 * n)))
JS_prod_session = qp(1, 8) * qp(7, 8) * qp(8, 8) * qp(6, 16) * qp(10, 16) / Q1
JS_prod_slater39 = qpn(3, 8) * qpn(5, 8) * qp(8, 8) / qp(2, 2)
ids['JS (session product)'] = (JS_sum, JS_prod_session)
ids['JS (Slater 39 product)'] = (JS_sum, JS_prod_slater39)

print(f'=== identities, exact to q^{M} ===')
for k, (a, b) in ids.items():
    d = next((i for i in range(M + 1) if a.c[i] != b.c[i]), None)
    print(f'  {k:24s}: {"VERIFIED" if d is None else "FAILS at q^%d" % d}; coefficients {a.c[:10]}')


# ---------------- Rocha-Caridi characters (independent of the product sides)
def rocha_caridi(p, pp, r, s):
    """q^(c/24 - h) chi_{r,s}(q) as a power series: (1/(q)_inf) sum_k [q^(A_k) - q^(B_k)] with exponents shifted so the
    leading term is q^0.  Convention: h = ((p r - pp s)^2 - (p - pp)^2)/(4 p pp)."""
    c = [0] * (M + 1)
    base = Fraction((p * r - pp * s) ** 2 - (p - pp) ** 2, 4 * p * pp)
    for k in range(-60, 61):
        for sgn, x in ((1, 2 * p * pp * k + p * r - pp * s), (-1, 2 * p * pp * k + p * r + pp * s)):
            ex = Fraction(x * x - (p - pp) ** 2, 4 * p * pp) - base
            assert ex.denominator == 1
            if 0 <= ex <= M:
                c[int(ex)] += sgn
    return PS(c) / Q1


print('\n=== relations claimed in the note ===')
# Jackson-Slater sum = q^(1/48) chi_0^Ising  (Rocha-Caridi M(3,4), (r,s) = (1,1), h = 0, p = 3, p' = 4)
ising0 = rocha_caridi(3, 4, 1, 1)
print('  JS sum = Ising vacuum (Rocha-Caridi M(3,4) chi_{1,1}) :', ising0.c == JS_sum.c)
# and = even part of (-q^(1/2); q)_inf : work in Q = q^(1/2): sum_n Q^(n^2)/(Q^2;Q^2)_n, keep even powers
Mold = M
M = 2 * Mold
QQ = ssum(lambda n: T(n * n, PS.one(), qp(2, 2, n)))
even = PS([QQ.c[2 * i] for i in range(Mold + 1)] + [0] * (M - Mold))
M = Mold
print('  JS sum = even part of (-q^(1/2); q)_inf = sum_n q^(n^2/2)/(q;q)_n          :', even.c[:M + 1] == JS_sum.c)
# M(2,7) characters: Rocha-Caridi (p, p') = (2, 7), s = 1, r = 1, 2, 3 ; compare theta_{a,7}/(q;q)
for r in (1, 2, 3):
    rc = rocha_caridi(2, 7, r, 1)
    hit = [a for a in (1, 2, 3) if (jtp(a, 7) / Q1).c == rc.c]
    print(f'  M(2,7) chi_(r={r},s=1) (Rocha-Caridi) = theta_(a,7)/(q;q) for a = {hit}')
# Rogers-Selberg sum = M(2,7) character x (q;q^2)_inf ; Rogers mod 14 = M(2,7) char at q^2 x (-q;q)_inf
odd = qp(1, 2)
for lab, a, r in (('RS A', 3, 3), ('RS B', 2, 2), ('RS C', 1, 1)):
    rc = rocha_caridi(2, 7, r, 1)
    print(f'  {lab} sum = chi^M(2,7)_(r={r}) x (q;q^2)_inf :', (rc * odd).c == ids[lab][0].c)
mneg = qpn(1, 1)
for lab, r in (('R14 (61)', 3), ('R14 (60)', 2), ('R14 (59)', 1)):
    rc = rocha_caridi(2, 7, r, 1).subs_power(2)
    print(f'  {lab} sum = chi^M(2,7)_(r={r})(q^2) x (-q;q)_inf :', (rc * mneg).c == ids[lab][0].c)
# Dyson_b(q) = Bailey_b(q^3) (q^9;q^9)/(q;q) at the level of SUM sides; and Dyson products = M(2,27) characters, 3 | a
fac = qp(9, 9) / Q1
for dl, bl in (('D27 A (S.93)', 'B9 A (S.42)'), ('D27 C (S.91)', 'B9 B (S.41)'), ('D27 D (S.90)', 'B9 C (S.40)')):
    lhs = ids[dl][0]; rhs = ids[bl][0].subs_power(3) * fac
    print(f'  {dl} sum = {bl} sum (q -> q^3) x (q^9;q^9)/(q;q) :', lhs.c == rhs.c)
print('  D27 B sum = (q^9;q^9)/(q;q) [trivial b = 3]:', ids['D27 B (S.92)'][0].c == fac.c)
for dl, a in (('D27 A (S.93)', 12), ('D27 B (S.92)', 9), ('D27 C (S.91)', 6), ('D27 D (S.90)', 3)):
    # M(2,27): p = 2, p' = 27, s = 1; theta_{a,27} with a = |2r - 27|/... find r with Rocha-Caridi = product
    hit = [r for r in range(1, 14) if rocha_caridi(2, 27, r, 1).c == ids[dl][1].c]
    print(f'  {dl} product = M(2,27) Rocha-Caridi character chi_(r,1) for r = {hit}  (a = {a})')

# ---------------- fifth-order mock theta f0
f0 = ssum(lambda n: T(n * n, PS.one(), qpn(1, 1, n)))
lhs = Q1 * f0
rhs = [0] * (M + 1)
for n in range(0, 60):
    for j in range(-n, n + 1):
        e1 = n * (5 * n + 1) // 2 - j * j
        for e_, sg in ((e1, 1), (e1 + 4 * n + 2, -1)):
            if 0 <= e_ <= M:
                rhs[e_] += sg * (-1) ** (j % 2)
print('\n  fifth order: (q;q)_inf f0(q) = sum_{n>=0,|j|<=n} (-1)^j q^(n(5n+1)/2 - j^2)(1 - q^(4n+2)) :', lhs.c == rhs)
print('  f0 coefficients:', f0.c[:16])
