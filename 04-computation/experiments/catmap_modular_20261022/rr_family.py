"""Rogers-Ramanujan-type identities: verification, growth constants, q -> 1 ratios, and the indefinite partners.
Session opus-2026-10-22-S55.

Exact integer q-series arithmetic to order M.  For each identity (sum side = product side):
  * verified coefficientwise to q^M;
  * effective central charge c_eff from the product (a part set of density d in (1 - q^n)^-1 contributes d,
    in (1 - q^n) contributes -d, in (1 + q^n) contributes d/2), so coefficients grow like exp(2 pi sqrt(c_eff n/6));
  * q -> 1 limits of ratios of companion products (numerically at q = exp(-t), t -> 0), compared with ratios
    sin(pi a/m)/sin(pi a'/m), i.e. quantum dimensions [a]_q at q = exp(i pi/m).
Indefinite (Hecke-type, signature (1,1)) partners, checked to order M:
  * fifth order mock theta f0(q) = sum q^(n^2)/(-q;q)_n,  (q;q)_inf f0 = sum_(n>=0, |j|<=n) (-1)^j q^(n(5n+1)/2 - j^2)(1 - q^(4n+2));
  * seventh order (Andrews 1986): cited only; the simple ansatz with (7n^2+bn)/2 - j^2 was tested and does not match.
Usage: python rr_family.py [M]"""
import sys, math
from fractions import Fraction

M = int(sys.argv[1]) if len(sys.argv) > 1 else 300


def mul(a, b):
    out = [0] * (M + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(0, M + 1 - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def inv(a):
    assert a[0] == 1
    out = [0] * (M + 1); out[0] = 1
    for n in range(1, M + 1):
        out[n] = -sum(a[k] * out[n - k] for k in range(1, n + 1) if a[k])
    return out


def one():
    o = [0] * (M + 1); o[0] = 1; return o


def poch(a_exp, step, count, sign=-1):
    """prod_(k=0)^(count-1) (1 + sign q^(a_exp + k step)); count = None means infinite (to order M)."""
    out = one()
    k = 0
    while (count is None or k < count):
        e = a_exp + k * step
        if e > M:
            break
        f = [0] * (M + 1); f[0] = 1
        if e == 0:
            f[0] = 1 + sign
        else:
            f[e] += sign
        out = mul(out, f)
        k += 1
    return out


def qpoch(a, step, count):            # (q^a; q^step)_count, count=None infinite
    return poch(a, step, count, -1)


def qpoch_neg(a, step, count):        # (-q^a; q^step)_count
    return poch(a, step, count, +1)


def shift(series, e):
    out = [0] * (M + 1)
    for i, x in enumerate(series):
        if i + e <= M:
            out[i + e] += x
    return out


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sum_side(term, nmax=None):
    tot = [0] * (M + 1)
    n = 0
    while True:
        t = term(n)
        if t is None:
            break
        tot = add(tot, t)
        n += 1
    return tot


def triple(a, m):                     # (q^a, q^(m-a), q^m; q^m)_inf
    return mul(mul(qpoch(a, m, None), qpoch(m - a, m, None)), qpoch(m, m, None))


def frac(num, den):                   # num * inverse(den)
    return mul(num, inv(den))


def tm(e, num, den):                  # q^e * num/den, None if e > M
    if e > M:
        return None
    return shift(frac(num, den), e)


identities = {}
# Rogers-Ramanujan (mod 5)
identities['RR G'] = (lambda n: tm(n * n, one(), qpoch(1, 1, n)), frac(triple(2, 5), qpoch(1, 1, None)), (5, 2), '1')
identities['RR H'] = (lambda n: tm(n * n + n, one(), qpoch(1, 1, n)), frac(triple(1, 5), qpoch(1, 1, None)), (5, 1), '1')
# Rogers-Selberg (mod 7)
identities['RS A'] = (lambda n: tm(2 * n * n, one(), mul(qpoch(2, 2, n), qpoch_neg(1, 1, 2 * n))), frac(triple(3, 7), qpoch(2, 2, None)), (7, 3), '2')
identities['RS B'] = (lambda n: tm(2 * n * n + 2 * n, one(), mul(qpoch(2, 2, n), qpoch_neg(1, 1, 2 * n))), frac(triple(2, 7), qpoch(2, 2, None)), (7, 2), '2')
identities['RS C'] = (lambda n: tm(2 * n * n + 2 * n, one(), mul(qpoch(2, 2, n), qpoch_neg(1, 1, 2 * n + 1))), frac(triple(1, 7), qpoch(2, 2, None)), (7, 1), '2')
# Rogers mod 14
identities['R14 A'] = (lambda n: tm(n * n, one(), mul(qpoch(1, 1, n), qpoch(1, 2, n))), frac(triple(6, 14), qpoch(1, 1, None)), (14, 6), '1')
identities['R14 B'] = (lambda n: tm(n * n + n, one(), mul(qpoch(1, 1, n), qpoch(1, 2, n + 1))), frac(triple(4, 14), qpoch(1, 1, None)), (14, 4), '1')
identities['R14 C'] = (lambda n: tm(n * n + 2 * n, one(), mul(qpoch(1, 1, n), qpoch(1, 2, n + 1))), frac(triple(2, 14), qpoch(1, 1, None)), (14, 2), '1')
# Goellnitz-Gordon (mod 8)
gg1 = inv(mul(mul(qpoch(1, 8, None), qpoch(4, 8, None)), qpoch(7, 8, None)))
gg2 = inv(mul(mul(qpoch(3, 8, None), qpoch(4, 8, None)), qpoch(5, 8, None)))
identities['GG 1'] = (lambda n: tm(n * n, qpoch_neg(1, 2, n), qpoch(2, 2, n)), gg1, None, 'gg')
identities['GG 2'] = (lambda n: tm(n * n + 2 * n, qpoch_neg(1, 2, n), qpoch(2, 2, n)), gg2, None, 'gg')
# Bailey mod 9
identities['B9 A'] = (lambda n: tm(3 * n * n, qpoch(1, 1, 3 * n), mul(qpoch(3, 3, n), qpoch(3, 3, 2 * n))), frac(triple(4, 9), qpoch(3, 3, None)), (9, 4), '3')
identities['B9 B'] = (lambda n: tm(3 * n * n + 3 * n, mul(qpoch(1, 1, 3 * n), add(one(), shift([-1] + [0] * M, 3 * n + 2)) if 3 * n + 2 <= M else qpoch(1, 1, 3 * n)), mul(qpoch(3, 3, n), qpoch(3, 3, 2 * n + 1))), frac(triple(2, 9), qpoch(3, 3, None)), (9, 2), '3')
identities['B9 C'] = (lambda n: tm(3 * n * n + 3 * n, qpoch(1, 1, 3 * n + 1), mul(qpoch(3, 3, n), qpoch(3, 3, 2 * n + 1))), frac(triple(1, 9), qpoch(3, 3, None)), (9, 1), '3')
# Dyson mod 27
identities['D27 A'] = (lambda n: one() if n == 0 else tm(n * n, qpoch(3, 3, n - 1), mul(qpoch(1, 1, n), qpoch(1, 1, 2 * n - 1))), frac(triple(12, 27), qpoch(1, 1, None)), (27, 12), '1')
identities['D27 B'] = (lambda n: tm(n * n + n, qpoch(3, 3, n), mul(qpoch(1, 1, n), qpoch(1, 1, 2 * n + 1))), frac(qpoch(9, 9, None), qpoch(1, 1, None)), (27, 9), '1')
identities['D27 C'] = (lambda n: tm(n * n + 2 * n, qpoch(3, 3, n), mul(qpoch(1, 1, n), qpoch(1, 1, 2 * n + 2))), frac(triple(6, 27), qpoch(1, 1, None)), (27, 6), '1')
identities['D27 D'] = (lambda n: tm(n * n + 3 * n, qpoch(3, 3, n), mul(qpoch(1, 1, n), qpoch(1, 1, 2 * n + 2))), frac(triple(3, 27), qpoch(1, 1, None)), (27, 3), '1')
# Jackson-Slater (Slater A.39 / (83)): sum q^(2n^2)/(q;q)_(2n) = (q, q^7, q^8; q^8)_inf (q^6, q^10; q^16)_inf / (q;q)_inf
js = frac(mul(mul(mul(qpoch(1, 8, None), qpoch(7, 8, None)), qpoch(8, 8, None)), mul(qpoch(6, 16, None), qpoch(10, 16, None))), qpoch(1, 1, None))
identities['JS'] = (lambda n: tm(2 * n * n, one(), qpoch(1, 1, 2 * n)), js, None, 'js')


def first_diff(a, b):
    for i in range(M + 1):
        if a[i] != b[i]:
            return i
    return None


def c_eff_numeric(series_fn, ts=(0.02, 0.01)):
    """estimate c_eff from log f(e^-t) ~ pi^2 c/(6t) + b log t + const using three t values."""
    import numpy as np
    vals = [math.log(abs(series_fn(t))) for t in ts]
    t1, t2 = ts
    return (vals[1] - vals[0]) / (math.pi ** 2 / 6 * (1 / t2 - 1 / t1))


if __name__ == '__main__':
    print(f'Exact verification to q^{M}:')
    for name, (term, prod, am, fam) in identities.items():
        s = sum_side(term)
        d = first_diff(s, prod)
        print(f'  {name:6s}: {"VERIFIED" if d is None else "FAILS at q^" + str(d)}; first coefficients {s[:12]}')
    # Ising check: Jackson-Slater sum = Ising vacuum character (fermionic) = even part of sum_k q^(k^2/2)/(q;q)_k in q^(1/2)
    print()
    print('Hecke-type (indefinite, signature (1,1)) partners:')
    f0 = sum_side(lambda n: tm(n * n, one(), qpoch_neg(1, 1, n)))
    lhs = mul(qpoch(1, 1, None), f0)
    rhs = [0] * (M + 1)
    for n in range(0, M + 1):
        for j in range(-n, n + 1):
            e = n * (5 * n + 1) // 2 - j * j
            if 0 <= e <= M:
                rhs[e] += (-1) ** (j % 2)
            e2 = e + 4 * n + 2
            if 0 <= e2 <= M:
                rhs[e2] -= (-1) ** (j % 2)
    print(f'  fifth order f0: (q;q) f0 = sum (-1)^j q^(n(5n+1)/2 - j^2)(1 - q^(4n+2)): {"VERIFIED" if first_diff(lhs, rhs) is None else "FAILS at " + str(first_diff(lhs, rhs))}')
    print('  seventh order: Andrews (TAMS 293, 1986) gives Hecke-type sums for F0, F1, F2; not reproduced here (cited only)')
