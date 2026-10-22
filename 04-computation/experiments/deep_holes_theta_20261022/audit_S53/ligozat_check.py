"""Audit check (S53): Ligozat orders and q-expansions of the eta quotients used for the Monster class names.

For an eta quotient prod_{delta | N} eta(delta tau)^{r_delta} of weight 0, Ligozat's formula gives the order at the
cusp c/d (d | N) in the local parameter of that cusp:
    v_d = (N/24) sum_delta gcd(d, delta)^2 r_delta / (gcd(d, N/d) d delta).
It is a Gamma_0(N)-function iff sum r = 0, sum delta r = 0 mod 24, sum (N/delta) r = 0 mod 24 and prod delta^r is a
rational square. W_N-invariance up to a constant: r_delta = r_{N/delta}. We print the orders and the first
coefficients of E - (constant), to compare with the twining series of THM-4759 (7), (10) and THM-4760 (4)."""
from fractions import Fraction
from math import gcd, isqrt

Q = {
    # THM-4760 (4)
    '4A  (4+)': (4, {1: -24, 2: 48, 4: -24}),
    '8A  (8+)': (8, {1: -8, 2: 8, 4: 8, 8: -8}),
    '12H (12+12)': (12, {1: -4, 3: 4, 4: 4, 12: -4}),
    '16C (16+)': (16, {1: -2, 2: -3, 4: 10, 8: -3, 16: -2}),
    '20F (20+20)': (20, {1: -2, 4: 2, 5: 2, 20: -2}),
    '24I (24+24)': (24, {1: -1, 2: -2, 3: -1, 4: 4, 6: 4, 8: -1, 12: -2, 24: -1}),
    '32A (32+)': (32, {1: -2, 2: 3, 4: -1, 8: -1, 16: 3, 32: -2}),
    # THM-4759 (7): squarefree composite Fricke levels
    '6B  (6+6)': (6, {1: -12, 2: 12, 3: 12, 6: -12}),
    '10D (10+10)': (10, {1: -6, 2: 6, 5: 6, 10: -6}),
    '15C (15+15)': (15, {1: -3, 3: 3, 5: 3, 15: -3}),
    '14C (14+14)': (14, {1: -4, 2: 4, 7: 4, 14: -4}),
    '26B (26+26)': (26, {1: -2, 2: 2, 13: 2, 26: -2}),
}
TARGET = {   # Z - 24 coefficients of q^1, q^2 (, q^3) claimed in the audited files
    '4A  (4+)': [276, 2048, 11202], '8A  (8+)': [36, 128, 386], '12H (12+12)': [14, 36], '16C (16+)': [8, 16],
    '20F (20+20)': [5, 10], '24I (24+24)': [4, 6], '32A (32+)': [2, 4], '6B  (6+6)': [78, 364, 1365],
    '10D (10+10)': [21, 62, 162], '15C (15+15)': [9, 19, 42], '14C (14+14)': [10, 24, 51], '26B (26+26)': [3, 6, 9],
}


def eta_series(r, L):
    """coefficients of prod_delta prod_n (1 - q^{delta n})^{r_delta} up to q^L."""
    c = [0] * (L + 1); c[0] = 1
    for delta, e in r.items():
        for n_ in range(1, L // delta + 1):
            step = delta * n_
            if e > 0:
                for _ in range(e):
                    for k in range(L, step - 1, -1):
                        c[k] -= c[k - step]
            else:
                for _ in range(-e):
                    for k in range(step, L + 1):
                        c[k] += c[k - step]
    return c


def issquare(fr):
    return fr >= 0 and isqrt(fr.numerator) ** 2 == fr.numerator and isqrt(fr.denominator) ** 2 == fr.denominator


allok = True
for name, (N, r) in Q.items():
    s0 = sum(r.values())
    s1 = sum(d * e for d, e in r.items())
    s2 = sum((N // d) * e for d, e in r.items())
    pr = Fraction(1)
    for d, e in r.items():
        pr *= Fraction(d) ** e
    lig = s0 == 0 and s1 % 24 == 0 and s2 % 24 == 0 and issquare(pr)
    orders = {}
    for d in [d for d in range(1, N + 1) if N % d == 0]:
        v = Fraction(N, 24) * sum(Fraction(gcd(d, dl) ** 2 * e, gcd(d, N // d) * d * dl) for dl, e in r.items())
        orders[d] = v
    sym = all(r.get(d, 0) == r.get(N // d, 0) for d in r)
    # q-expansion: q^{s1/24} prod(...)
    lead = Fraction(s1, 24)
    c = eta_series(r, 5)
    const = c[1]          # E = q^-1 (c0 + c1 q + c2 q^2 + ...): constant term is c1
    ser = c[2:5]
    ok = (lig and lead == -1 and orders[N] == -1 and orders[1] == -1 and all(v >= 0 for d, v in orders.items() if d not in (1, N))
          and sym and ser[:len(TARGET[name])] == TARGET[name][:3])
    allok &= ok
    print(f'{"ok  " if ok else "FAIL"} {name:12s} N={N:2d}: Ligozat {lig}; W_N-symmetric {sym}; orders at 1/d: '
          f'{ {d: str(v) for d, v in orders.items()} }; E = q^-1 + {const} + {ser} q..; claimed {TARGET[name]}')
print('ALL OK' if allok else 'SOME FAIL')
