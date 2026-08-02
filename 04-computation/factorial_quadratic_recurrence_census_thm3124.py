#!/usr/bin/env python3
"""Exact controls for THM-3124's quadratic factorial-moment recurrence."""

from math import factorial, isqrt

import sympy as sp


MAX_R = 200
PRIMES = (1_000_003, 1_000_033)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True


def trim(f):
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


def poly_divmod(f, g, p):
    f = trim(f[:])
    g = trim(g[:])
    if g == [0]:
        raise ZeroDivisionError
    inv_lead = pow(g[-1], p - 2, p)
    while len(f) >= len(g) and f != [0]:
        shift = len(f) - len(g)
        scale = f[-1] * inv_lead % p
        for j, gj in enumerate(g):
            f[j + shift] = (f[j + shift] - scale * gj) % p
        trim(f)
    return f


def poly_gcd(f, g, p):
    f, g = trim(f[:]), trim(g[:])
    while g != [0]:
        f, g = g, poly_divmod(f, g, p)
    inv_lead = pow(f[-1], p - 2, p)
    return [(x * inv_lead) % p for x in f]


def tables(p, max_n):
    choose = [[0] * (max_n + 1) for _ in range(max_n + 1)]
    for n in range(max_n + 1):
        choose[n][0] = choose[n][n] = 1
        for k in range(1, n):
            choose[n][k] = (choose[n - 1][k - 1] + choose[n - 1][k]) % p
    facts = [1] * (2 * max_n + 1)
    for n in range(1, len(facts)):
        facts[n] = facts[n - 1] * n % p
    return choose, facts


def moment_poly_mod(n, r, p, choose, facts):
    """Coefficients low-to-high in u for q=1-t/(r+2)+u*t^2."""
    b = -pow(r + 2, p - 2, p) % p
    bpowers = [1] * (n + 1)
    for ell in range(1, n + 1):
        bpowers[ell] = bpowers[ell - 1] * b % p
    coeffs = []
    for j in range(n + 1):
        inner = 0
        for ell in range(n - j + 1):
            inner += choose[n - j][ell] * bpowers[ell] * facts[2 * j + ell]
        coeffs.append(choose[n][j] * (inner % p) % p)
    return trim(coeffs)


def modular_census(p):
    choose, facts = tables(p, MAX_R + 1)
    flags = []
    preserved = True
    for r in range(1, MAX_R + 1):
        f = moment_poly_mod(r, r, p, choose, facts)
        g = moment_poly_mod(r + 1, r, p, choose, facts)
        preserved &= len(f) == r + 1 and f[-1] == facts[2 * r]
        preserved &= len(g) == r + 2 and g[-1] == facts[2 * r + 2]
        h = poly_gcd(f, g, p)
        if len(h) > 1:
            flags.append((r, len(h) - 1))
    return preserved, flags


def exact_moment(poly, t):
    expanded = sp.Poly(sp.expand(poly), t)
    return sp.expand(sum(coef * factorial(exp[0]) for exp, coef in expanded.terms()))


def verify_recurrence():
    t, a, b, c = sp.symbols("t a b c")
    q = a + b * t + c * t**2
    moments = [sp.Integer(1)] + [exact_moment(q**n, t) for n in range(1, 8)]
    checks = []
    for n in range(1, 7):
        rhs = (
            a**n * (a + (n + 1) * b)
            + 2 * (n + 1) * (2 * n + 1) * c * moments[n]
            + n * (n + 1) * (b**2 - 4 * a * c) * moments[n - 1]
        )
        checks.append(sp.expand(moments[n + 1] - rhs) == 0)
    return checks


def rational_gcd_degrees():
    c, t = sp.symbols("c t")
    answer = []
    for r in range(1, 8):
        b = sp.Rational(-1, r + 2)
        q = 1 + b * t + c * t**2
        mr = sp.Poly(exact_moment(q**r, t), c, domain=sp.QQ)
        mr1 = sp.Poly(exact_moment(q ** (r + 1), t), c, domain=sp.QQ)
        answer.append((r, sp.gcd(mr, mr1).degree()))
    return answer


def main():
    print("quadratic_recurrence_symbolic_n_1_through_6=", verify_recurrence(), sep="")
    print("rational_resonance_gcd_degrees_r_1_through_7=", rational_gcd_degrees(), sep="")
    for p in PRIMES:
        preserved, flags = modular_census(p)
        print(
            f"prime={p} isprime={is_prime(p)} candidate_windows=1..{MAX_R} "
            f"degree_preserved={preserved} common_factor_flags={flags}"
        )
    for p in PRIMES:
        positive = poly_gcd([1, 2, 1], [2, 3, 1], p)
        print(f"positive_control_prime={p} gcd_low_to_high={positive}")


if __name__ == "__main__":
    main()
