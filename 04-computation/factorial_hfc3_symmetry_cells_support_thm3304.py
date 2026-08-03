"""Exact/modular probes for omitted HFC(3) simplex symmetry cells.

The C3 calculation is homogeneous in barycentric variables.  All degree-dm
simplex moments have the same harmless denominator, so ``weighted_sum``
computes the numerator sum coeff(alpha) * alpha! modulo p.

The S3-sign calculation is over Q and treats the first new projective line,
degree five: Vandermonde * (A e1^2 + B e2).
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import factorial
from time import perf_counter

import sympy as sp


Exp = tuple[int, int, int]


def add(*polys, modulus=None):
    out = {}
    for poly in polys:
        for key, value in poly.items():
            value = value if modulus is None else value % modulus
            out[key] = out.get(key, 0) + value
            if modulus is not None:
                out[key] %= modulus
    return {key: value for key, value in out.items() if value}


def scale(poly, scalar, modulus=None):
    if modulus is None:
        return {key: value * scalar for key, value in poly.items()
                if value * scalar}
    return {key: value * scalar % modulus for key, value in poly.items()
            if value * scalar % modulus}


def mul(left, right, modulus=None):
    out = {}
    for a, av in left.items():
        for b, bv in right.items():
            key = (a[0] + b[0], a[1] + b[1], a[2] + b[2])
            out[key] = out.get(key, 0) + av * bv
            if modulus is not None:
                out[key] %= modulus
    return {key: value for key, value in out.items() if value}


def power(poly, exponent, modulus=None):
    result = {(0, 0, 0): 1}
    base = poly
    while exponent:
        if exponent & 1:
            result = mul(result, base, modulus)
        exponent >>= 1
        if exponent:
            base = mul(base, base, modulus)
    return result


def weak_compositions(total, length, prefix=()):
    if length == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from weak_compositions(total - value, length - 1,
                                     prefix + (value,))


def multinomial(parts):
    answer = factorial(sum(parts))
    for part in parts:
        answer //= factorial(part)
    return answer


def weighted_sum(poly, modulus=None):
    answer = 0
    for alpha, coefficient in poly.items():
        term = coefficient
        for exponent in alpha:
            term *= factorial(exponent)
        answer += term
        if modulus is not None:
            answer %= modulus
    return answer if modulus is None else answer % modulus


def moment_form(basis, exponent, modulus=None):
    """Coefficient dictionary of the moment form in basis coordinates."""
    powers = []
    for item in basis:
        powers.append([power(item, e, modulus) for e in range(exponent + 1)])
    out = {}
    for counts in weak_compositions(exponent, len(basis)):
        poly = {(0, 0, 0): 1}
        for index, count in enumerate(counts):
            poly = mul(poly, powers[index][count], modulus)
        coefficient = multinomial(counts) * weighted_sum(poly, modulus)
        if modulus is not None:
            coefficient %= modulus
        if coefficient:
            out[counts] = coefficient
    return out


def form_to_sympy(form, variables, modulus=None):
    if modulus is not None:
        return sp.Poly.from_dict(form, variables, modulus=modulus)
    return sp.Poly.from_dict(form, variables, domain=sp.QQ)


def primitive_root_of_unity_3(prime):
    for candidate in range(2, prime):
        root = pow(candidate, (prime - 1) // 3, prime)
        if root != 1 and (root * root + root + 1) % prime == 0:
            return root
    raise AssertionError("no cube root")


def cyclic_degree_four_basis(prime):
    """omega-eigen homogeneous degree-four basis on Delta_2.

    a=s1+omega*s2+omega^2*s3 has weight omega^2 and
    b=s1+omega^2*s2+omega*s3 has weight omega.  Multiplication by e1
    homogenizes without changing the restricted function.
    """
    omega = primitive_root_of_unity_3(prime)
    a = {(1, 0, 0): 1, (0, 1, 0): omega,
         (0, 0, 1): omega * omega % prime}
    b = {(1, 0, 0): 1, (0, 1, 0): omega * omega % prime,
         (0, 0, 1): omega}
    e1 = {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
    basis = [
        mul(b, power(e1, 3, prime), prime),
        mul(power(a, 2, prime), power(e1, 2, prime), prime),
        mul(mul(a, power(b, 2, prime), prime), e1, prime),
        power(b, 4, prime),
        mul(power(a, 3, prime), b, prime),
    ]
    return omega, basis


def vandermonde():
    x = {(1, 0, 0): Fraction(1)}
    y = {(0, 1, 0): Fraction(1)}
    z = {(0, 0, 1): Fraction(1)}
    return mul(mul(add(x, scale(y, -1)), add(y, scale(z, -1))),
               add(z, scale(x, -1)))


def sign_degree_five_probe():
    e1 = {(1, 0, 0): Fraction(1), (0, 1, 0): Fraction(1),
          (0, 0, 1): Fraction(1)}
    e2 = {(1, 1, 0): Fraction(1), (1, 0, 1): Fraction(1),
          (0, 1, 1): Fraction(1)}
    delta = vandermonde()
    basis = [mul(delta, power(e1, 2)), mul(delta, e2)]
    forms = {m: moment_form(basis, m) for m in (2, 4)}
    t = sp.symbols("t")
    polys = {}
    for m, form in forms.items():
        coeffs = [form.get((m - j, j), 0) for j in range(m + 1)]
        expression = sum(value * t ** j for j, value in enumerate(coeffs))
        polys[m] = sp.Poly(expression, t, domain=sp.QQ)
    gcd = sp.gcd(polys[2], polys[4])
    # The A=0 projective point is checked separately by its top coefficient.
    at_infinity = {m: forms[m].get((0, m), 0) for m in forms}
    return polys, gcd, at_infinity


def cyclic_degree_four_forms(prime=103, exponents=(3, 6, 9, 12, 15)):
    omega, basis = cyclic_degree_four_basis(prime)
    variables = sp.symbols("c0:5")
    forms = {}
    timings = {}
    for exponent in exponents:
        started = perf_counter()
        raw = moment_form(basis, exponent, prime)
        forms[exponent] = form_to_sympy(raw, variables, prime)
        timings[exponent] = perf_counter() - started
    return omega, variables, forms, timings


def cyclic_degree_four_forms_fast(prime=103, exponents=(3, 6, 9, 12, 15)):
    """Same restricted moment forms, using Fourier-monomial coordinates.

    On e1=1 the five basis functions are b, a^2, a*b^2, b^4,
    a^3*b.  A parameter monomial therefore needs only one cached Dirichlet
    moment of a^R*b^S.  This is much faster than rehomogenizing every product.
    """
    omega = primitive_root_of_unity_3(prime)
    a = {(1, 0, 0): 1, (0, 1, 0): omega,
         (0, 0, 1): omega * omega % prime}
    b = {(1, 0, 0): 1, (0, 1, 0): omega * omega % prime,
         (0, 0, 1): omega}
    exponent_pairs = ((0, 1), (2, 0), (1, 2), (0, 4), (3, 1))
    variables = sp.symbols("c0:5")
    @lru_cache(None)
    def dirichlet_moment(r, s):
        # The exponential/factorial generating function is
        #   product_j (1-u*a_j-v*b_j)^(-1)
        #     = 1/(1-u^3-v^3-3uv).
        # Hence the coefficient is a short trinomial sum.
        coefficient = 0
        for ell in range(min(r, s) + 1):
            if (r - ell) % 3 or (s - ell) % 3:
                continue
            i = (r - ell) // 3
            j = (s - ell) // 3
            coefficient += multinomial((i, j, ell)) * 3 ** ell
        numerator = 2 * factorial(r) * factorial(s) * coefficient % prime
        denominator = factorial(r + s + 2) % prime
        return numerator * pow(denominator, prime - 2, prime) % prime

    forms = {}
    timings = {}
    for exponent in exponents:
        started = perf_counter()
        raw = {}
        for counts in weak_compositions(exponent, 5):
            r = sum(count * pair[0]
                    for count, pair in zip(counts, exponent_pairs))
            s = sum(count * pair[1]
                    for count, pair in zip(counts, exponent_pairs))
            coefficient = multinomial(counts) * dirichlet_moment(r, s)
            coefficient %= prime
            if coefficient:
                raw[counts] = coefficient
        forms[exponent] = form_to_sympy(raw, variables, prime)
        timings[exponent] = perf_counter() - started
    return omega, variables, forms, timings


def chart_groebner(forms, variables, chart, used_exponents):
    remaining = variables[:chart] + variables[chart + 1:]
    equations = [forms[e].as_expr().subs(variables[chart], 1)
                 for e in used_exponents]
    started = perf_counter()
    basis = sp.groebner(equations, *remaining, modulus=forms[used_exponents[0]].get_modulus(),
                        order="grevlex")
    elapsed = perf_counter() - started
    contains_one = any(poly.as_expr() == 1 for poly in basis.polys)
    return contains_one, len(basis.polys), elapsed, basis


if __name__ == "__main__":
    print("THM-3304 support module; run the two declared companion scripts")
