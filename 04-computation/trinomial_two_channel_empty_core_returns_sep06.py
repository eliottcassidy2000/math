#!/usr/bin/env python3
"""Exact controls for the two-first-channel trinomial return theorem.

No external dependencies. Default support universe: primitive (-a,b,c),
1<=a,b<=40, b<c<=60; filter: exactly two channels at mass gcd(a+b,a+c).
An independent rational Laurent convolution checks named controls.
Run: python3 -B 04-computation/trinomial_two_channel_empty_core_returns_sep06.py
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, gcd, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def multinomial(v):
    return factorial(sum(v)) // prod(factorial(x) for x in v)


def channels(a, b, c, mass):
    result = []
    for y in range(mass + 1):
        numerator = a * mass - (a + c) * y
        if numerator >= 0 and numerator % (a + b) == 0:
            x = numerator // (a + b)
            if x + y <= mass:
                result.append((mass - x - y, x, y))
    return result


def shifted(v, d, j, scale=1):
    return tuple(scale * x + j * y for x, y in zip(v, d))


def check_support(a, b, c):
    g = gcd(a + b, a + c)
    first = channels(a, b, c, g)
    require(len(first) == 2, "wrong input stratum")
    A, B = (a + b) // g, (a + c) // g
    v, w = first
    n, s, t = v[0], w[1], v[2]
    h = B - A
    require(0 <= s < B and 0 <= t < A and h > 0, "residue bounds")
    require(a == A * B + A * s + B * t, "support parameterization")
    d = (h, -B, A)
    require(w == shifted(v, d, 1), "primitive channel step")
    left, right = 2 * t // A, 2 * s // B
    second = channels(a, b, c, 2 * g)
    predicted = [shifted(v, d, j, 2) for j in range(-left, 3 + right)]
    require(second == predicted, "full second slice, including carries")
    alpha, beta = multinomial(v), multinomial(w)
    C, D, E = [multinomial(shifted(v, d, j, 2)) for j in range(3)]
    delta = C * beta**2 - D * alpha * beta + E * alpha**2
    Nv = prod(comb(2 * x, x) for x in v)
    Nw = prod(comb(2 * x, x) for x in w)
    M = prod(comb(x + y, x) for x, y in zip(v, w))
    require(2 * M < Nv and 2 * M <= Nw, "half-budget inequalities")
    require(delta < 0, "internal determinant sign")
    z = -F(alpha, beta)
    normalized_second = sum(F(multinomial(shifted(v, d, j, 2))) * z**j
                            for j in range(-left, 3 + right))
    require(normalized_second < 0, "actual nonvanishing consequence")
    require(2 * g <= a + c, "span bound")
    return left, right, g, z, normalized_second


def rational_coefficients(A, B, z):
    # -B*u+A*v=1; then (1,z^u,z^v)^d=z, all coefficients rational/nonzero.
    for u in range(A + 1):
        if (1 + B * u) % A == 0:
            v = (1 + B * u) // A
            return F(1), z**u, z**v
    raise RuntimeError("Bezout search failed")


def direct_moments(charges, coefficients, limit):
    # Independent repeated Laurent multiplication, no channel parameterization.
    polynomial = {0: F(1)}
    result = []
    for _ in range(limit):
        product = defaultdict(F)
        for exponent, value in polynomial.items():
            for charge, coefficient in zip(charges, coefficients):
                product[exponent + charge] += value * coefficient
        polynomial = {k: v for k, v in product.items() if v}
        result.append(polynomial.get(0, F(0)))
    return result


def determinant(v, w):
    alpha, beta = multinomial(v), multinomial(w)
    return (multinomial(tuple(2 * x for x in v)) * beta**2
            - multinomial(tuple(x + y for x, y in zip(v, w))) * alpha * beta
            + multinomial(tuple(2 * x for x in w)) * alpha**2)


def main():
    print("FINITE-EXACT: two-first-channel trinomial return controls")
    print("universe primitive (-a,b,c), 1<=a,b<=40, b<c<=60")
    census = Counter()
    examples = {}
    for a in range(1, 41):
        for b in range(1, 41):
            for c in range(b + 1, 61):
                if gcd(gcd(a, b), c) != 1:
                    continue
                g = gcd(a + b, a + c)
                if len(channels(a, b, c, g)) != 2:
                    continue
                left, right, *_ = check_support(a, b, c)
                census[left, right] += 1
                examples.setdefault((left, right), (a, b, c))
    print("support_count", sum(census.values()))
    print("carry_profile_counts", sorted(census.items()))
    print("first_examples_by_carry", sorted(examples.items()))

    ratio_checks = 0
    for H in range(1, 81):
        for r in range(H):
            require(comb(2 * r + H, r) <= 2**(H - 1) * comb(2 * r, r),
                    "short-residue binomial lemma")
            ratio_checks += 1
    print("short_residue_binomial_controls", ratio_checks)

    controls = sorted(set(examples.values()) | {(13, 1, 8), (12, 27, 40), (3, 1, 9)})
    for a, b, c in controls:
        left, right, g, z, expected = check_support(a, b, c)
        A, B = (a + b) // g, (a + c) // g
        coefficients = rational_coefficients(A, B, z)
        moments = direct_moments((-a, b, c), coefficients, 2 * g)
        first_nonzero = next(i for i, value in enumerate(moments, 1) if value)
        require(first_nonzero == 2 * g, "independent tuned first-return replay")
        v = channels(a, b, c, g)[0]
        monomial = prod(coefficient**power for coefficient, power in zip(coefficients, v))
        require(moments[-1] / monomial**2 == expected, "independent normalized value")
        positive = direct_moments((-a, b, c), (F(1), F(1), F(1)), g)
        require(next(i for i, value in enumerate(positive, 1) if value) == g,
                "untuned positive first-return control")
        print("convolution_control", (-a, b, c), "carry", (left, right),
              "untuned_first", g, "tuned_first", first_nonzero,
              "normalized_second", str(expected))

    require(determinant((3, 4, 0), (4, 2, 1)) == 0, "zero boundary hostile")
    require(determinant((5, 4, 0), (6, 2, 1)) == 126309456,
            "positive boundary hostile")
    print("three_channel_hostiles (-4,3,10): Delta=0; (-4,5,14): Delta=126309456")

    sharp_count = 0
    for g in range(4, 101):
        if gcd(g, 3) != 1:
            continue
        q = -F(6, g - 2)
        computed = sum(F(multinomial((2 * g - 6 + j, 6 - 2 * j, j))) * q**(3 - j)
                       for j in range(4))
        formula = F(2 * g * (g - 1) * (2 * g - 1) * (23 * g*g - 47 * g + 20),
                    15 * (g - 2)**2)
        require(computed == formula and formula > 0, "width-three sharp family")
        check_support(3, g - 3, 2 * g - 3)
        sharp_count += 1
    print("width_three_sharp_family g=4..100 gcd(g,3)=1 controls", sharp_count)
    print("PASS: all exact checks; general theorem is analytic, not inferred from census")


if __name__ == "__main__":
    main()
