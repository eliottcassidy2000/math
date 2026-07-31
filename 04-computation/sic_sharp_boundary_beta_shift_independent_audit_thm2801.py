#!/usr/bin/env python3
"""Independent sparse/finite-difference audit for THM-2801.

This companion deliberately avoids SymPy.  Polynomials are dictionaries
over ``Fraction`` and the angular proof is checked through its alternating
finite-difference form rather than the endpoint-flat-jet form used in the
main proof.
"""

import ast
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def has_asserts(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Assert) for node in ast.walk(tree))


def clean(poly):
    return {powers: value for powers, value in poly.items() if value}


def add(left, right):
    result = dict(left)
    for powers, value in right.items():
        result[powers] = result.get(powers, Fraction(0)) + value
    return clean(result)


def scale(value, poly):
    return clean({powers: value * coefficient for powers, coefficient in poly.items()})


def multiply(left, right):
    result = {}
    for powers_left, value_left in left.items():
        for powers_right, value_right in right.items():
            powers = tuple(a + b for a, b in zip(powers_left, powers_right))
            result[powers] = (
                result.get(powers, Fraction(0)) + value_left * value_right
            )
    return clean(result)


def power(poly, exponent):
    variables = len(next(iter(poly)))
    result = {(0,) * variables: Fraction(1)}
    base = poly
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        exponent //= 2
        if exponent:
            base = multiply(base, base)
    return result


def monomial(powers, coefficient=1):
    return {tuple(powers): Fraction(coefficient)}


def falling(top, length):
    if length > top:
        return 0
    return factorial(top) // factorial(top - length)


def image_operator(poly, pairs):
    result = {}
    for powers, coefficient in poly.items():
        alpha = powers[:pairs]
        beta = powers[pairs:]
        if any(a > b for a, b in zip(alpha, beta)):
            continue
        output = tuple(b - a for a, b in zip(alpha, beta))
        value = coefficient
        for a, b in zip(alpha, beta):
            value *= falling(b, a)
        result[output] = result.get(output, Fraction(0)) + value
    return clean(result)


def odd_double_factorial(index):
    value = 1
    for term in range(1, index + 1, 2):
        value *= term
    return value


def two_pair_polynomial():
    # Variable order: xi1,xi2,z1,z2.
    R = add(monomial((1, 0, 1, 0)), monomial((0, 1, 0, 1)))
    Z = monomial((1, 0, 0, 1))
    W = monomial((0, 1, 1, 0), 2)
    T = add(monomial((1, 0, 1, 0)), monomial((0, 1, 0, 1), -1))
    first = add(R, Z)
    second = add(
        multiply(power(R, 2), W),
        scale(
            Fraction(-1, 2),
            multiply(add(scale(2, R), Z), power(T, 2)),
        ),
    )
    return multiply(first, second), Z


def audit_two_pair():
    F, Z = two_pair_polynomial()
    z2 = monomial((0, 0, 0, 1))
    current = monomial((0, 0, 0, 0))
    values = []
    for m in range(1, 6):
        current = multiply(current, F)
        expected = Fraction(
            factorial(4 * m + 2) * factorial(m),
            odd_double_factorial(2 * m + 1),
        )
        require(image_operator(current, 2) == {}, f"two-pair pure m={m}")
        require(
            image_operator(multiply(Z, current), 2) == {(0, 0): expected},
            f"two-pair scalar observer m={m}",
        )
        require(
            image_operator(multiply(z2, current), 2) == {(1, 0): expected},
            f"two-pair coordinate observer m={m}",
        )
        values.append(expected)
    return len(F), values


def audit_alternating_beta():
    checks = 0
    ladder_checks = 0
    for m in range(1, 81):
        def coefficient(index):
            if index < 0:
                return Fraction(0)
            return sum(
                Fraction(
                    (-1) ** j * comb(m, j) * comb(m + 2 * j, index),
                    2 * j + 1,
                )
                for j in range(m + 1)
                if index <= m + 2 * j
            )

        require(coefficient(m) == 0, f"alternating pure coefficient m={m}")
        beta = Fraction(
            2**m * factorial(m), odd_double_factorial(2 * m + 1)
        )
        require(coefficient(m - 1) == beta, f"alternating beta m={m}")
        for s in range(1, min(9, m + 4)):
            expected = beta * comb(m - 1, s - 1) if m >= s else 0
            require(coefficient(m - s) == expected,
                    f"alternating observer ladder m={m},s={s}")
            ladder_checks += 1
        checks += 1
    return checks, ladder_checks


def three_pair_polynomial(r):
    # Variable order: tau,w,v,t,z,y.
    tau = monomial((1, 0, 0, 0, 0, 0))
    w = monomial((0, 1, 0, 0, 0, 0))
    v = monomial((0, 0, 1, 0, 0, 0))
    t = monomial((0, 0, 0, 1, 0, 0))
    z = monomial((0, 0, 0, 0, 1, 0))
    y = monomial((0, 0, 0, 0, 0, 1))
    return multiply(
        power(tau, r),
        multiply(
            add(t, z),
            add(
                multiply(w, power(t, r)),
                scale(-1, multiply(v, multiply(y, power(add(t, y), r - 1)))),
            ),
        ),
    )


def expected_three_pair_output(r, m):
    result = {(1, 0, 0): Fraction(factorial(r * m + 1) * factorial(m))}
    if m == 1:
        result[(0, 1, 0)] = Fraction(factorial(r))
    if m == 2 and r > 1:
        result[(0, 0, 1)] = Fraction(-4 * (r - 1) * factorial(2 * r))
    return result


def audit_three_pair():
    z = monomial((0, 0, 0, 0, 1, 0))
    checks = 0
    support = {}
    for r in range(1, 5):
        f = three_pair_polynomial(r)
        support[r] = len(f)
        current = monomial((0, 0, 0, 0, 0, 0))
        for m in range(1, 5):
            current = multiply(current, f)
            require(image_operator(current, 3) == {},
                    f"three-pair pure r={r},m={m}")
            require(
                image_operator(multiply(z, current), 3)
                == expected_three_pair_output(r, m),
                f"three-pair marked r={r},m={m}",
            )
            checks += 1
    require(support[1] == 4, "r=1 is not four-term")
    require(support[3] == 8, "r=3 user support changed")
    return checks, support


def audit_bilinear_boundary():
    # Strictly upper triangular controls in n=2 and n=3.
    f2 = monomial((1, 0, 0, 1))
    f3 = add(
        monomial((1, 0, 0, 0, 1, 0)),
        monomial((0, 1, 0, 0, 0, 1)),
    )
    for f, pairs in ((f2, 2), (f3, 3)):
        for m in range(1, 9):
            require(image_operator(power(f, m), pairs) == {},
                    f"nilpotent bilinear moment pairs={pairs},m={m}")

    # A nonnilpotent diagonal control makes the determinant series active.
    diagonal = monomial((1, 0, 1, 0))
    require(image_operator(power(diagonal, 5), 2) == {(0, 0): Fraction(120)},
            "nonnilpotent bilinear positive control")
    return 17


def main():
    script = Path(__file__)
    require(not has_asserts(script), "truth-bearing Python assert found")
    support, values = audit_two_pair()
    beta_checks, ladder_checks = audit_alternating_beta()
    family_checks, family_support = audit_three_pair()
    bilinear_checks = audit_bilinear_boundary()

    print("THM-2801 INDEPENDENT SPARSE AND FINITE-DIFFERENCE AUDIT")
    print(f"two-pair sparse support: {support}")
    print(f"two-pair literal powers checked: {len(values)}")
    print(f"two-pair fifth marked value: {values[-1]}")
    print(f"alternating Beta identities checked: {beta_checks}")
    print(f"alternating observer-ladder identities checked: {ladder_checks}")
    print(f"general three-pair literal pairs checked: {family_checks}")
    print(f"three-pair supports r=1..4: {tuple(family_support[r] for r in range(1, 5))}")
    print(f"bilinear boundary controls checked: {bilinear_checks}")
    print("exceptional z/y channels and coordinate-only observer: PASS")


if __name__ == "__main__":
    main()
