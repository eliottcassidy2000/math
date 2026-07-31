#!/usr/bin/env python3
"""Exact companion for THM-2866.

The theorem proves that

    L(d_n V^3) / L(d_n V)

is strictly increasing in n for every nonzero nonnegative
adjacent-difference cone V.  The proof polarizes the cross product over
four cone atoms and splits the resulting cyclic kernel into three
separately nonnegative pieces.  This companion checks the algebraic
identities, the symbolic symmetric certificate, hostile boundaries, and
a bounded exact census.  Every truth-bearing gate raises explicitly, so
``python`` and ``python -O`` execute the same checks.
"""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations
from math import comb, factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def poly_add(left, right):
    out = [Fraction(0)] * max(len(left), len(right))
    for index in range(len(out)):
        out[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        )
    return out


def poly_mul(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            out[i + j] += first * second
    return out


def poly_scale(scale, poly):
    return [scale * coefficient for coefficient in poly]


def moment(poly):
    return sum(coefficient * factorial(index) for index, coefficient in enumerate(poly))


def f_poly(index):
    out = [Fraction(0)] * (index + 1)
    out[index] = Fraction(1, factorial(index))
    return out


def d_poly(index):
    return poly_add(f_poly(index + 1), poly_scale(-1, f_poly(index)))


def tensor(*indices):
    poly = [Fraction(1)]
    for index in indices:
        poly = poly_mul(poly, d_poly(index))
    return moment(poly)


def h(row, column):
    return comb(row + column, row)


def triple_constants(p, q, r):
    total = p + q + r
    e2 = p * q + p * r + q * r
    e3 = p * q * r
    b_value = 2 * (total + 1) ** 2 + total * e2 - e3
    triple = Fraction(
        factorial(total) * b_value,
        factorial(p + 1) * factorial(q + 1) * factorial(r + 1),
    )
    middle = Fraction(
        factorial(total),
        factorial(p) * factorial(q) * factorial(r),
    )
    high = Fraction(
        factorial(total + 3),
        factorial(p + 1) * factorial(q + 1) * factorial(r + 1),
    )
    return triple, middle, high


def triple_product_rhs(p, q, r):
    total = p + q + r
    triple, middle, high = triple_constants(p, q, r)
    return poly_add(
        poly_add(
            poly_scale(triple, f_poly(total + 1)),
            poly_scale(middle, d_poly(total)),
        ),
        poly_scale(high, d_poly(total + 2)),
    )


def four_tensor_formula(row, p, q, r):
    total = p + q + r
    triple, middle, high = triple_constants(p, q, r)
    return (
        triple * sum(h(row, index) for index in range(total + 1))
        + middle * h(row, total)
        + high * h(row, total + 2)
    )


def kernel_original(row, values):
    total = Fraction(0)
    for singled in range(4):
        column = values[singled]
        triple = values[:singled] + values[singled + 1 :]
        total += (
            four_tensor_formula(row + 1, *triple) * h(row, column)
            - four_tensor_formula(row, *triple) * h(row + 1, column)
        )
    return total


def kernel_components(row, values):
    prefix_part = Fraction(0)
    middle_part = Fraction(0)
    high_part = Fraction(0)
    for singled in range(4):
        column = values[singled]
        triple_indices = values[:singled] + values[singled + 1 :]
        total = sum(triple_indices)
        triple, middle, high = triple_constants(*triple_indices)
        prefix_part += Fraction(
            triple
            * h(row, column)
            * sum(h(row, index) * (index - column) for index in range(total + 1)),
            row + 1,
        )
        middle_part += Fraction(
            middle
            * h(row, column)
            * h(row, total)
            * (total - column),
            row + 1,
        )
        high_part += Fraction(
            high
            * h(row, column)
            * h(row, total + 2)
            * (total + 2 - column),
            row + 1,
        )
    return prefix_part, middle_part, high_part


def monomial_symmetric(variables, partition):
    padded = tuple(partition) + (0,) * (len(variables) - len(partition))
    exponent_vectors = set(permutations(padded))
    return sum(
        sp.prod(variable**exponent for variable, exponent in zip(variables, vector))
        for vector in exponent_vectors
    )


def symbolic_prefix_certificate():
    variables = sp.symbols("x0:4", nonnegative=True, integer=True)
    total = sum(variables)

    def b_value(triple):
        subtotal = sum(triple)
        e2 = sum(triple[i] * triple[j] for i, j in combinations(range(3), 2))
        e3 = sp.prod(triple)
        return 2 * (subtotal + 1) ** 2 + subtotal * e2 - e3

    left = sum(
        (variables[i] + 1)
        * b_value([variables[j] for j in range(4) if j != i])
        * (total - 3 * variables[i])
        for i in range(4)
    )

    m = lambda partition: monomial_symmetric(variables, partition)
    right = (
        2 * m((1,))
        + 8 * m((2,))
        + 4 * m((1, 1))
        + 6 * m((3,))
        + 4 * m((2, 1))
        + 12 * m((1, 1, 1))
        + 4 * (m((3, 1)) - m((2, 2)))
        + 4 * m((2, 1, 1))
        + 32 * m((1, 1, 1, 1))
        + 2 * (m((3, 1, 1)) - m((2, 2, 1)))
        + 8 * m((2, 1, 1, 1))
    )
    require(sp.expand(left - right) == 0, "prefix symmetric expansion")

    pair_square = sum(
        variables[i]
        * variables[j]
        * (variables[i] - variables[j]) ** 2
        for i, j in combinations(range(4), 2)
    )
    require(
        sp.expand(m((3, 1)) - 2 * m((2, 2)) - pair_square) == 0,
        "degree-four pair-square certificate",
    )

    triple_square = 0
    for omitted in range(4):
        triple = [variables[j] for j in range(4) if j != omitted]
        triple_square += sp.prod(triple) * (
            sum(value**2 for value in triple)
            - sum(triple[i] * triple[j] for i, j in combinations(range(3), 2))
        )
    require(
        sp.expand(m((3, 1, 1)) - m((2, 2, 1)) - triple_square) == 0,
        "degree-five triple-square certificate",
    )

    return len(sp.Poly(sp.expand(left), *variables).terms())


def endpoint_coefficient_certificate():
    g0, g1, g2 = sp.symbols("g0 g1 g2", nonzero=True)
    a0, a1, a2, a3, a4 = sp.symbols("a0 a1 a2 a3 a4")
    x, y = sp.symbols("x y")
    r0 = a0 / g0
    r2 = a4 / g2
    r1_left = (2 * a1 * g0 - a0 * g1) / g0**2
    r1_right = (2 * a3 * g2 - a4 * g1) / g2**2
    quadratic = g0 * x**2 + 2 * g1 * x * y + g2 * y**2
    quotient_left = r0 * x**2 + 2 * r1_left * x * y + r2 * y**2
    product_left = sp.Poly(sp.expand(quadratic * quotient_left), x, y)
    require(product_left.coeff_monomial(x**4) == a0, "left x4 coefficient")
    require(product_left.coeff_monomial(x**3 * y) == 4 * a1, "left x3y coefficient")
    require(product_left.coeff_monomial(y**4) == a4, "left y4 coefficient")
    require(
        sp.factor(
            product_left.coeff_monomial(x * y**3) - 4 * a3
            - 2 * g2 * (r1_left - r1_right)
        )
        == 0,
        "endpoint mismatch coefficient",
    )


def interval_face(a, b, c):
    u = poly_add(f_poly(b), poly_scale(-1, f_poly(a)))
    v = poly_add(f_poly(c), poly_scale(-1, f_poly(b)))
    u2 = poly_mul(u, u)
    v2 = poly_mul(v, v)
    g0 = moment(u2)
    g1 = moment(poly_mul(u, v))
    g2 = moment(v2)
    a0 = moment(poly_mul(u2, u2))
    a1 = moment(poly_mul(poly_mul(u2, u), v))
    a3 = moment(poly_mul(u, poly_mul(v2, v)))
    a4 = moment(poly_mul(v2, v2))
    left = 2 * a1 * g0 - a0 * g1
    right = 2 * a3 * g2 - a4 * g1
    mismatch = left * g2**2 - right * g0**2
    return left, right, mismatch


def main():
    endpoint_coefficient_certificate()
    symbolic_terms = symbolic_prefix_certificate()

    semiring_cells = 0
    for first in range(11):
        for second in range(11):
            require(
                poly_mul(f_poly(first), f_poly(second))
                == poly_scale(comb(first + second, first), f_poly(first + second)),
                "f times f multiplication",
            )
            fd_prefix = comb(first + second, first - 1) if first > 0 else 0
            require(
                poly_mul(f_poly(first), d_poly(second))
                == poly_add(
                    poly_scale(fd_prefix, f_poly(first + second)),
                    poly_scale(comb(first + second + 1, first), d_poly(first + second)),
                ),
                "f times d multiplication",
            )
            require(
                poly_mul(d_poly(first), d_poly(second))
                == poly_add(
                    poly_scale(comb(first + second, first), f_poly(first + second)),
                    poly_scale(
                        Fraction(
                            factorial(first + second + 2),
                            factorial(first + 1) * factorial(second + 1),
                        ),
                        d_poly(first + second + 1),
                    ),
                ),
                "d times d multiplication",
            )
            semiring_cells += 3

    product_cells = 0
    for p in range(7):
        for q in range(7):
            for r in range(7):
                direct = poly_mul(poly_mul(d_poly(p), d_poly(q)), d_poly(r))
                require(direct == triple_product_rhs(p, q, r), "triple product identity")
                triple, middle, high = triple_constants(p, q, r)
                require(triple == tensor(p, q, r), "triple constant")
                require(triple > 0 and middle > 0 and high > 0, "positive product constants")
                product_cells += 1

    tensor_cells = 0
    for row in range(7):
        for p in range(6):
            for q in range(6):
                for r in range(6):
                    require(
                        four_tensor_formula(row, p, q, r) == tensor(row, p, q, r),
                        "four-tensor product formula",
                    )
                    tensor_cells += 1

    kernel_cells = 0
    prefix_zero = 0
    middle_zero = 0
    for row in range(7):
        for values in combinations_with_replacement(range(13), 4):
            components = kernel_components(row, values)
            original = kernel_original(row, values)
            require(sum(components) == original, "cyclic kernel decomposition")
            require(components[0] >= 0, "prefix kernel sign")
            require(components[1] >= 0, "middle kernel sign")
            require(components[2] > 0, "high kernel strict sign")
            require(original > 0, "full cyclic kernel strict sign")
            require((components[0] == 0) == (sum(values) == 0), "prefix equality")
            require((components[1] == 0) == (sum(values) == 0), "middle equality")
            prefix_zero += int(components[0] == 0)
            middle_zero += int(components[1] == 0)
            kernel_cells += 1

    # The integration-by-parts step in the face application is deliberately
    # restricted to a>=1.  At a=0, U(0)=-1 contributes the missing boundary
    # terms in both the quadratic and quartic identities.
    u0 = poly_add(f_poly(1), poly_scale(-1, f_poly(0)))
    u0_prime = [coefficient * index for index, coefficient in enumerate(u0)][1:]
    require(u0[0] == -1, "boundary hostile value")
    require(
        moment(poly_mul(u0, u0))
        - 2 * moment(poly_mul(u0, u0_prime))
        == 1,
        "quadratic boundary correction",
    )
    require(
        moment(poly_mul(poly_mul(u0, u0), poly_mul(u0, u0)))
        - 4 * moment(poly_mul(poly_mul(poly_mul(u0, u0), u0), u0_prime))
        == 1,
        "quartic boundary correction",
    )

    face_cells = 0
    right_equalities = []
    for a in range(1, 24):
        for b in range(a + 1, 25):
            for c in range(b + 1, 26):
                left, right, mismatch = interval_face(a, b, c)
                require(left > 0, "top-face left endpoint sign")
                require(right <= 0, "top-face right endpoint sign")
                require(mismatch > 0, "top-face endpoint mismatch")
                if right == 0:
                    right_equalities.append((a, b, c))
                    require(
                        b == a + 1 and c == b + 1,
                        "right endpoint equality classification",
                    )
                face_cells += 1

    expected_equalities = [(a, a + 1, a + 2) for a in range(1, 24)]
    require(right_equalities == expected_equalities, "bounded equality ledger")

    # Hostile boundary for the cone hypothesis: allowing a negative
    # adjacent-difference coefficient reverses the response ladder.
    signed_w = poly_add(poly_scale(3, d_poly(0)), poly_scale(-1, d_poly(1)))
    signed_w3 = poly_mul(poly_mul(signed_w, signed_w), signed_w)
    signed_control = (
        moment(poly_mul(d_poly(0), signed_w)),
        moment(poly_mul(d_poly(1), signed_w)),
        moment(poly_mul(d_poly(0), signed_w3)),
        moment(poly_mul(d_poly(1), signed_w3)),
    )
    require(signed_control == (2, 1, 18, -27), "signed-cone hostile values")
    require(
        Fraction(signed_control[2], signed_control[0])
        > Fraction(signed_control[3], signed_control[1]),
        "signed-cone hostile reverses the ratio ladder",
    )

    print("THM2866 exact certificate")
    print(f"symbolic_prefix_terms={symbolic_terms}")
    print(f"positive_semiring_cells={semiring_cells}")
    print(f"triple_product_cells={product_cells}")
    print(f"four_tensor_cells={tensor_cells}")
    print(
        f"cyclic_kernel_cells={kernel_cells}; "
        f"prefix_zero={prefix_zero}; middle_zero={middle_zero}"
    )
    print(
        f"top_face_cells={face_cells}; "
        f"right_endpoint_equalities={len(right_equalities)}"
    )
    print("boundary_hostile=a0_requires_U(0)_corrections")
    print("signed_hostile=(D0,D1,C0,C1)=(2,1,18,-27)")
    print("THM2866_ALL_EXACT_CHECKS_PASSED")


if __name__ == "__main__":
    main()
