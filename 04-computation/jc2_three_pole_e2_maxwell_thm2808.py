#!/usr/bin/env python3
"""Exact finite controls for THM-2808.

The proof in THM-2808 is all-degree.  This companion checks the normalized
three-pole eliminant, its converse gates, and the noncrossing-chord count in
the explicitly declared finite universe below.
"""

import ast
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    return any(
        isinstance(node, ast.Assert)
        for node in ast.walk(ast.parse(Path(path).read_text(encoding="utf-8")))
    )


def reduce_mod(expr, variable, modulus):
    numerator, denominator = sp.cancel(expr).as_numer_denom()
    modulus_poly = sp.Poly(modulus, variable, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, variable, domain=sp.QQ)
    require(
        sp.gcd(denominator_poly, modulus_poly).degree() == 0,
        "nonunit denominator in quotient",
    )
    inverse = sp.invert(denominator_poly, modulus_poly)
    return sp.rem(
        sp.Poly(numerator, variable, domain=sp.QQ) * inverse,
        modulus_poly,
    ).as_expr()


def zero_mod_x(expr, x, parameter, modulus):
    return all(
        reduce_mod(coefficient, parameter, modulus) == 0
        for coefficient in sp.Poly(sp.expand(expr), x).all_coeffs()
    )


def coprime_to(expr, parameter, modulus):
    numerator = sp.cancel(expr).as_numer_denom()[0]
    return (
        sp.gcd(
            sp.Poly(numerator, parameter, domain=sp.QQ),
            sp.Poly(modulus, parameter, domain=sp.QQ),
        ).degree()
        == 0
    )


def normalized_data(a, b, c, x, parameter):
    degree = a + b + c
    critical_sum_numerator = (a + c) + (a + b) * parameter
    critical_product = sp.cancel(a * parameter / degree)
    critical_sum = sp.cancel(critical_sum_numerator / degree)
    critical_quadratic = sp.expand(
        x**2 - critical_sum * x + critical_product
    )
    collision_discriminant = sp.expand(
        critical_sum_numerator**2 - 4 * degree * a * parameter
    )

    denominator = sp.expand(
        x**a * (x - 1) ** b * (x - parameter) ** c
    )
    remainder = sp.rem(denominator, critical_quadratic, x)
    remainder_poly = sp.Poly(remainder, x)
    secant_coefficient = sp.expand(remainder_poly.coeff_monomial(x))
    critical_value = sp.expand(remainder_poly.coeff_monomial(1))
    maxwell = sp.Poly(
        sp.cancel(secant_coefficient / collision_discriminant),
        parameter,
        domain=sp.QQ,
    ).as_expr()

    numerator = sp.expand(denominator - critical_value)
    simple_factor, square_remainder = sp.div(
        numerator, critical_quadratic**2, x
    )
    pole_support = x * (x - 1) * (x - parameter)
    pole_factor = (
        x ** (a - 1)
        * (x - 1) ** (b - 1)
        * (x - parameter) ** (c - 1)
    )

    return {
        "degree": degree,
        "critical_sum_numerator": critical_sum_numerator,
        "critical_quadratic": critical_quadratic,
        "collision_discriminant": collision_discriminant,
        "denominator": denominator,
        "remainder": remainder,
        "secant_coefficient": secant_coefficient,
        "critical_value": critical_value,
        "maxwell": maxwell,
        "numerator": numerator,
        "simple_factor": simple_factor,
        "square_remainder": square_remainder,
        "pole_support": pole_support,
        "pole_factor": pole_factor,
    }


def check_partition(a, b, c, x, parameter):
    data = normalized_data(a, b, c, x, parameter)
    degree = data["degree"]
    maxwell = data["maxwell"]
    collision = data["collision_discriminant"]

    require(
        sp.expand(
            data["secant_coefficient"] - collision * maxwell
        )
        == 0,
        "collision factorization",
    )
    require(sp.degree(data["secant_coefficient"], parameter) == degree - 1,
            "raw secant degree")
    require(sp.degree(maxwell, parameter) == degree - 3, "Maxwell degree")
    require(sp.discriminant(maxwell, parameter) != 0, "Maxwell squarefree")
    require(
        all(
            coprime_to(gate, parameter, maxwell)
            for gate in (
                parameter,
                parameter - 1,
                collision,
                data["critical_value"],
            )
        ),
        "parameter/collision/value gate",
    )

    require(
        zero_mod_x(
            data["square_remainder"], x, parameter, maxwell
        ),
        "critical quadratic square division",
    )
    require(
        reduce_mod(
            sp.LC(sp.Poly(data["simple_factor"], x)),
            parameter,
            maxwell,
        )
        == 1,
        "simple factor monic",
    )
    require(sp.degree(data["simple_factor"], x) == degree - 4,
            "simple factor degree")

    derivative_identity = sp.expand(
        sp.diff(data["numerator"], x) * data["denominator"]
        - data["numerator"] * sp.diff(data["denominator"], x)
        - degree
        * data["critical_value"]
        * data["pole_factor"]
        * data["critical_quadratic"]
    )
    require(derivative_identity == 0, "cleared response derivative")

    if degree > 4:
        simple_discriminant = sp.discriminant(data["simple_factor"], x)
    else:
        simple_discriminant = sp.Integer(1)
    simple_resultant = sp.resultant(
        data["simple_factor"],
        data["critical_quadratic"] * data["pole_support"],
        x,
    )
    critical_resultant = sp.resultant(
        data["critical_quadratic"], data["pole_support"], x
    )
    require(
        all(
            coprime_to(gate, parameter, maxwell)
            for gate in (
                simple_discriminant,
                simple_resultant,
                critical_resultant,
            )
        ),
        "squarefree/disjoint converse gates",
    )

    # The marked noncrossing-chord atlas has one chart for each positive split
    # of one of the three pole cycles.
    marked_charts = sum(part - 1 for part in (a, b, c))
    require(marked_charts == degree - 3, "marked Nielsen count")
    chord_charts = []
    for split_index in range(3):
        chord_charts.extend(chord_trace((a, b, c), split_index))
    require(len(chord_charts) == marked_charts, "chord atlas size")

    return data


def proportional(left, right, variable):
    left_poly = sp.Poly(sp.expand(left), variable, domain=sp.QQ)
    right_poly = sp.Poly(sp.expand(right), variable, domain=sp.QQ)
    return sp.expand(
        left_poly.LC() * right_poly.as_expr()
        - right_poly.LC() * left_poly.as_expr()
    ) == 0


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def cycle_lengths(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        point = start
        length = 0
        while point in unseen:
            unseen.remove(point)
            point = permutation[point]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def chord_trace(parts, composite_index):
    degree = sum(parts)
    composite = parts[composite_index]
    left = parts[(composite_index + 1) % 3]
    right = parts[(composite_index + 2) % 3]
    full_cycle = tuple((index + 1) % degree for index in range(degree))
    traces = []
    for cut in range(1, composite):
        gaps = (left, cut, right, composite - cut)
        endpoints = (
            0,
            gaps[0],
            gaps[0] + gaps[1],
            gaps[0] + gaps[1] + gaps[2],
        )
        involution = list(range(degree))
        for first, second in (
            (endpoints[0], endpoints[1]),
            (endpoints[2], endpoints[3]),
        ):
            involution[first] = second
            involution[second] = first
        traced = cycle_lengths(compose(tuple(involution), full_cycle))
        require(traced == tuple(sorted(parts)), "noncrossing chord trace")
        traces.append((composite_index, cut, gaps))
    return traces


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    x, parameter = sp.symbols("x lambda")
    print("THM-2808 THREE-POLE E=2 MAXWELL AUDIT")
    print("universe: 4 <= N <= 10, all 1 <= a <= b <= c, a+b+c=N")
    print("N | partition | deg raw | deg Q | marked charts")

    cases = 0
    roots = 0
    for degree in range(4, 11):
        for a in range(1, degree - 1):
            for b in range(a, degree - a):
                c = degree - a - b
                if c < b:
                    continue
                data = check_partition(a, b, c, x, parameter)
                degree_q = int(sp.degree(data["maxwell"], parameter))
                print(
                    f"{degree:2d} | ({a},{b},{c})"
                    f" | {degree - 1:7d} | {degree_q:5d}"
                    f" | {degree_q:13d}"
                )
                cases += 1
                roots += degree_q

                # Swapping the poles at 0 and 1 sends lambda to 1-lambda.
                swapped = normalized_data(b, a, c, x, parameter)
                transformed = swapped["maxwell"].subs(
                    parameter, 1 - parameter
                )
                require(
                    proportional(
                        data["maxwell"], transformed, parameter
                    ),
                    "0/1 pole-swap covariance",
                )

    print(f"exact ordered-partition representatives checked: {cases}")
    print(f"simple Maxwell roots represented: {roots}")
    print("raw eliminant = collision discriminant * Maxwell polynomial")
    print("all bounded converse gates and marked chord counts: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
