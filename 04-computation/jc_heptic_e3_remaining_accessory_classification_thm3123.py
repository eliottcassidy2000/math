#!/usr/bin/env python3
"""Exact companion for THM-3123 (remaining heptic accessory strata).

This is the author-side exact companion.  It derives both accessory ideals,
checks the actual saturated geometric gates in their ordered quotient rings,
reconstructs the response identities, and supplies a separate permutation
census.  It uses exact SymPy arithmetic and no Python ``assert`` statements.
"""

from __future__ import annotations

from itertools import permutations

import sympy as sp


x, lam, mu, u, v = sp.symbols("x lambda mu u v")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def cycles(permutation: tuple[int, ...]) -> tuple[frozenset[int], ...]:
    seen: set[int] = set()
    result: list[frozenset[int]] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle: set[int] = set()
        point = start
        while point not in seen:
            seen.add(point)
            cycle.add(point)
            point = permutation[point]
        result.append(frozenset(cycle))
    return tuple(result)


def matchings(points: tuple[int, ...], count: int):
    if count == 0:
        yield ()
        return
    if len(points) < 2 * count:
        return
    first = points[0]
    yield from matchings(points[1:], count)
    for index in range(1, len(points)):
        second = points[index]
        remainder = points[1:index] + points[index + 1 :]
        for tail in matchings(remainder, count - 1):
            yield ((first, second),) + tail


def rotation_key(permutation: tuple[int, ...]) -> tuple[int, ...]:
    return min(
        tuple(
            (permutation[(point - shift) % 7] + shift) % 7
            for point in range(7)
        )
        for shift in range(7)
    )


def marked_key(
    permutation: tuple[int, ...], labelled_cycles: tuple[frozenset[int], ...]
):
    return min(
        (
            tuple(
                (permutation[(point - shift) % 7] + shift) % 7
                for point in range(7)
            ),
            tuple(
                tuple(sorted((point + shift) % 7 for point in cycle))
                for cycle in labelled_cycles
            ),
        )
        for shift in range(7)
    )


def group_closure(generators: list[tuple[int, ...]]) -> set[tuple[int, ...]]:
    identity = tuple(range(7))
    group = {identity}
    frontier = [identity]
    while frontier:
        left = frontier.pop()
        for right in generators:
            product = compose(left, right)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def permutation_power(
    permutation: tuple[int, ...], exponent: int
) -> tuple[int, ...]:
    result = tuple(range(len(permutation)))
    for _ in range(exponent):
        result = compose(permutation, result)
    return result


def explicit_s7_certificates() -> None:
    total_cycle = tuple((point + 1) % 7 for point in range(7))
    tau_4111 = (1, 0, 3, 2, 4, 6, 5)
    tau_3211 = (
        (0, 2, 1, 6, 5, 4, 3),
        (1, 0, 2, 6, 5, 4, 3),
        (5, 2, 1, 4, 3, 0, 6),
    )
    require(
        permutation_power(
            compose(tau_4111, permutation_power(total_cycle, 3)), 5
        )
        == (0, 4, 2, 3, 1, 5, 6),
        "the (4,1,1,1) transposition word differs",
    )
    expected = (
        (0, 1, 2, 5, 4, 3, 6),
        (0, 1, 2, 5, 4, 3, 6),
        (0, 1, 2, 3, 4, 6, 5),
    )
    require(
        tuple(
            permutation_power(compose(involution, total_cycle), 3)
            for involution in tau_3211
        )
        == expected,
        "a (3,2,1,1) transposition word differs",
    )


def quotient_unit(
    expression: sp.Expr,
    quotient: sp.GroebnerBasis,
    modulus: sp.Expr,
) -> bool:
    numerator = sp.cancel(expression).as_numer_denom()[0]
    remainder = sp.factor(quotient.reduce(numerator)[1])
    if remainder.has(mu):
        return False
    return sp.gcd(
        sp.Poly(remainder, lam, domain=sp.QQ),
        sp.Poly(modulus, lam, domain=sp.QQ),
    ).degree() == 0


def exact_accessory(passport: tuple[int, int, int, int]) -> dict[str, object]:
    a, b, c, d = passport
    require(c == d == 1 and sum(passport) == 7, "probe expects the two heptic lanes")
    divisor = sp.expand(x**a * (x - 1) ** b * (x - lam) * (x - mu))
    pole_product = sp.expand(x * (x - 1) * (x - lam) * (x - mu))
    critical = sp.expand(
        a * (x - 1) * (x - lam) * (x - mu)
        + b * x * (x - lam) * (x - mu)
        + x * (x - 1) * (x - mu)
        + x * (x - 1) * (x - lam)
    )
    energy = critical / 7
    require(
        sp.cancel(sp.diff(divisor, x) - 7 * divisor * energy / pole_product) == 0,
        "logarithmic derivative failed",
    )
    simple, remainder = sp.div(divisor, energy**2, x)
    quotient, derivative_remainder = sp.div(sp.diff(remainder, x), energy, x)
    require(derivative_remainder == 0, "E does not divide R'")
    coefficients = []
    for monomial in (x, 1):
        numerator = sp.together(sp.Poly(quotient, x).coeff_monomial(monomial)).as_numer_denom()[0]
        coefficients.append(sp.primitive(sp.Poly(numerator, lam, mu))[1].as_expr())
    symmetric = []
    for coefficient in coefficients:
        expression, residual, substitutions = sp.symmetrize(
            coefficient, [lam, mu], formal=True
        )
        require(residual == 0, "accessory equation is not symmetric")
        symmetric.append(
            sp.expand(
                expression.subs(
                    {substitutions[0][0]: u, substitutions[1][0]: v}
                )
            )
        )
    if passport == (4, 1, 1, 1):
        expected = [8 * u**2 + 9 * u - 7 * v + 8, -(25 * u**2 + 25 * u * v + 25 * u + 11 * v)]
        v_of_u = (8 * u**2 + 9 * u + 8) / 7
        cubic = 100 * u**3 + 244 * u**2 + 237 * u + 44
        expected_discriminant = -480200000
        expected_simple = x + 5 * (u + 1) / 7
        expected_accessory = 80 * v**2 * (u + 1) / 343
        expected_modulus = (
            2500 * lam**6
            + 6100 * lam**5
            + 10121 * lam**4
            + 12558 * lam**3
            + 10121 * lam**2
            + 6100 * lam
            + 2500
        )
        gate_polynomials = {
            "v": 8 * u**2 + 9 * u + 8,
            "one_minus_u_plus_v": 8 * u**2 + 2 * u + 15,
            "lambda_minus_mu_squared": 25 * u**2 + 36 * u + 32,
            "critical_discriminant": 7767 * u**2 + 12916 * u + 6192,
            "accessory_extra": u + 1,
            "simple_double_resultant": 8 * u**2 + 9 * u + 8,
        }
        expected_resultants = {
            "v": 3430000,
            "one_minus_u_plus_v": 68600000,
            "lambda_minus_mu_squared": 68600000,
            "critical_discriminant": 837402343750000,
            "accessory_extra": 49,
            "simple_double_resultant": 3430000,
        }
    else:
        expected = [24 * u**2 - 16 * u - 21 * v - 16, -(40 * u**2 + 50 * u * v - 32 * u - 61 * v)]
        v_of_u = (24 * u**2 - 16 * u - 16) / 21
        cubic = 75 * u**3 - 89 * u**2 - 31 * u + 61
        expected_discriminant = -149361408
        expected_simple = x + (5 * u - 4) / 7
        expected_accessory = 9 * v**2 * (5 * u - 4) / 343
        expected_modulus = (
            1875 * lam**6
            - 2225 * lam**5
            - 1967 * lam**4
            - 443 * lam**3
            + 2072 * lam**2
            + 1472 * lam
            + 512
        )
        gate_polynomials = {
            "v": 3 * u**2 - 2 * u - 2,
            "one_minus_u_plus_v": 24 * u**2 - 37 * u + 5,
            "lambda_minus_mu_squared": 75 * u**2 - 64 * u - 64,
            "critical_discriminant": 8817 * u**2 - 5332 * u - 4408,
            "accessory_extra": 5 * u - 4,
            "simple_double_resultant": 10 * u**2 + 5 * u - 23,
        }
        expected_resultants = {
            "v": 27783,
            "one_minus_u_plus_v": 36006768,
            "lambda_minus_mu_squared": 781396875,
            "critical_discriminant": 593373251953125,
            "accessory_extra": -2205,
            "simple_double_resultant": -20003760,
        }
    require(symmetric == expected, "symmetric accessory equations differ")
    require(
        sp.factor(expected[1].subs(v, v_of_u) / cubic)
        in (-sp.Rational(2, 7), -sp.Rational(16, 21)),
        "cubic elimination differs",
    )
    discriminant = sp.discriminant(cubic, u)
    require(discriminant == expected_discriminant, "cubic discriminant differs")
    resultants = {
        name: sp.resultant(cubic, polynomial, u)
        for name, polynomial in gate_polynomials.items()
    }
    require(resultants == expected_resultants, "symmetric boundary resultants differ")

    # Direct quotient-ring reconstruction in the original ordered coordinates.
    groebner = sp.groebner(coefficients, mu, lam, order="lex", domain=sp.QQ)
    require(len(groebner.polys) == 2, "ordered accessory basis is not triangular")
    require(
        sp.Poly(groebner.polys[-1].as_expr(), lam, domain=sp.QQ).monic()
        == sp.Poly(expected_modulus, lam, domain=sp.QQ).monic(),
        "ordered accessory sextic differs",
    )
    accessory = -remainder.subs(x, 0)
    response_constant = -7 * accessory

    def zero_mod(expression: sp.Expr) -> bool:
        numerator = sp.cancel(expression).as_numer_denom()[0]
        return all(
            groebner.reduce(coefficient)[1] == 0
            for coefficient in sp.Poly(sp.expand(numerator), x).all_coeffs()
        )

    numerator = divisor + accessory
    expected_simple_ordered = expected_simple.subs(u, lam + mu)
    expected_accessory_ordered = expected_accessory.subs(
        {u: lam + mu, v: lam * mu}
    )
    require(
        sp.expand(simple - expected_simple_ordered) == 0,
        "simple factor differs",
    )
    require(
        zero_mod(accessory - expected_accessory_ordered),
        "accessory constant differs",
    )
    require(zero_mod(remainder + accessory), "R is not constant")
    require(zero_mod(numerator - simple * energy**2), "D+A != S E^2")
    require(
        zero_mod(
            (sp.diff(numerator, x) * divisor - numerator * sp.diff(divisor, x))
            * pole_product
            - response_constant * energy * divisor
        ),
        "response derivative failed",
    )

    jacobian = sp.det(
        sp.Matrix(
            [
                [sp.diff(coefficients[0], lam), sp.diff(coefficients[0], mu)],
                [sp.diff(coefficients[1], lam), sp.diff(coefficients[1], mu)],
            ]
        )
    )
    actual_gates = {
        "pole_collision": lam * mu * (lam - 1) * (mu - 1) * (lam - mu),
        "critical_discriminant": sp.discriminant(energy, x),
        "nonzero_accessory": accessory,
        "simple_double_separation": sp.resultant(simple, energy, x),
        "jacobian": jacobian,
    }
    require(
        all(
            quotient_unit(gate, groebner, expected_modulus)
            for gate in actual_gates.values()
        ),
        "an actual geometric saturation gate is not a unit",
    )
    return {
        "passport": passport,
        "symmetric_equations": expected,
        "v_of_u": v_of_u,
        "cubic": cubic,
        "cubic_discriminant": discriminant,
        "boundary_resultants": resultants,
        "simple_factor": sp.factor(simple),
        "accessory_constant": expected_accessory,
        "actual_saturation_gates": "ALL UNITS",
        "groebner": [sp.factor(poly.as_expr()) for poly in groebner.polys],
    }


def positive_control_2221() -> None:
    divisor = sp.expand(x**2 * (x - 1) ** 2 * (x - lam) ** 2 * (x - mu))
    critical = sp.expand(
        2 * (x - 1) * (x - lam) * (x - mu)
        + 2 * x * (x - lam) * (x - mu)
        + 2 * x * (x - 1) * (x - mu)
        + x * (x - 1) * (x - lam)
    )
    energy = critical / 7
    _, remainder = sp.div(divisor, energy**2, x)
    quotient, derivative_remainder = sp.div(sp.diff(remainder, x), energy, x)
    require(derivative_remainder == 0, "2221 control lost E | R'")
    equations = []
    for monomial in (x, 1):
        numerator = sp.together(
            sp.Poly(quotient, x).coeff_monomial(monomial)
        ).as_numer_denom()[0]
        equations.append(sp.primitive(sp.Poly(numerator, lam, mu))[1].as_expr())
    groebner = sp.groebner(equations, mu, lam, order="lex", domain=sp.QQ)
    relation = (
        2 * lam**5
        - 5 * lam**4
        - 6 * lam**3
        + 14 * lam**2
        + 6 * lam
        + mu
        - 6
    )
    modulus = (lam**3 - 2 * lam**2 - lam + 1) * (
        lam**3 - lam**2 - 2 * lam + 1
    )
    require(
        groebner.reduce(relation)[1] == 0
        and groebner.reduce(modulus)[1] == 0,
        "THM-2840 positive control differs",
    )


def nielsen_control() -> dict[tuple[int, ...], tuple[int, int, list[int]]]:
    total_cycle = tuple((point + 1) % 7 for point in range(7))
    rows = []
    unique_matchings = {
        tuple(sorted(tuple(sorted(edge)) for edge in matching))
        for matching in matchings(tuple(range(7)), 3)
    }
    require(len(unique_matchings) == 105, "wrong matching universe")
    for matching in unique_matchings:
        involution_list = list(range(7))
        for left, right in matching:
            involution_list[left], involution_list[right] = right, left
        involution = tuple(involution_list)
        pole_cycles = cycles(inverse(compose(involution, total_cycle)))
        if len(pole_cycles) != 4:
            continue
        passport = tuple(sorted((len(cycle) for cycle in pole_cycles), reverse=True))
        rows.append((involution, pole_cycles, passport))
    require(len(rows) == 35, "wrong planar universe")
    histogram: dict[tuple[int, ...], int] = {}
    for _involution, _pole_cycles, passport in rows:
        histogram[passport] = histogram.get(passport, 0) + 1
    require(
        histogram
        == {(4, 1, 1, 1): 7, (3, 2, 1, 1): 21, (2, 2, 2, 1): 7},
        "wrong planar passport histogram",
    )
    result = {}
    for target in ((4, 1, 1, 1), (3, 2, 1, 1), (2, 2, 2, 1)):
        marked = set()
        unmarked: dict[tuple[int, ...], tuple[int, ...]] = {}
        for involution, pole_cycles, passport in rows:
            if passport != target:
                continue
            unmarked.setdefault(rotation_key(involution), involution)
            for labelled_cycles in permutations(pole_cycles):
                if tuple(len(cycle) for cycle in labelled_cycles) == target:
                    marked.add(marked_key(involution, labelled_cycles))
        orders = [
            len(group_closure([involution, total_cycle]))
            for involution in unmarked.values()
        ]
        result[target] = (len(marked), len(unmarked), sorted(orders))
    return result


def main() -> None:
    print("THM-3123 exact companion")
    positive_control_2221()
    print("THM-2840 positive control: PASS")
    for passport in ((4, 1, 1, 1), (3, 2, 1, 1)):
        record = exact_accessory(passport)
        for key, value in record.items():
            print(f"{key}: {value}")
    explicit_s7_certificates()
    print("explicit S7 transposition certificates: PASS")
    nielsen = nielsen_control()
    expected_nielsen = {
        (4, 1, 1, 1): (6, 1, [5040]),
        (3, 2, 1, 1): (6, 3, [5040, 5040, 5040]),
        (2, 2, 2, 1): (6, 1, [14]),
    }
    require(nielsen == expected_nielsen, "Nielsen/monodromy census differs")
    print("nielsen:", nielsen)
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
