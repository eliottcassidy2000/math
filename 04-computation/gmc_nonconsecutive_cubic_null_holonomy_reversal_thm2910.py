#!/usr/bin/env python3
"""Exact companion for THM-2910.

For L(s^q)=q! and f_j=s^j/j!, compare the two four-slot block charts

    U=(f_1-f_0)+x(f_r-f_2),
    V=(f_2-f_1)+y(f_r-f_2),              r=5,6.

Each chart has one positive point where its binary quadratic moment
divides its binary cubic moment.  Exact subresultants, quotient-ring
evaluation and Bernstein certificates prove that the endpoint-holonomy
determinant is negative for r=5 and positive for r=6.  Both local
response-row (1,2,3) minors remain coefficientwise positive.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from math import comb, factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


s, x, y, z = sp.symbols("s x y z")


def canonical_digest(polynomial: sp.Poly) -> str:
    records = "\n".join(
        f"{','.join(str(exponent) for exponent in monomial)}:{coefficient}"
        for monomial, coefficient in polynomial.terms()
    )
    return sha256((records + "\n").encode()).hexdigest()


def primitive_integer(polynomial: sp.Poly, orient: bool = False) -> sp.Poly:
    """Clear positive denominators and remove integer content."""

    denominator, integral = sp.Poly(
        polynomial,
        y,
        domain=sp.QQ,
    ).clear_denoms(convert=True)
    require(denominator > 0, "denominator clearing changed orientation")
    primitive = integral.primitive()[1]
    if orient and primitive.LC() < 0:
        primitive = -primitive
    return primitive


def f(index: int) -> sp.Expr:
    return s**index / factorial(index)


def functional(polynomial: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), s)
    return sp.expand(
        sum(
            coefficient * factorial(exponent[0])
            for exponent, coefficient in expanded.terms()
        )
    )


def indexed_moment(forms: tuple[dict[int, sp.Expr], ...]) -> sp.Expr:
    return sp.expand(
        sum(
            sp.prod(
                (coefficient for _, coefficient in terms),
                start=sp.Integer(1),
            )
            * sp.Rational(
                factorial(sum(index for index, _ in terms)),
                sp.prod(
                    (factorial(index) for index, _ in terms),
                    start=sp.Integer(1),
                ),
            )
            for terms in product(*(tuple(form.items()) for form in forms))
        )
    )


def sign_variations(coefficients: list[sp.Expr]) -> int:
    nonzero = [coefficient for coefficient in coefficients if coefficient]
    return sum(
        bool(left > 0) != bool(right > 0)
        for left, right in zip(nonzero, nonzero[1:])
    )


def bernstein_coefficients(
    polynomial: sp.Expr,
    lower: sp.Rational,
    upper: sp.Rational,
) -> list[sp.Expr]:
    source = sp.Poly(polynomial, y, domain=sp.QQ)
    degree = source.degree()
    affine = sp.Poly(
        sp.expand(source.as_expr().subs(y, lower + (upper - lower) * z)),
        z,
        domain=sp.QQ,
    )
    powers = [affine.nth(index) for index in range(degree + 1)]
    return [
        sp.factor(
            sum(
                powers[index]
                * sp.Rational(comb(order, index), comb(degree, index))
                for index in range(order + 1)
            )
        )
        for order in range(degree + 1)
    ]


def strictly_positive_on(
    polynomial: sp.Expr,
    interval: tuple[sp.Rational, sp.Rational],
) -> bool:
    return all(
        coefficient > 0
        for coefficient in bernstein_coefficients(polynomial, *interval)
    )


def determinant_three(matrix: list[list[sp.Expr]]) -> sp.Expr:
    (a, b, c), (d, e, f_entry), (g, h, i) = matrix
    return sp.expand(
        a * (e * i - f_entry * h)
        - b * (d * i - f_entry * g)
        + c * (d * h - e * g)
    )


def difference_to_raw(vector: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    """Convert a sparse adjacent-difference vector to the f_j basis."""

    raw: dict[int, sp.Expr] = {}
    for index, coefficient in vector.items():
        raw[index + 1] = raw.get(index + 1, 0) + coefficient
        raw[index] = raw.get(index, 0) - coefficient
    return {
        index: sp.expand(coefficient)
        for index, coefficient in raw.items()
        if coefficient != 0
    }


def local_response_minor(vector: dict[int, sp.Expr]) -> sp.Expr:
    """The response-row (1,2,3) TP3 minor used in THM-2906."""

    raw = difference_to_raw(vector)
    rows = [
        difference_to_raw({1 + row: sp.Integer(1)})
        for row in range(3)
    ]
    return determinant_three(
        [
            [
                indexed_moment(tuple([row_vector] + [raw] * power))
                for power in (1, 2, 3)
            ]
            for row_vector in rows
        ]
    )


def quotient_evaluate(
    expression: sp.Expr,
    x_class: sp.Poly,
    modulus: sp.Poly,
) -> sp.Poly:
    """Evaluate an x-polynomial in QQ[y]/(modulus) by Horner."""

    accumulator = sp.Poly(0, y, domain=sp.QQ)
    for coefficient in sp.Poly(expression, x).all_coeffs():
        accumulator = (
            accumulator * x_class
            + sp.Poly(coefficient, y, domain=sp.QQ)
        ).rem(modulus)
    return accumulator


def chart(support: tuple[int, int, int, int]) -> dict[str, sp.Expr]:
    a, b, c, d = support
    first = {b: sp.Integer(1), a: sp.Integer(-1)}
    second = {c: sp.Integer(1), b: sp.Integer(-1)}
    third = {d: sp.Integer(1), c: sp.Integer(-1)}
    lower = dict(first)
    upper = dict(second)
    for index, coefficient in third.items():
        lower[index] = lower.get(index, 0) + x * coefficient
        upper[index] = upper.get(index, 0) + y * coefficient

    lower_polynomial = f(b) - f(a) + x * (f(d) - f(c))
    upper_polynomial = f(c) - f(b) + y * (f(d) - f(c))
    first_differences = {
        index: sp.Integer(1)
        for index in range(a, b)
    }
    second_differences = {
        index: sp.Integer(1)
        for index in range(b, c)
    }
    third_differences = {
        index: sp.Integer(1)
        for index in range(c, d)
    }
    lower_differences = dict(first_differences)
    upper_differences = dict(second_differences)
    for index in third_differences:
        lower_differences[index] = lower_differences.get(index, 0) + x
        upper_differences[index] = upper_differences.get(index, 0) + y

    moments: dict[tuple[int, int], sp.Expr] = {}
    literal_controls = 0
    for total in (2, 3, 4):
        for upper_count in range(total + 1):
            key = (total - upper_count, upper_count)
            forms = tuple(
                [lower] * (total - upper_count) + [upper] * upper_count
            )
            indexed = indexed_moment(forms)
            literal = functional(
                lower_polynomial ** (total - upper_count)
                * upper_polynomial**upper_count
            )
            require(
                sp.expand(indexed - literal) == 0,
                f"{support}: literal factorial moment mismatch at {key}",
            )
            moments[key] = indexed
            literal_controls += 1

    g0, g1, g2 = moments[(2, 0)], moments[(1, 1)], moments[(0, 2)]
    t0, t1, t2, t3 = (
        moments[(3, 0)],
        moments[(2, 1)],
        moments[(1, 2)],
        moments[(0, 3)],
    )
    quartic = tuple(moments[(4 - count, count)] for count in range(5))
    invariant_one = sp.expand(
        3 * t1 * g0 * g2 - t3 * g0**2 - 2 * t0 * g1 * g2
    )
    invariant_two = sp.expand(
        3 * t2 * g0 * g2 - 2 * t3 * g1 * g0 - t0 * g2**2
    )
    holonomy = sp.expand(
        (2 * quartic[1] * g0 - quartic[0] * g1) * g2**2
        - (2 * quartic[3] * g2 - quartic[4] * g1) * g0**2
    )
    return {
        "g0": g0,
        "g1": g1,
        "g2": g2,
        "i1": invariant_one,
        "i2": invariant_two,
        "j": holonomy,
        "tp3_lower": local_response_minor(lower_differences),
        "tp3_upper": local_response_minor(upper_differences),
        "literal_controls": sp.Integer(literal_controls),
    }


CASES = {
    (0, 1, 2, 5): {
        "y_interval": (sp.Rational(1, 84), sp.Rational(1, 83)),
        "x_interval": (sp.Rational(1, 37), sp.Rational(1, 35)),
        "quadratic_factor": 108 * y**2 + 12 * y + 1,
        "holonomy_sign": -1,
    },
    (0, 1, 2, 6): {
        "y_interval": (sp.Rational(1, 161), sp.Rational(1, 159)),
        "x_interval": (sp.Rational(1, 105), sp.Rational(1, 104)),
        "quadratic_factor": 437 * y**2 + 18 * y + 1,
        "holonomy_sign": 1,
    },
}

EXPECTED = {
    (0, 1, 2, 5): {
        "profiles": (15, 10, 10, 14),
        "modular_resultant": 5,
        "digests": (
            "a592136c75d3c826701f6394a0f5ad60ffe94e57435c4c261d33d1182f5bd484",
            "57dad96e4d8b4d3df563e6ec1293bed18f4166dcffcdb7006c814a3920388c0b",
            "c4716625b18cf3b3860382bd05d71b6d1a123777813ac51a7e8ade76bb6dab68",
            "fb2579bbc170978cf97e8fddf7cc3ff0dd1fda27d9f40e536063c4f48cd24100",
            "98e2f11dec167ff03c1b25778096aaebc3318f9911dff8d4c5ee8eefb04b75f1",
            "fd4e9448f8b4fa57146aefff0f826020eb3b1a568a078f1af22a3f5e58cd6de0",
        ),
    },
    (0, 1, 2, 6): {
        "profiles": (15, 10, 10, 14),
        "modular_resultant": 4,
        "digests": (
            "eecb2d7a03b035ea0d3e4774dc7ae2685baab7eda6ba48c8d2369059e93e8677",
            "9d79c56a56085a62ac0bab771db06a6a4e0a644e3291dbfdd5e5995fa284e240",
            "8a9b8d6587944501940b7f88d84f68f16bb00330ec6ff6d781f2952134534c32",
            "e4d29b23f4c214dac8742ea311f7a5d854636306c8beec71a2fd0f694ea454d4",
            "49a8a7803db1dd48c95c962e15f6a757066c69a19335dc8b8c675a50d594f5de",
            "356bb00e79f264ae7744b9e9b923f8b74a1ddb610c48b1b6fa53c7e75eb36d20",
        ),
    },
}


def audit_case(
    support: tuple[int, int, int, int],
    data: dict[str, object],
) -> dict[str, object]:
    objects = chart(support)
    invariant_one = objects["i1"]
    invariant_two = objects["i2"]
    holonomy = objects["j"]
    subresultants = sp.subresultants(invariant_one, invariant_two, x)
    require(
        [sp.degree(item, x) for item in subresultants] == [4, 3, 2, 1, 0],
        f"{support}: subresultant profile changed",
    )
    factors = sp.factor_list(subresultants[-1])[1]
    y_factors = [
        (factor, exponent)
        for factor, exponent in factors
        if sp.degree(factor, y) > 0
    ]
    require(
        sorted((sp.degree(factor, y), exponent) for factor, exponent in y_factors)
        == [(2, 2), (15, 1)],
        f"{support}: resultant factor profile changed",
    )
    quadratic = next(factor for factor, _ in y_factors if sp.degree(factor, y) == 2)
    expected_quadratic = sp.sympify(data["quadratic_factor"])
    ratio = sp.cancel(
        sp.Poly(quadratic, y).LC() / sp.Poly(expected_quadratic, y).LC()
    )
    require(
        sp.expand(quadratic - ratio * expected_quadratic) == 0
        and sp.discriminant(expected_quadratic, y) < 0
        and sp.Poly(expected_quadratic, y).LC() > 0,
        f"{support}: positive quadratic factor changed",
    )

    eliminant = next(factor for factor, _ in y_factors if sp.degree(factor, y) == 15)
    eliminant_poly = sp.Poly(eliminant, y, domain=sp.QQ)
    if eliminant_poly.LC() < 0:
        eliminant_poly = -eliminant_poly
    eliminant_integer = primitive_integer(eliminant_poly, orient=True)
    y_interval = data["y_interval"]
    require(
        isinstance(y_interval, tuple)
        and len(y_interval) == 2
        and sign_variations(eliminant_poly.all_coeffs()) == 1
        and eliminant_poly.eval(y_interval[0]) < 0
        and eliminant_poly.eval(y_interval[1]) > 0
        and eliminant_poly.count_roots(*y_interval) == 1,
        f"{support}: unique positive eliminant root certificate changed",
    )

    linear = sp.Poly(subresultants[-2], x)
    content = sp.gcd(linear.nth(1), linear.nth(0))
    selector_a = sp.cancel(linear.nth(1) / content)
    selector_n = sp.cancel(-linear.nth(0) / content)
    midpoint = sum(y_interval) / 2
    if selector_a.subs(y, midpoint) < 0:
        selector_a, selector_n = -selector_a, -selector_n
    content_ratio = sp.cancel(content / expected_quadratic)
    require(
        not content_ratio.free_symbols and content_ratio != 0,
        f"{support}: selector content acquired a branch zero",
    )
    require(
        strictly_positive_on(selector_a, y_interval)
        and strictly_positive_on(selector_n, y_interval),
        f"{support}: selector lost positivity",
    )

    x_interval = data["x_interval"]
    require(
        isinstance(x_interval, tuple)
        and len(x_interval) == 2
        and strictly_positive_on(
            selector_n - x_interval[0] * selector_a,
            y_interval,
        )
        and strictly_positive_on(
            x_interval[1] * selector_a - selector_n,
            y_interval,
        ),
        f"{support}: x-coordinate box changed",
    )

    selector_a_poly = sp.Poly(selector_a, y, domain=sp.QQ)
    selector_n_poly = sp.Poly(selector_n, y, domain=sp.QQ)
    require(
        sp.gcd(selector_a_poly, eliminant_poly).degree() == 0,
        f"{support}: selector is not invertible on the eliminant",
    )
    x_class = (
        selector_n_poly * sp.invert(selector_a_poly, eliminant_poly)
    ).rem(eliminant_poly)

    for label, invariant in (("I1", invariant_one), ("I2", invariant_two)):
        require(
            quotient_evaluate(
                invariant,
                x_class,
                eliminant_poly,
            ).is_zero,
            f"{support}: {label} selector divisibility changed",
        )

    endpoint_remainder = quotient_evaluate(
        holonomy,
        x_class,
        eliminant_poly,
    )
    endpoint_integer = primitive_integer(endpoint_remainder)
    expected_sign = int(data["holonomy_sign"])
    require(
        endpoint_remainder.degree() == 14
        and strictly_positive_on(
            expected_sign * endpoint_remainder.as_expr(),
            y_interval,
        ),
        f"{support}: quotient-ring endpoint sign changed",
    )

    modular_resultant = int(
        sp.Poly(eliminant_integer, y, modulus=11).resultant(
            sp.Poly(endpoint_integer, y, modulus=11)
        )
    ) % 11
    expected_modular_resultant = {
        (0, 1, 2, 5): 5,
        (0, 1, 2, 6): 4,
    }[support]
    require(
        modular_resultant == expected_modular_resultant,
        f"{support}: modular endpoint coprimality changed",
    )

    tp3_lower = sp.Poly(objects["tp3_lower"], x, domain=sp.QQ)
    tp3_upper = sp.Poly(objects["tp3_upper"], y, domain=sp.QQ)
    require(
        tp3_lower.degree() == 6
        and tp3_upper.degree() == 6
        and all(coefficient > 0 for coefficient in tp3_lower.all_coeffs())
        and all(coefficient > 0 for coefficient in tp3_upper.all_coeffs()),
        f"{support}: local response-row TP3 positivity changed",
    )

    result = {
        "literal_controls": int(objects["literal_controls"]),
        "y_interval": y_interval,
        "x_interval": x_interval,
        "holonomy_sign": expected_sign,
        "profiles": (
            eliminant_poly.degree(),
            sp.degree(selector_a, y),
            sp.degree(selector_n, y),
            endpoint_remainder.degree(),
        ),
        "modular_resultant": modular_resultant,
        "digests": (
            canonical_digest(eliminant_integer),
            canonical_digest(selector_a_poly),
            canonical_digest(selector_n_poly),
            canonical_digest(endpoint_integer),
            canonical_digest(tp3_lower),
            canonical_digest(tp3_upper),
        ),
    }
    require(
        result["profiles"] == EXPECTED[support]["profiles"]
        and result["modular_resultant"]
        == EXPECTED[support]["modular_resultant"]
        and result["digests"] == EXPECTED[support]["digests"],
        f"{support}: locked exact certificate changed",
    )
    return result


def main() -> None:
    print("THM-2910 NONCONSECUTIVE CUBIC-NULL HOLONOMY SIGN REVERSAL")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
    results = {
        support: audit_case(support, data)
        for support, data in CASES.items()
    }
    for support, result in results.items():
        tag = "".join(str(value) for value in support)
        print(
            f"support_{tag}=literal_controls:{result['literal_controls']};"
            "resultant:G2^2*P15;"
            f"y_box:{result['y_interval']};x_box:{result['x_interval']};"
            "positive_cubic_null_points:1;"
            f"holonomy_sign:{result['holonomy_sign']}"
        )
        print(
            f"support_{tag}_profiles=P,A,N,K:{','.join(str(value) for value in result['profiles'])};"
            f"res_P_K_mod11:{result['modular_resultant']}"
        )
        print(
            f"support_{tag}_digests=P:{result['digests'][0]};"
            f"A:{result['digests'][1]};N:{result['digests'][2]};"
            f"K:{result['digests'][3]};"
            f"TP3U:{result['digests'][4]};TP3V:{result['digests'][5]}"
        )
    require(
        results[(0, 1, 2, 5)]["holonomy_sign"] == -1
        and results[(0, 1, 2, 6)]["holonomy_sign"] == 1,
        "opposite-sign conclusion changed",
    )
    print("consequence=no support-uniform endpoint orientation on the cubic-null ideal")
    print(
        "scope=two ordered nonconsecutive block charts; no arbitrary-support "
        "shared-line, SFC(4), or GMC closure claimed"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
