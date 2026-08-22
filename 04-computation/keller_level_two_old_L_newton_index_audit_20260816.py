#!/usr/bin/env python3
"""Exact depth-two old-L Newton/index companion for THM-3537.

The script reconstructs the literal second-preimage x-eliminant on the
transverse line (a,b,c)=(2/27+t,1,1), computes its t-adic Newton polygon and
regular residual polynomials, and compares the resulting tame normalization
discriminant exponent with the power-order discriminant exponent.

The all-level numerical monodromy scout is deliberately not imported.  This
file proves only the exact level-two local statement and records no proposed
recurrence at greater depths.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
import pickle
from pathlib import Path

import sympy as sp


EXPECTED_H_RAW_SHA256 = (
    "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
)
EXPECTED_SEMANTIC_SHA256 = (
    "02cbf6a236a84ca153bea6fa60fb3c1267593b9e0694c0cddc5b1d7dc5bd8f59"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def l_poly(a: sp.Expr, b: sp.Expr, c: sp.Expr) -> sp.Expr:
    return 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2


def t_poly(b: sp.Expr, c: sp.Expr) -> sp.Expr:
    return 4 - 3 * b * c


def order_at(poly: sp.Expr, variable: sp.Symbol) -> int:
    expanded = sp.Poly(sp.expand(poly), variable, domain=sp.QQ)
    return min(
        monomial[0]
        for monomial, coefficient in expanded.terms()
        if coefficient != 0
    )


def initial_coefficient(poly: sp.Expr, variable: sp.Symbol) -> sp.Rational:
    valuation = order_at(poly, variable)
    return sp.Poly(sp.expand(poly), variable, domain=sp.QQ).coeff_monomial(
        variable**valuation
    )


def degree_nine_eliminant(
    a: sp.Symbol,
    w: sp.Symbol,
    xi: sp.Symbol,
) -> sp.Poly:
    """Reconstruct the saturated x_2 eliminant before taking the slice."""

    b = sp.Integer(1)
    c = sp.Integer(1)
    y_denominator = 12 * a * w**2 - b**2 * w**2 + b * w + 2
    middle_y = b - 3 * a * w * (9 * a * c * w - b * w + 2) / y_denominator
    middle_z = (2 * w - 3 * w**2 * middle_y - c) / w**3
    inner_core = (
        l_poly(w, middle_y, middle_z) * xi**3
        + t_poly(middle_y, middle_z) * xi
        - 2 * middle_z
    )
    inner_numerator = sp.together(inner_core).as_numer_denom()[0]
    outer_core = l_poly(a, b, c) * w**3 + t_poly(b, c) * w - 2 * c
    raw_resultant = sp.resultant(
        sp.expand(outer_core), sp.expand(inner_numerator), w
    )

    as_xi = sp.Poly(raw_resultant, xi)
    coefficient_gcd = sp.Integer(0)
    for coefficient in as_xi.all_coeffs():
        if coefficient != 0:
            coefficient_gcd = sp.gcd(
                coefficient_gcd, sp.Poly(coefficient, a).as_expr()
            )
    primitive = sp.Poly(sp.cancel(raw_resultant / coefficient_gcd), xi)
    rational_content = sp.gcd_list(
        [
            sp.Rational(value)
            for coefficient in primitive.all_coeffs()
            if coefficient != 0
            for value in sp.Poly(coefficient, a).coeffs()
        ]
    )
    return sp.Poly(sp.expand(primitive.as_expr() / rational_content), xi)


def cross(
    left: tuple[int, int], middle: tuple[int, int], right: tuple[int, int]
) -> int:
    return (
        (middle[0] - left[0]) * (right[1] - left[1])
        - (middle[1] - left[1]) * (right[0] - left[0])
    )


def lower_hull(points: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    hull: list[tuple[int, int]] = []
    for point in sorted(points):
        while len(hull) >= 2 and cross(hull[-2], hull[-1], point) <= 0:
            hull.pop()
        hull.append(point)
    return tuple(hull)


def slope(left: tuple[int, int], right: tuple[int, int]) -> Fraction:
    return Fraction(right[1] - left[1], right[0] - left[0])


def residual_polynomial(
    polynomial: sp.Poly,
    parameter: sp.Symbol,
    variable: sp.Symbol,
    left: tuple[int, int],
    right: tuple[int, int],
    residual_variable: sp.Symbol,
) -> sp.Poly:
    segment_slope = slope(left, right)
    denominator = segment_slope.denominator
    require((right[0] - left[0]) % denominator == 0,
            (left, right, denominator))
    residual = sp.Integer(0)
    for index in range(left[0], right[0] + 1, denominator):
        coefficient = polynomial.coeff_monomial(variable**index)
        if coefficient == 0:
            continue
        expected_order = left[1] + segment_slope * (index - left[0])
        require(expected_order.denominator == 1, (index, expected_order))
        if order_at(coefficient, parameter) == expected_order.numerator:
            residual += initial_coefficient(coefficient, parameter) * (
                residual_variable ** ((index - left[0]) // denominator)
            )
    return sp.Poly(sp.expand(residual), residual_variable, domain=sp.QQ)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    a, w, xi, parameter, residual_variable = sp.symbols("a w xi t Y")

    target_l = sp.factor(l_poly(a, 1, 1))
    sliced_l = sp.factor(target_l.subs(a, sp.Rational(2, 27) + parameter))
    require(sliced_l == parameter * (27 * parameter + 2), sliced_l)
    require(sp.diff(sliced_l, parameter).subs(parameter, 0) == 2,
            "nontransverse L slice")

    h_path = root / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
    h_raw = h_path.read_bytes()
    require(sha256(h_raw).hexdigest() == EXPECTED_H_RAW_SHA256, "H bytes")
    h_poly = pickle.loads(h_raw)
    h_variables = sorted(h_poly.free_symbols, key=lambda symbol: symbol.name)
    h_at_point = sp.factor(
        h_poly.subs(dict(zip(h_variables, (sp.Rational(2, 27), 1, 1))))
    )
    require(h_at_point == sp.Rational(61815040, 27), h_at_point)

    eliminant_a = degree_nine_eliminant(a, w, xi)
    eliminant = sp.Poly(
        sp.expand(eliminant_a.as_expr().subs(a, sp.Rational(2, 27) + parameter)),
        xi,
    )
    require(eliminant.degree() == 9, eliminant.degree())
    require(eliminant.LC().subs(parameter, 0) != 0, "nonunit leading coefficient")

    valuation_rows: list[tuple[int, int | None, str]] = []
    points: list[tuple[int, int]] = []
    for index in range(10):
        coefficient = eliminant.coeff_monomial(xi**index)
        if coefficient == 0:
            valuation_rows.append((index, None, "0"))
            continue
        valuation = order_at(coefficient, parameter)
        initial = initial_coefficient(coefficient, parameter)
        valuation_rows.append((index, valuation, str(initial)))
        points.append((index, valuation))

    expected_rows = (
        (0, 2, "-1792"),
        (1, 2, "-58304"),
        (2, 1, "-3584"),
        (3, 1, "1664"),
        (4, 1, "248320/9"),
        (5, 1, "-9931808/27"),
        (6, 0, "-28672"),
        (7, 0, "-101376"),
        (8, None, "0"),
        (9, 0, "-61815040/27"),
    )
    require(tuple(valuation_rows) == expected_rows, valuation_rows)

    hull = lower_hull(points)
    expected_hull = ((0, 2), (2, 1), (6, 0), (9, 0))
    require(hull == expected_hull, hull)
    slopes = tuple(slope(left, right) for left, right in zip(hull, hull[1:]))
    require(slopes == (Fraction(-1, 2), Fraction(-1, 4), Fraction(0)), slopes)

    residuals = tuple(
        residual_polynomial(
            eliminant,
            parameter,
            xi,
            left,
            right,
            residual_variable,
        )
        for left, right in zip(hull, hull[1:])
    )
    expected_residuals = (
        -1792 * (2 * residual_variable + 1),
        -3584 * (8 * residual_variable + 1),
        -sp.Rational(256, 27)
        * (241465 * residual_variable**3 + 10692 * residual_variable + 3024),
    )
    for actual, expected in zip(residuals, expected_residuals):
        require(sp.expand(actual.as_expr() - expected) == 0,
                (actual.as_expr(), expected))
        require(sp.gcd(actual, actual.diff()).degree() == 0,
                ("inseparable residual", actual.as_expr()))
    unit_residual_discriminant = sp.factor(
        sp.discriminant(residuals[2].as_expr(), residual_variable)
    )
    require(
        unit_residual_discriminant
        == -sp.Rational(3398871051182646293954560, 27),
        unit_residual_discriminant,
    )

    root_valuations = (Fraction(1, 2), Fraction(1, 4), Fraction(0))
    root_multiplicities = (2, 4, 3)
    ramification_cycles = (2, 4, 1, 1, 1)
    normalization_exponent = sum(length - 1 for length in ramification_cycles)
    require(normalization_exponent == 4, normalization_exponent)

    discriminant = sp.discriminant(eliminant.as_expr(), xi)
    order_discriminant_exponent = order_at(discriminant, parameter)
    require(order_discriminant_exponent == 8, order_discriminant_exponent)
    index_length = (order_discriminant_exponent - normalization_exponent) // 2
    require(
        order_discriminant_exponent
        == normalization_exponent + 2 * index_length
        and index_length == 2,
        index_length,
    )

    record = {
        "point": ["2/27", "1", "1"],
        "L_slice": str(sliced_l),
        "H_at_point": str(h_at_point),
        "degree": eliminant.degree(),
        "valuation_rows": valuation_rows,
        "hull": hull,
        "slopes": [str(value) for value in slopes],
        "residuals": [str(sp.factor(poly.as_expr())) for poly in residuals],
        "unit_residual_discriminant": str(unit_residual_discriminant),
        "root_valuations": [str(value) for value in root_valuations],
        "root_multiplicities": root_multiplicities,
        "ramification_cycles": ramification_cycles,
        "normalization_exponent": normalization_exponent,
        "order_discriminant_exponent": order_discriminant_exponent,
        "index_length": index_length,
    }
    semantic_sha256 = digest_json(record)

    print("== fixed Keller level-two old-L Newton/index audit ==")
    print(f"point=(2/27,1,1);L_slice={sliced_l};H_at_point={h_at_point}")
    print(f"degree={eliminant.degree()};valuation_rows={tuple(valuation_rows)}")
    print(f"newton_hull={hull};slopes={tuple(str(value) for value in slopes)}")
    print(
        "residuals="
        + str(tuple(str(sp.factor(poly.as_expr())) for poly in residuals))
    )
    print(f"unit_residual_discriminant={unit_residual_discriminant}")
    print(
        f"root_valuations={tuple(str(value) for value in root_valuations)};"
        f"multiplicities={root_multiplicities}"
    )
    print(
        f"tame_cycles={ramification_cycles};"
        f"normalization_discriminant_exponent={normalization_exponent}"
    )
    print(
        f"literal_x2_order_discriminant_exponent={order_discriminant_exponent};"
        f"local_index_length={index_length}"
    )
    print(f"semantic_sha256={semantic_sha256}")
    print(
        "scope=exact depth-two old-L canonical transverse local model;"
        "no all-level old-prime recurrence,newest-prime x-index,arbitrary-map,or JC claim"
    )
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic hash")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
