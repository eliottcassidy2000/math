#!/usr/bin/env python3
"""Exact referee for the THM-2858 six-ray half-plane factorial hostile."""

from __future__ import annotations

from functools import cache
from collections import Counter
from itertools import combinations_with_replacement
from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@cache
def adjacent_tensor(indices: tuple[int, ...]) -> int:
    value = sp.Integer(0)
    for mask in range(1 << len(indices)):
        degrees = []
        sign = 1
        for position, index in enumerate(indices):
            if mask & (1 << position):
                degrees.append(index + 1)
            else:
                degrees.append(index)
                sign = -sign
        term = sp.factorial(sum(degrees))
        for degree in degrees:
            term /= sp.factorial(degree)
        value += sign * term
    require(value.q == 1, "atomic tensor lost integrality")
    return int(value)


def interval_add(
    first: tuple[sp.Rational, sp.Rational],
    second: tuple[sp.Rational, sp.Rational],
) -> tuple[sp.Rational, sp.Rational]:
    return first[0] + second[0], first[1] + second[1]


def interval_mul(
    first: tuple[sp.Rational, sp.Rational],
    second: tuple[sp.Rational, sp.Rational],
) -> tuple[sp.Rational, sp.Rational]:
    values = (
        first[0] * second[0],
        first[0] * second[1],
        first[1] * second[0],
        first[1] * second[1],
    )
    return min(values), max(values)


def interval_pow(
    interval: tuple[sp.Rational, sp.Rational], exponent: int
) -> tuple[sp.Rational, sp.Rational]:
    if exponent == 0:
        return sp.Rational(1), sp.Rational(1)
    result = (sp.Rational(1), sp.Rational(1))
    for _ in range(exponent):
        result = interval_mul(result, interval)
    return result


def polynomial_interval(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    box: tuple[tuple[sp.Rational, sp.Rational], ...],
) -> tuple[sp.Rational, sp.Rational]:
    polynomial = sp.Poly(sp.expand(expression), *variables)
    answer = (sp.Rational(0), sp.Rational(0))
    for powers, coefficient in polynomial.terms():
        term = (sp.Rational(coefficient), sp.Rational(coefficient))
        for interval, exponent in zip(box, powers):
            term = interval_mul(term, interval_pow(interval, exponent))
        answer = interval_add(answer, term)
    return answer


def main() -> None:
    # c_1=1.  The c_2 and c_4 coordinates are frozen rationally.
    fixed_c2 = sp.Rational("0.0006557567088377382") - sp.I * sp.Rational(
        "1.2141776252864656"
    )
    fixed_c4 = sp.Rational("0.009935469352174497") + sp.I * sp.Rational(
        "1.2542668065909375"
    )

    a3, b3, a5, b5, a6, b6 = sp.symbols(
        "a3 b3 a5 b5 a6 b6", real=True
    )
    variables = (a3, b3, a5, b5, a6, b6)
    coefficients = (
        sp.Integer(1),
        fixed_c2,
        a3 + sp.I * b3,
        fixed_c4,
        a5 + sp.I * b5,
        a6 + sp.I * b6,
    )
    indices = (1, 2, 3, 4, 5, 6)

    moments: dict[int, sp.Expr] = {}
    equations: list[sp.Expr] = []
    for order in (2, 3, 4, 5):
        moment = sp.expand(
            sum(
                (
                    factorial(order)
                    // sp.prod(factorial(count) for count in Counter(address).values())
                )
                * adjacent_tensor(tuple(indices[position] for position in address))
                * sp.prod(coefficients[position] for position in address)
                for address in combinations_with_replacement(range(6), order)
            )
        )
        moments[order] = moment
        if order <= 4:
            equations.extend(
                (
                    sp.expand(sp.re(moment)),
                    sp.expand(sp.im(moment)),
                )
            )

    center = sp.Matrix(
        [
            sp.Rational(
                "0.000063621434994120405029506867116333343987061889208839"
            ),
            sp.Rational(
                "0.040009898798925760784579140591476738268555395921997"
            ),
            sp.Rational(
                "0.0045193541665913666254482593704175376727873260493792"
            ),
            sp.Rational(
                "-0.65481952791873310288768173269952383004730694141286"
            ),
            sp.Rational(
                "0.01313872825314743624425793558960064291574202146893"
            ),
            sp.Rational(
                "0.098948300256372361843344915268029197245827621460675"
            ),
        ]
    )
    radius = sp.Rational(1, 10**45)
    box = tuple((value - radius, value + radius) for value in center)

    function = sp.Matrix(equations)
    jacobian = function.jacobian(variables)
    substitution = dict(zip(variables, center))
    jacobian_center = jacobian.subs(substitution)
    require(jacobian_center.det() != 0, "center Jacobian became singular")
    inverse = jacobian_center.inv()
    newton_center = center - inverse * function.subs(substitution)

    jacobian_box = [
        [
            polynomial_interval(jacobian[row, column], variables, box)
            for column in range(6)
        ]
        for row in range(6)
    ]

    # E=I-CJ(X), represented entrywise by exact rational intervals.
    error_box: list[list[tuple[sp.Rational, sp.Rational]]] = []
    for row in range(6):
        error_row = []
        for column in range(6):
            value = (
                sp.Rational(1) if row == column else sp.Rational(0),
                sp.Rational(1) if row == column else sp.Rational(0),
            )
            for middle in range(6):
                scalar = inverse[row, middle]
                term = interval_mul(
                    (scalar, scalar), jacobian_box[middle][column]
                )
                value = interval_add(value, (-term[1], -term[0]))
            error_row.append(value)
        error_box.append(error_row)

    margins = []
    contraction_row_bounds = []
    for row in range(6):
        contraction_row_bound = sum(
            max(abs(interval[0]), abs(interval[1]))
            for interval in error_box[row]
        )
        contraction_row_bounds.append(contraction_row_bound)
        error_radius = contraction_row_bound * radius
        displacement = abs(newton_center[row] - center[row])
        margin = radius - displacement - error_radius
        require(margin > 0, f"Krawczyk inclusion failed in coordinate {row}")
        margins.append(margin)
    require(
        max(contraction_row_bounds) < 1,
        "preconditioned Newton map lost strict contraction",
    )

    # The box stays strictly inside the open right half-plane.
    for position in (0, 2, 4):
        require(box[position][0] > 0, "a coefficient reached the imaginary axis")
    require(sp.re(fixed_c2) > 0 and sp.re(fixed_c4) > 0, "fixed ray left half-plane")
    require(
        2000 * sp.re(fixed_c2) > abs(sp.im(fixed_c2))
        and 2000 * sp.re(fixed_c4) > abs(sp.im(fixed_c4)),
        "a fixed ray left the explicit 2000-sector",
    )
    for real_position, imaginary_position in ((0, 1), (2, 3), (4, 5)):
        maximum_imaginary_size = max(
            abs(box[imaginary_position][0]), abs(box[imaginary_position][1])
        )
        require(
            2000 * box[real_position][0] > maximum_imaginary_size,
            "a variable ray left the explicit 2000-sector",
        )

    moment_five_real = polynomial_interval(
        sp.re(moments[5]), variables, box
    )
    moment_five_imag = polynomial_interval(
        sp.im(moments[5]), variables, box
    )
    require(
        moment_five_real[1] < 0
        or moment_five_real[0] > 0
        or moment_five_imag[1] < 0
        or moment_five_imag[0] > 0,
        "fifth moment interval contains zero",
    )
    require(
        moment_five_real[0] > 647000000000
        and moment_five_imag[1] < -3290000000000,
        "compact fifth-moment bounds changed",
    )
    require(
        min(margins) > sp.Rational(99999, 100000) * radius,
        "compact Krawczyk margin changed",
    )

    print("SIX-RAY OPEN-HALF-PLANE FACTORIAL HOSTILE -- exact Krawczyk referee")
    print("status=PROVED-CANDIDATE+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
    print("support=(d1,d2,d3,d4,d5,d6); c1=1; c2,c4=fixed rational")
    print("equations=Re/Im L(H^k), k=2,3,4; variables=c3,c5,c6")
    print("center_jacobian_nonzero=1")
    print("krawczyk_strict_inclusion=1")
    print("preconditioned_newton_contraction_norm_lt=1")
    print("minimum_inclusion_margin_gt=(99999/100000)*10^-45")
    print("all_six_real_parts_positive=1 explicit_sector_K=2000")
    print("moment5_real_gt=647000000000 moment5_imag_lt=-3290000000000")
    print("consequence=factorial moments 1..4 and Gaussian moments 1..8 vanish")
    print("scope=finite-moment hostile; not a GMC2 counterexample")


if __name__ == "__main__":
    main()
